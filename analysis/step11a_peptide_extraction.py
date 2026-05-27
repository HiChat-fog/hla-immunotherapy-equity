import csv, json, os, re, time, requests

import os
os.makedirs("/home/mw/project", exist_ok=True)
MUT_PATH = "/home/mw/project/recurrent_mutations.json"
OUT_PATH = "/home/mw/project/recurrent_peptides.csv"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{}.fasta"

with open(MUT_PATH, "r") as f:
    recurrent = json.load(f)

print(f"复发突变: {len(recurrent)}")

DRIVERS = {"TP53","ARID1A","KRAS","PIK3CA","CDH1","ERBB2","RHOA",
           "SMAD4","CTNNB1","APC","PTEN","FBXW7","ERBB3","MET"}

selected = [r for r in recurrent if r["gene"] in DRIVERS or r["n_patients"] >= 5]
print(f"入选: {len(selected)} (驱动基因+高复发)")

def parse_pc(pc_str):
    if not pc_str or not pc_str.startswith("p."): return None,None,None
    pc = pc_str[2:]
    aa3to1 = {"Ala":"A","Cys":"C","Asp":"D","Glu":"E","Phe":"F","Gly":"G","His":"H",
              "Ile":"I","Lys":"K","Leu":"L","Met":"M","Asn":"N","Pro":"P","Gln":"Q",
              "Arg":"R","Ser":"S","Thr":"T","Val":"V","Trp":"W","Tyr":"Y","Ter":"*"}
    m = re.match(r"([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", pc)
    if m:
        return int(m.group(2)), aa3to1.get(m.group(1),"X"), aa3to1.get(m.group(3),"X")
    m = re.match(r"([A-Z*])(\d+)([A-Z*])", pc)
    if m: return int(m.group(2)), m.group(1), m.group(3)
    return None,None,None

def get_seq(swiss):
    try:
        resp = requests.get(UNIPROT_URL.format(swiss), timeout=30)
        if resp.status_code == 200:
            return "".join(l.strip() for l in resp.text.split("\n") if not l.startswith(">"))
    except: pass
    return None

def extract(seq, pos, ref_aa, var_aa):
    idx = pos - 1
    if idx < 0 or idx >= len(seq): return []
    if ref_aa != "*" and seq[idx] != ref_aa: return []
    peptides = []
    for length in [9,10,11]:
        for offset in range(length):
            start = idx - offset
            if start < 0: continue
            end = start + length
            if end > len(seq): break
            wt = seq[start:end]
            rel = idx - start
            if 0 <= rel < len(wt):
                mt = wt[:rel] + var_aa + wt[rel+1:]
                if len(mt) == length and mt != wt:
                    peptides.append(mt)
    return list(set(peptides))

protein_cache = {}
peptide_list = []

for i, r in enumerate(selected):
    swiss = r.get("swiss","")
    if not swiss or swiss == "null": continue
    pos, ref_aa, var_aa = parse_pc(r["protein_change"])
    if pos is None: continue

    if swiss not in protein_cache:
        print(f"  UniProt {swiss}...", end=" ", flush=True)
        seq = get_seq(swiss)
        if seq:
            protein_cache[swiss] = seq
            print(f"{len(seq)}aa")
        else:
            print("FAIL")
            continue
        time.sleep(0.1)

    seq = protein_cache[swiss]
    peptides = extract(seq, pos, ref_aa, var_aa)
    for pep in peptides:
        peptide_list.append({"gene":r["gene"],"protein_change":r["protein_change"],
            "peptide":pep,"n_patients":r["n_patients"],"swiss_id":swiss})

with open(OUT_PATH, "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["gene","protein_change","peptide","n_patients","swiss_id"])
    w.writeheader()
    w.writerows(peptide_list)

print(f"\n肽段提取完成: {len(peptide_list)} 条, 文件: {OUT_PATH}")
