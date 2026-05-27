import csv, json, os
from collections import defaultdict

MUT_PATH = "/home/mw/input/TCGA_STADfireb49124912/TCGA_STAD_mutations_firebrowse.csv"
import os
os.makedirs("/home/mw/project", exist_ok=True)
OUT_PATH = "/home/mw/project/recurrent_mutations.json"

print("筛选复发突变 (>=2 患者)...")

recurrent = defaultdict(set)
details = {}

with open(MUT_PATH, "r", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        gene = row["Hugo_Symbol"]
        vc = row.get("Variant_Classification","")
        protein = row.get("Protein_Change","")
        swiss = row.get("SwissProt_entry_Id","")
        barcode = row["Tumor_Sample_Barcode"]
        patient = barcode[:12]

        if vc != "Missense_Mutation": continue
        if not protein or not protein.startswith("p."): continue
        if not swiss or swiss == "null" or swiss == "": continue

        key = (gene, protein)
        recurrent[key].add(patient)
        if key not in details:
            details[key] = {"swiss": swiss}

result = []
for (gene, protein), patients in recurrent.items():
    if len(patients) >= 2:
        result.append({
            "gene": gene,
            "protein_change": protein,
            "patients": sorted(patients),
            "n_patients": len(patients),
            "swiss": details[(gene, protein)]["swiss"]
        })

result.sort(key=lambda x: -x["n_patients"])

with open(OUT_PATH, "w", encoding="utf-8") as f:
    json.dump(result, f, indent=2, ensure_ascii=False)

print(f"复发突变: {len(result)} 个 (>=2 患者)")
print(f"已保存: {OUT_PATH}")
