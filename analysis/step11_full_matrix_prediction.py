"""
第11步：全矩阵新抗原-HLA结合预测
文件日志版 — 进度写入 step11_progress.txt
"""
import csv, os, time, requests, datetime

OUTPUT_DIR = os.path.dirname(__file__)
OUT_PATH = os.path.join(OUTPUT_DIR, "neoantigen_matrix_results.csv")
PROG_PATH = os.path.join(OUTPUT_DIR, "step11_progress.txt")
HLA_ALLELES = ["HLA-A*02:01","HLA-A*11:01","HLA-A*24:02","HLA-A*03:01"]
IEDB_URL = "http://tools-cluster-interface.iedb.org/tools_api/mhci/"

def log(msg):
    with open(PROG_PATH, "a", encoding="utf-8") as f:
        f.write(f"[{datetime.datetime.now().strftime('%H:%M:%S')}] {msg}\n")

def load_peptides():
    pep_path = os.path.join(OUTPUT_DIR, "recurrent_peptides.csv")
    plist = []
    with open(pep_path, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            plist.append(row)
    return plist

def load_done():
    done = set()
    if os.path.exists(OUT_PATH):
        with open(OUT_PATH, "r", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                done.add((row["peptide"], row["hla_allele"]))
    return done

def predict_one(peptide, allele):
    payload = {
        "method": "recommended",
        "sequence_text": peptide,
        "allele": allele,
        "length": [len(peptide)],
        "user_tool": "hla_equity_stad",
    }
    resp = requests.post(IEDB_URL, data=payload, timeout=60)
    if resp.status_code != 200 or not resp.text.strip():
        return None
    lines = resp.text.strip().split("\n")
    best_rank, best_ic50 = 100.0, 999999.0
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) >= 10:
            try:
                rank = float(parts[9])
                score = float(parts[8])
                if rank < best_rank:
                    best_rank = rank
                    best_ic50 = 50000 ** (1 - score) if 0 < score < 1 else 999999.0
            except (ValueError, IndexError):
                continue
    if best_rank < 0.5:   strength = "SB"
    elif best_rank < 2.0: strength = "WB"
    elif best_rank < 10.0: strength = "WB_edge"
    else:                 strength = "non"
    return {
        "peptide": peptide, "gene": gene_map.get(peptide, "?"),
        "hla_allele": allele, "percentile_rank": round(best_rank, 2),
        "ic50_nM": round(best_ic50, 1), "strength": strength,
    }

def save_results(results):
    with open(OUT_PATH, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "peptide","gene","hla_allele","percentile_rank","ic50_nM","strength"])
        writer.writeheader()
        writer.writerows(results)

# ======== MAIN ========
log("=== 第11步启动 ===")

peptide_list = load_peptides()
unique_peptides = sorted(set(p["peptide"] for p in peptide_list))

gene_map = {}
for p in peptide_list:
    g, pid = p["gene"], p["peptide"]
    cur = gene_map.get(pid, "")
    gene_map[pid] = cur if g in cur else (cur + "/" + g).strip("/")

done = load_done()
pending = [(p, a) for p in unique_peptides for a in HLA_ALLELES if (p, a) not in done]

log(f"肽段:{len(unique_peptides)} 已存:{len(done)} 待处理:{len(pending)}")

# 加载已有结果
existing = []
if os.path.exists(OUT_PATH):
    with open(OUT_PATH, "r", encoding="utf-8") as f:
        existing = list(csv.DictReader(f))

results = existing.copy()
t0 = time.time()
fails = 0
last_save = len(results)

for i, (peptide, allele) in enumerate(pending):
    r = predict_one(peptide, allele)
    if r:
        results.append(r)
    else:
        fails += 1

    # 每50次写进度
    if (i + 1) % 50 == 0:
        elapsed = time.time() - t0
        rate = (i + 1) / elapsed * 60 if elapsed > 0 else 0
        eta = (len(pending) - i - 1) / rate if rate > 0 else 0
        log(f"[{i+1}/{len(pending)}] 成功:{len(results)-last_save+50} 失败:{fails} 速率:{rate:.0f}/min ETA:{eta:.0f}min")

    # 每200条新结果保存
    if len(results) - last_save >= 200:
        save_results(results)
        last_save = len(results)
        log(f"  已保存 {len(results)} 条")

    time.sleep(0.15)

# 最终保存
save_results(results)
elapsed = time.time() - t0

# 汇总
hla_stats = {}
for a in HLA_ALLELES:
    s = a.replace("HLA-", "")
    ar = [r for r in results if r["hla_allele"] == a]
    sb = sum(1 for r in ar if r["strength"] == "SB")
    wb = sum(1 for r in ar if r["strength"] in ("SB","WB"))
    genes = set(r["gene"] for r in ar if r["strength"] in ("SB","WB"))
    hla_stats[s] = {"sb":sb, "wb":wb, "total":len(ar), "genes":len(genes)}

log(f"=== 完成! 耗时:{elapsed/60:.1f}min 总计:{len(results)}条 ===")
for a in ["A*02:01","A*11:01","A*24:02","A*03:01"]:
    d = hla_stats[a]
    log(f"  {a}: SB={d['sb']} WB={d['wb']} Total={d['total']} Genes={d['genes']}")

# 保存汇总
sp = os.path.join(OUTPUT_DIR, "hla_presentation_summary.csv")
with open(sp, "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f)
    w.writerow(["hla_allele","strong_binders","weak_binders","total","unique_genes"])
    for a in ["A*02:01","A*11:01","A*24:02","A*03:01"]:
        d = hla_stats[a]
        w.writerow([a, d["sb"], d["wb"], d["total"], d["genes"]])
log(f"汇总已保存: {sp}")
