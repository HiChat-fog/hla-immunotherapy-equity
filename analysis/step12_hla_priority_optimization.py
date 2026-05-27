"""
============================================================
  项目第12步：HLA 靶点优先级优化（贪心覆盖算法）
============================================================

  输入：第11步的全矩阵IEDB预测结果
  输出：如果给你N个药物研发名额，该优先靶向哪些HLA等位基因？

  核心指标：每个HLA等位基因的综合"价值得分"
    = 呈递效率 × 人群频率 × 疾病负担权重

  算法：贪心覆盖 ——
    1. 初始：所有胃癌患者都未被覆盖
    2. 每次选"覆盖最多剩余患者"的HLA
    3. 迭代直到覆盖饱和
    4. 输出优先级排序列表
"""
import csv, os, json
from collections import defaultdict

OUTPUT_DIR = os.path.dirname(__file__)

# ====== 第1步：加载全矩阵预测结果 ======
print("=" * 60)
print("第1步：加载 IEDB 全矩阵预测结果")
print("=" * 60)

pred_path = os.path.join(OUTPUT_DIR, "neoantigen_matrix_results.csv")
hla_peptides = defaultdict(set)      # HLA → {呈递的肽段}
hla_sb_peptides = defaultdict(set)   # HLA → {强结合肽段}
hla_genes = defaultdict(set)         # HLA → {覆盖的基因}

with open(pred_path, "r", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        allele = row["hla_allele"].replace("HLA-", "")
        peptide = row["peptide"]
        gene = row["gene"]
        strength = row["strength"]

        if strength in ("SB", "WB"):       # 可呈递（弱结合及以上）
            hla_peptides[allele].add(peptide)
            hla_genes[allele].add(gene)
        if strength == "SB":               # 强结合
            hla_sb_peptides[allele].add(peptide)

for a in sorted(hla_peptides.keys()):
    print(f"  {a}: {len(hla_peptides[a])} 肽段 ({len(hla_sb_peptides[a])} SB), "
          f"{len(hla_genes[a])} 个基因")

# ====== 第2步：加载 HLA 人群频率 ======
print("\n" + "=" * 60)
print("第2步：加载 HLA 人群频率数据")
print("=" * 60)

freq_path = os.path.join(OUTPUT_DIR, "hla_frequency_data.csv")
freq_data = {}  # (allele, population) → frequency_pct

with open(freq_path, "r", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        key = (row["allele"], row["population"])
        freq_data[key] = float(row["frequency_pct"])

# 关键人群的频率（用于加权）
# 东亚胃癌患者: Zhou 2019
# 东亚一般人群: NMDP East Asian
# 欧洲一般人群: NMDP European
POPULATIONS = {
    "chinese_gc": "Chinese GC patients",    # 中国胃癌患者实测 (Zhou 2019)
    "east_asian": "East Asian (NMDP)",       # 东亚一般人群 (NMDP)
    "european": "European (median)",         # 欧洲中位频率 (AFND)
}

for allele in ["A*02:01", "A*11:01", "A*24:02"]:
    for pop_key, pop_name in POPULATIONS.items():
        f = freq_data.get((allele, pop_name))
        if f:
            print(f"  {allele} @ {pop_key}: {f}%")

# ====== 第3步：计算综合价值得分 ======
print("\n" + "=" * 60)
print("第3步：计算 HLA 综合价值得分")
print("=" * 60)

HLA_LIST = sorted(hla_peptides.keys())

def calc_scores(hla_list, freq_col="chinese_gc"):
    """计算每个HLA的加权呈递得分"""
    scores = []
    pop_name = POPULATIONS[freq_col]
    for allele in hla_list:
        n_peptides = len(hla_peptides[allele])
        n_sb = len(hla_sb_peptides[allele])
        n_genes = len(hla_genes[allele])
        freq = freq_data.get((allele, pop_name), 0)

        # 综合得分 = 呈递的肽段数 × 人群频率(%)
        # 这里肽段数 = 可呈递的"公共新抗原"数量
        score_raw = n_peptides * freq / 100.0
        score_sb = n_sb * freq / 100.0

        scores.append({
            "allele": allele,
            "peptides_presented": n_peptides,
            "strong_binders": n_sb,
            "genes_covered": n_genes,
            "frequency_pct": freq,
            "population": pop_name,
            "raw_score": round(score_raw, 1),
            "sb_score": round(score_sb, 1),
        })
    return sorted(scores, key=lambda x: -x["raw_score"])

print("对中国胃癌患者的综合价值排序：")
scores_cn = calc_scores(HLA_LIST, "chinese_gc")
for i, s in enumerate(scores_cn):
    print(f"  {i+1}. {s['allele']}: 呈递{s['peptides_presented']}肽段 × "
          f"频率{s['frequency_pct']}% = 得分{s['raw_score']}")

print("\n对欧洲人群的综合价值排序：")
scores_eu = calc_scores(HLA_LIST, "european")
for i, s in enumerate(scores_eu):
    print(f"  {i+1}. {s['allele']}: 呈递{s['peptides_presented']}肽段 × "
          f"频率{s['frequency_pct']}% = 得分{s['raw_score']}")

# ====== 第4步：贪心覆盖优化 ======
print("\n" + "=" * 60)
print("第4步：贪心覆盖 - 优先开发哪些HLA能最快覆盖最多患者")
print("=" * 60)

# 假设场景：中国胃癌患者群体
# 每个HLA覆盖的患者比例 ∝ 该等位基因的频率
# （简化：不考虑连锁不平衡，假设等位基因独立）

def greedy_coverage(hla_list, freq_col, coverage_label):
    """贪心选择HLA靶点，最大化累积患者覆盖"""
    pop = POPULATIONS[freq_col]

    # 各HLA的频率（作为独立覆盖率的近似）
    available = {}
    for allele in hla_list:
        f = freq_data.get((allele, pop), 0)
        if f > 0:
            # 每个HLA的"覆盖价值" = 能呈递的肽段数 × 频率
            available[allele] = {
                "freq": f,
                "peptides": len(hla_peptides[allele]),
                "genes": len(hla_genes[allele]),
                "sb": len(hla_sb_peptides[allele]),
            }

    # 贪心迭代
    selected = []
    # 简化的覆盖率模型：P(至少一个HLA覆盖) = 1 - Π(1 - freq_i)
    # 这假设等位基因独立（实际上A*02:01和A*11:01确实在不同基因座）
    cumulative_cov = 0.0
    uncovered = 1.0

    while available and len(selected) < 8:
        best = None
        best_gain = 0
        for a, info in available.items():
            gain = uncovered * (info["freq"] / 100.0) * info["peptides"]
            if gain > best_gain:
                best_gain = gain
                best = a

        if best is None:
            break

        info = available.pop(best)
        new_cov = info["freq"] / 100.0 * uncovered
        cumulative_cov += new_cov
        uncovered *= (1 - info["freq"] / 100.0)

        selected.append({
            "rank": len(selected) + 1,
            "allele": best,
            "freq_pct": info["freq"],
            "peptides": info["peptides"],
            "genes": info["genes"],
            "sb": info["sb"],
            "marginal_cov_pct": round(new_cov * 100, 1),
            "cumulative_cov_pct": round(cumulative_cov * 100, 1),
        })

        print(f"  第{len(selected)}优先: {best} (频率{info['freq']}%) "
              f"→ 新增覆盖{new_cov*100:.1f}%, 累计{cumulative_cov*100:.1f}%")

    return selected

print("\n--- 对中国胃癌患者 ---")
result_cn = greedy_coverage(HLA_LIST, "chinese_gc", "中国胃癌患者")

print("\n--- 对欧洲人群 ---")
result_eu = greedy_coverage(HLA_LIST, "european", "欧洲人群")

# ====== 第5步：对比分析——"治疗公平性差距" ======
print("\n" + "=" * 60)
print("第5步：治疗公平性差距量化")
print("=" * 60)

# FDA 现状：只有 A*02:01
fda_allele = "A*02:01"
cn_freq = freq_data.get((fda_allele, POPULATIONS["chinese_gc"]), 0)
eu_freq = freq_data.get((fda_allele, POPULATIONS["european"]), 0)

# 如果按中国患者需求重新排优先级，应该优先开发的目标
top3_cn = [r["allele"] for r in result_cn[:3]]
top3_cov = result_cn[2]["cumulative_cov_pct"] if len(result_cn) >= 3 else result_cn[-1]["cumulative_cov_pct"]

print(f"""
  当前 FDA 批准的 HLA 限制型疗法：
    靶向: {fda_allele}
    中国胃癌患者覆盖: ~{cn_freq}%
    欧洲患者覆盖: ~{eu_freq}%
    公平性差距: {eu_freq/cn_freq:.1f}倍

  如果按胃癌新抗原呈递效率重新排优先级：
    中国患者优先开发: {', '.join(top3_cn)}
    前3位累计覆盖: {top3_cov:.1f}%
    → 是FDA现状({cn_freq}%)的 {top3_cov/cn_freq:.1f} 倍

  对比：欧洲患者优先级（FDA现状已覆盖最好的）
    优先开发: {', '.join(r['allele'] for r in result_eu[:3])}
""")

# ====== 第6步：保存结果 ======
print("=" * 60)
print("第6步：保存优化结果")
print("=" * 60)

# 综合得分
score_path = os.path.join(OUTPUT_DIR, "hla_priority_scores.csv")
with open(score_path, "w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=[
        "allele", "peptides_presented", "strong_binders", "genes_covered",
        "frequency_pct", "population", "raw_score", "sb_score"])
    writer.writeheader()
    writer.writerows(scores_cn)
print(f"已保存: {score_path}")

# 贪心排序
greedy_path = os.path.join(OUTPUT_DIR, "hla_greedy_priority.json")
with open(greedy_path, "w", encoding="utf-8") as f:
    json.dump({"chinese_gc": result_cn, "european": result_eu}, f, indent=2)
print(f"已保存: {greedy_path}")

print("\n" + "=" * 60)
print("第12步完成！")
print("=" * 60)
