"""
============================================================
  项目第13步：可视化整合
============================================================

  做什么：
  从已有数据文件生成报告所需的全部补充图表。
  - 图5：STAD 高频突变基因 Top 20
  - 图6：四 HLA 呈递效率对比
  - 图7：各 HLA 肽段 rank 分布
  - 图8：基因覆盖韦恩图
  - 图9：贪心覆盖对比（中国 vs 欧洲）
  - 图10：累积覆盖曲线
"""
import csv, os, sys
from collections import defaultdict, Counter

import numpy as np
from project_config import FIG_DIR, DRIVER_GENES

OUTPUT_DIR = os.path.dirname(__file__)

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    print("[OK] matplotlib")
except ImportError:
    print("[SKIP] matplotlib 未安装")
    sys.exit(0)

try:
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
except:
    pass

os.makedirs(FIG_DIR, exist_ok=True)

# ====== 加载数据 ======
mut_path = os.path.join(OUTPUT_DIR, "TCGA_STAD_mutations_firebrowse.csv")
pred_path = os.path.join(OUTPUT_DIR, "neoantigen_matrix_results.csv")
greedy_path = os.path.join(OUTPUT_DIR, "hla_greedy_priority.json")
freq_path = os.path.join(OUTPUT_DIR, "hla_frequency_data.csv")

# 图5：高频突变基因
print("生成图5：高频突变基因 Top 20...")
gene_nonsilent = Counter()
gene_patients = defaultdict(set)
non_silent = {'Missense_Mutation','Nonsense_Mutation','Frame_Shift_Del',
              'Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins',
              'Nonstop_Mutation','Splice_Site'}

with open(mut_path, "r", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        if row.get("Variant_Classification","") in non_silent:
            gene_nonsilent[row["Hugo_Symbol"]] += 1
            gene_patients[row["Hugo_Symbol"]].add(row["Tumor_Sample_Barcode"][:12])

top20 = gene_nonsilent.most_common(20)

fig, ax = plt.subplots(figsize=(12, 7))
genes = [g for g,_ in top20]
counts = [c for _,c in top20]
patients = [len(gene_patients[g]) for g in genes]
colors = ["#E74C3C" if g in DRIVER_GENES else "#3498DB" for g in genes]

bars = ax.bar(range(len(genes)), counts, color=colors, alpha=0.85)
for i, (g, c, p) in enumerate(zip(genes, counts, patients)):
    ax.text(i, c + 2, f"{c} ({p}人)", ha="center", fontsize=7, rotation=90)

ax.set_xticks(range(len(genes)))
ax.set_xticklabels(genes, rotation=45, ha="right", fontsize=9)
ax.set_ylabel("非沉默突变数", fontsize=12)
ax.set_title("TCGA-STAD 高频突变基因 Top 20（红=已知胃癌驱动基因）", fontsize=14, fontweight="bold")
legend = [Patch(color="#E74C3C", label="已知驱动基因"), Patch(color="#3498DB", label="其他")]
ax.legend(handles=legend, fontsize=10)
ax.grid(True, alpha=0.3, axis="y")
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "fig5_top_mutated_genes.png"), dpi=150)
plt.close()
print("  图5 → fig5_top_mutated_genes.png")

# 图6：四 HLA 呈递效率对比
print("生成图6：HLA 呈递效率对比...")
hla_stats = {}
with open(pred_path, "r", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        a = row["hla_allele"].replace("HLA-", "")
        if a not in hla_stats:
            hla_stats[a] = {"SB": 0, "WB": 0, "total": 0, "genes": set()}
        hla_stats[a]["total"] += 1
        if row["strength"] == "SB":
            hla_stats[a]["SB"] += 1
        if row["strength"] in ("SB", "WB"):
            hla_stats[a]["WB"] += 1
            hla_stats[a]["genes"].add(row["gene"])

alleles = ["A*02:01", "A*11:01", "A*24:02", "A*03:01"]
sb_vals = [hla_stats[a]["SB"] for a in alleles]
wb_vals = [hla_stats[a]["WB"] - hla_stats[a]["SB"] for a in alleles]
total_vals = [hla_stats[a]["total"] for a in alleles]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# 左：SB + WB 柱状图
x = np.arange(len(alleles))
ax1.bar(x, sb_vals, color="#E74C3C", label="强结合 (SB<0.5%)")
ax1.bar(x, wb_vals, bottom=sb_vals, color="#F39C12", label="弱结合 (WB<2%)")
ax1.set_xticks(x)
ax1.set_xticklabels(alleles, fontsize=11)
ax1.set_ylabel("呈递肽段数", fontsize=12)
ax1.set_title("各 HLA 呈递的胃癌新抗原数", fontsize=13, fontweight="bold")
ax1.legend(fontsize=9)

# 标注数字
for i, (sb, wb) in enumerate(zip(sb_vals, wb_vals)):
    ax1.text(i, sb + wb + 0.5, f"{sb+wb}", ha="center", fontsize=10, fontweight="bold")

# 右：覆盖基因数
gene_counts = [len(hla_stats[a]["genes"]) for a in alleles]
ax2.bar(alleles, gene_counts, color=["#E74C3C", "#F39C12", "#3498DB", "#2ECC71"], alpha=0.85)
for i, g in enumerate(gene_counts):
    ax2.text(i, g + 0.3, str(g), ha="center", fontsize=11, fontweight="bold")
ax2.set_ylabel("覆盖基因数", fontsize=12)
ax2.set_title("各 HLA 可呈递的独特基因数", fontsize=13, fontweight="bold")
ax2.grid(True, alpha=0.3, axis="y")
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "fig6_hla_presentation_comparison.png"), dpi=150)
plt.close()
print("  图6 → fig6_hla_presentation_comparison.png")

# 图7：各 HLA 肽段 rank 分布
print("生成图7：分布图...")
fig, ax = plt.subplots(figsize=(12, 6))
rank_data = {a: [] for a in alleles}
with open(pred_path, "r", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        a = row["hla_allele"].replace("HLA-", "")
        if a in rank_data:
            rank_data[a].append(float(row["percentile_rank"]))

positions = []
labels = []
for i, a in enumerate(alleles):
    ranks = rank_data[a]
    parts = ax.violinplot(ranks, [i], showmeans=True, showmedians=True)
    for pc in parts["bodies"]:
        pc.set_facecolor(["#E74C3C","#F39C12","#3498DB","#2ECC71"][i])
        pc.set_alpha(0.6)
    positions.append(i)
    labels.append(f"{a}\n(n={len(ranks)})")

ax.set_xticks(positions)
ax.set_xticklabels(labels, fontsize=10)
ax.set_ylabel("IEDB percentile rank (%)", fontsize=12)
ax.set_title("各 HLA 的肽段结合强度分布（越低表示结合越强）", fontsize=14, fontweight="bold")
ax.axhline(y=0.5, color="red", linestyle="--", alpha=0.5, label="SB 阈值 (0.5%)")
ax.axhline(y=2.0, color="orange", linestyle="--", alpha=0.5, label="WB 阈值 (2.0%)")
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3, axis="y")
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "fig7_rank_distribution.png"), dpi=150)
plt.close()
print("  图7 → fig7_rank_distribution.png")

# 图9：贪心覆盖对比
print("生成图9：贪心覆盖对比...")
import json
with open(greedy_path, "r") as f:
    greedy = json.load(f)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
for ax, (label, data) in zip([ax1, ax2], [("中国胃癌患者", greedy["chinese_gc"]), ("欧洲人群", greedy["european"])]):
    if not data:
        continue
    x = [d["rank"] for d in data]
    freq_only = [d["freq_pct"] for d in data]
    cum_cov = [d["cumulative_cov_pct"] for d in data]
    marginal = [d["marginal_cov_pct"] for d in data]

    ax.bar(x, freq_only, alpha=0.3, color="gray", label="该 HLA 频率")
    ax.bar(x, marginal, alpha=0.8, color=["#E74C3C","#F39C12","#3498DB","#2ECC71"][:len(x)],
           label="新增覆盖")
    ax.plot(x, cum_cov, "o-", color="black", linewidth=2, markersize=8, label="累计覆盖")
    for i, d in enumerate(data):
        ax.text(d["rank"], d["cumulative_cov_pct"] + 1, f"{d['cumulative_cov_pct']:.1f}%",
                ha="center", fontsize=10, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([d["allele"] for d in data])
    ax.set_ylabel("覆盖率 (%)", fontsize=12)
    ax.set_title(label, fontsize=13, fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis="y")

plt.suptitle("贪心覆盖优化：优先开发哪些 HLA 靶点？", fontsize=14, fontweight="bold")
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "fig9_greedy_coverage.png"), dpi=150)
plt.close()
print("  图9 → fig9_greedy_coverage.png")

# 图10：中欧最优覆盖对比汇总
print("生成图10：中欧对比汇总...")
fig, ax = plt.subplots(figsize=(10, 6))
cn_cov = greedy["chinese_gc"][-1]["cumulative_cov_pct"] if greedy["chinese_gc"] else 0
eu_cov = greedy["european"][-1]["cumulative_cov_pct"] if greedy["european"] else 0
fda_cn = 5.0
fda_eu = 27.1

x_labels = ["FDA 现状\n(A*02:01)", "贪心优化\n(前3 HLA)"]
cn_vals = [fda_cn, cn_cov]
eu_vals = [fda_eu, eu_cov]

x = np.arange(len(x_labels))
w = 0.35
ax.bar(x - w/2, cn_vals, w, color="#E74C3C", alpha=0.85, label="中国胃癌患者")
ax.bar(x + w/2, eu_vals, w, color="#3498DB", alpha=0.85, label="欧洲人群")

for i, (cv, ev) in enumerate(zip(cn_vals, eu_vals)):
    ax.text(i - w/2, cv + 1, f"{cv:.1f}%", ha="center", fontsize=11, fontweight="bold")
    ax.text(i + w/2, ev + 1, f"{ev:.1f}%", ha="center", fontsize=11, fontweight="bold")

ax.set_xticks(x)
ax.set_xticklabels(x_labels, fontsize=12)
ax.set_ylabel("患者覆盖率 (%)", fontsize=12)
ax.set_title("HLA 限制型免疫治疗的患者覆盖率：现状 vs 优化", fontsize=14, fontweight="bold")
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3, axis="y")
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "fig10_coverage_comparison.png"), dpi=150)
plt.close()
print("  图10 → fig10_coverage_comparison.png")

print(f"\n第13步完成！图表保存在 {FIG_DIR}/")
