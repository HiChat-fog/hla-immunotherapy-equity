"""
============================================================
  项目第5步：TME 免疫细胞组成分析
============================================================

  做什么：
  1. 加载亚洲和白人胃癌表达数据
  2. 用已知的免疫细胞 marker 基因，给每个样本算"免疫评分"
  3. 比较两组在各个免疫细胞类型上的差异
  4. 画图可视化

  原理：
  每种免疫细胞有自己独特的"特征基因"——这些基因在这种细胞里
  表达量特别高。我们取这些基因在样本中的平均表达值，作为这种
  免疫细胞在此样本中的"富集分数"。

  分数越高 → 这种免疫细胞在这个肿瘤中的占比可能越大。
============================================================
"""
import csv
import os
import math
from collections import Counter

OUTPUT_DIR = os.path.dirname(__file__)


# ====== 第1步：加载表达数据 ======
print("=" * 60)
print("第1步：加载表达矩阵")
print("=" * 60)

def load_expression(filepath):
    """
    加载 tsv 格式的表达矩阵。
    返回：
      genes: 基因名列表
      samples: 样本ID列表
      matrix: {gene_name: [value1, value2, ...]} 字典
    """
    with open(filepath, "r") as f:
        header = f.readline().strip().split("\t")
        samples = header[1:]   # 第一列是 'sample'（基因名列）

        genes = []
        matrix = {}
        for line in f:
            parts = line.strip().split("\t")
            gene = parts[0]
            values = [float(v) for v in parts[1:]]
            genes.append(gene)
            matrix[gene] = values

    return genes, samples, matrix


asian_genes, asian_samples, asian_matrix = load_expression(
    os.path.join(OUTPUT_DIR, "asian_expression.tsv")
)
white_genes, white_samples, white_matrix = load_expression(
    os.path.join(OUTPUT_DIR, "white_expression.tsv")
)

print(f"亚洲组: {len(asian_genes)} 基因 x {len(asian_samples)} 样本")
print(f"白人组: {len(white_genes)} 基因 x {len(white_samples)} 样本")

# 验证两组基因列表一致
if asian_genes != white_genes:
    # 找交集
    common = set(asian_genes) & set(white_genes)
    print(f"[注意] 基因列表不一致，取交集: {len(common)} 个共同基因")
else:
    print("[通过] 两组基因列表完全一致")


# ====== 第2步：定义免疫细胞 marker 基因集 ======
print("\n" + "=" * 60)
print("第2步：加载免疫细胞特征基因集")
print("=" * 60)

# 这些 marker 基因来自已发表的免疫学文献
# 来源: Bindea et al., Immunity (2013); Charoentong et al., Cell Reports (2017)
# 每个基因都是该免疫细胞类型中"表达量显著高于其他细胞"的标志基因
IMMUNE_MARKERS = {
    # 杀伤性 T 细胞（CTL）——免疫攻击的主力
    "CD8+ T cells": [
        "CD8A", "CD8B", "GZMA", "GZMB", "GZMK",
        "PRF1", "GNLY", "IFNG", "NKG7", "CST7",
    ],

    # 辅助 T 细胞 —— 协调免疫反应
    "CD4+ T cells": [
        "CD4", "IL7R", "CD3D", "CD3E", "CD2",
        "TRBC1", "TRBC2", "TRAC", "ICOS", "MAL",
    ],

    # 调节 T 细胞（Treg）—— 免疫抑制，通常预后不良
    "Treg cells": [
        "FOXP3", "IL2RA", "CTLA4", "TNFRSF18", "TNFRSF4",
        "TNFRSF9", "IKZF2", "IKZF4", "BATF", "CCR8",
    ],

    # B 细胞 —— 抗体产生，与免疫治疗响应有关
    "B cells": [
        "CD19", "CD79A", "CD79B", "MS4A1", "BLK",
        "PAX5", "BANK1", "CD22", "CD40", "TNFRSF17",
    ],

    # NK 细胞 —— 天然杀伤细胞，第一道免疫防线
    "NK cells": [
        "NKG7", "GNLY", "KLRD1", "KLRF1", "KLRB1",
        "NCR1", "NCR3", "PRF1", "GZMB", "CD160",
    ],

    # 巨噬细胞 M1 —— 促炎，抗肿瘤
    "M1 Macrophages": [
        "IL12A", "IL12B", "IL23A", "TNF", "IL1B",
        "CXCL9", "CXCL10", "CXCL11", "IRF5", "PTGS2",
    ],

    # 巨噬细胞 M2 —— 免疫抑制，促肿瘤
    "M2 Macrophages": [
        "CD163", "CD206", "MSR1", "IL10", "TGFB1",
        "CCL18", "CCL22", "ARG1", "CHI3L1", "F13A1",
    ],

    # 中性粒细胞 —— 与预后不良相关
    "Neutrophils": [
        "ELANE", "MPO", "CEACAM8", "FPR1", "FCGR3B",
        "CXCR1", "CXCR2", "ITGAM", "CSF3R", "IL8",
    ],

    # 树突状细胞 —— 抗原呈递，启动适应性免疫
    "Dendritic cells": [
        "CD1C", "CLEC4C", "NRP1", "ITGAX", "CD83",
        "CCR7", "LAMP3", "BATF3", "IRF8", "FLT3",
    ],
}

# 第3步：计算每个样本的免疫细胞评分
print("\n" + "=" * 60)
print("第3步：计算免疫细胞评分")
print("=" * 60)


def calculate_immune_scores(samples, matrix):
    """
    对每个样本，计算它在每种免疫细胞类型上的"特征基因平均表达值"。

    参数:
      samples: 样本ID列表
      matrix:  {gene_name: [value_across_samples]}

    返回:
      scores: {cell_type: [score_for_each_sample]}
    """
    scores = {}
    total_markers = 0
    found_markers = 0

    for cell_type, marker_genes in IMMUNE_MARKERS.items():
        # 检查这些 marker 基因在我们的数据里有多少
        present_genes = [g for g in marker_genes if g in matrix]
        found_markers += len(present_genes)
        total_markers += len(marker_genes)

        # 找出这些基因在矩阵中的位置
        gene_indices = []
        for g in present_genes:
            # matrix[g] 是一个列表，长度 = 样本数
            # 索引 i 对应 samples[i]
            gene_indices.append(g)

        # 对每个样本，计算这些基因的平均表达值
        cell_scores = []
        for sample_idx in range(len(samples)):
            gene_values = []
            for g in gene_indices:
                gene_values.append(matrix[g][sample_idx])
            if gene_values:
                # 用中位数而不是均值——因为个别极高表达的基因会拉偏均值
                sorted_vals = sorted(gene_values)
                n = len(sorted_vals)
                median = sorted_vals[n // 2] if n % 2 == 1 else \
                         (sorted_vals[n // 2 - 1] + sorted_vals[n // 2]) / 2
                cell_scores.append(median)
            else:
                cell_scores.append(0)

        scores[cell_type] = cell_scores

    print(f"Marker 基因命中率: {found_markers}/{total_markers}")
    return scores


asian_scores = calculate_immune_scores(asian_samples, asian_matrix)
white_scores = calculate_immune_scores(white_samples, white_matrix)


# ====== 第4步：比较两组的免疫细胞评分 ======
print("\n" + "=" * 60)
print("第4步：亚洲 vs 白人 免疫细胞评分比较")
print("=" * 60)


def compute_stats(values):
    """计算均值、标准差"""
    n = len(values)
    mean = sum(values) / n
    variance = sum((v - mean) ** 2 for v in values) / (n - 1) if n > 1 else 0
    std = math.sqrt(variance)
    return mean, std


def welch_ttest(values1, values2):
    """
    Welch t 检验——比较两组均值是否有显著差异。
    不需要假设两组方差相等。

    返回: (t统计量, p值近似值)
      这里用简化的 Welch 近似。
    """
    n1, n2 = len(values1), len(values2)
    m1, s1 = compute_stats(values1)
    m2, s2 = compute_stats(values2)

    # Welch t-statistic
    se = math.sqrt(s1**2 / n1 + s2**2 / n2)
    if se == 0:
        return 0, 1.0
    t_stat = (m1 - m2) / se

    # 自由度 (Welch-Satterthwaite)
    num = (s1**2 / n1 + s2**2 / n2) ** 2
    denom = (s1**2 / n1)**2 / (n1 - 1) + (s2**2 / n2)**2 / (n2 - 1)
    df = num / denom if denom > 0 else 1

    # p值近似（用正态逼近）
    # 真实项目中应使用 scipy.stats.t.sf
    import math as m
    # 简化: 用正态分布近似
    p_value = 2 * (1 - _norm_cdf(abs(t_stat)))

    return t_stat, p_value


def _norm_cdf(x):
    """标准正态分布 CDF 的近似（不需要 scipy）"""
    # Abramowitz and Stegun approximation
    if x < 0:
        return 1 - _norm_cdf(-x)
    # Constants
    b0 = 0.2316419
    b1 = 0.319381530
    b2 = -0.356563782
    b3 = 1.781477937
    b4 = -1.821255978
    b5 = 1.330274429
    t = 1 / (1 + b0 * x)
    phi = (1 / math.sqrt(2 * math.pi)) * math.exp(-x**2 / 2)
    return 1 - phi * (b1*t + b2*t**2 + b3*t**3 + b4*t**4 + b5*t**5)


# 计算每种免疫细胞的组间差异
print(f"\n{'细胞类型':<25} {'亚洲均值':>8} {'白人均值':>8} {'差异':>8} {'t值':>7} {'p值':>8} {'显著性':>6}")
print("-" * 80)

results = []
for cell_type in IMMUNE_MARKERS:
    a_vals = asian_scores[cell_type]
    w_vals = white_scores[cell_type]

    a_mean, a_std = compute_stats(a_vals)
    w_mean, w_std = compute_stats(w_vals)
    t_stat, p_val = welch_ttest(a_vals, w_vals)

    diff = a_mean - w_mean

    # 显著性标记
    sig = ""
    if p_val < 0.001:
        sig = "***"
    elif p_val < 0.01:
        sig = "**"
    elif p_val < 0.05:
        sig = "*"

    print(f"{cell_type:<25} {a_mean:8.3f} {w_mean:8.3f} {diff:+8.3f} {t_stat:+7.2f} {p_val:8.4f} {sig:>6}")

    results.append({
        "cell_type": cell_type,
        "asian_mean": a_mean,
        "white_mean": w_mean,
        "diff": diff,
        "t_stat": t_stat,
        "p_value": p_val,
        "significant": sig,
    })


# ====== 第5步：保存结果 ======
print("\n" + "=" * 60)
print("第5步：保存分析结果")
print("=" * 60)

result_path = os.path.join(OUTPUT_DIR, "immune_cell_scores.csv")
with open(result_path, "w", encoding="utf-8", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["cell_type", "asian_mean", "white_mean", "difference",
                     "t_statistic", "p_value", "significant"])
    for r in results:
        writer.writerow([
            r["cell_type"], r["asian_mean"], r["white_mean"],
            r["diff"], r["t_stat"], r["p_value"], r["significant"],
        ])

print(f"结果已保存: {result_path}")

# 保存每个样本的详细评分
sample_scores_path = os.path.join(OUTPUT_DIR, "sample_immune_scores.csv")
all_samples = list(asian_samples) + list(white_samples)
with open(sample_scores_path, "w", encoding="utf-8", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["sample_id", "group"] + list(IMMUNE_MARKERS.keys()))
    for i, sid in enumerate(asian_samples):
        row = [sid, "asian"]
        for ct in IMMUNE_MARKERS:
            row.append(asian_scores[ct][i])
        writer.writerow(row)
    for i, sid in enumerate(white_samples):
        row = [sid, "white"]
        for ct in IMMUNE_MARKERS:
            row.append(white_scores[ct][i])
        writer.writerow(row)

print(f"样本级评分已保存: {sample_scores_path}")


# ====== 第6步：文字解读 ======
print("\n" + "=" * 60)
print("初步解读")
print("=" * 60)

sig_results = [r for r in results if r["significant"]]
if sig_results:
    print(f"\n发现 {len(sig_results)} 个显著差异的免疫细胞类型：\n")
    for r in sorted(sig_results, key=lambda x: abs(x["diff"]), reverse=True):
        direction = "更高" if r["diff"] > 0 else "更低"
        print(f"  {r['cell_type']}: 亚洲组 {direction}")
        print(f"    亚洲: {r['asian_mean']:.3f}, 白人: {r['white_mean']:.3f}")
        print(f"    差异: {r['diff']:+.3f}, p={r['p_value']:.4f}")
        print()
else:
    print("\n未发现显著差异。可能需要：")
    print("  1. 更大的样本量")
    print("  2. 更精细的细胞亚群分类")
    print("  3. 考虑肿瘤分期等混杂因素")

print("=" * 60)
print("第5步完成！")
print("=" * 60)
print("""
下一步（第6步）：
  - 用 matplotlib 画对比箱线图
  - 画免疫评分热图
  - 让结果可视化
""")
