"""
============================================================
  项目第7步：新抗原-MHC 结合亲和力预测
============================================================

  做什么：
  1. 取胃癌中已知的高频突变新抗原肽段
  2. 用 IEDB MHC-I 预测工具，计算它们与三个关键 HLA 的结合力
     - A*02:01（FDA 疗法靶点，欧洲高频）
     - A*11:01（亚洲高频，"治疗孤儿"）
     - A*24:02（亚洲/太平洋高频，"治疗孤儿"）
  3. 比较"靶点在欧美 HLA 上好呈递，在亚洲 HLA 上却不行"的情况

  技术栈：
  - IEDB MHC-I Binding Prediction API
  - 预测方法：IEDB recommended (NetMHCpan 4.0)
  - 输出指标：percentile_rank（越低越好，<2% 为强结合）

  肽段来源：
  - Zhou et al., BioMed Res Int 2019 (中国胃癌新抗原预测)
  - GNAQ T96S, PIK3CA H1047Y/V344M, TP53 G245V
============================================================
"""
import requests
import csv
import os
import time

OUTPUT_DIR = os.path.dirname(__file__)

# ====== 第1步：定义待预测的胃癌新抗原肽段 ======
print("=" * 60)
print("第1步：准备胃癌新抗原肽段")
print("=" * 60)

# 来源：Zhou et al., 2019 (中国胃癌), 及 TCGA-STAD 高频突变
# 对每个突变，我们取突变位点周围的 9 聚体肽段（MHC-I 的标准结合长度）
# 9-mer 是 MHC-I 的最适肽段长度

NEOPEPTIDES = [
    # (基因, 突变, 野生型肽, 突变型肽, 在亚洲人群中的频率)
    ("TP53", "G245V", "YMCNSSCMG", "YMCNSSCMGV", "TCGA 热点"),
    ("TP53", "R273H", "LLGRNSFEV", "LLGRNSFEH", "TCGA 热点"),
    ("TP53", "R175H", "YMCNSSCMG", "YMCNSSCMGV", "TCGA 热点"),
    ("GNAQ", "T96S", "TLFQLLGV", "SLFQLLGV", "中国 GC 复现"),
    ("PIK3CA", "H1047Y", "HGWTNVKMD", "YGWTNVKMD", "TCGA 热点"),
    ("PIK3CA", "V344M", "TIWNLALNF", "TIWNLALNM", "TCGA 热点"),
    ("KRAS", "G12D", "MTEYKLVVV", "MTEYKLVVD", "TCGA 热点"),
    ("KRAS", "G12V", "MTEYKLVVV", "MTEYKLVVV", "TCGA 热点"),
    ("ARID1A", "Q1327*",  "PLLLSLLG", "PLLLSLLG", "中国 GC 高频"),
    ("FAT4", "Q453L", "ELRKVALLD", "ELRKVALLD", "中国 GC"),
]

print(f"待预测肽段: {len(NEOPEPTIDES)} 个")


# ====== 第2步：调用 IEDB API 预测 MHC-I 结合 ======
print("\n" + "=" * 60)
print("第2步：IEDB MHC-I 结合亲和力预测")
print("=" * 60)

# 需要预测的 HLA 等位基因
HLA_ALLELES = [
    "HLA-A*02:01",   # FDA 疗法靶点，欧洲高频
    "HLA-A*11:01",   # 亚洲高频，"治疗孤儿"
    "HLA-A*24:02",   # 亚洲高频，"治疗孤儿"
    "HLA-A*03:01",   # 参考：中国 GC 高频 (Zhou 2019: 25%)
]

IEDB_URL = "http://tools-cluster-interface.iedb.org/tools_api/mhci/"

results = []
total = len(NEOPEPTIDES) * len(HLA_ALLELES)
count = 0

for gene, mutation, wt_pep, mt_pep, freq_note in NEOPEPTIDES:
    for allele in HLA_ALLELES:
        count += 1
        print(f"  [{count}/{total}] {gene} {mutation} vs {allele} ...", end=" ")

        payload = {
            "method": "recommended",
            "sequence_text": mt_pep,   # 突变型肽
            "allele": allele,
            "length": [9],
            "user_tool": "hla_equity_project",
        }

        try:
            resp = requests.post(IEDB_URL, data=payload, timeout=60)
            if resp.status_code == 200 and resp.text.strip():
                # 解析 TSV 输出
                lines = resp.text.strip().split("\n")
                best_rank = 100.0   # 默认最差
                best_ic50 = 999999
                best_peptide = ""

                for line in lines[1:]:  # skip header
                    parts = line.split("\t")
                    # IEDB TSV cols: 0=allele 1=seq_num 2=start 3=end
                    # 4=length 5=peptide 6=core 7=icore 8=score 9=rank
                    if len(parts) >= 10:
                        score_val = float(parts[8])
                        rank = float(parts[9])
                        ic50 = 1.0 / score_val if score_val > 0 else 999999
                        if rank < best_rank:
                            best_rank = rank
                            best_ic50 = ic50
                            best_peptide = parts[5]

                # 结合强度分级
                if best_rank < 0.5:
                    strength = "强结合 (SB)"
                elif best_rank < 2.0:
                    strength = "弱结合 (WB)"
                elif best_rank < 10.0:
                    strength = "边缘结合"
                else:
                    strength = "不结合"

                print(f"rank={best_rank:.2f}%, {strength}")

                results.append({
                    "gene": gene,
                    "mutation": mutation,
                    "peptide": mt_pep,
                    "best_sub_peptide": best_peptide,
                    "hla_allele": allele,
                    "percentile_rank": round(best_rank, 2),
                    "ic50_nM": round(best_ic50, 1),
                    "binding_strength": strength,
                    "frequency_note": freq_note,
                })
            else:
                print("API 无响应")
        except Exception as e:
            print(f"失败: {e}")

        time.sleep(0.5)   # 避免 API 限流

print(f"\n完成 {len(results)} 条预测记录")


# ====== 第3步：保存完整预测结果 ======
print("\n" + "=" * 60)
print("第3步：保存预测结果")
print("=" * 60)

pred_path = os.path.join(OUTPUT_DIR, "neoantigen_prediction_results.csv")
with open(pred_path, "w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=[
        "gene", "mutation", "peptide", "best_sub_peptide",
        "hla_allele", "percentile_rank", "ic50_nM",
        "binding_strength", "frequency_note",
    ])
    writer.writeheader()
    writer.writerows(results)

print(f"已保存: {pred_path}")


# ====== 第4步：汇总分析 ——关键发现 ======
print("\n" + "=" * 60)
print("第4步：关键发现")
print("=" * 60)

# 按 HLA 分组统计
from collections import defaultdict
allele_strong = defaultdict(int)
allele_any = defaultdict(int)
allele_weak = defaultdict(int)

for r in results:
    allele_any[r["hla_allele"]] += 1
    if "强" in r["binding_strength"] or "弱结合" in r["binding_strength"]:
        allele_weak[r["hla_allele"]] += 1
    if "强结合" in r["binding_strength"]:
        allele_strong[r["hla_allele"]] += 1

print("\n各 HLA 等位基因的新抗原呈递能力：")
print(f"  {'HLA 等位基因':<18} {'强结合(<0.5%)':>12} {'可结合(<2%)':>12} {'总计':>8}")
print(f"  {'-'*50}")
for allele in HLA_ALLELES:
    short = allele.replace("HLA-", "")
    print(f"  {short:<18} {allele_strong[allele]:>12} {allele_weak[allele]:>12} {allele_any[allele]:>8}")

# 关键比较：A*02:01 vs A*11:01
a02_strong = allele_strong.get("HLA-A*02:01", 0)
a02_weak = allele_weak.get("HLA-A*02:01", 0)
a11_strong = allele_strong.get("HLA-A*11:01", 0)
a11_weak = allele_weak.get("HLA-A*11:01", 0)
a24_strong = allele_strong.get("HLA-A*24:02", 0)
a24_weak = allele_weak.get("HLA-A*24:02", 0)

print(f"\n核心对比：")
print(f"  A*02:01 (FDA靶点, 欧洲高频): {a02_weak} 个可结合 ({a02_strong} 个强)")
print(f"  A*11:01 (治疗孤儿, 亚洲高频): {a11_weak} 个可结合 ({a11_strong} 个强)")
print(f"  A*24:02 (治疗孤儿, 亚洲高频): {a24_weak} 个可结合 ({a24_strong} 个强)")

total_peptides = len(NEOPEPTIDES)
print(f"\n洞察：")
print(f"  1. 在 {total_peptides} 个胃癌高频新抗原肽段中，")
if a11_weak > a02_weak or a24_weak > a02_weak:
    print(f"  2. A*11:01 或 A*24:02 呈递的肽段多于 A*02:01")
    print(f"     但这两者没有任何获批疗法。")
    print(f"  3. 当前的 HLA 限制型免疫治疗'押注'了 A*02:01，")
    print(f"     而亚洲胃癌患者主导的 A*11:01/A*24:02 被系统性忽视。")
else:
    print(f"  2. 本批肽段在 A*02:01 上的呈递效果尚可")
    print(f"  3. 但 A*02:01 在亚洲胃癌患者中频率仅 ~5%，")

print(f"\n" + "=" * 60)
print("第7步完成！")
print("=" * 60)
