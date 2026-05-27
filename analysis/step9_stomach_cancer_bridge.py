"""
============================================================
  项目第9步：胃癌特异性数据桥接
============================================================

  做什么：
  用已发表的 GLOBOCAN 2022 和 GBD 2021 文献数据，
  将第8步的官方"肿瘤"大类锚定到"胃癌"。

  核心数据点：
  - 东亚贡献全球 53.8%~60.8% 的胃癌新发病例
  - 胃癌是全球第5高发癌症、第5大癌症死因

  引用来源：
  - GLOBOCAN 2022: Bray et al., CA Cancer J Clin 2024
  - Tan et al., Cancer Biol Med 2024
  - GBD 2021: JMIR Cancer 2025

  输出：胃癌全球分布的参考数据表 + 桥接论述
"""
import csv
import os

OUTPUT_DIR = os.path.dirname(__file__)

print("=" * 60)
print("第9步：胃癌特异性数据桥接")
print("=" * 60)

# ====== 第1步：GLOBOCAN 2022 胃癌全球关键数据 ======
print("\n-- GLOBOCAN 2022 全球胃癌关键统计 --")

# 来源: Bray et al. 2024, CA: A Cancer Journal for Clinicians
# https://pubmed.ncbi.nlm.nih.gov/38572751/
print("""
  GLOBOCAN 2022 (Bray et al., CA Cancer J Clin 2024):
    全球新发病例: ~968,350 (第5位, 占全部癌症的 4.9%)
    全球死亡: ~659,853 (第5位, 占全部癌症死亡的 6.8%)
    ASIR: 9.2/10万, ASMR: 6.1/10万
""")

# ====== 第2步：区域分布数据 ======
print("-- 胃癌区域分布 --")

# 来源: Tan et al., Cancer Biol Med 2024
# https://www.cancerbiomed.org/content/21/8/667
global_burden = [
    ("East Asia", 53.8, 48.2, "东亚承担全球半数以上的胃癌负担"),
    ("Europe", 13.6, 14.0, "欧洲是第二大负担区域"),
    ("South America", 5.3, 6.0, ""),
    ("South-Central Asia", 5.1, 6.4, ""),
    ("South-East Asia", 4.7, 5.6, ""),
    ("North America", 2.4, 1.8, "北美胃癌负担最低"),
    ("Africa", 4.2, 5.7, "非洲诊断不足，实际负担可能被低估"),
    ("Other", 10.9, 12.3, ""),
]

for region, inc_pct, mort_pct, note in global_burden:
    note_str = f"  ← {note}" if note else ""
    print(f"  {region:25s}: 新发病例 {inc_pct:5.1f}%, 死亡 {mort_pct:5.1f}%{note_str}")

# ====== 第3步：GBD 2021 验证数据 ======
print("\n-- GBD 2021 验证 (2021年数据) --")

# 来源:
# - Jiang et al., Clin Epidemiol 2024 (PMC11381218)
# - JMIR Cancer 2025 (PMC12334143)
# - BMC Public Health 2025 (PMC12326841)
print("""
  GBD 2021 (多个独立来源):
    东亚新发病例: 748,235 例 (占全球 ~60.8%)
    东亚死亡: 527,054 例 (占全球 ~55.2%)
    东亚 DALYs: ~12.1M (占全球 ~53.2%)

  GLOBOCAN 和 GBD 数据一致：
    东亚承担了全球胃癌负担的 50-60%
""")

# ====== 第4步：核心桥接论述 ======
print("\n" + "=" * 60)
print("第4步：桥接论述")
print("=" * 60)

print("""
  链式论证：

  (1) 第8步官方数据:
      东亚承担全球最高的肿瘤死亡绝对负担 (3.79M, 2023)
      肿瘤在东亚地区的死因排名持续上升

  (2) GLOBOCAN 2022 (Bray et al. 2024):
      胃癌是全球第5高发癌症，968,350 新发病例/年
      胃癌是全球第5大癌症死因，659,853 死亡/年

  (3) Tan et al. 2024, Cancer Biol Med:
      东亚贡献了全球 53.8% 的新发胃癌病例

  (4) GBD 2021 多中心研究一致验证:
      东亚胃癌 DALYs 超 1,200 万年
      高盐饮食、吸烟、幽门螺杆菌是东亚主要风险因素

  结论：
      东亚的肿瘤负担中，胃癌占比远高于全球平均水平。
      胃癌是东亚最具代表性的癌症类型之一，
      因此，研究胃癌患者的免疫治疗公平性问题，
      对东亚人群具有最高的公共卫生优先级。

  → 进入第10步：用 TCGA-STAD 突变数据量化
     "东亚胃癌患者的 HLA 等位基因与药物靶点匹配度"。
""")

# ====== 第5步：保存参考数据 ======
print("=" * 60)
print("第5步：保存胃癌全球参考数据")
print("=" * 60)

ref_path = os.path.join(OUTPUT_DIR, "stomach_cancer_global_burden.csv")
with open(ref_path, "w", encoding="utf-8", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["metric", "value", "source", "notes"])
    writer.writerow(["global_new_cases_2022", 968350, "GLOBOCAN 2022", "Bray et al., CA Cancer J Clin 2024"])
    writer.writerow(["global_deaths_2022", 659853, "GLOBOCAN 2022", ""])
    writer.writerow(["global_rank_incidence", 5, "GLOBOCAN 2022", "of 36 cancer types"])
    writer.writerow(["global_rank_mortality", 5, "GLOBOCAN 2022", "of 36 cancer types"])
    writer.writerow(["east_asia_share_incidence_pct", 53.8, "Tan et al. 2024", "Cancer Biol Med 2024;21(8):667-678"])
    writer.writerow(["east_asia_share_mortality_pct", 48.2, "Tan et al. 2024", ""])
    writer.writerow(["east_asia_cases_gbd2021", 748235, "GBD 2021", "Jiang et al., Clin Epidemiol 2024"])
    writer.writerow(["east_asia_deaths_gbd2021", 527054, "GBD 2021", ""])
    writer.writerow(["east_asia_dalys_millions", 12.1, "GBD 2021", "JMIR Cancer 2025"])
    writer.writerow(["china_asir_per_100k", 29.1, "GBD 2021", "age-standardized incidence rate"])
    writer.writerow(["mongolia_asir_per_100k", 36.8, "GBD 2021", "highest globally"])

print(f"已保存: {ref_path}")

print("\n" + "=" * 60)
print("第9步完成！")
print("=" * 60)
