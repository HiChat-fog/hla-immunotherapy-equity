"""
============================================================
  项目第4步（修正版）：用 GDC 官方数据重新匹配
============================================================
  GDC = NCI Genomic Data Commons，TCGA 数据的官方仓库
  临床数据完整性远超 cBioPortal
============================================================
"""
import csv
import requests
import os
from collections import Counter

OUTPUT_DIR = os.path.dirname(__file__)
EXP_FILE = os.path.join(OUTPUT_DIR, "STAD_expression.tsv")


# ====== 第1步：从 GDC 获取所有 STAD 病例的种族数据 ======
print("=" * 60)
print("第1步：从 GDC API 获取完整临床数据")
print("=" * 60)

query = {
    "filters": {
        "op": "in",
        "content": {
            "field": "cases.project.project_id",
            "value": ["TCGA-STAD"],
        },
    },
    "fields": (
        "case_id,"
        "submitter_id,"
        "demographic.race,"
        "demographic.ethnicity,"
        "demographic.gender,"
    ),
    "format": "JSON",
    "size": 1000,
}

resp = requests.post(
    "https://api.gdc.cancer.gov/cases",
    json=query,
    timeout=60,
)
cases = resp.json().get("data", {}).get("hits", [])
print(f"GDC 返回 {len(cases)} 个病例")

# 建立 submitter_id → race 的映射
# GDC 的 submitter_id 就是 TCGA 患者 ID（如 TCGA-FP-A4BE）
gdc_race = {}
for case in cases:
    submitter = case.get("submitter_id", "")
    demo = case.get("demographic") or {}
    race = demo.get("race", "").strip().lower()
    if submitter and race and race != "not reported":
        # 统一格式
        race_clean = race.upper()
        gdc_race[submitter] = race_clean

print(f"有种族标注的病例: {len(gdc_race)}")
race_dist = Counter(gdc_race.values())
print("种族分布:")
for race, cnt in race_dist.most_common():
    print(f"  {race}: {cnt}")


# ====== 第2步：从表达矩阵解析样本ID，匹配到GDC种族 ======
print("\n" + "=" * 60)
print("第2步：重新匹配表达样本到 GDC 种族")
print("=" * 60)


def sample_to_patient(sid):
    parts = sid.split("-")
    if len(parts) >= 3 and parts[0] == "TCGA":
        return "-".join(parts[:3])
    return None


with open(EXP_FILE, "r") as f:
    header = f.readline().strip().split("\t")

sample_ids = header[1:]
print(f"表达矩阵样本数: {len(sample_ids)}")

# 匹配
asian_samples = []
white_samples = []
unmatched = []
matched_other = []
patient_sample_count = {}  # 统计每个患者的样本数

for sid in sample_ids:
    pid = sample_to_patient(sid)
    if pid is None:
        unmatched.append((sid, "格式异常"))
        continue

    patient_sample_count[pid] = patient_sample_count.get(pid, 0) + 1
    race = gdc_race.get(pid, "")

    if race == "ASIAN":
        asian_samples.append(sid)
    elif race == "WHITE":
        white_samples.append(sid)
    elif race:
        matched_other.append((sid, pid, race))
    else:
        unmatched.append((sid, "GDC无种族数据"))

print(f"\n匹配结果:")
asian_patients = set(sample_to_patient(s) for s in asian_samples)
white_patients = set(sample_to_patient(s) for s in white_samples)
print(f"  亚洲样本: {len(asian_samples)} 个（来自 {len(asian_patients)} 个患者）")
print(f"  白人样本: {len(white_samples)} 个（来自 {len(white_patients)} 个患者）")
print(f"  其他种族: {len(matched_other)} 个")
print(f"  无法匹配: {len(unmatched)} 个")

if unmatched:
    # 分类未匹配的原因
    no_race = sum(1 for _, reason in unmatched if "无种族" in reason)
    bad_fmt = sum(1 for _, reason in unmatched if "格式" in reason)
    print(f"    其中: GDC无种族={no_race}, 格式异常={bad_fmt}")

# 对比之前 cBioPortal 的结果
print(f"\n  [对比] cBioPortal: 37亚洲/135白人 → GDC: {len(asian_samples)}亚洲/{len(white_samples)}白人")
print(f"  [提升] 亚洲样本增加了 {len(asian_samples) - 37} 个！")


# ====== 第3步：筛选 -- 只保留肿瘤样本 ======
print("\n" + "=" * 60)
print("第3步：筛选样本类型")
print("=" * 60)

# TCGA 样本类型编码: 01-09=肿瘤, 10-19=正常, 20-29=对照
# 我们只分析原发肿瘤 (01)
tumor_asian = []
tumor_white = []
normal_asian = []
normal_white = []

for sid in asian_samples:
    parts = sid.split("-")
    if len(parts) >= 4:
        sample_type = parts[3][:2]
        if sample_type == "01":
            tumor_asian.append(sid)
        else:
            normal_asian.append((sid, sample_type))

for sid in white_samples:
    parts = sid.split("-")
    if len(parts) >= 4:
        sample_type = parts[3][:2]
        if sample_type == "01":
            tumor_white.append(sid)
        else:
            normal_white.append((sid, sample_type))

print(f"亚洲组: {len(tumor_asian)} 肿瘤 + {len(normal_asian)} 正常/其他")
print(f"白人组: {len(tumor_white)} 肿瘤 + {len(normal_white)} 正常/其他")

if normal_asian:
    print(f"  亚洲非肿瘤样本: {normal_asian}")
if normal_white:
    print(f"  白人非肿瘤样本: {normal_white[:5]}...")


# ====== 第4步：生成最终表达子集 ======
print("\n" + "=" * 60)
print("第4步：生成最终数据文件")
print("=" * 60)

asian_set = set(tumor_asian)
white_set = set(tumor_white)

with open(EXP_FILE, "r") as f_in:
    all_lines = f_in.readlines()

header_parts = all_lines[0].strip().split("\t")
asian_cols = [i for i, sid in enumerate(header_parts) if sid in asian_set]
white_cols = [i for i, sid in enumerate(header_parts) if sid in white_set]

print(f"亚洲肿瘤样本: {len(asian_cols)} 列")
print(f"白人肿瘤样本: {len(white_cols)} 列")


def write_subset(lines, col_indices, out_path, label):
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        for line in lines:
            parts = line.strip().split("\t")
            row = [parts[0]] + [parts[i] for i in col_indices]
            writer.writerow(row)

    with open(out_path, "r") as f:
        n_genes = sum(1 for _ in f) - 1
    print(f"  {label}: {n_genes} 基因 x {len(col_indices)} 样本 → {out_path}")

write_subset(all_lines, asian_cols,
             os.path.join(OUTPUT_DIR, "asian_expression.tsv"), "亚洲肿瘤组")
write_subset(all_lines, white_cols,
             os.path.join(OUTPUT_DIR, "white_expression.tsv"), "白人肿瘤组")

# 第5步：更新标签文件
labels_path = os.path.join(OUTPUT_DIR, "sample_labels.csv")
with open(labels_path, "w", encoding="utf-8", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["sample_id", "patient_id", "race", "group"])
    for sid in tumor_asian:
        pid = sample_to_patient(sid)
        writer.writerow([sid, pid, "ASIAN", "asian"])
    for sid in tumor_white:
        pid = sample_to_patient(sid)
        writer.writerow([sid, pid, "WHITE", "white"])

print(f"\n标签文件: {labels_path}")
print(f"  亚洲肿瘤样本: {len(tumor_asian)}")
print(f"  白人肿瘤样本: {len(tumor_white)}")


# ====== 最终检查 ======
print("\n" + "=" * 60)
print("最终质量检查")
print("=" * 60)

print(f"[数据] 亚洲组 {len(tumor_asian)} vs 白人组 {len(tumor_white)}")
print(f"[来源] GDC 官方 API（NCI Genomic Data Commons）")
print(f"[日期] TCGA-STAD 项目，共 {len(cases)} 例")
print(f"[覆盖] 种族标注率 {len(gdc_race)}/{len(cases)} = {len(gdc_race)/len(cases)*100:.0f}%")

if len(set(tumor_asian) & set(tumor_white)) == 0:
    print("[通过] 无样本重叠")
else:
    print("[警告] 存在样本重叠！")

print(f"[通过] 所有样本均有明确种族标注")
print(f"[通过] 仅保留原发肿瘤样本（01型）")
