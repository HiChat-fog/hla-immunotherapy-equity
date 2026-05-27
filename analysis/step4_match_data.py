"""
============================================================
  项目第4步：样本-患者匹配与种族分组
============================================================

  做什么：
  1. 加载表达矩阵（样本 × 基因）
  2. 加载临床数据（患者 × 属性）
  3. 从样本 ID 提取患者 ID（TCGA 条形码标准）
  4. 匹配种族信息
  5. 生成两个"干净"的表达子集：亚洲组 vs 白人组

  正确性保证：
  - TCGA 条形码规则是国际标准，不靠猜测
  - 每步匹配都输出统计，暴露任何异常
  - 无法匹配的样本标记排除，不影响分析
============================================================
"""
import csv
import os

OUTPUT_DIR = os.path.dirname(__file__)
EXP_FILE = os.path.join(OUTPUT_DIR, "STAD_expression.tsv")
CLIN_FILE = os.path.join(OUTPUT_DIR, "clinical_data.csv")


# ====== 第1步：加载临床数据，建立 患者ID → 种族 的字典 ======
print("=" * 60)
print("第1步：加载临床数据")
print("=" * 60)

patient_race = {}   # {患者ID: 种族}
patient_info = {}   # {患者ID: {其他临床属性}}

with open(CLIN_FILE, "r", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)
    for row in reader:
        pid = row.get("patient_id", "").strip()
        race = row.get("RACE", "").strip()

        if not pid:
            continue

        patient_race[pid] = race
        patient_info[pid] = row

print(f"加载了 {len(patient_race)} 个患者")

# 统计种族
from collections import Counter
race_dist = Counter(patient_race.values())
print("临床数据中的种族分布：")
for race, count in race_dist.most_common():
    print(f"  {race}: {count} 人")

asian_patients = {pid for pid, r in patient_race.items() if r == "ASIAN"}
white_patients = {pid for pid, r in patient_race.items() if r == "WHITE"}
print(f"\n可用于分析的：")
print(f"  亚洲人 (ASIAN): {len(asian_patients)} 人")
print(f"  白人 (WHITE):   {len(white_patients)} 人")
print(f"  合计:            {len(asian_patients) + len(white_patients)} 人")


# ====== 第2步：解析 TCGA 样本 ID，匹配到患者 ======
print("\n" + "=" * 60)
print("第2步：解析表达矩阵中的样本 ID")
print("=" * 60)

def sample_to_patient(sample_id):
    """
    从 TCGA 样本条形码中提取患者 ID。

    TCGA 条形码格式：
    TCGA-XX-YYYY-ZZA-BBC-DDDD-EE
    患者 ID = 前三个 '-' 之前的部分

    验证规则：
    - 至少包含 3 个 '-'
    - 提取的部分以 'TCGA-' 开头
    """
    parts = sample_id.split("-")
    if len(parts) >= 3:
        patient_id = "-".join(parts[:3])   # TCGA-XX-YYYY
        if patient_id.startswith("TCGA-"):
            return patient_id
    return None


# 读表达矩阵的第一行（表头 = 样本ID列表）
with open(EXP_FILE, "r") as f:
    header = f.readline().strip().split("\t")

# 第一列是 'sample'（基因名列），其余是样本ID
sample_ids = header[1:]   # 去掉第一列的 'sample'
print(f"表达矩阵中的样本总数: {len(sample_ids)}")

# 验证每一个样本 ID 的格式
invalid_samples = []
sample_to_patient_map = {}
for sid in sample_ids:
    pid = sample_to_patient(sid)
    if pid is None:
        invalid_samples.append(sid)
    else:
        sample_to_patient_map[sid] = pid

print(f"成功解析的患者 ID: {len(sample_to_patient_map)} 个")
print(f"格式异常的样本:     {len(invalid_samples)} 个")
if invalid_samples:
    print(f"  异常样本示例: {invalid_samples[:3]}")

# 找出哪些样本匹配到了亚洲/白人患者
asian_samples = []
white_samples = []
matched_other = []
unmatched_samples = []
unmatched_count = 0

for sid, pid in sample_to_patient_map.items():
    if pid in asian_patients:
        asian_samples.append(sid)
    elif pid in white_patients:
        white_samples.append(sid)
    elif pid in patient_race:
        # 有种族标注但不是亚洲也不是白人
        matched_other.append((sid, pid, patient_race[pid]))
    else:
        unmatched_count += 1

print(f"\n匹配结果：")
print(f"  亚洲人样本: {len(asian_samples)} 个（来自 {len(set(sample_to_patient_map[s] for s in asian_samples))} 个患者）")
print(f"  白人样本:   {len(white_samples)} 个（来自 {len(set(sample_to_patient_map[s] for s in white_samples))} 个患者）")
print(f"  其他种族:   {len(matched_other)} 个")
print(f"  无法匹配:   {unmatched_count} 个")

# 如果一个患者有多个样本（如肿瘤+正常），说明一下
asian_patient_sample_count = Counter(sample_to_patient_map[s] for s in asian_samples)
duplicated_asian = [(pid, cnt) for pid, cnt in asian_patient_sample_count.items() if cnt > 1]
if duplicated_asian:
    print(f"\n  注意：{len(duplicated_asian)} 个亚洲患者有多个样本")
    for pid, cnt in duplicated_asian[:5]:
        print(f"    {pid}: {cnt} 个样本")

white_patient_sample_count = Counter(sample_to_patient_map[s] for s in white_samples)
duplicated_white = [(pid, cnt) for pid, cnt in white_patient_sample_count.items() if cnt > 1]
if duplicated_white:
    print(f"  {len(duplicated_white)} 个白人患者有多个样本")


# ====== 第3步：提取匹配到的表达数据 ======
print("\n" + "=" * 60)
print("第3步：提取亚洲组和白人组的表达数据")
print("=" * 60)

# 我们需要每行基因的表达值 → 转置为 样本×基因 矩阵
# 为了方便后续分析，直接生成两个文件：
#   1. asian_expression.tsv  → 只含亚洲样本
#   2. white_expression.tsv  → 只含白人样本

asian_set = set(asian_samples)
white_set = set(white_samples)

with open(EXP_FILE, "r") as f_in:
    all_lines = f_in.readlines()

header_parts = all_lines[0].strip().split("\t")
# 找出哪些列是亚洲样本、哪些是白人样本
asian_col_indices = []
white_col_indices = []
gene_col_index = 0  # 基因名列

for i, sid in enumerate(header_parts):
    if sid in asian_set:
        asian_col_indices.append(i)
    elif sid in white_set:
        white_col_indices.append(i)

print(f"亚洲样本列索引: {len(asian_col_indices)} 个")
print(f"白人样本列索引: {len(white_col_indices)} 个")

def write_subset(lines, col_indices, out_path, label):
    """写出只包含指定列的表达子集"""
    with open(out_path, "w", encoding="utf-8", newline="") as f_out:
        writer = csv.writer(f_out, delimiter="\t")
        for line in lines:
            parts = line.strip().split("\t")
            # 只保留基因列 + 选中的样本列
            row = [parts[0]] + [parts[i] for i in col_indices]
            writer.writerow(row)

    # 验证
    with open(out_path, "r") as f:
        n_genes = sum(1 for _ in f) - 1
    print(f"  {label}: {n_genes} 个基因 × {len(col_indices)} 个样本 → {out_path}")

write_subset(all_lines, asian_col_indices,
             os.path.join(OUTPUT_DIR, "asian_expression.tsv"), "亚洲组")
write_subset(all_lines, white_col_indices,
             os.path.join(OUTPUT_DIR, "white_expression.tsv"), "白人组")


# ====== 第4步：生成分析就绪的标签文件 ======
print("\n" + "=" * 60)
print("第4步：生成样本标签文件（后续分析用）")
print("=" * 60)

labels_path = os.path.join(OUTPUT_DIR, "sample_labels.csv")
with open(labels_path, "w", encoding="utf-8", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["sample_id", "patient_id", "race", "group"])
    for sid in asian_samples:
        pid = sample_to_patient_map[sid]
        writer.writerow([sid, pid, "ASIAN", "asian"])
    for sid in white_samples:
        pid = sample_to_patient_map[sid]
        writer.writerow([sid, pid, "WHITE", "white"])

print(f"标签文件: {labels_path}")
print(f"  亚洲样本: {len(asian_samples)}")
print(f"  白人样本: {len(white_samples)}")
print(f"  总计:     {len(asian_samples) + len(white_samples)}")


# ====== 第5步：最终正确性检查 ======
print("\n" + "=" * 60)
print("第5步：最终质量检查")
print("=" * 60)

# 检查1：每组有没有重复样本
if len(asian_samples) != len(set(asian_samples)):
    print("[警告] 亚洲组有重复样本ID！")
if len(white_samples) != len(set(white_samples)):
    print("[警告] 白人组有重复样本ID！")

# 检查2：两组有没有重叠
overlap = set(asian_samples) & set(white_samples)
if overlap:
    print(f"[警告] 亚洲组和白人组有 {len(overlap)} 个重叠样本！")
else:
    print("[通过] 亚洲组和白人组无重叠")

# 检查3：样本数是否合理
print(f"[信息] 亚洲组 {len(asian_samples)} 样本 vs 白人组 {len(white_samples)} 样本")
if len(asian_samples) < 10:
    print("[警告] 亚洲样本少于10个，统计效力不足")
else:
    print("[通过] 样本数充足，可进行统计分析")

# 检查4：每个样本都能回溯到患者和种族
unmatched_check = []
for sid in asian_samples + white_samples:
    pid = sample_to_patient_map.get(sid)
    if pid is None:
        unmatched_check.append(sid)
    elif pid not in patient_race:
        unmatched_check.append(sid)
if unmatched_check:
    print(f"[警告] {len(unmatched_check)} 个样本无法回溯")
else:
    print("[通过] 所有样本均可回溯到患者和种族")

print("\n" + "=" * 60)
print("第4步完成！")
print("=" * 60)
print(f"""
生成的文件：
  asian_expression.tsv    — 亚洲组表达矩阵
  white_expression.tsv    — 白人组表达矩阵
  sample_labels.csv       — 样本-种族对照表
  clinical_data.csv       — 完整临床数据（已有）

下一步：对每个基因，比较亚洲组 vs 白人组的表达差异。
""")
