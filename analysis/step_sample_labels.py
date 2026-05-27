import csv, os, json, requests
from collections import Counter

os.makedirs("/home/mw/project/tables", exist_ok=True)
OUT_PATH = "/home/mw/project/tables/sample_labels.csv"
EXP_PATH = "/home/mw/input/expression76957695/STAD_expression.tsv"

# ====== 第1步：从 GDC POST API 获取种族数据（跟本地 step4b 完全一样）======
print("第1步：GDC API 获取种族数据...")

query = {
    "filters": {
        "op": "in",
        "content": {"field": "cases.project.project_id", "value": ["TCGA-STAD"]}
    },
    "fields": "case_id,submitter_id,demographic.race,demographic.ethnicity,demographic.gender",
    "format": "JSON",
    "size": 1000
}

resp = requests.post("https://api.gdc.cancer.gov/cases", json=query, timeout=60)
cases = resp.json().get("data", {}).get("hits", [])
print(f"GDC 返回 {len(cases)} 个病例")

gdc_race = {}
for case in cases:
    submitter = case.get("submitter_id", "")
    demo = case.get("demographic") or {}
    race = (demo.get("race") or "").strip().lower()
    if submitter and race and race != "not reported":
        gdc_race[submitter] = race.upper()

print(f"有种族标注: {len(gdc_race)}")
race_dist = Counter(gdc_race.values())
for r, c in race_dist.most_common():
    print(f"  {r}: {c}")

# ====== 第2步：从表达矩阵提取样本 ID（跟本地 step4b 完全一样）======
print("\n第2步：匹配样本到种族...")

def sample_to_patient(sid):
    parts = sid.split("-")
    if len(parts) >= 3 and parts[0] == "TCGA":
        return "-".join(parts[:3])
    return None

with open(EXP_PATH, "r") as f:
    header = f.readline().strip().split("\t")
sample_ids = header[1:]
print(f"表达矩阵样本数: {len(sample_ids)}")

asian_samples = []
white_samples = []
other = []
unmatched = []

for sid in sorted(sample_ids):
    pid = sample_to_patient(sid)
    if pid is None:
        unmatched.append((sid, "bad_format"))
        continue
    race = gdc_race.get(pid, "")
    if race == "ASIAN":
        asian_samples.append(sid)
    elif race == "WHITE":
        white_samples.append(sid)
    elif race:
        other.append((sid, pid, race))
    else:
        unmatched.append((sid, "no_race"))

# ====== 第3步：只保留原发肿瘤 (01) ======
print("\n第3步：筛选原发肿瘤样本 (01型)...")

def is_tumor(sid):
    parts = sid.split("-")
    return len(parts) >= 4 and parts[3][:2] == "01"

tumor_asian = [s for s in asian_samples if is_tumor(s)]
tumor_white = [s for s in white_samples if is_tumor(s)]

print(f"亚洲: {len(tumor_asian)} 肿瘤 (原始 {len(asian_samples)})")
print(f"白人: {len(tumor_white)} 肿瘤 (原始 {len(white_samples)})")

# ====== 第4步：保存标签 ======
print("\n第4步：保存...")

with open(OUT_PATH, "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f)
    w.writerow(["sample_id", "patient_id", "race", "group"])
    for sid in tumor_asian:
        w.writerow([sid, sample_to_patient(sid), "ASIAN", "asian"])
    for sid in tumor_white:
        w.writerow([sid, sample_to_patient(sid), "WHITE", "white"])

print(f"亚裔: {len(tumor_asian)} | 白人: {len(tumor_white)} | 总计: {len(tumor_asian)+len(tumor_white)}")
print(f"已保存: {OUT_PATH}")
