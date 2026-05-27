"""
============================================================
  项目第10步：TCGA-STAD 突变数据下载 (MAF)
============================================================

  做什么：
  从 GDC API 下载 TCGA-STAD 的突变注释文件 (MAF)。
  MAF 包含每个样本的所有体细胞突变 ——
  氨基酸变化、基因组位置、突变类型等。
  这些信息是第11步提取突变肽段的前提。

  文件说明:
  - MAF格式基于 VCF 但更易读
  - 每行 = 一个突变
  - 关键列: Hugo_Symbol, HGVSp_Short, Tumor_Sample_Barcode
"""
import requests
import json
import os
import gzip

OUTPUT_DIR = os.path.dirname(__file__)

# ====== 第1步：通过 GDC API 找到 STAD 的 MAF 文件 ======
print("=" * 60)
print("第1步：在 GDC 中搜索 TCGA-STAD MAF 文件")
print("=" * 60)

# GDC API: 查询 open-access MAF 文件
# TCGA-STAD 的突变数据是开放获取的
files_endpoint = "https://api.gdc.cancer.gov/files"

params = {
    "filters": json.dumps({
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {
                    "field": "cases.project.project_id",
                    "value": ["TCGA-STAD"]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.data_type",
                    "value": ["Masked Somatic Mutation"]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.data_format",
                    "value": ["MAF"]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.access",
                    "value": ["open"]
                }
            }
        ]
    }),
    "size": "10",
    "fields": "file_id,file_name,file_size,data_type,data_format,cases.project.project_id"
}

resp = requests.get(files_endpoint, params=params)
data = resp.json()

print(f"API 查询: {resp.status_code}")
hits = data.get("data", {}).get("hits", [])
print(f"找到 {len(hits)} 个 MAF 文件")

for hit in hits:
    fsize_mb = hit.get("file_size", 0) / (1024 * 1024)
    print(f"  - {hit.get('file_name', 'N/A')}: {fsize_mb:.1f} MB (ID: {hit.get('file_id', 'N/A')})")

if not hits:
    print("\n未找到 open-access MAF。尝试放宽条件...")
    # Try without access filter
    params2 = {
        "filters": json.dumps({
            "op": "and",
            "content": [
                {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-STAD"]}},
                {"op": "in", "content": {"field": "files.data_type", "value": ["Masked Somatic Mutation"]}},
                {"op": "in", "content": {"field": "files.data_format", "value": ["MAF"]}}
            ]
        }),
        "size": "10",
        "fields": "file_id,file_name,file_size,access"
    }
    resp2 = requests.get(files_endpoint, params=params2)
    data2 = resp2.json()
    hits2 = data2.get("data", {}).get("hits", [])
    print(f"所有 MAF (含 controlled): {len(hits2)} 个")
    for hit in hits2:
        print(f"  - {hit.get('file_name')}: access={hit.get('access')}")


# ====== 第2步：下载 MAF 文件 ======
if hits:
    print("\n" + "=" * 60)
    print("第2步：下载 MAF 文件")
    print("=" * 60)

    file_id = hits[0]["file_id"]
    file_name = hits[0]["file_name"]
    file_size = hits[0]["file_size"]
    fsize_mb = file_size / (1024 * 1024)

    print(f"下载: {file_name} ({fsize_mb:.1f} MB)")

    # GDC data download endpoint
    download_url = f"https://api.gdc.cancer.gov/data/{file_id}"

    maf_path = os.path.join(OUTPUT_DIR, file_name)
    # 如果文件已存在，跳过
    if os.path.exists(maf_path) and os.path.getsize(maf_path) > 1000:
        print(f"文件已存在: {maf_path} ({os.path.getsize(maf_path)/1024/1024:.1f} MB)")
    else:
        dl_resp = requests.get(download_url, stream=True)
        print(f"  下载状态: {dl_resp.status_code}")

        if dl_resp.status_code == 200:
            total = 0
            with open(maf_path, "wb") as f:
                for chunk in dl_resp.iter_content(chunk_size=8192 * 1024):
                    f.write(chunk)
                    total += len(chunk)
                    pct = min(100, total / file_size * 100)
                    print(f"  [{total/1024/1024:.1f}/{fsize_mb:.1f} MB] {pct:.0f}%", end="\r")
            print(f"\n已保存: {maf_path}")
        else:
            print(f"下载失败! 状态码: {dl_resp.status_code}")
            print(f"可能需要 GDC 认证 token。")

# ====== 第3步：快速检查 MAF 内容 ======
maf_file = os.path.join(OUTPUT_DIR, "TCGA-STAD.mutect2_snv.tsv.gz")
# Find the actual downloaded file
actual_maf = None
for fname in os.listdir(OUTPUT_DIR):
    if fname.endswith(".maf") or fname.endswith(".maf.gz") or "mutect" in fname.lower():
        actual_maf = os.path.join(OUTPUT_DIR, fname)
        break

if actual_maf and os.path.getsize(actual_maf) > 1000:
    print("\n" + "=" * 60)
    print("第3步：MAF 文件内容预览")
    print("=" * 60)

    import csv

    # Handle gzipped files
    if actual_maf.endswith(".gz"):
        import gzip
        fh = gzip.open(actual_maf, "rt", encoding="utf-8", errors="replace")
    else:
        fh = open(actual_maf, "r", encoding="utf-8", errors="replace")

    reader = csv.reader(fh, delimiter="\t")
    header = next(reader)
    print(f"列数: {len(header)}")
    print(f"主要列: {[c for c in header[:30]]}")

    # 统计变异类型
    variant_types = {}
    total = 0
    samples = set()
    for row in reader:
        if len(row) < 10:
            continue
        total += 1
        # Variant_Classification is typically column 9
        vclass = row[9] if len(row) > 9 else "UNKNOWN"
        variant_types[vclass] = variant_types.get(vclass, 0) + 1
        # Tumor_Sample_Barcode is typically column 16
        if len(row) > 16:
            samples.add(row[16][:15])  # first 15 chars = patient ID

    print(f"\n总突变数: {total}")
    print(f"样本数: {len(samples)}")
    print(f"\n突变类型分布:")
    for vtype, count in sorted(variant_types.items(), key=lambda x: -x[1])[:10]:
        print(f"  {vtype}: {count} ({count/total*100:.1f}%)")

    fh.close()

print("\n" + "=" * 60)
print("第10步完成！")
print("=" * 60)
