"""
============================================================
  项目第8步：官方死亡原因数据分析
============================================================

  做什么：
  1. 加载24年全球死亡原因数据（比赛官方数据集）
  2. 聚焦：肿瘤 + 消化系统疾病
  3. 分析：时间趋势、区域差异、与经济水平的关系
  4. 产出：可视化图表

  定位：整个项目的基础层 —— "看到问题全景"
"""
import csv
import os
import sys

import numpy as np
from project_config import CN_TO_EN, REGIONS, FIG_DIR

OUTPUT_DIR = os.path.dirname(__file__)
DATA_DIR = os.path.join(OUTPUT_DIR, "死亡原因")

# ---- 常量配置 ----
YEAR_START = 2000
YEAR_END = 2024          # range 不包含此值
DISPLAY_YEARS = [2000, 2005, 2010, 2015, 2020, 2023]
COMPARE_START = 2000
COMPARE_END = 2023
GDP_YEAR = 2019           # GDP匹配用的年份
GDP_LOW = 3000            # 低收入阈值 (世界银行标准)
GDP_HIGH = 12000          # 高收入阈值 (世界银行标准)
KEY_COUNTRIES_EN = ["China", "United States", "Japan", "India", "Germany",
                     "Russian Federation", "Nigeria", "South Africa"]

# 用 matplotlib，如果没有则跳过可视化
try:
    import matplotlib
    matplotlib.use("Agg")  # 无 GUI 模式
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    print("[OK] matplotlib 就绪")
    HAS_PLT = True
except ImportError:
    print("[WARN] matplotlib 未安装，仅做数据汇总，不生成图")
    HAS_PLT = False

# 字体设置：支持中文
if HAS_PLT:
    try:
        plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "DejaVu Sans"]
        plt.rcParams["axes.unicode_minus"] = False
    except:
        pass


# ====== 第1步：加载全部24年数据 ======
print("=" * 60)
print("第1步：加载 2000-2023 年死亡原因数据")
print("=" * 60)

# 死亡人数存储：{年份: {国家: {原因: 死亡人数}}}
mortality = {}
# 国家-年份-原因的排名存储（线性插值时用）
ranking = {}
countries = set()
causes = set()

for year in range(YEAR_START, YEAR_END):
    fpath = os.path.join(DATA_DIR, f"{year}.csv")
    if not os.path.exists(fpath):
        print(f"[WARN] 跳过: {fpath} 不存在")
        continue
    mortality[year] = {}
    ranking[year] = {}
    with open(fpath, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        next(reader)  # 跳过表头
        for row in reader:
            if len(row) < 8:
                continue
            country = row[1]
            cause = row[5]
            measure = row[6]
            value_str = row[7]

            if measure == "死亡":
                try:
                    value = float(value_str)
                except (ValueError, TypeError):
                    continue

                countries.add(country)
                causes.add(cause)
                if country not in mortality[year]:
                    mortality[year][country] = {}
                mortality[year][country][cause] = value
            elif measure == "死亡排名":
                try:
                    rank = float(value_str)
                except (ValueError, TypeError):
                    rank = None
                if country not in ranking[year]:
                    ranking[year][country] = {}
                ranking[year][country][cause] = rank

print(f"加载完成: {len(countries)} 个国家, {len(causes)} 个死因类别, {24} 个年份")
print(f"死因类别: {sorted(causes)}")


# ====== 第2步：汇总肿瘤死亡数据 ======
print("\n" + "=" * 60)
print("第2步：肿瘤死亡负担 — 全球趋势与区域差异")
print("=" * 60)

# 全球肿瘤总死亡趋势
global_tumor_by_year = {}
global_digestive_by_year = {}
for year in range(YEAR_START, YEAR_END):
    total_tumor = 0
    total_digestive = 0
    for country in mortality[year]:
        total_tumor += mortality[year][country].get("肿瘤", 0)
        total_digestive += mortality[year][country].get("消化系统疾病", 0)
    global_tumor_by_year[year] = total_tumor
    global_digestive_by_year[year] = total_digestive

print("全球肿瘤死亡趋势:")
for year in DISPLAY_YEARS:
    if year in global_tumor_by_year:
        print(f"  {year}: {global_tumor_by_year[year]/1e6:.2f}M (肿瘤)",
              f"  {global_digestive_by_year[year]/1e6:.2f}M (消化系统)")

# 区域汇总
print("\n各区域肿瘤死亡趋势对比 (2023 vs 2000):")
for region_name, region_countries in REGIONS.items():
    tumor_2023 = 0
    tumor_2000 = 0
    for c in region_countries:
        if c in mortality[COMPARE_END]:
            tumor_2023 += mortality[COMPARE_END][c].get("肿瘤", 0)
        if c in mortality[COMPARE_START]:
            tumor_2000 += mortality[COMPARE_START][c].get("肿瘤", 0)
    if tumor_2000 > 0:
        change = (tumor_2023 - tumor_2000) / tumor_2000 * 100
        print(f"  {region_name:12s}: {tumor_2000/1e6:.2f}M → {tumor_2023/1e6:.2f}M ({change:+.1f}%)")


# ====== 第3步：肿瘤死亡排名分析 ======
print("\n" + "=" * 60)
print("第3步：肿瘤在各国死因中的排名变化")
print("=" * 60)

# 对比 2000 vs 2023 各国肿瘤死亡排名变化
rank_changes = []
for country in countries:
    rank_2000 = ranking.get(COMPARE_START, {}).get(country, {}).get("肿瘤", None)
    rank_2023 = ranking.get(COMPARE_END, {}).get(country, {}).get("肿瘤", None)
    if rank_2000 is not None and rank_2023 is not None:
        rank_changes.append((country, rank_2000, rank_2023, rank_2023 - rank_2000))

rank_changes.sort(key=lambda x: x[3])  # 按排名变化排序（负值=排名上升）

print("肿瘤死因排名上升最多的国家（排名变小 = 肿瘤死亡变得更突出）:")
for country, r2000, r2023, change in rank_changes[:10]:
    print(f"  {country:20s}: 第{r2000:.0f} → 第{r2023:.0f} ({change:+.0f})")

print("\n肿瘤死因排名下降最多的国家:")
for country, r2000, r2023, change in rank_changes[-10:]:
    print(f"  {country:20s}: 第{r2000:.0f} → 第{r2023:.0f} ({change:+.0f})")


# ====== 第4步：经济水平与肿瘤负担 ======
print("\n" + "=" * 60)
print("第4步：经济水平与肿瘤/消化系统疾病死亡的关系")
print("=" * 60)

# 加载经济数据（CN_TO_EN 从 project_config 导入）
econ_data = {}
pop_data = {}
econ_path = os.path.join(OUTPUT_DIR, "world_bank_economic_data.csv")
with open(econ_path, "r", encoding="utf-8-sig") as f:  # utf-8-sig handles BOM
    reader = csv.DictReader(f)
    for row in reader:
        country = row.get("country", "")
        year_str = row.get("year", "")
        gdp_str = row.get("gdp_per_capita", "")
        health_str = row.get("health_expenditure", "")
        pop_str = row.get("population", "")
        if country and year_str:
            try:
                yr = int(year_str)
                gdp = float(gdp_str) if gdp_str else None
                health = float(health_str) if health_str else None
                pop = float(pop_str) if pop_str else None
                if country not in econ_data:
                    econ_data[country] = {}
                econ_data[country][yr] = {"gdp": gdp, "health": health}
                if country not in pop_data:
                    pop_data[country] = {}
                pop_data[country][yr] = pop
            except (ValueError, TypeError):
                continue

# 匹配 GDP 与肿瘤死亡率（通过中文→英文映射）
matches_2019 = []
for country_cn in sorted(countries):
    tumor_death = mortality.get(GDP_YEAR, {}).get(country_cn, {}).get("肿瘤")
    if tumor_death is None or tumor_death == 0:
        continue

    country_en = CN_TO_EN.get(country_cn)
    if country_en is None:
        continue

    econ = econ_data.get(country_en, {}).get(GDP_YEAR)
    pop_val = pop_data.get(country_en, {}).get(GDP_YEAR)

    if econ and econ["gdp"] and econ["gdp"] > 0 and pop_val and pop_val > 0:
        rate_per_100k = tumor_death / pop_val * 100000
        matches_2019.append({
            "country_cn": country_cn,
            "country_en": country_en,
            "gdp": econ["gdp"],
            "health": econ["health"],
            "tumor_rate": rate_per_100k,
            "tumor_total": tumor_death,
            "population": pop_val,
        })

print(f"GDP-肿瘤负担匹配: {len(matches_2019)} 个国家")
matches_2019.sort(key=lambda x: x["gdp"])

# 分组对比
if matches_2019:
    low_income = [m for m in matches_2019 if m["gdp"] < GDP_LOW]
    mid_income = [m for m in matches_2019 if GDP_LOW <= m["gdp"] < GDP_HIGH]
    high_income = [m for m in matches_2019 if m["gdp"] >= GDP_HIGH]

    for label, group in [("Low income (<$3,000)", low_income),
                          ("Middle income ($3,000-12,000)", mid_income),
                          ("High income (>$12,000)", high_income)]:
        if group:
            avg_rate = sum(m["tumor_rate"] for m in group) / len(group)
            print(f"  {label}: {len(group)} countries, avg tumor mortality {avg_rate:.1f}/100k")

# Top/bottom highlights
if len(matches_2019) >= 5:
    print("\nHighest tumor mortality rate:")
    for m in sorted(matches_2019, key=lambda x: x["tumor_rate"], reverse=True)[:5]:
        print(f"  {m['country_cn']}: {m['tumor_rate']:.1f}/100k (GDP ${m['gdp']:.0f})")
    print("Lowest tumor mortality rate:")
    for m in sorted(matches_2019, key=lambda x: x["tumor_rate"])[:5]:
        print(f"  {m['country_cn']}: {m['tumor_rate']:.1f}/100k (GDP ${m['gdp']:.0f})")


# ====== 第5步：保存汇总数据 ======
print("\n" + "=" * 60)
print("第5步：保存分析结果")
print("=" * 60)

# 保存全球年度汇总
summary_path = os.path.join(OUTPUT_DIR, "global_cancer_burden_2000_2023.csv")
with open(summary_path, "w", encoding="utf-8", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["year", "tumor_deaths", "digestive_deaths",
                     "tumor_pct_change_from_2000", "digestive_pct_change_from_2000"])
    base_tumor = global_tumor_by_year.get(2000, 1)
    base_digestive = global_digestive_by_year.get(2000, 1)
    for year in range(YEAR_START, YEAR_END):
        t = global_tumor_by_year.get(year, 0)
        d = global_digestive_by_year.get(year, 0)
        writer.writerow([year, t, d,
                         round((t - base_tumor) / base_tumor * 100, 2),
                         round((d - base_digestive) / base_digestive * 100, 2)])

print(f"已保存: {summary_path}")

# 保存国家-年份-肿瘤死亡面板数据（给后续分析用）
panel_path = os.path.join(OUTPUT_DIR, "country_tumor_panel.csv")
with open(panel_path, "w", encoding="utf-8", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["country", "year", "tumor_deaths", "digestive_deaths", "tumor_rank"])
    for year in range(YEAR_START, YEAR_END):
        for country in sorted(countries):
            tumor = mortality[year].get(country, {}).get("肿瘤", 0)
            digestive = mortality[year].get(country, {}).get("消化系统疾病", 0)
            rank_val = ranking.get(year, {}).get(country, {}).get("肿瘤", "")
            writer.writerow([country, year, tumor, digestive, rank_val])

print(f"已保存: {panel_path}")

# 保存 GDP-肿瘤 匹配数据
gdp_tumor_path = os.path.join(OUTPUT_DIR, "gdp_tumor_matched.csv")
if matches_2019:
    with open(gdp_tumor_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["country_cn", "country_en", "gdp_per_capita", "health_expenditure",
                         "tumor_mortality_per_100k", "tumor_total_deaths", "population"])
        for m in matches_2019:
            writer.writerow([m["country_cn"], m["country_en"], m["gdp"],
                             m.get("health", ""), m["tumor_rate"],
                             m["tumor_total"], m["population"]])
    print(f"已保存: {gdp_tumor_path}")


# ====== 第6步：可视化 ======
if not HAS_PLT:
    print("\n[SKIP] 可视化需要 matplotlib，请安装后重新运行")
    print("=" * 60)
    print("第8步数据汇总完成！")
    sys.exit(0)

print("\n" + "=" * 60)
print("第6步：生成可视化图表")
print("=" * 60)

os.makedirs(FIG_DIR, exist_ok=True)

# ---- 图1：全球肿瘤及消化系统疾病死亡趋势 ({COMPARE_START}-{COMPARE_END}) ----
fig, ax1 = plt.subplots(figsize=(12, 6))

years = list(range(2000, 2024))
tumor_vals = [global_tumor_by_year[y] / 1e6 for y in years]
digestive_vals = [global_digestive_by_year[y] / 1e6 for y in years]

ax1.plot(years, tumor_vals, "o-", color="#E74C3C", linewidth=2, markersize=4, label="肿瘤")
ax1.plot(years, digestive_vals, "s-", color="#F39C12", linewidth=2, markersize=4, label="消化系统疾病")
ax1.set_xlabel("年份", fontsize=12)
ax1.set_ylabel("全球死亡人数（百万）", fontsize=12)
ax1.set_title("全球肿瘤及消化系统疾病死亡趋势（2000-2023）", fontsize=14, fontweight="bold")
ax1.legend(fontsize=11, loc="upper left")
ax1.grid(True, alpha=0.3)

ax1.annotate(f"{tumor_vals[-1]:.2f}M", xy=(2023, tumor_vals[-1]),
             xytext=(5, 5), textcoords="offset points", fontsize=9, color="#E74C3C")
ax1.annotate(f"{tumor_vals[0]:.2f}M", xy=(2000, tumor_vals[0]),
             xytext=(5, -15), textcoords="offset points", fontsize=9, color="#E74C3C")

plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "fig1_global_tumor_trend.png"), dpi=150)
plt.close()
print("  图1: 全球肿瘤死亡趋势")

# ---- 图2：区域肿瘤死亡趋势对比 ----
fig, ax = plt.subplots(figsize=(14, 7))

colors = ["#E74C3C", "#3498DB", "#2ECC71", "#9B59B6", "#F39C12", "#1ABC9C", "#E67E22"]
for idx, (region_name, region_countries) in enumerate(REGIONS.items()):
    region_vals = []
    for year in range(YEAR_START, YEAR_END):
        val = 0
        for c in region_countries:
            val += mortality[year].get(c, {}).get("肿瘤", 0)
        region_vals.append(val / 1e6)
    ax.plot(years, region_vals, "o-", color=colors[idx], linewidth=1.5, markersize=3,
            label=f"{region_name}（{len(region_countries)}国）")

ax.set_xlabel("年份", fontsize=12)
ax.set_ylabel("肿瘤死亡人数（百万）", fontsize=12)
ax.set_title("各区域肿瘤死亡趋势对比（2000-2023）", fontsize=14, fontweight="bold")
ax.legend(fontsize=9, loc="upper left", ncol=2)
ax.grid(True, alpha=0.3)
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "fig2_regional_tumor_trend.png"), dpi=150)
plt.close()
print("  图2: 区域肿瘤死亡趋势对比")

# ---- 图3：散点 + 收入组着色 ----
if matches_2019:
    fig, ax = plt.subplots(figsize=(12, 7))

    def income_group(gdp):
        if gdp < GDP_LOW: return "低收入 (<$3,000)", "#E74C3C"
        elif gdp < GDP_HIGH: return "中等收入 ($3,000-12,000)", "#F39C12"
        else: return "高收入 (>$12,000)", "#3498DB"

    groups = {}
    for m in matches_2019:
        label, color = income_group(m["gdp"])
        if label not in groups:
            groups[label] = {"color": color, "gdps": [], "rates": [], "names": [], "gdp_vals": []}
        groups[label]["gdps"].append(m["gdp"])
        groups[label]["rates"].append(m["tumor_rate"])
        groups[label]["names"].append(m["country_cn"])
        groups[label]["gdp_vals"].append(m["gdp"])

    for label, g in groups.items():
        ax.scatter(g["gdps"], g["rates"], s=45, c=g["color"], alpha=0.55,
                   edgecolors="white", linewidth=0.5, label=label)

    key_countries = KEY_COUNTRIES_EN
    for m in matches_2019:
        if m["country_en"] in key_countries:
            _, kc = income_group(m["gdp"])
            ax.scatter([m["gdp"]], [m["tumor_rate"]], s=60, c=kc, alpha=0.85,
                       edgecolors="black", linewidth=1.2, zorder=10)
            ax.annotate(m["country_cn"], (m["gdp"], m["tumor_rate"]),
                        fontsize=9, fontweight="bold",
                        xytext=(6, 6), textcoords="offset points")

    ax.set_xlabel("人均GDP（美元，2019年）", fontsize=12)
    ax.set_ylabel("肿瘤死亡率（每10万人，2019年）", fontsize=12)
    ax.set_title("各国肿瘤死亡率与人均GDP的关系（按收入分组）", fontsize=14, fontweight="bold")
    ax.set_xscale("log")
    ax.legend(fontsize=10, loc="upper left")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "fig3_income_gdp_vs_tumor.png"), dpi=150)
    plt.close()
    print("  图3: 散点+收入分组 → figures/fig3_income_gdp_vs_tumor.png")

# ---- 图4：肿瘤死因排名变化 ----
if rank_changes:
    fig, ax = plt.subplots(figsize=(12, 12))

    top20 = rank_changes[:20]
    countries_20 = [r[0][:20] for r in top20]
    ranks_2000 = [r[1] for r in top20]
    ranks_2023 = [r[2] for r in top20]

    y_pos = range(len(countries_20))
    ax.barh([y + 0.2 for y in y_pos], ranks_2000, height=0.4, color="#3498DB",
            alpha=0.7, label="2000年排名")
    ax.barh([y - 0.2 for y in y_pos], ranks_2023, height=0.4, color="#E74C3C",
            alpha=0.7, label="2023年排名")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(countries_20, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("死因排名（越小表示肿瘤死亡越突出）", fontsize=12)
    ax.set_title("肿瘤死因排名上升最多的国家（2000→2023）", fontsize=14, fontweight="bold")
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, axis="x")
    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "fig4_tumor_rank_change.png"), dpi=150)
    plt.close()
    print("  图4: 肿瘤死因排名变化")

print("\n" + "=" * 60)
print("第8步完成！")
print(f"  数据: {summary_path}")
print(f"  数据: {panel_path}")
print(f"  图表: {FIG_DIR}/ 目录下 4 张图")
print("=" * 60)
