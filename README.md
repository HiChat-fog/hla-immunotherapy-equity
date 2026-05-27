# HLA Immunotherapy Equity — 全球胃癌免疫治疗的公平性分析

[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Data: TCGA](https://img.shields.io/badge/data-TCGA--STAD-orange.svg)](https://portal.gdc.cancer.gov/)
[![DOI: pending](https://img.shields.io/badge/DOI-pending-lightgrey.svg)]()
[![Preprint](https://img.shields.io/badge/preprint-ChinaXiv-red.svg)]()

**FDA 批准的 HLA 限制型免疫疗法，系统性排除了胃癌负担最重的人群。**  
This repo contains the full computational pipeline and evidence behind this finding.

核心结论一句话：把 HLA 靶点从 A\*02:01 换成 A\*11:01 + A\*24:02，中国胃癌患者的理论覆盖率从 **5% → 60%**，翻 12 倍。

---

## Key Findings / 核心发现

| Metric | A\*02:01 (FDA-targeted) | A\*11:01 (top East Asian) |
|--------|--------------------------|---------------------------|
| Frequency in Chinese GC patients | 5.0% | 46.9% |
| Frequency in Europeans | 26.0% | 5.6% |
| Strong-binding neoantigens (SB) | 11 | **27** (2.5×) |
| KRAS SB peptides | **0** | **5** |
| Driver gene coverage | 4 (all shared) | 6 (incl. **FBXW7-exclusive**) |

- **A\*02:01 对胃癌突变肽段的强结合判定在统计上无法与随机肽库区分** — enrichment 仅 0.85×，A\*11:01 和 A\*03:01 都是 2.1×。
- **贪心优化**：按人群 HLA 频率重排研发优先级，中国患者覆盖率 **5.0% → 60.2%**（12.4×），欧洲仅 26.0% → 37.1%（1.4×）。差距不在研发水平，在运气 — A\*02:01 恰好是欧洲的最优靶点。
- **KRAS G12D/G12V 与 A\*02:01 结构不匹配**，现有疗法对最关键癌基因之一的 KRAS 完全无效 — 而 A\*11:01 能呈递。

<!--
**Visual abstract (below):** left = FDA status quo, right = optimized priority.
-->

---

## Repository Structure / 仓库结构

```
├── analysis/           # Computational pipeline (step 4–13)
│   ├── project_config.py               # Shared config, gene lists, country mappings
│   ├── step4_match_data.py             # TCGA sample ↔ clinical matching
│   ├── step4b_gdc_match.py             # GDC API ethnicity matching
│   ├── step5_immune_analysis.py        # Immune infiltration analysis
│   ├── step6_hla_frequency.py          # HLA frequency data integration
│   ├── step7_neoantigen_prediction.py  # Initial neoantigen prediction
│   ├── step8_official_data_analysis.py # GBD burden trend analysis
│   ├── step9_stomach_cancer_bridge.py  # GLOBOCAN gastric cancer data
│   ├── step10_download_mutations.py    # FireBrowse MAF download
│   ├── step11a0_recurrent_mutations.py # Recurrent mutation screening
│   ├── step11a_peptide_extraction.py   # UniProt peptide extraction (9-11mer)
│   ├── step11_full_matrix_prediction.py # IEDB 10,392-prediction pipeline
│   ├── step12_hla_priority_optimization.py  # Greedy coverage algorithm
│   ├── step13_visualization.py         # Figure generation
│   ├── step_hla_frequency_data.py      # HLA frequency data parsing
│   └── step_sample_labels.py           # Sample metadata processing
├── data/               # Input/output data files
│   ├── neoantigen_matrix_results.csv   # 10,333 IEDB predictions (core result)
│   ├── hla_frequency_data.csv          # Population-level HLA frequencies
│   ├── hla_greedy_priority.json        # Optimization results
│   ├── recurrent_peptides.csv          # 2,598 mutation-spanning peptides
│   └── ...                             # Additional data files
├── figures/            # All 12 figures (PNG, 150 dpi)
├── paper/              # LaTeX manuscript
│   ├── main.tex                        # Full LaTeX source
│   └── main.pdf                        # Compiled PDF
├── results/            # (placeholder for reproduction output)
├── .gitignore
├── requirements.txt    # Python dependencies
├── LICENSE             # MIT
└── README.md
```

## Quick Start / 快速复现

### Prerequisites / 前置条件

- Python 3.8+
- Internet connection (IEDB API, UniProt API are remote services)

### Install / 安装

```bash
git clone https://github.com/HiChat-fog/hla-immunotherapy-equity.git
cd hla-immunotherapy-equity
pip install -r requirements.txt
```

### Reproduce core results / 复现核心结果

The pipeline runs in order. Steps 1–3 are data preprocessing (hand-crafted), steps 4–13 are scripted.

**Full pipeline (from raw data):**

```bash
cd analysis

# 1. Download TCGA-STAD mutation data
python step10_download_mutations.py

# 2. Screen recurrent mutations (≥2 patients)
python step11a0_recurrent_mutations.py

# 3. Extract mutation-spanning peptides (UniProt API)
python step11a_peptide_extraction.py

# 4. Run full-matrix IEDB binding prediction (~4 hours)
python step11_full_matrix_prediction.py

# 5. Greedy coverage optimization
python step12_hla_priority_optimization.py

# 6. Generate figures
python step13_visualization.py
```

**Quick start with pre-computed data:**  
The `data/` directory already contains `neoantigen_matrix_results.csv` (10,333 IEDB predictions).  
Skip directly to the optimization and visualization:

```bash
cd analysis
python step12_hla_priority_optimization.py
python step13_visualization.py
# Figures appear in ../figures/
```

---

## Scientific Background / 科学背景

### The Problem

Tebentafusp (2022) and Afami-cel (2024) — the first FDA-approved TCR-T/CAR-T therapies for solid tumors — both require **HLA-A\*02:01** for antigen presentation. This allele is:

- **~27%** in European populations
- **~5%** in Chinese gastric cancer patients (Zhou et al., 2019)

Meanwhile, **A\*11:01** (46.9%) and **A\*24:02** (25.0%) together cover >70% of Chinese gastric cancer patients — but have **zero** corresponding drugs in development.

### What We Did

1. Downloaded **178,508 mutation records** from TCGA-STAD (FireBrowse)
2. Filtered to recurrent non-silent mutations → **2,598 unique mutation-spanning peptides** (9–11mer)
3. Submitted **10,392 IEDB MHC-I binding predictions** across 4 HLA alleles (A\*02:01, A\*11:01, A\*24:02, A\*03:01)
4. Ran a **greedy coverage algorithm** to determine optimal HLA-targeting priority by population

### What We Found

**A\*11:01 is not just more frequent — it is functionally superior for gastric cancer neoantigens.**

- SB peptides: A\*11:01 (27) = A\*03:01 (27) > A\*24:02 (13) > A\*02:01 (11)
- A\*02:01's 11 SB is *below* the random expectation of 13 (0.5% × 2,585 peptides) — **enrichment ratio 0.85×**
- A\*02:01 has **zero** strong-binding peptides for KRAS, CDH1, SMAD4, FBXW7, or ERBB3
- A\*11:01 has **exclusive** SB coverage of FBXW7 and co-covers KRAS with A\*03:01

**The coverage gap is extreme:**

| Scenario | Chinese patients | European patients |
|----------|-----------------|-------------------|
| FDA status quo | 5.0% | 26.0% |
| Optimized (top-2 HLA) | **60.2%** | 37.1% |
| Improvement | **12.4×** | 1.4× |

---

## Data Sources / 数据来源

| Data | Source | Description |
|------|--------|-------------|
| Global death causes (2000–2023) | GBD 2021 | 203 countries, cause-stratified |
| Gastric cancer burden | GLOBOCAN 2022 (IARC) | 185 countries, ASR |
| Gastric cancer mutations | TCGA-STAD (FireBrowse) | 178,508 records, 393 samples |
| Protein sequences | UniProt REST API | ~100 human proteins |
| HLA-peptide binding | IEDB MHC-I API (NetMHCpan-4.1) | 10,392 predictions |
| HLA population frequencies | AFND / NMDP / Zhou et al. 2019 | 25 records |
| Economic indicators | World Bank WDI | GDP, health expenditure |

---

## Citation / 引用

If you use this work, please cite:

```bibtex
@misc{hichat2026hla,
  author = {HiChat},
  title = {Global Equity Analysis of HLA-Restricted Immunotherapy},
  year = {2026},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/HiChat-fog/hla-immunotherapy-equity}}
}
```

## Author / 作者

**HiChat** — 本科独立研究，生信/免疫信息学方向。

GitHub: [@HiChat-fog](https://github.com/HiChat-fog)

---

## License / 许可

MIT License — see [LICENSE](LICENSE) for details.

---

*This project originated from the 2026 Chinese University Computer Design Competition and was later revised for public release.*
