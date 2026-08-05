# VarTriage ACMG Concordance Analysis

Reproducible analysis and manuscript source for:

> **Concordance of automated ACMG variant classification with expert-curated assertions: a systematic evaluation using the ClinGen Evidence Repository**
>
> Bolaji Fatai Oyeyemi, Muhammad Abiodun Sulaiman
>
> Target: BMC Genomics, 2026

## Repository Structure

```text
vartriage-acmg-concordance/
├── analysis/                      Analysis scripts
│   ├── 01_concordance.py              ClinGen eRepo concordance metrics
│   ├── 02_roc_calibration.py          REVEL/CADD ROC curves and threshold calibration
│   ├── 03_tag_contribution.py         Evidence tag frequency, co-occurrence, ablation
│   ├── 04_missing_data_ablation.py    Classification accuracy vs removed data sources
│   ├── 05_consequence_enrichment.py   Functional consequence distribution analysis
│   ├── 06_tool_comparison.py          Multi-tool sensitivity/specificity comparison
│   └── figures.py                     Publication figure generation
├── data/                          Input data (symlinked from vartriage-streaming-acmg)
│   ├── erepo/
│   │   ├── clingen_erepo.tsv
│   │   └── .api_cache.db
│   └── references/
│       └── revel_grch38.tsv
├── results/                       Analysis outputs (JSON)
├── figures/                       Generated publication figures (PDF)
├── manuscript/                    LaTeX manuscript source
├── reproduce.sh                   One-command reproduction
└── pyproject.toml                 Pinned dependencies
```

## Quick Start

```bash
# Set up environment
uv sync

# Symlink shared data from the validation repo (avoids duplication)
ln -sf ../vartriage-streaming-acmg/data/erepo data/erepo
ln -sf ../vartriage-streaming-acmg/data/references data/references

# Run all analyses
./reproduce.sh

# Or regenerate figures only
./reproduce.sh --quick
```

## Analyses

| # | Analysis | Script | Output |
| --- | ---------- | -------- | -------- |
| 1 | eRepo concordance | `01_concordance.py` | Sensitivity, PPV, Cohen's kappa |
| 2 | Predictor calibration | `02_roc_calibration.py` | ROC curves, optimal thresholds |
| 3 | Evidence tag contribution | `03_tag_contribution.py` | Tag frequency, co-occurrence, marginal contribution |
| 4 | Missing data ablation | `04_missing_data_ablation.py` | Accuracy degradation per removed source |
| 5 | Consequence enrichment | `05_consequence_enrichment.py` | Consequence distribution by expert assertion |
| 6 | Tool comparison | `06_tool_comparison.py` | Multi-tool sensitivity bar chart |

## Data

Uses the same 15,334 ClinGen eRepo Expert Panel variants from the initial eRepo validation. Data files are symlinked to avoid duplication. The `.api_cache.db` contains cached gnomAD and VEP responses from prior runs.

## Software Versions

| Component | Version |
| ----------- | --------- |
| Python | 3.11.14 |
| vartriage | 0.14.0 |
| scikit-learn | >=1.3 |
| Hardware | Apple M3 Pro, 36 GB |

## License

MIT
