# ACMG Concordance Analysis: VarTriage vs ClinGen Expert Panels

Reproducible analysis and manuscript source for:

> Bolaji Fatai Oyeyemi, Muhammad Abiodun Sulaiman. **Concordance of Automated ACMG Variant Classification with Expert-Curated Assertions: A Systematic Evaluation Using the ClinGen Evidence Repository.** Manuscript in preparation (2026).

---

## Key Findings

- **Pathogenic sensitivity:** 68.9% (95% CI: 67.8-70.1%), PPV 80.5%
- **Binary kappa:** 0.98 (near-perfect pathogenic/benign discrimination)
- **Largest accuracy driver:** REVEL scores + Bayesian relaxed combining rule (+36.9pp sensitivity)
- **Principal gap:** Benign sensitivity (16.9%) -- caused by missing BS3/BS4/BP6 evidence types
- **Consequence enrichment:** Nonsense OR = 160.25 for pathogenicity (chi-squared = 9,117)

---

## Reproduction

```bash
# Prerequisites: Python 3.11+, uv package manager
# The vartriage-streaming-acmg repo must be adjacent for shared data

uv sync

# Symlink shared data (run once)
ln -sf ../vartriage-streaming-acmg/data/erepo data/erepo
ln -sf ../vartriage-streaming-acmg/data/references data/references

# Run all analyses + generate figures
./reproduce.sh

# Figures only (from existing results)
./reproduce.sh --quick
```

---

## Analyses

| # | Script | What it produces |
| --- | -------- | ----------------- |
| 1 | `01_concordance.py` | Sensitivity, PPV, Cohen's kappa, confusion matrix |
| 2 | `02_roc_calibration.py` | REVEL ROC curve, AUC, ClinGen threshold validation |
| 3 | `03_tag_contribution.py` | Evidence tag frequency, relaxed vs strict combining impact |
| 4 | `04_missing_data_ablation.py` | Sensitivity per removed data source |
| 5 | `05_consequence_enrichment.py` | Chi-squared, odds ratios by consequence type |
| 6 | `06_tool_comparison.py` | VarTriage vs BIAS-2015 vs InterVar |
| 7 | `07_confidence_intervals.py` | Wilson CIs, weighted kappa |
| -- | `figures.py` | 6 publication-ready PDF figures |

---

## Repository Layout

```text
analysis/           Analysis scripts (numbered for execution order)
results/            JSON outputs from each analysis
  ablation/         Per-configuration ablation results
figures/            Publication figures (vector PDF)
manuscript/
  main.md           Markdown manuscript (source of truth)
  main.tex          LaTeX main file
  sections/         LaTeX section files (methods, results, discussion)
  references.bib    BibTeX (68 entries, all with DOIs)
  response_to_reviewers.md
data/               Symlinked from vartriage-streaming-acmg (gitignored)
reproduce.sh        One-command pipeline
pyproject.toml      Pinned dependencies
```

---

## Figures

| Figure | Content |
| -------- | --------- |
| Fig 1 | 5x5 confusion matrix heatmap (row-normalized %) |
| Fig 2 | REVEL ROC curve with ClinGen threshold markers |
| Fig 3 | Ablation waterfall (sensitivity change per source) |
| Fig 4 | Evidence tag frequency bar chart |
| Fig 5 | Consequence odds ratios (log-scale forest plot) |
| Fig 6 | Multi-tool sensitivity comparison |

---

## Dependencies

| Package | Version |
| --------- | --------- |
| Python | 3.11.14 |
| vartriage | 0.14.0 |
| numpy | >= 1.24 |
| scipy | >= 1.11 |
| scikit-learn | >= 1.3 |
| matplotlib | >= 3.8 |
| seaborn | >= 0.13 |

Install all via `uv sync`.

---

## Citation

```bibtex
@article{oyeyemi2026concordance,
  author={Oyeyemi, Bolaji Fatai and Sulaiman, Muhammad Abiodun},
  title={Concordance of Automated ACMG Variant Classification with Expert-Curated
         Assertions: A Systematic Evaluation Using the ClinGen Evidence Repository},
  year={2026},
  note={Manuscript in preparation}
}
```

---

## Related

- [vartriage](https://github.com/Behordeun/vartriage) -- The VarTriage library
- [vartriage-streaming-acmg](https://github.com/Behordeun/vartriage-streaming-acmg) -- VarTriage validation and benchmarking

## License

MIT
