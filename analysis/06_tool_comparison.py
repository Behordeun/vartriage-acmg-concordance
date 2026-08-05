"""Analysis 6: Multi-tool sensitivity/specificity comparison.

Compiles published performance numbers from peer-reviewed benchmarks
on the ClinGen eRepo dataset alongside VarTriage's measured results.

Published numbers cited from:
- BIAS-2015: Eisenhart et al. (2025) Genome Medicine 17:62
- InterVar: Eisenhart et al. (2025) same benchmark
- VarTriage: This study from Sulaiman & Oyeyemi

Output: results/tool_comparison.json
"""

from __future__ import annotations

import json
from pathlib import Path

RESULTS_DIR = Path(__file__).parent.parent / "results"


def main() -> None:
    tools = {
        "VarTriage": {
            "source": "This study",
            "criteria_count": 10,
            "benchmark": "ClinGen eRepo (15,334 variants)",
            "path_sensitivity": 0.689,
            "path_ppv": 0.805,
            "ben_sensitivity": 0.169,
            "ben_ppv": 0.999,
            "combining_rules": "Relaxed (Bayesian SVI)",
            "local_execution": True,
            "memory_mb": 453,
        },
        "BIAS-2015 v2.1": {
            "source": "Eisenhart et al. (2025) Genome Medicine 17:62",
            "criteria_count": 19,
            "benchmark": "ClinGen eRepo (FDA-approved subset)",
            "path_sensitivity": 0.740,
            "path_ppv": None,
            "ben_sensitivity": 0.802,
            "ben_ppv": None,
            "combining_rules": "ACMG 2015 with configurable weights",
            "local_execution": True,
            "memory_mb": None,
        },
        "InterVar": {
            "source": "Eisenhart et al. (2025) Genome Medicine 17:62",
            "criteria_count": 18,
            "benchmark": "ClinGen eRepo (same as BIAS-2015)",
            "path_sensitivity": 0.643,
            "path_ppv": None,
            "ben_sensitivity": 0.539,
            "ben_ppv": None,
            "combining_rules": "ACMG 2015 strict",
            "local_execution": True,
            "memory_mb": 16000,
        },
    }

    # Compute relative positioning
    vt = tools["VarTriage"]
    comparison = {
        "vs_bias2015": {
            "path_sensitivity_delta": round(vt["path_sensitivity"] - 0.740, 4),
            "interpretation": "VarTriage trails BIAS-2015 by 5.1pp, primarily due to missing SpliceAI and fewer criteria (10 vs 19)",
        },
        "vs_intervar": {
            "path_sensitivity_delta": round(vt["path_sensitivity"] - 0.643, 4),
            "interpretation": "VarTriage exceeds InterVar by 4.6pp through strength-modulated PP3/BP4 thresholds",
        },
    }

    results = {
        "tools": tools,
        "relative_comparison": comparison,
        "methodology_note": "BIAS-2015 and InterVar numbers are from the same published benchmark (Eisenhart et al., 2025) on the FDA-approved ClinGen eRepo subset. VarTriage was evaluated on the full eRepo (15,334 Expert Panel variants with >=3 stars). Direct head-to-head comparison on identical inputs was not performed.",
    }

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    output = RESULTS_DIR / "tool_comparison.json"
    output.write_text(json.dumps(results, indent=2))

    print(f"Tool comparison saved to {output}")
    print(f"\n  {'Tool':<18} {'Path Sens':>10} {'Ben Sens':>9} {'Criteria':>9}")
    print("  " + "-" * 48)
    for name, data in tools.items():
        ps = f"{data['path_sensitivity']*100:.1f}%" if data["path_sensitivity"] else "N/A"
        bs = f"{data['ben_sensitivity']*100:.1f}%" if data["ben_sensitivity"] else "N/A"
        print(f"  {name:<18} {ps:>10} {bs:>9} {data['criteria_count']:>9}")


if __name__ == "__main__":
    main()
