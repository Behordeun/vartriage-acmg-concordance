"""Analysis 5: Functional consequence enrichment.

Examines the distribution of variant consequences across expert classification
tiers, quantifying whether certain consequence types are enriched in pathogenic
vs benign expert assertions.

Analyses:
- Consequence distribution by expert tier (chi-squared test)
- Odds ratios for pathogenicity per consequence type
- Comparison of VEP consequence vs simplified annotation

Input: eRepo stratified results from the eRepo validation run
Output: results/consequence_enrichment.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import chi2_contingency, fisher_exact

RESULTS_DIR = Path(__file__).parent.parent / "results"
VALIDATION_RESULTS = Path(__file__).parent.parent / "data" / "validation_results"


def load_stratified() -> dict:
    """Load per-consequence stratified results."""
    p = VALIDATION_RESULTS / "erepo_stratified.json"
    if not p.exists():
        p = Path(__file__).parent.parent.parent / "vartriage-streaming-acmg" / "results" / "erepo_stratified.json"
    if not p.exists():
        print("ERROR: Cannot find erepo_stratified.json")
        sys.exit(1)
    return json.loads(p.read_text())


def _compute_or_entry(
    plp: int, blb: int, other_plp: int, other_blb: int, count: int
) -> dict:
    if blb == 0 or other_plp == 0:
        or_val: float = float("inf") if plp > 0 else 0.0
        p_val = None
    else:
        table = np.array([[plp, blb], [other_plp, other_blb]])
        denom = blb * other_plp
        or_val = (plp * other_blb) / denom if denom > 0 else float("inf")
        _, p_val = fisher_exact(table, alternative="two-sided")
    return {
        "expert_P_LP": plp,
        "expert_B_LB": blb,
        "total": count,
        "odds_ratio": round(or_val, 2) if or_val != float("inf") else "inf",
        "p_value": round(p_val, 6) if p_val is not None else None,
        "pathogenic_fraction": round(plp / count, 4) if count > 0 else 0,
    }


def compute_odds_ratios(stratified: dict) -> dict:
    """Compute odds ratio for pathogenicity per consequence type.

    OR = (P_LP in this consequence / B_LB in this consequence) /
         (P_LP in other consequences / B_LB in other consequences)
    """
    relaxed = stratified["relaxed"]
    total_plp = sum(d["expert_P_LP"] for d in relaxed.values())
    total_blb = sum(d["expert_B_LB"] for d in relaxed.values())
    return {
        csq: _compute_or_entry(
            data["expert_P_LP"],
            data["expert_B_LB"],
            total_plp - data["expert_P_LP"],
            total_blb - data["expert_B_LB"],
            data["count"],
        )
        for csq, data in relaxed.items()
    }


def compute_chi_squared(stratified: dict) -> dict:
    """Chi-squared test for independence between consequence type and expert assertion."""
    relaxed = stratified["relaxed"]
    consequences = list(relaxed.keys())

    # Contingency table: rows = consequences, cols = [P_LP, B_LB, VUS]
    table = []
    for csq in consequences:
        d = relaxed[csq]
        vus = d["count"] - d["expert_P_LP"] - d["expert_B_LB"]
        table.append([d["expert_P_LP"], d["expert_B_LB"], max(0, vus)])

    table = np.array(table)
    chi2, p_val, dof, _ = chi2_contingency(table)

    return {
        "chi2_statistic": round(float(chi2), 2),
        "p_value": float(p_val),
        "degrees_of_freedom": int(dof),
        "significant": bool(p_val < 0.001),
        "interpretation": "Consequence type and expert classification are not independent" if p_val < 0.001 else "No significant association",
    }


def main() -> None:
    stratified = load_stratified()

    odds_ratios = compute_odds_ratios(stratified)
    chi_sq = compute_chi_squared(stratified)

    # Consequence distribution summary
    distribution = stratified.get("consequence_distribution", {})

    results = {
        "consequence_distribution": distribution,
        "odds_ratios_for_pathogenicity": odds_ratios,
        "chi_squared_test": chi_sq,
        "key_findings": {
            "highest_pathogenic_enrichment": max(
                odds_ratios.keys(),
                key=lambda k: odds_ratios[k]["pathogenic_fraction"],
            ),
            "most_common_consequence": max(distribution, key=distribution.get) if distribution else None,
        },
    }

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    output = RESULTS_DIR / "consequence_enrichment.json"
    output.write_text(json.dumps(results, indent=2))

    print(f"Consequence enrichment saved to {output}")
    print(f"  Chi-squared: {chi_sq['chi2_statistic']} (p={chi_sq['p_value']:.2e})")
    print(f"\n  {'Consequence':<15} {'P/LP':>6} {'B/LB':>6} {'OR':>8} {'Path%':>7}")
    print("  " + "-" * 45)
    for csq, data in sorted(odds_ratios.items(), key=lambda x: x[1]["pathogenic_fraction"], reverse=True):
        print(f"  {csq:<15} {data['expert_P_LP']:>6} {data['expert_B_LB']:>6} {str(data['odds_ratio']):>8} {data['pathogenic_fraction']*100:>6.1f}%")


if __name__ == "__main__":
    main()
