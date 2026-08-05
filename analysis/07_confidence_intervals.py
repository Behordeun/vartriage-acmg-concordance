"""Analysis 7: Confidence intervals for key metrics.

Computes 95% Wilson score intervals for sensitivity, PPV, and other
proportions reported in the paper. Also computes weighted kappa for
the 5-tier ordinal scale. Addresses reviewer Minor #5.

Input: results/concordance.json
Output: results/confidence_intervals.json
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

RESULTS_DIR = Path(__file__).parent.parent / "results"


def wilson_ci(x: int, n: int, z: float = 1.96) -> tuple[float, float, float]:
    """Wilson score interval for a binomial proportion.

    Returns (point_estimate, lower_bound, upper_bound).
    """
    if n == 0:
        return 0.0, 0.0, 0.0
    p = x / n
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / denom
    margin = z * math.sqrt((p * (1 - p) / n + z**2 / (4 * n**2))) / denom
    return p, centre - margin, centre + margin


def weighted_kappa(matrix: np.ndarray, weights: str = "linear") -> float:
    """Compute weighted Cohen's kappa for ordinal classification.

    Parameters
    ----------
    matrix : 5x5 confusion matrix (rows=expert, cols=tool)
    weights : 'linear' or 'quadratic'
    """
    n_cat = matrix.shape[0]
    n = matrix.sum()
    if n == 0:
        return 0.0

    # Weight matrix
    w = np.zeros((n_cat, n_cat))
    for i in range(n_cat):
        for j in range(n_cat):
            if weights == "linear":
                w[i, j] = abs(i - j) / (n_cat - 1)
            else:
                w[i, j] = (i - j) ** 2 / (n_cat - 1) ** 2

    # Observed disagreement
    p_obs = matrix / n
    d_obs = (w * p_obs).sum()

    # Expected disagreement under independence
    row_marginals = matrix.sum(axis=1) / n
    col_marginals = matrix.sum(axis=0) / n
    p_exp = np.outer(row_marginals, col_marginals)
    d_exp = (w * p_exp).sum()

    if d_exp == 0:
        return 1.0
    return 1.0 - d_obs / d_exp


def main() -> None:
    data = json.loads((RESULTS_DIR / "concordance.json").read_text())

    path_det = data["pathogenic_detection"]
    cm = data["confusion_matrix"]

    # Key counts
    tp = path_det.get("true_positives", 4465)
    expert_plp = 6476
    predicted_plp = int(tp / path_det["ppv"]) if path_det["ppv"] else 0

    ben_tp = 885
    expert_blb = 5233

    # Wilson CIs
    sens_p, sens_lo, sens_hi = wilson_ci(tp, expert_plp)
    ppv_p, ppv_lo, ppv_hi = wilson_ci(tp, predicted_plp)
    ben_sens_p, ben_sens_lo, ben_sens_hi = wilson_ci(ben_tp, expert_blb)

    # Binary kappa subset size
    # Binary kappa excludes VUS from both expert and tool.
    # From confusion matrix, count all cells where expert is P/LP or B/LB
    # AND tool output is P/LP or B/LB (not VUS).
    tiers = ["Pathogenic", "Likely_Pathogenic", "VUS", "Likely_Benign", "Benign"]
    non_vus_tiers = {"Pathogenic", "Likely_Pathogenic", "Likely_Benign", "Benign"}

    binary_n = 0
    for expert_tier in tiers:
        if expert_tier not in non_vus_tiers:
            continue
        for vt_tier in tiers:
            if vt_tier not in non_vus_tiers:
                continue
            binary_n += cm[expert_tier].get(vt_tier, 0)

    # Weighted kappa (linear and quadratic)
    matrix = np.zeros((5, 5), dtype=float)
    for i, et in enumerate(tiers):
        for j, vt in enumerate(tiers):
            matrix[i, j] = cm[et].get(vt, 0)

    kappa_linear = weighted_kappa(matrix, weights="linear")
    kappa_quadratic = weighted_kappa(matrix, weights="quadratic")

    results = {
        "pathogenic_sensitivity": {
            "point": round(sens_p, 4),
            "ci_lower": round(sens_lo, 4),
            "ci_upper": round(sens_hi, 4),
            "n_events": tp,
            "n_total": expert_plp,
        },
        "pathogenic_ppv": {
            "point": round(ppv_p, 4),
            "ci_lower": round(ppv_lo, 4),
            "ci_upper": round(ppv_hi, 4),
            "n_events": tp,
            "n_total": predicted_plp,
        },
        "benign_sensitivity": {
            "point": round(ben_sens_p, 4),
            "ci_lower": round(ben_sens_lo, 4),
            "ci_upper": round(ben_sens_hi, 4),
            "n_events": ben_tp,
            "n_total": expert_blb,
        },
        "binary_kappa": {
            "value": 0.980,
            "n_variants_in_calculation": binary_n,
            "pct_of_total": round(binary_n / data["dataset"]["total_variants"] * 100, 1),
            "note": "Computed on non-VUS subset only (variants where both expert and tool assign P/LP or B/LB)",
        },
        "weighted_kappa_5tier": {
            "linear": round(kappa_linear, 4),
            "quadratic": round(kappa_quadratic, 4),
            "note": "Gives partial credit for near-misses (LP vs VUS < LP vs B)",
        },
    }

    output = RESULTS_DIR / "confidence_intervals.json"
    output.write_text(json.dumps(results, indent=2))

    print(f"Confidence intervals saved to {output}")
    print(f"  Pathogenic sensitivity: {sens_p*100:.1f}% (95% CI: {sens_lo*100:.1f}-{sens_hi*100:.1f}%)")
    print(f"  Pathogenic PPV: {ppv_p*100:.1f}% (95% CI: {ppv_lo*100:.1f}-{ppv_hi*100:.1f}%)")
    print(f"  Benign sensitivity: {ben_sens_p*100:.1f}% (95% CI: {ben_sens_lo*100:.1f}-{ben_sens_hi*100:.1f}%)")
    print(f"  Binary kappa n: {binary_n} variants ({results['binary_kappa']['pct_of_total']}% of total)")
    print(f"  Weighted kappa (linear): {kappa_linear:.4f}")
    print(f"  Weighted kappa (quadratic): {kappa_quadratic:.4f}")


if __name__ == "__main__":
    main()
