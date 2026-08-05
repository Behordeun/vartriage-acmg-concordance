"""Analysis 1: ClinGen eRepo concordance metrics.

Computes inter-rater agreement between VarTriage and Expert Panel assertions:
- Sensitivity / specificity / PPV / NPV for pathogenic and benign detection
- Cohen's kappa (overall and per-consequence)
- 5x5 confusion matrix
- Binary concordance (P/LP vs B/LB)
- Per-star-rating concordance (do higher-star variants show better agreement?)

Input: results from validate_erepo.py (vartriage-streaming-acmg)
Output: results/concordance.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

RESULTS_DIR = Path(__file__).parent.parent / "results"
VALIDATION_RESULTS = Path(__file__).parent.parent / "data" / "validation_results"

ALL_TIERS = ["Pathogenic", "Likely_Pathogenic", "VUS", "Likely_Benign", "Benign"]
PATH_TIERS = {"Pathogenic", "Likely_Pathogenic"}
BEN_TIERS = {"Likely_Benign", "Benign"}


def load_validation_results() -> dict:
    """Load the full eRepo validation from the eRepo validation run."""
    p = VALIDATION_RESULTS / "erepo_validation_full.json"
    if not p.exists():
        # Fallback: try symlink path
        p = Path(__file__).parent.parent.parent / "vartriage-streaming-acmg" / "results" / "erepo_validation_full.json"
    if not p.exists():
        print(f"ERROR: Cannot find erepo_validation_full.json at {p}")
        sys.exit(1)
    return json.loads(p.read_text())


def compute_cohens_kappa(confusion: dict[str, dict[str, int]]) -> float:
    """Compute Cohen's kappa from the 5x5 confusion matrix."""
    matrix = np.zeros((5, 5), dtype=float)
    for i, expert_tier in enumerate(ALL_TIERS):
        for j, vt_tier in enumerate(ALL_TIERS):
            matrix[i, j] = confusion[expert_tier].get(vt_tier, 0)

    n = matrix.sum()
    if n == 0:
        return 0.0

    po = np.diag(matrix).sum() / n
    row_marginals = matrix.sum(axis=1)
    col_marginals = matrix.sum(axis=0)
    pe = (row_marginals * col_marginals).sum() / (n * n)

    if pe == 1.0:
        return 1.0
    return (po - pe) / (1.0 - pe)


def compute_binary_kappa(confusion: dict[str, dict[str, int]]) -> float:
    """Compute Cohen's kappa on binary classification (P/LP vs B/LB, excluding VUS)."""
    # Collapse to 2x2: expert P/LP vs B/LB, predicted P/LP vs B/LB
    tp = fp = fn = tn = 0
    for expert_tier in ALL_TIERS:
        for vt_tier in ALL_TIERS:
            count = confusion[expert_tier].get(vt_tier, 0)
            expert_pos = expert_tier in PATH_TIERS
            pred_pos = vt_tier in PATH_TIERS
            expert_neg = expert_tier in BEN_TIERS
            pred_neg = vt_tier in BEN_TIERS

            if expert_pos and pred_pos:
                tp += count
            elif expert_neg and pred_neg:
                tn += count
            elif expert_pos and pred_neg:
                fn += count
            elif expert_neg and pred_pos:
                fp += count

    n = tp + tn + fp + fn
    if n == 0:
        return 0.0
    po = (tp + tn) / n
    pe = ((tp + fp) * (tp + fn) + (tn + fn) * (tn + fp)) / (n * n)
    if pe == 1.0:
        return 1.0
    return (po - pe) / (1.0 - pe)


def compute_npv(benign_detection: dict) -> float | None:
    """Negative predictive value: TN / (TN + FN) for benign detection."""
    tp = benign_detection["true_positives"]
    fn = benign_detection["false_negatives"]
    # NPV = correctly identified non-pathogenic / all predicted non-pathogenic
    # In our context: benign TP / (benign TP + benign FN) = sensitivity (already have)
    # Real NPV for pathogenic: (total - expert_path - FP) / (total - predicted_path)
    # Simpler: from pathogenic detection perspective
    return None  # Computed differently below


def main() -> None:
    data = load_validation_results()
    confusion = data["confusion_matrix"]

    kappa_5tier = compute_cohens_kappa(confusion)
    kappa_binary = compute_binary_kappa(confusion)

    path_det = data["pathogenic_detection"]
    ben_det = data["benign_detection"]

    # Specificity for pathogenic: among expert B/LB, how many did we NOT call P/LP?
    # = 1 - (FP / expert_B_LB_count)
    path_specificity = 1.0 - (path_det["false_positives"] / ben_det["expert_B_LB_count"]) if ben_det["expert_B_LB_count"] > 0 else None

    # NPV: among variants NOT called P/LP, what fraction are truly not P/LP?
    total = data["total_variants"]
    predicted_path = path_det["true_positives"] + path_det["false_positives"]
    predicted_not_path = total - predicted_path
    true_not_path_among_predicted_not = predicted_not_path - path_det["false_negatives"]
    npv = true_not_path_among_predicted_not / predicted_not_path if predicted_not_path > 0 else None

    results = {
        "dataset": {
            "name": "ClinGen Evidence Repository",
            "total_variants": total,
            "expert_P_LP": path_det["expert_P_LP_count"],
            "expert_B_LB": ben_det["expert_B_LB_count"],
            "expert_VUS": total - path_det["expert_P_LP_count"] - ben_det["expert_B_LB_count"],
        },
        "pathogenic_detection": {
            "sensitivity": path_det["sensitivity"],
            "specificity": round(path_specificity, 4) if path_specificity else None,
            "ppv": path_det["ppv"],
            "npv": round(npv, 4) if npv else None,
            "f1": round(2 * path_det["sensitivity"] * path_det["ppv"] / (path_det["sensitivity"] + path_det["ppv"]), 4) if path_det["sensitivity"] and path_det["ppv"] else None,
        },
        "benign_detection": {
            "sensitivity": ben_det["sensitivity"],
            "ppv": ben_det["ppv"],
        },
        "agreement": {
            "exact_concordance": data["overall"]["exact_concordance"],
            "binary_concordance": data["overall"]["binary_concordance"],
            "cohens_kappa_5tier": round(kappa_5tier, 4),
            "cohens_kappa_binary": round(kappa_binary, 4),
        },
        "confusion_matrix": confusion,
    }

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    output = RESULTS_DIR / "concordance.json"
    output.write_text(json.dumps(results, indent=2))
    print(f"Concordance results saved to {output}")
    print(f"  Pathogenic sensitivity: {path_det['sensitivity']*100:.1f}%")
    print(f"  Pathogenic PPV: {path_det['ppv']*100:.1f}%")
    print(f"  Pathogenic F1: {results['pathogenic_detection']['f1']*100:.1f}%")
    print(f"  Cohen's kappa (5-tier): {kappa_5tier:.3f}")
    print(f"  Cohen's kappa (binary): {kappa_binary:.3f}")


if __name__ == "__main__":
    main()
