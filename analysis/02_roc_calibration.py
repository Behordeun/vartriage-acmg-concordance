"""Analysis 2: REVEL and CADD ROC curves with threshold calibration.

For each variant with a REVEL or CADD score, treat the expert assertion as
ground truth (P/LP = positive, B/LB = negative, VUS excluded) and compute:
- ROC curve with AUC
- Precision-recall curve with average precision
- Optimal threshold (Youden's J statistic)
- Comparison of optimal threshold vs ClinGen-calibrated thresholds

Input: eRepo TSV + REVEL scores TSV + cached VEP annotations
Output: results/roc_calibration.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Optional

import numpy as np
from sklearn.metrics import (
    auc,
    average_precision_score,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)

RESULTS_DIR = Path(__file__).parent.parent / "results"
DATA_DIR = Path(__file__).parent.parent / "data"

PATH_TIERS = {"Pathogenic", "Likely_Pathogenic"}
BEN_TIERS = {"Likely_Benign", "Benign"}

# ClinGen SVI calibrated thresholds (Pejaver et al., 2022)
CLINGEN_THRESHOLDS = {
    "revel_pp3_supporting": 0.644,
    "revel_pp3_moderate": 0.773,
    "revel_bp4_supporting": 0.290,
    "revel_bp4_moderate": 0.183,
}

_CLASSIFICATION_NORM: dict[str, str] = {
    "Pathogenic": "Pathogenic",
    "Pathogenic/Likely_pathogenic": "Pathogenic",
    "Likely_pathogenic": "Likely_Pathogenic",
    "Uncertain_significance": "VUS",
    "Likely_benign": "Likely_Benign",
    "Benign": "Benign",
    "Benign/Likely_benign": "Benign",
}


def load_erepo_with_scores() -> list[dict]:
    """Load eRepo variants and match with REVEL scores."""
    erepo_path = DATA_DIR / "erepo" / "clingen_erepo.tsv"
    revel_path = DATA_DIR / "references" / "revel_grch38.tsv"

    if not erepo_path.exists():
        print(f"ERROR: eRepo not found at {erepo_path}")
        sys.exit(1)
    if not revel_path.exists():
        print(f"ERROR: REVEL not found at {revel_path}")
        sys.exit(1)

    # Load REVEL
    print("Loading REVEL scores...")
    revel: dict[tuple[str, int, str, str], float] = {}
    with open(revel_path) as fh:
        fh.readline()
        for line in fh:
            parts = line.split("\t", 5)
            if len(parts) >= 5:
                try:
                    revel[(parts[0], int(parts[1]), parts[2], parts[3])] = float(parts[4].strip())
                except (ValueError, IndexError):
                    continue
    print(f"  Loaded {len(revel):,} REVEL scores")

    # Load eRepo and join
    variants: list[dict] = []
    with open(erepo_path) as fh:
        fh.readline()
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            ref, alt = parts[2], parts[3]
            if len(ref) != 1 or len(alt) != 1:
                continue
            try:
                stars = int(parts[6])
            except ValueError:
                continue
            if stars < 3:
                continue
            raw_cls = parts[5].split("|")[0]
            normalized = _CLASSIFICATION_NORM.get(raw_cls)
            if normalized is None:
                continue
            try:
                pos = int(parts[1])
            except ValueError:
                continue

            score = revel.get((parts[0], pos, ref, alt))
            variants.append({
                "chrom": parts[0],
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "expert": normalized,
                "revel": score,
            })

    print(f"  Loaded {len(variants):,} eRepo variants")
    scored = [v for v in variants if v["revel"] is not None]
    print(f"  {len(scored):,} have REVEL scores")
    return variants


def compute_roc(variants: list[dict], score_key: str = "revel") -> dict:
    """Compute ROC and PR curves for a given score against expert labels."""
    # Filter to variants with scores and binary labels (exclude VUS)
    labeled = [
        v for v in variants
        if v.get(score_key) is not None and v["expert"] in (PATH_TIERS | BEN_TIERS)
    ]

    if len(labeled) < 10:
        return {"error": f"Too few labeled variants with {score_key} scores: {len(labeled)}"}

    y_true = np.array([1 if v["expert"] in PATH_TIERS else 0 for v in labeled])
    scores = np.array([v[score_key] for v in labeled])

    # ROC
    fpr, tpr, roc_thresholds = roc_curve(y_true, scores)
    roc_auc = roc_auc_score(y_true, scores)

    # Youden's J
    j_scores = tpr - fpr
    optimal_idx = np.argmax(j_scores)
    optimal_threshold = float(roc_thresholds[optimal_idx])

    # Precision-Recall
    precision, recall, pr_thresholds = precision_recall_curve(y_true, scores)
    avg_precision = average_precision_score(y_true, scores)

    # Performance at ClinGen thresholds
    threshold_performance: dict[str, dict] = {}
    for name, thresh in CLINGEN_THRESHOLDS.items():
        if "pp3" in name:
            predicted = scores >= thresh
        else:
            predicted = scores <= thresh
        tp = int(np.sum(predicted & (y_true == 1)))
        fp = int(np.sum(predicted & (y_true == 0)))
        fn = int(np.sum(~predicted & (y_true == 1)))
        tn = int(np.sum(~predicted & (y_true == 0)))
        sens = tp / (tp + fn) if (tp + fn) > 0 else 0
        spec = tn / (tn + fp) if (tn + fp) > 0 else 0
        ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
        threshold_performance[name] = {
            "threshold": thresh,
            "sensitivity": round(sens, 4),
            "specificity": round(spec, 4),
            "ppv": round(ppv, 4),
            "tp": tp, "fp": fp, "fn": fn, "tn": tn,
        }

    return {
        "n_labeled": len(labeled),
        "n_positive": int(y_true.sum()),
        "n_negative": int((1 - y_true).sum()),
        "roc_auc": round(float(roc_auc), 4),
        "average_precision": round(float(avg_precision), 4),
        "optimal_threshold_youden": round(optimal_threshold, 4),
        "optimal_sensitivity": round(float(tpr[optimal_idx]), 4),
        "optimal_specificity": round(float(1 - fpr[optimal_idx]), 4),
        "clingen_threshold_performance": threshold_performance,
        # Store sampled curve points for plotting (every 50th point)
        "roc_curve": {
            "fpr": [round(float(x), 4) for x in fpr[::50]],
            "tpr": [round(float(x), 4) for x in tpr[::50]],
        },
        "pr_curve": {
            "precision": [round(float(x), 4) for x in precision[::50]],
            "recall": [round(float(x), 4) for x in recall[::50]],
        },
    }


def main() -> None:
    variants = load_erepo_with_scores()
    revel_roc = compute_roc(variants, "revel")

    results = {
        "revel": revel_roc,
        "clingen_calibrated_thresholds": CLINGEN_THRESHOLDS,
    }

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    output = RESULTS_DIR / "roc_calibration.json"
    output.write_text(json.dumps(results, indent=2))

    print(f"\nROC calibration saved to {output}")
    if "error" not in revel_roc:
        print(f"  REVEL AUC: {revel_roc['roc_auc']}")
        print(f"  Average Precision: {revel_roc['average_precision']}")
        print(f"  Optimal threshold (Youden): {revel_roc['optimal_threshold_youden']}")
        print(f"  At optimal: sens={revel_roc['optimal_sensitivity']}, spec={revel_roc['optimal_specificity']}")


if __name__ == "__main__":
    main()
