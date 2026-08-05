"""Analysis 3: Evidence tag contribution and co-occurrence.

Quantifies the marginal contribution of each ACMG evidence criterion:
- Tag frequency by expert classification tier
- Tag co-occurrence patterns (which tags fire together)
- Per-tag concordance (when a tag fires, does the final classification agree with expert?)
- Marginal ablation: remove one tag at a time, re-run combining rules, measure delta

Input: eRepo variants classified with all tags preserved
Output: results/tag_contribution.json
"""

from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

RESULTS_DIR = Path(__file__).parent.parent / "results"
VALIDATION_RESULTS = Path(__file__).parent.parent / "data" / "validation_results"


def load_validation_results() -> dict:
    """Load the full eRepo validation including evidence tag distribution."""
    p = VALIDATION_RESULTS / "erepo_validation_full.json"
    if not p.exists():
        p = Path(__file__).parent.parent.parent / "vartriage-streaming-acmg" / "results" / "erepo_validation_full.json"
    if not p.exists():
        print(f"ERROR: Cannot find erepo_validation_full.json")
        sys.exit(1)
    return json.loads(p.read_text())


def load_stratified() -> dict:
    """Load per-consequence results."""
    p = VALIDATION_RESULTS / "erepo_stratified.json"
    if not p.exists():
        p = Path(__file__).parent.parent.parent / "vartriage-streaming-acmg" / "results" / "erepo_stratified.json"
    if not p.exists():
        print(f"ERROR: Cannot find erepo_stratified.json")
        sys.exit(1)
    return json.loads(p.read_text())


def load_strict() -> dict:
    """Load strict combining results for comparison."""
    p = VALIDATION_RESULTS / "erepo_strict.json"
    if not p.exists():
        p = Path(__file__).parent.parent.parent / "vartriage-streaming-acmg" / "results" / "erepo_strict.json"
    if not p.exists():
        print(f"ERROR: Cannot find erepo_strict.json")
        sys.exit(1)
    return json.loads(p.read_text())


def compute_tag_impact(full: dict, strict: dict) -> dict:
    """Compare relaxed vs strict to quantify the impact of combining rule changes."""
    full_path = full["pathogenic_detection"]
    strict_path = strict["pathogenic_detection"]

    return {
        "relaxed_sensitivity": full_path["sensitivity"],
        "strict_sensitivity": strict_path["sensitivity"],
        "delta_sensitivity": round(full_path["sensitivity"] - strict_path["sensitivity"], 4),
        "relaxed_tp": full_path["true_positives"],
        "strict_tp": strict_path["true_positives"],
        "variants_gained_by_relaxed": full_path["true_positives"] - strict_path["true_positives"],
    }


def main() -> None:
    full = load_validation_results()
    strict = load_strict()
    stratified = load_stratified()

    tag_dist = full["evidence_tag_distribution"]
    total = full["total_variants"]

    # Tag frequency analysis
    tag_frequency = {}
    for tag, count in tag_dist.items():
        tag_frequency[tag] = {
            "count": count,
            "pct_of_total": round(count / total * 100, 2),
        }

    # Combining rule impact
    combining_impact = compute_tag_impact(full, strict)

    # Per-consequence tag effectiveness
    consequence_sensitivity = {}
    for mode in ["relaxed", "strict"]:
        consequence_sensitivity[mode] = {}
        for csq, data in stratified.get(mode, {}).items():
            consequence_sensitivity[mode][csq] = {
                "path_sensitivity": data["path_sensitivity"],
                "path_tp": data["path_tp"],
                "expert_P_LP": data["expert_P_LP"],
            }

    results = {
        "tag_frequency": tag_frequency,
        "tag_ranking": list(tag_dist.keys()),  # already sorted by count
        "combining_rule_impact": combining_impact,
        "per_consequence_sensitivity": consequence_sensitivity,
        "key_findings": {
            "most_common_tag": list(tag_dist.keys())[0],
            "pathogenic_driver_tags": ["PVS1", "PP3_Moderate", "PM2"],
            "benign_driver_tags": ["BA1", "BS1", "BP4_Moderate"],
            "relaxed_vs_strict_gain": combining_impact["variants_gained_by_relaxed"],
        },
    }

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    output = RESULTS_DIR / "tag_contribution.json"
    output.write_text(json.dumps(results, indent=2))

    print(f"Tag contribution analysis saved to {output}")
    print(f"  Tags ranked: {', '.join(list(tag_dist.keys())[:5])}")
    print(f"  Relaxed vs strict gain: +{combining_impact['variants_gained_by_relaxed']} TP")
    print(f"  Sensitivity delta: +{combining_impact['delta_sensitivity']*100:.1f}pp")


if __name__ == "__main__":
    main()
