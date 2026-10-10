"""PM2 strength sensitivity: recompute eRepo classifications with PM2 at Supporting.

Deterministic replay of the committed per-variant evidence tags through the
engine's own point ladder (vartriage/classification/points.py), changing only
PM2 from Moderate (2 pt) to Supporting (1 pt). No engine run, no network: the
baseline evidence_tags are fixed inputs, so the only variable is PM2 weight.

Baseline (PM2=Moderate) classifications are read straight from the file and
checked against a replay to confirm the replay reproduces the engine exactly
before trusting the counterfactual.
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

CLASSIFICATIONS = (
    Path(__file__).parent.parent.parent
    / "vartriage-streaming-acmg"
    / "data"
    / "erepo"
    / "erepo_classifications.json"
)

# Strength -> magnitude, from points.py
STRENGTH_POINTS = {"Very_Strong": 8, "Strong": 4, "Moderate": 2, "Supporting": 1}

# Tag -> strength, from EVIDENCE_STRENGTH_MAP in models/variant.py.
BASELINE_STRENGTH = {
    "PVS1": "Very_Strong",
    "PVS1_STRONG": "Strong",
    "PS1": "Strong",
    "PM1": "Moderate",
    "PM2": "Moderate",
    "PM4": "Moderate",
    "PM5": "Moderate",
    "PP3": "Supporting",
    "PP3_MODERATE": "Moderate",
    "PP3_STRONG": "Strong",
    "PP5": "Supporting",
    "PP5_STRONG": "Strong",
    "BA1": "Standalone",
    "BS1": "Strong",
    "BS2": "Strong",
    "BP4": "Supporting",
    "BP4_MODERATE": "Moderate",
    "BP7": "Supporting",
}

BENIGN_TAGS = {"BA1", "BS1", "BS2", "BP4", "BP4_MODERATE", "BP7"}

PATHOGENIC_THRESHOLD = 10
LIKELY_PATHOGENIC_THRESHOLD = 6
LIKELY_BENIGN_THRESHOLD = -1
BENIGN_THRESHOLD = -7


def score(tags: list[str], strength_map: dict[str, str]) -> int:
    total = 0
    for tag in tags:
        if tag == "BA1":
            continue
        if tag not in strength_map:
            # Unknown tag: skip, but flag so we never silently miscount.
            raise KeyError(f"tag not in strength map: {tag}")
        magnitude = STRENGTH_POINTS[strength_map[tag]]
        if tag in BENIGN_TAGS:
            total -= magnitude
        else:
            total += magnitude
    return total


def classify(tags: list[str], strength_map: dict[str, str]) -> str:
    if "BA1" in tags:
        return "Benign"
    total = score(tags, strength_map)
    if total >= PATHOGENIC_THRESHOLD:
        return "Pathogenic"
    if total >= LIKELY_PATHOGENIC_THRESHOLD:
        return "Likely_Pathogenic"
    if total <= BENIGN_THRESHOLD:
        return "Benign"
    if total <= LIKELY_BENIGN_THRESHOLD:
        return "Likely_Benign"
    return "VUS"


def main() -> None:
    data = json.loads(CLASSIFICATIONS.read_text())
    variants = data["variants"]

    pm2_supporting = dict(BASELINE_STRENGTH)
    pm2_supporting["PM2"] = "Supporting"

    replay_mismatch = []
    carries_pm2 = 0
    moved = []
    base_dist: Counter[str] = Counter()
    new_dist: Counter[str] = Counter()
    unknown_tags: set[str] = set()

    for v in variants:
        tags = v["evidence_tags"]
        stored = v["acmg_classification"]
        try:
            base = classify(tags, BASELINE_STRENGTH)
            new = classify(tags, pm2_supporting)
        except KeyError as exc:
            unknown_tags.add(str(exc).split(": ")[-1])
            continue

        # Replay must reproduce the stored baseline, else our model is wrong.
        if base != stored:
            replay_mismatch.append((v["gene_name"], v["position"], tags, stored, base))

        base_dist[stored] += 1
        new_dist[new] += 1

        if "PM2" in tags:
            carries_pm2 += 1
        if new != base:
            moved.append(
                {
                    "gene": v["gene_name"],
                    "pos": v["position"],
                    "consequence": v["functional_consequence"],
                    "tags": tags,
                    "base_points": score(tags, BASELINE_STRENGTH),
                    "new_points": score(tags, pm2_supporting),
                    "from": base,
                    "to": new,
                    "clinvar": v["clinvar_assertion"],
                }
            )

    print(f"total variants: {len(variants)}")
    print(f"unknown tags encountered: {sorted(unknown_tags) or 'none'}")
    print(f"replay reproduces stored baseline: "
          f"{len(replay_mismatch)} mismatches "
          f"({'EXACT' if not replay_mismatch else 'MODEL DOES NOT MATCH ENGINE'})")
    if replay_mismatch[:5]:
        print("  first mismatches:")
        for m in replay_mismatch[:5]:
            print("   ", m)

    print(f"\nvariants carrying PM2: {carries_pm2}")
    print(f"variants that CHANGE classification under PM2=Supporting: {len(moved)}")

    move_pairs = Counter((m["from"], m["to"]) for m in moved)
    print("\nmovement (from -> to : count):")
    for (frm, to), n in sorted(move_pairs.items(), key=lambda x: -x[1]):
        print(f"  {frm:18s} -> {to:18s} : {n}")

    print("\nbaseline tier distribution (replayed):")
    for t in ["Pathogenic", "Likely_Pathogenic", "VUS", "Likely_Benign", "Benign"]:
        print(f"  {t:18s} {base_dist[t]:6d} -> {new_dist[t]:6d}  (delta {new_dist[t]-base_dist[t]:+d})")

    lp_to_vus = [m for m in moved if m["from"] == "Likely_Pathogenic" and m["to"] == "VUS"]
    print(f"\nLP -> VUS specifically: {len(lp_to_vus)}")

    out = {
        "total_variants": len(variants),
        "replay_exact": not replay_mismatch,
        "replay_mismatches": len(replay_mismatch),
        "pm2_carriers": carries_pm2,
        "changed": len(moved),
        "movement": {f"{k[0]}->{k[1]}": n for k, n in move_pairs.items()},
        "lp_to_vus": len(lp_to_vus),
        "base_dist": dict(base_dist),
        "new_dist": dict(new_dist),
        "moved_detail": moved,
    }
    outpath = Path(__file__).parent.parent / "results" / "pm2_supporting_sensitivity.json"
    outpath.write_text(json.dumps(out, indent=2))
    print(f"\nwrote {outpath}")


if __name__ == "__main__":
    sys.exit(main())
