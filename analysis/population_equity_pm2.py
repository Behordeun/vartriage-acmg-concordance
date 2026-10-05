"""Population-equity PM2 analysis (V4) — reference implementation.

Reconstructed from manuscript/population-equity/latex/sections/methods.tex and run
against the pinned vartriage v1.0.0 (19d658c) eRepo gnomAD r4 GraphQL cache
(data/erepo/.api_cache.db, source='gnomad'). The original generator was never
committed; this script implements the published method exactly and is labelled as
a reference implementation, not the original pipeline.

Method (per methods.tex):
  - Retain variants with a non-missing global AF AND >=1 non-missing per-ancestry AF.
  - Per-ancestry rarity threshold for PM2: AF <= 1e-4 applied independently to each
    AVAILABLE ancestry; per-population PM2 fires only when no reported ancestry exceeds it.
  - Global-only PM2: global AF <= 1e-4.
  - false-rare: global PM2-eligible AND >=1 reported ancestry AF > theta (common-in-pop),
    theta in {0.001 (primary), 0.005, 0.01}.
  - Ancestry enrichment at theta=1e-3: (ancestry share of false-rares) /
    (ancestry share of non-zero-AF reports across the retained set). Controls for
    gnomAD per-ancestry data volume.
  - absent_globally: global AF == 0.0.
"""
from __future__ import annotations

import json
import sqlite3
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
CACHE = REPO / "data" / "erepo" / ".api_cache.db"
OUT = REPO / "results" / "population_equity_pm2_v1.0.0.json"

POPS = ["afr", "amr", "asj", "eas", "fin", "nfe", "sas"]
PM2_RARE = 1e-4
COMMON_THRESHOLDS = [0.001, 0.005, 0.01]
PRIMARY_COMMON = 0.001


def load_variants() -> list[dict]:
    con = sqlite3.connect(f"file:{CACHE}?immutable=1", uri=True)
    rows = con.execute(
        "SELECT response_json FROM api_cache WHERE source='gnomad'"
    ).fetchall()
    con.close()
    out = []
    for (rj,) in rows:
        d = json.loads(rj)
        gaf = d.get("global_af")
        # per-ancestry AFs that are actually reported (non-null)
        pops = {p: d[p] for p in POPS if d.get(p) is not None}
        if gaf is None or not pops:
            continue  # need global AND >=1 per-ancestry
        out.append({"global_af": gaf, "pops": pops})
    return out


def main() -> None:
    variants = load_variants()
    n = len(variants)

    absent_globally = sum(1 for v in variants if v["global_af"] == 0.0)

    # PM2 eligibility
    global_pm2 = [v for v in variants if v["global_af"] <= PM2_RARE]
    per_pop_pm2 = [
        v for v in variants if all(af <= PM2_RARE for af in v["pops"].values())
    ]

    # false-rare at each theta (global PM2-eligible but common in >=1 ancestry)
    def false_rare_set(theta: float) -> list[dict]:
        return [
            v for v in global_pm2
            if any(af > theta for af in v["pops"].values())
        ]

    fr_primary = false_rare_set(PRIMARY_COMMON)
    # ancestry attribution for the primary threshold: count a false-rare under
    # EACH ancestry in which it is common
    by_ancestry: dict[str, int] = {p: 0 for p in POPS}
    for v in fr_primary:
        for p, af in v["pops"].items():
            if af > PRIMARY_COMMON:
                by_ancestry[p] += 1
    by_ancestry = {p: c for p, c in sorted(by_ancestry.items(), key=lambda x: -x[1]) if c}

    # false-common: global NOT rare but rare (<=PM2_RARE) in >=1 reported ancestry
    false_common = [
        v for v in variants
        if v["global_af"] > PM2_RARE
        and any(af <= PM2_RARE for af in v["pops"].values())
    ]

    # ancestry enrichment at primary threshold: (share of false-rare) /
    # (share of non-zero-AF reports across the retained set)
    nonzero_reports: dict[str, int] = {p: 0 for p in POPS}
    for v in variants:
        for p, af in v["pops"].items():
            if af > 0.0:
                nonzero_reports[p] += 1
    total_nonzero = sum(nonzero_reports.values())
    total_fr_attrib = sum(by_ancestry.values())
    enrichment = {}
    for p in POPS:
        fr_share = (by_ancestry.get(p, 0) / total_fr_attrib) if total_fr_attrib else 0.0
        rep_share = (nonzero_reports[p] / total_nonzero) if total_nonzero else 0.0
        enrichment[p] = round(fr_share / rep_share, 2) if rep_share > 0 else 0.0
    enrichment = dict(sorted(enrichment.items(), key=lambda x: -x[1]))

    threshold_sensitivity = {}
    for theta in COMMON_THRESHOLDS:
        fr = false_rare_set(theta)
        key = f"common_{theta:g}".replace("0.001", "1e-3").replace("0.005", "5e-3").replace("0.01", "1e-2")
        threshold_sensitivity[key] = {
            "false_rare": len(fr),
            "pct": round(100 * len(fr) / n, 2) if n else 0.0,
        }

    result = {
        "analysis": "per-population vs global-AF PM2 assignment (population equity)",
        "reference_implementation": True,
        "reference_note": (
            "Reconstructed from methods.tex and run against the pinned vartriage "
            "v1.0.0 (19d658c) eRepo gnomAD r4 GraphQL cache. The original generator "
            "was not committed; this script implements the published method exactly."
        ),
        "vartriage_pin": "19d658c",
        "vartriage_version": "1.0.0",
        "data_source": "eRepo gnomAD GraphQL API cache (.api_cache.db), gnomAD r4, genome-wide",
        "genome_build": "GRCh38",
        "n_variants_with_global_and_per_population_af": n,
        "absent_globally_af_zero": absent_globally,
        "pm2_thresholds": {
            "pm2_supporting_absent_or_rare_af_leq": PM2_RARE,
            "common_in_population_af_gt": PRIMARY_COMMON,
        },
        "global_only_pm2_eligible": {
            "count": len(global_pm2),
            "pct": round(100 * len(global_pm2) / n, 1) if n else 0.0,
        },
        "per_population_pm2_eligible_rare_in_all_reported_pops": {
            "count": len(per_pop_pm2),
            "pct": round(100 * len(per_pop_pm2) / n, 1) if n else 0.0,
        },
        "false_rare_global_rare_but_common_in_ge1_ancestry": {
            "count": len(fr_primary),
            "pct": round(100 * len(fr_primary) / n, 2) if n else 0.0,
            "by_ancestry": by_ancestry,
        },
        "false_common_global_notrare_but_rare_in_ge1_ancestry": {
            "count": len(false_common),
            "pct": round(100 * len(false_common) / n, 2) if n else 0.0,
        },
        "enrichment_at_common_1e-3": {
            "definition": "false-rare share in ancestry / that ancestry's share of nonzero-AF reports (controls for data volume)",
            **enrichment,
        },
        "threshold_sensitivity": threshold_sensitivity,
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2))
    print(f"WROTE {OUT}")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
