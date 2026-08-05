"""Analysis 4: Missing data ablation study.

Measures how classification accuracy degrades when data sources are removed.
Re-runs the eRepo validation with each source disabled:
- Baseline: all sources (gnomAD + VEP + REVEL)
- No REVEL: PM2 + PVS1 only (no PP3)
- No gnomAD: no PM2, no BA1/BS1 (frequency-based criteria disabled)
- No VEP: all variants treated as missense (consequence unknown)
- Minimal: gene annotation only (no scores, no frequencies)

Each configuration is run via validate_erepo.py with appropriate flags.
This script orchestrates the runs and compiles the ablation table.

Input: eRepo TSV + reference files + cached API data
Output: results/ablation.json
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

RESULTS_DIR = Path(__file__).parent.parent / "results"
VARTRIAGE_DIR = Path(__file__).parent.parent.parent / "vartriage-streaming-acmg"
VALIDATE_SCRIPT = VARTRIAGE_DIR / "analysis" / "validate_erepo.py"
EREPO = VARTRIAGE_DIR / "data" / "erepo" / "clingen_erepo.tsv"
REVEL = VARTRIAGE_DIR / "data" / "references" / "revel_grch38.tsv"
CACHE_DB = VARTRIAGE_DIR / "data" / "erepo" / ".api_cache.db"


def run_validation(config_name: str, extra_args: list[str], output_file: Path) -> dict | None:
    """Run validate_erepo.py with given arguments and return parsed results."""
    if output_file.exists():
        print(f"  [{config_name}] Using cached result: {output_file}")
        return json.loads(output_file.read_text())

    cmd = [
        sys.executable, str(VALIDATE_SCRIPT),
        "--erepo", str(EREPO),
        "--revel", str(REVEL),
        "--cache-db", str(CACHE_DB),
        "--output", str(output_file),
        "--relaxed-combining",
    ] + extra_args

    print(f"  [{config_name}] Running: {' '.join(cmd[-4:])}")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)

    if result.returncode != 0:
        print(f"  [{config_name}] FAILED: {result.stderr[:200]}")
        return None

    if output_file.exists():
        return json.loads(output_file.read_text())
    return None


def _extract_row(config: dict, data: dict) -> dict:
    """Extract metrics from a validation result into a flat row dict."""
    path_det = data.get("pathogenic_detection", {})
    ben_det = data.get("benign_detection", {})
    return {
        "description": config["description"],
        "path_sensitivity": path_det.get("sensitivity"),
        "path_ppv": path_det.get("ppv"),
        "path_tp": path_det.get("true_positives"),
        "ben_sensitivity": ben_det.get("sensitivity"),
        "ben_ppv": ben_det.get("ppv"),
        "binary_concordance": data.get("overall", {}).get("binary_concordance"),
    }


def _compute_deltas(results_table: dict[str, dict]) -> None:
    """Add delta_path_sensitivity relative to baseline, in-place."""
    baseline = results_table.get("baseline", {})
    base_sens = baseline.get("path_sensitivity")
    if not base_sens:
        return
    for name, row in results_table.items():
        if name == "baseline" or row.get("path_sensitivity") is None:
            continue
        row["delta_path_sensitivity"] = round(
            row["path_sensitivity"] - base_sens, 4
        )


def _print_table(results_table: dict[str, dict]) -> None:
    """Print a summary table of ablation results."""
    print(f"\n{'Config':<20} {'Path Sens':>10} {'Path PPV':>9} {'Delta':>8}")
    print("-" * 50)
    for name, row in results_table.items():
        if "error" in row:
            print(f"{name:<20} {'FAILED':>10}")
            continue
        sens = f"{row['path_sensitivity']*100:.1f}%" if row.get("path_sensitivity") else "N/A"
        ppv = f"{row['path_ppv']*100:.1f}%" if row.get("path_ppv") else "N/A"
        delta = (
            f"{row['delta_path_sensitivity']*100:+.1f}pp"
            if "delta_path_sensitivity" in row
            else "--"
        )
        print(f"{name:<20} {sens:>10} {ppv:>9} {delta:>8}")


def main() -> None:
    if not VALIDATE_SCRIPT.exists():
        print(f"ERROR: validate_erepo.py not found at {VALIDATE_SCRIPT}")
        print("Ensure the vartriage-streaming-acmg repo is adjacent to this repo.")
        sys.exit(1)

    if not EREPO.exists():
        print(f"ERROR: eRepo data not found at {EREPO}")
        sys.exit(1)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    ablation_dir = RESULTS_DIR / "ablation"
    ablation_dir.mkdir(exist_ok=True)

    configs = {
        "baseline": {
            "description": "All sources (gnomAD + VEP + REVEL)",
            "args": [],
        },
        "no_revel": {
            "description": "No REVEL scores (PP3/BP4 disabled for missense)",
            "args": ["--revel", "/dev/null"],  # empty REVEL file
        },
        "no_gnomad": {
            "description": "No gnomAD (PM2/BA1/BS1 disabled)",
            "args": ["--skip-gnomad"],
        },
        "no_vep": {
            "description": "No VEP (all variants treated as missense)",
            "args": ["--skip-vep"],
        },
        "no_gnomad_no_vep": {
            "description": "Minimal: REVEL scores only",
            "args": ["--skip-gnomad", "--skip-vep"],
        },
    }

    results_table: dict[str, dict] = {}

    print("Running ablation experiments...")
    for name, config in configs.items():
        output = ablation_dir / f"ablation_{name}.json"
        data = run_validation(name, config["args"], output)
        if data is None:
            results_table[name] = {"error": "run failed"}
        else:
            results_table[name] = _extract_row(config, data)

    _compute_deltas(results_table)

    output = RESULTS_DIR / "ablation.json"
    output.write_text(json.dumps(results_table, indent=2))

    print(f"\nAblation results saved to {output}")
    _print_table(results_table)


if __name__ == "__main__":
    main()
