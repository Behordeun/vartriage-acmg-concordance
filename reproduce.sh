#!/usr/bin/env bash
# reproduce.sh - Run all analyses for the ACMG concordance paper.
#
# Usage:
#   ./reproduce.sh           # Full analysis pipeline
#   ./reproduce.sh --quick   # Figures only (from existing results)

set -euo pipefail
cd "$(dirname "$0")"

# ============================================================================
# Quick mode: regenerate figures only
# ============================================================================

if [[ "${1:-}" == "--quick" ]]; then
    echo "=== Quick mode: regenerating figures ==="
    uv run --with matplotlib --with seaborn --with numpy --with scikit-learn python3 analysis/figures.py
    echo "Done. Figures in figures/"
    exit 0
fi

# ============================================================================
# Preflight: check VarTriage validation results are accessible
# ============================================================================

echo "=== Preflight checks ==="

VARTRIAGE_DIR="../vartriage-streaming-acmg"

if [[ ! -d "$VARTRIAGE_DIR/results" ]]; then
    echo "ERROR: VarTriage validation results not found at $VARTRIAGE_DIR/results"
    echo "Ensure vartriage-streaming-acmg repo is adjacent to this repo."
    exit 1
fi

# Symlink validation results for easy access
mkdir -p data/validation_results
for f in erepo_validation_full.json erepo_strict.json erepo_relaxed.json erepo_stratified.json; do
    if [[ ! -e "data/validation_results/$f" ]]; then
        ln -sf "../../$VARTRIAGE_DIR/results/$f" "data/validation_results/$f"
    fi
done

# Symlink data if not already present
if [[ ! -e "data/erepo" ]]; then
    ln -sf "../$VARTRIAGE_DIR/data/erepo" data/erepo
fi
if [[ ! -e "data/references" ]]; then
    ln -sf "../$VARTRIAGE_DIR/data/references" data/references
fi

echo "  Validation results: linked"
echo "  Data sources: linked"
echo ""

# ============================================================================
# Analysis 1: Concordance metrics
# ============================================================================

echo "=== Analysis 1/6: Concordance metrics ==="
uv run python3 analysis/01_concordance.py
echo ""

# ============================================================================
# Analysis 2: ROC curves and threshold calibration
# ============================================================================

echo "=== Analysis 2/6: ROC calibration ==="
uv run python3 analysis/02_roc_calibration.py
echo ""

# ============================================================================
# Analysis 3: Evidence tag contribution
# ============================================================================

echo "=== Analysis 3/6: Tag contribution ==="
uv run python3 analysis/03_tag_contribution.py
echo ""

# ============================================================================
# Analysis 4: Missing data ablation
# ============================================================================

echo "=== Analysis 4/6: Missing data ablation ==="
uv run python3 analysis/04_missing_data_ablation.py
echo ""

# ============================================================================
# Analysis 5: Consequence enrichment
# ============================================================================

echo "=== Analysis 5/6: Consequence enrichment ==="
uv run python3 analysis/05_consequence_enrichment.py
echo ""

# ============================================================================
# Analysis 6: Tool comparison
# ============================================================================

echo "=== Analysis 6/6: Tool comparison ==="
uv run python3 analysis/06_tool_comparison.py
echo ""

# ============================================================================
# Generate figures
# ============================================================================

echo "=== Generating figures ==="
uv run --with matplotlib --with seaborn --with numpy --with scikit-learn python3 analysis/figures.py
echo ""

# ============================================================================
# Summary
# ============================================================================

echo "============================================================"
echo "  ALL ANALYSES COMPLETE"
echo "============================================================"
echo ""
echo "  Results:"
for f in results/*.json; do
    echo "    $f"
done
echo ""
echo "  Figures:"
for f in figures/*.pdf; do
    echo "    $f"
done
echo ""
