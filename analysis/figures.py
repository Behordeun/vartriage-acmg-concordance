"""Publication figures for the ACMG concordance paper.

Generates 6 figures as PDF (vector) for journal submission:
  Fig 1: 5x5 confusion matrix heatmap
  Fig 2: REVEL ROC curve with ClinGen threshold annotations
  Fig 3: Ablation waterfall (sensitivity change per removed source)
  Fig 4: Evidence tag frequency bar chart
  Fig 5: Consequence-specific pathogenic odds ratios (forest-style)
  Fig 6: Multi-tool sensitivity comparison

Usage:
    uv run python3 analysis/figures.py
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import seaborn as sns

RESULTS_DIR = Path(__file__).parent.parent / "results"
FIGURES_DIR = Path(__file__).parent.parent / "figures"

# Publication style
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.1,
})

# Color palette: colorblind-safe
COLORS = {
    "primary": "#2166AC",
    "secondary": "#B2182B",
    "tertiary": "#4DAF4A",
    "neutral": "#636363",
    "highlight": "#FF7F00",
    "light_blue": "#92C5DE",
    "light_red": "#F4A582",
}


def load_results(filename: str) -> dict:
    path = RESULTS_DIR / filename
    return json.loads(path.read_text())


# Figure 1: Confusion Matrix Heatmap


def figure_confusion_matrix() -> None:
    """5x5 confusion matrix: Expert Panel (rows) vs VarTriage (columns)."""
    data = load_results("concordance.json")
    cm = data["confusion_matrix"]

    tiers = ["Pathogenic", "Likely_Pathogenic", "VUS", "Likely_Benign", "Benign"]
    labels = ["P", "LP", "VUS", "LB", "B"]

    matrix = np.zeros((5, 5), dtype=int)
    for i, expert_tier in enumerate(tiers):
        for j, vt_tier in enumerate(tiers):
            matrix[i, j] = cm[expert_tier].get(vt_tier, 0)

    # Normalize by row (expert assertion) to show % of each expert class
    row_sums = matrix.sum(axis=1, keepdims=True)
    matrix_pct = np.where(row_sums > 0, matrix / row_sums * 100, 0)

    fig, ax = plt.subplots(figsize=(4.5, 4.0))

    # Use diverging colormap centered on the diagonal
    sns.heatmap(
        matrix_pct,
        annot=True,
        fmt=".1f",
        cmap="Blues",
        xticklabels=labels,
        yticklabels=labels,
        ax=ax,
        cbar_kws={"label": "% of expert class"},
        linewidths=0.5,
        linecolor="white",
        vmin=0,
        vmax=100,
    )

    # Overlay raw counts in smaller font
    for i in range(5):
        for j in range(5):
            if matrix[i, j] > 0:
                ax.text(
                    j + 0.5, i + 0.78, f"n={matrix[i, j]}",
                    ha="center", va="center", fontsize=6,
                    color="gray" if matrix_pct[i, j] < 50 else "lightgray",
                )

    ax.set_xlabel("VarTriage classification")
    ax.set_ylabel("Expert Panel assertion")
    ax.set_title("Classification concordance (row-normalized %)")

    fig.savefig(FIGURES_DIR / "fig1_confusion_matrix.pdf")
    plt.close(fig)
    print("  Fig 1: Confusion matrix")


# Figure 2: ROC Curve


def figure_roc_curve() -> None:
    """REVEL ROC curve with ClinGen calibrated threshold markers."""
    data = load_results("roc_calibration.json")
    revel = data["revel"]
    roc = revel["roc_curve"]
    thresholds = revel["clingen_threshold_performance"]

    fpr = np.array(roc["fpr"])
    tpr = np.array(roc["tpr"])
    auc = revel["roc_auc"]

    fig, ax = plt.subplots(figsize=(4.5, 4.0))

    # ROC curve
    ax.plot(fpr, tpr, color=COLORS["primary"], linewidth=2,
            label=f"REVEL (AUC = {auc:.3f})")

    # Diagonal reference
    ax.plot([0, 1], [0, 1], color=COLORS["neutral"], linewidth=0.8,
            linestyle="--", label="Random classifier")

    # Mark ClinGen thresholds
    pp3_mod = thresholds["revel_pp3_moderate"]
    pp3_sup = thresholds["revel_pp3_supporting"]
    optimal = revel["optimal_threshold_youden"]

    # Optimal threshold (Youden)
    opt_sens = revel["optimal_sensitivity"]
    opt_spec = revel["optimal_specificity"]
    ax.plot(1 - opt_spec, opt_sens, "o", color=COLORS["highlight"],
            markersize=8, zorder=5, label=f"Optimal (t={optimal:.3f})")

    # PP3 moderate (0.773)
    ax.plot(1 - pp3_mod["specificity"], pp3_mod["sensitivity"], "^",
            color=COLORS["secondary"], markersize=7, zorder=5,
            label="PP3 Moderate (t=0.773)")

    # PP3 supporting (0.644)
    ax.plot(1 - pp3_sup["specificity"], pp3_sup["sensitivity"], "s",
            color=COLORS["tertiary"], markersize=6, zorder=5,
            label="PP3 Supporting (t=0.644)")

    ax.set_xlabel("False Positive Rate (1 - Specificity)")
    ax.set_ylabel("True Positive Rate (Sensitivity)")
    ax.set_title("REVEL pathogenicity prediction performance")
    ax.legend(loc="lower right", framealpha=0.9)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect("equal")

    fig.savefig(FIGURES_DIR / "fig2_roc_curve.pdf")
    plt.close(fig)
    print("  Fig 2: ROC curve")


# Figure 3: Ablation Waterfall


def figure_ablation() -> None:
    """Horizontal bar chart showing sensitivity change when each source is removed."""
    data = load_results("ablation.json")

    # Exclude baseline from the chart
    configs = [(k, v) for k, v in data.items() if k != "baseline"]
    configs.sort(key=lambda x: x[1].get("delta_path_sensitivity", 0))

    names = []
    deltas = []
    for name, row in configs:
        delta = row.get("delta_path_sensitivity", 0)
        names.append(row.get("description", name))
        deltas.append(delta * 100)  # Convert to percentage points

    fig, ax = plt.subplots(figsize=(5.5, 3.0))

    colors = [COLORS["secondary"] if d < 0 else COLORS["tertiary"] for d in deltas]
    bars = ax.barh(range(len(names)), deltas, color=colors, edgecolor="white",
                   linewidth=0.5, height=0.6)

    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=8)
    ax.set_xlabel("Change in pathogenic sensitivity (pp)")
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_title("Impact of data source removal on classification sensitivity")

    # Value labels
    for i, (bar, delta) in enumerate(zip(bars, deltas)):
        offset = 1.0 if delta < 0 else -1.0
        ha = "left" if delta < 0 else "right"
        ax.text(delta + offset, i, f"{delta:+.1f}pp", va="center", ha=ha,
                fontsize=7, color=COLORS["neutral"])

    # Baseline annotation
    baseline_sens = data["baseline"]["path_sensitivity"] * 100
    ax.text(0.98, 0.02, f"Baseline: {baseline_sens:.1f}%",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=7, color=COLORS["neutral"],
            bbox={"boxstyle": "round,pad=0.3", "facecolor": "lightyellow",
                  "edgecolor": "gray", "linewidth": 0.5})

    ax.set_xlim(min(deltas) - 5, max(deltas) + 8)

    fig.savefig(FIGURES_DIR / "fig3_ablation.pdf")
    plt.close(fig)
    print("  Fig 3: Ablation waterfall")


# Figure 4: Evidence Tag Frequency


def figure_tag_frequency() -> None:
    """Bar chart of evidence tag assignment frequency across all variants."""
    data = load_results("tag_contribution.json")
    tags = data["tag_frequency"]

    # Sort by count descending
    sorted_tags = sorted(tags.items(), key=lambda x: x[1]["count"], reverse=True)

    names = [t[0] for t in sorted_tags]
    counts = [t[1]["count"] for t in sorted_tags]

    # Color by evidence direction via tag prefix, so strength-tier and
    # constraint tags (PP3_Strong, PVS1_Strong, PM1, PS1, PM5) are classified
    # correctly rather than falling through to neutral.
    def _is_benign(tag: str) -> bool:
        return tag.startswith(("BA", "BS", "BP"))

    colors = [COLORS["primary"] if _is_benign(name) else COLORS["secondary"] for name in names]

    fig, ax = plt.subplots(figsize=(5.0, 3.5))

    bars = ax.bar(range(len(names)), counts, color=colors, edgecolor="white",
                  linewidth=0.5, width=0.7)

    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Variants with tag assigned")
    ax.set_title("ACMG evidence tag frequency (n=15,334)")

    # Legend
    path_patch = mpatches.Patch(color=COLORS["secondary"], label="Pathogenic evidence")
    ben_patch = mpatches.Patch(color=COLORS["primary"], label="Benign evidence")
    ax.legend(handles=[path_patch, ben_patch], loc="upper right")

    # Percentage labels on top of bars
    total = 15334
    for i, (bar, count) in enumerate(zip(bars, counts)):
        pct = count / total * 100
        if pct > 3:
            ax.text(i, count + 100, f"{pct:.0f}%", ha="center", va="bottom",
                    fontsize=6, color=COLORS["neutral"])

    fig.savefig(FIGURES_DIR / "fig4_tag_frequency.pdf")
    plt.close(fig)
    print("  Fig 4: Tag frequency")


# Figure 5: Consequence Odds Ratios


def figure_consequence_odds() -> None:
    """Forest-style plot of odds ratios for pathogenicity by consequence type."""
    data = load_results("consequence_enrichment.json")
    odds = data["odds_ratios_for_pathogenicity"]

    # Sort by OR descending, handle "inf"
    items = []
    for csq, vals in odds.items():
        or_val = vals["odds_ratio"]
        if or_val == "inf":
            or_val = 200.0  # Cap for display
        items.append((csq, float(or_val), vals["expert_P_LP"], vals["expert_B_LB"]))

    items.sort(key=lambda x: x[1], reverse=True)

    names = [x[0] for x in items]
    ors = [x[1] for x in items]
    _ = [x[2] for x in items]
    _ = [x[3] for x in items]

    fig, ax = plt.subplots(figsize=(5.0, 3.0))

    y_pos = range(len(names))
    colors = [COLORS["secondary"] if or_val > 1 else COLORS["primary"] for or_val in ors]

    ax.barh(y_pos, ors, color=colors, edgecolor="white", linewidth=0.5, height=0.6)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(names, fontsize=8)
    ax.set_xlabel("Odds ratio for pathogenicity")
    ax.set_xscale("log")
    ax.axvline(1.0, color="black", linewidth=0.8, linestyle="--")
    ax.set_title("Pathogenic enrichment by functional consequence")

    # Annotations: sample sizes
    for i, (name, or_val, plp, blb) in enumerate(items):
        label = f"P/LP={plp}, B/LB={blb}"
        ax.text(max(ors) * 1.1, i, label, va="center", fontsize=6,
                color=COLORS["neutral"])

    ax.set_xlim(0.005, max(ors) * 30)

    fig.savefig(FIGURES_DIR / "fig5_consequence_odds.pdf")
    plt.close(fig)
    print("  Fig 5: Consequence odds ratios")


# Figure 6: Tool Comparison


def figure_tool_comparison() -> None:
    """Grouped bar chart comparing pathogenic/benign sensitivity across tools."""
    data = load_results("tool_comparison.json")
    tools = data["tools"]

    tool_names = list(tools.keys())
    path_sens = [tools[t].get("path_sensitivity", 0) or 0 for t in tool_names]
    ben_sens = [tools[t].get("ben_sensitivity", 0) or 0 for t in tool_names]
    criteria = [tools[t].get("criteria_count", 0) for t in tool_names]

    x = np.arange(len(tool_names))
    width = 0.35

    fig, ax = plt.subplots(figsize=(5.0, 3.5))

    bars1 = ax.bar(x - width / 2, [s * 100 for s in path_sens], width,
                   color=COLORS["secondary"], edgecolor="white", linewidth=0.5,
                   label="Pathogenic sensitivity")
    bars2 = ax.bar(x + width / 2, [s * 100 for s in ben_sens], width,
                   color=COLORS["primary"], edgecolor="white", linewidth=0.5,
                   label="Benign sensitivity")

    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"{name}\n({c} criteria)" for name, c in zip(tool_names, criteria)],
        fontsize=8,
    )
    ax.set_ylabel("Sensitivity (%)")
    ax.set_title("Multi-tool classification performance on ClinGen eRepo")
    ax.legend(loc="upper right")
    ax.set_ylim(0, 100)

    # Value labels
    for bar in bars1:
        h = bar.get_height()
        if h > 0:
            ax.text(bar.get_x() + bar.get_width() / 2, h + 1, f"{h:.1f}%",
                    ha="center", va="bottom", fontsize=7)
    for bar in bars2:
        h = bar.get_height()
        if h > 0:
            ax.text(bar.get_x() + bar.get_width() / 2, h + 1, f"{h:.1f}%",
                    ha="center", va="bottom", fontsize=7)

    fig.savefig(FIGURES_DIR / "fig6_tool_comparison.pdf")
    plt.close(fig)
    print("  Fig 6: Tool comparison")


# Main


def main() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    print("Generating publication figures...")
    figure_confusion_matrix()
    figure_roc_curve()
    figure_ablation()
    figure_tag_frequency()
    figure_consequence_odds()
    figure_tool_comparison()

    print(f"\nAll figures saved to {FIGURES_DIR}/")
    for f in sorted(FIGURES_DIR.glob("*.pdf")):
        print(f"  {f.name}")


if __name__ == "__main__":
    main()
