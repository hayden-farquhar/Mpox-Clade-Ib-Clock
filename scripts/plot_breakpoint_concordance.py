#!/usr/bin/env python3
"""Figure 1 — cross-detector breakpoint concordance for the positive control.

Reconstructs the breakpoint-concordance schematic for the inter-clade Ib/IIb
recombinant (OZ375330.1). Three detector rows are aligned over an NC_063383
reference-genome bar:

  * 3SEQ smallest-P triplet  (P = 1.18e-32; breakpoint intervals 93,744-93,894 / 166,694-166,842)
  * 3SEQ alternative triplet (breakpoints 94,134 / 166,694)
  * GARD insurance 2-bp model (converged breakpoints 94,133 / 166,780; cAIC 505,381.68)

The PRIMARY (smallest-P) 3SEQ triplet is the headline concordance: GARD's first
breakpoint (94,133) sits ~240 bp from the primary interval (93,744-93,894) and its
second (166,780) falls inside the primary interval (166,694-166,842) — coincident.
The exact 1-bp match is a property of 3SEQ's ALTERNATIVE triplet (94,134) only and is
shown as a secondary annotation, not the headline, to avoid overstating concordance.
A shaded band marks the intra-clade Ib recombination hotspot reported by Feehley
et al. (~166,442-166,579), adjacent to the second breakpoint.

All values are from RUNLOG.md sec.5.4 (3SEQ) and the insurance-GARD final block.
The legend is placed below the axis so it no longer overlaps the genome bar.

Output: outputs/figures/breakpoint_concordance.png
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

try:
    import sys as _sys
    _sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
    import _integrity_guard  # noqa: F401
except ImportError:
    pass

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "outputs" / "figures" / "breakpoint_concordance.png"

GENOME_LEN = 196_967
BP1 = 94_133
BP2 = 166_780

# 3SEQ smallest-P triplet breakpoint intervals (primary call)
SMALLEST_BP1 = (93_744, 93_894)
SMALLEST_BP2 = (166_694, 166_842)
# 3SEQ alternative triplet breakpoints (single positions)
ALT_BP1 = 94_134
ALT_BP2 = 166_694
# Feehley et al. intra-clade Ib recombination hotspot adjacent to breakpoint 2
FEEHLEY_HOTSPOT = (166_442, 166_579)

C_SMALLEST = "#b3202c"   # dark red
C_ALT = "#f4a6a6"        # light red / pink
C_GARD = "#2a6094"       # blue
C_FEEHLEY = "#8a6d3b"    # muted brown for hotspot band

# row y-positions (top to bottom)
Y_SMALLEST = 3.0
Y_ALT = 2.2
Y_BAR = 1.0

ROW_LABELS = [
    (Y_SMALLEST, "3SEQ (primary, smallest-P triplet)"),
    (Y_ALT, "3SEQ (alternative triplet)"),
    (Y_BAR, "GARD (insurance, 2-bp model)"),
]


def main():
    fig, ax = plt.subplots(figsize=(13.0, 5.0))

    # --- reference genome bar ---
    ax.plot([0, GENOME_LEN], [Y_BAR, Y_BAR], color="#c9c9c9", lw=10,
            solid_capstyle="butt", zorder=1)

    # --- Feehley et al. intra-clade Ib hotspot band (spans all rows) ---
    ax.axvspan(FEEHLEY_HOTSPOT[0], FEEHLEY_HOTSPOT[1], ymin=0.05, ymax=0.78,
               color=C_FEEHLEY, alpha=0.16, zorder=0)

    # --- 3SEQ smallest-P triplet breakpoint INTERVALS (dark red spans) ---
    for lo, hi in (SMALLEST_BP1, SMALLEST_BP2):
        ax.hlines(Y_SMALLEST, lo, hi, color=C_SMALLEST, lw=7,
                  capstyle="butt", zorder=3)
        for x in (lo, hi):
            ax.vlines(x, Y_SMALLEST - 0.16, Y_SMALLEST + 0.16,
                      color=C_SMALLEST, lw=2, zorder=3)
    # --- 3SEQ alternative triplet ticks (pink) ---
    for x in (ALT_BP1, ALT_BP2):
        ax.vlines(x, Y_ALT - 0.18, Y_ALT + 0.18, color=C_ALT, lw=3, zorder=3)

    # --- GARD converged breakpoints (blue triangles on the bar) ---
    for x, lab in ((BP1, f"{BP1:,}"), (BP2, f"{BP2:,}")):
        ax.plot(x, Y_BAR, marker="v", markersize=16, color=C_GARD,
                markeredgecolor="#1d466d", zorder=4)
        ax.text(x, Y_BAR - 0.32, lab, ha="center", va="top", fontsize=11,
                fontweight="bold", color=C_GARD)

    # --- left-hand row labels ---
    for y, lab in ROW_LABELS:
        ax.text(-3_500, y, lab, ha="right", va="center", fontsize=12,
                fontweight="bold")

    # --- centred genome-bar caption ---
    ax.text(GENOME_LEN / 2, Y_BAR - 0.62,
            f"NC_063383 reference genome ({GENOME_LEN:,} bp)",
            ha="center", va="top", fontsize=12, fontweight="bold")

    # --- top concordance annotations with arrows to the GARD triangles ---
    # Headline = PRIMARY (smallest-P) triplet; the 1-bp match is a secondary note.
    ax.annotate(
        "GARD 94,133 ↔ 3SEQ primary 93,744–93,894\n"
        "(within ~240 bp; 1 bp vs alternative triplet 94,134)",
        xy=(BP1, Y_BAR + 0.12), xytext=(48_000, 4.30),
        ha="center", va="bottom", fontsize=11, fontweight="bold", color=C_GARD,
        arrowprops=dict(arrowstyle="->", color=C_GARD, lw=1.4),
    )
    ax.annotate(
        "GARD 166,780 inside 3SEQ primary 166,694–166,842\n"
        "(coincident point estimates)",
        xy=(BP2, Y_BAR + 0.12), xytext=(150_000, 4.30),
        ha="center", va="bottom", fontsize=11, fontweight="bold", color=C_GARD,
        arrowprops=dict(arrowstyle="->", color=C_GARD, lw=1.4),
    )
    # --- Feehley hotspot label ---
    ax.text(sum(FEEHLEY_HOTSPOT) / 2, Y_ALT - 0.55,
            "Feehley et al.\nintra-clade Ib\nhotspot",
            ha="center", va="top", fontsize=8.5, style="italic",
            color=C_FEEHLEY)

    # --- legend in the empty upper-left quadrant ---
    handles = [
        Patch(facecolor=C_SMALLEST, edgecolor="none",
              label="3SEQ primary triplet breakpoint interval (P = 1.18 × 10⁻³²)"),
        Patch(facecolor=C_ALT, edgecolor="none",
              label="3SEQ alternative triplet"),
        Line2D([0], [0], marker="v", color="none", markerfacecolor=C_GARD,
               markeredgecolor=C_GARD, markersize=12,
               label="GARD 2-bp converged breakpoint (cAIC 505,381.68)"),
        Patch(facecolor=C_FEEHLEY, edgecolor="none", alpha=0.16,
              label="Feehley et al. intra-clade Ib hotspot (~166,442–166,579)"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.16),
              ncol=2, fontsize=10.5, framealpha=0.95, borderpad=0.7,
              handletextpad=0.6, columnspacing=1.8, frameon=True)

    # --- title ---
    ax.set_title(
        "Cross-detector breakpoint concordance for the inter-clade Ib/IIb "
        "recombinant (positive control OZ375330.1)",
        fontsize=14, fontweight="bold", pad=16)

    # --- axes cosmetics ---
    ax.set_xlim(-30_000, 207_000)
    ax.set_ylim(0.0, 4.9)
    ax.set_xlabel("Position (NC_063383 reference frame)", fontsize=13,
                  fontweight="bold")
    ax.set_xticks(range(0, 200_001, 20_000))
    ax.set_xticklabels([f"{i} kb" for i in range(0, 201, 20)], fontsize=11)
    ax.set_yticks([])
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)

    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=300, bbox_inches="tight")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
