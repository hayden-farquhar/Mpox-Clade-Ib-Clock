#!/usr/bin/env python3
"""Recombination over-call figure for the restructured manuscript.

Two panels showing why a single-detector (3SEQ) scan over-calls on partial
genomes and how a population-level second detector (PhiPack) separates the one
genuine recombinant from the false positives:

  A. Genome length vs PhiPack Phi-permutation P for every 3SEQ candidate,
     coloured by partial-genome status, positive control highlighted.
  B. Genome-length distribution of the candidates relative to the reference.

Input: data/processed/freeze_20260522/recomb_3seq/candidate_annotation_table.tsv
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

try:
    import sys as _sys
    _sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
    import _integrity_guard  # noqa: F401
except ImportError:
    pass

ROOT = Path(__file__).resolve().parent.parent
FREEZE = ROOT / "data" / "processed" / "freeze_20260522"
REF_LEN = 196967
ALPHA = 0.05


def resolve_partial(df: pd.DataFrame) -> pd.Series:
    """Completeness from the NCBI structured COMPLETE/PARTIAL field.

    The headline completeness counts standardise on the NCBI ``Completeness``
    field (the same field used in Table S6 and committed to in Methods). For the
    three clade-IIb outgroup genomes that are absent from the clade-I metadata
    freeze, fall back to the description-derived ``is_partial_genome`` flag.
    """
    meta = pd.read_csv(FREEZE / "metadata.tsv", sep="\t")[["accession", "Completeness"]]
    merged = df.merge(meta, on="accession", how="left")
    is_partial = merged["Completeness"].eq("PARTIAL")
    missing = merged["Completeness"].isna()
    is_partial = is_partial.where(~missing, merged["is_partial_genome"] == True)  # noqa: E712
    is_partial.index = df.index
    return is_partial


def main():
    df = pd.read_csv(FREEZE / "recomb_3seq" / "candidate_annotation_table.tsv", sep="\t")
    df["Length"] = pd.to_numeric(df["Length"], errors="coerce")
    df["phi_perm_p"] = pd.to_numeric(df["phi_perm_p"], errors="coerce")
    df = df.dropna(subset=["Length", "phi_perm_p"]).copy()
    df["completeness_partial"] = resolve_partial(df)

    pc = df[df["is_positive_control"] == True]  # noqa: E712
    rest = df[df["is_positive_control"] != True]  # noqa: E712
    partial = rest[rest["completeness_partial"] == True]  # noqa: E712
    complete = rest[rest["completeness_partial"] != True]  # noqa: E712

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.6))

    # Panel A: length vs Phi P
    axA.scatter(complete["Length"] / 1000, complete["phi_perm_p"], s=26,
                c="#3b6ea5", alpha=0.65, edgecolors="white", linewidths=0.4,
                label=f"Complete, unconfirmed (n={len(complete)})")
    axA.scatter(partial["Length"] / 1000, partial["phi_perm_p"], s=26,
                c="#d98a2b", alpha=0.7, edgecolors="white", linewidths=0.4,
                label=f"Partial genome, unconfirmed (n={len(partial)})")
    axA.scatter(pc["Length"] / 1000, pc["phi_perm_p"], s=150, marker="*",
                c="#b3202c", edgecolors="black", linewidths=0.6, zorder=5,
                label="Confirmed recombinant (positive control)")
    axA.axhline(ALPHA, color="grey", ls="--", lw=1)
    axA.text(190.3, ALPHA + 0.015, "Phi P = 0.05", fontsize=8, color="grey")
    axA.set_xlabel("Genome length (kb)")
    axA.set_ylabel("PhiPack Phi-permutation P")
    axA.set_ylim(-0.03, 1.0)
    axA.legend(fontsize=7.5, loc="center left", framealpha=0.9)
    axA.set_title("A. Single-detector candidates vs second-detector confirmation", fontsize=10)

    # Panel B: length distribution
    axB.hist(rest["Length"] / 1000, bins=24, color="#9aa7b5", edgecolor="white")
    axB.axvline(REF_LEN / 1000, color="#b3202c", ls="--", lw=1.2)
    axB.text(REF_LEN / 1000 - 0.4, axB.get_ylim()[1] * 0.92, "reference\n197 kb",
             fontsize=8, ha="right", color="#b3202c")
    n_partial = int((rest["completeness_partial"] == True).sum())  # noqa: E712
    axB.set_xlabel("Genome length (kb)")
    axB.set_ylabel("3SEQ candidates")
    axB.set_title(f"B. {n_partial}/{len(rest)} unconfirmed candidates are partial genomes",
                  fontsize=10)

    fig.tight_layout()
    out = ROOT / "outputs" / "figures" / "overcall_partial_genomes.png"
    fig.savefig(out, dpi=300)
    print(f"-> {out}  (complete={len(complete)}, partial={len(partial)}, "
          f"pos_control_phi_p={pc['phi_perm_p'].iloc[0] if len(pc) else 'NA'})")


if __name__ == "__main__":
    main()
