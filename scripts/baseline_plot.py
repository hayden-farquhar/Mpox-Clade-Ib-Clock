#!/usr/bin/env python3
"""
Phase 2 baseline diagnostic plot: APOBEC3 SNV count vs collection date,
coloured by Nextclade clade. Replicates the O'Toole-et-al. signature plot.

Also writes outputs/tables/baseline.csv with the per-clade headline numbers
needed for the manuscript Results section.
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # non-interactive backend for headless / Snakemake runs
import matplotlib.pyplot as plt
import pandas as pd


def parse_date_loose(s):
    """Parse loose date strings: YYYY, YYYY-MM, YYYY-MM-DD."""
    if not isinstance(s, str) or not s.strip():
        return pd.NaT
    s = s.strip()
    for fmt, take in (("%Y-%m-%d", 10), ("%Y-%m", 7), ("%Y", 4)):
        if len(s) < take:
            continue
        try:
            return pd.Timestamp(datetime.strptime(s[:take], fmt))
        except ValueError:
            continue
    return pd.NaT


def main() -> int:
    ap = argparse.ArgumentParser(description="Phase 2 baseline diagnostic plot.")
    ap.add_argument("--in-counts", required=True, type=Path)
    ap.add_argument("--in-metadata", required=True, type=Path)
    ap.add_argument("--out-png", required=True, type=Path)
    ap.add_argument("--out-baseline", required=True, type=Path)
    args = ap.parse_args()

    counts = pd.read_csv(args.in_counts, sep="\t")
    meta = pd.read_csv(args.in_metadata, sep="\t", dtype=str)

    # Find the collection-date column (NCBI Datasets uses 'Isolate Collection date')
    date_col = next(
        (c for c in ("Isolate Collection date", "isolate-collection-date", "collection_date")
         if c in meta.columns),
        None,
    )
    if date_col is None:
        sys.exit(f"No collection-date column in metadata; cols = {list(meta.columns)}")

    meta["collection_date_parsed"] = meta[date_col].apply(parse_date_loose)

    df = counts.merge(
        meta[["accession", "clade", "collection_date_parsed"]] if "clade" in meta.columns
        else meta[["accession", "collection_date_parsed"]],
        on="accession",
        how="left",
    )

    df = df.dropna(subset=["collection_date_parsed"])
    print(f"[plot] {len(df):,} genomes with parseable collection dates", file=sys.stderr)

    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))
    palette = {"Ia": "#1f77b4", "Ib": "#d62728", "IIa": "#7f7f7f",
               "IIb": "#2ca02c", "Ib/IIb": "#9467bd"}
    if "clade" in df.columns:
        for clade, sub in df.groupby(df["clade"].fillna("unassigned")):
            ax.scatter(
                sub["collection_date_parsed"],
                sub["apobec3_snvs"],
                s=12,
                alpha=0.55,
                label=f"{clade} (n={len(sub):,})",
                color=palette.get(clade, "#888888"),
            )
        ax.legend(loc="upper left", fontsize=8, framealpha=0.95)
    else:
        ax.scatter(df["collection_date_parsed"], df["apobec3_snvs"], s=12, alpha=0.5)

    ax.set_xlabel("Collection date")
    ax.set_ylabel("APOBEC3-context SNVs per genome")
    ax.set_title("APOBEC3-context SNV accumulation by collection date")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    args.out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out_png, dpi=300)
    print(f"[plot] wrote {args.out_png}", file=sys.stderr)

    # Headline baseline table — per-clade and overall
    rows = []
    if "clade" in df.columns:
        for clade, sub in df.groupby(df["clade"].fillna("unassigned")):
            rows.append({
                "clade": clade,
                "n_genomes": len(sub),
                "median_total_snvs": float(sub["total_snvs"].median()),
                "median_apobec3_snvs": float(sub["apobec3_snvs"].median()),
                "median_apobec3_fraction": float(sub["apobec3_fraction"].median()),
                "pooled_apobec3_fraction": float(
                    sub["apobec3_snvs"].sum() / sub["total_snvs"].sum()
                    if sub["total_snvs"].sum() else 0.0
                ),
            })
    rows.append({
        "clade": "ALL",
        "n_genomes": len(df),
        "median_total_snvs": float(df["total_snvs"].median()),
        "median_apobec3_snvs": float(df["apobec3_snvs"].median()),
        "median_apobec3_fraction": float(df["apobec3_fraction"].median()),
        "pooled_apobec3_fraction": float(
            df["apobec3_snvs"].sum() / df["total_snvs"].sum()
            if df["total_snvs"].sum() else 0.0
        ),
    })

    baseline_df = pd.DataFrame(rows)
    baseline_df.to_csv(args.out_baseline, index=False)
    print(f"[plot] wrote {args.out_baseline}", file=sys.stderr)
    print(baseline_df.to_string(index=False), file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
