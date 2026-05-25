#!/usr/bin/env python3
"""
Branch-quantity diagnostic plot: per-tip branch-quantity APOBEC3 SNV count vs
collection date, for the H1ʹ analysis (Amendment 02). Companion to the
tip-vs-reference diagnostic in scripts/baseline_plot.py.

Only clade-Ib tips are present in the input counts file.
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def parse_date_loose(s):
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
    ap = argparse.ArgumentParser(description="Branch-quantity APOBEC3 diagnostic plot.")
    ap.add_argument("--in-counts", required=True, type=Path)
    ap.add_argument("--in-metadata", required=True, type=Path)
    ap.add_argument("--out-png", required=True, type=Path)
    args = ap.parse_args()

    counts = pd.read_csv(args.in_counts, sep="\t")
    meta = pd.read_csv(args.in_metadata, sep="\t", dtype=str)

    date_col = next(
        (c for c in ("Isolate Collection date", "isolate-collection-date", "collection_date")
         if c in meta.columns), None)
    if date_col is None:
        sys.exit(f"No collection-date column in metadata; cols = {list(meta.columns)}")
    meta["collection_date_parsed"] = meta[date_col].apply(parse_date_loose)

    df = counts.merge(meta[["accession", "collection_date_parsed"]],
                      on="accession", how="left")
    df = df.dropna(subset=["collection_date_parsed"])
    print(f"[plot] {len(df):,} Ib tips with parseable collection dates", file=sys.stderr)

    import matplotlib.dates as mdates
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    import numpy as np
    rng = np.random.default_rng(42)
    jitter_apo = rng.uniform(-0.2, 0.2, size=len(df))
    jitter_non = rng.uniform(-0.2, 0.2, size=len(df))
    ax1.scatter(df["collection_date_parsed"], df["apobec3_snvs"] + jitter_apo,
                s=22, alpha=0.55, color="#d62728", edgecolor="white", linewidth=0.4,
                label=f"APOBEC3-context (Ib, n={len(df):,})")
    ax1.scatter(df["collection_date_parsed"], df["non_apobec3_snvs"] + jitter_non,
                s=22, alpha=0.45, color="#1f77b4", edgecolor="white", linewidth=0.4,
                label="Non-APOBEC3 (same tips)")
    ax1.set_xlabel("Collection date")
    ax1.set_ylabel("Branch SNVs per tip (MRCA -> tip)")
    ax1.set_title("A. Branch-quantity SNV accumulation (clade Ib)")
    ax1.legend(loc="upper left", fontsize=9, framealpha=0.95)
    ax1.grid(True, alpha=0.3)
    ax1.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1, 7]))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    plt.setp(ax1.get_xticklabels(), rotation=45, ha="right")

    ax2.scatter(df["collection_date_parsed"], df["apobec3_fraction"],
                s=14, alpha=0.6, color="#d62728")
    ax2.axhline(0.70, color="black", linestyle="--", linewidth=1,
                label="H1prime null p0 = 0.70")
    pooled = df["apobec3_snvs"].sum() / df["total_snvs"].sum() if df["total_snvs"].sum() else 0.0
    ax2.axhline(pooled, color="#d62728", linestyle=":", linewidth=1.5,
                label=f"Pooled observed = {pooled:.3f}")
    ax2.set_xlabel("Collection date")
    ax2.set_ylabel("Per-tip APOBEC3 fraction")
    ax2.set_title("B. Per-tip APOBEC3 fraction vs date")
    ax2.set_ylim(0, 1.05)
    ax2.legend(loc="lower right", fontsize=9, framealpha=0.95)
    ax2.grid(True, alpha=0.3)
    ax2.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1, 7]))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    plt.setp(ax2.get_xticklabels(), rotation=45, ha="right")

    fig.tight_layout()
    args.out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out_png, dpi=300)
    print(f"[plot] wrote {args.out_png}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
