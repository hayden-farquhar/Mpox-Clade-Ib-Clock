#!/usr/bin/env python3
"""
H1 test (pre-registration §1.4 H1, §5.2): one-sided exact binomial test of the
clade-wide APOBEC3 SNV fraction against the null p0 = 0.70 at alpha = 0.05.
H1 is supported if the lower bound of the 95% Wilson interval exceeds 0.70.

Failure mode (§5.7 stopping table):
- Observed fraction <0.50  → emergency halt; downstream analysis stops.
- Observed fraction in [0.50, 0.70 lower bound]  → H1 not supported but report
  transparently and proceed.

Restricts to clade Ib by default (the population the H1 test is about);
the metadata must carry a `clade` column (e.g. from Nextclade).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import binomtest


def wilson_interval(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion."""
    if n == 0:
        return (0.0, 0.0)
    from scipy.stats import norm
    z = norm.ppf(1 - alpha / 2)
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5) / denom
    return max(0.0, centre - half), min(1.0, centre + half)


def main() -> int:
    ap = argparse.ArgumentParser(description="H1 test: clade-Ib APOBEC3 fraction vs p0 = 0.70.")
    ap.add_argument("--in-counts", required=True, type=Path,
                    help="TSV produced by count_apobec3.py")
    ap.add_argument("--in-metadata", required=True, type=Path,
                    help="Per-genome metadata (joined on accession; must have 'clade' column)")
    ap.add_argument("--restrict-clade", default="Ib",
                    help="Clade prefix to restrict to (default: Ib)")
    ap.add_argument("--null-p", type=float, default=0.70)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--halt-threshold", type=float, default=0.50)
    ap.add_argument("--out-result", required=True, type=Path)
    args = ap.parse_args()

    counts = pd.read_csv(args.in_counts, sep="\t")
    meta = pd.read_csv(args.in_metadata, sep="\t", dtype=str)

    df = counts.merge(meta[["accession", "clade"]] if "clade" in meta.columns else meta[["accession"]],
                       on="accession", how="left")

    n_total = len(df)
    if args.restrict_clade and "clade" in df.columns:
        before = len(df)
        df = df[df["clade"].fillna("").str.startswith(args.restrict_clade)]
        print(f"[restrict] {before:,} -> {len(df):,} genomes after clade prefix '{args.restrict_clade}'",
              file=sys.stderr)
    else:
        print(f"[restrict] no clade restriction applied", file=sys.stderr)

    if len(df) == 0:
        sys.exit("No genomes match the clade restriction; cannot run H1.")

    k = int(df["apobec3_snvs"].sum())
    n = int(df["total_snvs"].sum())
    if n == 0:
        sys.exit("Zero SNVs across the restricted corpus; cannot run H1.")
    observed = k / n

    # Wilson 95% interval
    lo, hi = wilson_interval(k, n, alpha=args.alpha)

    # One-sided exact binomial test (alternative: observed > p0)
    test = binomtest(k=k, n=n, p=args.null_p, alternative="greater")
    p_value = test.pvalue

    supported = lo > args.null_p
    emergency_halt = observed < args.halt_threshold

    result = {
        "n_genomes_in_test": int(len(df)),
        "n_genomes_in_corpus": n_total,
        "restrict_clade_prefix": args.restrict_clade,
        "total_snvs": n,
        "apobec3_snvs": k,
        "non_apobec3_snvs": n - k,
        "observed_apobec3_fraction": observed,
        "wilson_95_lower": lo,
        "wilson_95_upper": hi,
        "alpha": args.alpha,
        "null_p": args.null_p,
        "alternative": "greater",
        "exact_binomial_p_value": p_value,
        "h1_supported": bool(supported),
        "h1_supported_rule": (
            f"Wilson lower 95% bound ({lo:.4f}) > null p0 ({args.null_p})"
        ),
        "emergency_halt_triggered": bool(emergency_halt),
        "emergency_halt_threshold": args.halt_threshold,
        "decision": (
            "EMERGENCY HALT — downstream analysis stops; document investigation"
            if emergency_halt
            else "H1 supported — proceed"
            if supported
            else "H1 not supported but observed >= halt threshold — proceed, report transparently"
        ),
    }

    args.out_result.parent.mkdir(parents=True, exist_ok=True)
    args.out_result.write_text(json.dumps(result, indent=2))

    # Print headline result to stderr for the snakemake log
    print(
        f"\n[H1] observed = {observed:.4f}  (k={k:,} of n={n:,})\n"
        f"[H1] 95% Wilson interval = [{lo:.4f}, {hi:.4f}]\n"
        f"[H1] one-sided exact-binomial P = {p_value:.3e}\n"
        f"[H1] decision: {result['decision']}\n",
        file=sys.stderr,
    )

    if emergency_halt:
        return 2  # distinct exit code for the halt trigger
    return 0


if __name__ == "__main__":
    sys.exit(main())
