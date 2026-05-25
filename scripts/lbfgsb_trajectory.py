#!/usr/bin/env python3
"""
Capture the L-BFGS-B convergence trajectory of the 25 random starts used by
fit_saturation_dating.py on the §5.3 saturation-aware Poisson model.

Output: a CSV with columns
    start_idx, t0_init_ordinal, lam_init, t0_init_year,
    t0_converged_ordinal, lam_converged, t0_converged_year,
    final_nll, converged_status, basin

where 'basin' classifies the converged point as ~2019 (within 1 year of
2019-01-09) or ~2023 (within 1 year of 2023-10-01) or 'other'.

Supports §6.4 Amendment 04 (post-hoc local-minimum-trap diagnostic).
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize


SEED = 42


def mu_saturation(t_days, L, lam):
    return L * (1.0 - np.exp(-lam * t_days))


def neg_log_lik(params, k, dates_ord, L):
    t0, lam = params
    ages = dates_ord - t0
    if (ages <= 0).any():
        return 1e18
    mu = mu_saturation(ages, L, lam)
    if (mu <= 0).any():
        return 1e18
    ll = np.sum(k * np.log(mu) - mu)
    return -ll


def parse_date(s):
    for fmt, take in (("%Y-%m-%d", 10), ("%Y-%m", 7), ("%Y", 4)):
        if len(s) < take:
            continue
        try:
            return datetime.strptime(s[:take], fmt)
        except ValueError:
            continue
    return None


def ord_to_year(o):
    d = datetime.fromordinal(int(o))
    start = datetime(d.year, 1, 1).toordinal()
    end = datetime(d.year + 1, 1, 1).toordinal()
    return d.year + (o - start) / (end - start)


def classify_basin(t0_year):
    if 2018.0 <= t0_year <= 2020.0:
        return "~2019 local minimum"
    if 2022.5 <= t0_year <= 2024.5:
        return "~2023 global minimum"
    return "other"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-cluster-members", required=True, type=Path)
    ap.add_argument("--in-branch-counts", required=True, type=Path)
    ap.add_argument("--in-results", required=True, type=Path)
    ap.add_argument("--n-starts", type=int, default=25)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--out-csv", required=True, type=Path)
    args = ap.parse_args()

    members = pd.read_csv(args.in_cluster_members, sep="\t")
    counts = pd.read_csv(args.in_branch_counts, sep="\t")
    results = json.loads(args.in_results.read_text())
    L_sites = int(results["L"])

    merged = members.merge(counts, on="accession", how="inner")
    merged["dt"] = merged["collection_date"].astype(str).map(parse_date)
    merged = merged.dropna(subset=["dt"])
    merged["t_ordinal"] = merged["dt"].map(lambda d: d.toordinal())
    k_arr = merged["apobec3_snvs"].to_numpy(dtype=float)
    dates_arr = merged["t_ordinal"].to_numpy(dtype=float)

    rng = np.random.default_rng(args.seed)
    earliest = dates_arr.min()
    freeze_ord = datetime(2026, 5, 22).toordinal()
    bounds = [
        (datetime(1990, 1, 1).toordinal(), freeze_ord - 30),
        (1e-10, 1e-3),
    ]

    rows = []
    for i in range(args.n_starts):
        t0_init = rng.uniform(bounds[0][0], min(bounds[0][1], earliest - 1))
        lam_init = 10 ** rng.uniform(-10, -3)
        try:
            r = minimize(
                neg_log_lik,
                x0=[t0_init, lam_init],
                args=(k_arr, dates_arr, L_sites),
                method="L-BFGS-B",
                bounds=bounds,
            )
            t0_conv = float(r.x[0])
            lam_conv = float(r.x[1])
            nll_conv = float(r.fun)
            converged = bool(r.success)
        except Exception as e:
            t0_conv = lam_conv = float("nan")
            nll_conv = float("inf")
            converged = False

        t0_init_year = ord_to_year(t0_init)
        t0_conv_year = ord_to_year(t0_conv) if np.isfinite(t0_conv) else float("nan")
        basin = classify_basin(t0_conv_year) if np.isfinite(t0_conv_year) else "non-converged"

        rows.append({
            "start_idx": i,
            "t0_init_ordinal": int(t0_init),
            "lam_init": lam_init,
            "t0_init_year": t0_init_year,
            "t0_converged_ordinal": int(t0_conv) if np.isfinite(t0_conv) else None,
            "t0_converged_date": datetime.fromordinal(int(t0_conv)).strftime("%Y-%m-%d") if np.isfinite(t0_conv) else "non-converged",
            "lam_converged": lam_conv,
            "t0_converged_year": t0_conv_year,
            "final_nll": nll_conv,
            "converged_status": converged,
            "basin": basin,
        })

    df = pd.DataFrame(rows)
    df.to_csv(args.out_csv, index=False)
    print(f"[done] wrote {args.out_csv} ({len(df)} rows)")

    # Summary
    print(f"\n=== Summary across {args.n_starts} random starts ===")
    print(df.groupby("basin").size().to_string())
    print()
    print("Best (lowest NLL) per basin:")
    for basin, sub in df.groupby("basin"):
        best = sub.loc[sub["final_nll"].idxmin()]
        print(f"  {basin:30} best NLL={best['final_nll']:.2f}  t0={best['t0_converged_date']}  lam={best['lam_converged']:.3e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
