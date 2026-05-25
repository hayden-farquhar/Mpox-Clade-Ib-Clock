#!/usr/bin/env python3
"""
Regenerate the (λ, t₀) likelihood-landscape figure for the §5.3 saturation
Poisson dating analysis, with title aligned to the local-minimum-trap framing.

The figure shows the absolute log-likelihood of the saturation-Poisson model
on a (λ, t₀) grid, with the registered L-BFGS-B converged primary fit (~2019
local minimum) and the global ML minimum (~2023) annotated. The ΔLL scale
shows the registered primary is ~100 LL units worse than the global minimum.
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def neg_loglik(t0_ordinal, lam, t_ages_days, k_counts, L_sites):
    if (t_ages_days <= 0).any():
        return np.inf
    mu = L_sites * (1.0 - np.exp(-lam * t_ages_days))
    mu = np.clip(mu, 1e-300, None)
    # Poisson NLL up to additive constants
    return float(np.sum(mu - k_counts * np.log(mu)))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-cluster-members", required=True, type=Path)
    ap.add_argument("--in-branch-counts", required=True, type=Path,
                    help="apobec3_counts_branch.tsv from build_branch_apobec3_counts.py")
    ap.add_argument("--in-results", required=True, type=Path,
                    help="saturation_dating_results.json")
    ap.add_argument("--out-fig", required=True, type=Path)
    ap.add_argument("--grid-n-t", type=int, default=200)
    ap.add_argument("--grid-n-lam", type=int, default=200)
    ap.add_argument("--t-min-year", type=float, default=2017.0)
    ap.add_argument("--t-max-year", type=float, default=2024.0)
    ap.add_argument("--lam-min", type=float, default=1e-8)
    ap.add_argument("--lam-max", type=float, default=1e-5)
    args = ap.parse_args()

    members = pd.read_csv(args.in_cluster_members, sep="\t")
    counts = pd.read_csv(args.in_branch_counts, sep="\t")
    results = json.loads(args.in_results.read_text())

    cluster = results["clusters"][0]
    L_sites = int(results["L"])

    def parse_date(s: str) -> datetime | None:
        for fmt, take in (("%Y-%m-%d", 10), ("%Y-%m", 7), ("%Y", 4)):
            if len(s) < take:
                continue
            try:
                return datetime.strptime(s[:take], fmt)
            except ValueError:
                continue
        return None

    merged = members.merge(counts, on="accession", how="inner")
    merged["dt"] = merged["collection_date"].astype(str).map(parse_date)
    merged = merged.dropna(subset=["dt"])
    merged["t_ordinal"] = merged["dt"].map(lambda d: d.toordinal())
    k_counts = merged["apobec3_snvs"].to_numpy(dtype=float)
    t_ordinal = merged["t_ordinal"].to_numpy(dtype=float)

    # Build grid
    t_min_ord = datetime(int(args.t_min_year), 1, 1).toordinal()
    t_max_ord = datetime(int(args.t_max_year), 1, 1).toordinal()
    t0_grid_ord = np.linspace(t_min_ord, t_max_ord, args.grid_n_t)
    lam_grid = np.logspace(np.log10(args.lam_min), np.log10(args.lam_max), args.grid_n_lam)

    nll = np.full((args.grid_n_lam, args.grid_n_t), np.inf)
    for i, lam in enumerate(lam_grid):
        for j, t0_ord in enumerate(t0_grid_ord):
            ages = t_ordinal - t0_ord
            if (ages <= 0).any():
                continue
            nll[i, j] = neg_loglik(t0_ord, lam, ages, k_counts, L_sites)

    nll_min = float(np.nanmin(nll))
    delta_nll = nll - nll_min  # delta-NLL = delta-LL (positive)
    j_min, i_min = np.unravel_index(np.argmin(nll), nll.shape)
    print(f"[grid] min ΔLL=0 at lam={lam_grid[j_min]:.3e}, t0_ord={t0_grid_ord[i_min]:.1f} ({datetime.fromordinal(int(t0_grid_ord[i_min])).strftime('%Y-%m-%d')})")

    # Registered primary fit
    t0_primary_ord = cluster["saturation"]["t0_hat_ordinal"]
    lam_primary = cluster["saturation"]["lambda_hat_per_site_per_day"]
    # ΔLL of primary vs global grid minimum
    ages_p = t_ordinal - t0_primary_ord
    nll_primary = neg_loglik(t0_primary_ord, lam_primary, ages_p, k_counts, L_sites)
    delta_primary = nll_primary - nll_min
    print(f"[primary] L-BFGS-B primary fit: t0={cluster['saturation']['t0_hat']}, lam={lam_primary:.3e}, ΔLL from grid global min = {delta_primary:.1f}")

    # Compute best-λ-for-each-t0 ridge line
    best_lam_for_t0 = np.array([lam_grid[np.argmin(nll[:, j])] for j in range(args.grid_n_t)])

    # Convert t0 ordinals to fractional years for plotting
    def ord_to_year(o):
        d = datetime.fromordinal(int(o))
        start = datetime(d.year, 1, 1).toordinal()
        end = datetime(d.year + 1, 1, 1).toordinal()
        return d.year + (o - start) / (end - start)

    t_years = np.array([ord_to_year(o) for o in t0_grid_ord])
    t_primary_year = ord_to_year(t0_primary_ord)
    t_global_year = ord_to_year(t0_grid_ord[i_min])

    # Plot
    fig, ax = plt.subplots(figsize=(11, 6.5))

    levels = [0, 2, 5, 12, 35, 100, 300]
    cmap = plt.get_cmap("viridis_r")
    cs = ax.contourf(t_years, lam_grid, delta_nll, levels=levels, cmap=cmap, extend="max")

    cb = fig.colorbar(cs, ax=ax, label=r"$\Delta$LL from global ML minimum (lower = better)", extend="max")
    cb.set_ticks([0, 2, 5, 12, 35, 100, 300])

    # Best-λ ridge line
    ax.plot(t_years, best_lam_for_t0, ls="--", color="white", lw=1.2, alpha=0.7,
            label=r"Best $\lambda$ for each $t_0$ (ridge trace)")

    # Annotate registered primary
    ax.plot(t_primary_year, lam_primary, marker="o", markersize=14,
            markerfacecolor="red", markeredgecolor="white", markeredgewidth=2,
            label=f"Registered L-BFGS-B primary (~2019; ΔLL ≈ {delta_primary:.0f})")

    # Annotate global minimum
    ax.plot(t_global_year, lam_grid[j_min], marker="s", markersize=14,
            markerfacecolor="gold", markeredgecolor="black", markeredgewidth=2,
            label=f"Global ML minimum on grid (~{datetime.fromordinal(int(t0_grid_ord[i_min])).strftime('%Y-%m')}; ΔLL = 0)")

    # Annotate TreeTime cross-check
    treetime_year = 2023 + (datetime(2023, 6, 2).toordinal() - datetime(2023, 1, 1).toordinal()) / 365.25
    # Find TreeTime-implied lambda from the ridge
    j_treetime = np.argmin(np.abs(t_years - treetime_year))
    treetime_lam_implied = best_lam_for_t0[j_treetime]
    ax.plot(treetime_year, treetime_lam_implied, marker="^", markersize=14,
            markerfacecolor="cyan", markeredgecolor="black", markeredgewidth=1.5,
            label="TreeTime strict-clock cross-check (2023-06-02)")

    ax.set_yscale("log")
    ax.set_xlabel(r"Cluster MRCA date $t_0$", fontsize=12)
    ax.set_ylabel(r"APOBEC3 substitution rate $\lambda$ (per site per day)", fontsize=12)
    ax.set_title(
        r"Saturation-Poisson $(\lambda, t_0)$ likelihood landscape:" "\n"
        "registered primary is a local minimum; data prefer ~2023 (agreeing with TreeTime)",
        fontsize=11,
    )
    ax.grid(True, ls="--", alpha=0.3, color="white")
    ax.legend(loc="lower left", fontsize=9, framealpha=0.9)

    args.out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(args.out_fig, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[done] -> {args.out_fig}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
