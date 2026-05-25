#!/usr/bin/env python3
"""
Fit the §5.3 saturation-aware Poisson dating model to a clade-Ib cluster.

Models
------
Saturation-aware (primary):
    E[edited sites at t] = L · (1 − exp(−λ · t))
Uncorrected comparator (linearisation, valid as λ·t → 0):
    E[edits at t]        = L · λ · t

Per-genome counts are treated as Poisson around the model expectation.
Genomes whose collection date is before t₀ are penalised: optimiser cannot
choose t₀ > earliest collection date in the cluster.

Likelihood
----------
L(λ, t₀ | data) = ∏_i  Poisson(k_i; μ_i(λ, t₀))

where k_i is the per-genome APOBEC3 SNV count (branch-quantity per Amendment
02 H1ʹ) and μ_i is the model expectation at age (collection_date_i − t₀).

Fit
---
scipy.optimize.minimize, method = L-BFGS-B, with bounds:
    t₀ ∈ [1990-01-01, freeze_date − 30 days]
    λ  ∈ [1e-10, 1e-3] per site per day
25 random starting points, best-LL kept.

Uncertainty
-----------
Parametric bootstrap with 2000 resamples. For each iteration, per-genome
counts are re-drawn from the fitted-model's expectation and the corresponding
model is re-fit. Reported 95% interval = 2.5th–97.5th percentiles of t₀.

Identifiability check
---------------------
Profile log-likelihood for t₀ is inspected. Cluster is flagged "poorly
identified" if profile is flatter than 2 LL units across ±1 year of t₀_hat.

H2 test
-------
Δ = |t₀_corrected − t₀_uncorrected| in days. H2 supported if Δ > 60 days on
the cluster (or max-over-clusters Δ > 60 days when multiple clusters present).
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import poisson

ACGT = set("ACGT")
SEED = 42


def parse_loose(s):
    if not isinstance(s, str) or not s.strip():
        return None
    s = s.strip()
    for fmt, take in (("%Y-%m-%d", 10), ("%Y-%m", 7), ("%Y", 4)):
        if len(s) < take:
            continue
        try:
            return datetime.strptime(s[:take], fmt)
        except ValueError:
            continue
    return None


def mu_saturation(t_days, L, lam):
    """E[edited sites at t] under saturation model."""
    t = np.maximum(t_days, 0.0)
    return L * (1.0 - np.exp(-lam * t))


def mu_linear(t_days, L, lam):
    """E[edits at t] under uncorrected linearisation."""
    t = np.maximum(t_days, 0.0)
    return L * lam * t


def neg_log_lik(params, k, dates, L, model):
    """
    params = (t0_ordinal_days_since_epoch, lam)
    Returns -log-likelihood. Returns +inf if any genome has collection_date < t0
    (incompatible with the model — penalised via inf).
    """
    t0, lam = params
    ages = dates - t0  # days
    if (ages < 0).any():
        return 1e18 + (-ages.clip(max=0)).sum() * 1e6  # smooth penalty
    mu = model(ages, L, lam) if callable(model) else mu_saturation(ages, L, lam)
    if (mu <= 0).any():
        return 1e18
    # Poisson log-likelihood
    ll = np.sum(k * np.log(mu) - mu)  # constant terms dropped
    return -ll


def fit_one(k_arr, dates_arr, L, model, freeze_date_ordinal, n_starts=25, seed=SEED):
    """Fit a single model with `n_starts` random restarts. Returns (best_params, best_nll)."""
    rng = np.random.default_rng(seed)
    earliest = dates_arr.min()
    bounds = [(datetime(1990, 1, 1).toordinal(), freeze_date_ordinal - 30),
              (1e-10, 1e-3)]
    best = None
    for _ in range(n_starts):
        t0_init = rng.uniform(bounds[0][0], min(bounds[0][1], earliest - 1))
        lam_init = 10 ** rng.uniform(-10, -3)
        try:
            r = minimize(neg_log_lik, x0=[t0_init, lam_init],
                          args=(k_arr, dates_arr, L, model),
                          method="L-BFGS-B", bounds=bounds)
            if best is None or r.fun < best.fun:
                best = r
        except Exception:
            continue
    return best


def parametric_bootstrap(k_arr, dates_arr, L, fitted_t0, fitted_lam, model,
                          freeze_date_ordinal, n_iter=2000, seed=SEED):
    """
    Parametric bootstrap: redraw counts from fitted-model expectation and refit.
    Returns array of bootstrapped t0 estimates.
    """
    rng = np.random.default_rng(seed)
    ages = dates_arr - fitted_t0
    mu = model(ages, L, fitted_lam)
    t0_boot = []
    for it in range(n_iter):
        k_resampled = rng.poisson(mu)
        r = fit_one(k_resampled, dates_arr, L, model, freeze_date_ordinal,
                     n_starts=5, seed=seed + it)
        if r is not None and np.isfinite(r.fun):
            t0_boot.append(r.x[0])
    return np.array(t0_boot)


def profile_ll(k_arr, dates_arr, L, fitted_t0, model, freeze_date_ordinal,
                halfwidth_days=365):
    """
    Profile log-likelihood across t0 ∈ [fitted_t0 − halfwidth, fitted_t0 + halfwidth].
    For each candidate t0, optimise λ alone. Returns (t0_grid, profile_nll).
    Used for the identifiability check (§5.3: poorly identified if profile is
    flatter than 2 LL units across ±1 year).
    """
    grid = np.linspace(fitted_t0 - halfwidth_days, fitted_t0 + halfwidth_days, 41)
    earliest = dates_arr.min()
    grid = grid[grid < earliest]  # respect t0 < earliest constraint
    prof = []
    for t0 in grid:
        def nll_lam(lam):
            ages = dates_arr - t0
            if (ages < 0).any():
                return 1e18
            mu = model(ages, L, lam[0])
            if (mu <= 0).any():
                return 1e18
            return -np.sum(k_arr * np.log(mu) - mu)
        r = minimize(nll_lam, x0=[1e-6], method="L-BFGS-B", bounds=[(1e-10, 1e-3)])
        prof.append(r.fun if r.success else np.nan)
    return grid, np.array(prof)


def main() -> int:
    ap = argparse.ArgumentParser(description="§5.3 saturation-aware Poisson dating.")
    ap.add_argument("--in-counts", required=True, type=Path,
                    help="Per-tip branch-quantity APOBEC3 counts TSV "
                         "(from build_branch_apobec3_counts.py)")
    ap.add_argument("--in-cluster-members", required=True, type=Path,
                    help="cluster_members.tsv from build_dating_clusters.py")
    ap.add_argument("--in-L", required=True, type=Path,
                    help="L_eligible_sites.json from count_apobec3.py")
    ap.add_argument("--freeze-date", required=True,
                    help="Freeze date ISO (YYYY-MM-DD); used as upper bound on t0")
    ap.add_argument("--n-bootstrap", type=int, default=2000)
    ap.add_argument("--n-starts", type=int, default=25)
    ap.add_argument("--out-results", required=True, type=Path)
    args = ap.parse_args()

    counts = pd.read_csv(args.in_counts, sep="\t")
    members = pd.read_csv(args.in_cluster_members, sep="\t", dtype=str)
    L_json = json.loads(args.in_L.read_text())
    L = int(L_json["L_total"])
    freeze_ord = datetime.strptime(args.freeze_date, "%Y-%m-%d").toordinal()

    print(f"[fit] L = {L:,} eligible APOBEC3 target sites", file=sys.stderr)
    print(f"[fit] freeze date = {args.freeze_date} (ord {freeze_ord})", file=sys.stderr)
    print(f"[fit] bootstrap = {args.n_bootstrap}, random starts = {args.n_starts}", file=sys.stderr)

    df = counts.merge(members[["accession", "cluster_id", "collection_date"]],
                       on="accession", how="inner")
    df["cd_parsed"] = df["collection_date"].apply(parse_loose)
    df = df.dropna(subset=["cd_parsed"])
    df["cd_ord"] = df["cd_parsed"].apply(lambda d: d.toordinal())

    results = {"L": L, "freeze_date": args.freeze_date, "clusters": []}

    for cid, sub in df.groupby("cluster_id"):
        if len(sub) < 10:
            print(f"[skip] cluster {cid}: n={len(sub)} < 10", file=sys.stderr)
            continue
        k_arr = sub["apobec3_snvs"].astype(int).to_numpy()
        dates_arr = sub["cd_ord"].astype(float).to_numpy()
        print(f"[cluster] {cid}: n={len(sub)}, "
              f"earliest={datetime.fromordinal(int(dates_arr.min())).date()}, "
              f"latest={datetime.fromordinal(int(dates_arr.max())).date()}, "
              f"k range={k_arr.min()}-{k_arr.max()}", file=sys.stderr)

        out_cluster = {"cluster_id": cid, "n_genomes": int(len(sub))}

        for name, model in (("saturation", mu_saturation), ("uncorrected", mu_linear)):
            print(f"  fitting {name} model ...", file=sys.stderr)
            best = fit_one(k_arr, dates_arr, L, model, freeze_ord,
                            n_starts=args.n_starts, seed=SEED)
            if best is None:
                out_cluster[name] = {"status": "fit_failed"}
                continue
            t0_hat, lam_hat = best.x
            print(f"    t0_hat = {datetime.fromordinal(int(t0_hat)).date()}, "
                  f"lam_hat = {lam_hat:.3e} per site per day", file=sys.stderr)

            # Bootstrap
            print(f"    parametric bootstrap n={args.n_bootstrap} ...", file=sys.stderr)
            boot = parametric_bootstrap(k_arr, dates_arr, L, t0_hat, lam_hat, model,
                                          freeze_ord, n_iter=args.n_bootstrap, seed=SEED)
            if len(boot) < 10:
                out_cluster[name] = {"status": "bootstrap_failed",
                                       "t0_hat": datetime.fromordinal(int(t0_hat)).isoformat(),
                                       "lambda_hat": lam_hat}
                continue
            ci_lo, ci_hi = np.percentile(boot, [2.5, 97.5])

            # Profile-LL identifiability check
            grid, prof = profile_ll(k_arr, dates_arr, L, t0_hat, model, freeze_ord)
            in_window = (grid >= t0_hat - 365) & (grid <= t0_hat + 365)
            profile_range = (prof[in_window].max() - prof[in_window].min()) if in_window.any() else float("nan")
            poorly_identified = bool(profile_range < 2.0) if np.isfinite(profile_range) else False

            out_cluster[name] = {
                "status": "fit_ok",
                "t0_hat": datetime.fromordinal(int(t0_hat)).isoformat(),
                "t0_hat_ordinal": float(t0_hat),
                "lambda_hat_per_site_per_day": float(lam_hat),
                "ci95_lo": datetime.fromordinal(int(ci_lo)).isoformat(),
                "ci95_hi": datetime.fromordinal(int(ci_hi)).isoformat(),
                "n_bootstrap_successful": int(len(boot)),
                "profile_range_LL_units": float(profile_range) if np.isfinite(profile_range) else None,
                "poorly_identified": poorly_identified,
            }
            print(f"    95% CI: {out_cluster[name]['ci95_lo']} -- {out_cluster[name]['ci95_hi']}", file=sys.stderr)
            print(f"    profile-LL range (±1y) = {profile_range:.2f} LL units; "
                  f"{'POORLY IDENTIFIED' if poorly_identified else 'identifiable'}", file=sys.stderr)

        # H2 statistic for this cluster
        if (out_cluster.get("saturation", {}).get("status") == "fit_ok"
            and out_cluster.get("uncorrected", {}).get("status") == "fit_ok"):
            t_sat = out_cluster["saturation"]["t0_hat_ordinal"]
            t_unc = out_cluster["uncorrected"]["t0_hat_ordinal"]
            out_cluster["delta_days"] = float(abs(t_sat - t_unc))
            out_cluster["delta_signed_days"] = float(t_sat - t_unc)
            out_cluster["h2_supported_this_cluster"] = bool(out_cluster["delta_days"] > 60)

        results["clusters"].append(out_cluster)

    # H2 max over clusters
    deltas = [c["delta_days"] for c in results["clusters"] if "delta_days" in c]
    results["h2_max_delta_days"] = float(max(deltas)) if deltas else None
    results["h2_supported"] = bool(results["h2_max_delta_days"] is not None
                                    and results["h2_max_delta_days"] > 60)

    args.out_results.parent.mkdir(parents=True, exist_ok=True)
    args.out_results.write_text(json.dumps(results, indent=2))
    print(f"\n[done] results -> {args.out_results}", file=sys.stderr)
    print(f"[H2] max-over-clusters Δ = {results['h2_max_delta_days']} days; "
          f"H2 supported? {results['h2_supported']}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
