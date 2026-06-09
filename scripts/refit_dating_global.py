#!/usr/bin/env python3
"""Global-minimum re-fit of the saturation-aware Poisson MRCA model.

The registered fit (scripts/fit_saturation_dating.py) draws random starts that
never seed the high-rate near-earliest region, so its best-of-25 rule settles in
a ~2019 local minimum ~260 LL units worse than the global ML. This script
locates the global ML by a dense (t0, log10 lambda) grid plus local L-BFGS-B
polish, then bootstraps *within the global basin* so the reported CI reflects the
correct mode. Same data, same L, same likelihood as the registered fit.
"""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize

FREEZE = "2026-05-22"
SEED = 42
ROOT = Path(__file__).resolve().parent.parent
FREEZE_DIR = ROOT / "data" / "processed" / "freeze_20260522"


def parse_loose(s):
    if not isinstance(s, str) or not s.strip():
        return None
    s = s.strip()
    for fmt, take in (("%Y-%m-%d", 10), ("%Y-%m", 7), ("%Y", 4)):
        if len(s) >= take:
            try:
                return datetime.strptime(s[:take], fmt)
            except ValueError:
                pass
    return None


def mu_sat(ages, L, lam):
    return L * (1.0 - np.exp(-lam * np.maximum(ages, 0.0)))


def mu_lin(ages, L, lam):
    return L * lam * np.maximum(ages, 0.0)


def nll(params, k, dates, L, model):
    t0, lam = params
    ages = dates - t0
    if (ages < 0).any():
        return 1e18 + (-ages.clip(max=0)).sum() * 1e6
    mu = model(ages, L, lam)
    if (mu <= 0).any():
        return 1e18
    return -np.sum(k * np.log(mu) - mu)


def global_fit(k, dates, L, model, freeze_ord):
    """Dense grid + local polish. Returns (t0, lam, nll)."""
    earliest = dates.min()
    t0_hi = min(freeze_ord - 30, earliest - 1)
    t0_grid = np.linspace(datetime(1990, 1, 1).toordinal(), t0_hi, 200)
    lam_grid = 10 ** np.linspace(-10, -3, 200)
    best = (None, None, np.inf)
    for t0 in t0_grid:
        for lam in lam_grid:
            f = nll((t0, lam), k, dates, L, model)
            if f < best[2]:
                best = (t0, lam, f)
    bounds = [(datetime(1990, 1, 1).toordinal(), t0_hi), (1e-10, 1e-3)]
    r = minimize(nll, x0=[best[0], best[1]], args=(k, dates, L, model),
                 method="L-BFGS-B", bounds=bounds)
    if r.fun < best[2]:
        best = (r.x[0], r.x[1], r.fun)
    return best


def basin_bootstrap(k, dates, L, t0, lam, model, freeze_ord, n=2000, seed=SEED):
    """Parametric bootstrap; each refit is polished from the global-ML point so
    the CI characterises the global basin rather than the optimiser's pathology."""
    rng = np.random.default_rng(seed)
    mu = model(dates - t0, L, lam)
    earliest = dates.min()
    t0_hi = min(freeze_ord - 30, earliest - 1)
    bounds = [(datetime(1990, 1, 1).toordinal(), t0_hi), (1e-10, 1e-3)]
    out = []
    for _ in range(n):
        kk = rng.poisson(mu)
        r = minimize(nll, x0=[t0, lam], args=(kk, dates, L, model),
                     method="L-BFGS-B", bounds=bounds)
        if r.success and np.isfinite(r.fun):
            out.append(r.x[0])
    return np.array(out)


def main():
    counts = pd.read_csv(FREEZE_DIR / "apobec3_counts_branch.tsv", sep="\t")
    members = pd.read_csv(FREEZE_DIR / "dating" / "cluster_members.tsv",
                          sep="\t", dtype=str)
    L = int(json.loads((FREEZE_DIR / "L_eligible_sites_v2.json").read_text())["L_total"])
    freeze_ord = datetime.strptime(FREEZE, "%Y-%m-%d").toordinal()

    df = counts.merge(members[["accession", "cluster_id", "collection_date"]],
                      on="accession", how="inner")
    df["cd"] = df["collection_date"].apply(parse_loose)
    df = df.dropna(subset=["cd"])
    df["ord"] = df["cd"].apply(lambda d: d.toordinal())

    k = df["apobec3_snvs"].astype(int).to_numpy()
    dates = df["ord"].astype(float).to_numpy()
    print(f"n={len(df)}  L={L:,}  earliest={datetime.fromordinal(int(dates.min())).date()}"
          f"  latest={datetime.fromordinal(int(dates.max())).date()}")

    out = {"L": L, "n_genomes": int(len(df)), "freeze_date": FREEZE}
    for name, model in (("saturation", mu_sat), ("uncorrected", mu_lin)):
        t0, lam, f = global_fit(k, dates, L, model, freeze_ord)
        boot = basin_bootstrap(k, dates, L, t0, lam, model, freeze_ord)
        lo, hi = np.percentile(boot, [2.5, 97.5])
        out[name] = {
            "t0_hat": datetime.fromordinal(int(round(t0))).isoformat()[:10],
            "lambda_per_site_per_day": float(lam),
            "nll": float(f),
            "ci95_lo": datetime.fromordinal(int(round(lo))).isoformat()[:10],
            "ci95_hi": datetime.fromordinal(int(round(hi))).isoformat()[:10],
            "n_boot": int(len(boot)),
        }
        print(f"[{name}] t0={out[name]['t0_hat']}  lam={lam:.3e}  NLL={f:.2f}  "
              f"CI {out[name]['ci95_lo']}..{out[name]['ci95_hi']}")

    if "saturation" in out and "uncorrected" in out:
        ts = datetime.strptime(out["saturation"]["t0_hat"], "%Y-%m-%d").toordinal()
        tu = datetime.strptime(out["uncorrected"]["t0_hat"], "%Y-%m-%d").toordinal()
        out["delta_days_global_basin"] = abs(ts - tu)
        print(f"[H2] global-basin Delta = {out['delta_days_global_basin']} days")

    dest = FREEZE_DIR / "dating" / "saturation_dating_global_refit.json"
    dest.write_text(json.dumps(out, indent=2))
    print(f"-> {dest}")


if __name__ == "__main__":
    main()
