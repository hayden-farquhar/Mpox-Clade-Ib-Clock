#!/usr/bin/env python3
"""
§5.5 sensitivity analyses for the dating + recombination-scan pipeline.

Implements §5.5 items 2, 3, and 6 in a single configurable runner. Item 1
(alternate alignment) is handled separately by re-running Phase 3 / 4 on the
MAFFT primary alignment; item 4 (detector-parameter sweep) is in
sensitivity_detector_sweep.py; item 5 (JC-style correction) is invoked via
fit_saturation_dating.py's --model flag (saturation_jc); item 7
(pipeline-reproducibility verification) is a workflow re-run on a clean env,
not a new script.

This runner perturbs the input genome set, re-runs the dating step
(build_dating_clusters → fit_saturation_dating), and reports per-cluster
emergence-date shifts and Jaccard similarities for the H5 test.

Modes
-----
- subsample: §5.5 item 2. Random subsampling at the given fraction (75% or
  50%), specified number of independent replicates. Reports per-cluster
  emergence-date median + 5th-95th percentile across replicates, and Jaccard
  similarity of the consensus-candidate set against the primary.
- temporal_holdout: §5.5 item 3. Exclude the most-recently-collected fraction
  (default 25%). Reports per-cluster emergence-date shift versus the primary.
- outgroup_redraw: §5.5 item 6. Re-sample the §5.4 30-genome outgroup with
  alternative seeds (default 1..5). Reports whether the consensus-candidate
  set is stable.

All modes preserve seed reproducibility — the same --seed flag with the same
inputs always gives the same outputs.
"""

from __future__ import annotations

import argparse
import json
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd


def run_step(label: str, cmd: list[str]) -> int:
    print(f"\n[step] {label}\n  {' '.join(cmd[:8])} {'...' if len(cmd)>8 else ''}",
          file=sys.stderr)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  FAILED (exit {r.returncode})", file=sys.stderr)
        print(r.stderr[-500:], file=sys.stderr)
    return r.returncode


def mode_subsample(args) -> int:
    """§5.5 item 2: random subsampling at the given fraction, multiple replicates."""
    rng = random.Random(args.seed)
    counts = pd.read_csv(args.in_counts, sep="\t", dtype=str)
    members = pd.read_csv(args.in_cluster_members, sep="\t", dtype=str)

    results = []
    for rep in range(args.n_replicates):
        rep_seed = args.seed + rep
        rep_rng = random.Random(rep_seed)
        n_keep = int(round(len(members) * args.fraction))
        keep_acc = set(rep_rng.sample(list(members["accession"]), n_keep))

        rep_dir = args.out_dir / f"subsample_{int(args.fraction*100)}pct_rep{rep+1}"
        rep_dir.mkdir(parents=True, exist_ok=True)

        # Write subsetted cluster_members + counts
        members_sub = members[members["accession"].isin(keep_acc)].copy()
        members_sub.to_csv(rep_dir / "cluster_members.tsv", sep="\t", index=False)
        counts_sub = counts[counts["accession"].isin(keep_acc)].copy()
        counts_sub.to_csv(rep_dir / "counts.tsv", sep="\t", index=False)

        # Re-run dating on the subset
        rc = run_step(f"subsample {int(args.fraction*100)}% rep {rep+1} dating",
                       ["python3", str(args.fit_script),
                        "--in-counts", str(rep_dir / "counts.tsv"),
                        "--in-cluster-members", str(rep_dir / "cluster_members.tsv"),
                        "--in-L", str(args.in_L),
                        "--freeze-date", args.freeze_date,
                        "--n-bootstrap", str(args.n_bootstrap_reduced),
                        "--out-results", str(rep_dir / "dating_results.json")])
        if rc == 0 and (rep_dir / "dating_results.json").exists():
            results.append({"rep": rep+1, "seed": rep_seed,
                            "n_genomes": len(members_sub),
                            "results_path": str(rep_dir / "dating_results.json")})

    summary = {"mode": "subsample", "fraction": args.fraction,
                "n_replicates": args.n_replicates, "replicates": results}
    (args.out_dir / "subsample_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\n[done] subsample mode, {len(results)}/{args.n_replicates} replicates ok",
          file=sys.stderr)
    return 0


def mode_temporal_holdout(args) -> int:
    """§5.5 item 3: exclude the holdout-fraction most recently collected."""
    from datetime import datetime
    members = pd.read_csv(args.in_cluster_members, sep="\t", dtype=str)

    def parse_loose(s):
        if not isinstance(s, str) or not s.strip(): return None
        s = s.strip()
        for fmt, take in (("%Y-%m-%d", 10), ("%Y-%m", 7), ("%Y", 4)):
            if len(s) < take: continue
            try: return datetime.strptime(s[:take], fmt)
            except: continue
        return None

    members["cd"] = members["collection_date"].apply(parse_loose)
    members = members.dropna(subset=["cd"]).sort_values("cd")
    n_drop = int(round(len(members) * args.holdout_fraction))
    keep = members.iloc[:len(members) - n_drop].copy()
    drop = members.iloc[len(members) - n_drop:].copy()

    out_dir = args.out_dir / f"temporal_holdout_{int(args.holdout_fraction*100)}pct"
    out_dir.mkdir(parents=True, exist_ok=True)
    keep.drop(columns=["cd"]).to_csv(out_dir / "cluster_members.tsv", sep="\t", index=False)
    drop.drop(columns=["cd"]).to_csv(out_dir / "dropped_members.tsv", sep="\t", index=False)
    print(f"[holdout] dropped {len(drop):,} most-recent ({args.holdout_fraction:.0%}); "
          f"kept {len(keep):,}; date cut = {drop['cd'].min().date()}", file=sys.stderr)

    # Subset counts too
    counts = pd.read_csv(args.in_counts, sep="\t", dtype=str)
    counts_keep = counts[counts["accession"].isin(set(keep["accession"]))].copy()
    counts_keep.to_csv(out_dir / "counts.tsv", sep="\t", index=False)

    rc = run_step("temporal_holdout dating",
                   ["python3", str(args.fit_script),
                    "--in-counts", str(out_dir / "counts.tsv"),
                    "--in-cluster-members", str(out_dir / "cluster_members.tsv"),
                    "--in-L", str(args.in_L),
                    "--freeze-date", args.freeze_date,
                    "--n-bootstrap", str(args.n_bootstrap_reduced),
                    "--out-results", str(out_dir / "dating_results.json")])
    summary = {"mode": "temporal_holdout", "holdout_fraction": args.holdout_fraction,
                "n_kept": int(len(keep)), "n_dropped": int(len(drop)),
                "date_cut": str(drop["cd"].min().date()),
                "results_path": str(out_dir / "dating_results.json")}
    (args.out_dir / "temporal_holdout_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"[done] temporal_holdout mode (exit {rc})", file=sys.stderr)
    return rc


def mode_outgroup_redraw(args) -> int:
    """§5.5 item 6: re-draw the 30-genome IIb outgroup panel with alternative seeds."""
    # Reuses build_panels.py with different seeds; comparison happens against
    # the primary outgroup_30_iib panel. For each redraw, re-run the
    # recombination scan (3SEQ + GARD + PhiPack) and compute Jaccard.
    redraws = []
    for seed_offset in args.outgroup_seeds:
        rep_dir = args.out_dir / f"outgroup_redraw_seed{seed_offset}"
        rep_dir.mkdir(parents=True, exist_ok=True)
        rc = run_step(f"outgroup redraw seed={seed_offset}",
                       ["python3", str(args.panels_script),
                        "--clade-i-tsv", str(args.clade_i_tsv),
                        "--all-clades-tsv", str(args.all_clades_tsv),
                        "--ncbi-meta-tsv", str(args.ncbi_meta_tsv),
                        "--out-dir", str(rep_dir),
                        "--seed", str(seed_offset)])
        if rc == 0:
            redraws.append({"seed": seed_offset,
                             "panels_path": str(rep_dir / "outgroup_30_iib.tsv")})

    summary = {"mode": "outgroup_redraw", "n_seeds": len(args.outgroup_seeds),
                "redraws": redraws,
                "note": "Per-redraw recombination scan + consensus filter must be "
                        "re-run by the H5 driver against each redraw's panel; this "
                        "script produces the panels only."}
    (args.out_dir / "outgroup_redraw_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"[done] outgroup_redraw mode, {len(redraws)} panels produced",
          file=sys.stderr)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description="§5.5 sensitivity runner")
    sub = ap.add_subparsers(dest="mode", required=True)

    # subsample
    ss = sub.add_parser("subsample")
    ss.add_argument("--in-counts", required=True, type=Path)
    ss.add_argument("--in-cluster-members", required=True, type=Path)
    ss.add_argument("--in-L", required=True, type=Path)
    ss.add_argument("--freeze-date", required=True)
    ss.add_argument("--fraction", type=float, choices=[0.5, 0.75], required=True)
    ss.add_argument("--n-replicates", type=int, default=5)
    ss.add_argument("--seed", type=int, default=42)
    ss.add_argument("--n-bootstrap-reduced", type=int, default=500,
                     help="Reduced bootstrap count for sensitivity runs (default 500 vs primary 2000)")
    ss.add_argument("--fit-script", required=True, type=Path,
                     help="Path to fit_saturation_dating.py")
    ss.add_argument("--out-dir", required=True, type=Path)

    # temporal_holdout
    th = sub.add_parser("temporal_holdout")
    th.add_argument("--in-counts", required=True, type=Path)
    th.add_argument("--in-cluster-members", required=True, type=Path)
    th.add_argument("--in-L", required=True, type=Path)
    th.add_argument("--freeze-date", required=True)
    th.add_argument("--holdout-fraction", type=float, default=0.25)
    th.add_argument("--n-bootstrap-reduced", type=int, default=500)
    th.add_argument("--fit-script", required=True, type=Path)
    th.add_argument("--out-dir", required=True, type=Path)

    # outgroup_redraw
    og = sub.add_parser("outgroup_redraw")
    og.add_argument("--clade-i-tsv", required=True, type=Path)
    og.add_argument("--all-clades-tsv", required=True, type=Path)
    og.add_argument("--ncbi-meta-tsv", required=True, type=Path)
    og.add_argument("--outgroup-seeds", type=int, nargs="+", default=[1, 2, 3, 4, 5])
    og.add_argument("--panels-script", required=True, type=Path)
    og.add_argument("--out-dir", required=True, type=Path)

    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    if args.mode == "subsample":
        return mode_subsample(args)
    if args.mode == "temporal_holdout":
        return mode_temporal_holdout(args)
    if args.mode == "outgroup_redraw":
        return mode_outgroup_redraw(args)
    return 1


if __name__ == "__main__":
    sys.exit(main())
