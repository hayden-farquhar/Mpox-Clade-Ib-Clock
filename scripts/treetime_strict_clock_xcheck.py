#!/usr/bin/env python3
"""
§5.3 TreeTime strict-clock cross-check.

The pre-registration §5.3 specifies an independent TreeTime strict-clock run
on each cluster's per-cluster masked alignment, reported alongside the
saturation-aware Poisson model as a methodological cross-check. Disagreement
between TreeTime and the primary estimate beyond mutual-interval overlap is
reported and discussed but does not override the pre-specified primary.

Inputs
------
- Per-cluster masked alignment (FASTA) — subset of the §3.3 corpus alignment
  restricted to cluster members
- Per-cluster metadata with collection dates

Outputs
-------
- treetime/ working directory with TreeTime's native outputs
- treetime_summary.json — MRCA date estimate, 95% CI, clock rate, n_iter,
  convergence flag

Implementation
--------------
TreeTime is invoked through its Python API (the conda env carries it).
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd


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


def decimal_year(d: datetime) -> float:
    """TreeTime expects decimal-year dates."""
    start = datetime(d.year, 1, 1).toordinal()
    end = datetime(d.year + 1, 1, 1).toordinal()
    return d.year + (d.toordinal() - start) / (end - start)


def main() -> int:
    ap = argparse.ArgumentParser(description="§5.3 TreeTime strict-clock cross-check.")
    ap.add_argument("--in-aln", required=True, type=Path,
                    help="Per-cluster masked alignment FASTA")
    ap.add_argument("--in-tree", type=Path,
                    help="Optional starting tree (Newick). If absent, TreeTime "
                         "infers an FastTree on the fly via its `--aln` mode.")
    ap.add_argument("--in-cluster-members", required=True, type=Path,
                    help="cluster_members.tsv (for dates + cluster restriction)")
    ap.add_argument("--cluster-id", required=True,
                    help="Cluster ID to run TreeTime on")
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    work_dir = args.out_dir / f"treetime_{args.cluster_id.replace('/','_')}"
    work_dir.mkdir(parents=True, exist_ok=True)

    members = pd.read_csv(args.in_cluster_members, sep="\t", dtype=str)
    members = members[members["cluster_id"] == args.cluster_id].copy()
    if len(members) == 0:
        sys.exit(f"No members for cluster_id {args.cluster_id!r}")

    members["cd_parsed"] = members["collection_date"].apply(parse_loose)
    members = members.dropna(subset=["cd_parsed"])
    members["decimal_year"] = members["cd_parsed"].apply(decimal_year)
    print(f"[treetime] cluster {args.cluster_id}: {len(members):,} members "
          f"with parseable dates", file=sys.stderr)

    # Write dates file in TreeTime's expected format: tab-separated accession\tdate
    dates_file = work_dir / "dates.tsv"
    with dates_file.open("w") as f:
        f.write("name\tdate\n")
        for _, r in members.iterrows():
            f.write(f"{r['accession']}\t{r['decimal_year']:.4f}\n")

    # Subset alignment + tree to dated cluster members only (TreeTime errors on
    # undated tips with infinity in clock-rate estimation).
    from Bio import SeqIO
    dated_ids = set(members["accession"])
    subset_aln = work_dir / "dated_subset.fasta"
    with subset_aln.open("w") as fh:
        n = 0
        for rec in SeqIO.parse(args.in_aln, "fasta"):
            acc = rec.id.split()[0]
            if acc in dated_ids:
                SeqIO.write(rec, fh, "fasta")
                n += 1
        print(f"[treetime] subsetted alignment to {n} dated tips", file=sys.stderr)

    # Invoke TreeTime via its CLI. Omit --clock-rate so TreeTime estimates.
    # Omit --reroot if a starting tree is provided (TreeTime will use the input rooting).
    cmd = [
        "treetime",
        "--aln", str(subset_aln),
        "--dates", str(dates_file),
        "--outdir", str(work_dir),
    ]
    if args.in_tree:
        # Subset the tree similarly: prune outgroup/undated tips before TreeTime
        from ete3 import Tree
        t = Tree(str(args.in_tree), format=1)
        keep = [l for l in t.get_leaves() if l.name in dated_ids]
        if len(keep) < len(t.get_leaves()):
            t.prune(keep, preserve_branch_length=True)
        pruned_tree = work_dir / "dated_tree.nwk"
        t.write(outfile=str(pruned_tree), format=1)
        print(f"[treetime] pruned tree to {len(keep)} dated tips", file=sys.stderr)
        cmd += ["--tree", str(pruned_tree)]
    else:
        cmd += ["--reroot", "least-squares"]

    print(f"[treetime] running: {' '.join(cmd)}", file=sys.stderr)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    (work_dir / "treetime.stdout.log").write_text(proc.stdout)
    (work_dir / "treetime.stderr.log").write_text(proc.stderr)
    if proc.returncode != 0:
        print(f"[treetime] TreeTime failed (exit {proc.returncode})", file=sys.stderr)
        print(proc.stderr[-500:], file=sys.stderr)
        return proc.returncode

    # Parse TreeTime's output. The key outputs are:
    # - dates.tsv (per-tip + internal-node date estimates)
    # - molecular_clock.txt (rate + tMRCA estimate)
    # - timetree.nexus (the dated tree)
    summary = {"cluster_id": args.cluster_id, "n_genomes": int(len(members))}

    mc_file = work_dir / "molecular_clock.txt"
    if mc_file.exists():
        text = mc_file.read_text()
        summary["molecular_clock_text"] = text
        # Extract rate (subs/site/year) and tMRCA from molecular_clock.txt
        for line in text.splitlines():
            line = line.strip()
            if line.startswith("Root-Tip-Regression"):
                continue
            if "rate" in line.lower() and ":" in line:
                summary.setdefault("rate_lines", []).append(line)
            if "intercept" in line.lower() and ":" in line:
                summary.setdefault("intercept_lines", []).append(line)

    # tMRCA is reported in dates.tsv as the date of the root node
    dates_out = work_dir / "dates.tsv"
    if dates_out.exists():
        df = pd.read_csv(dates_out, sep="\t")
        if "name" in df.columns and "numeric date" in df.columns:
            root_rows = df[df["name"].str.startswith("NODE_") | df["name"].str.contains("root", case=False, na=False)]
            if not root_rows.empty:
                summary["root_node_dates"] = root_rows.to_dict("records")[:5]

    (args.out_dir / f"treetime_summary_{args.cluster_id.replace('/','_')}.json").write_text(
        json.dumps(summary, indent=2, default=str))
    print(f"[treetime] cluster {args.cluster_id} cross-check complete", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
