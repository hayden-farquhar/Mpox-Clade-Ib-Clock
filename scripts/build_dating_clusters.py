#!/usr/bin/env python3
"""
Define dating clusters per OSF Pre-Reg CASR2 §3.5 + §5.0 step 7.

Inputs
------
- Per-genome metadata + Nextclade clade/lineage assignments (clade-i Nextclade TSV)
- §3.3 corpus accession list (restricts to the 198 Ib survivors)
- Optional: recombinant-exclusion list (consensus candidates from §5.4) — these
  are removed from cluster membership per §5.0 step 6

Cluster eligibility (§3.5)
--------------------------
- min_genomes ≥ 10
- date span ≥ 6 months

Cluster definition (§5.3)
-------------------------
"Define clusters by a threshold on patristic distance, or use Nextstrain's
lineage assignments."

This script supports two modes:
- mode = "nextstrain_lineage": group by Nextclade lineage tag. Returns one
  cluster per unique non-null lineage; Ib genomes without finer lineage tags
  (the current Nextclade clade-i dataset has no Ib sub-lineages) collapse to
  a single "Ib-global" cluster.
- mode = "patristic": cut an IQ-TREE tree at a chosen height, produce
  sub-clusters from the resulting partition. Requires --in-tree.

Outputs
-------
- clusters_manifest.tsv — one row per eligible cluster: cluster_id, mode,
  n_genomes, earliest, latest, span_days, eligible (always True in output)
- cluster_members.tsv — accession → cluster_id mapping
- cluster_definition.json — mode + parameters used, freeze date, eligibility
  summary
"""

from __future__ import annotations

import argparse
import json
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


def define_by_lineage(meta: pd.DataFrame, nc: pd.DataFrame) -> pd.DataFrame:
    """One cluster per unique Nextclade lineage; NaN lineages collapse to 'Ib-global'."""
    # If meta already carries a 'lineage' column (e.g. when the merged metadata.tsv
    # from Phase 1 is passed in directly), use it; otherwise pull lineage from nc.
    if "lineage" in meta.columns:
        merged = meta[["accession", "lineage"]].copy()
    else:
        merged = meta.merge(nc[["accession", "lineage"]], on="accession", how="left")
    merged["cluster_id"] = merged["lineage"].fillna("Ib-global")
    return merged[["accession", "cluster_id"]]


def define_by_patristic(meta: pd.DataFrame, tree_path: Path,
                         threshold: float) -> pd.DataFrame:
    """
    Cut a Newick tree at the given height, produce sub-clusters.

    Implementation: load tree via ete3; for each tip, walk up until
    cumulative branch-length from tip exceeds `threshold` (or root reached);
    the node label at that point defines the cluster_id.
    """
    from ete3 import Tree
    tree = Tree(str(tree_path), format=1)
    cluster_for: dict[str, str] = {}
    for leaf in tree.get_leaves():
        cum = 0.0
        node = leaf
        while node.up is not None and cum + (node.dist or 0.0) < threshold:
            cum += node.dist or 0.0
            node = node.up
        cluster_for[leaf.name] = node.name or f"node_{id(node)}"
    out = pd.DataFrame({"accession": list(cluster_for.keys()),
                        "cluster_id": list(cluster_for.values())})
    return out[out["accession"].isin(meta["accession"])]


def main() -> int:
    ap = argparse.ArgumentParser(description="§3.5 dating-cluster definition.")
    ap.add_argument("--mode", choices=["nextstrain_lineage", "patristic"],
                    default="nextstrain_lineage")
    ap.add_argument("--in-metadata", required=True, type=Path,
                    help="Per-genome metadata with collection_date + clade columns")
    ap.add_argument("--in-nextclade-tsv", required=True, type=Path,
                    help="Clade-i Nextclade TSV (for lineage column)")
    ap.add_argument("--in-tree", type=Path,
                    help="IQ-TREE Newick (required for patristic mode)")
    ap.add_argument("--patristic-threshold", type=float, default=0.0005,
                    help="Cumulative branch-length threshold for patristic mode "
                         "(default 0.0005 substitutions/site)")
    ap.add_argument("--restrict-clade", default="Ib")
    ap.add_argument("--exclude-recombinants", type=Path,
                    help="Optional text file with one accession per line to "
                         "exclude from clustering (consensus recombinant candidates)")
    ap.add_argument("--min-genomes", type=int, default=10)
    ap.add_argument("--min-date-span-months", type=int, default=6)
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(args.in_metadata, sep="\t", dtype=str)
    # Standardise date column
    date_col = next((c for c in ("Isolate Collection date", "collection_date")
                     if c in meta.columns), None)
    if date_col is None:
        sys.exit(f"No date column in metadata; cols = {list(meta.columns)[:20]}")
    meta["collection_date"] = meta[date_col]

    # Restrict to clade
    if args.restrict_clade and "clade" in meta.columns:
        meta = meta[meta["clade"].fillna("").str.startswith(args.restrict_clade)]
    print(f"[input] {len(meta):,} {args.restrict_clade} genomes from metadata", file=sys.stderr)

    # Exclude recombinants if provided
    if args.exclude_recombinants and args.exclude_recombinants.exists():
        with args.exclude_recombinants.open() as fh:
            excl = {line.strip().split()[0] for line in fh if line.strip()}
        before = len(meta)
        meta = meta[~meta["accession"].isin(excl)]
        print(f"[exclude] {before - len(meta):,} recombinant accessions removed; "
              f"{len(meta):,} remain", file=sys.stderr)

    # Cluster definition
    nc = pd.read_csv(args.in_nextclade_tsv, sep="\t", dtype=str,
                     usecols=["seqName", "clade", "lineage"])
    nc["accession"] = nc["seqName"].str.split().str[0]

    if args.mode == "nextstrain_lineage":
        members = define_by_lineage(meta, nc)
    else:
        if not args.in_tree:
            sys.exit("--in-tree is required for patristic mode")
        members = define_by_patristic(meta, args.in_tree, args.patristic_threshold)

    # Attach dates + apply eligibility
    members = members.merge(meta[["accession", "collection_date"]], on="accession", how="left")
    members["cd_parsed"] = members["collection_date"].apply(parse_loose)

    rows = []
    eligible = []
    for cid, sub in members.groupby("cluster_id"):
        n = len(sub)
        dates = sub["cd_parsed"].dropna()
        if len(dates) < 2:
            early, late, span = "", "", 0
        else:
            early = dates.min().strftime("%Y-%m-%d")
            late = dates.max().strftime("%Y-%m-%d")
            span = (dates.max() - dates.min()).days
        elig = (n >= args.min_genomes) and (span >= args.min_date_span_months * 30)
        rows.append({"cluster_id": cid, "mode": args.mode, "n_genomes": n,
                     "earliest": early, "latest": late, "span_days": span,
                     "eligible": elig})
        if elig:
            eligible.append(cid)

    manifest = pd.DataFrame(rows)
    manifest.to_csv(args.out_dir / "clusters_manifest.tsv", sep="\t", index=False)
    print(f"[clusters] {len(rows)} candidate clusters; {len(eligible)} eligible "
          f"under §3.5 (n≥{args.min_genomes}, span≥{args.min_date_span_months} months)",
          file=sys.stderr)

    # Write per-genome cluster mapping for eligible clusters only
    members_elig = members[members["cluster_id"].isin(eligible)].copy()
    members_elig[["accession", "cluster_id", "collection_date"]].to_csv(
        args.out_dir / "cluster_members.tsv", sep="\t", index=False)
    print(f"[members] {len(members_elig):,} genomes assigned to eligible clusters",
          file=sys.stderr)

    # Summary JSON
    summary = {
        "mode": args.mode,
        "patristic_threshold": args.patristic_threshold if args.mode == "patristic" else None,
        "restrict_clade": args.restrict_clade,
        "min_genomes": args.min_genomes,
        "min_date_span_months": args.min_date_span_months,
        "n_input_genomes": int(len(meta)),
        "n_candidate_clusters": int(len(rows)),
        "n_eligible_clusters": int(len(eligible)),
        "eligible_cluster_ids": eligible,
    }
    (args.out_dir / "cluster_definition.json").write_text(json.dumps(summary, indent=2))

    return 0


if __name__ == "__main__":
    sys.exit(main())
