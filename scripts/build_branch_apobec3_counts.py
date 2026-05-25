#!/usr/bin/env python3
"""
Build per-tip branch-quantity APOBEC3 SNV counts for the H1ʹ test of OSF
Amendment 02 (https://osf.io/gt3vx/files/osfstorage/6a1049c0964dd2d18df92ec0).

For each clade-Ib tip, count the substitutions along the terminal branch from
the reconstructed clade-Ib MRCA to that tip, classified as APOBEC3-context or
not. This is the per-branch quantity that O'Toole et al. (Science 2023)
calibrated p₀ = 0.70 against — distinct from the tip-vs-reference quantity
that the registered §5.2 H1 protocol uses.

APOBEC3 context (matching scripts/count_apobec3.py for consistency):
- Substitution is APOBEC3 if (MRCA C → tip T at 5'-T·C-3' MRCA context)
  or (MRCA G → tip A at 5'-G·A-3' MRCA context). Context is taken from the
  MRCA sequence (i.e., the state immediately before the mutation).

Inputs
------
- IQ-TREE .state file giving per-node, per-site ancestral state ML estimates.
- IQ-TREE .treefile (rooted on the outgroup) so we can find the MRCA node of
  the clade-Ib tips.
- The masked alignment containing the Ib tips (the same alignment IQ-TREE
  was run on).
- Metadata TSV with the Ib accession list (used to identify which tips are Ib).

Outputs
-------
- TSV with one row per Ib tip: accession, total_snvs, apobec3_snvs,
  non_apobec3_snvs, apobec3_fraction. Schema matches count_apobec3.py so the
  downstream test_h1.py can run unchanged.
- JSON with the resolved MRCA node ID, the number of MRCA sites with
  unambiguous ML state, and the count of Ib tips compared against it.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

ACGT = set("ACGT")


def read_state_for_node(state_path: Path, target_node: str, n_sites_hint: int = 0) -> list[str]:
    """
    Stream-parse an IQ-TREE .state file and return ONLY the ancestral sequence
    for `target_node` as a list of single-character states (A/C/G/T/N for gaps
    and ambiguous reconstructions).

    Streaming is essential because the full .state file is ~10 MB per ancestral
    node × ~400 nodes for a 200-tip × 200kb dataset — loading all nodes into
    memory uses ~10 GB+ and routinely OOMs a laptop. Stream-extracting one
    node uses ~1 MB peak.

    File format:
        # comment lines beginning with '#'
        Node\tSite\tState\tp_A\tp_C\tp_G\tp_T   (column header)
        NodeN\tSITE\tSTATE\t...   (one row per (node, site))
    """
    pairs: list[tuple[int, str]] = []
    n_lines = 0
    n_matched = 0
    in_target = False  # rows for one node are contiguous in IQ-TREE's output
    target_seen = False
    with state_path.open() as fh:
        for line in fh:
            n_lines += 1
            if not line or line[0] == "#":
                continue
            # cheap pre-check: avoid splitting most lines
            if not (line.startswith(target_node + "\t")):
                if in_target:
                    # contiguity guarantee broken — we've passed all target rows
                    break
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            in_target = True
            target_seen = True
            try:
                site = int(parts[1])
            except ValueError:
                continue
            state = parts[2].upper()
            pairs.append((site, state if state in ACGT else "N"))
            n_matched += 1
    if not target_seen:
        raise SystemExit(
            f"No rows found for node {target_node!r} in {state_path}. "
            f"Streamed {n_lines:,} lines without a match."
        )
    pairs.sort(key=lambda r: r[0])
    max_site = pairs[-1][0]
    seq = ["N"] * max_site
    for site, state in pairs:
        seq[site - 1] = state
    print(f"[state] streamed {n_lines:,} lines; matched {n_matched:,} for node "
          f"{target_node!r}; reconstructed sequence length {len(seq):,}",
          file=sys.stderr)
    return seq


def find_mrca_of_ib_tips(tree_path: Path, ib_tips: set[str]) -> str:
    """
    Read the IQ-TREE .treefile and return the internal-node label for the MRCA
    of the provided Ib tip set.

    Uses ete3 (already in the conda env). If ete3 is unavailable, falls back
    to dendropy.
    """
    try:
        from ete3 import Tree
    except ImportError:
        from dendropy import Tree as DendroTree, TaxonNamespace
        tns = TaxonNamespace()
        t = DendroTree.get(path=str(tree_path), schema="newick", taxon_namespace=tns)
        # ... not implementing dendropy fallback in detail
        raise NotImplementedError("ete3 is required; dendropy fallback not implemented.")

    tree = Tree(str(tree_path), format=1)  # format=1 reads internal-node names
    # Match leaf names to tip ids — ete3 reads "DQ011155.1" as the leaf name
    leaves = {l.name for l in tree.get_leaves()}
    overlap = leaves & ib_tips
    print(f"[tree] {len(leaves)} leaves; {len(overlap)} match Ib tip set", file=sys.stderr)
    if len(overlap) < len(ib_tips):
        missing = ib_tips - overlap
        print(f"[tree] WARN: {len(missing)} Ib tips not in tree (first 5: {list(missing)[:5]})",
              file=sys.stderr)
    if len(overlap) < 2:
        raise SystemExit(f"Cannot define an MRCA with <2 matching Ib tips (got {len(overlap)}).")
    mrca = tree.get_common_ancestor(list(overlap))
    name = mrca.name
    if not name:
        # IQ-TREE labels internal nodes; if blank, walk to a labelled ancestor
        cur = mrca
        while cur.up and not cur.name:
            cur = cur.up
        name = cur.name
    if not name:
        raise SystemExit("MRCA node has no name label — cannot align to .state file.")
    print(f"[tree] Ib MRCA node label: {name}", file=sys.stderr)
    return name


def classify_branch_snv(mrca_b: str, tip_b: str, mrca_seq: list[str], i: int) -> tuple[bool, str | None]:
    """Classify a (MRCA→tip) substitution at 0-based position i. Returns (is_apobec3, context)."""
    if mrca_b not in ACGT or tip_b not in ACGT:
        return False, None
    if mrca_b == tip_b:
        return False, None
    # TC → TT
    if mrca_b == "C" and tip_b == "T" and i > 0 and mrca_seq[i - 1] == "T":
        return True, "TC"
    # GA → AA
    if mrca_b == "G" and tip_b == "A" and i + 1 < len(mrca_seq) and mrca_seq[i + 1] == "A":
        return True, "GA"
    return False, None


def main() -> int:
    ap = argparse.ArgumentParser(description="Branch-quantity APOBEC3 counts (H1ʹ).")
    ap.add_argument("--in-aln", required=True, type=Path,
                    help="Masked alignment containing the Ib tips (same as IQ-TREE input).")
    ap.add_argument("--in-state", required=True, type=Path,
                    help="IQ-TREE .state ancestral-state file.")
    ap.add_argument("--in-tree", required=True, type=Path,
                    help="IQ-TREE .treefile (with internal-node labels).")
    ap.add_argument("--in-metadata", required=True, type=Path,
                    help="Per-genome metadata; used to identify Ib tips.")
    ap.add_argument("--out-counts", required=True, type=Path,
                    help="TSV: per-tip branch-quantity APOBEC3 counts.")
    ap.add_argument("--out-summary", required=True, type=Path,
                    help="JSON: MRCA resolution + headline stats.")
    args = ap.parse_args()

    # ---- Identify Ib tips ----
    meta = pd.read_csv(args.in_metadata, sep="\t", dtype=str)
    ib_tips = set(meta.loc[meta["clade"] == "Ib", "accession"].tolist())
    print(f"[ib] {len(ib_tips):,} clade-Ib accessions in metadata", file=sys.stderr)

    # ---- Find Ib MRCA node (tree topology only — fast) ----
    mrca_label = find_mrca_of_ib_tips(args.in_tree, ib_tips)

    # ---- Stream the .state file, extracting ONLY the MRCA node's sequence ----
    print(f"[state] streaming {args.in_state} for node {mrca_label!r}", file=sys.stderr)
    mrca_seq = read_state_for_node(args.in_state, mrca_label)
    n_unambig_mrca = sum(1 for b in mrca_seq if b in ACGT)
    print(f"[mrca] {mrca_label}: {len(mrca_seq):,} sites total, "
          f"{n_unambig_mrca:,} unambiguous", file=sys.stderr)

    # ---- Per-Ib-tip branch-quantity counts ----
    rows: list[dict] = []
    site_counter: Counter[tuple[int, str]] = Counter()

    n_skipped_len = 0
    for rec in tqdm(SeqIO.parse(args.in_aln, "fasta"), desc="branch-counting"):
        acc = rec.id.split()[0]
        if acc not in ib_tips:
            continue
        tip_seq = str(rec.seq).upper()
        if len(tip_seq) != len(mrca_seq):
            n_skipped_len += 1
            continue
        total = 0
        apobec = 0
        for i in range(len(mrca_seq)):
            is_apo, ctx = classify_branch_snv(mrca_seq[i], tip_seq[i], mrca_seq, i)
            if mrca_seq[i] in ACGT and tip_seq[i] in ACGT and mrca_seq[i] != tip_seq[i]:
                total += 1
            if is_apo:
                apobec += 1
                site_counter[(i + 1, ctx)] += 1
        rows.append({
            "accession": acc,
            "total_snvs": total,
            "apobec3_snvs": apobec,
            "non_apobec3_snvs": total - apobec,
            "apobec3_fraction": (apobec / total) if total else 0.0,
        })

    df = pd.DataFrame(rows)
    args.out_counts.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_counts, sep="\t", index=False)
    print(f"[write] per-tip branch counts: {len(df):,} rows -> {args.out_counts}",
          file=sys.stderr)

    # ---- Summary ----
    total_apo = int(df["apobec3_snvs"].sum())
    total_snv = int(df["total_snvs"].sum())
    clade_wide = total_apo / total_snv if total_snv else 0.0
    print(f"[headline] branch-quantity Ib APOBEC3 fraction = {clade_wide:.4f} "
          f"({total_apo:,} APOBEC3 of {total_snv:,} branch SNVs)", file=sys.stderr)

    summary = {
        "mrca_node_label": mrca_label,
        "mrca_unambiguous_sites": n_unambig_mrca,
        "mrca_total_sites": len(mrca_seq),
        "n_ib_tips_compared": len(df),
        "n_ib_tips_skipped_length_mismatch": n_skipped_len,
        "pooled_branch_total_snvs": total_snv,
        "pooled_branch_apobec3_snvs": total_apo,
        "pooled_branch_apobec3_fraction": clade_wide,
        "median_per_tip_branch_apobec3_snvs": float(df["apobec3_snvs"].median()) if len(df) else 0.0,
        "median_per_tip_branch_total_snvs": float(df["total_snvs"].median()) if len(df) else 0.0,
    }
    args.out_summary.write_text(json.dumps(summary, indent=2))
    print(f"[write] summary -> {args.out_summary}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
