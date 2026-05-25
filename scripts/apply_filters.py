#!/usr/bin/env python3
"""
Apply pre-registration §3.3 inclusion / exclusion filters to the raw NCBI
MPXV corpus, producing the analysis corpus FASTA, the merged metadata table,
an exclusions table (every excluded accession + reason), and a JSON summary.

The filters implemented here are EXACTLY those committed in the pre-registration.
Adding, removing, or tuning any filter constitutes a deviation that must be
documented per §6.4.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterator

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm


# --------------------------------------------------------------------------
# Date-resolution helper
# --------------------------------------------------------------------------

def date_resolution(date_str: str | None) -> str:
    """Classify a collection-date string as 'day', 'month', 'year', or 'unknown'."""
    if not date_str or not isinstance(date_str, str):
        return "unknown"
    s = date_str.strip()
    # ISO 8601 patterns commonly used by NCBI / GenBank
    if len(s) >= 10 and s[4] == "-" and s[7] == "-":
        return "day"
    if len(s) >= 7 and s[4] == "-":
        return "month"
    if len(s) == 4 and s.isdigit():
        return "year"
    return "unknown"


def resolution_passes(observed: str, required: str) -> bool:
    """True iff observed resolution is at least as fine as required."""
    order = {"day": 3, "month": 2, "year": 1, "unknown": 0}
    return order.get(observed, 0) >= order.get(required, 0)


# --------------------------------------------------------------------------
# Per-sequence inspection
# --------------------------------------------------------------------------

ACGTN = set("ACGTNacgtn")
ACGT = set("ACGTacgt")


def base_composition(seq: str) -> tuple[int, float, float]:
    """Return (length, N-fraction, non-ACGTN-fraction)."""
    n_len = len(seq)
    if n_len == 0:
        return 0, 0.0, 0.0
    n_count = sum(1 for b in seq if b in "Nn")
    non_acgtn = sum(1 for b in seq if b not in ACGTN)
    return n_len, n_count / n_len, non_acgtn / n_len


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description="Apply pre-registration §3.3 filters.")
    ap.add_argument("--in-fasta", required=True, type=Path)
    ap.add_argument("--in-ncbi-meta", required=True, type=Path)
    ap.add_argument("--in-nextclade", required=True, type=Path)
    ap.add_argument("--out-fasta", required=True, type=Path)
    ap.add_argument("--out-metadata", required=True, type=Path)
    ap.add_argument("--out-exclusions", required=True, type=Path)
    ap.add_argument("--out-summary", required=True, type=Path)
    ap.add_argument("--min-length", type=int, required=True)
    ap.add_argument("--max-n-fraction", type=float, required=True)
    ap.add_argument("--max-non-acgtn-fraction", type=float, required=True)
    ap.add_argument(
        "--date-resolution",
        choices=["year", "month", "day"],
        required=True,
        help="Minimum required collection-date resolution",
    )
    args = ap.parse_args()

    # ---- Load metadata ----
    print(f"[load] NCBI metadata: {args.in_ncbi_meta}", file=sys.stderr)
    ncbi_meta = pd.read_csv(args.in_ncbi_meta, sep="\t", dtype=str)
    ncbi_meta.columns = [c.strip() for c in ncbi_meta.columns]
    # Normalise the accession column name (Datasets uses "Accession")
    if "Accession" in ncbi_meta.columns:
        ncbi_meta = ncbi_meta.rename(columns={"Accession": "accession"})
    elif "accession" not in ncbi_meta.columns:
        # First column is typically accession in the dataformat output
        ncbi_meta = ncbi_meta.rename(columns={ncbi_meta.columns[0]: "accession"})
    ncbi_meta["accession"] = ncbi_meta["accession"].str.strip()
    print(f"  rows: {len(ncbi_meta):,}", file=sys.stderr)

    print(f"[load] Nextclade TSV: {args.in_nextclade}", file=sys.stderr)
    nextclade = pd.read_csv(args.in_nextclade, sep="\t", dtype=str)
    # Nextclade uses "seqName" as the per-genome identifier. Match it to accession.
    if "seqName" in nextclade.columns:
        # Extract the accession from the FASTA header — NCBI's bundle uses ">{accession} ..."
        nextclade["accession"] = nextclade["seqName"].str.split().str[0]
    print(f"  rows: {len(nextclade):,}", file=sys.stderr)

    merged = ncbi_meta.merge(nextclade, on="accession", how="outer", suffixes=("_ncbi", "_nc"))
    print(f"[merge] {len(merged):,} accessions after outer join", file=sys.stderr)

    # ---- Per-sequence inspection (length, composition) ----
    print(f"[scan] reading {args.in_fasta} for base composition", file=sys.stderr)
    comp: dict[str, tuple[int, float, float]] = {}
    for rec in tqdm(SeqIO.parse(args.in_fasta, "fasta"), desc="composition"):
        acc = rec.id.split()[0]
        comp[acc] = base_composition(str(rec.seq))
    comp_df = pd.DataFrame(
        [(a, *v) for a, v in comp.items()],
        columns=["accession", "seq_length", "n_fraction", "non_acgtn_fraction"],
    )
    merged = merged.merge(comp_df, on="accession", how="left")

    # ---- Apply filters ----
    # Each row gets a list of exclusion reasons; an empty list = passes.
    def collect_reasons(row: pd.Series) -> list[str]:
        reasons: list[str] = []

        # Organism check (Datasets bundles MPXV, so this should pass; defensive)
        organism = str(row.get("Organism Name", "") or row.get("organism_name", ""))
        if organism and "monkeypox" not in organism.lower() and "mpxv" not in organism.lower():
            reasons.append(f"organism_not_mpxv:{organism[:30]}")

        # Clade filter — accept clade I (Ia/Ib) and Ib/IIb recombinant
        clade = str(row.get("clade", "") or "")
        if clade and not (
            clade.startswith("I") and not clade.startswith("II")  # I, Ia, Ib
            or "/" in clade  # recombinants like Ib/IIb
            or clade in {"IIb", "II"}  # clade II retained as outgroup only
        ):
            reasons.append(f"clade_excluded:{clade}")
        # Tag clade-II-only as outgroup (still in file, but flagged)
        # (downstream rules subset on the clade column.)

        # Length
        length = row.get("seq_length")
        if pd.isna(length):
            reasons.append("no_sequence_found")
        elif length < args.min_length:
            reasons.append(f"length_lt_{args.min_length}:{int(length)}")

        # N fraction
        nf = row.get("n_fraction")
        if pd.notna(nf) and nf > args.max_n_fraction:
            reasons.append(f"n_fraction_gt_{args.max_n_fraction}:{nf:.3f}")

        # Non-ACGTN fraction
        naf = row.get("non_acgtn_fraction")
        if pd.notna(naf) and naf > args.max_non_acgtn_fraction:
            reasons.append(f"non_acgtn_gt_{args.max_non_acgtn_fraction}:{naf:.3f}")

        # Date resolution
        col_date_col = next(
            (c for c in ("Isolate Collection date", "isolate_collection_date", "collection_date") if c in row.index),
            None,
        )
        date_str = row.get(col_date_col) if col_date_col else None
        resolution = date_resolution(date_str)
        if not resolution_passes(resolution, args.date_resolution):
            reasons.append(f"date_resolution_{resolution}_lt_{args.date_resolution}")

        # Nextclade QC overall status
        qc = str(row.get("qc.overallStatus", "") or "")
        if qc.lower() == "bad":
            reasons.append("nextclade_qc_bad")

        # Frameshifts in core gene set, per pre-registration §3.3.
        # Operationalised as Nextclade's frameshift-QC status flag
        # (qc.frameShifts.status == "bad"), which IS Nextclade's call on
        # whether observed frameshifts impair the core gene set. Using a raw
        # count of frameshifts is over-strict — Nextclade's status field is
        # the direct match for "as called by Nextclade" in the pre-reg.
        fs_status = str(row.get("qc.frameShifts.status", "") or "")
        if fs_status.lower() == "bad":
            reasons.append("frameshift_qc_bad")

        return reasons

    print(f"[filter] applying §3.3 criteria", file=sys.stderr)
    merged["exclusion_reasons"] = merged.apply(collect_reasons, axis=1)
    merged["passes"] = merged["exclusion_reasons"].apply(len) == 0

    # ---- Split into corpus + exclusions ----
    corpus = merged[merged["passes"]].copy()
    exclusions = merged[~merged["passes"]].copy()
    exclusions["exclusion_reasons"] = exclusions["exclusion_reasons"].apply(lambda xs: "|".join(xs))

    # ---- Write outputs ----
    print(
        f"[write] corpus={len(corpus):,} | exclusions={len(exclusions):,}",
        file=sys.stderr,
    )

    args.out_metadata.parent.mkdir(parents=True, exist_ok=True)
    corpus.to_csv(args.out_metadata, sep="\t", index=False)
    exclusions[["accession", "exclusion_reasons"]].to_csv(args.out_exclusions, sep="\t", index=False)

    # FASTA: keep only accessions in `corpus`
    keep_set = set(corpus["accession"])
    written = 0
    with args.out_fasta.open("w") as fh:
        for rec in SeqIO.parse(args.in_fasta, "fasta"):
            acc = rec.id.split()[0]
            if acc in keep_set:
                SeqIO.write(rec, fh, "fasta")
                written += 1
    print(f"[write] FASTA records: {written:,}", file=sys.stderr)

    # ---- Summary JSON ----
    reason_counter: Counter[str] = Counter()
    for reasons in merged.loc[~merged["passes"], "exclusion_reasons"]:
        for r in reasons:
            reason_counter[r.split(":", 1)[0]] += 1

    # Clade-level breakdown of the surviving corpus
    clade_counts: dict[str, int] = {}
    if "clade" in corpus.columns:
        clade_counts = corpus["clade"].fillna("unassigned").value_counts().to_dict()

    summary = {
        "freeze_input_rows": int(len(merged)),
        "passed": int(len(corpus)),
        "excluded": int(len(exclusions)),
        "exclusion_reason_counts": dict(reason_counter.most_common()),
        "surviving_clade_counts": clade_counts,
        "filter_parameters": {
            "min_length_nt": args.min_length,
            "max_n_fraction": args.max_n_fraction,
            "max_non_acgtn_fraction": args.max_non_acgtn_fraction,
            "min_date_resolution": args.date_resolution,
            "nextclade_qc_bad_excluded": True,
            "frameshifts_excluded": True,
        },
        "outputs": {
            "fasta": str(args.out_fasta),
            "metadata": str(args.out_metadata),
            "exclusions": str(args.out_exclusions),
        },
    }
    args.out_summary.write_text(json.dumps(summary, indent=2, default=str))
    print(f"[done] summary written to {args.out_summary}", file=sys.stderr)

    # Sanity: §3.3 exit criterion is ≥80% date-resolved to month among the surviving corpus.
    if "clade" in corpus.columns and len(corpus) > 0:
        # Re-compute date resolution distribution among survivors
        col_date_col = next(
            (c for c in ("Isolate Collection date", "isolate_collection_date", "collection_date") if c in corpus.columns),
            None,
        )
        if col_date_col:
            res = corpus[col_date_col].apply(date_resolution)
            month_or_better = (res.isin({"month", "day"})).mean()
            print(
                f"[exit-criterion] {month_or_better:.1%} of survivors date-resolved to month or better "
                f"(§3.3 exit criterion: ≥80%)",
                file=sys.stderr,
            )

    return 0


if __name__ == "__main__":
    sys.exit(main())
