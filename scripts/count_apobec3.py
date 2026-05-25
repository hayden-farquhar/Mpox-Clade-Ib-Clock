#!/usr/bin/env python3
"""
Count APOBEC3-context SNVs per genome relative to the clade-I reference.

Per pre-registration §4.2:
- APOBEC3 context flag is TRUE if (C→T at 5'-TC context) OR (G→A at 5'-GA context).
- TC context: reference base at position i is C and the base at position i-1 is T.
- GA context: reference base at position i is G and the base at position i+1 is A.
  (These are the same physical edits on opposite strands.)

Per pre-registration §5.3, the eligible-target-sites count L is the number of
positions on the masked reference where a TC or GA dinucleotide context is
present and both bases are unambiguous and unmasked. L is written separately
as JSON for the dating step.

Outputs
-------
- TSV with one row per genome: accession, total_snvs, apobec3_snvs, apobec3_fraction
- BED with one row per APOBEC3-edited site: chrom, start (0-based), end, "TC" or "GA"
- JSON with L = count of eligible APOBEC3 target sites on the masked reference
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

ACGT = set("ACGT")


def eligible_apobec3_sites(reference: str) -> list[tuple[int, str]]:
    """Find every (position, context) where reference has TC or GA dinucleotide.

    Returns list of (1-based position of the C or G, 'TC' or 'GA').
    Positions inside masked regions (ref base = N) or with ambiguous flanks
    are excluded.
    """
    sites: list[tuple[int, str]] = []
    n = len(reference)
    for i, b in enumerate(reference):
        if b not in ACGT:
            continue
        # TC context: position i = C, position i-1 = T
        if b == "C" and i > 0 and reference[i - 1] == "T":
            sites.append((i + 1, "TC"))  # report as 1-based
        # GA context: position i = G, position i+1 = A
        if b == "G" and i + 1 < n and reference[i + 1] == "A":
            sites.append((i + 1, "GA"))
    return sites


def classify_snv(ref: str, alt: str, ref_seq: str, i: int) -> tuple[bool, str | None]:
    """
    Classify a (ref, alt) substitution at 0-based position i in ref_seq.
    Returns (is_apobec3, context_label_or_None).
    """
    if ref not in ACGT or alt not in ACGT:
        return False, None
    if ref == alt:
        return False, None
    # TC → TT
    if ref == "C" and alt == "T" and i > 0 and ref_seq[i - 1] == "T":
        return True, "TC"
    # GA → AA (the GA→AA equivalent of TC→TT on the opposite strand)
    if ref == "G" and alt == "A" and i + 1 < len(ref_seq) and ref_seq[i + 1] == "A":
        return True, "GA"
    return False, None


def main() -> int:
    ap = argparse.ArgumentParser(description="Count APOBEC3-context SNVs per genome.")
    ap.add_argument("--in-aln", required=True, type=Path)
    ap.add_argument("--in-reference", required=True, type=Path,
                    help="Reference FASTA (used only to validate identity of the first MSA record).")
    ap.add_argument("--in-metadata", required=True, type=Path)
    ap.add_argument("--out-counts", required=True, type=Path)
    ap.add_argument("--out-sites", required=True, type=Path)
    ap.add_argument("--out-L", required=True, type=Path)
    args = ap.parse_args()

    # ---- Load reference + alignment ----
    print(f"[load] reference: {args.in_reference}", file=sys.stderr)
    ref_records = list(SeqIO.parse(args.in_reference, "fasta"))
    if not ref_records:
        sys.exit("No records in reference FASTA.")
    ref_id = ref_records[0].id

    print(f"[load] alignment: {args.in_aln}", file=sys.stderr)
    aln_records = list(SeqIO.parse(args.in_aln, "fasta"))
    if not aln_records:
        sys.exit("No records in alignment.")

    # The first record in the alignment should be the reference (rule align_mafft
    # puts the reference at the top via cat <ref> <corpus> | mafft --reorder).
    # mafft --reorder may move it; find it by ID.
    aln_ref_idx = next(
        (i for i, r in enumerate(aln_records) if r.id.startswith(ref_id)),
        None,
    )
    if aln_ref_idx is None:
        sys.exit(f"Reference {ref_id} not found in alignment.")
    ref_aln = str(aln_records[aln_ref_idx].seq).upper()
    print(f"[ref] {ref_id}: {len(ref_aln):,} aligned columns "
          f"({sum(1 for b in ref_aln if b != '-'):,} non-gap)", file=sys.stderr)

    # ---- L: count eligible APOBEC3 target sites on the masked reference ----
    eligible = eligible_apobec3_sites(ref_aln)
    n_tc = sum(1 for _, c in eligible if c == "TC")
    n_ga = sum(1 for _, c in eligible if c == "GA")
    L = len(eligible)
    print(f"[L] eligible APOBEC3 target sites on masked reference: {L:,} "
          f"({n_tc:,} TC + {n_ga:,} GA)", file=sys.stderr)
    args.out_L.write_text(json.dumps({
        "L_total": L,
        "L_tc": n_tc,
        "L_ga": n_ga,
        "reference_alignment_length": len(ref_aln),
        "reference_id": ref_id,
    }, indent=2))

    # ---- Per-genome SNV counting ----
    rows: list[dict] = []
    site_counter: Counter[tuple[int, str]] = Counter()
    for idx, rec in enumerate(tqdm(aln_records, desc="counting")):
        if idx == aln_ref_idx:
            continue  # skip reference
        seq = str(rec.seq).upper()
        if len(seq) != len(ref_aln):
            print(f"[warn] {rec.id}: length mismatch ({len(seq)} vs {len(ref_aln)}); skipping",
                  file=sys.stderr)
            continue
        total = 0
        apobec = 0
        for i in range(len(ref_aln)):
            ref_b = ref_aln[i]
            alt_b = seq[i]
            if ref_b not in ACGT or alt_b not in ACGT:
                continue
            if ref_b == alt_b:
                continue
            total += 1
            is_apobec, ctx = classify_snv(ref_b, alt_b, ref_aln, i)
            if is_apobec:
                apobec += 1
                site_counter[(i + 1, ctx)] += 1
        accession = rec.id.split()[0]
        rows.append({
            "accession": accession,
            "total_snvs": total,
            "apobec3_snvs": apobec,
            "non_apobec3_snvs": total - apobec,
            "apobec3_fraction": (apobec / total) if total else 0.0,
        })

    # ---- Write outputs ----
    df = pd.DataFrame(rows)
    args.out_counts.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_counts, sep="\t", index=False)
    print(f"[write] per-genome counts: {len(df):,} rows -> {args.out_counts}", file=sys.stderr)

    # APOBEC3 sites BED
    with args.out_sites.open("w") as fh:
        fh.write("# chrom\tstart\tend\tcontext\toccupancy\n")
        for (pos, ctx), n in sorted(site_counter.items()):
            fh.write(f"{ref_id}\t{pos-1}\t{pos}\t{ctx}\t{n}\n")
    print(f"[write] APOBEC3 sites: {len(site_counter):,} unique sites edited", file=sys.stderr)

    # ---- Headline stats ----
    total_apo = int(df["apobec3_snvs"].sum())
    total_snv = int(df["total_snvs"].sum())
    clade_wide = total_apo / total_snv if total_snv else 0.0
    print(
        f"[headline] clade-wide APOBEC3 fraction = {clade_wide:.4f} "
        f"({total_apo:,} APOBEC3 of {total_snv:,} total SNVs)",
        file=sys.stderr,
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
