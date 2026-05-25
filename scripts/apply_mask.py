#!/usr/bin/env python3
"""
Apply a BED file of masked intervals to every sequence in a multiple-sequence
alignment (FASTA). Masked positions become N in every record.

BED format: 0-based half-open intervals (start, end). One interval per line.
Lines beginning with '#' or 'track' are skipped.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import MutableSeq
from tqdm import tqdm


def parse_bed(bed_path: Path) -> list[tuple[int, int]]:
    intervals: list[tuple[int, int]] = []
    with bed_path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue  # malformed line
            try:
                start, end = int(parts[1]), int(parts[2])
            except ValueError:
                continue  # column-header row or other non-numeric line
            if end > start:
                intervals.append((start, end))
    return intervals


def main() -> int:
    ap = argparse.ArgumentParser(description="Apply a BED mask to a FASTA alignment.")
    ap.add_argument("--in-aln", required=True, type=Path)
    ap.add_argument("--bed", required=True, type=Path)
    ap.add_argument("--out-aln", required=True, type=Path)
    args = ap.parse_args()

    intervals = parse_bed(args.bed)
    total_mask_bp = sum(e - s for s, e in intervals)
    print(f"[mask] {len(intervals):,} intervals, {total_mask_bp:,} bp total to mask", file=sys.stderr)

    args.out_aln.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    aln_length: int | None = None
    with args.out_aln.open("w") as out_fh:
        for rec in tqdm(SeqIO.parse(args.in_aln, "fasta"), desc="masking"):
            seq = MutableSeq(str(rec.seq))
            if aln_length is None:
                aln_length = len(seq)
            for start, end in intervals:
                lo = max(0, start)
                hi = min(len(seq), end)
                if hi > lo:
                    seq[lo:hi] = "N" * (hi - lo)
            rec.seq = seq.toseq() if hasattr(seq, "toseq") else seq
            SeqIO.write(rec, out_fh, "fasta")
            written += 1

    print(
        f"[mask] wrote {written:,} sequences (alignment length {aln_length:,} bp; "
        f"mask covers ~{total_mask_bp/aln_length:.1%} of each genome)",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
