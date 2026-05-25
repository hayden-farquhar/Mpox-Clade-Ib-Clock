#!/usr/bin/env python3
"""
Fallback mask builder: produces a BED file covering only the inverted terminal
repeats (ITRs) of the MPXV reference. Used when no Nextstrain mpox masking
BED is available locally and the user hasn't supplied one via MASK_BED.

Heuristic: find the longest perfect reverse-complement match between the 5'
end and the 3' end of the reference (the canonical poxvirus ITR signature).
Reports its coordinates as two BED intervals (one per ITR).

This is a deliberately conservative fallback — it masks only ITRs and does
NOT mask hyper-variable regions. For a publication-grade analysis, supply
the curated mask BED from the Nextstrain mpox build via MASK_BED in
workflow/config.yaml.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from Bio import SeqIO

# ITR sizes for MPXV cluster around 6-10 kb. Cap the search to avoid spurious
# long internal matches. Minimum reportable ITR is 1 kb.
MAX_ITR_BP = 12_000
MIN_ITR_BP = 1_000


def reverse_complement(s: str) -> str:
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return s.translate(table)[::-1]


def find_itr_length(seq: str, max_bp: int) -> int:
    """
    Find the longest k such that seq[:k] == reverse_complement(seq[-k:]).
    Linear scan from large k downward with early stop.
    """
    n = len(seq)
    upper = min(max_bp, n // 2)
    # Sweep large-to-small, return at first match (which is the longest).
    # For speed, compare slices in chunks: this is O(n^2) worst case but
    # for n ~ 200 kb and upper ~ 12 kb, ~2.4 GB comparisons worst case —
    # still seconds in Python because string comparison is C-level.
    seq = seq.upper()
    head = seq[:upper]
    tail_rc = reverse_complement(seq[-upper:])
    # Find the longest matching prefix between head and tail_rc.
    k_max = 0
    for k in range(upper, MIN_ITR_BP - 1, -1):
        if head[:k] == tail_rc[:k]:
            k_max = k
            break
    return k_max


def main() -> int:
    ap = argparse.ArgumentParser(description="Build a fallback ITR-only mask BED.")
    ap.add_argument("--reference", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--max-itr-bp", type=int, default=MAX_ITR_BP)
    args = ap.parse_args()

    rec = next(SeqIO.parse(args.reference, "fasta"))
    seq = str(rec.seq).upper()
    n = len(seq)

    itr_k = find_itr_length(seq, args.max_itr_bp)
    if itr_k == 0:
        print(
            f"[mask] no ITR ≥{MIN_ITR_BP} bp detected on {rec.id} "
            f"(genome length {n:,}); writing EMPTY mask BED",
            file=sys.stderr,
        )
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text("# fallback mask: no ITR detected; nothing masked\n")
        return 0

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        fh.write("# Fallback mask: ITRs only, auto-derived from the reference.\n")
        fh.write("# For publication-grade analyses, supply a curated mask BED.\n")
        fh.write(f"{rec.id}\t0\t{itr_k}\tITR_5prime\n")
        fh.write(f"{rec.id}\t{n - itr_k}\t{n}\tITR_3prime\n")

    print(
        f"[mask] {rec.id}: detected {itr_k:,} bp ITR; "
        f"masked 0–{itr_k:,} and {n - itr_k:,}–{n:,} (genome length {n:,})",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
