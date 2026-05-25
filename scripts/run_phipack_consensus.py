#!/usr/bin/env python3
"""
Per-candidate PhiPack runs for the §5.4 consensus filter (Amendment 03).

For each candidate sequence flagged by 3SEQ as a recombinant, run PhiPack on
the (candidate + parental panel + positive control) subset alignment. Capture
the three PhiPack statistics (Φ permutation, NSS, Max Chi²) and apply the
Amendment 03 internal-ensemble rule: ≥2 of 3 statistics rejecting the null at
P < 0.05 → the candidate is a "PhiPack call".

Inputs
------
- §5.4 masked alignment (FASTA)
- 3SEQ longRec file (one candidate accession per line, possibly with description)
- Parental panel TSVs (Ia, Ib, IIb)
- Positive control accession (default OZ375330.1)

Outputs
-------
- One JSON per candidate at out_dir/per_candidate/<accession>.json
- Aggregated summary at out_dir/phipack_consensus_summary.tsv:
    accession, n_subset_sequences, phi_perm_p, phi_normal_p, nss_p, maxchi_p,
    n_stats_significant_at_p05, phipack_call (True/False)
- Progress log at out_dir/phipack_run.log
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
from Bio import SeqIO


PVAL_PATTERNS = {
    "phi_perm_p":    re.compile(r"PHI\s*\(Permutation\):\s*([\d.eE+-]+)"),
    "phi_normal_p":  re.compile(r"PHI\s*\(Normal\):\s*([\d.eE+-]+)"),
    "nss_p":         re.compile(r"NSS:\s*([\d.eE+-]+)"),
    "maxchi_p":      re.compile(r"Max\s*Chi\^2:\s*([\d.eE+-]+)"),
}


def parse_phi_output(text: str) -> dict[str, float | None]:
    out: dict[str, float | None] = {k: None for k in PVAL_PATTERNS}
    for k, pat in PVAL_PATTERNS.items():
        m = pat.search(text)
        if m:
            try:
                out[k] = float(m.group(1))
            except ValueError:
                out[k] = None
    return out


def load_candidates(longrec_path: Path, limit: int | None = None) -> list[str]:
    """One accession per line — first whitespace-token."""
    accs = []
    with longrec_path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            acc = line.split()[0]
            if acc:
                accs.append(acc)
    return accs[:limit] if limit else accs


def load_panel_accessions(panel_paths: list[Path]) -> set[str]:
    panel: set[str] = set()
    for p in panel_paths:
        df = pd.read_csv(p, sep="\t", dtype=str)
        panel.update(df["accession"].dropna().tolist())
    return panel


def main() -> int:
    ap = argparse.ArgumentParser(description="§5.4 per-candidate PhiPack consensus filter.")
    ap.add_argument("--in-aln", required=True, type=Path,
                    help="§5.4 masked alignment FASTA (e.g. alignment_phase4.masked.fasta)")
    ap.add_argument("--in-longrec", required=True, type=Path,
                    help="3SEQ longRec file with candidate accessions")
    ap.add_argument("--panel-ia", required=True, type=Path)
    ap.add_argument("--panel-ib", required=True, type=Path)
    ap.add_argument("--panel-iib", required=True, type=Path)
    ap.add_argument("--positive-control", default="OZ375330.1")
    ap.add_argument("--phi-binary", default="Phi")
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--n-permutations", type=int, default=1000)
    ap.add_argument("--limit", type=int, default=0,
                    help="Stop after N candidates (0 = all)")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    per_cand_dir = args.out_dir / "per_candidate"
    per_cand_dir.mkdir(parents=True, exist_ok=True)
    log_path = args.out_dir / "phipack_run.log"
    log_fh = log_path.open("w")

    def log(msg: str):
        line = f"[{time.strftime('%H:%M:%S')}] {msg}"
        print(line, file=sys.stderr)
        log_fh.write(line + "\n")
        log_fh.flush()

    # --- Load inputs ---
    log(f"loading alignment {args.in_aln}")
    aln_records = {rec.id.split()[0]: rec for rec in SeqIO.parse(args.in_aln, "fasta")}
    log(f"  {len(aln_records):,} sequences in §5.4 alignment")

    panel = load_panel_accessions([args.panel_ia, args.panel_ib, args.panel_iib])
    log(f"  parental panel: {len(panel)} accessions (Ia + Ib + IIb)")

    cands = load_candidates(args.in_longrec, limit=args.limit if args.limit else None)
    log(f"  candidates to test: {len(cands)}")

    pos = args.positive_control

    # --- Per-candidate runs ---
    rows = []
    skipped = []
    t0 = time.time()
    for i, cand in enumerate(cands, 1):
        cand_json = per_cand_dir / f"{cand}.json"
        if cand_json.exists():
            try:
                rows.append(json.loads(cand_json.read_text()))
                continue
            except Exception:
                pass  # re-run

        if cand not in aln_records:
            log(f"[{i}/{len(cands)}] SKIP {cand} (not in §5.4 alignment)")
            skipped.append({"accession": cand, "reason": "missing_from_alignment"})
            continue

        subset_ids = {cand, pos} | panel
        subset_seqs = [aln_records[a] for a in subset_ids if a in aln_records]
        if len(subset_seqs) < 4:
            log(f"[{i}/{len(cands)}] SKIP {cand} (subset size {len(subset_seqs)} < 4)")
            skipped.append({"accession": cand, "reason": f"subset_too_small_{len(subset_seqs)}"})
            continue

        # Write subset alignment to a tmpfile
        subset_path = per_cand_dir / f"{cand}.subset.fasta"
        SeqIO.write(subset_seqs, subset_path, "fasta")

        # Run Phi
        cmd = [args.phi_binary, "-f", str(subset_path), "-p", str(args.n_permutations), "-o"]
        cwd = per_cand_dir / cand
        cwd.mkdir(parents=True, exist_ok=True)
        t_run = time.time()
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=600, cwd=str(cwd))
        except subprocess.TimeoutExpired:
            log(f"[{i}/{len(cands)}] TIMEOUT {cand} (>10 min)")
            skipped.append({"accession": cand, "reason": "phi_timeout"})
            continue

        pvals = parse_phi_output(r.stdout + "\n" + r.stderr)
        n_sig = sum(1 for k in ("phi_perm_p", "nss_p", "maxchi_p")
                     if pvals.get(k) is not None and pvals[k] < 0.05)
        phipack_call = n_sig >= 2

        rec = {
            "accession": cand,
            "n_subset_sequences": len(subset_seqs),
            "subset_includes_positive_control": pos in subset_ids,
            "phi_perm_p": pvals["phi_perm_p"],
            "phi_normal_p": pvals["phi_normal_p"],
            "nss_p": pvals["nss_p"],
            "maxchi_p": pvals["maxchi_p"],
            "n_stats_significant_at_p05": n_sig,
            "phipack_call": bool(phipack_call),
            "phi_runtime_sec": round(time.time() - t_run, 2),
        }
        cand_json.write_text(json.dumps(rec, indent=2))
        rows.append(rec)

        # Clean up the subset FASTA + per-candidate cwd (Phi temp files) to save disk
        subset_path.unlink(missing_ok=True)
        for f in cwd.glob("*"):
            f.unlink(missing_ok=True)
        cwd.rmdir()

        elapsed = time.time() - t0
        eta_min = (elapsed / i) * (len(cands) - i) / 60
        log(f"[{i}/{len(cands)}] {cand}: phi_perm={pvals['phi_perm_p']!r} "
            f"nss={pvals['nss_p']!r} maxchi={pvals['maxchi_p']!r} "
            f"call={phipack_call}  ({rec['phi_runtime_sec']}s; ETA {eta_min:.0f} min)")

    # --- Aggregate ---
    if rows:
        df = pd.DataFrame(rows)
        summary_tsv = args.out_dir / "phipack_consensus_summary.tsv"
        df.to_csv(summary_tsv, sep="\t", index=False)
        n_calls = int(df["phipack_call"].sum())
        log(f"DONE: {len(rows)} candidates tested, {n_calls} PhiPack calls (≥2/3 stats sig at P<0.05)")
        log(f"  summary: {summary_tsv}")

    if skipped:
        skip_tsv = args.out_dir / "phipack_skipped.tsv"
        pd.DataFrame(skipped).to_csv(skip_tsv, sep="\t", index=False)
        log(f"  skipped: {len(skipped)} -> {skip_tsv}")

    log(f"total wall-clock: {(time.time()-t0)/60:.1f} min")
    log_fh.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
