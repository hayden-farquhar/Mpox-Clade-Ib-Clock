#!/usr/bin/env python3
"""
APOBEC3-context masking experiment for the P86 clade-Ib recombination scan.

Decisive test for the reviewer's collinearity critique: the published enrichment
test (flagged-Ib APOBEC3 load vs unflagged-Ia load) cannot isolate APOBEC3 as the
*cause* of the 3SEQ flag, because flag status is 100% collinear with clade. This
script manufactures the missing within-clade contrast.

Three 3SEQ conditions, identical except for which columns are masked:
  baseline  : reordered alignment, no extra masking (must reproduce ~198/198 Ib)
  apobec3   : APOBEC3-context homoplasy columns (TC/GA with a cognate edit) -> N
  random    : the SAME NUMBER of variable non-APOBEC3 columns -> N (seeded control),
              replicated over ten seeds (anchor + RANDOM_REPLICATE_SEEDS) to report
              the control as a distribution (Table S7) rather than a single draw

Children are restricted to the 198 clade-Ib sequences (-b1 -e198); all 926
sequences remain as parents. Measure: number of distinct Ib accessions appearing
as the recombinant child (C_name) in each condition's rec.csv. If the Ib flag
rate collapses under `apobec3` but survives under `random`, causation is pinned on
APOBEC3 homoplasy directly.

Also computes the within-clade-Ib Spearman correlation between each Ib child's
minimum corrected 3SEQ p-value (baseline) and its APOBEC3 fraction -- the
alternative within-clade test the reviewer proposed.

Resumable: a condition whose rec.csv already exists with a completed log is skipped.
"""
from __future__ import annotations

import csv
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

# ---------------------------------------------------------------------------
# Paths (no hardcoded /Users absolute paths: resolve relative to this file)
# ---------------------------------------------------------------------------
PROJ = Path(__file__).resolve().parents[1]
ALN = PROJ / "data/interim/freeze_20260522/alignment_phase4.masked.fasta"
META = PROJ / "data/processed/freeze_20260522/metadata.tsv"
A3COUNTS = PROJ / "data/processed/freeze_20260522/apobec3_counts_v2.tsv"
OUTDIR = PROJ / "outputs/apobec3_mask_experiment"

TOOLS = Path.home() / "tools"
THREESEQ = TOOLS / "3seq/3seq"
PTABLE = TOOLS / "3seq/3seq_ptable_300"

REF_ACC = "DQ011155.1"        # reference row (first record)
N_IB = 198                    # clade-Ib children
SEED = 20260608

CONDITIONS = ("baseline", "apobec3", "random")

# The single `random` condition above (seed=SEED) is the anchor of a replicate
# distribution: these additional integer seeds each draw an independent set of k
# variable non-APOBEC3 columns, turning the one-shot control into the Table S7
# distribution. Runs are resumable (a seed whose completed rec.csv already exists
# is parsed, not recomputed), so the deposited run_random_seed* outputs reproduce
# the table exactly on this machine; on a fresh machine the loop reproduces the
# methodology (the 3SEQ triplet heuristic is the same; per-seed counts are draw-
# and platform-dependent, which is itself the point the distribution makes).
RANDOM_REPLICATE_SEEDS = (1, 2, 3, 4, 5, 6, 7, 8, 9)


# ---------------------------------------------------------------------------
# IO helpers
# ---------------------------------------------------------------------------
def log(msg: str) -> None:
    print(f"[apobec3_mask] {msg}", flush=True)


def accession_of(header: str) -> str:
    return header.split()[0]


def read_alignment(path: Path):
    """Return (accessions list in file order, numpy uint8 array shape (n, L))."""
    accs, rows, cur = [], [], []
    L = None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    s = "".join(cur).encode()
                    if L is None:
                        L = len(s)
                    elif len(s) != L:
                        raise SystemExit("alignment is not rectangular")
                    rows.append(np.frombuffer(s, dtype=np.uint8))
                    cur = []
                accs.append(accession_of(line[1:].strip()))
            else:
                cur.append(line.strip())
    if cur:
        s = "".join(cur).encode()
        rows.append(np.frombuffer(s, dtype=np.uint8))
    arr = np.vstack(rows)
    return accs, arr


def load_clade_sets():
    ib, ia = set(), set()
    with open(META) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            c = row.get("clade", "")
            if c == "Ib":
                ib.add(row["accession"])
            elif c == "Ia":
                ia.add(row["accession"])
    return ib, ia


def load_a3_fraction():
    frac = {}
    with open(A3COUNTS) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            try:
                frac[row["accession"]] = float(row["apobec3_fraction"])
            except (KeyError, ValueError):
                pass
    return frac


# ---------------------------------------------------------------------------
# Column selection
# ---------------------------------------------------------------------------
ORD = {b: ord(b) for b in "ACGTN-"}


def apobec3_homoplasy_columns(arr: np.ndarray, ref_idx: int):
    """Columns where the reference is in TC or GA context AND at least one
    sequence carries the cognate APOBEC3 edit (C->T for TC, G->A for GA)."""
    ref = arr[ref_idx]
    L = ref.shape[0]
    prev_ref = np.empty_like(ref); prev_ref[1:] = ref[:-1]; prev_ref[0] = 0
    next_ref = np.empty_like(ref); next_ref[:-1] = ref[1:]; next_ref[-1] = 0

    tc_ctx = (ref == ORD["C"]) & (prev_ref == ORD["T"])   # edited C is at column i
    ga_ctx = (ref == ORD["G"]) & (next_ref == ORD["A"])   # edited G is at column i

    cognate_T = (arr == ORD["T"]).any(axis=0)             # some seq has T at column
    cognate_A = (arr == ORD["A"]).any(axis=0)

    tc_cols = tc_ctx & cognate_T
    ga_cols = ga_ctx & cognate_A
    cols = np.where(tc_cols | ga_cols)[0]
    return cols, int(tc_cols.sum()), int(ga_cols.sum())


def variable_columns(arr: np.ndarray):
    """Columns with >=2 distinct unambiguous ACGT bases across all sequences."""
    distinct = np.zeros(arr.shape[1], dtype=np.int8)
    for b in "ACGT":
        distinct += (arr == ORD[b]).any(axis=0).astype(np.int8)
    return np.where(distinct >= 2)[0]


# ---------------------------------------------------------------------------
# FASTA writing (accession-only headers => clean rec.csv parsing)
# ---------------------------------------------------------------------------
def write_condition_fasta(path: Path, accs, arr, order_idx, mask_cols):
    work = arr
    if mask_cols is not None and len(mask_cols):
        work = arr.copy()
        work[:, mask_cols] = ORD["N"]
    with open(path, "w") as fh:
        for i in order_idx:
            fh.write(f">{accs[i]}\n")
            fh.write(work[i].tobytes().decode("ascii"))
            fh.write("\n")


# ---------------------------------------------------------------------------
# 3SEQ run + parse
# ---------------------------------------------------------------------------
def log_complete(workdir: Path, run_id: str) -> bool:
    lg = workdir / f"{run_id}.3s.log"
    if not lg.exists():
        return False
    txt = lg.read_text()
    return "100.000%" in txt or "Number of triples tested" in txt


def run_3seq_children(fasta: Path, workdir: Path, run_id: str, lo: int, hi: int):
    workdir.mkdir(parents=True, exist_ok=True)
    rec = workdir / f"{run_id}.3s.rec.csv"
    if rec.exists() and log_complete(workdir, run_id):
        log(f"  [{run_id}] resume: completed rec.csv present, skipping run")
        return rec
    # this 3SEQ build uses 1-based, ATTACHED -b/-e flags restricting the CHILD set
    cmd = [str(THREESEQ), "-full", str(fasta), "-ptable", str(PTABLE),
           "-id", run_id, f"-b{lo}", f"-e{hi}"]
    log(f"  [{run_id}] running 3SEQ children {lo}-{hi} (this takes ~1-2 h) ...")
    subprocess.run(cmd, cwd=workdir, capture_output=True, text=True, input="Y\n")
    if not rec.exists():
        alt = workdir / "3s.rec.csv"
        if alt.exists():
            rec = alt
    return rec


def parse_rec(rec: Path, ib_set):
    """Return {accession: min_corrected_p} for distinct Ib children (C_name).

    Headers were written as bare accessions, so each line splits cleanly on
    commas: fields[2] = C_name; fields[9], fields[10] = the two Dunn-Sidak p's.
    The trailing breakpoints field may itself contain commas, but we only read
    fixed low-index fields, so that is irrelevant.
    """
    out = {}
    if not rec or not rec.exists():
        return out
    with open(rec) as fh:
        header = fh.readline()  # discard
        for line in fh:
            parts = line.rstrip("\n").split(",")
            if len(parts) < 11:
                continue
            child = parts[2].strip()
            if child not in ib_set:
                continue
            ps = []
            for j in (9, 10):
                try:
                    v = float(parts[j])
                    if 0 < v <= 1:
                        ps.append(v)
                except ValueError:
                    pass
            if not ps:
                continue
            p = min(ps)
            out[child] = min(out[child], p) if child in out else p
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    for b in (THREESEQ, PTABLE):
        if not b.exists():
            raise SystemExit(f"missing dependency: {b}")
    OUTDIR.mkdir(parents=True, exist_ok=True)

    log("loading alignment ...")
    accs, arr = read_alignment(ALN)
    n, L = arr.shape
    log(f"  {n} sequences x {L} columns")
    if REF_ACC not in accs:
        raise SystemExit(f"{REF_ACC} not found in alignment")
    ref_idx = accs.index(REF_ACC)

    ib_set, ia_set = load_clade_sets()
    log(f"  clade Ib={len(ib_set)}  Ia={len(ia_set)}")
    a3frac = load_a3_fraction()

    # order: 198 clade-Ib children FIRST (file order preserved within group), rest after
    ib_idx = [i for i, a in enumerate(accs) if a in ib_set]
    rest_idx = [i for i, a in enumerate(accs) if a not in ib_set]
    if len(ib_idx) != N_IB:
        log(f"  WARNING: expected {N_IB} Ib, found {len(ib_idx)}")
    order_idx = ib_idx + rest_idx
    ib_children = [accs[i] for i in ib_idx]

    # column sets
    log("identifying APOBEC3-context homoplasy columns ...")
    a3_cols, n_tc, n_ga = apobec3_homoplasy_columns(arr, ref_idx)
    log(f"  APOBEC3 homoplasy columns: {len(a3_cols)} (TC={n_tc}, GA={n_ga})")

    log("identifying variable non-APOBEC3 columns for matched random control ...")
    var_cols = variable_columns(arr)
    a3_set = set(a3_cols.tolist())
    candidate = np.array([c for c in var_cols if c not in a3_set])
    rng = np.random.default_rng(SEED)
    k = min(len(a3_cols), len(candidate))
    rand_cols = np.sort(rng.choice(candidate, size=k, replace=False))
    log(f"  variable cols={len(var_cols)}  candidate (non-A3)={len(candidate)}  "
        f"sampled random control={len(rand_cols)} (seed={SEED})")

    mask_for = {"baseline": None, "apobec3": a3_cols, "random": rand_cols}

    # write + run each condition
    results = {}
    for cond in CONDITIONS:
        fasta = OUTDIR / f"aln_{cond}.fasta"
        if not fasta.exists():
            log(f"writing {fasta.name} ...")
            write_condition_fasta(fasta, accs, arr, order_idx, mask_for[cond])
        wd = OUTDIR / f"run_{cond}"
        rec = run_3seq_children(fasta, wd, cond, 1, N_IB)
        child_p = parse_rec(rec, ib_set)
        results[cond] = child_p
        log(f"  [{cond}] distinct Ib children flagged: {len(child_p)} / {N_IB}")
        # incremental save
        (OUTDIR / "results_partial.json").write_text(json.dumps(
            {c: {"n_flagged": len(results[c])} for c in results}, indent=2))

    # ---- matched-random control replicate distribution (Table S7) ----
    log("running matched-random control across replicate seeds ...")
    rand_counts = {"anchor_20260608": len(results.get("random", {}))}
    for s in RANDOM_REPLICATE_SEEDS:
        rng_s = np.random.default_rng(int(s))
        cols_s = np.sort(rng_s.choice(candidate, size=k, replace=False))
        fasta_s = OUTDIR / f"aln_random_seed{s}.fasta"
        if not fasta_s.exists():
            log(f"writing {fasta_s.name} ...")
            write_condition_fasta(fasta_s, accs, arr, order_idx, cols_s)
        wd_s = OUTDIR / f"run_random_seed{s}"
        rec_s = run_3seq_children(fasta_s, wd_s, f"random_seed{s}", 1, N_IB)
        rand_counts[f"seed{s}"] = len(parse_rec(rec_s, ib_set))
        log(f"  [random seed {s}] flagged: {rand_counts[f'seed{s}']} / {N_IB}")
    rvals = [int(v) for v in rand_counts.values()]
    rand_dist = dict(
        per_seed=rand_counts, n=len(rvals),
        min=min(rvals), max=max(rvals),
        mean=round(float(np.mean(rvals)), 1),
        median=float(np.median(rvals)),
        sd=round(float(np.std(rvals, ddof=1)), 1),
    )
    (OUTDIR / "random_distribution_results.json").write_text(
        json.dumps(rand_dist, indent=2))
    log("random-control distribution: " + json.dumps(rand_dist))

    # within-clade Spearman on baseline
    base = results.get("baseline", {})
    pairs = [(a3frac[a], base[a]) for a in ib_children if a in base and a in a3frac]
    if len(pairs) >= 3:
        x = [p[0] for p in pairs]
        y = [p[1] for p in pairs]
        rho, pval = spearmanr(x, y)
        within = dict(n=len(pairs), spearman_rho=float(rho), spearman_p=float(pval),
                      note="x=apobec3_fraction, y=min corrected 3SEQ p (baseline); "
                           "negative rho => higher APOBEC3 load tracks stronger "
                           "(smaller-p) recombination signal within clade Ib")
    else:
        within = dict(n=len(pairs), spearman_rho=None, spearman_p=None,
                      note="insufficient paired data")

    summary = dict(
        n_sequences=int(n), n_columns=int(L),
        n_ib_children=len(ib_children),
        apobec3_columns=dict(total=int(len(a3_cols)), tc=n_tc, ga=n_ga),
        random_control_columns=int(len(rand_cols)),
        seed=SEED,
        flagged_ib=dict(
            baseline=len(results.get("baseline", {})),
            apobec3=len(results.get("apobec3", {})),
            random=len(results.get("random", {})),
        ),
        within_clade_correlation=within,
        random_control_distribution=rand_dist,
    )
    (OUTDIR / "results_summary.json").write_text(json.dumps(summary, indent=2))

    # tidy TSV of per-child min-p across conditions
    with open(OUTDIR / "per_ib_child_minp.tsv", "w") as fh:
        fh.write("accession\tapobec3_fraction\tbaseline_minp\tapobec3_minp\trandom_minp\n")
        for a in ib_children:
            fh.write("\t".join(str(v) for v in (
                a,
                a3frac.get(a, ""),
                results.get("baseline", {}).get(a, ""),
                results.get("apobec3", {}).get(a, ""),
                results.get("random", {}).get(a, ""),
            )) + "\n")

    log("DONE")
    log(json.dumps(summary, indent=2))


if __name__ == "__main__":
    sys.exit(main())
