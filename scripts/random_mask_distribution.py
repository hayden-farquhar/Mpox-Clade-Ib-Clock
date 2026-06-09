#!/usr/bin/env python3
"""Multi-seed distribution for the matched random-column masking control (P86).

The primary masking experiment (apobec3_mask_recomb.py) runs one random-column
control under a single seed (20260608 -> 114/198 clade-Ib children retained). A
single draw cannot show whether that count is stable or a lucky pick of which
845 variable non-APOBEC3 columns happen to be masked. This driver repeats the
random condition across several seeds and reports the distribution of retained
flag counts, so the control can be stated as a mean and range rather than a
point.

It reuses the column-selection and 3SEQ-run helpers from apobec3_mask_recomb so
the per-seed procedure is byte-identical to the primary experiment, differing
only in the RNG seed. The canonical seed 20260608 reuses the already-completed
run_random directory (no re-run).

Robustness: resumable (a seed whose rec.csv exists with a completed log is
skipped), incremental JSON save after every seed, no hardcoded absolute paths.

Each 3SEQ run is ~1.3 h, so a full sweep is multi-hour; run in the background.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from apobec3_mask_recomb import (  # noqa: E402
    N_IB,
    OUTDIR,
    PTABLE,
    REF_ACC,
    THREESEQ,
    apobec3_homoplasy_columns,
    log,
    load_clade_sets,
    parse_rec,
    read_alignment,
    run_3seq_children,
    variable_columns,
    write_condition_fasta,
)
from apobec3_mask_recomb import ALN  # noqa: E402

# Seed 20260608 is the canonical primary-experiment seed; it reuses run_random.
CANONICAL_SEED = 20260608
EXTRA_SEEDS = [1, 2, 3, 4, 5, 6, 7, 8, 9]
ALL_SEEDS = [CANONICAL_SEED] + EXTRA_SEEDS

DIST_JSON = OUTDIR / "random_distribution.json"


def load_existing() -> dict:
    if DIST_JSON.exists():
        return json.loads(DIST_JSON.read_text())
    return {"seeds": {}, "note": "per-seed clade-Ib children retained under "
            "matched random-column masking (845 cols -> N)"}


def summarise(counts: list[int]) -> dict:
    a = np.array(counts, dtype=float)
    return dict(
        n_seeds=len(counts),
        counts=sorted(int(c) for c in counts),
        mean=float(a.mean()),
        sd=float(a.std(ddof=1)) if len(a) > 1 else 0.0,
        min=int(a.min()),
        max=int(a.max()),
        median=float(np.median(a)),
    )


def main() -> int:
    for b in (THREESEQ, PTABLE):
        if not b.exists():
            raise SystemExit(f"missing dependency: {b}")

    log("loading alignment ...")
    accs, arr = read_alignment(ALN)
    ref_idx = accs.index(REF_ACC)
    ib_set, _ = load_clade_sets()

    ib_idx = [i for i, a in enumerate(accs) if a in ib_set]
    rest_idx = [i for i, a in enumerate(accs) if a not in ib_set]
    order_idx = ib_idx + rest_idx

    a3_cols, _, _ = apobec3_homoplasy_columns(arr, ref_idx)
    var_cols = variable_columns(arr)
    a3_set = set(a3_cols.tolist())
    candidate = np.array([c for c in var_cols if c not in a3_set])
    k = min(len(a3_cols), len(candidate))
    log(f"  random-control pool: {len(candidate)} non-APOBEC3 variable cols; "
        f"sampling k={k} per seed")

    state = load_existing()

    for seed in ALL_SEEDS:
        skey = str(seed)
        if skey in state["seeds"] and state["seeds"][skey].get("n_flagged") is not None:
            log(f"  [seed {seed}] resume: already recorded "
                f"({state['seeds'][skey]['n_flagged']}/{N_IB}), skipping")
            continue

        if seed == CANONICAL_SEED:
            # reuse the primary experiment's completed run_random
            wd = OUTDIR / "run_random"
            rec = wd / "random.3s.rec.csv"
            if not rec.exists():
                rng = np.random.default_rng(seed)
                rand_cols = np.sort(rng.choice(candidate, size=k, replace=False))
                fasta = OUTDIR / "aln_random.fasta"
                if not fasta.exists():
                    write_condition_fasta(fasta, accs, arr, order_idx, rand_cols)
                rec = run_3seq_children(fasta, wd, "random", 1, N_IB)
        else:
            rng = np.random.default_rng(seed)
            rand_cols = np.sort(rng.choice(candidate, size=k, replace=False))
            fasta = OUTDIR / f"aln_random_seed{seed}.fasta"
            if not fasta.exists():
                log(f"  [seed {seed}] writing masked alignment ...")
                write_condition_fasta(fasta, accs, arr, order_idx, rand_cols)
            wd = OUTDIR / f"run_random_seed{seed}"
            rec = run_3seq_children(fasta, wd, f"random_seed{seed}", 1, N_IB)

        child_p = parse_rec(rec, ib_set)
        n_flagged = len(child_p)
        state["seeds"][skey] = {"n_flagged": n_flagged}
        log(f"  [seed {seed}] clade-Ib children flagged: {n_flagged}/{N_IB}")

        counts = [v["n_flagged"] for v in state["seeds"].values()
                  if v.get("n_flagged") is not None]
        state["summary"] = summarise(counts)
        DIST_JSON.write_text(json.dumps(state, indent=2))

    log("DONE")
    log(json.dumps(state.get("summary", {}), indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
