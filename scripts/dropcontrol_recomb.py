#!/usr/bin/env python3
"""Parent-panel sensitivity test for the P86 clade-Ib recombination scan.

Reviewer's last substantive question on paragraph [80]: the primary scan
reconstructs all 198 clade-Ib children from the same best parent pair, a
clade-Ia genome and the clade-Ib/IIb positive control OZ375330.1. The 100%
over-call is therefore conditioned on the known recombinant being in the
candidate-parent panel; it is the universal donor for the IIb-derived tracts
that diverged clade-Ib genomes share. Does the flag count survive removing the
control from the parental panel?

This re-runs 3SEQ on the SAME primary (nextstrain-masked) alignment used for the
baseline scan, with OZ375330.1 deleted from the alignment so it can no longer be
a candidate parent. The 198 clade-Ib children are still tested as recipients
(`-b1 -e198`); the 30 clade-IIb outgroup genomes remain available as alternative
inter-clade donors. Two informative outcomes:

  * flag rate stays ~198/198  -> the ~200-flag burden is intrinsic to clade-Ib
    divergence (3SEQ falls back to a IIb outgroup donor); the headline magnitude
    generalises to corpora lacking a known recombinant.
  * flag rate falls           -> a share of the over-call is induced by having
    the recombinant in the panel; a pre-discovery surveillance corpus would see
    fewer flags, and the magnitude is partly self-induced by the control.

Input is the deposited reordered baseline FASTA (198 Ib children first), so the
child block is identical to the baseline run; only OZ375330.1 (row 910, in the
non-child remainder) is dropped. Heavy IO is staged to /tmp per portfolio policy
(Documents/ is iCloud-synced). Resumable: a completed rec.csv is parsed, not
recomputed.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

PROJ = Path(__file__).resolve().parents[1]
BASELINE_FASTA = PROJ / "outputs/apobec3_mask_experiment/aln_baseline.fasta"
META = PROJ / "data/processed/freeze_20260522/metadata.tsv"
OUTDIR = PROJ / "outputs/apobec3_mask_experiment/run_dropcontrol"

TOOLS = Path.home() / "tools"
THREESEQ = TOOLS / "3seq/3seq"
PTABLE = TOOLS / "3seq/3seq_ptable_300"

CONTROL = "OZ375330.1"   # positive control / universal donor to drop from parents
N_IB = 198               # clade-Ib children (rows 1..198 of the reordered FASTA)
RUN_ID = "dropcontrol"
STAGE = Path("/tmp") / "p86_dropcontrol"


def log(msg: str) -> None:
    print(f"[dropcontrol] {msg}", flush=True)


def accession_of(header: str) -> str:
    return header.split()[0]


def read_fasta(path: Path):
    """Return list of (accession, full_header, sequence) in file order."""
    recs, header, cur = [], None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    recs.append((accession_of(header), header, "".join(cur)))
                header = line[1:].strip()
                cur = []
            else:
                cur.append(line.strip())
    if header is not None:
        recs.append((accession_of(header), header, "".join(cur)))
    return recs


def load_ib_set():
    ib = set()
    with open(META) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        ci = head.index("clade")
        ai = head.index("accession")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > max(ci, ai) and parts[ci] == "Ib":
                ib.add(parts[ai])
    return ib


def parse_rec(rec: Path, ib_set):
    """{child: min corrected p} for distinct clade-Ib children, plus parent pairs."""
    minp, parents = {}, {}
    if not rec or not rec.exists():
        return minp, parents
    with open(rec) as fh:
        fh.readline()
        for line in fh:
            parts = line.rstrip("\n").split(",")
            if len(parts) < 11:
                continue
            p_name, q_name, child = parts[0].strip(), parts[1].strip(), parts[2].strip()
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
            if child not in minp or p < minp[child]:
                minp[child] = p
                parents[child] = (p_name, q_name)
    return minp, parents


def log_complete(workdir: Path) -> bool:
    lg = workdir / f"{RUN_ID}.3s.log"
    if not lg.exists():
        return False
    txt = lg.read_text()
    return "100.000%" in txt or "Number of triples tested" in txt


def main():
    for b in (THREESEQ, PTABLE, BASELINE_FASTA, META):
        if not b.exists():
            raise SystemExit(f"missing input: {b}")
    OUTDIR.mkdir(parents=True, exist_ok=True)
    STAGE.mkdir(parents=True, exist_ok=True)

    ib_set = load_ib_set()
    log(f"clade-Ib accessions: {len(ib_set)}")

    recs = read_fasta(BASELINE_FASTA)
    log(f"baseline alignment: {len(recs)} records")
    first198 = [a for a, _, _ in recs[:N_IB]]
    n_ib_in_block = sum(1 for a in first198 if a in ib_set)
    if n_ib_in_block != N_IB:
        raise SystemExit(
            f"child block integrity failed: only {n_ib_in_block}/{N_IB} of the "
            f"first {N_IB} records are clade Ib; aborting before 3SEQ")
    log(f"child-block check: {n_ib_in_block}/{N_IB} first records are clade Ib (OK)")

    kept = [r for r in recs if r[0] != CONTROL]
    dropped = len(recs) - len(kept)
    if dropped != 1:
        raise SystemExit(f"expected to drop exactly 1 record ({CONTROL}); dropped {dropped}")
    if any(a == CONTROL for a, _, _ in kept):
        raise SystemExit(f"{CONTROL} still present after filtering")
    log(f"dropped {CONTROL} from candidate parents; {len(kept)} records remain")

    staged_fasta = STAGE / "aln_dropcontrol.fasta"
    with open(staged_fasta, "w") as fh:
        for acc, _, seq in kept:
            fh.write(f">{acc}\n{seq}\n")
    log(f"staged child-first alignment to {staged_fasta}")

    rec = OUTDIR / f"{RUN_ID}.3s.rec.csv"
    if rec.exists() and log_complete(OUTDIR):
        log("resume: completed rec.csv present, skipping 3SEQ run")
    else:
        # stage the run in /tmp, then copy artefacts back into the repo
        cmd = [str(THREESEQ), "-full", str(staged_fasta), "-ptable", str(PTABLE),
               "-id", RUN_ID, f"-b1", f"-e{N_IB}"]
        log(f"running 3SEQ children 1-{N_IB} (no {CONTROL} parent); ~1-2 h ...")
        subprocess.run(cmd, cwd=STAGE, capture_output=True, text=True, input="Y\n")
        for suf in ("3s.rec.csv", "3s.log", "3s.longRec", "3s.pvalHist"):
            src = STAGE / f"{RUN_ID}.{suf}"
            if src.exists():
                shutil.copy2(src, OUTDIR / src.name)
        if not rec.exists():
            alt = STAGE / "3s.rec.csv"
            if alt.exists():
                shutil.copy2(alt, rec)

    minp, parents = parse_rec(rec, ib_set)
    n_flagged = len(minp)

    # classify the new donor (the non-Ia parent) for each child
    ia_set = set()
    with open(META) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        ci, ai = head.index("clade"), head.index("accession")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) > max(ci, ai) and p[ci] == "Ia":
                ia_set.add(p[ai])

    pair_counts = Counter()
    donor_counts = Counter()
    for child, (pn, qn) in parents.items():
        pair_counts[f"{pn} + {qn}"] += 1
        # the donor is whichever parent is NOT clade Ia
        donor = pn if pn not in ia_set else qn
        if pn not in ia_set and qn not in ia_set:
            donor = f"{pn}|{qn}(neither-Ia)"
        donor_counts[donor] += 1

    summary = dict(
        test="drop positive control OZ375330.1 from candidate parents",
        baseline_flagged_ib=198,
        baseline_parent="all 198 = clade-Ia parent + OZ375330.1 (positive control)",
        dropcontrol_flagged_ib=n_flagged,
        n_ib_children=N_IB,
        control_present_as_parent=False,
        top_parent_pairs=pair_counts.most_common(10),
        donor_composition=donor_counts.most_common(10),
    )
    (OUTDIR / "dropcontrol_summary.json").write_text(json.dumps(summary, indent=2))
    log("SUMMARY:\n" + json.dumps(summary, indent=2))
    log(f"RESULT: {n_flagged}/{N_IB} clade-Ib children still flagged without the control as parent")


if __name__ == "__main__":
    sys.exit(main())
