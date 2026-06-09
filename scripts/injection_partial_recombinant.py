"""Reviewer point #3: would the two-detector rule discard a GENUINE inter-clade
recombinant if it arrived on a partial (truncated) genome?

We take the one confirmed inter-clade recombinant (OZ375330.1, the UK Ib/IIb
positive control) inside the frozen GARD-subset alignment, progressively truncate
it to gaps from both ends (simulating an incomplete assembly that captured only a
central window), and at each truncation level re-run BOTH detectors against the
frozen parental panel:

  - 3SEQ (-full, single child) under the ptable used in the main scan
  - PhiPack (Phi -p 1000 -o): PHI(perm), NSS, Max Chi^2; "PhiPack call" = >=2/3 at P<0.05

The point of confirmation is the regime where 3SEQ still flags the recombinant but
PhiPack drops below the 2/3 rule: there, the consensus protocol DISCARDS a real
inter-clade event purely because the genome is partial. Read-only on the freeze;
writes one TSV + JSON to outputs/tables/.
"""
from pathlib import Path
import json
import re
import shutil
import subprocess
import tempfile

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
ALN = ROOT / "data/interim/freeze_20260522/alignment_phase4_gard_subset.masked.fasta"
MANIFEST = ROOT / "data/processed/freeze_20260522/panels/gard_subset_manifest.tsv"
OUT = ROOT / "outputs/tables"
OUT.mkdir(parents=True, exist_ok=True)

THREESEQ = Path.home() / "tools/3seq/3seq"
PTABLE = Path.home() / "tools/3seq/3seq_ptable_300"
PHI = Path.home() / "tools/PhiPack/Phi"

CONTROL = "OZ375330.1"
N_PER_CLADE = 10           # parental representatives per clade (~31-seq subset)
RETAIN = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.15, 0.1]
ALPHA = 0.05

PHI_PERM = re.compile(r"PHI\s*\(Permutation\):\s*([\d.eE+\-]+|--)")
PHI_NORM = re.compile(r"PHI\s*\(Normal\):\s*([\d.eE+\-]+|--)")
NSS_P = re.compile(r"NSS:\s*([\d.eE+\-]+|--)")
MAXCHI_P = re.compile(r"Max\s*Chi\^2:\s*([\d.eE+\-]+|--)")


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def write_fasta(path, seqs):
    with open(path, "w") as fh:
        for k, v in seqs.items():
            fh.write(f">{k}\n")
            for i in range(0, len(v), 80):
                fh.write(v[i:i + 80] + "\n")


def truncate_central(seq, retain):
    """Keep a central retain-fraction window; mask the flanks to gaps."""
    L = len(seq)
    keep = int(round(retain * L))
    if keep >= L:
        return seq
    start = (L - keep) // 2
    end = start + keep
    return "-" * start + seq[start:end] + "-" * (L - end)


def parse_float(tok):
    if tok is None or tok == "--":
        return None
    try:
        return float(tok)
    except ValueError:
        return None


def run_phi(subset_path, workdir):
    res = subprocess.run([str(PHI), "-f", str(subset_path), "-p", "1000", "-o"],
                         cwd=workdir, capture_output=True, text=True)
    out = res.stdout + "\n" + res.stderr
    def grab(rx):
        m = rx.search(out)
        return parse_float(m.group(1)) if m else None
    return dict(phi_perm_p=grab(PHI_PERM), phi_norm_p=grab(PHI_NORM),
                nss_p=grab(NSS_P), maxchi_p=grab(MAXCHI_P), raw=out)


def run_3seq(subset_path, child_idx, workdir, run_id):
    # this 3SEQ build uses 1-based, ATTACHED -b/-e flags (-b1 = first sequence)
    subprocess.run([str(THREESEQ), "-full", str(subset_path),
                    "-ptable", str(PTABLE), "-id", run_id,
                    f"-b{child_idx}", f"-e{child_idx}"],
                   cwd=workdir, capture_output=True, text=True, input="Y\n")
    rec = Path(workdir) / f"{run_id}.3s.rec.csv"
    if not rec.exists():
        rec = Path(workdir) / "3s.rec.csv"
    if not rec.exists():
        return dict(flagged=False, min_corr_p=None, n_triples=0)
    txt = rec.read_text()
    lines = [ln for ln in txt.splitlines() if CONTROL in ln]
    if not lines:
        return dict(flagged=False, min_corr_p=None, n_triples=0)
    # corrected (Dunn-Sidak) p-value is the last numeric-looking field that is a
    # probability; scan all fields and take the smallest value in (0,1].
    min_p = None
    for ln in lines:
        for tok in re.split(r"[,\t]", ln):
            v = parse_float(tok.strip())
            if v is not None and 0 < v <= 1:
                min_p = v if min_p is None else min(min_p, v)
    return dict(flagged=True, min_corr_p=min_p, n_triples=len(lines))


def main():
    for binp in (THREESEQ, PTABLE, PHI):
        if not binp.exists():
            raise SystemExit(f"missing dependency: {binp}")
    aln = read_fasta(ALN)
    man = pd.read_csv(MANIFEST, sep="\t")
    panel = {}
    for clade, role in [("Ia", "ia_representative"), ("Ib", "ib_representative"),
                        ("IIb", "outgroup_iib")]:
        accs = [a for a in man.loc[man.role == role, "accession"] if a in aln][:N_PER_CLADE]
        for a in accs:
            panel[a] = aln[a]
    assert CONTROL in aln, f"{CONTROL} not in alignment"
    L = len(aln[CONTROL])

    rows = []
    for retain in RETAIN:
        with tempfile.TemporaryDirectory() as td:
            trunc = truncate_central(aln[CONTROL], retain)
            retained_bp = len(trunc.replace("-", "").replace("N", "").replace("n", ""))
            # subset: truncated control FIRST (child index 1) + parental panel
            subset = {CONTROL: trunc}
            subset.update(panel)
            sp = Path(td) / "subset.fasta"
            write_fasta(sp, subset)
            t3 = run_3seq(sp, 1, td, f"inj_{int(retain*100)}")
            phi = run_phi(sp, td)
            calls = [phi["phi_perm_p"], phi["nss_p"], phi["maxchi_p"]]
            n_sig = sum(1 for p in calls if p is not None and p < ALPHA)
            phipack_call = n_sig >= 2
            rows.append(dict(
                retain_frac=retain,
                retained_nongap_bp=retained_bp,
                threeseq_flag=t3["flagged"],
                threeseq_min_corr_p=t3["min_corr_p"],
                threeseq_n_triples=t3["n_triples"],
                phi_perm_p=phi["phi_perm_p"],
                nss_p=phi["nss_p"],
                maxchi_p=phi["maxchi_p"],
                phipack_n_sig=n_sig,
                phipack_call=phipack_call,
                two_detector_consensus=bool(t3["flagged"] and phipack_call),
            ))
            print(f"retain={retain:>4}  bp={retained_bp:>7}  "
                  f"3SEQ={'Y' if t3['flagged'] else 'n'}  "
                  f"Phi(perm/nss/maxchi)={phi['phi_perm_p']}/{phi['nss_p']}/{phi['maxchi_p']}  "
                  f"PhiPack={'Y' if phipack_call else 'n'}  "
                  f"consensus={'Y' if rows[-1]['two_detector_consensus'] else 'n'}")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "injection_partial_recombinant.tsv", sep="\t", index=False)

    # locate the discordance regime: 3SEQ flags but PhiPack does not
    disc = df[(df.threeseq_flag) & (~df.phipack_call)]
    summary = dict(
        control=CONTROL,
        alignment_len=L,
        panel_size=len(panel),
        panel_composition=dict(man.loc[man.accession.isin(panel), "clade"].value_counts()),
        full_length_phipack_call=bool(df.iloc[0]["phipack_call"]),
        full_length_threeseq_flag=bool(df.iloc[0]["threeseq_flag"]),
        discordance_first_retain=(float(disc.retain_frac.max()) if not disc.empty else None),
        discordance_first_bp=(int(disc.loc[disc.retain_frac.idxmax(), "retained_nongap_bp"])
                              if not disc.empty else None),
        discordance_levels=[dict(retain=float(r.retain_frac),
                                 retained_bp=int(r.retained_nongap_bp),
                                 phi_perm_p=r.phi_perm_p, nss_p=r.nss_p, maxchi_p=r.maxchi_p)
                            for r in disc.itertuples()],
    )
    (OUT / "injection_partial_recombinant.json").write_text(json.dumps(summary, indent=2, default=str))
    print("\n=== summary ===")
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
