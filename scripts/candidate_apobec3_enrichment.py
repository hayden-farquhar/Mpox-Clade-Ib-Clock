"""Characterise the 3SEQ over-call by clade and APOBEC3 load.

Reviewer point #2: a 3SEQ recombination flag on a COMPLETE genome cannot be a
truncation artefact, so what drives it?

The data answer the mechanism question directly: 3SEQ flags 100% of clade Ib
(198/198) and 0% of clade Ia (0/694), and completeness does NOT discriminate
(0/336 Ia partials flagged vs 104/104 Ib partials). The over-call tracks clade,
not genome length. This script quantifies (a) the clade-stratified flag rate by
completeness, (b) the APOBEC3 load of the flagged Ib panel vs the unflagged Ia
background, and (c) where the confirmed inter-clade recombinant sits.

Read-only on the freeze; writes one summary table + JSON. No detector re-run.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

FREEZE = Path(__file__).resolve().parent.parent / "data/processed/freeze_20260522"
OUT = Path(__file__).resolve().parent.parent / "outputs/tables"
OUT.mkdir(parents=True, exist_ok=True)

cand = pd.read_csv(FREEZE / "recomb_3seq/candidate_annotation_table.tsv", sep="\t")
ap = pd.read_csv(FREEZE / "apobec3_counts_v2.tsv", sep="\t")
meta = pd.read_csv(FREEZE / "metadata.tsv", sep="\t", low_memory=False)

cand = cand.merge(ap, on="accession", how="left")
control = cand[cand.is_positive_control == True].iloc[0]
nonctrl_acc = set(cand.loc[cand.is_positive_control != True, "accession"])

meta = meta.merge(ap, on="accession", how="left")
meta = meta.assign(flagged=meta.accession.isin(nonctrl_acc))

# (a) clade x completeness flag rate (metadata covers the 892 clade-I genomes).
flagtab = (meta.groupby(["clade", "Completeness"])["flagged"]
           .agg(flagged="sum", n="count").reset_index())
flagtab = flagtab.assign(rate=(flagtab.flagged / flagtab.n).round(4))

# (b) APOBEC3 load: flagged Ib vs unflagged Ia (and complete-only sub-comparison).
ib = meta[meta.clade == "Ib"].copy()
ia = meta[meta.clade == "Ia"].copy()
ib_complete = ib[ib.Completeness == "COMPLETE"].copy()
ia_complete = ia[ia.Completeness == "COMPLETE"].copy()


def desc(s):
    s = pd.to_numeric(s, errors="coerce").dropna()
    return dict(n=int(s.size), median=round(float(np.median(s)), 4),
                q1=round(float(np.percentile(s, 25)), 4),
                q3=round(float(np.percentile(s, 75)), 4))


def mwu(a, b):
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    u, p = mannwhitneyu(a, b, alternative="greater")
    return dict(U=float(u), p_greater=float(p))


summary = pd.DataFrame([
    dict(group="Ib_flagged_all", **desc(ib.apobec3_fraction)),
    dict(group="Ib_flagged_complete", **desc(ib_complete.apobec3_fraction)),
    dict(group="Ia_unflagged_all", **desc(ia.apobec3_fraction)),
    dict(group="Ia_unflagged_complete", **desc(ia_complete.apobec3_fraction)),
])

ctrl_frac = pd.to_numeric(pd.Series([control.apobec3_fraction]), errors="coerce").iloc[0]
ctrl_cnt = control.apobec3_snvs
ib_frac = pd.to_numeric(ib.apobec3_fraction, errors="coerce").dropna()
ctrl_in_counts = pd.notna(ctrl_frac)

result = dict(
    flag_rate_overall=dict(
        Ia=f"{int(meta[(meta.clade=='Ia')].flagged.sum())}/{int((meta.clade=='Ia').sum())}",
        Ib=f"{int(meta[(meta.clade=='Ib')].flagged.sum())}/{int((meta.clade=='Ib').sum())}"),
    apobec3_fraction_Ib_vs_Ia_all=mwu(ib.apobec3_fraction, ia.apobec3_fraction),
    apobec3_fraction_Ib_vs_Ia_complete=mwu(ib_complete.apobec3_fraction, ia_complete.apobec3_fraction),
    apobec3_count_Ib_vs_Ia_all=mwu(ib.apobec3_snvs, ia.apobec3_snvs),
    control_accession=str(control.accession),
    control_apobec3_in_cladeI_counts_file=bool(ctrl_in_counts),
    control_apobec3_fraction=(round(float(ctrl_frac), 4) if ctrl_in_counts else None),
    control_apobec3_snvs=(int(ctrl_cnt) if ctrl_in_counts else None),
    control_percentile_within_Ib=(round(float((ib_frac < ctrl_frac).mean() * 100), 1)
                                  if ctrl_in_counts else None),
)

summary.to_csv(OUT / "candidate_apobec3_enrichment.csv", index=False)
flagtab.to_csv(OUT / "flag_rate_by_clade_completeness.csv", index=False)
(OUT / "candidate_apobec3_enrichment.json").write_text(json.dumps(result, indent=2))

pd.set_option("display.width", 200, "display.max_columns", 30)
print("=== flag rate by clade x completeness ===")
print(flagtab.to_string(index=False))
print("\n=== APOBEC3 fraction by group ===")
print(summary.to_string(index=False))
print("\n=== tests + control position ===")
print(json.dumps(result, indent=2))
