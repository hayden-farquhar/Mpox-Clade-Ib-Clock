# Mpox-Clade-Ib-Clock

Reproducibility code and result artefacts for:

**Isolating inter-clade recombination in mpox genomic surveillance: a two-detector protocol that separates one genuine event from intra-clade and partial-genome signal**

Hayden Farquhar MBBS MPHTM &nbsp; · &nbsp; Independent Researcher, Finley, New South Wales, Australia &nbsp; · &nbsp; ORCID: [0009-0002-6226-440X](https://orcid.org/0009-0002-6226-440X)

- **Preprint:** Zenodo deposit, DOI to be minted at first tagged release
- **Pre-registration:** [OSF Pre-Registration CASR2](https://doi.org/10.17605/OSF.IO/CASR2)
- **Project deposit (pre-registration + amendments + result artefacts):** [osf.io/gt3vx](https://osf.io/gt3vx)

This repository will be archived to Zenodo automatically on first tagged release; the resulting Zenodo DOI cites the exact code state at the time of preprint upload.

## Overview

This repository contains the full reproducibility chain for a pre-registered analysis of the public mpox virus (MPXV) genome corpus at a 2026-05-22 freeze. The headline contribution is a recombination-surveillance result: when a single genome-wide detector (3SEQ) is run over the clade-I corpus it flags every clade-Ib genome (198/198) and no clade-Ia genome, independent of genome completeness, whereas requiring confirmation by a second, independent detector (PhiPack) isolates the one genuine inter-clade Ib/IIb recombinant (the WHO/Pullan positive control OZ375330.1) from the clade-Ib flag pool. The analysis quantifies the specificity that the two-detector rule buys and the partial-genome sensitivity it costs, and supplies two sensitivity tests behind that claim: a within-clade APOBEC3-column / random-column masking experiment, and a parent-panel drop-control test (removing the positive control as a candidate parent collapses the flag count 198/198 → 10/198, showing the over-call magnitude is largely induced by the known recombinant acting as a universal donor).

A pre-registered saturation-aware APOBEC3 molecular-clock analysis of clade Ib is reported as a secondary methods note: it recovers the expected APOBEC3 human-transmission signature, dates the clade-Ib MRCA to mid-to-late 2023 (TreeTime strict-clock cross-check, consistent with the Kamituga emergence) with wide uncertainty, and records two transferable cautions for groups applying the O'Toole APOBEC3-clock framework to a newly emerging clade: a counting-protocol/null calibration mismatch, and a **local-minimum trap** in the saturation-Poisson likelihood landscape (the registered L-BFGS-B fit settles ~260 log-likelihood units worse than the global ML optimum found by a direct (λ, t₀) grid search). The L-BFGS-B random-start convergence trajectory data is included so the latter is directly reproducible.

## Data sources

| Source | URL | Access | Licence |
|---|---|---|---|
| NCBI Virus (MPXV genomes) | <https://www.ncbi.nlm.nih.gov/labs/virus/vssi/> | Free, public | Public domain |
| Nextstrain mpox build (clade-i, all-clades) | <https://nextstrain.org/mpox> / <https://github.com/nextstrain/mpox> | Free, public | MIT (code), CC-BY (analyses) |
| Nextclade mpox clade-i dataset | <https://github.com/nextstrain/nextclade_data/tree/master/data/nextstrain/mpox/clade-i> | Free, public | MIT |
| Pullan et al. 2025 positive control (OZ375330.1) | <https://www.ncbi.nlm.nih.gov/nuccore/OZ375330.1> | Free, public | Public domain |

All raw sequence data and metadata are publicly accessible at the freeze date (2026-05-22). The corpus FASTA file produced by the pipeline (`data/processed/corpus.fasta`, ~360 MB, 893 sequences) is not committed to this repository; instructions for regenerating it from NCBI are in `data/raw/README.md`. GISAID is excluded by design.

## Requirements

- Python 3.11
- conda/mamba (recommended) or pip
- External bioinformatics tools (must be installed separately — see `env/README.md` for guidance): Nextclade v3.21.2, MAFFT v7.526, IQ-TREE v3.1.2, TreeTime v0.12.1, HyPhy v2.5.99 (for GARD), 3SEQ v1.9.1, PhiPack v1.1, NCBI Datasets CLI v18.28.0

### Conda / mamba (recommended)

```bash
mamba env create -f env/environment.yml
mamba activate clade-ib-clock
```

### Pip (fallback for Python-only dependencies; external bioinformatics tools must be installed separately)

```bash
pip install -r env/requirements.txt
```

The full pinned-version list used in the analysis is documented in §2.11 of the accompanying manuscript.

## Reproduction

The analysis is orchestrated by a Snakemake workflow at `workflow/Snakefile`. The pipeline is broken into three phases that mirror the registered §5.0 execution order; each phase has its own `.smk` rule module.

```bash
# Phase 1: data freeze (NCBI fetch + Nextclade QC + §3.3 filters)
snakemake --cores all --config FREEZE_DATE=20260522 freeze_complete

# Phase 2: alignment + mask
snakemake --cores all align_and_mask

# Phase 3: APOBEC3 counting + H1/H1' tests
snakemake --cores all phase3_complete

# Phase 4: recombination scan (3SEQ + PhiPack + GARD)
snakemake --cores all phase4_complete

# Phase 3B: saturation-aware dating + sensitivity
snakemake --cores all dating_complete sensitivity_complete
```

Phases 4 and 5 also include cloud-compute steps (GARD on a 31-sequence subset; full-iterative MAFFT for the §5.5 alignment-method sensitivity); see `docs/` for the analysis log of those steps.

For a single-command exploratory reproduction of the headline statistical results (H1, H1', H2, H3, H4 + the L-BFGS-B trajectory + Spearman ρ diagnostic) from the deposited intermediate artefacts, see the section "Reproducing headline numbers" below.

## Script descriptions

| Script | Description | Phase |
|---|---|---|
| `build_freeze_manifest.py` | Build the freeze manifest with provenance, SHA-256 hashes, and tool versions | 1 |
| `apply_filters.py` | Apply pre-registration §3.3 inclusion / exclusion filters to the raw NCBI corpus | 1 |
| `build_default_mask.py` | Fallback mask builder (ITR-only); the production mask is the upstream nextstrain/mpox BED, see Amendment 01 | 2 |
| `apply_mask.py` | Apply a BED file of masked intervals to every sequence in a multiple-sequence alignment | 2 |
| `count_apobec3.py` | Count APOBEC3-context SNVs per genome relative to the clade-I reference (tip-vs-reference protocol; §5.2) | 3 |
| `build_branch_apobec3_counts.py` | Build per-tip branch-quantity APOBEC3 SNV counts for the H1' test (Amendment 02) | 3 |
| `test_h1.py` | H1 / H1' one-sided exact-binomial test against the IIb-calibrated null *p*₀ = 0.70 | 3 |
| `baseline_plot.py` | Supplementary Figure S1: tip-vs-reference APOBEC3 SNV accumulation, clade Ia vs Ib (`snv_by_date_v2.png`) | 3 |
| `branch_baseline_plot.py` | Diagnostic: per-tip branch-quantity APOBEC3 evidence for H1' (`branch_snv_by_date.png`) | 3 |
| `build_panels.py` | Build the §5.4 outgroup panel and three parental reference panels (stratified random sampling, seed 42) | 4 |
| `run_phipack_consensus.py` | Per-candidate PhiPack runs for the §5.4 consensus filter (Amendment 03; PhiPack substituted for RDP5) | 4 |
| `dropcontrol_recomb.py` | Parent-panel sensitivity test: re-run 3SEQ on the 198 clade-Ib children with the positive control OZ375330.1 removed from the candidate-parent panel (198/198 → 10/198) | 4 |
| `apobec3_mask_recomb.py` | Within-clade APOBEC3-context vs random-column vs baseline masking experiment: isolates whether APOBEC3 homoplasy (not clade per se) drives the 3SEQ flag (Table S6/S7) | 4 |
| `random_mask_distribution.py` | Multi-seed (n=10) distribution for the matched random-column masking control, so the control is reported as a median and range rather than one draw (Table S7) | 4 |
| `candidate_apobec3_enrichment.py` | Characterise the 3SEQ over-call by clade and genome completeness, and the APOBEC3 load of the flagged Ib panel vs the unflagged Ia background | 4 |
| `injection_partial_recombinant.py` | Partial-genome sensitivity: progressively truncate the confirmed inter-clade recombinant and re-run both detectors to test what the two-detector rule discards on incomplete assemblies | 4 |
| `plot_breakpoint_concordance.py` | Figure 1: cross-detector breakpoint concordance for the positive control (3SEQ primary/alternative triplets + GARD two-breakpoint model) | 4 |
| `plot_overcall_partial_genomes.py` | Figure 2: genome length vs PhiPack Phi-permutation *P* for every 3SEQ candidate + candidate genome-length distribution | 4 |
| `build_dating_clusters.py` | Define dating clusters per §3.5 + §5.0 step 7 (Nextstrain-lineage mode) | 3B |
| `fit_saturation_dating.py` | Fit the §5.3 saturation-aware Poisson dating model (primary L-BFGS-B with 25 random starts) | 3B |
| `refit_dating_global.py` | Global-ML re-fit of the saturation-Poisson MRCA model via a dense (t₀, log₁₀λ) grid + L-BFGS-B polish, with within-basin bootstrap — locates the global optimum the registered fit's random starts miss | 3B |
| `treetime_strict_clock_xcheck.py` | §5.3 TreeTime strict-clock cross-check on the dating cluster | 3B |
| `lbfgsb_trajectory.py` | §6.4 Amendment 04: capture the L-BFGS-B convergence trajectory of the 25 random starts | 3B |
| `plot_likelihood_landscape.py` | Diagnostic: (λ, t₀) likelihood-landscape figure with the registered primary + global ML minimum + TreeTime cross-check annotated (`likelihood_landscape.png`) | 3B |
| `sensitivity_runner.py` | §5.5 sensitivity analyses (items 2 random subsampling, 3 temporal hold-out, 6 outgroup-composition perturbation) | 5 |
| `sensitivity_detector_sweep.py` | §5.5 item 4: detector-parameter sweep (3SEQ alpha-scale, PhiPack internal-agree, GARD ΔAIC) | 5 |

## Outputs

### Figures (publication, 300 dpi)

| File | Paper reference | What it shows | Generated by |
|---|---|---|---|
| `outputs/figures/breakpoint_concordance.png` | Figure 1 | Cross-detector breakpoint concordance for the positive control OZ375330.1 | `plot_breakpoint_concordance.py` |
| `outputs/figures/overcall_partial_genomes.png` | Figure 2 | Genome length vs PhiPack Phi-permutation *P* for every 3SEQ candidate + candidate length distribution | `plot_overcall_partial_genomes.py` |
| `outputs/figures/snv_by_date_v2.png` | Supplementary Figure S1 | Tip-vs-reference APOBEC3 SNV accumulation, clade Ia vs Ib | `baseline_plot.py` |
| `outputs/figures/branch_snv_by_date.png` | Diagnostic (methods note) | Per-tip branch-quantity APOBEC3 evidence for H1' | `branch_baseline_plot.py` |
| `outputs/figures/likelihood_landscape.png` | Diagnostic (methods note) | Saturation-Poisson (λ, t₀) likelihood landscape with primary fit + global ML + TreeTime annotated | `plot_likelihood_landscape.py` |
| `outputs/figures/sensitivity_mrca_bimodality.png` | Diagnostic (methods note) | Saturation-Poisson MRCA bimodality across the §5.5 sensitivity replicates | `sensitivity_runner.py` |

### Supplementary tables

| Paper reference | Backing file(s) | Content |
|---|---|---|
| Supplementary Table S1 | `outputs/tables/Table_S1_candidate_annotation.tsv` | Per-candidate annotation of the sequences flagged by 3SEQ as recombinant children |
| Supplementary Table S2 | `outputs/tables/Table_S2_sensitivity_summary.tsv` | MRCA fit by subsample (the saturation-versus-uncorrected difference as a basin artefact) |
| Supplementary Table S3 | `outputs/tables/Table_S3_detector_sweep.csv` | Detector-parameter sweep summary (§5.5 item 4) |
| Supplementary Table S4 | `outputs/tables/Table_S4_lbfgsb_trajectory.csv` | L-BFGS-B random-start convergence trajectory for the saturation-Poisson dating fit |
| Supplementary Table S5 | `outputs/tables/injection_partial_recombinant.tsv`, `.json` | Partial-genome injection test: detector retention versus retained genome length |
| Supplementary Table S6 | `outputs/tables/candidate_apobec3_enrichment.csv`, `.json`, `outputs/tables/flag_rate_by_clade_completeness.csv` | Single-detector flag rate by clade and completeness, with APOBEC3 load |
| Supplementary Table S7 | `outputs/apobec3_mask_experiment/results_summary.json`, `random_distribution_results.json` | APOBEC3-column masking experiment (the triplet signal is not specifically APOBEC3-driven) |
| Main-text Table 1 | `outputs/tables/baseline_v2.csv` | Baseline 3SEQ flag table underlying the headline 198/198 clade-Ib flag rate |

### APOBEC3-column masking experiment and parent-panel drop-control (Supplementary Table S7 + sensitivity)

The masking experiment runs 3SEQ under three conditions that differ only in which alignment columns are masked (baseline, APOBEC3-context homoplasy columns, and a matched random-column control over ten seeds), plus a parent-panel drop-control. The large per-condition input alignments (`aln_*.fasta`, ~174 MB each) are **not** committed; they are regenerated by `apobec3_mask_recomb.py` from the freeze alignment. The small 3SEQ run artefacts and summaries are deposited:

| File | Content |
|---|---|
| `outputs/apobec3_mask_experiment/results_summary.json` | Baseline vs APOBEC3-masked vs single-seed random-mask flag counts |
| `outputs/apobec3_mask_experiment/random_distribution_results.json` | Ten-seed random-column control distribution (n=10; median 119/198, range 64–197, mean 126.4, SD 54.0) |
| `outputs/apobec3_mask_experiment/run_baseline/`, `run_apobec3/`, `run_random{,_seed1..9}/` | 3SEQ logs, rec.csv, longRec, pvalHist for every masking condition/seed |
| `outputs/apobec3_mask_experiment/per_ib_child_minp.tsv` | Per clade-Ib child minimum corrected *P* across conditions |
| `outputs/apobec3_mask_experiment/run_dropcontrol/dropcontrol_summary.json` | Parent-panel drop-control: removing OZ375330.1 as a candidate parent collapses the flag count 198/198 → 10/198, with donor composition |
| `outputs/apobec3_mask_experiment/run_dropcontrol/dropcontrol.3s.{rec.csv,log,pvalHist}` | 3SEQ artefacts for the drop-control run |

### Deposited result JSONs

| File | Content |
|---|---|
| `data/processed/h1prime_result.json` | H1' test result: Ib branch-quantity APOBEC3 fraction 0.8358 (Wilson 95% 0.8131–0.8563), *P* = 4.5 × 10⁻²⁶ |
| `data/processed/saturation_dating_results.json` | Saturation-aware Poisson dating fit: registered primary *t*₀ = 2019-01-09, λ = 9.99 × 10⁻⁸ /site/day, profile-LL 3,030.77 units |
| `data/processed/lbfgsb_trajectory.csv` | L-BFGS-B random-start trajectory: 2 of 25 starts in the 2019 local-minimum basin, 0 of 25 in the 2023 global-ML basin |
| `data/processed/freeze_manifest.json` | Freeze manifest: SHA-256 hashes of input artefacts at the 2026-05-22 freeze |

## Reproducing headline numbers

From the deposited intermediate artefacts in `data/processed/`, the headline statistical results can be reproduced with a few lines of Python:

```python
import json, pandas as pd
from scipy.stats import binomtest

# H1' test (Amendment 02 branch-quantity protocol)
counts = pd.read_csv("data/processed/apobec3_counts_branch.tsv", sep="\t")  # generated by build_branch_apobec3_counts.py
k = int(counts["apobec3_snvs"].sum())
n = int(counts["total_snvs"].sum())
print(f"Ib branch-quantity APOBEC3 fraction = {k/n:.4f} ({k:,}/{n:,})")
print(f"One-sided exact binomial P (k >= obs | p0=0.70) = {binomtest(k, n, 0.70, 'greater').pvalue:.2e}")

# L-BFGS-B trajectory basin counts (Amendment 04)
traj = pd.read_csv("data/processed/lbfgsb_trajectory.csv")
print(traj.groupby("basin").size())
```

## Pre-registration and §6.4 deviation amendments

The full analysis plan was pre-registered at the Open Science Framework before any sequence data was retrieved: OSF Pre-Registration CASR2 (<https://doi.org/10.17605/OSF.IO/CASR2>). Four deviation amendments were filed under §6.4 of the pre-registration during execution; the full text of each is in `docs/`:

- `docs/osf_amendment_01_mask_deviation.md` (2026-05-22) — Mask source substitution
- `docs/osf_amendment_02_h1_calibration.md` (2026-05-22) — Counting-protocol calibration (branch-quantity H1' exploratory analysis)
- `docs/osf_amendment_03_phipack_substitution.md` (2026-05-23) — Recombination-detector substitution (PhiPack for RDP5)
- `docs/osf_amendment_04_spearman_alignment_diagnostic.md` (2026-05-25) — Post-hoc per-tip Spearman ρ alignment-noise diagnostic + direct (λ, t₀) grid-search

Each amendment is also deposited at OSF under the parent project [gt3vx](https://osf.io/gt3vx) with server-issued timestamps preserved.

## Citation

If you use this code or any of its outputs, please cite both the preprint and the Zenodo deposit of this repository. After Zenodo mints the DOI on first tagged release, the canonical citation is:

```
Farquhar, H. (2026). Mpox-Clade-Ib-Clock: Reproducibility code for "Isolating
inter-clade recombination in mpox genomic surveillance: a two-detector protocol
that separates one genuine event from intra-clade and partial-genome signal".
Zenodo. https://doi.org/[Zenodo DOI on release]
```

A `CITATION.cff` file is included at the repository root; GitHub will surface a "Cite this repository" button once the repo is public.

## Licence

Code: **MIT License** (see `LICENSE`)
Data, documentation, and figures: **Creative Commons Attribution 4.0 International (CC-BY 4.0)** (see `LICENSE`)

## Reporting issues

Bug reports, reproduction issues, and methodological questions: please open a GitHub Issue at <https://github.com/hayden-farquhar/Mpox-Clade-Ib-Clock/issues>. For scientific correspondence, contact the author at <hayden.farquhar@icloud.com>.
