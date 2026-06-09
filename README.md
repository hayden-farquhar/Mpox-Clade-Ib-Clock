# Mpox-Clade-Ib-Clock

Reproducibility code and result artefacts for:

**A local-minimum trap in saturation-aware APOBEC3-clock dating: pre-registered analysis of mpox clade Ib from the public corpus**

Hayden Farquhar MBBS MPHTM &nbsp; · &nbsp; Independent Researcher, Finley, New South Wales, Australia &nbsp; · &nbsp; ORCID: [0009-0002-6226-440X](https://orcid.org/0009-0002-6226-440X)

- **Preprint:** Zenodo deposit, DOI to be minted at first tagged release
- **Pre-registration:** [OSF Pre-Registration CASR2](https://doi.org/10.17605/OSF.IO/CASR2)
- **Project deposit (pre-registration + amendments + result artefacts):** [osf.io/gt3vx](https://osf.io/gt3vx)

This repository will be archived to Zenodo automatically on first tagged release; the resulting Zenodo DOI cites the exact code state at the time of preprint upload.

## Overview

This repository contains the full reproducibility chain for a pre-registered analysis of the public mpox virus (MPXV) genome corpus at a 2026-05-22 freeze. The analysis combines (a) a saturation-aware Poisson APOBEC3-clock dating of clade Ib, (b) a pre-registered three-detector recombination scan for inter-clade Ib/IIb recombinants under a ≥2/3 consensus rule, and (c) three methodological observations about counting-protocol calibration, alignment-method sensitivity of branch-quantity counts, and a per-tip Spearman ρ diagnostic for separating alignment-noise from biological signal.

The headline methodological finding is a **local-minimum trap** in the saturation-Poisson likelihood landscape: the registered L-BFGS-B fitting procedure converged to a local minimum approximately 260 log-likelihood units worse than the global ML minimum identified by a direct (λ, t₀) grid search. The repository includes the L-BFGS-B random-start convergence trajectory data that makes this finding directly reproducible.

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
| `baseline_plot.py` | Figure 2: tip-vs-reference APOBEC3 SNV accumulation, clade Ia vs Ib | 3 |
| `branch_baseline_plot.py` | Figure 1: per-tip branch-quantity APOBEC3 evidence for H1' | 3 |
| `build_panels.py` | Build the §5.4 outgroup panel and three parental reference panels (stratified random sampling, seed 42) | 4 |
| `run_phipack_consensus.py` | Per-candidate PhiPack runs for the §5.4 consensus filter (Amendment 03; PhiPack substituted for RDP5) | 4 |
| `dropcontrol_recomb.py` | Parent-panel sensitivity test: re-run 3SEQ on the 198 clade-Ib children with the positive control OZ375330.1 removed from the candidate-parent panel | 4 |
| `build_dating_clusters.py` | Define dating clusters per §3.5 + §5.0 step 7 (Nextstrain-lineage mode) | 3B |
| `fit_saturation_dating.py` | Fit the §5.3 saturation-aware Poisson dating model (primary L-BFGS-B with 25 random starts) | 3B |
| `treetime_strict_clock_xcheck.py` | §5.3 TreeTime strict-clock cross-check on the dating cluster | 3B |
| `lbfgsb_trajectory.py` | §6.4 Amendment 04: capture the L-BFGS-B convergence trajectory of the 25 random starts | 3B |
| `plot_likelihood_landscape.py` | Figure 4: regenerate the (λ, t₀) likelihood-landscape figure with the registered primary + global ML minimum + TreeTime cross-check annotated | 3B |
| `sensitivity_runner.py` | §5.5 sensitivity analyses (items 2 random subsampling, 3 temporal hold-out, 6 outgroup-composition perturbation) | 5 |
| `sensitivity_detector_sweep.py` | §5.5 item 4: detector-parameter sweep (3SEQ alpha-scale, PhiPack internal-agree, GARD ΔAIC) | 5 |

## Outputs

### Figures (publication, 300 dpi)

| File | Paper reference | What it shows |
|---|---|---|
| `outputs/figures/Figure_1_branch_snv_by_date.png` | Figure 1 | Per-tip branch-quantity APOBEC3 evidence for H1' |
| `outputs/figures/Figure_2_snv_by_date.png` | Figure 2 | Tip-vs-reference APOBEC3 SNV accumulation, clade Ia vs Ib |
| `outputs/figures/Figure_3_breakpoint_concordance.png` | Figure 3 | Cross-detector breakpoint concordance for the positive control OZ375330.1 |
| `outputs/figures/Figure_4_likelihood_landscape.png` | Figure 4 | Saturation-Poisson (λ, t₀) likelihood landscape with primary fit + global ML + TreeTime annotated |
| `outputs/figures/Figure_5_sensitivity_bimodality.png` | Figure 5 | Saturation-Poisson MRCA bimodality across the 11 §5.5 sensitivity replicates |

### Supplementary tables

| File | Paper reference | Content |
|---|---|---|
| `outputs/tables/Table_S1_candidate_annotation.tsv` | Supplementary Table S1 | Per-candidate annotation of the 202 sequences flagged by 3SEQ as recombinant children |
| `outputs/tables/Table_S2_sensitivity_summary.tsv` | Supplementary Table S2 | Per-replicate sensitivity summary across 11 §5.5 replicates |
| `outputs/tables/Table_S3_detector_sweep.csv` | Supplementary Table S3 | Detector-parameter sweep summary (§5.5 item 4) |
| `outputs/tables/Table_S4_lbfgsb_trajectory.csv` | Supplementary Table S4 | L-BFGS-B random-start convergence trajectory (25 starts × 7 columns) for the saturation-Poisson dating fit |

### Recombination parent-panel sensitivity (drop-control test)

| File | Content |
|---|---|
| `outputs/apobec3_mask_experiment/run_dropcontrol/dropcontrol_summary.json` | Result of removing the positive control OZ375330.1 from the candidate-parent panel: clade-Ib flag count falls 198/198 → 10/198, with donor composition |
| `outputs/apobec3_mask_experiment/run_dropcontrol/dropcontrol.3s.rec.csv` | 3SEQ recombinant-triple records for the drop-control run |
| `outputs/apobec3_mask_experiment/run_dropcontrol/dropcontrol.3s.log` | 3SEQ run log (triples tested, runtime) |
| `outputs/apobec3_mask_experiment/run_dropcontrol/dropcontrol.3s.pvalHist` | 3SEQ corrected-*P* histogram for the drop-control run |

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
Farquhar, H. (2026). Mpox-Clade-Ib-Clock: Reproducibility code for "A local-minimum
trap in saturation-aware APOBEC3-clock dating: pre-registered analysis of mpox
clade Ib from the public corpus". Zenodo. https://doi.org/[Zenodo DOI on release]
```

A `CITATION.cff` file is included at the repository root; GitHub will surface a "Cite this repository" button once the repo is public.

## Licence

Code: **MIT License** (see `LICENSE`)
Data, documentation, and figures: **Creative Commons Attribution 4.0 International (CC-BY 4.0)** (see `LICENSE`)

## Reporting issues

Bug reports, reproduction issues, and methodological questions: please open a GitHub Issue at <https://github.com/hayden-farquhar/Mpox-Clade-Ib-Clock/issues>. For scientific correspondence, contact the author at <hayden.farquhar@icloud.com>.
