# OSF Pre-Registration: Quantitative APOBEC3-Signature Dating and Recombination Scan of Mpox Clade Ib

**Template:** OSF Standard Pre-Data Collection Registration
**Registration date:** 2026-05-22
**Authors:** Hayden Farquhar MBBS MPHTM (Independent Researcher, Finley, NSW, Australia; ORCID 0009-0002-6226-440X)
**Status at registration:** No clade Ib sequence data has been downloaded, aligned, or analysed by the author. The environment and reference Nextclade dataset have been verified as reachable; the data freeze has not yet been taken. The OSF registration timestamp on this document serves as the priority anchor for the analysis plan recorded herein.

---

## Study at a glance

| Hypothesis | Quantity tested | Test | Pre-specified threshold for support |
|---|---|---|---|
| H1 | Clade-Ib APOBEC3 SNV fraction | One-sided exact binomial vs p₀ = 0.70 | Wilson lower 95% bound > 0.70 |
| H2 | Per-cluster Δ between corrected and uncorrected emergence date | Δ in days, max across clusters | max Δ > 60 days |
| H3 | Count of consensus candidate recombinants (excl. WHO Dec 2025 case) | Count at ≥2-of-3 detector agreement (default thresholds) | Count ≥ 1 |
| H4 | Recovery of the WHO Dec 2025 positive-control recombinant | ≥2 of the three detectors call the case (full 3/3 vs 2/3 also reported) | ≥ 2 / 3 |
| H5 | Jaccard similarity of consensus-candidate sets under alignment and subsampling perturbations | J = \|A ∩ B\| / \|A ∪ B\| | J ≥ 0.80 (each perturbation independently) |

Primary data source: NCBI Virus (anonymous Entrez API). Reference genome: NC_063383 (clade I). Detectors: 3SEQ, RDP5 ensemble, GARD. GISAID is excluded by design.

---

## Glossary of key terms

- **APOBEC3 context:** The dinucleotide sequence context (5′-T·C-3′ on one strand, equivalently 5′-G·A-3′ on the opposite strand at the same physical site) targeted by the APOBEC3 cytidine deaminase enzymes. Edits at these sites (C→T or G→A respectively) accumulate during sustained human-to-human transmission of MPXV.
- **MRCA:** Most recent common ancestor of a cluster of sequences; equivalently, the inferred date at which the cluster originated.
- **Saturation:** The condition in which a meaningful fraction of available APOBEC3 target sites have been edited more than once, so that the observed edit count under-estimates the true mutational accumulation.
- **Consensus candidate (this study):** A genome called as a recombinant by ≥2 of the three pre-specified recombination detectors at their default significance thresholds.
- **Positive control (this study):** The WHO Disease Outbreak News mpox recombinant case reported 14 February 2026 and reclassified in the Nextclade mpox dataset on 16 December 2025.

---

## 1. Study Information

### 1.1 Title

Quantitative APOBEC3-Signature Dating of Mpox Clade Ib Transmission Clusters and a Systematic Recombination Scan of the Public Clade I + Ib/IIb GenBank Corpus

### 1.2 Description

Monkeypox virus (MPXV) clade IIb sustained human-to-human transmission since 2022 has been shown to accumulate a strong APOBEC3 cytidine-deamination signature (TC→TT on one strand, equivalently GA→AA on the opposite strand at the same physical site), and this signature has been used as a quantitative molecular clock (O'Toole et al., *Science* 2023). The 2024–2026 expansion of clade Ib transmission in central and eastern Africa, the first imported clade Ia case in China, and the WHO Disease Outbreak News report (14 February 2026) of an Ib/IIb recombinant lineage have established a need for the same quantitative dating framework, applied to clade Ib, with explicit accounting for the limited pool of TC/GA contexts (saturation) at current sample sizes.

The most recent published work on MPXV clade Ib genomics has approached recombination from complementary but non-overlapping angles. Feehley et al. (*Genetics*, 2026) develop a haploid linkage-disequilibrium model on the 2023–2025 Ib corpus and report that the majority of recombinant SNP pairs show strong evidence of historical recombination, but do not perform breakpoint-resolved scans or molecular-clock dating. Alfonsi et al. (*J Mol Biol*, 2026) extend their RecombinHunt framework to MPXV and other viruses with a publicly-hosted automated tool, but the framework is designed for cross-virus monitoring rather than clade-Ib-specific findings and does not couple to a clock. The CDC *Emerging Infectious Diseases* review of clade Ib emergence (Aug 2025) notes APOBEC3 mutational signatures as evidence of sustained transmission but explicitly does not perform dating or recombination scans, and calls for such analyses. No peer-reviewed paper or preprint identified at the time of registration combines (a) a saturation-corrected APOBEC3 molecular clock on clade Ib clusters with (b) a breakpoint-resolved consensus recombination scan of the public clade I + Ib/IIb corpus, with the WHO-reported Ib/IIb case as a positive control.

This study has two confirmatory aims: (i) date the major clade Ib transmission clusters using an APOBEC3-context-aware Poisson model that explicitly accounts for site saturation, and (ii) perform a pre-specified three-method recombination scan (3SEQ, RDP5 ensemble, GARD) of the public clade I + Ib/IIb corpus to detect additional inter-clade recombinant lineages beyond the WHO-reported case. All sequence data and analytical inputs are drawn from open NCBI, Nextstrain, and Nextclade resources; the Virological.org public forum is consulted for documentation purposes only and does not influence detector calls. GISAID is excluded by design.

### 1.3 Research questions

1. **RQ1 — Dating.** What are the estimated emergence dates (median + 95% credible interval) of the major clade Ib transmission clusters when APOBEC3 edit counts are modelled with an explicit per-site saturation correction?
2. **RQ2 — Saturation.** At current clade Ib sample sizes, what fraction of inferred APOBEC3 sites are likely to be multiply hit, and how much do uncorrected emergence-date estimates differ from the corrected estimates?
3. **RQ3 — Recombination.** Are there detectable additional Ib/IIb (or Ia-involving) recombinant genomes in the public corpus beyond the WHO Dec 2025 / Nextclade 16-Dec-2025 index case? What are their breakpoint distributions and parental clade assignments?
4. **RQ4 — Sensitivity.** How robust are recombinant calls and emergence-date estimates to (a) alignment method, (b) recombination-detector choice, (c) subsampling of the input corpus, and (d) the choice of saturation-correction model?

### 1.4 Hypotheses

All hypotheses are pre-specified before any sequence is downloaded. Each is paired with a specific statistical test and a quantitative decision threshold in §5.

**H1 (APOBEC3 dominance — positive-control sanity check).** The proportion of single-nucleotide variants on the clade Ib alignment that fall in the APOBEC3 dinucleotide context (TC→TT, with the equivalent GA→AA on the opposite strand) will exceed 0.70, consistent with the value reported for clade IIb by O'Toole et al. (2023). *Test:* one-sided one-sample exact binomial test against the null p₀ = 0.70 at α = 0.05; H1 is supported if the lower bound of the 95% Wilson interval exceeds 0.70. *Effect size:* the observed proportion with its 95% Wilson interval. *Power:* with the expected ≥500 SNVs across a clade-Ib corpus of ≥30 genomes, the test has >99% power to reject the null if the true value matches the IIb-reported ≈0.93, and >80% power if the true value is ≥0.78. *Failure mode:* an observed value <0.50 halts downstream analysis and triggers a documented investigation (misalignment, contamination, or fundamentally different mutational process).

**H2 (Saturation bias is non-trivial).** For at least one eligible clade-Ib cluster (§3.5), the absolute difference between the saturation-corrected and the uncorrected Poisson emergence-date estimate will exceed 60 days. *Test:* per-cluster Δ = |t₀(corrected) − t₀(uncorrected)| in days, with H2 supported if max-over-clusters Δ > 60 days. *Effect size:* the per-cluster Δ in days, with its bootstrap 95% interval. *Threshold rationale:* 60 days is the smallest difference that would meaningfully change a public-health narrative about whether a transmission cluster predates a documented introduction event; it is also approximately one mpox serial-interval doubling window. *Power note:* with ≥10 genomes per eligible cluster and ~50 APOBEC3 SNVs per genome, the parametric-bootstrap CI on per-cluster Δ is expected to be in the ~60–120 day range; the test is therefore powered to detect Δ values exceeding the threshold for clusters with substantial saturation but may fail to support H2 for clusters with very mild saturation even if a meaningful systematic bias exists. The signed per-cluster Δ is reported regardless of the test outcome. *Falsifiable:* if max Δ ≤ 60 days for all clusters, we conclude that the simpler uncorrected estimator is adequate at current sample sizes.

**H3 (Additional recombinants exist in the public corpus).** At least one genome other than the WHO Dec 2025 / Nextclade 16-Dec-2025 index case will be called as a candidate inter-clade recombinant by ≥2 of the three pre-specified detectors (3SEQ, RDP5 ensemble, GARD) at detector-default significance thresholds. *Test:* count of consensus candidates excluding the positive control; H3 supported if count ≥1. *Effect size:* the integer count *k*, together with a Wilson 95% interval for the population proportion *k*/*N* (where *N* is the corpus size after exclusions); the proportion is descriptive only — the H3 test itself is on the count. *Falsifiable:* a null result (zero additional consensus candidates) is equally publishable and would constitute evidence that, at the resolution of the current corpus, the WHO-reported case is an isolated event. *Upper-bound caution:* if the count exceeds 10, we will treat this as a possible detector-failure signal (saturation of the multi-method consensus by alignment artefact) and report it only after a documented manual inspection of each candidate alignment block.

**H4 (Positive control — known recombinant is recovered).** The WHO Dec 2025 / Nextclade 16-Dec-2025 index recombinant genome will be called as a recombinant by at least two of the three detectors (≥2/3 consensus), matching the threshold used for the H3 consensus-candidate definition. *Test:* binary, ≥2/3 required. *Effect size:* the per-detector significance value reported for the positive control. *Rationale for ≥2/3 rather than 3/3:* whole-genome GARD is known to have lower sensitivity than 3SEQ and RDP5 on long viral genomes with sparse breakpoints; a single-detector negative on the positive control therefore does not by itself indicate methodology failure, but two negatives (any pair) would. *Failure mode:* if fewer than 2/3 detectors call the positive control, H3 reporting is halted until the affected detector configurations are re-verified against their published documentation; if the failure persists with documented defaults, RQ3 is downgraded to exploratory and so reported. The full 3/3 vs 2/3 result is reported for transparency regardless of outcome.

**H5 (Robustness to alignment and subsampling).** The set of consensus candidate recombinants identified in H3 will be ≥80% concordant (Jaccard similarity) between (i) the primary MAFFT-based alignment and an alternative Nextalign-based alignment, and (ii) the full corpus and the union of five 75% random-subsample re-runs. *Test:* per-perturbation Jaccard similarity J = |A ∩ B| / |A ∪ B|; H5 supported separately for (i) and (ii) if J ≥ 0.80. *Effect size:* J with bootstrap 95% interval over the candidate set. *Threshold rationale:* 0.80 is adopted as a stringent stability threshold; the result is reported with its bootstrap CI regardless of whether the threshold is met. *Falsifiable consequence:* candidates lost under either perturbation are downgraded to exploratory and reported in the supplement only.

---

## 2. Design Plan

### 2.1 Study type

Retrospective phylogenetic / molecular epidemiology analysis of publicly available viral genome sequences. No human participants. No private or controlled-access data.

### 2.2 Study design

Four-phase, fully confirmatory design (every test is pre-specified in §5):

1. **Baseline replication (tests H1).** Re-implement the O'Toole et al. APOBEC3 edit-counting framework on the clade Ib alignment and confirm signature dominance.
2. **Recombination scan and positive-control verification (tests H3 and H4; answers RQ3).** Run all three pre-specified detectors with frozen default parameters; apply the ≥2/3 consensus rule; verify positive-control recovery.
3. **Quantitative dating with saturation correction (tests H2; answers RQ1 and RQ2).** Fit the pre-specified context-aware Poisson model and a competing uncorrected Poisson model to each cluster (after excluding the recombinants identified in phase 2); report emergence-date estimates with parametric-bootstrap uncertainty.
4. **Sensitivity (tests H5; answers RQ4).** Repeat phases 2 and 3 under alternative alignment, subsampling, and saturation-model choices.

Phase ordering matches the pipeline order in §5.0. Phase 2 is executed before phase 3 so that recombinants can be removed from cluster membership before dating models are fit.

### 2.3 Exploratory analyses

This pre-registration commits to a fully confirmatory analysis; no analyses are designated as exploratory at the time of registration. Any analysis added after registration that responds to features of the data observed during the analysis itself will be clearly labelled as exploratory in the manuscript and excluded from the confirmatory hypothesis tests of §1.4. The single permitted exception, anticipated *a priori*, is a follow-up qualitative description of any consensus candidate recombinants identified in §5.4 (clinical / epidemiological context per accession, breakpoint visualisation, parental-clade tree placement) — these descriptions are not hypothesis tests and are reported as descriptive supplements to the confirmatory H3 count.

### 2.4 Blinding

Not applicable. All analyses are deterministic functions of the input alignment and metadata; outcome scoring is fully automated.

### 2.5 Reporting guidelines

No formal community-endorsed reporting checklist exists for APOBEC3-clock dating or breakpoint-resolved recombination scanning of viral genomes. Methods reporting will follow the data and code transparency principles of the Nextstrain pathogen-build documentation; recombination analyses will follow the RDP5 user-guide and 3SEQ documentation defaults verbatim; and the manuscript will adhere to the general principles of the STROBE statement for observational studies where applicable to study design and reporting of inclusion/exclusion criteria.

---

## 3. Sampling Plan

### 3.1 Data sources

| Source | Access | Role |
|---|---|---|
| NCBI Virus / GenBank | Free, Entrez API | Primary sequence + metadata corpus |
| Nextstrain mpox open builds | Free, anonymous S3 (`data.nextstrain.org`) | Cross-check clade assignments; supplementary tree topology |
| Nextclade datasets (mpox, clade-i and all-clades builds, GitHub `nextstrain/nextclade_data`) | Free, public GitHub | Reference genome, clade definitions, QC framework |
| Virological.org public forum (category: MPXV) | Free, public web | Documentation of any community-reported breakpoint hints used as a priori candidate set; will not influence detector calls |

GISAID is **deliberately excluded** to remove terms-of-use risk for an unaffiliated researcher and because the open NCBI corpus is sufficient for clade Ib and the recombinant analysis.

### 3.2 Reference genome and coordinate system

NC_063383 (clade I MPXV reference) is the alignment reference for all coordinate work in this study. All breakpoint coordinates, masked-region coordinates, and reported variant positions are expressed in NC_063383 coordinates.

### 3.3 Corpus construction

A single dated data freeze will be taken from NCBI Virus on or after the OSF registration timestamp. The freeze will be stored with its date and a SHA-256 hash of the concatenated FASTA in the project repository, and the freeze date will be reported in the manuscript.

**Inclusion criteria** (all must be satisfied):
- Organism: *Monkeypox virus*.
- Nextclade-assigned clade I (sub-clade Ia or Ib) OR Nextclade-assigned Ib/IIb recombinant; clade-II-only genomes are retained as outgroup reference only and are not used for the dating analysis.
- Genome length ≥190,000 nucleotides after Nextclade QC.
- ≤10% N bases and ≤5% non-ACGTN bases across the full genome.
- Collection date resolved to at least the calendar month (year-only dates are excluded from the dating analysis; they may be retained for recombination analyses where appropriate).
- Public unrestricted release.

**Exclusion criteria** (any one is sufficient for exclusion, applied independently of the per-base inclusion thresholds above):
- Genomes flagged by Nextclade QC as `bad` for any of: missing data, mixed sites, private mutations exceeding 4 standard deviations of the dataset distribution.
- Genomes with frameshifts in the core gene set as called by Nextclade.
- Genomes whose source publication restricts redistribution.
- Sequencing-platform-suspected contigs (assembled from <30× coverage where coverage metadata is available).

### 3.4 Sample size

The clade Ib + Ib/IIb publicly available corpus is the population, not a sample. At the time of registration the author has not queried the NCBI Virus interface for current counts in order to avoid any pre-analytic look at the data. A pre-specified minimum of 30 date-resolved clade-Ib genomes after the inclusion / exclusion filters of §3.3 is required to proceed with the dating analysis (RQ1). If the freeze contains fewer than 30 date-resolved Ib genomes, RQ1 is dropped and the analysis proceeds with RQ3/RQ4 only; this contingency is pre-specified and is not a post-hoc scope change. The threshold of 30 is chosen to ensure stable maximum-likelihood estimation of the two-parameter Poisson dating model under parametric bootstrap, and is consistent with common practice for viral molecular-clock studies relying on within-clade temporal signal.

### 3.5 Clusters for the dating analysis

Clusters will be defined by Nextstrain's lineage assignments on the clade-i build at the time of the data freeze, with the additional constraint that a cluster must contain ≥10 date-resolved genomes spanning ≥6 months of collection dates to be eligible for independent dating. Singletons and small clusters are pooled into a "minor lineages" group and dated jointly. The list of clusters considered and their cluster boundaries will be locked at the time of the data freeze and reported in the manuscript.

---

## 4. Variables

### 4.1 Sequence-level variables (per genome)

- Accession; collection date (with resolution flag); country; host (reported descriptively — no analysis stratifies on host, since clade-Ib genomes are overwhelmingly human-derived); Nextclade clade; Nextclade lineage; assembly length; N fraction; Nextclade QC overall status; cluster assignment (§3.5).

### 4.2 Variant-level variables

For each variant call on the clade-I alignment relative to NC_063383:
- Reference and alternative allele.
- Tri-nucleotide context (one base 5′, one base 3′).
- APOBEC3 context flag: TRUE if (C→T at TC context) or (G→A at GA context); FALSE otherwise.
- Position; gene; synonymous / non-synonymous flag.

### 4.3 Cluster-level variables

For each cluster (§3.5):
- APOBEC3 SNV count per genome.
- Non-APOBEC3 SNV count per genome.
- Estimated emergence date (median + 95% bootstrap interval) under (a) the context-aware Poisson model with saturation correction, (b) the uncorrected Poisson model, and (c) a TreeTime strict-clock cross-check.
- Per-cluster estimated saturation fraction (proportion of inferred APOBEC3 sites estimated to be multiply hit).

### 4.4 Recombinant-call-level variables

For every genome called by at least one detector (full per-detector record, reported in supplement) and for every consensus candidate (called by ≥2/3 detectors, reported in main text):
- Detector(s) that called it, and which detectors did not.
- Significance value reported by each detector at default thresholds.
- Estimated breakpoint coordinates in NC_063383 reference frame.
- Inferred parental clades / lineages for each segment.
- Concordance flag across alignment and subsampling perturbations (§5.5).
- Consensus-candidate status (TRUE if ≥2/3, FALSE if exactly 1/3).

---

## 5. Analysis Plan

### 5.0 Pipeline order

The analysis pipeline runs as a single deterministic Snakemake DAG in the following order. Every step's output is a frozen artefact that downstream steps depend on; no step is re-litigated after a downstream result.

1. **Fetch and metadata harmonisation.** Pull the data freeze from NCBI Virus; harmonise collection-date resolution flags.
2. **Quality control.** Run Nextclade against the mpox clade-i dataset to obtain clade/lineage assignments, QC flags, and frameshift calls. Apply the inclusion / exclusion filters of §3.3.
3. **Primary alignment.** MAFFT `--auto` against NC_063383; hard-mask ITRs and known hyper-variable regions (§5.1).
4. **APOBEC3 edit counting and H1 test.** Compute the clade-wide APOBEC3 SNV fraction; test H1 (§5.2).
5. **Recombination scan and H3/H4 tests.** Run 3SEQ, RDP5 ensemble, and GARD with frozen default parameters; apply the ≥2/3 consensus rule; verify positive-control recovery; test H3 and H4 (§5.4).
6. **Recombinant exclusion for dating.** Genomes called as consensus candidate recombinants in step 5 (including the positive control) are removed from cluster membership for the subsequent dating analysis.
7. **Cluster definition.** Apply Nextstrain lineage assignments and the cluster-eligibility constraints of §3.5 to the recombinant-excluded genome set; lock the cluster list.
8. **Dating models and H2 test.** Fit the saturation-aware Poisson model (primary) and the uncorrected linearisation (comparator) per locked cluster; compute per-cluster Δ and test H2 (§5.3).
9. **Sensitivity analyses and H5 test.** Re-run steps 3–8 under the perturbations of §5.5; compute Jaccard similarities and test H5.

This ordering resolves the dependency between cluster definition and recombinant identification — the recombination scan completes first, then dating proceeds on a recombinant-free cluster set. No iterative re-fit occurs.

### 5.1 Alignment and masking

- **Primary alignment:** MAFFT v7.5xx `--auto` against NC_063383; thread count `-1` (all available).
- **Masking:** Hard-mask the two inverted terminal repeats (ITRs) and the known mpox hyper-variable regions as specified in the Nextstrain mpox masking BED file at the time of the freeze. The exact BED file (commit hash) will be recorded in the analysis log.
- **Alternative alignment for sensitivity (§5.5):** Nextalign v3 against NC_063383, using the mpox clade-i Nextclade dataset's reference and gene-map. Nextalign is chosen as the alternative aligner because (a) it produces a coordinate-frame-compatible MSA aligned to the same reference as the primary MAFFT alignment, (b) it uses a different alignment algorithm (banded seed-and-extend) than MAFFT's progressive/iterative approach, providing a genuine methodological perturbation, and (c) it is already part of the Nextclade tool stack used in §5.0 step 2. The same masking BED is applied to the Nextalign output.

### 5.2 APOBEC3 edit counting (tests H1)

For each genome, the count of variants matching the APOBEC3 context flag (§4.2) is computed against the masked alignment. The clade-wide APOBEC3 fraction is the sum of APOBEC3 SNVs divided by the sum of all SNVs across the clade Ib corpus. **H1 test:** one-sided one-sample exact binomial test of the clade-wide APOBEC3 fraction against the null p₀ = 0.70 at α = 0.05; H1 is supported if the lower bound of the 95% Wilson interval exceeds 0.70. The observed proportion, the exact-binomial P-value, and the Wilson interval are all reported regardless of outcome.

### 5.3 Saturation-aware Poisson dating model (answers RQ1; tests H2)

Let *L* be the number of eligible APOBEC3 target sites on the masked clade-I reference, defined as the count of positions where the reference and one flanking base together form a TC (5′-T·C-3′) or GA (5′-G·A-3′) dinucleotide context and where neither base is masked. Positions inside a masked region (ITRs, hyper-variable regions per §5.1) are excluded from *L*. Ambiguous reference bases (N, R, Y, W, S, K, M, etc.) and gap characters are excluded from *L*. The exact value of *L*, computed on the masked NC_063383 reference, is reported in the manuscript.

Under the assumption that APOBEC3 edits arise independently at each eligible site at a per-site, per-time rate λ, the expected number of *distinct* edited sites observed in a genome sampled *t* time units after the cluster most-recent-common-ancestor (MRCA) is

E[edited sites at *t*] = L · (1 − exp(−λ · t)).

This is the **saturation-aware** (Poisson-with-context) model. The corresponding **uncorrected model** is its linearisation E[edits at *t*] = L · λ · t (valid as λ·t → 0), which is the form applied by O'Toole et al. (2023) to clade IIb.

Both models are fitted per cluster by maximum likelihood on the per-genome APOBEC3 SNV count, treating the per-genome count as Poisson-distributed around its expectation. The free parameters are the cluster MRCA date *t*₀ and the per-site rate λ. The fitting is implemented in Python with `scipy.optimize.minimize` (method = L-BFGS-B), with bounds *t*₀ ∈ [1990-01-01, freeze-date − 30 days] and λ ∈ [10⁻¹⁰, 10⁻³] per site per day. Optimisation is repeated from 25 random starting points and the best-likelihood fit retained.

Genomes whose collection date is incompatible with the model (collection date < *t*₀ for the candidate fit) are penalised by the likelihood; the optimisation cannot select a *t*₀ later than the earliest collection date in the cluster.

**Cluster emergence date estimate:** the maximum-likelihood estimate of *t*₀, expressed as a calendar date.
**Uncertainty:** parametric bootstrap with 2,000 resamples — for each bootstrap iteration, per-genome APOBEC3 SNV counts are re-drawn from the *fitted model's* expectation. For the saturation-aware primary model, that expectation per genome *i* is *L* · (1 − exp(−λ̂ · (collection_date_*i* − *t̂*₀))); for the uncorrected comparator it is λ̂ · *L* · (collection_date_*i* − *t̂*₀). In each iteration the corresponding model is refit on the resampled counts. The reported 95% interval is the 2.5th–97.5th percentile of the bootstrap distribution of *t*₀.
**Identifiability check:** before reporting any cluster estimate, the profile log-likelihood for *t*₀ is inspected; clusters whose profile is flatter than 2 log-likelihood units across a ±1-year window are flagged as poorly identified, their bootstrap intervals widened accordingly, and they are reported with an explicit identifiability caveat. The criterion is applied symmetrically across all clusters.
**Cross-check:** an independent TreeTime strict-clock run on the same per-cluster masked alignment, reporting its MRCA date with the program's default confidence interval, is reported alongside but is not the primary estimate. Disagreement between TreeTime and the primary estimate beyond mutual-interval overlap is reported and discussed but does not override the pre-specified primary.

**Exclusion of consensus candidate recombinants from dating.** Per the pipeline order (§5.0 steps 5–6), the recombination scan completes before cluster definition. The WHO Dec 2025 / Nextclade 16-Dec-2025 recombinant genome and any additional consensus candidate recombinants identified in §5.4 are removed from cluster membership before the dating models are fit, because their mosaic ancestry violates the constant-rate molecular-clock assumption that the APOBEC3 model relies upon. The list of excluded genomes is reported in the supplement.

**H2 test:** for each cluster, compute Δ = |corrected *t*₀ − uncorrected *t*₀| in days. H2 is supported if max-over-clusters Δ > 60 days. The signed difference and magnitude per cluster are reported regardless, alongside the per-cluster estimated saturation fraction (the fraction of *observed* edited APOBEC3 sites that are multiply hit, derived from the fitted λ and the median per-genome age *t̄* in the cluster as

(1 − exp(−λ·*t̄*) − λ·*t̄*·exp(−λ·*t̄*)) ÷ (1 − exp(−λ·*t̄*))

which is the conditional probability P(site hit ≥ 2 times | site hit ≥ 1 time) under the Poisson model with rate λ over duration *t̄*).

### 5.4 Recombination scan (answers RQ3; tests H3 and H4)

The three pre-specified detectors are run on the masked clade-I + Ib/IIb alignment supplemented with a curated reference outgroup of clade-IIb genomes (composition described below).

**Outgroup panel.** The outgroup is a stratified random sample of 30 high-quality clade-IIb genomes drawn from the data freeze, stratified by collection year (2022, 2023, 2024, 2025+) and by WHO region, with the constraint that each (year × region) stratum contributes at least one and at most 10 genomes. The random seed for stratified sampling is 42. The exact accession list is locked at the data freeze and reported in the supplement.

**Parental reference panels.** For breakpoint polarisation, three reference panels of ≤10 genomes each (clade Ia, clade Ib, clade IIb) are constructed by the same stratified-random procedure at the freeze. The seed is 42 and the accession lists are reported in the supplement.

**Positive-control accession lock.** The WHO Dec 2025 / Nextclade 16-Dec-2025 index recombinant is described qualitatively in this pre-registration but is committed to a specific GenBank accession at the time of the data freeze. The accession is identified by cross-referencing the Nextclade `nextstrain/nextclade_data` mpox-dataset commit dated 16 December 2025 (which reclassified the case as an Ib/IIb recombinant) with the corresponding NCBI Virus record. The locked accession is recorded in the analysis log and reported in the manuscript; all subsequent H3 / H4 tests reference that accession unambiguously.

| Detector | Version policy | Significance threshold | Multiple-comparisons handling |
|---|---|---|---|
| 3SEQ | Latest official release at the freeze date; version recorded in the analysis log | Detector default (P-value with the Dunn–Šidák correction as implemented in 3SEQ) | Internal to 3SEQ |
| RDP5 ensemble (RDP, GENECONV, BootScan, MaxChi, Chimaera, SiScan, 3SEQ-internal) | Latest official release at the freeze date; version recorded in the analysis log | Detector default (P < 0.05 with the program's internal Bonferroni correction); ≥3 of the 7 internal RDP5 methods must agree for an internal "ensemble call" to be reported | Internal to RDP5 |
| GARD (HyPhy) | Latest official release at the freeze date; version recorded in the analysis log | Default ΔAIC threshold for breakpoint acceptance; whole-genome GARD plus a partitioned GARD run using the CDS partitions defined by the Nextclade mpox clade-i gene map | Internal to GARD |

A genome is a **consensus candidate** recombinant if ≥2 of these three detectors call it.

**Ambiguous-base policy.** Positions where ≥10% of corpus genomes carry an N or gap are treated as missing for the purposes of all three detectors and the polarisation step.

**Expected runtime.** On a single modern laptop CPU, 3SEQ on the masked clade-I + Ib/IIb + 30-genome outgroup alignment is expected to complete within minutes; RDP5 ensemble within 1–2 hours; GARD (whole-genome plus gene-partitioned) within 6–48 hours depending on corpus size. Total recombination-scan wall-clock is expected to be in the 12–72 hour range and is documented in the analysis log for the actual run.

**H3 test:** count consensus candidates other than the locked positive-control accession (see "Positive-control accession lock" below). H3 is supported if count ≥1. The reported effect size is the integer count *k*, together with a Wilson 95% interval for the descriptive proportion *k*/*N* (where *N* is the corpus size after exclusions); the H3 test itself is on the count, not the proportion.

**H4 test:** verify that ≥2 of the three detectors call the locked positive-control accession (see "Positive-control accession lock" above) as recombinant. H4 is supported if the case is called by ≥2/3 detectors at default significance, matching the H3 consensus-candidate definition for fairness. The per-detector significance value for the positive control is reported regardless, and the full 3/3 vs 2/3 outcome is reported. Failure (fewer than 2/3) halts H3 reporting until configurations are re-verified per documented defaults; if it still fails, RQ3 is downgraded to exploratory.

**Breakpoint and parental assignment.** For each consensus candidate, breakpoint coordinates (in NC_063383 frame) are taken as the median of the per-detector breakpoint estimates among the calling detectors. Parental clade assignments for each segment are derived by SNV polarisation: for each segment, count the number of clade-Ia, clade-Ib, and clade-IIb diagnostic SNVs from the parental reference panels and assign the segment to the clade with the highest count, provided that count exceeds the next-highest count by ≥3 diagnostic SNVs; segments without a clear majority are flagged as ambiguous.

### 5.5 Sensitivity analyses (tests H5; answers RQ4)

All sensitivity analyses are pre-specified and are run on the same data freeze. Random seeds are fixed at 42 throughout.

1. **Alternate alignment.** Repeat §5.3 and §5.4 on the Nextalign-based alternative alignment (§5.1). Report Jaccard similarity of the consensus-candidate set and the per-cluster shift in emergence date.
2. **Random subsampling.** Repeat §5.3 and §5.4 on five independent random subsamples at 75% of the eligible genomes and five at 50% (10 sub-runs total). For each consensus candidate, report the proportion of the five 75% sub-runs in which it is retained. For each cluster emergence date, report the median across sub-runs with the 5th–95th percentile range.
3. **Temporal hold-out.** Re-fit the cluster dating after excluding the 25% most-recently-collected genomes (which contribute disproportionately to clock-rate identifiability under saturation). Report the per-cluster shift in emergence date. This addresses whether emergence-date estimates are dominated by the most recent samples.
4. **Detector-parameter sweep.** Re-run 3SEQ at default and at ×2 and ÷2 of the default significance penalty; re-run RDP5 with the internal-method agreement threshold lowered to 2/7 and raised to 4/7; re-run GARD with the ΔAIC threshold ±25%. Report the resulting consensus-candidate counts at each setting in a single table. The default setting is the primary analysis; the sweep is purely descriptive.
5. **Saturation-model alternative.** Refit §5.3 with a Jukes–Cantor-style per-site multiply-hit correction in addition to the primary Poisson-context formulation. The exact functional form of the JC-style correction adapted to the two-context APOBEC3 setting (a derivation analogous to the standard `−(3/4)·ln(1 − (4/3)·d)` correction, restricted to the TC/GA target alphabet) is reported in the manuscript Methods; both estimators use the same fitted λ for fairness. Report the per-cluster difference in emergence date and discuss model choice in the limitations.
6. **Outgroup-composition perturbation.** Re-run §5.4 with the clade-IIb outgroup panel replaced by five independent stratified-random redraws (seeds 1–5). Report whether the consensus-candidate set is stable across redraws.
7. **Pipeline-reproducibility verification.** Re-run the complete Snakemake workflow on a second machine (or in a clean conda environment on the same machine) starting from the freeze FASTA + metadata and confirm that the consensus-candidate set is identical to the primary run and that per-cluster emergence-date estimates agree within ±1 day (a tolerance chosen to absorb floating-point non-determinism in `scipy.optimize` across architectures while still detecting any substantive workflow divergence). Any disagreement beyond this tolerance is documented in the analysis log and investigated.

**H5 test:** Jaccard similarity of consensus-candidate sets between primary and (1) is computed and reported with its bootstrap 95% CI; H5 is supported separately for (1) if J ≥ 0.80. The same threshold is applied to the union of consensus calls across the five 75% subsamples in (2). Both decisions are reported regardless of outcome.

### 5.6 Multiple-comparisons handling

H1–H5 are five distinct null hypotheses addressing five distinct quantities; no family-wise correction is applied across them. Within the recombination scan, multiple-comparisons correction is performed internally by each detector at its default setting and is not adjusted further. H2 is a single max-over-clusters statistic and therefore does not require a per-cluster family-wise correction; the per-cluster Δ values are reported descriptively alongside the test outcome.

### 5.7 Pre-specified stopping / pivot conditions

| Trigger | Pre-specified response |
|---|---|
| Corpus contains <30 date-resolved clade-Ib genomes | RQ1 dropped; H1, H2 not tested; analysis proceeds with RQ3/RQ4 only and is reported as a recombination-scan-only study |
| H1 not supported (Wilson lower 95% bound ≤ 0.70) but observed fraction ≥ 0.50 | Report the H1 outcome and observed fraction transparently; proceed with downstream analyses; discuss the lower-than-expected APOBEC3 dominance in the manuscript |
| H1 emergency-halt trigger (observed APOBEC3 fraction <0.50) | Halt downstream analysis; document the diagnostic investigation in the manuscript (misalignment, contamination, or fundamentally different mutational process); do not report cluster dates |
| H4 fails (WHO Dec 2025 case recovered by fewer than 2/3 detectors) | Halt H3 reporting; repeat configurations from documented defaults; if it still fails, downgrade RQ3 to exploratory |
| A peer-reviewed paper publishing the same APOBEC3-dating + recombination-scan combination on clade Ib appears between registration and submission | Re-position the manuscript as an independent verification; cite the prior work; cite the OSF registration DOI and timestamp in the cover letter and manuscript as evidence of independent priority; do not change the pre-specified analysis plan |

---

## 6. Other

### 6.1 Software and reproducibility

- **Language:** Python 3.11.
- **Pinned packages:** biopython, ete3, dendropy, scipy, statsmodels, numpy, pandas. Exact version pins are recorded in the public code repository's conda environment specification and pip requirements files at the time of the freeze.
- **External tools and pinned versions (recorded in the analysis log at the freeze date):** Nextclade CLI v3.x (which includes the Nextalign v3 aligner used in §5.1), MAFFT v7.5xx, IQ-TREE 2.x, TreeTime ≥0.11, 3SEQ (latest official), RDP5 (latest official), HyPhy / GARD (latest official).
- **Pipeline:** A Snakemake workflow will orchestrate the pipeline of §5.0 (fetch → QC → alignment → masking → APOBEC3 counting → recombination scan → recombinant exclusion → cluster definition → dating → sensitivity). The Snakefile and rule modules will be deposited in the public code repository at the time of preprint upload; they do not exist at the time of this registration and are not part of the OSF deposit.
- **Random seeds:** All resampling routines use seed = 42, set at the top of each script.
- **Data freeze:** A single dated NCBI Virus pull plus a Nextstrain build snapshot, both stored with SHA-256 hashes and pull timestamps.
- **Code release:** The full analysis repository (workflow + scripts + environment specification + frozen reference panel manifests) will be deposited in a public GitHub repository at the time of preprint upload, with the same release tagged to a Zenodo DOI.

### 6.2 Data availability

All input sequences and metadata are publicly available from NCBI Virus and Nextstrain. The data freeze used for the analysis (FASTA + metadata + SHA-256 manifest) will be deposited as a public Zenodo dataset at the time of preprint upload.

### 6.3 Known limitations (stated a priori)

1. **GenBank metadata quality.** Collection dates are sometimes year-only; these are excluded from the dating analysis per §3.3. Country and host fields are not validated against primary clinical records.
2. **Recombination detectors are noisy.** The ≥2-of-3 consensus rule is the standard mitigation; single-method-only calls are reported as exploratory (supplement only) and are not part of the H3 test.
3. **Sample-size limits on dating precision.** Small clusters yield wide bootstrap intervals; per-cluster sample sizes will be reported alongside every emergence-date estimate.
4. **Clade Ia is excluded from the dating analysis.** Repeated zoonotic introductions violate the constant-rate molecular-clock assumption that the APOBEC3 model relies upon. Clade Ia genomes appear in the analysis only as parental panel references for the recombination scan.
5. **The corpus is open NCBI only.** Genomes deposited only in GISAID are not analysed. This may shift cluster boundaries and emergence-date estimates relative to a GISAID-inclusive analysis; the limitation is disclosed.
6. **Detector default thresholds are taken as fixed.** Re-tuning detector thresholds against the data would be a degree of post-hoc freedom; the parameter sweep in §5.5 quantifies sensitivity to this choice but does not select among them.
7. **The Poisson dating model assumes per-site independence.** APOBEC3 deaminase activity is known to have penta- and heptamer local-sequence preferences beyond the dinucleotide TC/GA context used here, which can induce correlated edits across nearby sites and violate the strict independence assumed in §5.3. The site-level Poisson model is therefore a tractable approximation rather than a fully mechanistic description; the saturation-model-alternative sensitivity analysis (§5.5 item 5) provides one robustness check, and any cluster whose emergence-date estimate is destabilised by it is flagged in the manuscript.

### 6.4 Deviations from this pre-registration

Any deviation from this protocol will be:
1. Documented in the manuscript Methods with a one-sentence justification.
2. Clearly labelled as exploratory (not confirmatory) in any reported result that depends on the deviation.
3. Compared against the pre-registered analysis as the primary reference.

Foreseeable deviations: exact corpus size will not be known until the freeze; the cluster list (§3.5) is fixed at the freeze and is not re-litigated thereafter; pinned tool versions are determined by the latest official release at the freeze date.

### 6.5 Ethical approval

Not required. The study analyses publicly available viral genome sequences and the associated NCBI / Nextstrain metadata. No human participants, no identifiable patient information, no animal subjects.

### 6.6 Funding

No external funding. All compute is laptop-scale CPU; no paid services are used.

### 6.7 Competing interests

The author declares no competing interests.

### 6.8 Files included in this OSF registration

At the time of registration (2026-05-22), the OSF deposit consists of a single file: this pre-registration document. No data have been downloaded or analysed; no scripts, workflow modules, evaluation manifests, or frozen reference panels yet exist.

The following artefacts will be added to a separate public repository (GitHub plus a Zenodo-DOI-tagged release) at the time of preprint upload, and are listed here so reviewers can verify that they were created post-registration:

- Snakemake workflow (Snakefile and rule modules) implementing the pipeline of §5.0.
- The data freeze: FASTA, metadata table, masking BED file, and SHA-256 manifest.
- The locked outgroup-panel and parental-reference-panel accession lists (§5.4).
- The locked cluster list (§3.5).
- Conda environment specification and pip requirements pinning every dependency.
- All resampling seeds and the random-state log.
- The analysis log recording tool versions, run dates, and runtime per stage.

---

## References

- Aksamentov I, Roemer C, Hodcroft EB, Neher RA. Nextclade: clade assignment, mutation calling and quality control for viral genomes. *J Open Source Softw*. 2021;6(67):3773.
- Alfonsi T, Topcuoglu YS, Chiara M, Bernasconi A. OpenRecombinHunt: Automatic detection of recombination in publicly available viral sequences. *J Mol Biol*. 2026; doi:10.1016/j.jmb.2026.169811.
- Feehley PJ, Feehley MC, Hsieh ZY, Poyer AT, Contreras GP, Yeh TY. A new mathematical model for genetic linkage and recombination of haploid virus populations. *Genetics*. 2026; doi:10.1093/genetics/iyag118.
- Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Mol Biol Evol*. 2013;30(4):772–780.
- Kosakovsky Pond SL, Posada D, Gravenor MB, Woelk CH, Frost SDW. GARD: a genetic algorithm for recombination detection. *Bioinformatics*. 2006;22(24):3096–3098.
- Lam HM, Ratmann O, Boni MF. Improved algorithmic complexity for the 3SEQ recombination detection algorithm. *Mol Biol Evol*. 2018;35(1):247–251.
- Martin DP, Varsani A, Roumagnac P, Botha G, Maslamoney S, Schwab T, Kelz Z, Kumar V, Murrell B. RDP5: a computer program for analyzing recombination in, and removing signals of recombination from, nucleotide sequence datasets. *Virus Evol*. 2021;7(1):veaa087.
- O'Toole Á, Neher RA, Ndodo N, et al. APOBEC3 deaminase editing in mpox virus as evidence for sustained human transmission since at least 2016. *Science*. 2023;382(6670):595–600.
- Sagulenko P, Puller V, Neher RA. TreeTime: maximum-likelihood phylodynamic analysis. *Virus Evol*. 2018;4(1):vex042.
- US Centers for Disease Control and Prevention. Emergence of Clade Ib Monkeypox Virus — Current State of Evidence. *Emerg Infect Dis*. 2025;31(8). Available at wwwnc.cdc.gov/eid/article/31/8/24-1551_article.
- Nextstrain mpox open build documentation, nextstrain.org/mpox.
- WHO Disease Outbreak News. Mpox — Multi-country outbreak: situation update on recombinant lineage. 14 February 2026.
