# OSF Pre-Registration CASR2 — Deviation Amendment 03: Recombination-detector substitution (RDP5 → PhiPack)

**Parent registration:** OSF Pre-Registration `casr2` ([https://doi.org/10.17605/OSF.IO/CASR2](https://doi.org/10.17605/OSF.IO/CASR2))
**Parent OSF project:** [https://osf.io/gt3vx/](https://osf.io/gt3vx/) (node id `gt3vx`)
**Amendment filed under:** §6.4 of the parent registration ("Deviations from this pre-registration")
**Amendment date:** 2026-05-23
**Amendment scope:** §5.4 ("Recombination scan"), §6.1 ("Software and reproducibility") — substitutes PhiPack (Bruen, Philippe & Bryant 2006, *Genetics* 172:2665) for the registered RDP5 ensemble as the third detector in the ≥2/3 consensus rule. The consensus rule is preserved at ≥2/3; only the identity of the third detector changes.
**Author:** Hayden Farquhar (ORCID `0009-0002-6226-440X`)

---

## 1. What the registered protocol specified

The parent registration's §5.4 ("Recombination scan") specifies three detectors:

- **3SEQ** — exhaustive triplet exact-binomial test (Lam, Ratmann & Boni 2018);
- **RDP5 ensemble** — seven internal methods (RDP, GENECONV, BootScan, MaxChi, Chimaera, SiScan, 3SEQ-internal) with the internal-ensemble call defined as ≥3-of-7 internal methods agreeing at default thresholds (Martin et al. 2021, *Virus Evolution* 7:veaa087);
- **GARD** — phylogenetic-partition + cAIC breakpoint detection (Kosakovsky Pond et al. 2006, *Mol Biol Evol* 23:1891).

A genome is a "consensus candidate recombinant" if ≥2 of these three detectors call it at default thresholds. §6.1 lists RDP5 in the pinned-software table with the version policy "Latest official release at the freeze date; version recorded in the analysis log."

## 2. What changed and why

**RDP5 is no longer publicly distributed from its registered author's primary distribution site as of the data freeze date (2026-05-22).** Specifically:

- The host previously distributing the RDP5 installer (a personal-domain academic site at the University of Cape Town) does not resolve. Last cached snapshot in the Internet Archive Wayback Machine is dated 2025-09-04; the host has been offline thereafter.
- The Wayback Machine does not preserve the installer binary itself (large `.exe` files are not part of standard web-archive captures).
- A systematic search across GitHub repositories, SourceForge, Bitbucket, Bioconda, conda-forge, and anaconda.org personal channels returned no public mirror of the RDP5 installer or source code at the freeze date.
- The registered author remains contactable by email, but direct-request distribution does not satisfy the "Latest official release at the freeze date" requirement of §6.1, which presumes a stable, citation-able distribution URL.

A secondary infrastructure constraint: the Windows-only RDP5 binary requires the Wine compatibility layer to run on macOS. Wine on macOS was formally deprecated by Homebrew on 2026-05-21 with a sunset date of 2026-09-01 (the cask declares: "Deprecated because it does not pass the macOS Gatekeeper check"). This deprecation does not block running RDP5 *now* but means any reproducibility verification more than three months after the freeze date would face a vanishing compatibility layer in addition to the missing installer.

Together these constraints make the registered RDP5 component unrecoverable through standard channels at the freeze date. A substitution under §6.4 is therefore required before the recombination scan can complete with the registered ≥2/3 consensus rule intact.

## 3. The substitute: PhiPack

PhiPack (Bruen, Philippe & Bryant 2006, *Genetics* 172:2665, >2,000 citations) is selected as the replacement third detector. Reasons:

1. **Methodologically orthogonal to both 3SEQ and GARD.** 3SEQ tests individual triplets via an exact-binomial framework on informative-site patterns. GARD partitions the alignment and compares phylogenies on partitions via cAIC. PhiPack implements a third paradigm — population-level pairwise-homoplasy / incompatibility tests — that does not depend on either triplet enumeration or phylogenetic partitioning. The three detectors therefore retain the orthogonality property that motivated the registered ≥2/3 consensus rule.

2. **Open-source, multi-platform, peer-reviewed.** Available as the `phipack` package in the Bioconda channel for `osx-64` and `linux-64` at the freeze date. Single command-line binary, no GUI dependency, no Wine layer.

3. **Three internal statistics analogous to RDP5's internal ensemble.** PhiPack reports Φ (pairwise homoplasy index), NSS (neighbour similarity score), and Max Chi² in a single run. This is structurally similar to the registered RDP5 internal-ensemble construction, allowing a directly comparable "internal-ensemble call" definition (see §5 below).

4. **Standard alternative to RDP5 in the published recombination-detection literature.** PhiPack is the most widely used non-RDP5 single-paper recombination test in viral phylogenetics and is the de-facto third detector in most published consensus-rule applications.

PhiPack's primary statistic is Φ. The Bruen et al. 2006 paper presents Φ as a robust test of the null hypothesis of clonal evolution (no recombination) against the alternative of any recombination signal. NSS and Max Chi² are reported alongside Φ as cross-validating statistics drawing on different aspects of the homoplasy structure.

## 4. Updated software-pinning entry

§6.1's pinned-software table is updated to replace the RDP5 entry with:

| Tool | Version policy |
|---|---|
| PhiPack | Latest stable release in the Bioconda `osx-64` channel at the freeze date; version and Bioconda build string recorded in the analysis log |

The 3SEQ and GARD entries are unchanged. All other pinned tools (Nextclade, MAFFT, IQ-TREE, TreeTime, HyPhy) are unchanged.

## 5. Updated §5.4 detector composition and consensus rule

The detector table in §5.4 is updated:

| Detector | Version policy | Significance threshold | Multiple-comparisons handling |
|---|---|---|---|
| 3SEQ | Latest official release at the freeze date; version recorded in the analysis log | Detector default (P-value with the Dunn–Šidák correction as implemented in 3SEQ) | Internal to 3SEQ |
| **PhiPack** | Latest stable release in the Bioconda `osx-64` channel at the freeze date; version + Bioconda build string recorded in the analysis log | Detector default: Φ test rejects the null of clonal evolution at P < 0.05 (the Bruen et al. 2006 default) on the candidate-plus-parental-panel subset alignment, with NSS and Max Chi² reported alongside | Internal to PhiPack |
| GARD (HyPhy) | Latest official release at the freeze date; version recorded in the analysis log | Default ΔAIC threshold for breakpoint acceptance; whole-genome GARD plus a partitioned GARD run using the CDS partitions defined by the Nextclade mpox clade-i gene map | Internal to GARD |

**A genome is a consensus candidate recombinant if ≥2 of these three detectors call it at default thresholds.** This preserves the registered ≥2/3 consensus rule unchanged.

### Per-candidate operationalisation of PhiPack calls

PhiPack's primary statistic (Φ) tests the alignment as a whole rather than producing per-sequence calls in the manner of 3SEQ or RDP5. To produce per-sequence calls compatible with the registered ≥2/3 consensus rule, PhiPack is operationalised as follows:

- For each candidate sequence (defined as any sequence called by at least one of 3SEQ or GARD at default thresholds), PhiPack is run on a subset alignment containing: the candidate, the three parental reference panels (≤10 clade-Ia, ≤10 clade-Ib, ≤10 clade-IIb), and the locked positive-control accession (also serving as a "known-recombinant" reference). The candidate is flagged as a PhiPack call if Φ rejects the null at P < 0.05 on this subset alignment.
- NSS and Max Chi² P-values are reported alongside Φ for every per-candidate run, regardless of outcome.
- An "internal-ensemble call" is defined as: ≥2 of PhiPack's 3 internal statistics (Φ, NSS, Max Chi²) reject the null at default thresholds on the candidate-plus-parental-panel subset alignment. The ≥2-of-3 threshold approximates the registered RDP5 ≥3-of-7 ratio (~43%) while accounting for PhiPack's smaller internal-method count. This is the threshold used for the consensus-rule call.

This operationalisation maps PhiPack onto the per-sequence call format that the registered consensus rule assumes, while preserving PhiPack's default Φ + NSS + Max Chi² statistics at their published thresholds.

## 6. Effect on the registered hypotheses

- **H3 (count of consensus candidates other than the positive control).** Unchanged in statement and threshold; the candidate set is determined by the ≥2/3 consensus across the new three-detector composition. H3 is supported if count ≥1.
- **H4 (recovery of the positive control).** Unchanged in statement and threshold; the positive control must be called by ≥2 of the three detectors (now 3SEQ + GARD + PhiPack) at default significance.
- **H5 (alignment-perturbation sensitivity; Jaccard ≥0.80).** Unchanged.

The §5.7 stopping conditions are unchanged. Failure of H4 (positive control recovered by fewer than 2 of the three detectors) still halts H3 reporting.

## 7. Optional later addition (clearly out of scope of this amendment)

If the RDP5 author makes the installer available through a stable distribution channel before manuscript submission, the registered protocol may be additionally run with RDP5 as a fourth, supplementary detector reported only as a §5.5 sensitivity analysis. The registered ≥2/3 consensus rule and the H3 / H4 outcomes would continue to be determined by 3SEQ + GARD + PhiPack as specified in this amendment; RDP5 sensitivity would be a transparency check, not a re-test. This is named as an option only to clarify what would and would not constitute an amendment to the registered protocol; no commitment is made here.

## 8. Scope of this amendment

This amendment changes only the third detector in the recombination scan. No other element of the registered analysis plan is modified:

- The reference genome (NC_063383), the masking BED (corrected under Amendment 01 at [https://osf.io/gt3vx/files/osfstorage/6a1049bdd97d989484f9305c](https://osf.io/gt3vx/files/osfstorage/6a1049bdd97d989484f9305c)), and the §5.4 alignment composition are unchanged.
- 3SEQ and GARD detector configurations are unchanged.
- The outgroup panel (30 stratified clade-IIb genomes) and the three parental reference panels (≤10 each clade-Ia, clade-Ib, clade-IIb) are unchanged.
- The positive-control accession lock is unchanged.
- The H1ʹ branch-quantity exploratory analysis added under Amendment 02 at [https://osf.io/gt3vx/files/osfstorage/6a1049c0964dd2d18df92ec0](https://osf.io/gt3vx/files/osfstorage/6a1049c0964dd2d18df92ec0) is unchanged.
- The §5.3 saturation-aware dating model and the §5.5 sensitivity analyses are unchanged.

## 9. Priority statement

This amendment is filed before any PhiPack run on the §5.4 alignment has been executed. The 3SEQ result (described and deposited in the Phase 3 results document at [https://osf.io/gt3vx/files/osfstorage/6a10ce8e07dad27cfdf6c2e2](https://osf.io/gt3vx/files/osfstorage/6a10ce8e07dad27cfdf6c2e2)) was generated under the as-registered protocol before this amendment was drafted; the GARD run is also underway under the as-registered protocol at amendment-drafting time. PhiPack's role in the consensus rule is defined here in advance of its first invocation on the analytical alignment, satisfying the standard pre-specification requirement for any detector added or substituted via the §6.4 mechanism.
