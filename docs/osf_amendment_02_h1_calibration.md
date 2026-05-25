# OSF Pre-Registration CASR2 — Deviation Amendment 02: H1 calibration framing

**Parent registration:** OSF Pre-Registration `casr2` ([https://doi.org/10.17605/OSF.IO/CASR2](https://doi.org/10.17605/OSF.IO/CASR2))
**Parent OSF project:** [https://osf.io/gt3vx/](https://osf.io/gt3vx/) (node id `gt3vx`)
**Amendment filed under:** §6.4 of the parent registration ("Deviations from this pre-registration")
**Amendment date:** 2026-05-22
**Amendment scope:** §2.3 ("Exploratory analyses") and §5.2 ("APOBEC3 edit counting (tests H1)") — adds a registered exploratory analysis with a biologically aligned framing of H1. The original §5.2 H1 test is retained and reported transparently; this amendment does not retract it.
**Author:** Hayden Farquhar (ORCID `0009-0002-6226-440X`)

---

## 1. What triggered this amendment

The registered H1 test (§5.2) was executed against the 2026-05-22 freeze Phase 3a outputs and triggered the §5.7-row-3 emergency-halt condition (observed APOBEC3 fraction < 0.50). Headline result:

| Quantity | Value |
|---|---|
| n_genomes_in_test (clade Ib only) | 198 |
| Total SNVs (Ib pooled) | 20,698 |
| APOBEC3 SNVs (Ib pooled) | 2,494 |
| Observed APOBEC3 fraction | 0.1205 |
| Wilson 95% CI | (0.1161, 0.1250) |
| Registered null p₀ | 0.70 |
| One-sided exact-binomial P | 1.0000 |
| §5.7-row-3 halt threshold | 0.50 |

Under §5.7-row-3, the halt requires diagnostic investigation before any downstream analysis can proceed. The diagnostic (RUNLOG, 2026-05-22 22:00–22:30 UTC) is summarised below; this amendment is the documented corrective action.

## 2. What the diagnostic found

The 12.05% observed APOBEC3 fraction is **not** evidence of contamination, misalignment, or a fundamentally different mutational process in clade Ib. It is the consequence of a calibration mismatch between the registered counting protocol and the registered null:

1. **Registered counting protocol (§4.2 + §5.2):** SNVs are computed per genome relative to the masked NC_063383 / DQ011155.1 reference. NC_063383 / DQ011155.1 is a 1979 clade-Ia isolate (Zaire-1979-005). Every SNV the tip carries against that reference is counted, including SNVs that arose during the deep clade-Ib-vs-clade-Ia divergence (pre-2022 zoonotic-era evolution).

2. **Registered null p₀ = 0.70:** Taken from O'Toole et al. (*Science* 2023). That figure is the APOBEC3 fraction along the **post-zoonotic clade-IIb human-transmission branch** — i.e., a per-branch quantity calculated on mutations that arose since the clade-IIb MRCA, restricted to the human-to-human expansion phase.

3. **Mismatch.** For clade-Ib tips counted vs a clade-Ia reference, the per-tip total-SNV denominator is dominated by clade-Ib-vs-clade-Ia divergence accumulated in the zoonotic phase (mostly random base substitutions, only marginally APOBEC3-enriched). The clock-era (post-2022 Kamituga-expansion) APOBEC3 edits in the numerator sit on top of that floor. The dilution drives the registered tip-vs-reference fraction below the per-branch null. The observed 12% Ib fraction is what is expected given (a) ~12 transmission-acquired APOBEC3 edits per clade-Ib tip plus (b) ~90 deep-divergence SNVs in the denominator that are not part of the recent transmission process the O'Toole null is calibrated to.

## 3. Visual confirmation

The diagnostic scatter at `outputs/figures/snv_by_date.png` plots APOBEC3-context SNV count per genome against collection date for both clades. Two qualitative observations are diagnostic of the dilution interpretation:

- **Clade Ia (n=694, 1971–2025):** flat 2–10 APOBEC3 SNVs per genome with no time trend. This is the noise floor of clade-Ia-vs-clade-Ia-reference background — exactly what is expected when counting tips against a reference of the same sub-clade.
- **Clade Ib (n=198, 2022–2026):** clear positive time trend, from ~7 APOBEC3 SNVs in 2022 to ~22 SNVs in 2026. This is the O'Toole-style accumulation signature, intact.

The biological clock signature is present in the data. Only the H1 statistical-test framing was misaligned.

## 4. Corrective action under §6.4 + §2.3

This amendment **adds a registered exploratory analysis** to the parent registration, under §2.3 ("Exploratory analyses"), and **retains the original §5.2 H1 test** for transparent reporting alongside the new analysis. The new analysis is committed in advance of execution.

### 4.1 New registered exploratory analysis: branch-quantity APOBEC3 fraction (H1ʹ)

**Statement.** Compute the APOBEC3 fraction along terminal branches from a reconstructed clade-Ib most-recent-common-ancestor (MRCA), and test that **branch-quantity** APOBEC3 fraction against the same null p₀ = 0.70 that O'Toole et al. (2023) calibrated for the per-branch clade-IIb quantity.

**Procedure.**

1. Subset the masked primary MAFFT alignment to the 198 clade-Ib genomes.
2. Infer a maximum-likelihood tree on the clade-Ib subset with IQ-TREE (substitution model selected by ModelFinder; bootstrap is not required because no internal-node interval is reported, only the root state).
3. Reconstruct the ancestral sequence at the root (clade-Ib MRCA) with IQ-TREE's `--ancestral` option.
4. For each clade-Ib tip, compute the count of substitutions along its **terminal branch from the reconstructed MRCA to the tip**, classified as APOBEC3-context (TC→TT or GA→AA in the reconstructed MRCA context) or non-APOBEC3.
5. Pool the per-tip branch-quantity counts to obtain the clade-Ib **branch-quantity APOBEC3 fraction**.
6. Test against p₀ = 0.70 with a one-sided exact binomial test at α = 0.05; the new test is denoted H1ʹ and is supported if the lower bound of the 95% Wilson interval exceeds 0.70.

**Stopping rule for H1ʹ.** The §5.7-row-3 emergency-halt threshold (observed < 0.50) applies to H1ʹ as well, with the same diagnostic-investigation requirement before downstream analyses proceed. The H1ʹ stopping rule is independent of the original §5.2 H1 stopping rule; the registered §5.2 stopping rule does not need to be re-fired by H1ʹ.

**Implementation artefact (committed in advance).** A new script `scripts/build_branch_apobec3_counts.py` will be written before execution; the script will be deposited in the public code repository at the time of preprint upload, alongside the existing `scripts/count_apobec3.py` (which implements the tip-vs-reference counting of the registered §5.2 protocol).

### 4.2 Reporting plan for the original §5.2 H1 test

The original H1 test result obtained at the 2026-05-22 freeze is reported in the manuscript Methods + Results sections with the following explicit framing:

- The exact pre-registered statistical statement, point estimate, Wilson CI, and one-sided P-value as already computed.
- The §5.7-row-3 emergency halt and the diagnostic investigation that followed.
- The calibration mismatch identified by the diagnostic.
- This amendment as the corrective action and the H1ʹ result as the biologically interpretable comparator.

This satisfies the §5.7-row-3 protocol's requirement to "document the diagnostic investigation in the manuscript", and satisfies §6.4's general requirement to disclose every deviation transparently.

### 4.3 Lifting of the §5.7-row-3 halt

The §5.7-row-3 halt blocks downstream cluster dating ("do not report cluster dates"). On the strength of:

- the diagnostic establishing that the halt-trigger is a registered-protocol calibration mismatch rather than a biological signal corruption,
- the scatter-plot visual confirmation that the APOBEC3 clock signature is intact in clade Ib,
- the addition of H1ʹ as a biologically aligned hypothesis test that the data are expected to support, and
- the mask correction in Amendment 01,

the §5.7-row-3 halt is **lifted conditional on**:

1. Filing of this amendment and Amendment 01 at OSF before downstream re-entry (this document is the filing).
2. Successful execution of H1ʹ on the corrected Phase 2 / Phase 3 outputs and reporting of its result transparently regardless of outcome.
3. If H1ʹ itself triggers the §5.7-row-3 threshold (observed branch-quantity APOBEC3 fraction < 0.50 in clade Ib), the halt re-fires and downstream cluster dating is not reported.

## 5. Scope of this amendment

This amendment introduces one new registered exploratory hypothesis (H1ʹ) and one new script (`scripts/build_branch_apobec3_counts.py`) under §2.3 + §5.2. No other element of the registered analysis plan is modified:

- The reference genome, alignment method, masking protocol (corrected under Amendment 01), and clade-Ib genome set are unchanged.
- The recombination-scan protocol (§5.4) and its H3/H4 hypotheses are unchanged.
- The dating-model protocol (§5.3) and its H2 hypothesis are unchanged. H2 remains tested on the registered counting; H1ʹ does not propagate into the H2 saturation model, which continues to use the per-genome counts as specified in §5.3.
- The sensitivity-analyses protocol (§5.5) and its H5 hypothesis are unchanged.

## 6. Priority statement

This amendment is filed before any of the new analyses it specifies are executed. The diagnostic that motivated it (RUNLOG 2026-05-22 22:00–22:30 UTC) was performed under the registered §5.7-row-3 protocol, which mandates diagnostic investigation when the halt threshold is crossed. No analysis outside that registered diagnostic has been run. The branch-quantity counting, the H1ʹ test, the Phase 2 re-run on the corrected mask, and the Phase 3 re-run on the MAFFT primary alignment are all to-be-executed under the corrected protocol after this amendment is filed.
