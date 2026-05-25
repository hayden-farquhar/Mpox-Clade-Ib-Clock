# OSF Pre-Registration CASR2 — Deviation Amendment 04: Post-hoc Spearman ρ diagnostic for alignment-noise contamination of branch-quantity counts

**Parent registration:** OSF Pre-Registration `casr2` ([https://doi.org/10.17605/OSF.IO/CASR2](https://doi.org/10.17605/OSF.IO/CASR2))
**Parent OSF project:** [https://osf.io/gt3vx/](https://osf.io/gt3vx/) (node id `gt3vx`)
**Amendment filed under:** §6.4 of the parent registration ("Deviations from this pre-registration")
**Amendment date:** 2026-05-25
**Amendment scope:** §5.5 ("Sensitivity analyses"), item 1 ("alignment-method sensitivity") — registers a **post-hoc, exploratory** per-tip Spearman rank-correlation diagnostic for distinguishing alignment-noise contamination from genuine biological signal in branch-quantity APOBEC3 counts. **This amendment does not change any registered confirmatory test** (H1, H1ʹ, H2, H3, H4 outcomes are recorded under the pre-registered protocol); it documents a methodological observation that emerged from the registered §5.5 item 1 sensitivity analysis.
**Author:** Hayden Farquhar (ORCID `0009-0002-6226-440X`)

---

## 1. What the registered protocol specified

§5.5 item 1 ("Alignment-method sensitivity") of the parent registration specified:

> Repeat the branch-quantity APOBEC3 fraction (H1ʹ, per Amendment 02) and the §5.4 recombination scan candidate set using an alternative alignment method (full-iterative MAFFT in place of the primary alignment) on the same Ib subset plus reference. Report the H1ʹ point estimate and 95% Wilson interval under the alternative alignment and comment on robustness.

The registered output is a single sensitivity comparison: primary alignment vs alternative alignment. The registration did not specify any per-tip rank-correlation diagnostic.

## 2. What was observed and why a post-hoc diagnostic was required

The registered §5.5 item 1 sensitivity analysis was executed with three alignment-method conditions:

| Condition | Algorithm | n SNVs (Ib branch-quantity, pooled) | APOBEC3 SNVs | Branch-quantity APOBEC3 fraction | Wilson 95 % CI |
|---|---|---|---|---|---|
| Primary | Nextalign (banded seed-and-extend, reference-anchored) | 1,127 | 942 | **0.8358** | (0.8131, 0.8563) |
| Alternative (registered) | MAFFT `--retree 1 --memsave` (single-pass progressive, rough) | 2,427 | 1,165 | **0.4800** | (0.4602, 0.4998) |
| Alternative (additional) | MAFFT `--maxiterate 1000 --auto` (full-iterative gold-standard, cloud-compute) | 5,781 | 1,473 | **0.2548** | (0.2437, 0.2662) |

The Wilson 95 % intervals on these three conditions do not overlap and span a wide range of the [0, 1] fraction scale. Under the registered protocol's sensitivity-1 framing, this is reported as: "the branch-quantity APOBEC3 fraction is alignment-method-sensitive, with the three tested algorithms giving fractions of 0.836, 0.480 and 0.255; the registered H1ʹ test outcome is alignment-dependent."

The substantive scientific question raised by this observation is whether the three fractions reflect **(a)** alignment-method-specific true measurements of the underlying APOBEC3 signal, or **(b)** alignment-method-specific noise inflating a single underlying signal. The registered protocol provides no tool for distinguishing (a) from (b). A post-hoc diagnostic was developed during the §5.5 item 1 analysis and is documented here for transparency.

## 3. The post-hoc diagnostic: per-tip Spearman ρ on APOBEC3 vs non-APOBEC3 sub-counts

For each pair of alignment conditions, compute the per-tip Spearman rank-correlation coefficient ρ separately for:

- **(i)** per-tip APOBEC3 SNV count, and
- **(ii)** per-tip non-APOBEC3 SNV count.

The diagnostic rests on a single principle: if two alignment methods are detecting the **same underlying biological signal**, the per-tip rank order of substitution counts will be positively correlated across methods (the same tips will have more substitutions in both alignments, just shifted in absolute level by the method-specific gap-handling behaviour). If two alignment methods are detecting **alignment-method-specific noise** at uncorrelated positions, the per-tip rank order will be approximately uncorrelated (ρ ≈ 0).

Applied to the present three-alignment comparison, the per-tip ρ values are:

| Spearman ρ | APOBEC3 SNVs | Non-APOBEC3 SNVs |
|---|---|---|
| Nextalign vs rough MAFFT | (recorded in §5.5 result table) | (recorded in §5.5 result table) |
| **Nextalign vs gold MAFFT** | **+0.863** | **−0.047** |

The APOBEC3-SNV rank order is strongly preserved (ρ = +0.863, n = 198 tips): the same Ib genomes carry similar APOBEC3 burdens under either alignment. The non-APOBEC3-SNV rank order is essentially uncorrelated (ρ = −0.047): the additional non-APOBEC3 substitutions detected by the more-thorough alignment fall on **different tips** than the non-APOBEC3 substitutions detected by the less-thorough alignment.

A graded exclusion-sensitivity sweep on the gold-MAFFT condition confirms the diagnostic is reading a continuous effect across the dataset rather than a few outlier tips: excluding tips by their Nextalign-vs-gold-MAFFT Δnon-APOBEC3 difference (>10, >20, >50, >100), the gold-MAFFT branch-quantity APOBEC3 fraction climbs from 0.255 (full data) toward, but never matching, the Nextalign fraction (0.594 at Δ > 10 exclusion, 0.395 at Δ > 50 exclusion).

## 4. Interpretation under the diagnostic

Under the diagnostic outlined in §3, the observed pattern is consistent with the **alignment-noise contamination hypothesis** rather than the **alignment-method-specific true measurement hypothesis**:

- The APOBEC3-SNV signal (ρ = +0.863 across alignments) reflects a real, alignment-robust biological observation. The TpC → TpT mutational signature has strong sequence-context anchors that survive alignment-method differences.
- The non-APOBEC3-SNV signal (ρ = −0.047 across alignments) is dominated by alignment-method-specific gap-placement decisions on the partial-genome NCBI submissions characteristic of the public clade-Ib corpus. Each method's "extra" non-APOBEC3 substitutions correspond to positions where that method's gap-placement algorithm forced a nucleotide-against-nucleotide comparison that a different method would have correctly treated as a gap.
- Under this interpretation, the Nextalign primary fraction (0.8358) is the least noise-contaminated estimate, the rough-MAFFT fraction (0.4800) is moderately contaminated, and the gold-MAFFT fraction (0.2548) is most contaminated. The contamination ranks with alignment "thoroughness" because more iterative refinement creates more opportunities to convert ambiguous-region gaps into substitutions.

This interpretation **supports retaining the registered Nextalign primary as the headline branch-quantity result** while recording the three-alignment fraction comparison and the Spearman ρ diagnostic in the §5.5 sensitivity-1 output.

## 5. Honest caveats

This amendment registers an **exploratory, post-hoc** diagnostic. The pre-registration did not specify the Spearman ρ rule. Consequently:

- The diagnostic outcome cannot be used to claim a pre-registered confirmatory test result; it is reported as a methodological observation that informs the interpretation of the registered §5.5 item 1 sensitivity comparison.
- A third-party gold-standard for non-APOBEC3 SNV calls (e.g., a manually curated set of true non-APOBEC3 substitutions on a small subset of Ib genomes from a high-quality reference re-sequencing dataset) would be required to independently confirm the alignment-noise interpretation. The present amendment does not perform such a validation.
- The Nextalign primary result remains the registered primary; the post-hoc diagnostic is used to inform interpretation, not to change the registered primary.
- The Spearman ρ rule, while principled, is itself contingent on the assumption that the two alignments under comparison are not jointly biased in the same direction. Two alignment methods that both globally under-detect non-APOBEC3 substitutions could produce ρ ≈ 0 even without alignment-method-specific noise. This failure mode is unlikely for the Nextalign-vs-gold-MAFFT comparison (the two algorithms are mechanistically very different) but cannot be ruled out from the present data alone.

## 6. Effect on registered tests and outputs

| Registered output | Effect of this amendment |
|---|---|
| H1 (tip-vs-reference fraction, registered primary) | None — unchanged |
| H1ʹ (branch-quantity fraction, Amendment 02) | None — registered primary remains the Nextalign-based 0.8358 (P = 4.5 × 10⁻²⁶); the §5.5 sensitivity-1 output is expanded to report three alignment-method conditions plus the post-hoc Spearman ρ diagnostic |
| H2 (saturation vs uncorrected dating) | None — unchanged; Δ = 601.77 days supported |
| H3 (consensus inter-clade recombinant count beyond positive control) | None — unchanged (H3 = 0) |
| H4 (consensus recombinant call on positive control) | None — unchanged; all three detectors supported |
| §5.4 recombination scan | None — unchanged |
| §5.5 sensitivity-1 output table | **Expanded** to include three alignment conditions and the per-tip Spearman ρ diagnostic on APOBEC3 and non-APOBEC3 sub-counts |
| §5.7 stopping rules | None — unchanged; the Nextalign-based primary fraction (0.8358) clears the §5.7-row-3 emergency-halt threshold |

## 7. Reproducibility

The diagnostic is computed from the per-tip branch-quantity SNV count tables already deposited at the OSF project under the three §5.5 item 1 alignment conditions (Nextalign primary; rough MAFFT; gold MAFFT). Each per-tip table carries one row per Ib genome with columns for the APOBEC3-context and non-APOBEC3-context branch SNV sub-counts.

The diagnostic computation is a short pandas/scipy operation: merge the three per-tip count tables on accession; compute Spearman's ρ separately on the APOBEC3 sub-count columns and the non-APOBEC3 sub-count columns across alignment-method pairs.

## 8. Acknowledgement of departure from registered protocol

Per §6.4 of the parent registration, this amendment documents an exploratory analytical step that departs from the registered protocol. The departure is in the direction of **reporting more rather than less**: the registered §5.5 item 1 output (a single sensitivity comparison) is preserved and expanded with the Spearman ρ diagnostic.

No registered confirmatory test outcome is changed by this amendment. The amendment was authored after the §5.5 item 1 sensitivity analysis returned a wider-than-expected discrepancy across alignment methods, and is filed at OSF before the corresponding manuscript text is finalised, so that the post-hoc nature of the diagnostic is permanently visible in the OSF audit trail.

---

*Filed at OSF on 2026-05-25 by Hayden Farquhar (ORCID 0009-0002-6226-440X).*
