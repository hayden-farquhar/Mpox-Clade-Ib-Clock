# OSF Pre-Registration CASR2 — Deviation Amendment 01: Masking BED source

**Parent registration:** OSF Pre-Registration `casr2` ([https://doi.org/10.17605/OSF.IO/CASR2](https://doi.org/10.17605/OSF.IO/CASR2))
**Parent OSF project:** [https://osf.io/gt3vx/](https://osf.io/gt3vx/) (node id `gt3vx`)
**Amendment filed under:** §6.4 of the parent registration ("Deviations from this pre-registration")
**Amendment date:** 2026-05-22
**Amendment scope:** §5.1 ("Alignment and masking") — the masking BED source.
**Author:** Hayden Farquhar (ORCID `0009-0002-6226-440X`)

---

## 1. What the registered protocol specified

Pre-registration §5.1 states:

> "Masking: Hard-mask the two inverted terminal repeats (ITRs) and the known mpox hyper-variable regions as specified in the Nextstrain mpox masking BED file at the time of the freeze. The exact BED file (commit hash) will be recorded in the analysis log."

## 2. What was actually executed at the 2026-05-22 freeze

The Snakemake `resolve_mask_bed` rule used a three-priority resolution: (1) user-supplied `MASK_BED` config, (2) `<nextclade_dataset>/mask_overview.bed`, (3) built-in fallback masking only the ITRs.

At the 2026-05-22 freeze:

- No user-supplied `MASK_BED` was set in `workflow/config.yaml`.
- The `nextstrain/mpox/clade-i` Nextclade dataset distributed via `nextclade dataset get` did **not** ship a `mask_overview.bed`. Verified by directory listing: the dataset contains `reference.fasta`, `sequences.fasta`, `pathogen.json`, `genome_annotation.gff3`, `tree.json`, `CHANGELOG.md`, `README.md` — no BED file.
- Consequently the `resolve_mask_bed` rule fell through to the built-in fallback, which produced an ITR-only mask (`DQ011155.1 0 6389 ITR_5prime` + `DQ011155.1 190578 196967 ITR_3prime`).

The hyper-variable-region component of the §5.1 mask specification was therefore silently omitted from the Phase 2 masked alignment used as input to Phase 3a APOBEC3 counting.

## 3. Why this matters

Hyper-variable regions and known indel-prone loci contribute disproportionately to per-genome SNV counts without carrying APOBEC3 clock signal. Omitting them from the mask inflates the `total_snvs` denominator in the H1 test and dilutes the APOBEC3 fraction. The mask deviation is a contributor — though not the primary cause — of the H1 emergency-halt triggered at the 2026-05-22 22:00 UTC run (see also Amendment 02). Restoring the registered hyper-variable mask is required before any downstream analysis can be re-entered under §5.7.

## 4. Corrective action taken

A new masking BED has been sourced directly from the upstream `nextstrain/mpox` GitHub repository, pinned to commit `1ef5ef2ae1f92cf07dc659c1dfb05e6d5b7abdf1` (master tip at fetch time `2026-05-22T12:14:28Z`).

Four candidate masks were fetched from the repository for comparison:

| File | Lines | SHA-256 (truncated) |
|---|---|---|
| `nextclade/resources/clade-i/mask.bed` | 33 | `3dcd9187…3ac` |
| `phylogenetic/defaults/clade-i/mask.bed` | 25 | `f75d944a…51c` |
| `phylogenetic/defaults/mask.bed` | 11 | `5a4acac3…01d` |
| `phylogenetic/defaults/mask_overview.bed` | 16 | `9b98f5cf…a30` |

The clade-I-specific Nextclade mask (`nextclade/resources/clade-i/mask.bed`, SHA-256 `3dcd918726959d56617718b7631fcaf99347e4b182da150bb6051c22c43443ac`) is selected. Reasons:

1. **Same dataset family as the rest of the pipeline.** The Phase 1 QC ran through the `nextstrain/mpox/clade-i` Nextclade dataset; the matching mask is sourced from the same dataset's source tree.
2. **Most comprehensive coverage.** 33 intervals vs 11–25 in the alternatives.
3. **Both ITRs covered as named intervals** (`0–1500 mask from beginning` and `190885–196967 mask ITR from end`) in the same coordinate frame as the rest of the mask.
4. **Explicit Ib-specific homoplasy coverage.** 15 of the 33 intervals are tagged `homoplasic in Ib` — sites with parallel substitutions across Ib lineages that inflate non-APOBEC3 background SNV counts when included.

All four candidate masks are retained in the freeze under `data/raw/freeze_20260522/nextstrain_mask/` with their SHA-256 hashes and a `PROVENANCE.md` recording the choice rationale. The chosen mask is committed as `CHOSEN_mask.bed`. The Snakemake `workflow/config.yaml` has been updated:

```yaml
MASK_BED: data/raw/freeze_20260522/nextstrain_mask/CHOSEN_mask.bed
```

This causes the existing `resolve_mask_bed` rule's first-priority branch to fire, bypassing the Nextclade-dataset and built-in-fallback paths.

## 5. Coordinate-frame note

The upstream BED uses generic `chr` as the chromosome name. `scripts/apply_mask.py` is chrom-agnostic (uses only `start, end` columns), so no rewrite of the chrom column is required.

The upstream intervals are in reference-coordinate (NC_063383 / DQ011155.1 base-pair) space. `scripts/apply_mask.py` applies intervals to alignment-column positions directly. For the 2026-05-22 freeze's MAFFT primary alignment, the reference row carries 109 gap characters out of 196,967 columns despite the `--keeplength` invocation; this artefact is being investigated at re-run. If the corrected re-run confirms a gap-free reference row, reference-coordinate → alignment-column mapping is 1:1 and no translation is needed. If gap chars persist, a translation utility will be added that builds the mapping from the alignment's reference row before `apply_mask.py` is invoked. Either path will be documented in the analysis log at re-run time.

## 6. Effect on existing Phase 2 and Phase 3a artefacts

The Phase 2 masked alignments (`alignment_mafft.masked.fasta`, `alignment_nextalign.masked.fasta`) and Phase 3a APOBEC3 counts (`apobec3_counts_nextalign.tsv`, `h1_result_nextalign.json`) produced under the ITR-only mask are **superseded** by this amendment. They are retained on disk for audit and reproducibility but are not the analytical artefacts that will be reported in the manuscript. Phase 2 will be re-run with the upstream mask and Phase 3 re-run on the corrected alignment.

## 7. Scope of this amendment

This amendment changes only the masking BED source. No other element of the registered analysis plan is modified:

- The reference genome remains NC_063383 / DQ011155.1.
- The MAFFT `--auto` primary aligner and Nextalign sensitivity aligner are unchanged.
- The §3.3 inclusion/exclusion filters and the §3.4 sample-size stopping rule are unchanged.
- The hypothesis tests H1–H5 and the §5.7 stopping conditions are unchanged (Amendment 02 addresses a separate calibration issue with H1).

## 8. Priority statement

This amendment is filed before any downstream analysis (Phase 3 re-run, recombination scan, dating model fitting) has been executed against the corrected mask. The 2026-05-22 22:00 UTC H1 emergency-halt diagnostic — which identified the mask deviation and is the basis for this amendment — was conducted under the registered §5.7-row-3 protocol. No exploratory analyses outside the registered plan have been run; the diagnostic was purely descriptive (counts, scatter plot, narrative interpretation).
