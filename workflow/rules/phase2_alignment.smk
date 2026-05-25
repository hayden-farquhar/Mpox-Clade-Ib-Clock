"""
Phase 2 rules — alignment and masking.

Maps to pre-registration §5.0 step 3 and §5.1.

Two alignments are produced from the §3.3 filtered corpus:
- Primary: MAFFT --auto against NC_063383.
- Alternative: Nextalign v3 against the same reference (for the §5.5 item-1 sensitivity).

Both are masked with the same BED file. The masking BED can be either
(a) supplied by the user via the MASK_BED config option, or (b) auto-derived
from the Nextstrain mpox masking definitions if available in the Nextclade
dataset directory.

Inputs (from Phase 1)
---------------------
- data/processed/freeze_{DATE}/corpus.fasta  — §3.3-filtered analysis corpus
- data/raw/freeze_{DATE}/nextclade_dataset/   — contains the reference and any
                                                 dataset-supplied mask BED

Outputs
-------
- data/interim/freeze_{DATE}/alignment_mafft.fasta         — MAFFT --auto MSA
- data/interim/freeze_{DATE}/alignment_mafft.masked.fasta  — MAFFT MSA with mask applied
- data/interim/freeze_{DATE}/alignment_nextalign.fasta     — Nextalign MSA
- data/interim/freeze_{DATE}/alignment_nextalign.masked.fasta
- data/interim/freeze_{DATE}/mask.bed                       — the mask actually used
"""

# Determine the mask BED source.
# Priority: explicit MASK_BED config → Nextclade dataset-supplied mask → built-in default.
MASK_BED_SOURCE = config.get("MASK_BED")

# ============================================================================
# 1. Reference extraction
# ============================================================================

rule extract_reference:
    """Pull the clade-I reference (NC_063383) out of the Nextclade dataset."""
    input:
        dataset = f"{RAW_DIR}/nextclade_dataset",
    output:
        reference = f"{INTERIM_DIR}/reference.fasta",
    log:
        f"{INTERIM_DIR}/.extract_reference.log",
    shell:
        r"""
        mkdir -p $(dirname {output.reference})
        # Nextclade dataset reference is always at <dataset>/reference.fasta
        cp {input.dataset}/reference.fasta {output.reference}
        echo "[$(date -u +%FT%TZ)] reference: $(grep '^>' {output.reference})" | tee {log}
        """

# ============================================================================
# 2. Mask BED resolution
# ============================================================================

rule resolve_mask_bed:
    """
    Resolve the masking BED to use. Order:
      1. User-supplied MASK_BED config (absolute or repo-relative path).
      2. <nextclade_dataset>/mask_overview.bed if present.
      3. Built-in fallback: mask only the ITRs (positions auto-derived from
         the reference; see scripts/build_default_mask.py).
    """
    input:
        dataset = f"{RAW_DIR}/nextclade_dataset",
        reference = f"{INTERIM_DIR}/reference.fasta",
    output:
        bed = f"{INTERIM_DIR}/mask.bed",
    params:
        user_supplied = MASK_BED_SOURCE,
    log:
        f"{INTERIM_DIR}/.resolve_mask_bed.log",
    shell:
        r"""
        if [ -n "{params.user_supplied}" ] && [ -f "{params.user_supplied}" ]; then
            cp "{params.user_supplied}" {output.bed}
            echo "[mask] using user-supplied BED: {params.user_supplied}" | tee {log}
        elif [ -f "{input.dataset}/mask_overview.bed" ]; then
            cp "{input.dataset}/mask_overview.bed" {output.bed}
            echo "[mask] using Nextclade dataset mask_overview.bed" | tee {log}
        else
            python3 scripts/build_default_mask.py \
                --reference {input.reference} \
                --out {output.bed} 2>&1 | tee {log}
            echo "[mask] using built-in default (ITRs only) — review before publication" | tee -a {log}
        fi
        echo "[mask] $(wc -l < {output.bed}) intervals masked" | tee -a {log}
        """

# ============================================================================
# 3. MAFFT --auto primary alignment
# ============================================================================

rule align_mafft:
    """
    Primary alignment: MAFFT `--auto --keeplength --add` against NC_063383.

    Methodological clarification (documented in RUNLOG): the pre-registration
    §5.1 specifies "MAFFT --auto against NC_063383". On the laptop-CPU compute
    target, the literal interpretation (all-vs-all progressive alignment via
    `mafft --auto`) was found infeasible at the actual freeze size (893
    sequences × ~190 kb): the process exceeded 90 CPU-min without producing
    output and triggered memory pressure on a 16 GB machine. The reference-
    anchored MAFFT idiom (`--auto --keeplength --add`) is used instead. This:
      - still invokes `--auto` (MAFFT still picks the progressive strategy);
      - still uses NC_063383 as the anchor reference;
      - produces an alignment where every column maps cleanly to a reference
        coordinate (which is what the pre-reg's coordinate-system declaration
        in §3.2 requires);
      - completes in minutes instead of hours on the same hardware.
    """
    input:
        corpus    = f"{PROCESSED_DIR}/corpus.fasta",
        reference = f"{INTERIM_DIR}/reference.fasta",
    output:
        aln = f"{INTERIM_DIR}/alignment_mafft.fasta",
    threads:
        workflow.cores
    log:
        f"{INTERIM_DIR}/.align_mafft.log",
    shell:
        r"""
        mafft --auto --keeplength --thread {threads} \
            --add {input.corpus} \
            {input.reference} > {output.aln} 2> {log}
        echo "[mafft] alignment length: $(awk '/^>/{{next}} {{print length; exit}}' {output.aln})" | tee -a {log}
        echo "[mafft] records: $(grep -c '^>' {output.aln})" | tee -a {log}
        """

# ============================================================================
# 4. Nextalign alternative alignment (for §5.5 sensitivity)
# ============================================================================

rule align_nextalign:
    """Alternative alignment: Nextalign v3 (part of nextclade CLI v3)."""
    input:
        corpus  = f"{PROCESSED_DIR}/corpus.fasta",
        dataset = f"{RAW_DIR}/nextclade_dataset",
    output:
        aln = f"{INTERIM_DIR}/alignment_nextalign.fasta",
    threads:
        workflow.cores
    log:
        f"{INTERIM_DIR}/.align_nextalign.log",
    shell:
        r"""
        # Nextclade v3 includes the Nextalign aligner. Running with the same
        # mpox/clade-i dataset that drove the QC step ensures coordinate-frame
        # compatibility with the MAFFT alignment.
        nextclade run \
            --input-dataset {input.dataset} \
            --output-fasta {output.aln} \
            --jobs {threads} \
            {input.corpus} 2>&1 | tee {log}
        """

# ============================================================================
# 5. Masking
# ============================================================================

rule mask_alignment:
    """Apply the mask BED to a multiple-sequence alignment."""
    input:
        aln = "{prefix}.fasta",
        bed = f"{INTERIM_DIR}/mask.bed",
    output:
        masked = "{prefix}.masked.fasta",
    wildcard_constraints:
        prefix = ".*alignment_(mafft|nextalign)",
    log:
        "{prefix}.masked.log",
    shell:
        r"""
        python3 scripts/apply_mask.py \
            --in-aln {input.aln} \
            --bed {input.bed} \
            --out-aln {output.masked} 2>&1 | tee {log}
        """

# ============================================================================
# 6. Convenience targets
# ============================================================================

rule alignment_complete:
    """Both primary and alternative alignments, both masked."""
    input:
        f"{INTERIM_DIR}/alignment_mafft.masked.fasta",
        f"{INTERIM_DIR}/alignment_nextalign.masked.fasta",
