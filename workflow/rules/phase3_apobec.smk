"""
Phase 3 rules — APOBEC3 edit counting and the H1 test.

Maps to pre-registration §5.0 step 4 and §5.2.

Inputs (from Phase 2)
---------------------
- data/interim/freeze_{DATE}/alignment_mafft.masked.fasta  — primary masked MSA
- data/interim/freeze_{DATE}/reference.fasta                — NC_063383

Outputs
-------
- data/processed/freeze_{DATE}/apobec3_counts.tsv   — per-genome APOBEC3 vs non-APOBEC3 SNV counts
- data/processed/freeze_{DATE}/apobec3_sites.bed    — per-site occupancy on the reference
- data/processed/freeze_{DATE}/h1_result.json       — Wilson interval + exact-binomial P-value
- outputs/figures/snv_by_date.png                    — diagnostic plot (SNV count vs collection date)
- outputs/tables/baseline.csv                        — headline numbers for the manuscript
"""

# ============================================================================
# 1. APOBEC3 SNV counting per genome
# ============================================================================

rule count_apobec3:
    """
    Count APOBEC3-context SNVs (TC→TT and the GA→AA equivalent) per genome,
    relative to NC_063383, using the masked primary alignment.
    """
    input:
        aln       = f"{INTERIM_DIR}/alignment_mafft.masked.fasta",
        reference = f"{INTERIM_DIR}/reference.fasta",
        metadata  = f"{PROCESSED_DIR}/metadata.tsv",
    output:
        counts = f"{PROCESSED_DIR}/apobec3_counts.tsv",
        sites  = f"{PROCESSED_DIR}/apobec3_sites.bed",
        L_json = f"{PROCESSED_DIR}/L_eligible_sites.json",
    log:
        f"{PROCESSED_DIR}/.count_apobec3.log",
    shell:
        r"""
        python3 scripts/count_apobec3.py \
            --in-aln {input.aln} \
            --in-reference {input.reference} \
            --in-metadata {input.metadata} \
            --out-counts {output.counts} \
            --out-sites {output.sites} \
            --out-L {output.L_json} 2>&1 | tee {log}
        """

# ============================================================================
# 2. H1 test — clade-wide APOBEC3 fraction vs p0 = 0.70
# ============================================================================

rule test_h1:
    """
    H1: one-sided exact binomial against p0 = 0.70; supported if Wilson lower
    95% bound > 0.70. Emergency-halt trigger if observed < 0.50.
    """
    input:
        counts   = f"{PROCESSED_DIR}/apobec3_counts.tsv",
        metadata = f"{PROCESSED_DIR}/metadata.tsv",
    output:
        result = f"{PROCESSED_DIR}/h1_result.json",
    log:
        f"{PROCESSED_DIR}/.test_h1.log",
    shell:
        r"""
        python3 scripts/test_h1.py \
            --in-counts {input.counts} \
            --in-metadata {input.metadata} \
            --restrict-clade Ib \
            --null-p 0.70 \
            --alpha 0.05 \
            --halt-threshold 0.50 \
            --out-result {output.result} 2>&1 | tee {log}
        """

# ============================================================================
# 3. Baseline diagnostic plots and tables (replicates O'Toole-style SNV-vs-date)
# ============================================================================

rule baseline_plot:
    """SNV count vs collection date — diagnostic sanity check."""
    input:
        counts   = f"{PROCESSED_DIR}/apobec3_counts.tsv",
        metadata = f"{PROCESSED_DIR}/metadata.tsv",
    output:
        png      = "outputs/figures/snv_by_date.png",
        baseline = "outputs/tables/baseline.csv",
    log:
        f"{PROCESSED_DIR}/.baseline_plot.log",
    shell:
        r"""
        mkdir -p outputs/figures outputs/tables
        python3 scripts/baseline_plot.py \
            --in-counts {input.counts} \
            --in-metadata {input.metadata} \
            --out-png {output.png} \
            --out-baseline {output.baseline} 2>&1 | tee {log}
        """

# ============================================================================
# 4. Convenience targets
# ============================================================================

rule apobec_complete:
    """All Phase-3 outputs produced."""
    input:
        f"{PROCESSED_DIR}/apobec3_counts.tsv",
        f"{PROCESSED_DIR}/h1_result.json",
        "outputs/figures/snv_by_date.png",
        "outputs/tables/baseline.csv",
