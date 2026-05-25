"""
Phase 1 rules — data acquisition and quality control.

Maps to pre-registration §3.3 (corpus construction), §5.0 steps 1 (fetch) and 2 (QC).

Inputs
------
- NCBI Virus / GenBank: all public MPXV genomes via the ncbi-datasets-cli.
- Nextstrain mpox open builds: clade-i and all-clades JSON snapshots.
- Nextclade mpox clade-i dataset: reference, gene map, QC framework.

Outputs (under data/raw/, data/interim/, data/processed/)
--------------------------------------------------------
- data/raw/freeze_{DATE}/mpox_genomes.zip          — NCBI Datasets bundle
- data/raw/freeze_{DATE}/mpox_genomes.fasta        — extracted FASTA
- data/raw/freeze_{DATE}/mpox_metadata.tsv         — NCBI raw metadata
- data/raw/freeze_{DATE}/nextstrain_clade_i.json   — Nextstrain build snapshot
- data/raw/freeze_{DATE}/nextstrain_all_clades.json — cross-clade snapshot
- data/raw/freeze_{DATE}/nextclade_dataset/        — Nextclade reference + tree
- data/interim/freeze_{DATE}/nextclade.tsv         — Nextclade QC + clade calls
- data/processed/freeze_{DATE}/corpus.fasta        — post-filter analysis corpus
- data/processed/freeze_{DATE}/metadata.tsv        — merged metadata for the corpus
- data/processed/freeze_{DATE}/exclusions.tsv      — every excluded accession + reason
- data/processed/freeze_{DATE}/freeze_manifest.json — provenance, SHA-256 hashes, tool versions
"""

# ============================================================================
# 1. Fetch — NCBI Virus genome bundle via ncbi-datasets-cli
# ============================================================================

rule fetch_ncbi_genomes:
    """
    Pull all public Monkeypox virus genomes from NCBI Datasets.
    Uses the modern ncbi-datasets-cli (NCBI Datasets v2 API).
    """
    output:
        zipfile = f"{RAW_DIR}/mpox_genomes.zip",
    params:
        email = config["NCBI_EMAIL"],
    log:
        f"{RAW_DIR}/.fetch_ncbi_genomes.log",
    shell:
        r"""
        mkdir -p $(dirname {output.zipfile})
        echo "[$(date -u +%FT%TZ)] datasets download virus genome taxon 'Monkeypox virus'" | tee {log}
        datasets download virus genome taxon "Monkeypox virus" \
            --filename {output.zipfile} \
            --include genome,annotation 2>&1 | tee -a {log}
        echo "[$(date -u +%FT%TZ)] SHA256: $(shasum -a 256 {output.zipfile} | awk '{{print $1}}')" | tee -a {log}
        """

rule extract_ncbi_bundle:
    """Unpack the NCBI Datasets zip into FASTA + metadata TSV."""
    input:
        zipfile = f"{RAW_DIR}/mpox_genomes.zip",
    output:
        fasta    = f"{RAW_DIR}/mpox_genomes.fasta",
        metadata = f"{RAW_DIR}/mpox_metadata.tsv",
    log:
        f"{RAW_DIR}/.extract_ncbi_bundle.log",
    shell:
        r"""
        TMPDIR=$(mktemp -d)
        unzip -o {input.zipfile} -d $TMPDIR 2>&1 | tee {log}
        # NCBI Datasets bundles place sequences at ncbi_dataset/data/genomic.fna
        # and metadata via the included data_report.jsonl.
        find $TMPDIR -name 'genomic.fna' -exec cat {{}} + > {output.fasta}
        # Convert data_report.jsonl to a TSV with the fields we need.
        REPORT=$(find $TMPDIR -name 'data_report.jsonl' | head -1)
        if [ -z "$REPORT" ]; then
            echo "ERROR: data_report.jsonl not found in NCBI bundle." | tee -a {log}
            exit 1
        fi
        # Field names per NCBI dataformat v18 (see `dataformat tsv virus-genome --help`).
        dataformat tsv virus-genome \
            --fields accession,length,isolate-collection-date,geo-location,host-name,completeness,sourcedb,isolate-lineage,virus-tax-id,release-date \
            --inputfile $REPORT > {output.metadata}
        echo "[$(date -u +%FT%TZ)] Extracted $(wc -l < {output.metadata}) metadata rows." | tee -a {log}
        rm -rf $TMPDIR
        """

# ============================================================================
# 2. Fetch — Nextstrain build snapshots + Nextclade dataset
# ============================================================================

rule fetch_nextstrain_builds:
    """Pull the two Nextstrain mpox build snapshots used in this analysis."""
    output:
        clade_i    = f"{RAW_DIR}/nextstrain_clade_i.json",
        all_clades = f"{RAW_DIR}/nextstrain_all_clades.json",
    log:
        f"{RAW_DIR}/.fetch_nextstrain_builds.log",
    shell:
        r"""
        mkdir -p $(dirname {output.clade_i})
        echo "[$(date -u +%FT%TZ)] fetching Nextstrain mpox JSONs" | tee {log}
        curl -fL --retry 3 --retry-delay 5 \
            -o {output.clade_i} \
            https://data.nextstrain.org/mpox_clade-I.json 2>&1 | tee -a {log}
        curl -fL --retry 3 --retry-delay 5 \
            -o {output.all_clades} \
            https://data.nextstrain.org/mpox_all-clades.json 2>&1 | tee -a {log}
        echo "[$(date -u +%FT%TZ)] SHA256 clade-i: $(shasum -a 256 {output.clade_i} | awk '{{print $1}}')" | tee -a {log}
        echo "[$(date -u +%FT%TZ)] SHA256 all-clades: $(shasum -a 256 {output.all_clades} | awk '{{print $1}}')" | tee -a {log}
        """

rule fetch_nextclade_dataset:
    """Pull the Nextclade mpox/clade-i dataset (reference + clade tree + gene map)."""
    output:
        directory(f"{RAW_DIR}/nextclade_dataset"),
    params:
        dataset = config["NEXTCLADE_DATASET_PRIMARY"],
    log:
        f"{RAW_DIR}/.fetch_nextclade_dataset.log",
    shell:
        r"""
        nextclade dataset get --name {params.dataset} --output-dir {output} 2>&1 | tee {log}
        # Record the commit / version of the dataset that was pulled.
        if [ -f {output}/pathogen.json ]; then
            echo "[$(date -u +%FT%TZ)] dataset version:" | tee -a {log}
            python3 -c "import json; d=json.load(open('{output}/pathogen.json')); print(json.dumps({{k: d.get(k) for k in ('attributes','version','meta')}}, indent=2))" 2>&1 | tee -a {log}
        fi
        """

# ============================================================================
# 3. Quality control — Nextclade
# ============================================================================

rule run_nextclade:
    """
    Run Nextclade on the raw NCBI MPXV FASTA against the mpox/clade-i dataset.
    Produces a TSV with clade assignment, QC flags, frameshift calls, and
    per-genome diagnostics needed for the §3.3 inclusion / exclusion filter.
    """
    input:
        fasta   = f"{RAW_DIR}/mpox_genomes.fasta",
        dataset = f"{RAW_DIR}/nextclade_dataset",
    output:
        tsv = f"{INTERIM_DIR}/nextclade.tsv",
    threads:
        workflow.cores
    log:
        f"{INTERIM_DIR}/.run_nextclade.log",
    shell:
        r"""
        mkdir -p $(dirname {output.tsv})
        nextclade run \
            --input-dataset {input.dataset} \
            --output-tsv {output.tsv} \
            --jobs {threads} \
            {input.fasta} 2>&1 | tee {log}
        echo "[$(date -u +%FT%TZ)] Nextclade rows: $(wc -l < {output.tsv})" | tee -a {log}
        """

# ============================================================================
# 4. Apply §3.3 inclusion / exclusion filters
# ============================================================================

rule apply_filters:
    """
    Apply pre-registration §3.3 inclusion and exclusion criteria.
    Produces (a) the analysis corpus FASTA + metadata, (b) an exclusions TSV
    documenting every excluded accession with the reason.
    """
    input:
        fasta     = f"{RAW_DIR}/mpox_genomes.fasta",
        ncbi_meta = f"{RAW_DIR}/mpox_metadata.tsv",
        nextclade = f"{INTERIM_DIR}/nextclade.tsv",
    output:
        corpus_fasta = f"{PROCESSED_DIR}/corpus.fasta",
        metadata     = f"{PROCESSED_DIR}/metadata.tsv",
        exclusions   = f"{PROCESSED_DIR}/exclusions.tsv",
        summary      = f"{PROCESSED_DIR}/filter_summary.json",
    params:
        filters = config["FILTERS"],
    log:
        f"{PROCESSED_DIR}/.apply_filters.log",
    shell:
        r"""
        mkdir -p $(dirname {output.corpus_fasta})
        python3 scripts/apply_filters.py \
            --in-fasta {input.fasta} \
            --in-ncbi-meta {input.ncbi_meta} \
            --in-nextclade {input.nextclade} \
            --out-fasta {output.corpus_fasta} \
            --out-metadata {output.metadata} \
            --out-exclusions {output.exclusions} \
            --out-summary {output.summary} \
            --min-length {params.filters[min_genome_length_nt]} \
            --max-n-fraction {params.filters[max_n_fraction]} \
            --max-non-acgtn-fraction {params.filters[max_non_acgtn_fraction]} \
            --date-resolution {params.filters[require_collection_date_resolution]} \
            2>&1 | tee {log}
        """

# ============================================================================
# 5. Build the freeze manifest — SHA-256 hashes, tool versions, provenance
# ============================================================================

rule build_freeze_manifest:
    """
    Assemble the freeze manifest: provenance, SHA-256 hashes of every freeze
    file, tool versions, NCBI / Nextstrain / Nextclade dataset versions.
    This is the artefact that backs the reproducibility claims in the pre-reg.
    """
    input:
        raw_fasta     = f"{RAW_DIR}/mpox_genomes.fasta",
        raw_meta      = f"{RAW_DIR}/mpox_metadata.tsv",
        ns_clade_i    = f"{RAW_DIR}/nextstrain_clade_i.json",
        ns_all_clades = f"{RAW_DIR}/nextstrain_all_clades.json",
        nc_dataset    = f"{RAW_DIR}/nextclade_dataset",
        nc_tsv        = f"{INTERIM_DIR}/nextclade.tsv",
        corpus_fasta  = f"{PROCESSED_DIR}/corpus.fasta",
        corpus_meta   = f"{PROCESSED_DIR}/metadata.tsv",
        exclusions    = f"{PROCESSED_DIR}/exclusions.tsv",
        filter_summary = f"{PROCESSED_DIR}/filter_summary.json",
    output:
        manifest = f"{PROCESSED_DIR}/freeze_manifest.json",
        sha256   = f"{PROCESSED_DIR}/SHA256SUMS",
    params:
        freeze_date = FREEZE_DATE,
    log:
        f"{PROCESSED_DIR}/.build_freeze_manifest.log",
    shell:
        r"""
        python3 scripts/build_freeze_manifest.py \
            --freeze-date {params.freeze_date} \
            --raw-dir {RAW_DIR} \
            --interim-dir {INTERIM_DIR} \
            --processed-dir {PROCESSED_DIR} \
            --out-manifest {output.manifest} \
            --out-sha256 {output.sha256} \
            2>&1 | tee {log}
        """
