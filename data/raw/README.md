# Raw data acquisition

The raw NCBI sequence corpus and Nextstrain build snapshots used by this analysis are not committed to this repository because (a) the corpus FASTA is large (~360 MB for 893 sequences), (b) the data is freely retrievable from public sources at any time, and (c) the freeze manifest (`data/processed/freeze_manifest.json`) records the SHA-256 hashes of every input artefact so an independent retrieval at the same freeze date can be verified bit-for-bit identical.

This README documents how to retrieve the inputs.

## 1. Mpox genome corpus from NCBI Virus

The 2026-05-22 freeze used the NCBI Datasets command-line interface (v18.28.0) with the taxonomy term `Monkeypox virus`. The first-fetch command is implemented by the Snakemake `freeze_complete` target; the equivalent manual command is:

```bash
datasets download virus genome taxon "Monkeypox virus" \
  --released-after 1970-01-01 \
  --released-before 2026-05-22 \
  --include genome,annotation
unzip ncbi_dataset.zip -d data/raw/ncbi_dataset/
```

For the locked positive-control accession (Pullan et al. 2025), retrieve separately:

```bash
datasets download virus genome accession OZ375330.1 --include genome
```

## 2. Nextstrain mpox build snapshots

The clade-i and all-clades JSON build snapshots were retrieved from the Nextstrain server on the same day:

```bash
curl -L -o data/raw/nextstrain_mpox_clade-i.json \
  https://data.nextstrain.org/mpox_clade-i.json
curl -L -o data/raw/nextstrain_mpox_all-clades.json \
  https://data.nextstrain.org/mpox_all-clades.json
```

## 3. Nextclade mpox clade-i reference dataset

The Nextclade clade-i dataset (CHANGELOG date 2026-04-14) was fetched via the Nextclade CLI:

```bash
nextclade dataset get \
  --name nextstrain/mpox/clade-i \
  --output-dir data/raw/nextclade_dataset_clade-i/
```

## 4. Mask BED file

The masking BED used in the production pipeline (§6.4 Amendment 01) was vendored from the `nextstrain/mpox` repository at the freeze fetch time `2026-05-22T12:14:28Z`. The specific commit identifier is recorded in `data/processed/freeze_manifest.json`. To retrieve the same BED:

```bash
git clone https://github.com/nextstrain/mpox.git /tmp/nextstrain_mpox
cd /tmp/nextstrain_mpox && git checkout <commit-id-from-freeze-manifest>
cp /tmp/nextstrain_mpox/clade-i/defaults/mask.bed data/raw/mask.bed
```

## Verification

After retrieval, the SHA-256 hashes of the raw artefacts should match those in `data/processed/freeze_manifest.json`. The Snakemake `freeze_complete` target verifies this automatically on first run and refuses to proceed if any hash does not match the expected manifest entry.
