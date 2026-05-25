# Data dictionary

This dictionary documents every dataset produced by the analysis pipeline and consumed by the deposited result artefacts in `data/processed/`. For raw NCBI / Nextstrain source data, see `data/raw/README.md`.

## `data/processed/freeze_manifest.json`

Provenance manifest for the 2026-05-22 data freeze.

| Field | Type | Description |
|---|---|---|
| `freeze_date` | string (YYYY-MM-DD) | The freeze date — `2026-05-22` |
| `inputs` | object | Per-input artefact record: source URL, SHA-256 hash, byte count, retrieval timestamp |
| `tools` | object | Pinned tool versions: Nextclade, MAFFT, IQ-TREE, TreeTime, HyPhy, 3SEQ, PhiPack, NCBI Datasets CLI |
| `python_env` | object | Pinned Python package versions (numpy, pandas, scipy, biopython, dendropy, ete3, statsmodels, snakemake) |

## `data/processed/h1prime_result.json`

Result of the H1' test (Amendment 02 branch-quantity protocol).

| Field | Type | Description |
|---|---|---|
| `n_genomes_in_test` | int | Number of clade-Ib genomes contributing to the pooled fraction (198) |
| `restrict_clade_prefix` | string | Clade prefix used to restrict the test (`Ib`) |
| `total_snvs` | int | Pooled total branch SNVs across all 198 tips |
| `apobec3_snvs` | int | Pooled APOBEC3-context branch SNVs across all 198 tips |
| `non_apobec3_snvs` | int | Pooled non-APOBEC3 branch SNVs |
| `observed_apobec3_fraction` | float | Pooled fraction — 0.8358 |
| `wilson_95_lower` / `wilson_95_upper` | float | Wilson 95 % confidence interval bounds |
| `null_p` | float | Pre-registered null *p*₀ — 0.70 |
| `alternative` | string | One-sided test direction — `greater` |
| `exact_binomial_p_value` | float | P-value against H₀: p = 0.70 |
| `h1_supported` | bool | Decision: TRUE if Wilson lower bound exceeds null *p*₀ |

## `data/processed/saturation_dating_results.json`

Result of the §5.3 saturation-aware Poisson dating fit on the clade-Ib cluster.

| Field | Type | Description |
|---|---|---|
| `L` | int | Number of eligible APOBEC3 target sites on the masked reference (23,140) |
| `freeze_date` | string | Data freeze date |
| `clusters` | array | Per-cluster fit record |
| `clusters[i].cluster_id` | string | Cluster identifier (`Ib-global` for the single 198-genome Ib cluster) |
| `clusters[i].n_genomes` | int | Genome count in the cluster |
| `clusters[i].saturation.status` | string | Fit status (`fit_ok` / `not_identifiable` / `fit_failed`) |
| `clusters[i].saturation.t0_hat` | string (ISO datetime) | Saturation-aware Poisson MRCA point estimate (registered primary) |
| `clusters[i].saturation.lambda_hat_per_site_per_day` | float | Fitted per-site rate λ |
| `clusters[i].saturation.ci95_lo` / `ci95_hi` | string | Parametric bootstrap 95 % CI bounds for *t*₀ |
| `clusters[i].saturation.profile_range_LL_units` | float | Profile log-likelihood range across ±1 year (identifiability check) |
| `clusters[i].uncorrected.*` | mixed | Same fields for the uncorrected linearisation comparator |
| `clusters[i].delta_days` | float | Days between saturation-aware and uncorrected estimates (H2 metric) |
| `clusters[i].h2_supported` | bool | Decision: TRUE if Δ > 60 days |

## `data/processed/lbfgsb_trajectory.csv`

Per-start convergence trajectory of the 25 random L-BFGS-B starts used by the registered saturation-aware Poisson fit (§6.4 Amendment 04).

| Column | Type | Description |
|---|---|---|
| `start_idx` | int | Start index (0–24) |
| `t0_init_ordinal` | int | Initial *t*₀ as proleptic-Gregorian ordinal day |
| `lam_init` | float | Initial λ (log-uniform sample on [10⁻¹⁰, 10⁻³]) |
| `t0_init_year` | float | Initial *t*₀ expressed as fractional year |
| `t0_converged_ordinal` | int | Converged *t*₀ as ordinal day |
| `t0_converged_date` | string (YYYY-MM-DD) | Converged *t*₀ as date |
| `lam_converged` | float | Converged λ |
| `t0_converged_year` | float | Converged *t*₀ as fractional year |
| `final_nll` | float | Final negative log-likelihood |
| `converged_status` | bool | SciPy L-BFGS-B `success` flag |
| `basin` | string | Classification: `~2019 local minimum`, `~2023 global minimum`, or `other` (based on converged *t*₀ year only) |

## `outputs/tables/Table_S1_candidate_annotation.tsv`

Per-candidate annotation of the 202 sequences flagged by 3SEQ as recombinant children in the §5.4 scan.

| Column | Type | Description |
|---|---|---|
| `accession` | string | NCBI GenBank accession |
| `is_positive_control` | bool | TRUE for OZ375330.1 only |
| `phipack_call` | bool | PhiPack consensus call (≥2 of 3 internal statistics significant at default P < 0.05) |
| `phi_perm_p` | float | Φ permutation P-value (1000 permutations) |
| `nss_p` | float | NSS P-value |
| `maxchi_p` | float | Max Chi² P-value |
| `n_stats_significant_at_p05` | int | Count of internal statistics significant at P < 0.05 |
| `clade` | string | Nextclade clade assignment (Ib for all 202 in the present scan) |
| `Length` | int | Genome length in bp |
| `Isolate Collection date` | string | NCBI collection date (resolution varies) |
| `Geographic Location` | string | NCBI geographic location field |
| `is_unverified` | bool | NCBI UNVERIFIED prefix in the title |
| `is_partial_genome` | bool | NCBI partial-genome flag in the title |
| `ncbi_description` | string | Full NCBI sequence description |

## `outputs/tables/Table_S2_sensitivity_summary.tsv`

Per-replicate summary across the 11 §5.5 sensitivity replicates (5 × 75 %-subsample + 5 × 50 %-subsample + 1 temporal-holdout).

| Column | Type | Description |
|---|---|---|
| `analysis` | string | Analysis type (`subsample_50`, `subsample_75`, `temporal_holdout`) |
| `replicate` | string | Replicate identifier |
| `n_genomes` | int | Number of genomes in the subsample |
| `sat_t0` | string | Saturation-aware Poisson MRCA point estimate |
| `sat_ci_lo` / `sat_ci_hi` | string | Bootstrap 95 % CI bounds |
| `unc_t0` | string | Uncorrected linearisation MRCA |
| `delta_days` | float | Days between saturation-aware and uncorrected estimates |
| `h2_supported` | bool | TRUE if Δ > 60 days on this replicate |

## `outputs/tables/Table_S3_detector_sweep.csv`

Detector-parameter sweep summary (§5.5 item 4).

| Column | Type | Description |
|---|---|---|
| `detector` | string | `3SEQ` / `PhiPack` / `GARD` |
| `setting` | string | Parameter setting (e.g., `alpha×0.5`, `internal≥2/3`, `ΔAIC×1.25`) |
| `n_candidates` | int | Count of consensus candidates at this setting |
| `positive_control_in_set` | bool | TRUE if the locked positive control OZ375330.1 passes at this setting |

## `outputs/tables/Table_S4_lbfgsb_trajectory.csv`

Mirror of `data/processed/lbfgsb_trajectory.csv` — see column definitions above.
