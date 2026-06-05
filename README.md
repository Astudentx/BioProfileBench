# BioProfileBench

[中文说明](README-CN.md)

BioProfileBench is a Python toolkit for benchmarking predicted microbiome abundance profiles against truth abundance profiles. It supports one-to-one and batch comparisons, taxonomic-rank aggregation, prediction-only filtering, stable filter identifiers, and R-friendly result tables.

Current version: `v1.0.0`

## Key Principles

- The truth profile is treated as the fixed reference and is **never abundance-filtered**.
- Filtering is applied only to the prediction profile.
- Truth and prediction taxa are aligned using their union before metrics are calculated.
- Each sample and each taxonomic rank are benchmarked independently.
- Multiple filters are calculated from the same unfiltered prediction matrix and combined using logical **AND**, so filtering is order-independent.

## Features

- One-to-one benchmarking with the `run` command.
- Many-to-many batch benchmarking with a manifest TSV.
- Taxonomic aggregation across kingdom, phylum, class, order, family, genus, and species.
- Prediction-side cell-level, taxon-level, prevalence, top-N, and cumulative-abundance filters.
- Grid-based and data-driven dynamic threshold generation.
- Stable short filter IDs such as `F000` and `FBCB70CDD`.
- Presence/absence, abundance-weighted, correlation, and microbiome distance metrics.
- Overall and per-taxonomic-rank summary tables.

## Dependencies

- Python 3.9+
- pandas
- numpy
- scipy

## Input Abundance Table

Input files must be tab-separated. The first column may be `ID`, `TaxID`, `SpeciesID`, `FeatureID`, or `GeneID`. Sample columns contain abundance values. `Taxonomy` and/or `Lineage` are used for taxonomic aggregation.

```text
ID	Sample1	Sample2	Taxonomy	Lineage
gene1	10	0	Escherichia_coli	k__Bacteria;p__Proteobacteria;g__Escherichia;s__Escherichia_coli
gene2	0	5	Klebsiella_pneumoniae	k__Bacteria;p__Proteobacteria;g__Klebsiella;s__Klebsiella_pneumoniae
```

The truth and prediction tables must contain the same sample columns.

## Quick Start

### One-to-One Benchmark

```bash
python BioProfileBench.py run \
  --truth truth.tsv \
  --pred prediction.tsv \
  --config config.yaml \
  --out results \
  --kind Bacteria
```

`--kind` is a general profile label. It may be `ARGs`, `VFs`, `Bac`, `Bacteria`, or another user-defined category. The legacy `--gene-kind` option remains available as a compatibility alias.

### Batch Benchmark

```bash
python BioProfileBench.py batch \
  --manifest manifest.tsv \
  --config config.yaml \
  --threads 4 \
  --out results
```

Recommended manifest format:

```text
profile_id	dataset	method	kind	prediction_path	truth_path
run_001	Dataset1	Prediction	ARGs	prediction_ARG.tsv	truth_ARG.tsv
run_002	Dataset1	Prediction	Bacteria	prediction_Bac.tsv	truth_Bac.tsv
```

The legacy `gene_kind` manifest column is still accepted, but `kind` is preferred.

## Minimal Configuration

```yaml
taxonomy:
  source: Lineage
  levels: [kingdom, phylum, class, order, family, genus, species]
  exclude_unclassified: false

abundance:
  normalize: relative
  presence_threshold: 0

benchmark:
  digits: 6

filtering:
  modes: [raw, common_filter]
  common_filters:
    pred_min_relative_abundance: [1e-06]
```

## Filtering Workflow

For each requested taxonomic rank, BioProfileBench performs the following steps:

1. Aggregate truth and prediction features to the selected taxonomic rank.
2. Apply `abundance.normalize` independently to truth and prediction.
3. Build prediction filter configurations for that rank.
4. Calculate every filter mask independently from the same prediction matrix.
5. Combine all enabled masks using logical AND.
6. Set failed cell-level entries to zero.
7. Remove taxa that fail row-level filters.
8. Remove taxa whose final abundance is zero across every sample.
9. Align filtered prediction and unfiltered truth using the taxa union.
10. Calculate benchmark metrics for each sample.

Because masks are calculated independently before combination, changing filter order does not change the result. However, combining several strict filters may still remove most or all predicted taxa.

### Filtering Modes

```yaml
filtering:
  modes: [raw, common_filter, optimized_filter]
```

- `raw`: no prediction filtering; filter ID is `F000`.
- `common_filter`: apply the configured filter grid for direct comparison.
- `optimized_filter`: apply the configured filter grid and mark results for downstream best-filter selection.

Any non-`raw` mode uses the resolved filter combinations. The mode label itself does not change the mathematical filter behavior and is not included in the stable filter ID.

## Supported Prediction Filters

### Cell-Level Filters

Cell-level filters evaluate each `taxon × sample` abundance independently. Failed cells are replaced with zero.

#### `pred_min_abundance`

```yaml
pred_min_abundance: [0, 10, 100]
```

Keeps a cell only when:

```text
prediction abundance > threshold
```

This filter uses the matrix after `abundance.normalize`. Therefore, when `normalize: relative`, this threshold is also interpreted on a relative-abundance scale.

#### `pred_min_relative_abundance`

```yaml
pred_min_relative_abundance: [1e-06, 1e-05, 1e-04]
```

For each sample, prediction abundances are converted to relative abundance. A cell is kept only when:

```text
relative abundance > threshold
```

This filter always evaluates relative abundance, regardless of whether the main benchmark uses raw or relative abundance.

### Taxon-Level Filters

Taxon-level filters evaluate one taxon across all samples. A failed taxon row is removed entirely.

#### `pred_min_total_abundance`

```yaml
pred_min_total_abundance: [10, 100]
```

Keeps a taxon when its abundance summed across all samples is at least the threshold.

#### `pred_min_mean_abundance`

```yaml
pred_min_mean_abundance: [0.001, 0.01]
```

Keeps a taxon when its mean abundance across samples is at least the threshold.

#### `pred_min_max_abundance`

```yaml
pred_min_max_abundance: [0.01, 0.1]
```

Keeps a taxon when its maximum abundance among all samples is at least the threshold.

#### `pred_min_prevalence`

```yaml
pred_min_prevalence: [0.1, 0.2, 3]
```

Prevalence is the number of samples in which prediction abundance is greater than `abundance.presence_threshold`.

- `value <= 0`: no effective prevalence restriction.
- `0 < value < 1`: minimum sample proportion. For example, `0.2` means present in at least 20% of samples.
- `value >= 1`: minimum sample count. For example, `3` means present in at least three samples.

Fractional sample requirements are rounded upward.

### Top-Abundance Filters

Top-abundance filters support two modes:

```yaml
top_mode: global       # default
top_mode: per_sample
```

- `global`: rank taxa by abundance summed across all samples and remove entire taxon rows.
- `per_sample`: rank taxa separately within each sample and set failed cells to zero.

#### `pred_top_k`

```yaml
pred_top_k: [50, 100]
top_mode: [global]
```

Keeps the top K most abundant predicted taxa.

#### `pred_top_percent`

```yaml
pred_top_percent: [0.1, 0.2]
top_mode: [global]
```

- Values between 0 and 1 are interpreted as fractions of the taxa count.
- Values greater than 1 are interpreted as an explicit number of taxa.

For example, `0.1` keeps approximately the top 10% of taxa, while `100` keeps approximately the top 100 taxa.

#### `pred_cumulative_abundance`

```yaml
pred_cumulative_abundance: [0.90, 0.95, 0.99]
top_mode: [global]
```

Sorts taxa by descending abundance and keeps taxa until the selected cumulative abundance fraction is reached. The first taxon that crosses the threshold is also retained, preventing an empty selection.

## Combining Filters

Filter values listed under `common_filters` are expanded as a Cartesian product.

```yaml
filtering:
  modes: [raw, common_filter]
  common_filters:
    pred_min_relative_abundance: [1e-06, 1e-05]
    pred_min_prevalence: [0.1, 0.2]
```

This produces four effective non-raw combinations:

```text
1e-06 + prevalence 0.1
1e-06 + prevalence 0.2
1e-05 + prevalence 0.1
1e-05 + prevalence 0.2
```

Within one combination, all filters must pass. Cell-level masks and taxon-level masks are independently calculated from the original prediction matrix, then combined using AND. This avoids filter-order dependence.

## Dynamic Thresholds

Dynamic filters generate threshold values from the prediction matrix separately for each taxonomic rank.

### Quantile Thresholds

```yaml
filtering:
  modes: [raw, optimized_filter]
  dynamic_filters:
    pred_min_relative_abundance:
      strategy: quantile
      quantiles: [0.5, 0.75, 0.9, 0.95]
```

Uses selected quantiles of non-zero prediction abundances.

### Log-Spaced Thresholds

```yaml
filtering:
  dynamic_filters:
    pred_min_abundance:
      strategy: logspace
      min: 1e-06
      max: 1e-02
      count: 5
```

Generates logarithmically spaced thresholds.

### Target Taxa Count

```yaml
filtering:
  dynamic_filters:
    pred_min_total_abundance:
      strategy: target_taxa_count
      targets: [50, 100, 200]
```

Derives total-abundance thresholds intended to retain approximately the selected taxa counts.

### Target Retained Abundance

```yaml
filtering:
  dynamic_filters:
    pred_min_total_abundance:
      strategy: target_retained_abundance
      targets: [0.90, 0.95, 0.99]
```

Derives thresholds from cumulative total abundance.

### Prevalence Sweep

```yaml
filtering:
  dynamic_filters:
    pred_min_prevalence:
      strategy: prevalence_sweep
```

Generates every integer prevalence threshold from zero to the total sample count.

Dynamic thresholds are data-dependent. Their stable filter IDs are based on the final resolved numeric thresholds, not only the dynamic strategy name.

## Stable Filter IDs

Every effective filter combination receives a deterministic short identifier.

- `F000`: no effective filtering.
- Other IDs, such as `FBCB70CDD`, represent a specific resolved filter parameter combination.
- Equivalent values such as `0`, `0.0`, or omitted no-op thresholds produce the same ID.
- IDs ignore file paths, batch name, folder, taxonomic rank, and analysis mode.
- The same effective filter parameters produce the same ID across different runs.

Actual filter parameters are recorded in `filter_trace.tsv`, so compact IDs can be used in benchmark tables without losing reproducibility.

## Filter Diagnostics

`filter_diagnostics.tsv` records the effect of each filter:

- `filter_name`: filter being evaluated.
- `threshold`: configured threshold.
- `taxa_before`: number of taxa before filtering.
- `taxa_pass`: taxa passing that individual filter.
- `taxa_fail`: taxa failing that individual filter.
- `unique_fail_taxa`: taxa removed only by that filter while passing the others.
- `pred_retained_abundance`: fraction of prediction abundance retained.
- `FINAL`: final result after all masks are combined.

Use this table to identify redundant thresholds, overly strict filter combinations, and cases where nearly all prediction taxa are removed.

## Benchmark Metrics

BioProfileBench reports:

- Presence/absence: TP, FP, FN, TN, Precision, Recall, F1, Jaccard, Accuracy, Specificity, MCC, AUPRC.
- Abundance overlap: AbundancePrecision, AbundanceRecall, AbundanceF1, WeightedJaccard.
- Correlation/similarity: Spearman, Pearson, Cosine.
- Distance/error: L1, L2, Bray-Curtis (`BC`), rooted Jensen-Shannon distance (`rJSD`), MAE, RMSE.

## Main Outputs

- `benchmark.tsv`: one row per sample × taxonomic rank × file pair × filter ID.
- `benchmark_long.tsv`: detailed internal benchmark result table.
- `benchmark_summary.tsv`: summary by `Kind`, truth file, prediction file, filter ID, and taxonomic rank.
- `benchmark_summary.<level>.tsv`: separate summary for each taxonomic rank.
- `filter_trace.tsv`: stable filter ID and resolved filter parameters.
- `filter_diagnostics.tsv`: effect of individual filters and final combination.
- `pair_filter_ranking.tsv`: ranking across all taxonomic ranks.
- `pair_filter_level_ranking.tsv`: ranking separately for each taxonomic rank.
- `failed_or_skipped_runs.tsv`: empty or skipped benchmark records.

## Notes and Limitations

- Prediction filtering occurs after taxonomic aggregation, so the same threshold may behave differently at species and genus levels.
- Dynamic thresholds are recalculated for each prediction profile and taxonomic rank.
- Strict AND combinations may produce empty prediction profiles; these are reported rather than silently discarded.
- `Taxonomy` without prefixed rank information is treated as a species label. Use `Lineage` for multi-rank benchmarking.
