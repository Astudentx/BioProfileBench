# BioProfileBench

BioProfileBench is a lightweight Python toolkit for benchmarking microbiome abundance profiles against truth abundance tables.

Current version: `v0.0.0`

## Features

- Compare one truth abundance table against one prediction abundance table.
- Run batch benchmarks from a manifest TSV.
- Aggregate profiles across taxonomic ranks: kingdom, phylum, class, order, family, genus, and species.
- Apply prediction-side filters with stable short filter IDs.
- Export R-friendly benchmark tables and per-rank summary tables.
- Compute presence/absence, abundance-weighted, correlation, and microbiome distance metrics.

## Main Outputs

- `benchmark.tsv`: sample-level benchmark results.
- `benchmark_summary.tsv`: summary by truth file, prediction file, filter ID, and taxonomic rank.
- `benchmark_summary.<level>.tsv`: per-rank summaries.
- `filter_trace.tsv`: filter ID to filter-parameter mapping.
- `pair_filter_ranking.tsv`: ranking of truth/prediction/filter combinations.

## Quick Start

```bash
python BioProfileBench.py run \
  --truth truth.tsv \
  --pred prediction.tsv \
  --config config.yaml \
  --out results \
  --gene-kind ARGs
```

Batch mode:

```bash
python BioProfileBench.py batch \
  --manifest manifest.tsv \
  --config config.yaml \
  --threads 4 \
  --out results
```

## Manifest Format

```text
profile_id	dataset	method	gene_kind	prediction_path	truth_path
run_001	Dataset1	Prediction	ARGs	prediction.tsv	truth.tsv
```

## Minimal Config

```yaml
taxonomy:
  source: Lineage
  levels: [kingdom, phylum, class, order, family, genus, species]
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

## Dependencies

- Python 3.9+
- pandas
- numpy
- scipy

## Status

This is an initial development release. The API and output format may still change before a stable release.
