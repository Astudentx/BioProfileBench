"""Profile pair benchmarking."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from .filters import apply_prediction_filters, build_filter_config_grid, filter_config_table
from .metrics import CORE_COLUMNS, EXTRA_COLUMNS, compute_metrics
from .profile import Profile
from .taxonomy import aggregate_by_rank, choose_taxonomy_source, normalize_levels


METADATA_COLUMNS = [
    "profile_id",
    "dataset",
    "method",
    "gene_kind",
    "parameter_tag",
    "analysis_mode",
    "filter_config_id",
    "Status",
    "Reason",
    "threshold_selection_mode",
    "taxa_before",
    "taxa_after",
    "pred_retained_abundance",
]

BENCHMARK_COLUMNS = METADATA_COLUMNS + CORE_COLUMNS + EXTRA_COLUMNS


def order_benchmark_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Keep core columns first while preserving extra manifest/BK metadata."""
    ordered = [col for col in BENCHMARK_COLUMNS if col in df.columns]
    extras = [col for col in df.columns if col not in ordered]
    return df[ordered + extras]


def config_get(config: dict[str, Any], section: str, key: str, default: Any = None) -> Any:
    return config.get(section, {}).get(key, default)


def validate_samples(truth: Profile, prediction: Profile) -> list[str]:
    truth_samples = truth.sample_cols
    pred_samples = prediction.sample_cols
    if set(truth_samples) != set(pred_samples):
        missing_pred = sorted(set(truth_samples) - set(pred_samples))
        missing_truth = sorted(set(pred_samples) - set(truth_samples))
        details = []
        if missing_pred:
            details.append(f"missing in prediction: {missing_pred}")
        if missing_truth:
            details.append(f"missing in truth: {missing_truth}")
        raise ValueError("Sample columns are inconsistent: " + "; ".join(details))
    return truth_samples


def format_value(value: Any) -> str:
    if value is None:
        return "null"
    if isinstance(value, float):
        return f"{value:g}"
    return str(value)


def make_info(
    normalize: str,
    taxonomy_source: str,
    presence_threshold: float,
    filter_config: dict[str, Any],
    diagnostics: dict[str, Any],
) -> str:
    parts = [
        f"normalize={normalize}",
        f"taxonomy_source={taxonomy_source}",
        f"presence_threshold={presence_threshold:g}",
        f"analysis_mode={filter_config['analysis_mode']}",
        f"filter_config_id={filter_config['filter_config_id']}",
        f"taxa_before={diagnostics.get('taxa_before', 'NA')}",
        f"taxa_after={diagnostics.get('taxa_after', 'NA')}",
        f"pred_retained_abundance={format_value(diagnostics.get('pred_retained_abundance'))}",
    ]
    for key, value in sorted(filter_config.get("filters", {}).items()):
        parts.append(f"{key}={format_value(value)}")
    return "; ".join(parts)


def benchmark_sample(
    truth_matrix: pd.DataFrame,
    pred_matrix: pd.DataFrame,
    sample: str,
    level: str,
    presence_threshold: float,
    info: str,
) -> dict[str, Any]:
    union = truth_matrix.index.union(pred_matrix.index).sort_values()
    true_vec = truth_matrix.reindex(union, fill_value=0.0)[sample]
    pred_vec = pred_matrix.reindex(union, fill_value=0.0)[sample]
    row = {
        "Sample": sample,
        "Taxonomic Level": level,
        "Info": info,
    }
    row.update(compute_metrics(true_vec, pred_vec, presence_threshold))
    return row


def benchmark_profile_pair(
    truth: Profile,
    prediction: Profile,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Benchmark one truth/prediction profile pair."""
    samples = validate_samples(truth, prediction)
    taxonomy_config = config.get("taxonomy", {})
    abundance_config = config.get("abundance", {})

    levels = normalize_levels(taxonomy_config.get("levels"))
    taxonomy_source = choose_taxonomy_source(truth, prediction, taxonomy_config.get("source", "auto"))
    normalize = abundance_config.get("normalize", "none")
    presence_threshold = float(abundance_config.get("presence_threshold", 0.0))
    exclude_unclassified = bool(taxonomy_config.get("exclude_unclassified", False))

    benchmark_rows: list[dict[str, Any]] = []
    trace_tables: list[pd.DataFrame] = []
    config_tables: list[pd.DataFrame] = []
    skipped_rows: list[dict[str, Any]] = []
    decomposition_rows: list[dict[str, Any]] = []
    write_decomposition = bool(config.get("diagnostics", {}).get("fp_fn_decomposition", False))

    for level in levels:
        truth_rank = aggregate_by_rank(truth, level, taxonomy_source, normalize, exclude_unclassified)
        pred_rank = aggregate_by_rank(prediction, level, taxonomy_source, normalize, exclude_unclassified)
        filter_configs = build_filter_config_grid(config, pred_rank)
        resolved = filter_config_table(filter_configs)
        resolved.insert(0, "profile_id", prediction.profile_id)
        resolved.insert(1, "Taxonomic Level", level)
        config_tables.append(resolved)

        for filter_config in filter_configs:
            pred_filtered, trace, diagnostics = apply_prediction_filters(
                pred_rank,
                filter_config,
                presence_threshold=presence_threshold,
            )
            trace.insert(0, "profile_id", prediction.profile_id)
            trace.insert(1, "dataset", prediction.dataset)
            trace.insert(2, "method", prediction.method)
            trace.insert(3, "gene_kind", prediction.gene_kind)
            trace.insert(4, "parameter_tag", prediction.parameter_tag)
            trace.insert(5, "analysis_mode", filter_config["analysis_mode"])
            trace.insert(6, "filter_config_id", filter_config["filter_config_id"])
            trace.insert(7, "Taxonomic Level", level)
            trace_tables.append(trace)

            info = make_info(normalize, taxonomy_source, presence_threshold, filter_config, diagnostics)
            threshold_selection_mode = (
                "per_method_optimized"
                if filter_config["analysis_mode"] == "optimized_filter"
                else filter_config["analysis_mode"]
            )

            union = truth_rank.index.union(pred_filtered.index)
            if len(union) == 0:
                skipped_rows.append(
                    {
                        **prediction.metadata,
                        "analysis_mode": filter_config["analysis_mode"],
                        "filter_config_id": filter_config["filter_config_id"],
                        "Taxonomic Level": level,
                        "Status": "SKIPPED",
                        "Reason": "no_taxa_after_filter_and_no_truth_taxa",
                    }
                )
                continue

            for sample in samples:
                true_total = float(truth_rank[sample].sum()) if sample in truth_rank else 0.0
                pred_total = float(pred_filtered[sample].sum()) if sample in pred_filtered else 0.0
                if pred_total <= 0 and true_total > 0:
                    status = "EMPTY_PREDICTION"
                    reason = "prediction_total_zero_after_filter"
                elif pred_total <= 0 and true_total <= 0:
                    status = "EMPTY_BOTH"
                    reason = "truth_and_prediction_total_zero"
                else:
                    status = "OK"
                    reason = ""

                metrics = benchmark_sample(
                    truth_rank,
                    pred_filtered,
                    sample,
                    level,
                    presence_threshold,
                    info,
                )
                row = {
                    **prediction.metadata,
                    "analysis_mode": filter_config["analysis_mode"],
                    "filter_config_id": filter_config["filter_config_id"],
                    "Status": status,
                    "Reason": reason,
                    "threshold_selection_mode": threshold_selection_mode,
                    "taxa_before": diagnostics["taxa_before"],
                    "taxa_after": diagnostics["taxa_after"],
                    "pred_retained_abundance": diagnostics["pred_retained_abundance"],
                }
                row.update(metrics)
                benchmark_rows.append(row)

                if write_decomposition:
                    union = truth_rank.index.union(pred_filtered.index).sort_values()
                    true_vec = truth_rank.reindex(union, fill_value=0.0)[sample]
                    pred_vec = pred_filtered.reindex(union, fill_value=0.0)[sample]
                    true_present = true_vec > presence_threshold
                    pred_present = pred_vec > presence_threshold
                    for taxon in union:
                        if true_present.loc[taxon] and pred_present.loc[taxon]:
                            error_type = "TP"
                        elif not true_present.loc[taxon] and pred_present.loc[taxon]:
                            error_type = "FP"
                        elif true_present.loc[taxon] and not pred_present.loc[taxon]:
                            error_type = "FN"
                        else:
                            error_type = "TN"
                        decomposition_rows.append(
                            {
                                **prediction.metadata,
                                "analysis_mode": filter_config["analysis_mode"],
                                "filter_config_id": filter_config["filter_config_id"],
                                "Taxonomic Level": level,
                                "Sample": sample,
                                "taxon": taxon,
                                "error_type": error_type,
                                "true_abundance": true_vec.loc[taxon],
                                "pred_abundance": pred_vec.loc[taxon],
                            }
                        )

    benchmark_df = order_benchmark_columns(pd.DataFrame(benchmark_rows)) if benchmark_rows else pd.DataFrame()
    trace_df = pd.concat(trace_tables, ignore_index=True) if trace_tables else pd.DataFrame()
    config_df = pd.concat(config_tables, ignore_index=True).drop_duplicates() if config_tables else pd.DataFrame()
    skipped_df = pd.DataFrame(skipped_rows)
    decomposition_df = pd.DataFrame(decomposition_rows)
    return benchmark_df, trace_df, config_df, skipped_df, decomposition_df


def summarize_benchmark(benchmark: pd.DataFrame) -> pd.DataFrame:
    """Summarize benchmark rows across samples."""
    if benchmark.empty:
        return pd.DataFrame()
    group_cols = [
        "dataset",
        "method",
        "gene_kind",
        "parameter_tag",
        "analysis_mode",
        "filter_config_id",
        "Taxonomic Level",
    ]
    rows = []
    for keys, group in benchmark.groupby(group_cols, dropna=False):
        key_data = dict(zip(group_cols, keys))
        ok = group[group["Status"].isin(["OK", "EMPTY_PREDICTION", "EMPTY_BOTH"])]
        row = {
            **key_data,
            "mean_F1": ok["F1"].mean(),
            "median_F1": ok["F1"].median(),
            "mean_AbundanceF1": ok["AbundanceF1"].mean(),
            "mean_BC": ok["BC"].mean(),
            "mean_rJSD": ok["rJSD"].mean(),
            "mean_Spearman": ok["Spearman"].mean(),
            "sample_count": int(len(ok)),
            "taxa_before_mean": ok["taxa_before"].mean(),
            "taxa_after_mean": ok["taxa_after"].mean(),
            "pred_retained_abundance_mean": ok["pred_retained_abundance"].mean(),
            "skipped_count": int((group["Status"] == "SKIPPED").sum()),
            "empty_prediction_count": int((group["Status"] == "EMPTY_PREDICTION").sum()),
        }
        rows.append(row)
    return pd.DataFrame(rows)


def parse_objective(objective: str, default_direction: str = "maximize") -> tuple[str | None, str, str]:
    """Parse objective strings such as genus:F1 or family:BC:min."""
    parts = objective.split(":")
    if len(parts) == 1:
        return None, parts[0], default_direction
    if len(parts) == 2:
        return parts[0], parts[1], default_direction
    direction = "minimize" if parts[2].lower() in {"min", "minimize"} else "maximize"
    return parts[0], parts[1], direction


def select_best_filters(summary: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    """Select best optimized filters by configured objective."""
    if summary.empty:
        return pd.DataFrame()
    opt = config.get("optimization", {})
    level, metric, direction = parse_objective(
        opt.get("objective", "genus:F1"),
        opt.get("direction", "maximize"),
    )
    candidates = summary[summary["analysis_mode"] == "optimized_filter"].copy()
    if level:
        candidates = candidates[candidates["Taxonomic Level"] == level]
    metric_col = metric if metric in candidates.columns else f"mean_{metric}"
    if candidates.empty or metric_col not in candidates.columns:
        return pd.DataFrame()

    group_by = opt.get("group_by", ["dataset", "method", "gene_kind", "Taxonomic Level"])
    group_by = [col for col in group_by if col in candidates.columns]
    rows = []
    ascending = direction == "minimize"
    for keys, group in candidates.groupby(group_by, dropna=False):
        group = group.dropna(subset=[metric_col]).copy()
        if group.empty:
            continue
        group = group.sort_values(
            by=[metric_col, "taxa_after_mean", "pred_retained_abundance_mean"],
            ascending=[ascending, False, False],
        )
        best = group.iloc[0].to_dict()
        best["objective"] = opt.get("objective", "genus:F1")
        best["objective_metric"] = metric_col
        best["objective_direction"] = direction
        best["threshold_selection_mode"] = "per_method_optimized"
        rows.append(best)
    return pd.DataFrame(rows)


def build_method_comparison(summary: pd.DataFrame) -> pd.DataFrame:
    if summary.empty:
        return pd.DataFrame()
    group_cols = ["dataset", "method", "gene_kind", "analysis_mode", "Taxonomic Level"]
    return (
        summary.groupby(group_cols, dropna=False)
        .agg(
            mean_F1=("mean_F1", "mean"),
            mean_AbundanceF1=("mean_AbundanceF1", "mean"),
            mean_BC=("mean_BC", "mean"),
            mean_rJSD=("mean_rJSD", "mean"),
            config_count=("filter_config_id", "nunique"),
        )
        .reset_index()
    )
