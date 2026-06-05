"""Simplified output tables for BioProfileBench."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from .metrics import CORE_COLUMNS, EXTRA_COLUMNS


PUBLIC_DROP_COLUMNS = {"dataset"}
PATH_COLUMNS = {"prediction_path", "truth_path", "source_path"}
TAXONOMIC_LEVELS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]


SIMPLE_BENCHMARK_COLUMNS = CORE_COLUMNS + EXTRA_COLUMNS + [
    "File",
    "TruthFile",
    "PredFile",
    "GeneKind",
    "Status",
]


def basename_or_empty(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value)
    if not text:
        return ""
    return Path(text).name




def prepare_public_table(df: pd.DataFrame, drop_dataset: bool = True) -> pd.DataFrame:
    """Clean output tables for downstream use.

    Public result files should not expose absolute paths by default. Path-like
    metadata columns are kept, but values are reduced to file names. The dataset
    column is dropped because batch pairing is already encoded by truth/pred file
    names and filter IDs.
    """
    if df.empty:
        return df.copy()
    out = df.copy()
    if drop_dataset:
        out = out.drop(columns=[col for col in PUBLIC_DROP_COLUMNS if col in out.columns])
    for col in PATH_COLUMNS:
        if col in out.columns:
            out[col] = out[col].map(basename_or_empty)
    return out


def compact_filter_params(row: pd.Series) -> str:
    """Compact non-empty filter columns into a stable key=value string."""
    skip = {
        "Info",
        "profile_id",
        "Taxonomic Level",
        "analysis_mode",
        "filter_config_id",
        "threshold_source",
        "FilterParams",
    }
    parts = []
    for key in sorted(col for col in row.index if col not in skip):
        value = row[key]
        if pd.isna(value) or value == "":
            continue
        if isinstance(value, float):
            value_text = f"{value:g}"
        else:
            value_text = str(value)
        parts.append(f"{key}={value_text}")
    return "; ".join(parts) if parts else "none"


def build_simple_benchmark_table(benchmark: pd.DataFrame) -> pd.DataFrame:
    """Return the compact R-friendly sample x rank benchmark table."""
    if benchmark.empty:
        return pd.DataFrame(columns=SIMPLE_BENCHMARK_COLUMNS)

    out = benchmark.copy()
    out["Info"] = out.get("filter_config_id", out.get("Info", "raw")).astype(str)
    if "Taxonomic Level" in out.columns:
        out["Taxonomic Level"] = out["Taxonomic Level"].map(lambda value: str(value).capitalize())
    out["PredFile"] = out.get("prediction_path", "").map(basename_or_empty)
    out["TruthFile"] = out.get("truth_path", "").map(basename_or_empty)
    out["File"] = out["PredFile"]
    out["GeneKind"] = out.get("gene_kind", "")
    out["Status"] = out.get("Status", "")

    columns = [col for col in SIMPLE_BENCHMARK_COLUMNS if col in out.columns]
    return out[columns]


def build_filter_trace_table(filter_configs: pd.DataFrame) -> pd.DataFrame:
    """Map compact Info/filter IDs to actual filtering parameter values."""
    if filter_configs.empty:
        return pd.DataFrame(
            columns=[
                "Info",
                "Taxonomic Level",
                "analysis_mode",
                "filter_config_id",
                "threshold_source",
                "FilterParams",
            ]
        )

    out = filter_configs.copy().drop_duplicates()
    out["Info"] = out["filter_config_id"].astype(str)
    if "Taxonomic Level" in out.columns:
        out["Taxonomic Level"] = out["Taxonomic Level"].map(lambda value: str(value).capitalize())
    out["FilterParams"] = out.apply(compact_filter_params, axis=1)

    preferred = [
        "Info",
        "Taxonomic Level",
        "analysis_mode",
        "filter_config_id",
        "threshold_source",
        "FilterParams",
    ]
    filter_cols = [col for col in out.columns if col not in preferred and col != "profile_id"]
    ordered = [col for col in preferred if col in out.columns] + sorted(filter_cols)
    return out[ordered].drop_duplicates()


def _metric_mean(group: pd.DataFrame, column: str) -> float:
    return group[column].mean() if column in group.columns else pd.NA


def _build_ranking(
    simple_benchmark: pd.DataFrame,
    filter_trace: pd.DataFrame,
    group_cols: list[str],
) -> pd.DataFrame:
    if simple_benchmark.empty:
        return pd.DataFrame()

    usable = simple_benchmark[simple_benchmark.get("Status", "OK").isin(["OK", "EMPTY_PREDICTION", "EMPTY_BOTH", ""])].copy()
    if usable.empty:
        usable = simple_benchmark.copy()

    rows = []
    for keys, group in usable.groupby(group_cols, dropna=False):
        key_data = dict(zip(group_cols, keys if isinstance(keys, tuple) else (keys,)))
        rows.append(
            {
                **key_data,
                "mean_F1": _metric_mean(group, "F1"),
                "median_F1": group["F1"].median() if "F1" in group.columns else pd.NA,
                "mean_AbundanceF1": _metric_mean(group, "AbundanceF1"),
                "mean_Spearman": _metric_mean(group, "Spearman"),
                "mean_BC": _metric_mean(group, "BC"),
                "mean_rJSD": _metric_mean(group, "rJSD"),
                "mean_Precision": _metric_mean(group, "Precision"),
                "mean_Recall": _metric_mean(group, "Recall"),
                "mean_Jaccard": _metric_mean(group, "Jaccard"),
                "mean_WeightedJaccard": _metric_mean(group, "WeightedJaccard"),
                "mean_Pearson": _metric_mean(group, "Pearson"),
                "mean_Cosine": _metric_mean(group, "Cosine"),
                "mean_AUPRC": _metric_mean(group, "AUPRC"),
                "mean_Accuracy": _metric_mean(group, "Accuracy"),
                "mean_MatthewsCorrelationCoefficient": _metric_mean(group, "MatthewsCorrelationCoefficient"),
                "mean_MeanAbsoluteError": _metric_mean(group, "MeanAbsoluteError"),
                "mean_RootMeanSquaredError": _metric_mean(group, "RootMeanSquaredError"),
                "sample_count": int(group["Sample"].nunique()) if "Sample" in group.columns else int(len(group)),
                "row_count": int(len(group)),
            }
        )
    ranked = pd.DataFrame(rows)
    if ranked.empty:
        return ranked

    trace_cols = [col for col in ["Info", "FilterParams", "analysis_mode"] if col in filter_trace.columns]
    if trace_cols:
        trace = filter_trace[trace_cols].drop_duplicates(subset=["Info"])
        ranked = ranked.merge(trace, on="Info", how="left")

    sort_cols = [col for col in ["mean_F1", "mean_AbundanceF1", "mean_BC", "mean_rJSD"] if col in ranked.columns]
    ascending = [False, False, True, True][: len(sort_cols)]
    if sort_cols:
        ranked = ranked.sort_values(sort_cols, ascending=ascending, na_position="last")
    return ranked.reset_index(drop=True)




def build_benchmark_summary(simple_benchmark: pd.DataFrame, filter_trace: pd.DataFrame) -> pd.DataFrame:
    """Summarize benchmark results by truth file, prediction file, filter, and rank."""
    group_cols = [
        col
        for col in ["GeneKind", "TruthFile", "PredFile", "Info", "Taxonomic Level"]
        if col in simple_benchmark.columns
    ]
    summary = _build_ranking(simple_benchmark, filter_trace, group_cols)
    if summary.empty:
        return summary
    if "filter_config_id" not in summary.columns and "Info" in summary.columns:
        insert_at = summary.columns.get_loc("Info") + 1
        summary.insert(insert_at, "filter_config_id", summary["Info"])
    preferred = [
        "GeneKind",
        "TruthFile",
        "PredFile",
        "Info",
        "filter_config_id",
        "FilterParams",
        "analysis_mode",
        "Taxonomic Level",
    ]
    ordered = [col for col in preferred if col in summary.columns]
    ordered += [col for col in summary.columns if col not in ordered]
    return summary[ordered]


def split_summary_by_level(summary: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Split benchmark summary into standard taxonomic-rank tables."""
    if summary.empty or "Taxonomic Level" not in summary.columns:
        return {level: pd.DataFrame(columns=summary.columns) for level in TAXONOMIC_LEVELS}
    levels = summary["Taxonomic Level"].astype(str).str.lower()
    return {level: summary.loc[levels == level].copy() for level in TAXONOMIC_LEVELS}


def build_pair_filter_ranking(simple_benchmark: pd.DataFrame, filter_trace: pd.DataFrame) -> pd.DataFrame:
    """Rank each truth/pred/filter combination across all samples and ranks."""
    group_cols = [col for col in ["GeneKind", "TruthFile", "PredFile", "Info"] if col in simple_benchmark.columns]
    return _build_ranking(simple_benchmark, filter_trace, group_cols)


def build_pair_filter_level_ranking(simple_benchmark: pd.DataFrame, filter_trace: pd.DataFrame) -> pd.DataFrame:
    """Rank each truth/pred/filter combination separately for each taxonomic rank."""
    group_cols = [
        col
        for col in ["GeneKind", "TruthFile", "PredFile", "Info", "Taxonomic Level"]
        if col in simple_benchmark.columns
    ]
    return _build_ranking(simple_benchmark, filter_trace, group_cols)
