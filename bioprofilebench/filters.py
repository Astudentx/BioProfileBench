"""Prediction-only filtering engine."""

from __future__ import annotations

import itertools
import math
from numbers import Number
from typing import Any

import numpy as np
import pandas as pd

from .io import stable_hash


FILTER_KEYS = [
    "pred_min_abundance",
    "pred_min_relative_abundance",
    "pred_min_total_abundance",
    "pred_min_mean_abundance",
    "pred_min_max_abundance",
    "pred_min_prevalence",
    "pred_top_k",
    "pred_top_percent",
    "pred_cumulative_abundance",
    "top_mode",
]

ZERO_NOOP_FILTERS = {
    "pred_min_abundance",
    "pred_min_relative_abundance",
    "pred_min_total_abundance",
    "pred_min_mean_abundance",
    "pred_min_max_abundance",
    "pred_min_prevalence",
}




def normalize_filter_value(value: Any) -> str:
    """Normalize filter values so equivalent numeric thresholds share one code."""
    if value is None or value == "":
        return ""
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, Number):
        return f"{float(value):.12g}"
    text = str(value).strip()
    try:
        return f"{float(text):.12g}"
    except ValueError:
        return text.lower() if text.lower() in {"global", "per_sample"} else text


def canonical_filter_payload(filters: dict[str, Any]) -> dict[str, str]:
    """Return semantic filter parameters used for stable filter codes.

    The code intentionally ignores analysis mode, file path, taxonomic level,
    folder, and batch. It also removes no-op defaults, so equivalent filters such
    as 0, 0.0, omitted top_mode, and top_mode=global share the same code.
    """
    payload: dict[str, str] = {}
    top_filters = {"pred_top_k", "pred_top_percent", "pred_cumulative_abundance"}
    for key in FILTER_KEYS:
        if key == "top_mode" and not any(filters.get(name) is not None for name in top_filters):
            continue
        value = filters.get(key)
        if value is None or value == "":
            continue
        normalized = normalize_filter_value(value)
        if key in ZERO_NOOP_FILTERS:
            try:
                if float(normalized) <= 0:
                    continue
            except ValueError:
                pass
        payload[key] = normalized
    if any(filters.get(name) is not None for name in top_filters):
        payload.setdefault("top_mode", "global")
    return payload


def make_filter_config_id(filters: dict[str, Any]) -> str:
    """Build a short deterministic filter ID from semantic filter parameters."""
    payload = canonical_filter_payload(filters)
    if not payload:
        return "F000"
    return f"F{stable_hash(payload, length=8).upper()}"


def as_list(value: Any) -> list[Any]:
    if value is None:
        return [None]
    if isinstance(value, list):
        return value
    return [value]


def generate_dynamic_values(name: str, spec: dict[str, Any], matrix: pd.DataFrame) -> list[Any]:
    """Generate data-driven threshold values for one filter."""
    if name == "pred_min_relative_abundance":
        column_totals = matrix.sum(axis=0)
        matrix = matrix.copy().astype(float)
        nonzero_columns = column_totals > 0
        matrix.loc[:, nonzero_columns] = matrix.loc[:, nonzero_columns].div(
            column_totals[nonzero_columns],
            axis=1,
        )
    values = matrix.to_numpy(dtype=float).ravel()
    nonzero = values[values > 0]
    if len(nonzero) == 0:
        return [0]

    strategy = spec.get("strategy", "quantile")
    if strategy == "quantile":
        quantiles = as_list(spec.get("quantiles", [0.5, 0.75, 0.9, 0.95]))
        return sorted(set(float(np.quantile(nonzero, q)) for q in quantiles))
    if strategy == "logspace":
        count = int(spec.get("count", 5))
        lo = float(spec.get("min", nonzero.min()))
        hi = float(spec.get("max", nonzero.max()))
        if lo <= 0 or hi <= 0 or np.isclose(lo, hi):
            return sorted(set([lo, hi]))
        return [float(x) for x in np.logspace(np.log10(lo), np.log10(hi), count)]
    if strategy == "target_taxa_count":
        targets = [int(x) for x in as_list(spec.get("targets", [10, 50, 100]))]
        totals = matrix.sum(axis=1).sort_values(ascending=False)
        out = []
        for target in targets:
            if target <= 0:
                out.append(float(totals.max()) + 1.0)
            elif target >= len(totals):
                out.append(0.0)
            else:
                out.append(float(totals.iloc[target]))
        return sorted(set(out))
    if strategy == "target_retained_abundance":
        targets = [float(x) for x in as_list(spec.get("targets", [0.95, 0.99]))]
        totals = matrix.sum(axis=1).sort_values(ascending=False)
        total = totals.sum()
        if total <= 0:
            return [0]
        cumulative = totals.cumsum() / total
        out = []
        for target in targets:
            idx = int(np.searchsorted(cumulative.to_numpy(), target, side="left"))
            idx = min(idx, len(totals) - 1)
            out.append(float(totals.iloc[idx]))
        return sorted(set(out))
    if strategy == "prevalence_sweep":
        max_prev = matrix.shape[1]
        return list(range(0, max_prev + 1))
    raise ValueError(f"Unsupported dynamic threshold strategy for {name}: {strategy}")


def build_filter_config_grid(
    config: dict[str, Any],
    prediction_matrix: pd.DataFrame | None = None,
) -> list[dict[str, Any]]:
    """Build filter configurations for raw/common/optimized analyses."""
    filtering = config.get("filtering", {})
    modes = as_list(filtering.get("modes", ["raw"]))
    common_filters = dict(filtering.get("common_filters", {}))

    if prediction_matrix is not None:
        for name, spec in dict(filtering.get("dynamic_filters", {})).items():
            if not isinstance(spec, dict):
                raise ValueError(f"dynamic_filters.{name} must be a mapping")
            common_filters[name] = generate_dynamic_values(name, spec, prediction_matrix)

    filter_names = list(common_filters.keys())
    filter_values = [as_list(common_filters[name]) for name in filter_names]
    if not filter_names:
        filter_combinations = [{}]
    else:
        filter_combinations = [
            {name: value for name, value in zip(filter_names, values)}
            for values in itertools.product(*filter_values)
        ]

    configs: list[dict[str, Any]] = []
    for mode in modes:
        if mode == "raw":
            filters = {}
            configs.append(
                {
                    "analysis_mode": "raw",
                    "filter_config_id": make_filter_config_id(filters),
                    "filters": filters,
                    "threshold_source": "none",
                }
            )
            continue
        for filters in filter_combinations:
            configs.append(
                {
                    "analysis_mode": mode,
                    "filter_config_id": make_filter_config_id(filters),
                    "filters": filters,
                    "threshold_source": "dynamic_or_grid" if filtering.get("dynamic_filters") else "grid",
                }
            )
    deduplicated: list[dict[str, Any]] = []
    seen_ids: set[str] = set()
    for filter_config in configs:
        filter_id = filter_config["filter_config_id"]
        if filter_id in seen_ids:
            continue
        seen_ids.add(filter_id)
        deduplicated.append(filter_config)
    return deduplicated


def prevalence_to_count(value: float, sample_count: int) -> int:
    if value <= 0:
        return 0
    if value < 1:
        return int(math.ceil(value * sample_count))
    return int(math.ceil(value))


def top_k_cell_mask(matrix: pd.DataFrame, k: int) -> pd.DataFrame:
    mask = pd.DataFrame(False, index=matrix.index, columns=matrix.columns)
    if k is None or k <= 0:
        return pd.DataFrame(True, index=matrix.index, columns=matrix.columns)
    for sample in matrix.columns:
        top_index = matrix[sample].sort_values(ascending=False).head(k).index
        mask.loc[top_index, sample] = True
    return mask


def top_percent_count(percent: float, taxa_count: int) -> int:
    if percent is None:
        return taxa_count
    if percent <= 0:
        return 0
    if percent <= 1:
        return max(1, int(math.ceil(percent * taxa_count)))
    return max(1, int(math.ceil(percent)))


def cumulative_cell_mask(matrix: pd.DataFrame, fraction: float) -> pd.DataFrame:
    mask = pd.DataFrame(False, index=matrix.index, columns=matrix.columns)
    if fraction is None:
        return pd.DataFrame(True, index=matrix.index, columns=matrix.columns)
    for sample in matrix.columns:
        values = matrix[sample].sort_values(ascending=False)
        total = values.sum()
        if total <= 0:
            continue
        cumulative = values.cumsum() / total
        keep = cumulative <= fraction
        if not keep.any() and len(keep):
            keep.iloc[0] = True
        elif (~keep).any():
            first_over = keep[~keep].index[0]
            keep.loc[first_over] = True
        mask.loc[values.index[keep], sample] = True
    return mask


def apply_prediction_filters(
    pred_matrix: pd.DataFrame,
    filter_config: dict[str, Any],
    presence_threshold: float = 0.0,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Apply prediction-only filters using independently built masks."""
    original = pred_matrix.copy().astype(float)
    column_totals = original.sum(axis=0)
    relative_original = original.copy()
    nonzero_columns = column_totals > 0
    relative_original.loc[:, nonzero_columns] = relative_original.loc[:, nonzero_columns].div(
        column_totals[nonzero_columns],
        axis=1,
    )
    filters = filter_config.get("filters", {})
    top_mode = filters.get("top_mode", "global") or "global"
    total_abundance = float(original.sum().sum())

    cell_masks: dict[str, pd.DataFrame] = {}
    row_masks: dict[str, pd.Series] = {}

    threshold = filters.get("pred_min_abundance")
    if threshold is not None:
        cell_masks["pred_min_abundance"] = original > float(threshold)

    threshold = filters.get("pred_min_relative_abundance")
    if threshold is not None:
        cell_masks["pred_min_relative_abundance"] = relative_original > float(threshold)

    threshold = filters.get("pred_min_total_abundance")
    if threshold is not None:
        row_masks["pred_min_total_abundance"] = original.sum(axis=1) >= float(threshold)

    threshold = filters.get("pred_min_mean_abundance")
    if threshold is not None:
        row_masks["pred_min_mean_abundance"] = original.mean(axis=1) >= float(threshold)

    threshold = filters.get("pred_min_max_abundance")
    if threshold is not None:
        row_masks["pred_min_max_abundance"] = original.max(axis=1) >= float(threshold)

    threshold = filters.get("pred_min_prevalence")
    if threshold is not None:
        min_count = prevalence_to_count(float(threshold), original.shape[1])
        row_masks["pred_min_prevalence"] = (original > presence_threshold).sum(axis=1) >= min_count

    top_k = filters.get("pred_top_k")
    if top_k is not None:
        top_k = int(top_k)
        if top_mode == "per_sample":
            cell_masks["pred_top_k"] = top_k_cell_mask(original, top_k)
        else:
            keep = original.sum(axis=1).sort_values(ascending=False).head(top_k).index
            row_masks["pred_top_k"] = pd.Series(original.index.isin(keep), index=original.index)

    top_percent = filters.get("pred_top_percent")
    if top_percent is not None:
        count = top_percent_count(float(top_percent), len(original.index))
        if top_mode == "per_sample":
            cell_masks["pred_top_percent"] = top_k_cell_mask(original, count)
        else:
            keep = original.sum(axis=1).sort_values(ascending=False).head(count).index
            row_masks["pred_top_percent"] = pd.Series(original.index.isin(keep), index=original.index)

    cumulative = filters.get("pred_cumulative_abundance")
    if cumulative is not None:
        cumulative = float(cumulative)
        if top_mode == "per_sample":
            cell_masks["pred_cumulative_abundance"] = cumulative_cell_mask(original, cumulative)
        else:
            totals = original.sum(axis=1).sort_values(ascending=False)
            total = totals.sum()
            keep = pd.Series(False, index=original.index)
            if total > 0:
                cumulative_fraction = totals.cumsum() / total
                selected = cumulative_fraction <= cumulative
                if not selected.any() and len(selected):
                    selected.iloc[0] = True
                elif (~selected).any():
                    selected.loc[selected[~selected].index[0]] = True
                keep.loc[totals.index[selected]] = True
            row_masks["pred_cumulative_abundance"] = keep

    final_cell = pd.DataFrame(True, index=original.index, columns=original.columns)
    for mask in cell_masks.values():
        final_cell &= mask

    final_row = pd.Series(True, index=original.index)
    for mask in row_masks.values():
        final_row &= mask

    filtered = original.where(final_cell, 0.0)
    filtered = filtered.loc[final_row]
    filtered = filtered.loc[filtered.sum(axis=1) > 0]

    row_pass_masks: dict[str, pd.Series] = {}
    for name, mask in row_masks.items():
        row_pass_masks[name] = mask.reindex(original.index, fill_value=False)
    for name, mask in cell_masks.items():
        row_pass_masks[name] = mask.any(axis=1).reindex(original.index, fill_value=False)

    trace_rows = []
    for name, row_pass in row_pass_masks.items():
        other_masks = [mask for other, mask in row_pass_masks.items() if other != name]
        unique_fail = ~row_pass
        for other_mask in other_masks:
            unique_fail &= other_mask

        if name in row_masks:
            retained = original.loc[row_pass].sum().sum()
        else:
            retained = original.where(cell_masks[name], 0.0).sum().sum()
        trace_rows.append(
            {
                "filter_name": name,
                "threshold": filters.get(name),
                "taxa_before": len(original),
                "taxa_pass": int(row_pass.sum()),
                "taxa_fail": int((~row_pass).sum()),
                "unique_fail_taxa": int(unique_fail.sum()),
                "pred_retained_abundance": np.nan if total_abundance <= 0 else float(retained / total_abundance),
            }
        )

    final_row_pass = pd.Series(original.index.isin(filtered.index), index=original.index)
    trace_rows.append(
        {
            "filter_name": "FINAL",
            "threshold": "",
            "taxa_before": len(original),
            "taxa_pass": int(final_row_pass.sum()),
            "taxa_fail": int((~final_row_pass).sum()),
            "unique_fail_taxa": int((~final_row_pass).sum()),
            "pred_retained_abundance": np.nan
            if total_abundance <= 0
            else float(filtered.sum().sum() / total_abundance),
        }
    )

    diagnostics = {
        "taxa_before": len(original),
        "taxa_after": len(filtered),
        "pred_total_before": total_abundance,
        "pred_total_after": float(filtered.sum().sum()),
        "pred_retained_abundance": np.nan
        if total_abundance <= 0
        else float(filtered.sum().sum() / total_abundance),
    }
    return filtered, pd.DataFrame(trace_rows), diagnostics


def filter_config_table(filter_configs: list[dict[str, Any]]) -> pd.DataFrame:
    rows = []
    for config in filter_configs:
        row = {
            "analysis_mode": config["analysis_mode"],
            "filter_config_id": config["filter_config_id"],
            "threshold_source": config.get("threshold_source", "grid"),
        }
        row.update(config.get("filters", {}))
        rows.append(row)
    return pd.DataFrame(rows)
