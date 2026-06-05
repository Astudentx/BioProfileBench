"""BK project metadata parsing and R-compatible output helpers."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pandas as pd


BK_FILE_COLUMNS = [
    "File",
    "GeneKind",
    "GeneIdentify",
    "GeneLength",
    "GeneEvalue",
    "TrueValueIdentify",
    "TrueValueCoverage",
    "TrueValueType",
    "filter1",
    "filter2",
    "filter3",
    "filter4",
]

R_COMPAT_METRIC_COLUMNS = [
    "Sample",
    "Level",
    "Info",
    "TP",
    "FP",
    "FN",
    "TN",
    "Precision",
    "Recall",
    "F1",
    "Spearman",
    "AbundancePrecision",
    "AbundanceRecall",
    "AbundanceF1",
    "L1",
    "L2",
    "BC",
    "rJSD",
]

R_COMPAT_METADATA_COLUMNS = [
    "File",
    "GeneKind",
    "GeneIdentify",
    "GeneLength",
    "GeneEvalue",
    "TrueValueIdentify",
    "TrueValueCoverage",
    "TrueValueType",
    "filter1",
    "filter2",
    "filter3",
    "filter4",
    "Method",
    "Group",
    "rep",
    "SpeciesNumbers",
    "StrainNumbers",
]

R_COMPAT_COLUMNS = R_COMPAT_METRIC_COLUMNS + R_COMPAT_METADATA_COLUMNS


FILE_PATTERN = re.compile(
    r"^(?P<GeneKind>[A-Za-z]+)"
    r"\.i(?P<GeneIdentify>\d+)"
    r"\.l(?P<GeneLength>\d+)"
    r"\.(?P<GeneEvalue>[^_]+)"
    r"_i(?P<TrueValueIdentify>\d+)c(?P<TrueValueCoverage>\d+)"
    r"_?(?P<TrueValueType>[A-Za-z0-9]*)"
    r"(?P<FilterPart>(?:\.f\d+_[A-Za-z0-9eE+\-]+(?:\.\d+)?)*)"
)


SAMPLE_PATTERN = re.compile(r"^(?P<Group>[A-Za-z]+\d+)-(?P<rep0>\d+)$")


def parse_bk_file_metadata(path_or_name: str | Path) -> dict[str, str]:
    """Parse BK-style file names into plotting metadata."""
    name = Path(path_or_name).name
    parsed_name = name[:-3] if name.endswith(".BK") else name
    match = FILE_PATTERN.match(parsed_name)
    data = {column: "" for column in BK_FILE_COLUMNS}
    data["File"] = name
    if not match:
        return data
    matched = {key: value or "" for key, value in match.groupdict().items()}
    filter_part = matched.pop("FilterPart", "")
    data.update(matched)
    filter_map = {
        "f4": "filter1",
        "f5": "filter2",
        "f1": "filter3",
        "f6": "filter4",
    }
    for token in re.findall(r"(f\d+_[A-Za-z0-9eE+\-]+(?:\.\d+)?)", filter_part):
        prefix = token.split("_", 1)[0]
        column = filter_map.get(prefix)
        if column:
            data[column] = token
    return data


def filter_info_from_metadata(metadata: dict[str, Any]) -> str:
    """Build legacy Info such as f4=0.001_f5=0_f1=0_f6=0 from parsed filters."""
    parts = []
    for key in ["filter1", "filter2", "filter3", "filter4"]:
        value = metadata.get(key, "")
        if not value:
            continue
        text = str(value).rstrip(".")
        if "_" in text:
            name, threshold = text.split("_", 1)
            parts.append(f"{name}={threshold}")
        else:
            parts.append(text)
    return "_".join(parts)


def parse_sample_metadata(sample: str, group_size_map: dict[str, Any] | None = None) -> dict[str, Any]:
    """Parse sample names like T3-0 into Group=T3 and rep=1."""
    group_size_map = group_size_map or {}
    match = SAMPLE_PATTERN.match(str(sample))
    if match:
        group = match.group("Group")
        rep = int(match.group("rep0")) + 1
    else:
        group = ""
        rep = ""

    group_sizes = group_size_map.get(group, {}) if isinstance(group_size_map, dict) else {}
    return {
        "Group": group,
        "rep": rep,
        "SpeciesNumbers": group_sizes.get("SpeciesNumbers", ""),
        "StrainNumbers": group_sizes.get("StrainNumbers", ""),
    }


def capitalized_level(level: Any) -> str:
    text = str(level)
    return text[:1].upper() + text[1:].lower() if text else text


def build_r_compatible_table(
    benchmark: pd.DataFrame,
    group_size_map: dict[str, Any] | None = None,
) -> pd.DataFrame:
    """Convert benchmark_long into the legacy AllBK_filter-like plotting table."""
    if benchmark.empty:
        return pd.DataFrame(columns=R_COMPAT_COLUMNS)

    rows = []
    for _, row in benchmark.iterrows():
        metadata = {column: row.get(column, "") for column in BK_FILE_COLUMNS}
        if not metadata.get("File"):
            metadata.update(parse_bk_file_metadata(row.get("prediction_path", row.get("source_path", row.get("profile_id", "")))))
        for column in BK_FILE_COLUMNS:
            if not metadata.get(column) and column in row:
                metadata[column] = row.get(column, "")

        sample_meta = parse_sample_metadata(row.get("Sample", ""), group_size_map)
        legacy_info = filter_info_from_metadata(metadata) or row.get("Info", "")
        out = {
            "Sample": row.get("Sample", ""),
            "Level": capitalized_level(row.get("Taxonomic Level", row.get("Level", ""))),
            "Info": legacy_info,
            "TP": row.get("TP", ""),
            "FP": row.get("FP", ""),
            "FN": row.get("FN", ""),
            "TN": row.get("TN", ""),
            "Precision": row.get("Precision", ""),
            "Recall": row.get("Recall", ""),
            "F1": row.get("F1", ""),
            "Spearman": row.get("Spearman", ""),
            "AbundancePrecision": row.get("AbundancePrecision", ""),
            "AbundanceRecall": row.get("AbundanceRecall", ""),
            "AbundanceF1": row.get("AbundanceF1", ""),
            "L1": row.get("L1", ""),
            "L2": row.get("L2", ""),
            "BC": row.get("BC", ""),
            "rJSD": row.get("rJSD", ""),
            "Method": row.get("Method", row.get("method", "")),
            **metadata,
            **sample_meta,
        }
        rows.append(out)

    return pd.DataFrame(rows, columns=R_COMPAT_COLUMNS)
