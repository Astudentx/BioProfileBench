"""Utilities for existing BK benchmark result files."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from .bk_metadata import build_r_compatible_table, parse_bk_file_metadata
from .io import write_table


METHOD_BY_DIR = {
    "02.BP2": "BP-Tracer",
    "03.Reads": "Reads",
    "04.Contigs": "Contigs",
    "05.MAGs": "MAGs",
    "06.Cor": "Correlation",
}


def infer_method_from_path(path: str | Path) -> str:
    parts = Path(path).parts
    for part in parts:
        if part in METHOD_BY_DIR:
            return METHOD_BY_DIR[part]
    return ""


def collect_bk_files(inputs: list[str]) -> list[Path]:
    files: list[Path] = []
    for item in inputs:
        path = Path(item)
        if path.is_dir():
            files.extend(sorted(path.rglob("*.BK")))
        elif path.is_file():
            files.append(path)
        else:
            raise FileNotFoundError(f"Input path does not exist: {path}")
    return files


def load_existing_bk_results(paths: list[Path]) -> pd.DataFrame:
    tables = []
    for path in paths:
        df = pd.read_csv(path, sep="\t")
        if "Taxonomic Level" not in df.columns and "Level" in df.columns:
            df = df.rename(columns={"Level": "Taxonomic Level"})
        metadata = parse_bk_file_metadata(path.name)
        metadata["prediction_path"] = str(path)
        metadata["source_path"] = str(path)
        method = infer_method_from_path(path)
        metadata["method"] = method
        metadata["Method"] = method
        for key, value in metadata.items():
            df[key] = value
        tables.append(df)
    return pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()


def format_existing_bk_results(args: Any, config: dict[str, Any] | None = None) -> Path:
    """Convert old .BK benchmark outputs into AllBK_filter-compatible TSV."""
    config = config or {}
    report_config = config.get("report", {})
    sample_config = config.get("sample_metadata", {})
    group_size_map = report_config.get("group_size_map", sample_config.get("group_size_map", {}))
    files = collect_bk_files(args.input)
    benchmark = load_existing_bk_results(files)
    output = build_r_compatible_table(benchmark, group_size_map)
    digits = int(config.get("benchmark", {}).get("digits", 6))
    write_table(output, args.out, digits)
    return Path(args.out)
