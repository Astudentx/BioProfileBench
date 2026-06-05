"""Input/output helpers for BioProfileBench."""

from __future__ import annotations

import ast
import hashlib
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd

from .bk_metadata import parse_bk_file_metadata
from .profile import Profile


ID_COLUMN = "ID"
ID_ALIASES = {ID_COLUMN, "TaxID", "SpeciesID", "FeatureID", "GeneID"}
TAXONOMY_COLUMN = "Taxonomy"
LINEAGE_COLUMN = "Lineage"
TAXID_COLUMN = "TaxID"
RESERVED_COLUMNS = {ID_COLUMN, TAXONOMY_COLUMN, LINEAGE_COLUMN, TAXID_COLUMN}


def stable_hash(data: Any, length: int = 12) -> str:
    """Return a stable short hash for config dictionaries and run IDs."""
    payload = json.dumps(data, sort_keys=True, ensure_ascii=False, default=str)
    return hashlib.md5(payload.encode("utf-8")).hexdigest()[:length]


def parse_scalar(value: str) -> Any:
    """Parse scalar values from a small YAML subset."""
    value = value.strip()
    lower = value.lower()
    if lower in {"null", "none", "na", ""}:
        return None
    if lower == "true":
        return True
    if lower == "false":
        return False
    if value.startswith("[") and value.endswith("]"):
        inner = value[1:-1].strip()
        if not inner:
            return []
        parts = [part.strip() for part in inner.split(",")]
        return [parse_scalar(part) for part in parts]
    try:
        return ast.literal_eval(value)
    except Exception:
        return value.strip("'\"")


def strip_comment(line: str) -> str:
    """Strip comments from simple YAML lines."""
    in_single = False
    in_double = False
    for i, char in enumerate(line):
        if char == "'" and not in_double:
            in_single = not in_single
        elif char == '"' and not in_single:
            in_double = not in_double
        elif char == "#" and not in_single and not in_double:
            return line[:i]
    return line


def load_simple_yaml(path: str | Path) -> dict[str, Any]:
    """Load a lightweight YAML subset without requiring PyYAML.

    Supported syntax is enough for this tool's configs:
    nested mappings by indentation and inline lists such as [a, b, 1e-5].
    """
    path = Path(path)
    root: dict[str, Any] = {}
    stack: list[tuple[int, dict[str, Any]]] = [(-1, root)]

    for line_no, raw_line in enumerate(path.read_text().splitlines(), start=1):
        line = strip_comment(raw_line).rstrip()
        if not line.strip():
            continue
        if re.match(r"^\s*-", line):
            raise ValueError(
                f"{path}:{line_no}: block lists are not supported by the built-in parser; "
                "use inline lists like [a, b, c]."
            )
        indent = len(line) - len(line.lstrip(" "))
        text = line.strip()
        if ":" not in text:
            raise ValueError(f"{path}:{line_no}: expected 'key: value' syntax")
        key, value = text.split(":", 1)
        key = key.strip()
        value = value.strip()

        while stack and indent <= stack[-1][0]:
            stack.pop()
        current = stack[-1][1]
        if not value:
            new_dict: dict[str, Any] = {}
            current[key] = new_dict
            stack.append((indent, new_dict))
        else:
            current[key] = parse_scalar(value)

    return root


def load_config(path: str | Path) -> dict[str, Any]:
    """Load JSON or the built-in YAML subset."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config file does not exist: {path}")
    if path.suffix.lower() == ".json":
        return json.loads(path.read_text())
    return load_simple_yaml(path)


def read_abundance_table(path: str | Path) -> pd.DataFrame:
    """Read and validate a tab-separated abundance table."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Input file does not exist: {path}")
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if df.empty:
        raise ValueError(f"Input table is empty: {path}")
    first_column = df.columns[0]
    if first_column not in ID_ALIASES:
        raise ValueError(
            f"First column must be one of {sorted(ID_ALIASES)} in {path}, found '{first_column}'"
        )
    if first_column != ID_COLUMN:
        df = df.rename(columns={first_column: ID_COLUMN})
    if df.columns.duplicated().any():
        duplicated = sorted(set(df.columns[df.columns.duplicated()].tolist()))
        raise ValueError(f"Duplicated columns in {path}: {duplicated}")
    return df


def detect_sample_columns(df: pd.DataFrame) -> list[str]:
    """Detect abundance sample columns."""
    sample_cols = [col for col in df.columns if col not in RESERVED_COLUMNS]
    if not sample_cols:
        raise ValueError("No sample abundance columns detected")
    return sample_cols


def coerce_abundance(df: pd.DataFrame, sample_cols: list[str], path: str | Path) -> pd.DataFrame:
    """Convert abundance columns to numeric values."""
    df = df.copy()
    for col in sample_cols:
        converted = pd.to_numeric(df[col], errors="coerce")
        invalid = converted.isna() & ~df[col].astype(str).str.strip().isin(["", "NA", "nan", "None"])
        if invalid.any():
            print(
                f"WARNING: {path}: column '{col}' has {int(invalid.sum())} non-numeric values; converted to 0.",
                flush=True,
            )
        converted = converted.fillna(0.0).astype(float)
        if (converted < 0).any():
            print(
                f"WARNING: {path}: column '{col}' has negative values; clipped to 0.",
                flush=True,
            )
            converted = converted.clip(lower=0.0)
        df[col] = converted
    return df


def build_profile(
    path: str | Path,
    profile_id: str,
    dataset: str = "",
    method: str = "",
    gene_kind: str = "",
    parameter_tag: str = "",
    extra_metadata: dict[str, Any] | None = None,
) -> Profile:
    """Read a table and convert it into a Profile."""
    path = Path(path)
    df = read_abundance_table(path)
    sample_cols = detect_sample_columns(df)
    df = coerce_abundance(df, sample_cols, path)

    abundance = df[[ID_COLUMN] + sample_cols].set_index(ID_COLUMN)
    taxonomy_cols = [col for col in [ID_COLUMN, TAXONOMY_COLUMN, LINEAGE_COLUMN, TAXID_COLUMN] if col in df.columns]
    taxonomy = df[taxonomy_cols].set_index(ID_COLUMN, drop=False)

    feature_metadata = df.drop(columns=sample_cols).set_index(ID_COLUMN, drop=False)
    merged_extra = parse_bk_file_metadata(path)
    if extra_metadata:
        merged_extra.update({key: value for key, value in extra_metadata.items() if value != ""})
    merged_extra.setdefault("source_path", str(path))
    merged_extra.setdefault("prediction_path", str(path))
    merged_extra.setdefault("Method", method)
    return Profile(
        profile_id=profile_id,
        dataset=dataset,
        method=method,
        gene_kind=gene_kind,
        parameter_tag=parameter_tag,
        abundance_matrix=abundance,
        taxonomy_table=taxonomy,
        feature_metadata=feature_metadata,
        source_path=path,
        extra_metadata=merged_extra,
    )


def load_manifest(path: str | Path) -> pd.DataFrame:
    """Load batch manifest and resolve relative paths."""
    path = Path(path)
    manifest = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = ["profile_id", "dataset", "method", "gene_kind", "prediction_path", "truth_path"]
    missing = [col for col in required if col not in manifest.columns]
    if missing:
        raise ValueError(f"Manifest missing required columns: {missing}")
    if manifest["profile_id"].duplicated().any():
        duplicated = manifest.loc[manifest["profile_id"].duplicated(), "profile_id"].tolist()
        raise ValueError(f"Manifest profile_id values must be unique; duplicated: {duplicated}")

    base = path.parent
    for col in ["prediction_path", "truth_path"]:
        manifest[col] = manifest[col].map(lambda value: str((base / value).resolve()) if value and not Path(value).is_absolute() else value)
        missing_paths = [value for value in manifest[col] if not Path(value).exists()]
        if missing_paths:
            raise FileNotFoundError(f"Manifest column '{col}' has missing files: {missing_paths[:5]}")
    return manifest


def write_table(df: pd.DataFrame, path: str | Path, digits: int = 6) -> None:
    """Write a TSV table with consistent NA and float formatting."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, na_rep="NA", float_format=f"%.{digits}f")
