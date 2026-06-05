"""Taxonomy parsing and rank aggregation."""

from __future__ import annotations

import re
from typing import Any

import pandas as pd

from .io import LINEAGE_COLUMN, TAXONOMY_COLUMN
from .profile import Profile


UNCLASSIFIED = "Unclassified"
LEVEL_ORDER = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
LEVEL_ALIASES = {
    "k": "kingdom",
    "kingdom": "kingdom",
    "domain": "kingdom",
    "d": "kingdom",
    "p": "phylum",
    "phylum": "phylum",
    "c": "class",
    "class": "class",
    "o": "order",
    "order": "order",
    "f": "family",
    "family": "family",
    "g": "genus",
    "genus": "genus",
    "s": "species",
    "species": "species",
}
PREFIX_TO_LEVEL = {
    "k": "kingdom",
    "d": "kingdom",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}
MISSING_TOKENS = {"", "na", "n/a", "nan", "none", "null", "-", "unclassified"}


def is_missing(value: Any) -> bool:
    if value is None:
        return True
    try:
        if pd.isna(value):
            return True
    except TypeError:
        pass
    return str(value).strip().lower() in MISSING_TOKENS


def clean_taxon(value: Any, prefix_numeric_taxid: bool = False) -> str:
    if is_missing(value):
        return UNCLASSIFIED
    text = re.sub(r"\s+", " ", str(value).strip())
    if is_missing(text):
        return UNCLASSIFIED
    if prefix_numeric_taxid and re.fullmatch(r"\d+", text):
        return f"taxid_{text}"
    return text


def normalize_level(level: str) -> str:
    key = str(level).strip().lower()
    if key not in LEVEL_ALIASES:
        valid = ", ".join(LEVEL_ORDER + ["k", "p", "c", "o", "f", "g", "s"])
        raise ValueError(f"Unsupported taxonomic level '{level}'. Valid: {valid}")
    return LEVEL_ALIASES[key]


def normalize_levels(levels: list[str] | None) -> list[str]:
    if not levels:
        return list(LEVEL_ORDER)
    out: list[str] = []
    seen = set()
    for level in levels:
        canonical = normalize_level(level)
        if canonical not in seen:
            out.append(canonical)
            seen.add(canonical)
    return out


def parse_prefixed_taxonomy(text: Any) -> dict[str, str]:
    if is_missing(text):
        return {}
    parsed: dict[str, str] = {}
    for part in re.split(r"[;|]", str(text)):
        token = part.strip()
        if not token:
            continue
        short = re.match(r"^([A-Za-z])\s*(?:__|:)\s*(.*)$", token)
        long = re.match(
            r"^(kingdom|domain|phylum|class|order|family|genus|species)\s*(?:__|:)\s*(.*)$",
            token,
            flags=re.IGNORECASE,
        )
        if short:
            level = PREFIX_TO_LEVEL.get(short.group(1).lower())
            name = short.group(2)
        elif long:
            raw_level = long.group(1).lower()
            level = "kingdom" if raw_level == "domain" else raw_level
            name = long.group(2)
        else:
            continue
        if level:
            parsed[level] = clean_taxon(name)
    return parsed


def parse_lineage(lineage: Any, level: str) -> str:
    """Parse prefixed lineage or plain taxid lineage."""
    level = normalize_level(level)
    if is_missing(lineage):
        return UNCLASSIFIED
    prefixed = parse_prefixed_taxonomy(lineage)
    if prefixed:
        return prefixed.get(level, UNCLASSIFIED)

    parts = [part.strip() for part in re.split(r"[;|]", str(lineage))]
    idx = LEVEL_ORDER.index(level)
    if idx >= len(parts):
        return UNCLASSIFIED
    return clean_taxon(parts[idx], prefix_numeric_taxid=True)


def parse_taxonomy(taxonomy: Any, level: str) -> str:
    """Parse Taxonomy as prefixed ranks or species-level label."""
    level = normalize_level(level)
    if is_missing(taxonomy):
        return UNCLASSIFIED
    prefixed = parse_prefixed_taxonomy(taxonomy)
    if prefixed:
        return prefixed.get(level, UNCLASSIFIED)
    return clean_taxon(taxonomy) if level == "species" else UNCLASSIFIED


def column_has_values(profile: Profile, column: str) -> bool:
    return column in profile.taxonomy_table.columns and profile.taxonomy_table[column].map(lambda x: not is_missing(x)).any()


def choose_taxonomy_source(truth: Profile, prediction: Profile, requested: str = "auto") -> str:
    requested = requested or "auto"
    if requested not in {"auto", LINEAGE_COLUMN, TAXONOMY_COLUMN}:
        raise ValueError("taxonomy.source must be one of: auto, Lineage, Taxonomy")
    shared_lineage = column_has_values(truth, LINEAGE_COLUMN) and column_has_values(prediction, LINEAGE_COLUMN)
    shared_taxonomy = column_has_values(truth, TAXONOMY_COLUMN) and column_has_values(prediction, TAXONOMY_COLUMN)
    if requested == LINEAGE_COLUMN:
        if not shared_lineage:
            raise ValueError("Requested Lineage, but both profiles do not contain usable Lineage values")
        return LINEAGE_COLUMN
    if requested == TAXONOMY_COLUMN:
        if not shared_taxonomy:
            raise ValueError("Requested Taxonomy, but both profiles do not contain usable Taxonomy values")
        return TAXONOMY_COLUMN
    if shared_lineage:
        return LINEAGE_COLUMN
    if shared_taxonomy:
        return TAXONOMY_COLUMN
    raise ValueError("No shared usable taxonomy source found: need Lineage or Taxonomy in both profiles")


def normalize_abundance(matrix: pd.DataFrame, method: str) -> pd.DataFrame:
    """Normalize columns to relative abundance when requested."""
    if method == "none":
        return matrix.copy()
    if method != "relative":
        raise ValueError("abundance.normalize must be 'none' or 'relative'")
    out = matrix.copy().astype(float)
    totals = out.sum(axis=0)
    nonzero = totals > 0
    out.loc[:, nonzero] = out.loc[:, nonzero].div(totals[nonzero], axis=1)
    return out


def aggregate_by_rank(
    profile: Profile,
    level: str,
    source: str,
    normalize: str = "none",
    exclude_unclassified: bool = False,
) -> pd.DataFrame:
    """Aggregate feature abundances to one taxonomic rank."""
    level = normalize_level(level)
    if source == LINEAGE_COLUMN:
        taxa = profile.taxonomy_table[source].map(lambda value: parse_lineage(value, level))
    elif source == TAXONOMY_COLUMN:
        taxa = profile.taxonomy_table[source].map(lambda value: parse_taxonomy(value, level))
    else:
        raise ValueError(f"Unsupported taxonomy source: {source}")

    matrix = profile.abundance_matrix.copy()
    matrix["__taxon__"] = taxa.reindex(matrix.index).fillna(UNCLASSIFIED).map(clean_taxon)
    grouped = matrix.groupby("__taxon__", dropna=False).sum()
    grouped.index = grouped.index.map(clean_taxon)
    grouped = grouped.groupby(level=0).sum()
    if exclude_unclassified and UNCLASSIFIED in grouped.index:
        grouped = grouped.drop(index=UNCLASSIFIED)
    grouped.index.name = "Taxon"
    return normalize_abundance(grouped.sort_index(), normalize)
