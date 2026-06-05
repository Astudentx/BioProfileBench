"""Profile data model and profile construction helpers."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd


@dataclass
class Profile:
    """Unified abundance profile used by all benchmark components."""

    profile_id: str
    dataset: str
    method: str
    gene_kind: str
    parameter_tag: str
    abundance_matrix: pd.DataFrame
    taxonomy_table: pd.DataFrame
    sample_metadata: dict[str, Any] = field(default_factory=dict)
    feature_metadata: pd.DataFrame | None = None
    source_path: Path | None = None
    extra_metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def sample_cols(self) -> list[str]:
        return list(self.abundance_matrix.columns)

    @property
    def feature_ids(self) -> list[str]:
        return list(self.abundance_matrix.index)

    @property
    def metadata(self) -> dict[str, str]:
        metadata = {
            "profile_id": self.profile_id,
            "dataset": self.dataset,
            "method": self.method,
            "gene_kind": self.gene_kind,
            "parameter_tag": self.parameter_tag,
        }
        metadata.update(self.extra_metadata)
        return metadata
