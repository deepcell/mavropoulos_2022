"""metrics.py for scRNA-seq analysis workflow."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Iterable, Sequence

import pandas as pd


@dataclass(frozen=True)
class Metric:
    """Container for metric metadata.

    Args:
        name: metric name
        dtype: metric data type
        convert: function for converting specified values
        source: optional original metric name

    Attributes:
        name: metric name
        dtype: metric data type
        convert: function for converting specified values
        source: original metric name

    """

    name: str
    dtype: Any
    convert: Callable[[str], Any] | None = None
    source: str | None = None

    def __post_init__(self) -> None:
        """Finish initialization."""
        if not self.name:
            raise ValueError("name cannot be blank or missing")
        if not self.source:
            object.__setattr__(self, "source", self.name)


class Metrics(frozenset[Metric]):
    """Frozen set of Metrics metadata with some useful operations."""

    def __or__(self, other: Metrics) -> Metrics:  # type: ignore[override]
        """Union with another set of Metrics."""
        return Metrics(super().__or__(other))

    def build_map(self, level: str) -> dict[str, Metric]:
        """Build 1-to-1 mapping of Metric attributes to Metric objects.

        Args:
            metrics: iterable of Metrics
            level: attribute level to index by (source|name)

        Returns:
            Mapping between Metric attributes and Metric objects.

        """
        return build_metrics_map(self, level)


def build_metrics_map(metrics: Iterable[Metric], level: str) -> dict[str, Metric]:
    """Build 1-to-1 mapping of Metric attributes to Metric objects.

    Args:
        metrics: iterable of Metrics
        level: attribute level to index by (source|name)

    Returns:
        Mapping between Metric attributes and Metric objects.

    """
    return {getattr(metric, level): metric for metric in metrics}


def merge(dfs: Sequence[pd.DataFrame], axis: str | int = 0) -> pd.DataFrame:
    """Merge DataFrames (dfs) through concatenation by axis.

    If merging by index (axis=0), indexes are required to be disjoint and column names the same across dfs.
    If merging by columns (axis=1), indexes are required to be the same and column names disjoint across dfs.

    Args:
        dfs:  DataFrame instances to merge
        axis: the axis to concatenate on

    Returns:
        Concatenated DataFrame, or None

    Raises:
        ValueError: if merging by index with non-disjoint indexes
        ValueError: if merging by index with differing column names
        ValueError: if merging by columns with differing indexes
        ValueError: if merging by columns with non-disjoint column names

    """
    if axis in (0, "columns"):
        first_index = set(dfs[0].index)
        if any(first_index != set(df.index) for df in dfs[1:]):
            raise ValueError("indices must be equal")
        all_columns: list[str] = sum((df.columns.tolist() for df in dfs), [])
        if len(set(all_columns)) != len(all_columns):
            raise ValueError("columns must be disjoint")
    elif axis in (1, "index"):
        first_columns = set(dfs[0].columns)
        if any(first_columns != set(df.columns) for df in dfs[1:]):
            raise ValueError("columns must be equal")
        all_indices: list[int] = sum((df.index.tolist() for df in dfs), [])
        if len(set(all_indices)) != len(all_indices):
            raise ValueError("indices must be disjoint")
    else:
        raise ValueError("Invalid axis specified.")

    return pd.concat(dfs, axis=axis, verify_integrity=True)


def update_names(metrics: pd.DataFrame, mapping: dict[str, Metric]) -> pd.DataFrame:
    """Change column names to those designated in mapping.

    Drops column names not specified in mapping.

    Args:
        metrics: data subjected to name updates
        mapping: dictionary mapping column names to Metric definitions

    Returns:
        New DataFrame with new column names.

    """
    name_mapping = {c: mapping[c].name for c in metrics.columns if c in mapping}
    return metrics.rename(columns=name_mapping)  # [name_mapping.keys()]


def update_types(metrics: pd.DataFrame, mapping: dict[str, Metric]) -> pd.DataFrame:
    """Update DataFrame column types to those designated in mapping.

    Requires that all DataFrame columns be specified in mapping.

    Args:
        metrics: data subjected to type updates
        mapping: dictionary mapping column names to Metric definitions

    Returns:
        New DataFrame.

    """
    return metrics.astype({c: mapping[c].dtype for c in metrics.columns})


def update_metrics(metrics: pd.DataFrame, mapping: dict[str, Metric]) -> pd.DataFrame:
    """Update DataFrame column names types to those designated in mapping.

    Args:
        metrics: data subjected to type updates
        mapping: dictionary mapping column names to Metric definitions

    Returns:
        New DataFrame

    """
    return update_types(update_names(metrics, mapping), mapping)
