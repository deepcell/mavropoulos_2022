"""Unit tests for rnaseq.py."""

from __future__ import annotations

import pandas as pd

from pandas.testing import assert_frame_equal
from pytest import raises

from askcell.summarize import get_wgs_cov_bias, get_wgs_scaled_coverage, merge_bias_cov
from askcell.utils.pathutils import Path


test_datadir = Path("tests/data/test_summarize")


def test_get_wgs_cov_bias() -> None:
    """Test get_wgs_base_stats."""
    base_line = pd.read_csv(test_datadir / "normalization_baseline.csv", header=None)
    observed = get_wgs_cov_bias(base_line)
    observed.sort_values(["chrom", "pos_mb"], inplace=True, ignore_index=True)
    expected = pd.read_csv(test_datadir / "expected_cov_bias.csv")
    assert observed.shape == (39, 5)
    assert_frame_equal(observed, expected)


def test_get_wgs_scaled_coverage() -> None:
    """Test get_wgs_scaled_coverage."""
    df = pd.read_csv(test_datadir / "genome_mb_coverage.csv", header=None)
    observed = get_wgs_scaled_coverage(df)
    assert_frame_equal(observed, pd.read_csv(test_datadir / "expected_scaled_cov.csv", index_col=0))


def test_merge_bias_cov() -> None:
    """Test merge_bias_cov."""
    cov_bias = pd.read_csv(test_datadir / "expected_cov_bias.csv", header=0)
    scaled_cov = pd.read_csv(test_datadir / "expected_scaled_cov.csv", index_col=0, header=0)
    observed = merge_bias_cov(cov_bias, scaled_cov)
    assert_frame_equal(observed, pd.read_csv(test_datadir / "expected_merged.csv", index_col=0, header=0))

    with raises(ValueError, match="No data left in the merged dataframe after filtering."):
        observed = merge_bias_cov(cov_bias, scaled_cov, include_chromosome="chr1")
