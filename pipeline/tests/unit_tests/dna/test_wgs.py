"""Unit tests for wgs.py."""

from __future__ import annotations

import unittest
import unittest.mock as mock

from collections import defaultdict

import pandas as pd
import pysam

from pandas.testing import assert_frame_equal
from pytest import mark

from askcell.dna.wgs import (
    count_gc_content,
    get_galignment_tally,
    get_gc_tally,
    load_config,
    parse_bam,
    parse_samtools_stats,
)
from askcell.utils.fileutils import load_toml
from askcell.utils.interpolate import interpolate
from askcell.utils.merge import merge
from askcell.utils.pathutils import Path
from askcell.utils.types import MutableStrMapping


testdata_dir = Path("tests", "data", "test_wgs", "220407_MN01265_0135_A000H3LYGL")
alignment = testdata_dir / "DLS-62321-Histiocytes_S9_L001.bam"


@unittest.skip("need access to bucket to access bam files")
@mark.output_validated
def test_parse_bam() -> None:
    """Test parse_bam."""
    observed_alignment, observed_gc, observed_duplicate_stats = parse_bam(
        "DLS-62321-Histiocytes_S9",
        alignment,
    )
    expected_alignment = pd.read_csv(
        testdata_dir / "expected_mb_coverage.csv",
        sep=",",
        names=["cov", "ref_name", "ref_chunk", "sample_id"],
        header=0,
    )
    expected_gc = pd.read_csv(
        testdata_dir / "gc_histogram.csv",
        sep=",",
        names=["count", "range", "sample_id"],
    )
    assert_frame_equal(observed_alignment, expected_alignment)
    assert_frame_equal(observed_gc, expected_gc)
    assert_frame_equal(
        observed_duplicate_stats,
        pd.DataFrame({"smample_id": ["DLS-62321-Histiocytes_S9"], "total": [1722010], "unique": [1487802]}),
    )


@mark.output_validated
def test_count_gc_content() -> None:
    """Test count_gc_content."""
    sequence = (
        "TAATGATACGGCGACCACCGAGATCTACACCGAC"
        "TCTCACACTCTTTCCCTACACGACGCTCTTC"
        "CGATCTGGATTAAACCAAACCCAACTACGCA"
        "AAATCTTAGCATACTCCTCAATTACCCACATAGGATGAATAACAGCAGTTCTGAC"
    )
    observed = count_gc_content(sequence)
    expected = 0.456
    assert observed == expected


@unittest.skip("need access to bucket to access bam files")
@mark.output_validated
def test_get_gc_tally() -> None:
    """Test get_gc_tally."""
    with pysam.AlignmentFile(str(alignment), "r", threads=4) as f:
        for x in f.fetch(multiple_iterators=True):
            observed = get_gc_tally(x)
    assert observed[0.0] == 8
    assert observed[0.006] == 2


@unittest.skip("need access to bucket to access bam files")
@mark.output_validated
def test_get_galignment_tally() -> None:
    """Test get_galignment_tally."""
    prev = None
    unique = 0
    observed: defaultdict[tuple[str | None, int], int] = defaultdict(int)
    with pysam.AlignmentFile(str(alignment), "r", threads=4) as f:
        for x in f.fetch(multiple_iterators=False):
            observed, prev, unique = get_galignment_tally(x, existing_tally=observed, prev=prev, unique=unique)
    assert observed[("chr1", 0)] == 64
    assert len(list(observed.keys())) == 3195
    assert unique == 1487802


@mark.output_validated
def test_get_mismatch_stats() -> None:
    """Test get_mismatch_stats."""
    samtools_stats = testdata_dir / "samtools_stats.txt"
    error_rate, gc_stats = parse_samtools_stats(
        "x",
        samtools_stats,
    )
    assert_frame_equal(
        error_rate,
        pd.DataFrame(
            [
                ["x", 4.241944e-03],
            ]
        ),
    )
    assert_frame_equal(
        gc_stats.iloc[:1],
        pd.DataFrame(
            [
                ["x", 0.0, 0.002, 0.000, 0.000, 0.000, 0.000, 0.000],
            ]
        ),
    )


@mark.parametrize(
    "config, set_config, extra",
    [
        (
            Path("tests/data/test_wgs/load_config/user_config.toml"),
            ["readcounter.options.quality=30"],
            {
                "merge_bias_cv": {"options": {"bias_cv_threshold": 0.3}},
            },
        ),
    ],
)
def test_load_config(
    config: Path | None,
    set_config: list[str] | None,
    extra: MutableStrMapping,
    tmp_path: Path,
) -> None:
    """Test load_config.

    Args:
        config: Optionally override default package configs built from schema using data files
        set_config: Nested key assignment expressions of the form key[.key]*=value used to override data values
        extra: extra config values to update
        tmp_path: pytest temp directory fixture

    """
    test_data = Path("tests/data/test_wgs/load_config")

    params = {
        "params": {
            "input_dir": "mock_input_dir",
            "output_dir": tmp_path,
        },
    }

    with mock.patch("subprocess.check_call", new=tmp_path):
        load_config(
            config=config,
            input_dir="mock_input_dir",
            runid="mock_run_id",
            output_dir=tmp_path,
            set_config=set_config,
        )

    expected = interpolate(merge([load_toml(test_data / "expected_config.toml"), params, extra]))
    assert expected == load_toml(tmp_path / "config.toml")
