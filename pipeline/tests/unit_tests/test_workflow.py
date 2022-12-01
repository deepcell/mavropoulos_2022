"""Test workflow.py, used to run scRNA-seq workflow."""

from __future__ import annotations

import argparse
import math

from types import MappingProxyType
from typing import Any

import pandas as pd

from pandas.testing import assert_frame_equal
from pytest import mark

from askcell import Path, parse_path
from askcell.scrna.request import AnalysisRequest
from askcell.utils.types import ImmutableStrMapping
from askcell.workflow import (
    bin_files_directories_from_config,
    convert_pct_to_fraction,
    get_local_paths_settings,
    get_scrna_settings,
    remove_thousand_separator,
    write_data,
)


test_data_dir = Path("tests/data/test_scrna")


@mark.parametrize(
    "command, expected",
    [
        ("analyze", ["input.input_dir=/inputs"]),
        ("count", ["input.reads_dir=/reads", "workflow.sampleid=s1"]),
        ("demux", ["input.run_dir=/220309_VL00213_12_AAAT3NCM5"]),
    ],
)
def test_get_scrna_settings(command: str, expected: list[str]) -> None:
    """Test scrna settings."""
    args = argparse.Namespace()

    args.command = command
    args.request = "request.xlsx"
    args.metadata = "samples.tsv"
    args.local_temp_dir = "/tmp"
    args.output_dir = "/outputs"
    args.input_dir = "/inputs"
    args.reads_dir = "/reads"
    args.run_dir = "/220309_VL00213_12_AAAT3NCM5"
    args.sid = "s1"
    args.keep_temp = None
    args.set_config = None

    observed = get_scrna_settings(
        args,
        AnalysisRequest(test_data_dir / "request.xlsx"),
    )

    expected_base = [
        "input.request=request.xlsx",
        "input.metadata=samples.tsv",
        "workflow.runid=220309_VL00213_12_AAAT3NCM5",
        "workflow.fcid=AAAT3NCM5",
        "workflow.requestid=Test",
        "output.output_dir=/outputs",
        "scratch.scratch_dir=/tmp",
        "params.max_mito_pct_allowed=20",
        "params.max_ribo_pct_allowed=5",
        "params.regress_out_cellcycle=False",
        "params.n_pcs=50",
        "params.n_neighbors=20",
    ]
    assert observed == expected_base + expected


@mark.parametrize(
    "config, local_dir, expected",
    [
        (
            {
                "input": {
                    "file1": parse_path("s3://a/b/c"),
                    "file2": parse_path("/tmp/a/c/"),
                    "file3": parse_path("gs://e/f/g"),
                    "notfile": 100,
                },
                "output": {
                    "file1": parse_path("s3://a/b/c"),
                    "file2": parse_path("/tmp/a/c/"),
                    "file3": parse_path("gs://e/f/g"),
                    "notfile": 100,
                },
            },
            Path("/tmp"),
            [
                "input.file1=/tmp/b/c",
                "input.file3=/tmp/f/g",
                "output.file1=/tmp/b/c",
                "output.file3=/tmp/f/g",
            ],
        )
    ],
)
def test_get_local_paths_settings(
    config: ImmutableStrMapping,
    local_dir: Path,
    expected: list[str],
) -> None:
    """Test get_local_paths_settings."""
    observed = list(get_local_paths_settings(config, local_dir))
    assert sorted(observed) == sorted(expected)


@mark.parametrize(
    "filename, is_parquet, delimiter",
    [
        ("test.tsv", False, "\t"),
        ("test.csv", False, ","),
        ("test.pq", True, None),
    ],
)
def test_write_data(
    tmp_dirname: Path,
    filename: Path,
    is_parquet: bool,
    delimiter: str | None,
) -> None:
    """Test write_data."""
    data = pd.DataFrame(
        {
            "A": [1, 2, 3],
            "B": [4, 5, 6],
        }
    )

    filename = tmp_dirname / filename
    write_data(filename, data)
    assert filename.exists()

    if is_parquet:
        observed = pd.read_parquet(str(filename))
    else:
        observed = pd.read_csv(str(filename), sep=delimiter)

    assert_frame_equal(data, observed)


@mark.parametrize(
    "section, files, dirs, expected_files, expected_dirs",
    [
        (
            MappingProxyType(
                {
                    "v1": Path("/test/test.txt"),
                    "v1_dir": Path("/test"),
                }
            ),
            [],
            [],
            [Path("/test/test.txt")],
            [Path("/test")],
        )
    ],
)
def test_bin_files_directories_from_config(
    section: Any,
    files: list[Path],
    dirs: list[Path],
    expected_files: list[Path],
    expected_dirs: list[Path],
) -> None:
    """Test bin_files_directories_from_config."""
    bin_files_directories_from_config(section, files, dirs)
    assert files == expected_files
    assert dirs == expected_dirs


@mark.parametrize(
    "value, expected",
    [
        ("1.23%", 0.0123),
        ("1.23 %", 0.0123),
        ("1%", 0.01),
        ("1 %", 0.01),
        ("0%", 0),
        ("0 %", 0),
    ],
)
def test_convert_to_fraction(value: Any, expected: float) -> None:
    """Test convert_to_fraction."""
    observed = convert_pct_to_fraction(value)

    assert isinstance(observed, float)
    assert math.isclose(expected, observed)


@mark.parametrize(
    "value, expected",
    [
        ("1,234,567", 1234567),
        ("123", 123),
    ],
)
def test_remove_thousand_separator(value: Any, expected: int) -> None:
    """Test remove_thousand_separator."""
    observed = remove_thousand_separator(value)

    assert isinstance(observed, int)
    assert expected == observed
