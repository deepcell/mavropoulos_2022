"""Unit tests for progenity_app_demultiplexer/stat_reader.py."""

from __future__ import annotations

import math

from pathlib import Path
from typing import Any

from pytest import mark

from askcell import PathLike
from askcell.stat_reader import read_conversionstat, read_demuxstat


testdata_dir = Path("tests", "data")


@mark.parametrize(
    "rundir, exp_demuxstats",
    [
        (
            "H5VLNDMXX",
            {
                "LANE_1": {
                    "X1390060": {
                        "CONTROL": None,
                        "INDEX": ["CTGACGTC"],
                        "PF_Q30_BASES_PCT": 83.29,
                        "INDEX_READS_ONE_MISMATCH_PCT": 1.67,
                        "NUM_READS": 4630664,
                        "PF_PCT": 99.49,
                        "YIELD_MBASES": 472,
                        "PF_QUALITY_SCORE_MEAN": 34.21,
                        "CLUSTERS_RAW": 15371,
                        "CLUSTERS_PF": 15292,
                        "CLUSTERS_RAW_PER_LANE_PCT": 0.0005,
                        "INDEX_READS_PERFECT_PCT": 98.33,
                    }
                },
                "LANE_2": {
                    "X1390060": {
                        "CONTROL": None,
                        "INDEX": ["CTGACGTC"],
                        "PF_Q30_BASES_PCT": 82.79,
                        "INDEX_READS_ONE_MISMATCH_PCT": 1.72,
                        "NUM_READS": 4689199,
                        "PF_PCT": 100.00,
                        "YIELD_MBASES": 478,
                        "PF_QUALITY_SCORE_MEAN": 34.14,
                        "CLUSTERS_RAW": 8975,
                        "CLUSTERS_PF": 8975,
                        "CLUSTERS_RAW_PER_LANE_PCT": 0.0003,
                        "INDEX_READS_PERFECT_PCT": 98.28,
                    }
                },
            },
        ),
    ],
)
def test_demuxstat_reader(
    rundir: PathLike,
    exp_demuxstats: dict[str, Any],
) -> None:
    """Test read_demuxstat given an Illumina Stats.json and ConversionStats.xml file."""
    conv_filename = testdata_dir / rundir / "Stats" / "ConversionStats.xml"
    stats_filename = testdata_dir / rundir / "Stats" / "Stats.json"
    demuxstats = read_demuxstat(stats_filename, read_conversionstat(conv_filename))

    for lane in exp_demuxstats:
        for sid in exp_demuxstats[lane]:
            # Verify keys are identical
            assert set(exp_demuxstats[lane][sid]) == set(demuxstats[lane][sid])

            # Verify values are close enough
            for key, expected in exp_demuxstats[lane][sid].items():
                observed = demuxstats[lane][sid][key]
                assert is_close_enough(expected, observed)


def is_close_enough(a: Any, b: Any) -> bool:
    """Compare values a and b and return True if values a and b are close enough or equal, else return False.

    Definition of "close enough" : if objects a and b are instances of int or float,
    difference between a and b should be either <= absolute difference of 0.01
    or <= default relative difference of 1e-09.

    Args:
       a: 1st object.
       b: 2nd object.

    Returns:
       True if values are "close enough", otherwise False.

    """
    a = convert_value(a)
    b = convert_value(b)

    if isinstance(a, float) or isinstance(b, float):
        return math.isclose(a, b, rel_tol=0.001, abs_tol=0.01)
    else:
        return bool(a == b)


def convert_value(value: Any) -> Any:
    """Convert a value into a sensible type.

    Return original value, if value is instance of int, float, list, str, or None.
    else try to convert value to int or float and return converted value.

    Args:
        value: object to be converted.

    Returns:
        original object or converted object.

    """
    if isinstance(value, (int, float)):
        return value
    elif isinstance(value, list):
        return [convert_value(v) for v in value]
    elif isinstance(value, tuple):
        return tuple(convert_value(v) for v in value)
    elif not value:
        return value

    try:
        return int(value)
    except ValueError:
        pass

    try:
        return float(value)
    except ValueError:
        return value
