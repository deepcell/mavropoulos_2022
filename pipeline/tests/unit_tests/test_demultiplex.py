"""Unit tests for askcell.demultplex.py."""

from pathlib import Path

import pandas as pd

from interop.py_interop_run import xml_file_not_found_exception
from pandas.testing import assert_frame_equal
from pytest import mark, raises

from askcell.demultiplex import get_runsummary


testdata_dir = Path("tests/data/test_demultiplex")


@mark.output_validated
def test_get_runsummary() -> None:
    """Test get_runsummary."""
    with raises(xml_file_not_found_exception, match="cannot open file nonexistent_dir/RunInfo.xml"):
        get_runsummary("nonexistent_dir")
    data = testdata_dir / "220308_MN01265_0132_A000H3M3CL"
    assert_frame_equal(
        get_runsummary(data),
        pd.DataFrame(
            {
                "total_yield": [1.433197],
                "bases_gt_q30": [0.966619],
                "clusters_frac_pf": [0.931771],
                "cluster_density": [69.766667],
            }
        ),
    )
