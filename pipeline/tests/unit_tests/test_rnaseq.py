"""Unit tests for rnaseq.py."""

from __future__ import annotations

import numpy as np
import pandas as pd

from pandas.testing import assert_frame_equal

from askcell.rna.rnaseq import (
    DTYPES_ALIGNMENT_STATS,
    get_aln_stats,
    get_gene_count_stats,
    get_rRNA_stats,
    get_trim_stats,
    summarize_stats,
)
from askcell.utils.pathutils import Path


testdata_dir = Path("tests", "data", "test_rnaseq")


def test_get_rRNA_stats() -> None:
    """Test get rRNA stats."""
    filename = testdata_dir / "sample1.rRNA.alignments.json"
    sampleid = "sample1"

    observed = get_rRNA_stats(filename, sampleid)
    expected = pd.DataFrame(
        {
            "secondary": pd.Series([0], dtype=np.uint64),
            "mapped": pd.Series([12219], dtype=np.uint64),
        },
    ).set_index(pd.Index([sampleid]))

    assert_frame_equal(observed, expected)


def test_get_aln_stats() -> None:
    """Test get alignment stats."""
    filename = testdata_dir / "sample1.log.final.out"
    sampleid = "sample1"

    observed = get_aln_stats(filename, sampleid)
    expected = pd.DataFrame(
        {
            "Uniquely mapped reads number": pd.Series([514994], dtype=np.uint64),
            "Number of input reads": pd.Series([581476], dtype=np.uint64),
            "Number of reads mapped to multiple loci": pd.Series([61482], dtype=np.uint64),
            "Number of reads mapped to too many loci": pd.Series([498], dtype=np.uint64),
            "Number of reads unmapped: too short": pd.Series([4262], dtype=np.uint64),
            "Number of reads unmapped: other": pd.Series([240], dtype=np.uint64),
            "Mismatch rate per base, %": pd.Series(["0.31%"], dtype=str),
        },
    ).set_index(pd.Index([sampleid]))
    expected.columns.name = "Name"

    assert_frame_equal(observed[DTYPES_ALIGNMENT_STATS], expected)


def test_get_gene_count_stats() -> None:
    """Test get gene count stats."""
    filename = testdata_dir / "sample1.ReadsPerGene.out.tab"
    sampleid = "sample1"

    observed = get_gene_count_stats(filename, sampleid)
    expected = pd.DataFrame(
        {
            "N_FWD": pd.Series([416], dtype=np.uint64),
            "N_REV": pd.Series([118], dtype=np.uint64),
        },
    ).set_index(pd.Index([sampleid]))

    assert_frame_equal(observed, expected)


def test_get_trim_stats() -> None:
    """Test get trim stats."""
    filename = testdata_dir / "sample1.trimmed.info.txt"
    sampleid = "sampleid"

    observed = get_trim_stats(filename, sampleid)
    expected = pd.DataFrame(
        {
            "MeanFragLen": pd.Series([67], dtype=np.float64),
            "MeanAdapterLen": pd.Series([46], dtype=np.float64),
            "PCTAdapter": pd.Series([70], dtype=np.float64),
            "PCT_A": pd.Series([96.5834], dtype=np.float64),
            "PCT_G": pd.Series([0.93168], dtype=np.float64),
        },
    ).set_index(pd.Index([sampleid]))

    assert_frame_equal(observed, expected, check_exact=False)


def test_summarize_stats() -> None:
    """Summarize stats."""
    sampleid = "sampleid"
    observed = summarize_stats(
        pd.concat(
            [
                get_rRNA_stats(testdata_dir / "sample1.rRNA.alignments.json", sampleid),
                get_aln_stats(testdata_dir / "sample1.log.final.out", sampleid),
                get_gene_count_stats(testdata_dir / "sample1.ReadsPerGene.out.tab", sampleid),
                get_trim_stats(testdata_dir / "sample1.trimmed.info.txt", sampleid),
            ],
            axis=1,
        ),
    )

    expected = pd.DataFrame(
        {
            "NumFrags": pd.Series([581476], dtype=np.uint64),
            "MismatchRate": pd.Series(["0.31%"], dtype=str),
            "MeanFragLen": pd.Series([67], dtype=np.float64),
            "MeanAdapterLen": pd.Series([46], dtype=np.float64),
            "PCTAdapter": pd.Series([70], dtype=np.float64),
            "PCT_A": pd.Series([96.5834], dtype=np.float64),
            "PCT_G": pd.Series([0.93168], dtype=np.float64),
            "PCT_UniqMapped": pd.Series([88.566682], dtype=np.float64),
            "PCT_MultiMapped": pd.Series([10.573437], dtype=np.float64),
            "PCT_Unmapped": pd.Series([0.8598807], dtype=np.float64),
            "PCT_FwdStranded": pd.Series([0.0715421], dtype=np.float64),
            "PCT_RevStranded": pd.Series([0.0202932], dtype=np.float64),
            "PCT_rRNA": pd.Series([2.1013765], dtype=np.float64),
        },
    ).set_index(pd.Index([sampleid], name="SampleName"))

    assert_frame_equal(observed, expected, check_exact=False)
