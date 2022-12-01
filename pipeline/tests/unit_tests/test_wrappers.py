"""Unit tests for wrappers.py."""

from __future__ import annotations

from askcell.utils.pathutils import parse_path
from askcell.utils.shellutils import Command
from askcell.wrappers import make_command_read_counter, make_command_run_ichorcna


def test_make_command_read_counter() -> None:
    """Test make_command_read_counter."""
    observed = make_command_read_counter(
        aligned_read="tests/data/mock.bam",
        output="tests/data/mock.wig",
    )

    chromosome = (
        "chr1,chr2,chr3,chr4,chr5,chr6,chr7,"
        "chr8,chr9,chr10,chr11,chr12,chr13,"
        "chr14,chr15,chr16,chr17,chr18,chr19,"
        "chr20,chr21,chr22,chrX,chrY"
    )

    expected = Command(
        "readCounter",
        "--window",
        1000000,
        "--quality",
        20,
        "--chromosome",
        chromosome,
        "tests/data/mock.bam",
    ).redirect_output(stdout="tests/data/mock.wig")

    assert observed == expected


def test_make_command_run_ichorcna() -> None:
    """Test make_command_run_ichorcna."""
    sampleid = "mock_id"
    bigwig = "tests/data/mock.wig"
    output = "tests/data/output"

    observed = make_command_run_ichorcna(
        sampleid=sampleid,
        bigwig=bigwig,
        output=output,
        estimate_normal=True,
    )

    expected = Command(
        "Rscript",
        parse_path("askcell/dna/runIchorCNA.R"),
        "--id",
        sampleid,
        "--WIG",
        bigwig,
        "--ploidy",
        "c(1.75,2,2.25)",
        "--normal",
        "c(0.01,0.1,0.25,0.5,0.75,0.9,0.99)",
        "--gcWig",
        parse_path("askcell/data/gc_hg38_1000kb.wig"),
        "--mapWig",
        parse_path("askcell/data/map_hg38_1000kb.wig"),
        "--centromere",
        parse_path("askcell/data/GRCh38.GCA_000001405.2_centromere_acen.txt"),
        "--normalPanel",
        parse_path("askcell/data/cnv.pon_median.rds"),
        "--chrs",
        'c(1:22,"X")',
        "--chrTrain",
        "c(1:22)",
        "--txnE",
        0.9999,
        "--txnStrength",
        10000,
        "--outDir",
        parse_path(output),
        "--estimateNormal",
        True,
    )

    assert observed == expected
