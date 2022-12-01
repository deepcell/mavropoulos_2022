"""Demultiplex analysis task for scRNA-seq.

Demultiplex the reads for a run.
This is currently focused around 10x platform data processing
using cellranger mkfastq.

These are the steps in this task.

    * Create samplesheet for cellranger mkfastq and execute
    * Collect and write metrics
    * Place the reads in a single dir to sample-level dirs

"""

from __future__ import annotations

import csv

from typing import Any

import numpy as np
import pandas as pd

from .. import PathLike, only_one
from ..stat_reader import read_conversionstat, read_demuxstat, summarize_demuxstat
from ..utils.config import PackageConfig
from ..utils.rowcontainer import RowContainer
from ..utils.shellutils import run_cmd
from ..utils.types import ImmutableStrMapping
from ..workflow import write_data

# from ..pathutils import set_path_readonly
from ..wrappers import make_command_cellranger_mkfastq
from .metrics import Metric, Metrics, update_names  # check_metrics
from .request import AnalysisRequest, ChromiumIndex


# TODO: Update the list
DEMUX_METRICS = Metrics(
    [
        Metric("SC_DEMUX_TOTAL_READS", int, source="SUMMARYSTAT.READS_TOTAL"),
        Metric("SC_DEMUX_FRACTION_Q30", float, source="SUMMARYSTAT.CLUSTERS_PF_TOTAL,SUMMARYSTAT.CLUSTERS_RAW_TOTAL"),
        Metric("SC_DEMUX_MEDIAN_READS", float, source="SUMMARYSTAT.READS_PER_SAMPLE_MEDIAN"),
    ]
)

DEMUX_MAPPINGS = DEMUX_METRICS.build_map("source")


def run_cellranger_mkfastq(
    inputs: ImmutableStrMapping,
    outputs: ImmutableStrMapping,
    scratch: ImmutableStrMapping,
    params: ImmutableStrMapping,
    request: AnalysisRequest,
    sample_configs: dict[str, tuple[PackageConfig, PackageConfig]],
    **kwargs: Any,
) -> None:
    """Run cellranger mkfastq as a step."""
    # TODO: add option to create samplesheet from old seq run samplesheet
    # create_samplesheet_from_old_samplesheet(
    #     ss_old=input[""],
    #     chromium_index=inputs["index"],
    #     ss=scratch["ss"],
    # )
    create_samplesheet(
        request=request,
        chromium_index=inputs["index"],
        ss=scratch["ss"],
    )
    run_cmd(make_command_cellranger_mkfastq(**params))


def run_sequencing_qc(
    inputs: ImmutableStrMapping,
    outputs: ImmutableStrMapping,
    scratch: ImmutableStrMapping,
    params: ImmutableStrMapping,
    request: AnalysisRequest,
    sample_configs: dict[str, tuple[PackageConfig, PackageConfig]],
    **kwargs: Any,
) -> None:
    """Run sequencing metrics QC as a step."""
    metrics = collect_metrics(
        inputs["conversion_stats"],
        inputs["bcl2fastq_stats"],
    )

    qc_metrics = update_names(check_metrics(params, metrics), DEMUX_MAPPINGS)

    write_data(
        scratch["metrics_result"],
        qc_metrics,
        use_index=True,
    )
    write_data(
        scratch["metrics_result_pq"],
        qc_metrics,
        use_index=True,
    )


def create_samplesheet(
    request: AnalysisRequest,
    chromium_index: PathLike,
    ss: PathLike,
    index_version: str = "A",
) -> None:
    """Create samplesheet to use during cellranger mkfastq.

    Workflow A = Illumina Forward Strand Sequencing Workflow
    Workflow B = Illumina Reverse Complement Sequencing Workflow

    Args:
        request: sequencing run and sample information
        chromium_index: index sequence to index identifier mapping file
        ss: 10x samplesheet
        index_version: A for fwd strand sequencing and B for rev. complement

    """
    with open(chromium_index) as csvfile:
        indices = RowContainer(
            ChromiumIndex,
            (ChromiumIndex(row) for row in csv.DictReader(csvfile, delimiter=",")),
        )
    index_maps = indices.build_index(*["index1", f"index2_{index_version}"])

    with open(ss, "w") as fout:
        writer = csv.writer(fout, delimiter=",")
        writer.writerow(["Lane", "Sample", "Index"])
        for s in request.samples:
            if s.index1 in index_maps and s.index2 in index_maps[s.index1]:
                writer.writerow(["1", s.sid, only_one(index_maps[s.index1][s.index2]).name])
            else:
                raise ValueError(f"Index {s.index1}, {s.index2}: not match found.")


def create_samplesheet_from_old_samplesheet(
    ss_old: PathLike,
    chromium_index: PathLike,
    ss: PathLike,
    index_version: str = "A",
) -> None:
    """Load sequencing run samplesheet to use during cellranger mkfastq.

    Use the samplesheet created during the sequencing run to
    create a new samplesheet for 10x recommended for cellranger mkfastq.

    Workflow A = Illumina Forward Strand Sequencing Workflow
    Workflow B = Illumina Reverse Complement Sequencing Workflow

    Args:
        ss_old: samplesheet generated during the sequencing run
        chromium_index: index sequence to index identifier mapping file
        ss: 10x samplesheet
        index_version: A for fwd strand sequencing and B for rev. complement

    """
    is_index_section = False

    with open(chromium_index) as csvfile:
        indices = RowContainer(
            ChromiumIndex,
            (ChromiumIndex(row) for row in csv.DictReader(csvfile, delimiter=",")),
        )
    index_maps = indices.build_index(*["index1", f"index2_{index_version}"])

    with open(ss_old, newline="") as fin, open(ss, "w") as fout:
        reader = csv.reader(fin, delimiter=",")
        writer = csv.writer(fout, delimiter=",")
        writer.writerow(["Lane", "Sample", "Index"])

        for row in reader:
            if row[0] == "[BCLConvert_Data]":
                is_index_section = True
                header = next(reader)
            elif is_index_section and len(row) == 0:
                is_index_section = False

            if is_index_section:
                s = dict(zip(header, row))
                idx = only_one(index_maps[s["Index"]][s["Index2"]])
                writer.writerow(["1", s["Sample_ID"], idx.name])


def collect_metrics(
    conversion_stats: PathLike,
    stats: PathLike,
) -> pd.DataFrame:
    """Collect sequencing run metrics and demultiplexed stats.

    Args:
        conversion_stats: conversion stats in xml
        stats: stats output from demultiplexer in JSON

    Returns:
        metrics data on sequencing run performance

    """
    convstat = read_conversionstat(conversion_stats)
    demuxstat = read_demuxstat(stats, convstat)

    stat = {
        "DEMUXSTAT": demuxstat,
        "SUMMARYSTAT": summarize_demuxstat(demuxstat),
    }

    return pd.json_normalize(stat)


def check_metrics(
    conditions: ImmutableStrMapping,
    data: pd.DataFrame,
) -> pd.DataFrame:
    """Check the metrics against the acceptance criteria.

    Args:
        conditions: metrics and passing criteria information
        data: metrics data with metrics name and value

    Returns:
        metrics data with metrics name, value, acceptance criteria, pass/fail flags

    """
    # call check_qc on all metrics to be checked

    # FIXME: automate the comparison using config
    # FIXME: populate more metrics checks
    metric_name = f'{DEMUX_MAPPINGS["SUMMARYSTAT.READS_TOTAL"].name}_QC'
    data[metric_name] = np.where(
        data["SUMMARYSTAT.READS_TOTAL"] > conditions["min_total_reads"],
        "Pass",
        "Fail",
    )

    metric_name = f'{DEMUX_MAPPINGS["SUMMARYSTAT.CLUSTERS_PF_TOTAL,SUMMARYSTAT.CLUSTERS_RAW_TOTAL"].name}_QC'
    data["SUMMARYSTAT.CLUSTERS_PF_TOTAL,SUMMARYSTAT.CLUSTERS_RAW_TOTAL"] = (
        data["SUMMARYSTAT.CLUSTERS_PF_TOTAL"] / data["SUMMARYSTAT.CLUSTERS_RAW_TOTAL"]
    )
    data[metric_name] = np.where(
        data["SUMMARYSTAT.CLUSTERS_PF_TOTAL,SUMMARYSTAT.CLUSTERS_RAW_TOTAL"] > conditions["min_fraction_q30"],
        "Pass",
        "Fail",
    )

    metric_name = f'{DEMUX_MAPPINGS["SUMMARYSTAT.READS_PER_SAMPLE_MEDIAN"].name}_QC'
    data[metric_name] = np.where(
        data["SUMMARYSTAT.READS_PER_SAMPLE_MEDIAN"] > conditions["min_median_reads"],
        "Pass",
        "Fail",
    )

    return data
