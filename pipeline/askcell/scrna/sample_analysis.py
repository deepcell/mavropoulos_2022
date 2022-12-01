"""Sample-level analysis task for scRNA-seq.

Align and count the reads in a single sample.
This is currently very focused around 10x platform data.

Using the raw reads in fastqs, it utilizes the cellranger tool to
generate alignment and count data.

These are the steps in this task.

    * Create command for cellranger count and execute
    * Collect and write metrics

"""

from __future__ import annotations

import shutil

from typing import Any, Callable

import numpy as np
import pandas as pd

from .. import PathLike, only_one, parse_path
from ..utils.config import PackageConfig
from ..utils.shellutils import run_cmd
from ..utils.types import ImmutableStrMapping
from ..workflow import (  # check_qc,
    convert_pct_to_fraction,
    remove_thousand_separator,
    write_data,
)
from ..wrappers import make_command_cellranger_count
from .metrics import Metric, Metrics, update_names
from .request import AnalysisRequest, Sample


CELLRANGER_METRICS = Metrics(
    [
        Metric("SC_COUNT_CELLS", int, remove_thousand_separator, source="Estimated Number of Cells"),
        Metric("SC_MEAN_READS_PER_CELL", int, remove_thousand_separator, source="Mean Reads per Cell"),
        Metric("SC_MEDIAN_GENES_PER_CELL", int, remove_thousand_separator, source="Median Genes per Cell"),
        Metric("SC_COUNT_READS", int, remove_thousand_separator, source="Number of Reads"),
        Metric("SC_FRACTION_BARCODES_IN_WHITELIST", float, convert_pct_to_fraction, source="Valid Barcodes"),
        Metric("SC_FRACTION_READS_IN_UMI", float, convert_pct_to_fraction, source="Sequencing Saturation"),
        Metric("SC_FRACTION_Q30BASES_IN_BARCODE", float, convert_pct_to_fraction, source="Q30 Bases in Barcode"),
        Metric("SC_FRACTION_Q30BASES_IN_READ", float, convert_pct_to_fraction, source="Q30 Bases in RNA Read"),
        Metric("SC_FRACTION_Q30BASES_IN_UMI", float, convert_pct_to_fraction, source="Q30 Bases in UMI"),
        Metric("SC_FRACTION_READS_MAPPED", float, convert_pct_to_fraction, source="Reads Mapped to Genome"),
        Metric(
            "SC_FRACTION_READS_MAPPED_UNIQUE",
            float,
            convert_pct_to_fraction,
            source="Reads Mapped Confidently to Genome",
        ),
        Metric(
            "SC_FRACTION_READS_MAPPED_INTERGENIC",
            float,
            convert_pct_to_fraction,
            source="Reads Mapped Confidently to Intergenic Regions",
        ),
        Metric(
            "SC_FRACTION_READS_MAPPED_INTRONIC",
            float,
            convert_pct_to_fraction,
            source="Reads Mapped Confidently to Intronic Regions",
        ),
        Metric(
            "SC_FRACTION_READS_MAPPED_EXON",
            float,
            convert_pct_to_fraction,
            source="Reads Mapped Confidently to Exonic Regions",
        ),
        Metric(
            "SC_FRACTION_READS_MAPPED_TRANSCRIPTOME",
            float,
            convert_pct_to_fraction,
            source="Reads Mapped Confidently to Transcriptome",
        ),
        Metric(
            "SC_FRACTION_READS_MAPPED_ANTISENSE",
            float,
            convert_pct_to_fraction,
            source="Reads Mapped Antisense to Gene",
        ),
        Metric("SC_FRACTION_READS_IN_CELLS", float, convert_pct_to_fraction, source="Fraction Reads in Cells"),
        Metric("SC_COUNT_GENES", int, remove_thousand_separator, source="Total Genes Detected"),
        Metric("SC_MEDIAN_UMIS_PER_CELL", int, remove_thousand_separator, source="Median UMI Counts per Cell"),
        Metric("SC_COUNT_CELLS_QC", str),
        Metric("SC_MEAN_READS_PER_CELL_QC", str),
        Metric("SC_MEDIAN_GENES_PER_CELL_QC", str),
        Metric("SC_MEDIAN_UMIS_PER_CELL_QC", str),
        Metric("SC_COUNT_GENES_QC", str),
    ]
)

CELLRANGER_MAPPINGS = CELLRANGER_METRICS.build_map("source")


def run_cellranger_count(
    inputs: ImmutableStrMapping,
    outputs: ImmutableStrMapping,
    scratch: ImmutableStrMapping,
    params: ImmutableStrMapping,
    request: AnalysisRequest,
    sample_configs: dict[str, tuple[PackageConfig, PackageConfig]],
    **kwargs: Any,
) -> None:
    """Run cellranger count as a step."""
    shutil.rmtree(scratch["10x_dir"])
    run_cmd(make_command_cellranger_count(**params))

    sid = params["id"]
    write_data(
        scratch["metadata2"],
        make_sample_metadata(
            sid,
            only_one(request.samples.build_index("sid")[sid]),
            outputs["raw_matrix_h5"],
        ),
    )


def run_qc(
    inputs: ImmutableStrMapping,
    outputs: ImmutableStrMapping,
    scratch: ImmutableStrMapping,
    params: ImmutableStrMapping,
    request: AnalysisRequest,
    sample_configs: dict[str, tuple[PackageConfig, PackageConfig]],
    **kwargs: Any,
) -> None:
    """Run metrics QC as a step."""
    metrics = collect_metrics(
        inputs["sid"],
        filename=inputs["metrics"],
        metrics=CELLRANGER_METRICS,
    )

    qc_metrics = update_names(check_metrics(params, metrics), CELLRANGER_MAPPINGS)

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


def collect_metrics(
    sid: str,
    filename: PathLike,
    metrics: Metrics | None = None,
) -> pd.DataFrame:
    """Collect align/count metrics.

    FIXME: Separate out collect and metrics check into 2 functions.

    Args:
        sid: name of the sample
        filename: cellranger metrics summary
        metrics: collection of metrics metadata

    Returns:
        metrics data on sample alignment and count performance

    """
    converter_maps: dict[str, Callable[[str], Any]] = {}

    if metrics is not None:
        converter_maps = {
            metric.source: metric.convert
            for metric in metrics
            if metric.source is not None and metric.convert is not None
        }

    data = pd.read_csv(
        parse_path(filename),
        converters=converter_maps,
        index_col=False,
    )
    return data.set_index([pd.Index([sid], name="SID")])


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
    metric_name = f'{CELLRANGER_MAPPINGS["Estimated Number of Cells"].name}_QC'
    data[metric_name] = np.where(
        data["Estimated Number of Cells"] > conditions["min_number_cells"],
        "Pass",
        "Fail",
    )

    metric_name = f'{CELLRANGER_MAPPINGS["Mean Reads per Cell"].name}_QC'
    data[metric_name] = np.where(
        data["Mean Reads per Cell"] > conditions["min_reads_per_cell_rate"],
        "Pass",
        "Fail",
    )

    metric_name = f'{CELLRANGER_MAPPINGS["Median Genes per Cell"].name}_QC'
    data[metric_name] = np.where(
        data["Median Genes per Cell"] > conditions["min_median_genes_per_cell"],
        "Pass",
        "Fail",
    )

    metric_name = f'{CELLRANGER_MAPPINGS["Median UMI Counts per Cell"].name}_QC'
    data[metric_name] = np.where(
        data["Median UMI Counts per Cell"] > conditions["min_umi_count_per_cell"],
        "Pass",
        "Fail",
    )

    metric_name = f'{CELLRANGER_MAPPINGS["Total Genes Detected"].name}_QC'
    data[metric_name] = np.where(
        data["Total Genes Detected"] > conditions["min_total_genes_detected"],
        "Pass",
        "Fail",
    )
    return data


def make_sample_metadata(
    sid: str,
    sample: Sample,
    count: PathLike,
) -> pd.DataFrame:
    """Make new metadata with count path."""
    return pd.DataFrame(
        {
            "SampleName": [sample.name],
            "RunID": [sample.runid],
            "Instrument": [sample.instrument.name if sample.instrument is not None else "."],
            "Index1": [sample.index1],
            "Index2": [sample.index2],
            "Count": [str(count)],
            "SampleType": [sample.label if sample.label is not None else "."],
        }
    )
