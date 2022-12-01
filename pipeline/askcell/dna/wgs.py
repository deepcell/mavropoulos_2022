"""wga.py for wga analysis. TODO: add later."""

from __future__ import annotations

import argparse
import csv
import re
import tempfile

from collections import defaultdict
from io import StringIO
from typing import Any

import pandas as pd
import pysam

from askcell.utils.fileutils import save_toml

from .. import only_one
from ..static_data import GENOME_INDEX, WGS_BASELINE
from ..summarize import get_wgs_cov_bias, get_wgs_scaled_coverage, merge_bias_cov
from ..utils.config import (
    PackageResource,
    load_config_data,
    load_package_data,
    load_package_schema,
)
from ..utils.pathutils import Path, PathLike, SomePathLikes, parse_path
from ..utils.shellutils import Command, run_cmd
from ..utils.types import MutableStrMapping, unparse_data
from ..wrappers import (
    make_command_bwa_mem,
    make_command_picard_markduplicate,
    make_command_read_counter,
    make_command_run_ichorcna,
    make_command_samtools_index,
    make_command_samtools_stats,
    make_command_samtools_view,
)


def parse_bam(
    sample_id: str,
    alignment: PathLike,
    min_mapq: int = 50,
    alignment_bin_size: int = 1000000,
    gc_bin_size: int = 1000,
    threads: int = 4,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Parse a single bam file to get gc content and good alignment stats.

    Args:
        sample_id: sample ID.
        alignment: alignment .bam file.
        min_mapq: threshold value for good coverage.
        alignment_bin_size: default bin size for alignment.
        gc_bin_size: default bin size for gc histogram.
        threads: number of threads to use for concurrent processing.

    Returns:
        Dataframe containing "cov", "chrom", "pos_mb", and "sample_name".
        Dataframe containing "bin index", "count".

    """
    galignment_tally: defaultdict[tuple[str | None, int], int] = defaultdict(int)
    gc_tally: defaultdict[float, int] = defaultdict(int)
    prev: tuple[str | None, int] | None = None
    measurements: int = 0
    unique: int = 0

    with pysam.AlignmentFile(str(parse_path(alignment)), "r", threads=threads) as f:
        for segment in f.fetch(multiple_iterators=True):
            if segment.mapping_quality > min_mapq:
                measurements += 1
                galignment_tally, prev, unique = get_galignment_tally(
                    segment,
                    existing_tally=galignment_tally,
                    prev=prev,
                    unique=unique,
                    bin_size=alignment_bin_size,
                )

                gc_tally = get_gc_tally(
                    segment,
                    existing_tally=gc_tally,
                    bin_size=gc_bin_size,
                )

    galignment_tally_df = pd.DataFrame(
        data={"cov": list(galignment_tally.values())},
        index=pd.MultiIndex.from_tuples(tuples=galignment_tally.keys(), names=["ref_name", "ref_chunk"]),
    )
    galignment_tally_df["sample_id"] = sample_id
    galignment_tally_df.reset_index(inplace=True)
    galignment_tally_df = galignment_tally_df[["cov", "ref_name", "ref_chunk", "sample_id"]]

    gc_tally_df = pd.DataFrame.from_dict(
        {key: [gc_tally[key]] for key in sorted(gc_tally.keys())},
        orient="index",
    )
    gc_tally_df["sample_id"] = sample_id
    gc_tally_df.reset_index(inplace=True)
    gc_tally_df.set_axis(["range", "count", "sample_id"], axis=1, inplace=True)
    gc_tally_df = gc_tally_df[["count", "range", "sample_id"]]

    duplicate_stats = pd.DataFrame({"smample_id": [sample_id], "total": [measurements], "unique": [unique]})

    return galignment_tally_df, gc_tally_df, duplicate_stats


def count_gc_content(
    sequence: str,
    bin_size: int = 1000,
) -> float:
    """Get G / C content in a given sequence.

    Args:
        sequence: sequence to be counted.
        bin_size: default bin size for gc histogram.

    Returns:
        fraction of g or c content.

    """
    return int(len(re.findall(r"[GCgc]", sequence)) / len(sequence) * bin_size) / bin_size


def get_gc_tally(
    segment: pysam.libcalignedsegment.AlignedSegment,
    existing_tally: defaultdict[float, int] | None = None,
    bin_size: int = 1000,
) -> defaultdict[float, int]:
    """Get tally for GC per file count.

    Args:
        segment: pysam alignment segment
        existing_tally: A tally of counts in each G / C fraction bin
        bin_size: default bin size for gc histogram.
        min_mapq: threshold value for good coverage.

    Returns:
        An updated tally of counts in each G / C fraction bin

    """
    if existing_tally is None:
        existing_tally = defaultdict(int)

    existing_tally[
        count_gc_content(
            str(segment.query_sequence),
            bin_size=bin_size,
        )
    ] += 1
    return existing_tally


def get_galignment_tally(
    segment: pysam.libcalignedsegment.AlignedSegment,
    existing_tally: defaultdict[tuple[str | None, int], int] | None = None,
    prev: tuple[str | None, int] | None = None,
    unique: int = 0,
    bin_size: int = 1000000,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Get tally for good alignment.

    Args:
        segment: pysam alignment segment.
        existing_tally: A tally of counts for good alignments.
        prev: (reference_name, reference_start + 1) from the previous line.
        unique: amount of unique measurements.
        bin_size: default bin size for alignment.

    Returns:
        An updated tally of counts in each G / C fraction bin

    """
    if existing_tally is None:
        existing_tally = defaultdict(int)

    if (segment.reference_name, segment.reference_start + 1) != prev:
        unique += 1
        prev = (segment.reference_name, segment.reference_start + 1)
        existing_tally[(segment.reference_name, (segment.reference_start + 1) // bin_size)] += 1
    return existing_tally, prev, unique


def parse_samtools_stats(
    sample_id: str,
    samtool_stats: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Parse file generated by samtools stats for error rate and gc stats.

    Args:
        sample_id: sample ID.
        samtool_stats: samtool stats file created from `samtools stats <>`.

    Returns:
        Error rate dataframe and gcd stats dataframe.

    """
    gcd_stats: pd.DataFrame = pd.DataFrame()

    with open(samtool_stats) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("SN	error rate:"):
                er = float(line.split("\t")[2])
                error_rate = pd.DataFrame([[sample_id, er]])

            if line.startswith("GCD"):
                gcd = pd.read_csv(StringIO(line), sep="\t", header=None)
                gcd[0] = sample_id
                gcd_stats = pd.concat([gcd_stats, gcd], ignore_index=True)
    return error_rate, gcd_stats


def load_config(
    *,
    input_dir: PathLike | None = None,
    output_dir: PathLike | None = None,
    sampleid: str | None = None,
    runid: str | None = None,
    config: SomePathLikes | None = None,
    set_config: list[str] | None = None,
) -> MutableStrMapping:
    """Return a unified configuration for this command given command line arguments.

    Combines CLI options with schema and configuration to produce an immutable, unified, and validated configuration
    object for this command using the config utility module.
    The configuration schema for this command is located within the
    pipline package in `askcell.config.schema` sub-package.

    Args:
        input_dir: Optionally override path for input directory
        output_dir: Optionally override path for output directory
        sampleid: Optionally override sample ID
        runid: Optionaly override run ID
        config: Optionally override default package configs built from schema using data files
        set_config: Nested key assignment expressions of the form key[.key]*=value used to override data values

    Returns:
        Mutable nested dictionary-like objected

    """
    params = {
        "params": {
            "input_dir": input_dir,
            "output_dir": output_dir,
        },
        "sample": {
            "sampleid": sampleid,
            "runid": runid,
        },
    }

    # Build setting data in reverse precedence order (low to high)
    data: list[MutableStrMapping] = list(load_config_data(config)) if config is not None else []

    # Load and validate concentration command configuration from data stored within the package
    package_schema = load_package_schema(PackageResource("askcell.config.schema", "wgs.toml"))
    package_data = load_package_data(
        package_schema,
        data=data + [params],
        settings=set_config,
    )

    save_toml(
        Path(package_data["output"]["output_dir"]) / "config.toml",
        unparse_data(package_schema, package_data),
    )

    return package_data


def run_wgs_samples(
    args: argparse.Namespace,
) -> None:
    """Run wga workflow.

    Args:
        args: Parameters parsed by argparse

    """
    # FIXME: Run 1 sample at a time, refactor later for multi-processing

    args.output_dir.mkdir(exist_ok=True, parents=True)

    config = load_config(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        config=args.config,
        set_config=args.set_config,
    )
    output_dir = parse_path(config["output"]["output_dir"])

    complete_genome_coverage = pd.DataFrame()

    with tempfile.TemporaryDirectory() as tmpdir, open(args.samplesheet) as ss:
        reader = csv.DictReader(ss, delimiter="\t")
        base_input_dir = parse_path(config["input"]["input_dir"])

        for row in reader:
            sampleid = row["SampleName"]
            config["sample"]["sampleid"] = sampleid
            runid = row["SequencingRun"]
            config["sample"]["run_id"] = runid
            config["scratch"]["output_dir"] = parse_path(tmpdir) / sampleid
            config["sample"]["input_dir"] = base_input_dir / runid

            results = run_wgs_sample(args, config)
            norm_base_line = pd.read_csv(
                parse_path(WGS_BASELINE),
                header=None,
            )

            # Replacement for R script (wgs_postprocessing_titan2.R)
            cov_bias = get_wgs_cov_bias(
                norm_base_line,
                normalization_samples=config["get_wgs_cov_bias"]["options"]["norm_samples"],
            )

            wgs_scaled_coverage = get_wgs_scaled_coverage(
                results["genome_mb_coverage"],
            )

            merged_genome_coverage = merge_bias_cov(
                cov_bias,
                wgs_scaled_coverage,
                bias_cv_threshold=config["merge_bias_cv"]["options"]["bias_cv_threshold"],
                include_chromosome=config["merge_bias_cv"]["options"]["include_chromosome"],
            )

            complete_genome_coverage = pd.concat([complete_genome_coverage, merged_genome_coverage])

        (output_dir / "total").mkdir(exist_ok=True, parents=True)
        complete_genome_coverage.to_csv(
            output_dir / "total" / "complete_genome_coverage.csv",
            index=False,
        )


def run_wgs_sample(
    args: argparse.Namespace,
    config: MutableStrMapping | None = None,
) -> dict[str, Any]:
    """Process a single whole genome sample analysis."""
    args.output_dir.mkdir(exist_ok=True, parents=True)

    if config is None:  # entry point is from "askcell wgs-sample"
        config = load_config(
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            runid=args.runid,
            sampleid=args.sample_id,
            config=args.config,
            set_config=args.set_config,
        )

    input_dir = parse_path(config["sample"]["input_dir"])
    output_dir = parse_path(config["output"]["output_dir"])
    scratch_dir = parse_path(config["scratch"]["output_dir"] or parse_path(tempfile.mkdtemp()))

    sampleid = config["sample"]["sampleid"]

    alignment_dir = output_dir / "alignments"
    logs = output_dir / "log"
    metrics_dir = output_dir / "metrics"
    metrics_out = metrics_dir / f"{sampleid}"
    ichorcna_dir = output_dir / "ichorcna"
    ichorcna_out = ichorcna_dir / f"{sampleid}"

    scratch_dir.mkdir(exist_ok=True, parents=True)
    alignment_dir.mkdir(exist_ok=True, parents=True)
    logs.mkdir(exist_ok=True, parents=True)
    metrics_dir.mkdir(exist_ok=True, parents=True)
    metrics_out.mkdir(exist_ok=True, parents=True)
    ichorcna_dir.mkdir(exist_ok=True, parents=True)
    ichorcna_out.mkdir(exist_ok=True, parents=True)

    # reads
    read1 = only_one(input_dir.glob(f"Fastq/{sampleid}*_R1_001.fastq.gz"))
    read2 = only_one(input_dir.glob(f"Fastq/{sampleid}*_R2_001.fastq.gz"))

    # Aligned reads
    aligned_read = alignment_dir / f"{sampleid}.so.bam"
    aligned_log = logs / f"{sampleid}.alignments.log"

    # Metrics
    dedup_metrics = metrics_out / f"{sampleid}.dedup_metrics.txt"
    samtool_stats = metrics_out / f"{sampleid}_samtool_stats.txt"

    # readCounter stats
    counter_stats = ichorcna_out / f"{sampleid}.wig"

    # construct list of commands
    cmds = make_commands_sample_analysis(
        sampleid=sampleid,
        read1=read1,
        read2=read2,
        dedup_metrics=dedup_metrics,
        aligned_read=aligned_read,
        aligned_log=aligned_log,
        scratch_dir=scratch_dir,
        threads=args.threads,
        min_mapq=config["parse_bam"]["options"]["min_mapq"],
        counter_stats=counter_stats,
        chromosome=config["readcounter"]["options"]["chromosome"],
        alignment_bsize=config["parse_bam"]["options"]["alignment_bsize"],
        min_counterquality=config["readcounter"]["options"]["quality"],
    )

    cmds.append(
        make_command_samtools_stats(
            aligned_read=aligned_read,
            output_stats=samtool_stats,
            threads=args.threads,
        )
    )

    cmds.append(
        make_command_run_ichorcna(
            sampleid=sampleid,
            bigwig=counter_stats,
            output=ichorcna_out,
            ploidy=config["ichorcna"]["options"]["ploidy"],
            normal_fraction=config["ichorcna"]["options"]["normal_fraction"],
            include_homd=config["ichorcna"]["flags"]["include_homd"],
            chrs=config["ichorcna"]["options"]["chrs"],
            chr_train=config["ichorcna"]["options"]["chr_train"],
            estimate_normal=config["ichorcna"]["flags"]["estimate_normal"],
            estimate_ploidy=config["ichorcna"]["flags"]["estimate_ploidy"],
            estimate_scprevalence=config["ichorcna"]["flags"]["estimate_scprevalence"],
            txne=config["ichorcna"]["options"]["txne"],
            txn_strength=config["ichorcna"]["options"]["txn_strength"],
            logfile=logs / f"{sampleid}.ichorcna.log",
        )
    )

    for cmd in cmds:
        if args.dryrun:
            print(cmd)
        else:
            run_cmd(cmd)

    galignment, gc, duplicates = parse_bam(
        sampleid,
        aligned_read,
        min_mapq=config["parse_bam"]["options"]["min_mapq"],
        alignment_bin_size=config["parse_bam"]["options"]["alignment_bsize"],
        gc_bin_size=config["parse_bam"]["options"]["gc_bsize"],
        threads=args.threads,
    )

    error_rate, gcd_stats = parse_samtools_stats(
        sampleid,
        metrics_out / f"{sampleid}_samtool_stats.txt",
    )

    galignment.to_csv(metrics_out / "genome_mb_coverage.csv", index=False)
    gc.to_csv(metrics_out / "gc_histogram.csv", index=False)
    duplicates.to_csv(metrics_out / "duplicate_stats.csv", index=False)
    error_rate.to_csv(metrics_out / "mismatch_stats.csv", index=False)
    gcd_stats.to_csv(metrics_out / "gc_stats.csv", index=False)

    return {
        "genome_mb_coverage": galignment,
    }


def make_commands_sample_analysis(
    sampleid: str,
    *,
    reference_genome_index: PathLike = GENOME_INDEX,
    read1: PathLike,
    read2: PathLike,
    dedup_metrics: PathLike,
    aligned_read: PathLike,
    aligned_log: PathLike,
    scratch_dir: PathLike,
    threads: int,
    min_mapq: int,
    counter_stats: PathLike,
    chromosome: str,
    alignment_bsize: int,
    min_counterquality: int,
) -> list[Command]:
    """Make commands to process a single sample for whole-genome-sequencing analysis."""
    temp_dir = parse_path(scratch_dir)
    readgroup = f"@RG\\tID:{sampleid}\\tSM:{sampleid}\\tLB:WES\\tPL:Illumina\\tPU:run"
    aligned_q20_read = temp_dir / f"{sampleid}.q20.bam"
    dedup_read = temp_dir / f"{sampleid}.dedup.bam"

    return [
        make_command_bwa_mem(
            reference=reference_genome_index,
            read1=read1,
            read2=read2,
            aligned_read=aligned_read,
            readgroup=readgroup,
            temp_dir=temp_dir,
            logfile=aligned_log,
            min_mapq=min_mapq,
            threads=threads,
        ),
        make_command_samtools_view(
            aligned_read=aligned_read,
            output_read=aligned_q20_read,
            min_mapq=min_mapq,
        ),
        make_command_samtools_index(aligned_read=aligned_read),
        make_command_samtools_index(aligned_read=aligned_q20_read),
        make_command_picard_markduplicate(
            aligned_read=aligned_q20_read,
            duped_read=dedup_read,
            metrics=dedup_metrics,
        ),
        make_command_read_counter(
            aligned_read=aligned_read,
            output=counter_stats,
            chromosome=chromosome,
            window=alignment_bsize,
            quality=min_counterquality,
        ),
    ]
