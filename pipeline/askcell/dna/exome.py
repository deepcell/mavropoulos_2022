"""Exome sequencing analysis workflow."""

from __future__ import annotations

import argparse
import csv
import tempfile

from .. import only_one
from ..static_data import (
    GENOME,
    GENOME_DICT,
    GENOME_INDEX,
    GENOME_SDF,
    GNOMAD,
    REFERENCE_CONTROLS_CALLS,
    REFERENCE_CONTROLS_CONFIDENT_REGIONS,
    TWIST_EXOME_TARGET,
)
from ..utils.pathutils import PathLike, parse_path
from ..utils.shellutils import Command, run_cmd
from ..wrappers import (
    make_command_bwa_mem,
    make_command_gatk_haplotypecaller,
    make_command_gatk_recalibrate,
    make_command_picard_collecthsmetrics,
    make_command_picard_markduplicate,
    make_command_rtg_vcfeval,
    make_command_samtools_index,
    make_command_samtools_view,
)


def run_exome_samples(args: argparse.Namespace) -> None:
    """Run exome sequencing workflow.

    Currently, Twist exome 2.0 is used.

    """
    # FIXME: Run 1 sample at a time, refactor later for multi-processing
    args.output_dir.mkdir(exist_ok=False, parents=True)

    with tempfile.TemporaryDirectory() as tmpdir, open(args.samplesheet) as ss:
        reader = csv.DictReader(ss, delimiter="\t")
        for row in reader:
            args.sampleid = row["SampleName"]
            args.runid = row["SequencingRun"]
            args.control_name = row["ControlName"]
            args.scratch_dir = parse_path(tmpdir) / args.sampleid

            run_exome_sample(args)


def run_exome_sample(args: argparse.Namespace) -> None:
    """Process a single exome sample analysis."""
    input_dir = args.input_dir
    output_dir = args.output_dir
    scratch_dir = args.scratch_dir or parse_path(tempfile.mkdtemp())

    sampleid = args.sampleid
    control_name = args.control_name

    is_control = True if control_name in REFERENCE_CONTROLS_CALLS else False

    alignment_dir = output_dir / "alignments"
    recal_dir = output_dir / "recal"
    var_dir = output_dir / "calls"
    metrics_dir = output_dir / "metrics"
    metrics_out = metrics_dir / f"{sampleid}"
    logs = args.logs_dir or output_dir / "log"

    scratch_dir.mkdir(exist_ok=True, parents=True)
    alignment_dir.mkdir(exist_ok=True, parents=True)
    recal_dir.mkdir(exist_ok=True, parents=True)
    var_dir.mkdir(exist_ok=True, parents=True)
    metrics_out.mkdir(exist_ok=True, parents=True)
    logs.mkdir(exist_ok=True, parents=True)

    if is_control:
        control_dir = output_dir / f"control_metrics/{sampleid}/giab"
        # rtg vcfeval throws an error, if output directory already exists
        # control_dir.mkdir(exist_ok=True, parents=True)

    # reads
    read1 = only_one(input_dir.glob(f"**/Fastq/{sampleid}*_R1_001.fastq.gz"))
    read2 = only_one(input_dir.glob(f"**/Fastq/{sampleid}*_R2_001.fastq.gz"))

    # Aligned reads
    aligned_read = alignment_dir / f"{sampleid}.so.bam"
    recalib_read = alignment_dir / f"{sampleid}.recal.bam"
    aligned_log = logs / f"{sampleid}.alignments.log"

    # Metrics
    recal_table = metrics_out / f"{sampleid}.recal_data.table"
    dedup_metrics = metrics_out / f"{sampleid}.dedup_metrics.txt"
    hsmetrics = metrics_out / f"{sampleid}.hsmetrics.txt"
    base_coverage = metrics_out / f"{sampleid}.base.txt"
    target_coverage = metrics_out / f"{sampleid}.target.txt"

    # Calls
    calls = var_dir / f"{sampleid}.vcf.gz"

    # construct list of commands
    cmds = make_commands_sample_analysis(
        sampleid=sampleid,
        reference=GENOME,
        reference_index=GENOME_INDEX,
        reference_dict=GENOME_DICT,
        target=TWIST_EXOME_TARGET,
        known_sites=GNOMAD,
        read1=read1,
        read2=read2,
        aligned_read=aligned_read,
        aligned_log=aligned_log,
        recalib_read=recalib_read,
        recalib_table=recal_table,
        dedup_metrics=dedup_metrics,
        hsmetrics=hsmetrics,
        base_coverage=base_coverage,
        target_coverage=target_coverage,
        calls=calls,
        scratch_dir=scratch_dir,
        threads=args.threads,
    )

    if is_control:
        cmds += [
            make_command_rtg_vcfeval(
                baseline=REFERENCE_CONTROLS_CALLS[control_name],
                called_variants=calls,
                reference=GENOME_SDF,
                outdir=control_dir,
                regions=REFERENCE_CONTROLS_CONFIDENT_REGIONS[control_name],
                evaluation_regions=TWIST_EXOME_TARGET,
            ),
        ]

    for cmd in cmds:
        if args.dryrun:
            print(cmd)
        else:
            run_cmd(cmd)


def make_commands_sample_analysis(
    sampleid: str,
    *,
    reference: PathLike,
    reference_index: PathLike,
    reference_dict: PathLike,
    target: PathLike,
    known_sites: PathLike,
    read1: PathLike,
    read2: PathLike,
    aligned_read: PathLike,
    aligned_log: PathLike,
    recalib_read: PathLike,
    recalib_table: PathLike,
    dedup_metrics: PathLike,
    hsmetrics: PathLike,
    base_coverage: PathLike,
    target_coverage: PathLike,
    calls: PathLike,
    scratch_dir: PathLike,
    threads: int,
) -> list[Command]:
    """Make commands to process a single sample for exome-sequencing analysis."""
    temp_dir = parse_path(scratch_dir)
    readgroup = f"@RG\\tID:{sampleid}\\tSM:{sampleid}\\tLB:WES\\tPL:Illumina\\tPU:run"
    aligned_q20_read = temp_dir / f"{sampleid}.q20.bam"
    dedup_read = temp_dir / f"{sampleid}.dedup.bam"

    return [
        make_command_bwa_mem(
            reference=reference_index,
            read1=read1,
            read2=read2,
            aligned_read=aligned_read,
            readgroup=readgroup,
            temp_dir=temp_dir,
            logfile=aligned_log,
            threads=threads,
        ),
        make_command_samtools_view(
            aligned_read=aligned_read,
            output_read=aligned_q20_read,
            min_mapq=20,
        ),
        make_command_samtools_index(aligned_read=aligned_q20_read),
        make_command_picard_markduplicate(
            aligned_read=aligned_q20_read,
            duped_read=dedup_read,
            metrics=dedup_metrics,
        ),
        make_command_samtools_index(aligned_read=dedup_read),
        make_command_picard_collecthsmetrics(
            reference=reference,
            reference_dict=reference_dict,
            aligned_read=dedup_read,
            metrics=hsmetrics,
            bait_bed=target,
            target_bed=target,
            tmp_dir=temp_dir,
            base_coverage=base_coverage,
            target_coverage=target_coverage,
        ),
        make_command_gatk_recalibrate(
            reference=reference,
            known_sites=known_sites,
            aligned_read=dedup_read,
            calibrated_read=recalib_read,
            calibrated_table=recalib_table,
            intervals=target,
        ),
        make_command_gatk_haplotypecaller(
            reference=reference,
            aligned_read=recalib_read,
            calls=calls,
            intervals=target,
            no_soft_clipped_bases=True,
        ),
    ]
