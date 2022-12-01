"""Amplicon target (such as swift panel) DNA analysis workflow."""

from __future__ import annotations

import argparse
import csv
import tempfile

from .. import only_one
from ..static_data import (
    GENOME,
    GENOME_DICT,
    GENOME_INDEX,
    GNOMAD,
    PANELS_AMPLICON_KEY,
    PANELS_PRIMER_KEY,
    PANELS_TARGET_KEY,
    R_SUMMARIZE,
    SWIFT_SID_SNP_DB,
)
from ..utils.pathutils import PathLike, parse_path
from ..utils.shellutils import Command, run_cmd
from ..wrappers import (
    make_command_aggregate,
    make_command_atropos_trim,
    make_command_bcftools_call,
    make_command_bcftools_index,
    make_command_bcftools_mpileup,
    make_command_bwa_mem,
    make_command_call_germline,
    make_command_call_mutect2,
    make_command_check_swap_metrics,
    make_command_collect_metrics,
    make_command_filter_mutect2,
    make_command_picard_collecthsmetrics,
    make_command_samtools_ampliconclip,
    make_command_samtools_index,
)


# FIXME: lacks logging functionality, lots of log outputs from askcell.sh is missing


def run_tgs_samples(args: argparse.Namespace) -> None:
    """Run tga workflow for oncology analysis.

    Currently, supported panels include 'swift_72g', 'swift_lung', and 'swift_brca'.
    Type definitions: S = post-sort, P = pre-sort, B = bulk

    """
    # FIXME: Run 1 sample at a time, refactor later for multi-processing
    args.output_dir.mkdir(exist_ok=False, parents=True)

    with tempfile.TemporaryDirectory() as tmpdir, open(args.samplesheet) as ss:
        reader = csv.DictReader(ss, delimiter="\t")
        base_input_dir = args.input_dir
        matched_seen = set()

        for row in reader:
            args.sampleid = row["SampleName"]
            args.runid = row["SequencingRun"]
            args.scratch_dir = parse_path(tmpdir) / args.sampleid
            args.input_dir = only_one(base_input_dir.glob(f"**/{args.runid}"))

            args.matchedid = row["MatchedNormalSampleName"]
            args.matched_runid = row["MatchedNormalSequencingRun"]
            args.matched_dir = only_one(base_input_dir.glob(f"**/{args.matched_runid}"))
            args.matched_exist = args.matchedid in matched_seen

            args.panel = row["Panel"]

            run_tgs_sample(args)
            matched_seen.add(args.matchedid)

    run_cmd(
        make_command_aggregate(
            script=R_SUMMARIZE,  # FIXME: no hard-coded filepath
            output_dir=args.output_dir,
            samplesheet=args.samplesheet,
            bulk_db=SWIFT_SID_SNP_DB,
            workflow="tgs",
        ),
    )


def run_tgs_sample(args: argparse.Namespace) -> None:
    """Process a single paired sample analysis."""
    # FIXME: refactor into initialization step for directory creation
    # FIXME: more robust input check, to ensure that fastqs are uniquely available
    # FIXME: also needs a smarter way to find panel, i.e. swift-sid fails, while swift_sid is recognized
    if args.panel not in PANELS_TARGET_KEY:
        # FIXME: panel check is limited, that it doesn't recognize unless available from static_data.py
        raise ValueError(f"{args.panel} is not supported")

    args.scratch_dir = args.scratch_dir or parse_path(tempfile.mkdtemp())

    # tumor and matched sample (normal) runs, default non-paired portion is same as sid-sample analysis
    args1 = argparse.Namespace(**vars(args))
    args2 = argparse.Namespace(**vars(args))
    args2.sampleid = args.matchedid
    args2.input_dir = args.matched_dir

    # FIXME: this runs extra germline call, not part of tgs workflow, address this in future PR
    run_sid_sample(args1)
    if "matched_exist" not in vars(args) or not args.matched_exist:
        run_sid_sample(args2)

    # somatic calls & swap checks that requires both tumor and normal data
    run_paired_sample(args)


def run_sid_samples(args: argparse.Namespace) -> None:
    """Run sid workflow for genotyping analysis.

    Currently, supported panels include 'swift_sid' for sid workflow.
    Type definitions: S = post-sort, P = pre-sort, B = bulk

    """
    # FIXME: Run 1 sample, refactor later for multi-processing
    args.output_dir.mkdir(exist_ok=False, parents=True)

    with tempfile.TemporaryDirectory() as tmpdir, open(args.samplesheet) as ss:
        reader = csv.DictReader(ss, delimiter="\t")
        for row in reader:
            args.sampleid = row["SampleName"]
            args.runid = row["SequencingRun"]
            args.scratch_dir = parse_path(tmpdir)
            args.panel = row["Panel"]

            run_sid_sample(args)

    run_cmd(
        make_command_aggregate(
            script=R_SUMMARIZE,  # FIXME: no hard-coded filepath
            output_dir=args.output_dir,
            samplesheet=args.samplesheet,
            bulk_db=SWIFT_SID_SNP_DB,
            workflow="sid",
        ),
    )


def run_sid_sample(args: argparse.Namespace) -> None:
    """Process a single sample for default amplicon analysis."""
    # FIXME: refactor into initialization step for directory creation, input check
    if args.panel not in PANELS_TARGET_KEY:
        raise ValueError(f"{args.panel} is not supported")

    input_dir = args.input_dir
    output_dir = args.output_dir
    scratch_dir = args.scratch_dir or parse_path(tempfile.mkdtemp())

    sampleid = args.sampleid

    trim_dir = output_dir / "trimmed"
    alignment_dir = output_dir / "alignments"
    var_dir = output_dir / "germlines"
    metrics_dir = output_dir / "metrics"
    metrics_out = metrics_dir / f"{sampleid}"
    plots_dir = output_dir / "plots"
    logs = args.logs_dir or output_dir / "log"

    scratch_dir.mkdir(exist_ok=True, parents=True)
    trim_dir.mkdir(exist_ok=True, parents=True)
    alignment_dir.mkdir(exist_ok=True, parents=True)
    var_dir.mkdir(exist_ok=True, parents=True)
    metrics_dir.mkdir(exist_ok=True, parents=True)
    metrics_out.mkdir(exist_ok=True, parents=True)
    plots_dir.mkdir(exist_ok=True, parents=True)
    logs.mkdir(exist_ok=True, parents=True)

    # FIXME: all paths defined in configuration
    # input

    # TODO: if fastqs from multiple lines go into separate file,
    #       thus valid to have multiple fastq outputs for a single sample
    #       add concat function to put multiple fastqs into a single one
    read1 = only_one(input_dir.glob(f"**/Fastq/{sampleid.replace('_', '-')}*_R1_001.fastq.gz"))
    read2 = only_one(input_dir.glob(f"**/Fastq/{sampleid.replace('_', '-')}*_R2_001.fastq.gz"))

    # Trimmed reads
    trim_read1 = trim_dir / f"{sampleid}.R1.trimmed.fastq.gz"
    trim_read2 = trim_dir / f"{sampleid}.R2.trimmed.fastq.gz"
    trim_info = trim_dir / f"{sampleid}.trimmed.info.txt"
    trim_report = trim_dir / f"{sampleid}.trimmed.stats.txt"
    trim_log = logs / f"{sampleid}.adapter_trim.log"

    # Aligned reads
    aligned_read = alignment_dir / f"{sampleid}.so.bam"
    aligned_log = logs / f"{sampleid}.alignments.log"

    # Naive caller
    counts = var_dir / f"{sampleid}.swift_sid.cnt.csv"
    # FIXME: calls output is currently determined by replacing cnt.csv to snp.tsv from naive caller
    # calls = var_dir / f"{sampleid}.swift_sid.snp.tsv"

    # FIXME: metrics, currently metrics_dir is input to collect_metrics
    # metrics = metrics_dir / f"{sampleid}.metrics.csv"
    # metrics_log = logs / f"{sampleid}.collect_dna_metrics.log"
    hsmetrics = metrics_out / f"{sampleid}.hsmetrics.txt"
    base_coverage = metrics_out / f"{sampleid}.base_coverage.txt"
    target_coverage = metrics_out / f"{sampleid}.target_coverage.txt"

    # Construct list of commands
    cmds = make_commands_sid_analysis(
        sampleid=sampleid,
        target=PANELS_TARGET_KEY[args.panel],
        amplicon_target=PANELS_AMPLICON_KEY[args.panel],
        primer_target=PANELS_PRIMER_KEY[args.panel],
        panel=args.panel,
        read1=read1,
        read2=read2,
        trim_read1=trim_read1,
        trim_read2=trim_read2,
        trim_info=trim_info,
        trim_report=trim_report,
        trim_log=trim_log,
        aligned_read=aligned_read,
        aligned_log=aligned_log,
        calls=counts,
        metrics_dir=metrics_dir,
        scratch_dir=scratch_dir,
        hsmetrics=hsmetrics,
        base_coverage=base_coverage,
        target_coverage=target_coverage,
        adapter1=args.adapter1,
        adapter2=args.adapter2,
        trim_algorithm="insert",
        threads=args.threads,
    )

    for cmd in cmds:
        if args.dryrun:
            print(cmd)
        else:
            run_cmd(cmd)


def run_paired_sample(args: argparse.Namespace) -> None:
    """Process a paired sample analysis, do somatic calling and swap check."""
    # FIXME: refactor into initialization step for directory creation, input check
    if args.panel not in PANELS_TARGET_KEY:
        raise ValueError(f"{args.panel} is not supported")

    # input_dir = args.input_dir
    output_dir = args.output_dir

    sampleid = args.sampleid
    matchedid = args.matchedid

    var_dir = output_dir / "somatics" / sampleid
    metrics_dir = output_dir / "metrics"
    swapcheck_dir = output_dir / "swap_checker" / f"{sampleid}.vs.{matchedid}"
    alignment_dir = output_dir / "alignments"

    var_dir.mkdir(exist_ok=True, parents=True)
    swapcheck_dir.mkdir(exist_ok=True, parents=True)
    metrics_dir.mkdir(exist_ok=True, parents=True)

    # FIXME: all paths defined in configuration,
    # some output definitions are intentionallly commented out, instead being deleted, for config use later

    # input
    normal_aligned = alignment_dir / f"{matchedid}.so.bam"
    tumor_aligned = alignment_dir / f"{sampleid}.so.bam"

    # swap-check output
    swap_calls = var_dir / f"{sampleid}.vs.{matchedid}.vcf.gz"
    swap_metrics = var_dir / f"{sampleid}.vs.{matchedid}.sample_swap_check.metrics.tsv"
    # swap_log = logs / f"{sampleid}.vs.{matchedid}.sample_swap_check.log"

    # mutect2
    calls = var_dir / f"{sampleid}.somatics.raw.vcf"
    evidence_bam = var_dir / f"{sampleid}.somatics.bamout"
    filtered_calls = var_dir / f"{sampleid}.somatics.filtered.vcf"
    call_stats = var_dir / f"{sampleid}.somatics.raw.vcf.stats"
    # call_log = logs / f"{sampleid}.vs.{matchedid}.somatics.log"

    # Construct list of commands
    cmds = make_commands_paired_analysis(
        sampleid=sampleid,
        matchedid=matchedid,
        target=PANELS_TARGET_KEY[args.panel],
        normal_aligned_read=normal_aligned,
        tumor_aligned_read=tumor_aligned,
        somatic_calls=calls,
        filtered_somatic_calls=filtered_calls,
        somatic_call_stats=call_stats,
        evidence_bam=evidence_bam,
        swap_calls=swap_calls,
        swap_metrics=swap_metrics,
        threads=args.threads,
    )

    for cmd in cmds:
        if args.dryrun:
            print(cmd)
        else:
            run_cmd(cmd)


def make_commands_sid_analysis(
    sampleid: str,
    *,
    reference_genome: PathLike = GENOME,
    reference_genome_index: PathLike = GENOME_INDEX,
    reference_dict: PathLike = GENOME_DICT,
    target: PathLike,
    amplicon_target: PathLike,
    primer_target: PathLike,
    panel: str,
    read1: PathLike,
    read2: PathLike,
    trim_read1: PathLike,
    trim_read2: PathLike,
    trim_info: PathLike,
    trim_report: PathLike,
    trim_log: PathLike,
    aligned_read: PathLike,
    aligned_log: PathLike,
    calls: PathLike,
    metrics_dir: PathLike,
    scratch_dir: PathLike,
    hsmetrics: PathLike,
    base_coverage: PathLike,
    target_coverage: PathLike,
    adapter1: str,
    adapter2: str,
    trim_algorithm: str,
    threads: int,
) -> list[Command]:
    """Make commands to process a single sample for swift-sid analysis."""
    readgroup = f"@RG\\tID:{sampleid}\\tSM:{sampleid}\\tLB:{panel}\\tPL:Illumina\\tPU:run"
    aligned_read_with_primer = parse_path(scratch_dir) / "aligned_w_primer.bam"

    return [
        make_command_atropos_trim(
            algorithm=trim_algorithm,
            adapter1=adapter1,
            adapter2=adapter2,
            read1=read1,
            read2=read2,
            trimmed_read1=trim_read1,
            trimmed_read2=trim_read2,
            trim_info=trim_info,
            trim_report=trim_report,
            logfile=trim_log,
            nextseq_trim_qscore=20,
            sampleid=sampleid,
            threads=threads,
        ),
        make_command_bwa_mem(
            reference=reference_genome_index,
            read1=trim_read1,
            read2=trim_read2,
            aligned_read=aligned_read_with_primer,
            readgroup=readgroup,
            temp_dir=scratch_dir,
            logfile=aligned_log,
            threads=threads,
        ),
        make_command_samtools_ampliconclip(
            reference=reference_genome,
            aligned_read=aligned_read_with_primer,
            clipped_read=aligned_read,
            primer=primer_target,
            tmp_dir=scratch_dir,
        ),
        make_command_samtools_index(aligned_read=aligned_read),
        make_command_picard_collecthsmetrics(
            reference=reference_genome,
            reference_dict=reference_dict,
            aligned_read=aligned_read,
            metrics=hsmetrics,
            bait_bed=target,
            target_bed=target,
            tmp_dir=scratch_dir,
            base_coverage=base_coverage,
            target_coverage=target_coverage,
        ),
        make_command_call_germline(
            reference=reference_genome,
            aligned_read=aligned_read,
            targets=target,
            calls=calls,
            sid=sampleid,
        ),
        make_command_collect_metrics(
            input_type="dna",
            reference=reference_genome,
            aligned_read=aligned_read,
            metrics_dir=metrics_dir,
            sid=sampleid,
            target=amplicon_target,
        ),
    ]


def make_commands_paired_analysis(
    sampleid: str,
    matchedid: str,
    *,
    reference_genome: PathLike = GENOME,
    genome_db: PathLike = GNOMAD,
    target: PathLike,
    normal_aligned_read: PathLike,
    tumor_aligned_read: PathLike,
    somatic_calls: PathLike,
    filtered_somatic_calls: PathLike,
    somatic_call_stats: PathLike,
    evidence_bam: PathLike,
    swap_calls: PathLike,
    swap_metrics: PathLike,
    threads: int,
) -> list[Command]:
    """Make commands to process a single sample, swift-tga paired analysis."""
    aligned_reads = (normal_aligned_read, tumor_aligned_read)

    cmd_pileup = make_command_bcftools_mpileup(
        reference=reference_genome,
        aligned_reads=aligned_reads,
        threads=threads,
        target=target,
        annotate="FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR",
        output_type="u",
        include_flags=2,
    )

    cmd_call = make_command_bcftools_call(
        output=swap_calls,
        output_type="z",
        use_multiallelic=True,
        out_variants_only=True,
    )

    return [
        cmd_pileup | cmd_call,
        make_command_bcftools_index(swap_calls),
        make_command_check_swap_metrics(swap_calls, swap_metrics),
        make_command_call_mutect2(
            reference=reference_genome,
            aligned_reads=aligned_reads,
            calls=somatic_calls,
            intervals=target,
            make_haplotype_graph=True,
            dont_use_soft_clipped_bases=True,
            germline_resource=genome_db,
            evidence_bam=evidence_bam,
            normal_sample_ids=matchedid,
            tumor_sample_ids=sampleid,
            filters="FragmentLengthReadFilter",
        ),
        make_command_filter_mutect2(
            reference=reference_genome,
            variants=somatic_calls,
            filtered_calls=filtered_somatic_calls,
            stats=somatic_call_stats,
            min_unique_alt_cnt=20,
            min_reads_per_strand=2,
            max_events_in_region=2,
            min_allele_fraction=0.01,
        ),
    ]
