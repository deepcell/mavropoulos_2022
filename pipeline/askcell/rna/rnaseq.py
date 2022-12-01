"""RNA-seq analysis workflow."""

from __future__ import annotations

import argparse
import csv
import tempfile

import numpy as np
import pandas as pd

from .. import only_one
from ..static_data import (
    AG_TRIM,
    GENCODE_GENE_BED,
    P5,
    R_SUMMARIZE,
    STAR_GENOME_INDEX,
    SWIFT_SID_SNP_DB,
    rRNA_FA,
)
from ..utils.pathutils import PathLike, parse_path
from ..utils.shellutils import Command, run_cmd
from ..wrappers import (
    make_command_aggregate,
    make_command_atropos_trim,
    make_command_bwa_mem,
    make_command_samtools_flagstat,
    make_command_samtools_view,
    make_command_star_alignreads,
)


DTYPES_ALIGNMENT_STATS = {
    "Uniquely mapped reads number": np.uint64,
    "Number of input reads": np.uint64,
    "Number of reads mapped to multiple loci": np.uint64,
    "Number of reads mapped to too many loci": np.uint64,
    "Number of reads unmapped: too short": np.uint64,
    "Number of reads unmapped: other": np.uint64,
    "Mismatch rate per base, %": str,
}

DTYPES_TRIMINFO = {
    "name": str,
    "num_errors": np.uint8,
    "adapter_start": np.uint32,
    "adapter_stop": np.uint32,
    "left_seq": str,
    "adapter_seq": str,
    "right_seq": str,
    "name_adapter": str,
    "q-left": str,
    "q-adapter": str,
    "q-right": str,
    "is_reverse": np.bool_,
}

DTYPES_NOTRIMINFO = {
    "name": str,
    "trim_flag": str,
    "read": str,
    "quality": str,
}

DTYPES_GENE_COUNTS = {
    "name": str,
    "N_FR": np.uint64,
    "N_FWD": np.uint64,
    "N_REV": np.uint64,
}


def run_rnaseq_samples(args: argparse.Namespace) -> None:
    """Run rnaseq workflow.

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

            run_rnaseq_sample(args)

    # aggregated gene counts
    expression_dir = args.output_dir / "expression"
    expression_dir.mkdir(exist_ok=True, parents=True)

    # counts = expression_dir / "aggregate.gene.counts.tsv"

    run_cmd(
        make_command_aggregate(
            script=R_SUMMARIZE,
            output_dir=args.output_dir,
            samplesheet=args.samplesheet,
            bulk_db=SWIFT_SID_SNP_DB,
            workflow="rna",
            genes=GENCODE_GENE_BED,
        ),
    )


def run_rnaseq_sample(args: argparse.Namespace) -> None:
    """Run rnaseq workflow for a single sample."""
    # FIXME: refactor into initialization step
    input_dir = args.input_dir
    output_dir = args.output_dir
    scratch_dir = args.scratch_dir or parse_path(tempfile.mkdtemp())

    sampleid = args.sampleid

    trim_dir = output_dir / "trimmed"
    alignment_dir = output_dir / "alignments" / sampleid
    metrics_dir = output_dir / "metrics"

    contamination_dir = metrics_dir / "rRNA_contam_alignments"
    logs = args.logs_dir or output_dir / "log"

    scratch_dir.mkdir(exist_ok=True, parents=True)
    trim_dir.mkdir(exist_ok=True, parents=True)
    alignment_dir.mkdir(exist_ok=True, parents=True)
    metrics_dir.mkdir(exist_ok=True, parents=True)
    logs.mkdir(exist_ok=True, parents=True)
    contamination_dir.mkdir(exist_ok=True, parents=True)

    # FIXME: all paths defined in configuration
    # input
    # TODO: if it is valid to have multiple fastqs (i.e. fq per lane),
    #       add concat function to put multiple fastq into a single one
    read = only_one(input_dir.glob(f"**/Fastq/{sampleid.replace('_', '-')}*.1.fastq.gz"))

    # Trimmed reads
    trim_read = trim_dir / f"{sampleid}.R1.trimmed.fastq.gz"
    trim_info = trim_dir / f"{sampleid}.trimmed.info.txt"
    trim_report = trim_dir / f"{sampleid}.trimmed.stats.txt"
    trim_log = logs / f"{sampleid}.adapter_trim.log"

    # Aligned reads
    aligned_read = alignment_dir / "Aligned.sortedByCoord.out.bam"
    # aligned_read_index = alignment_dir / "Aligned.sortedByCoord.out.bam.bai"
    # aligned_log = logs / f"{sampleid}.alignments.log"
    aligned_log_final = alignment_dir / "Log.final.out"
    # aligned_log_out = alignment_dir / "Log.out"
    # aligned_log_progress = alignment_dir / "Log.progress.out"
    reads_per_gene = alignment_dir / "ReadsPerGene.out.tab"

    # rRNA alignment
    # unmapped_fq = alignment_dir / "Unmapped.out.mate1"
    rRNA_aligned_read = contamination_dir / f"{sampleid}.unmapped_to_rRNA.bam"
    rRNA_aligned_log = logs / f"{sampleid}.rRNA.alignments.json"

    # FIXME: metrics, currently metrics_dir is input to collect_metrics
    metrics = metrics_dir / f"{sampleid}.metrics.csv"
    # metrics_log = logs / f"{sampleid}.collect_rna_metrics.log"

    # Construct list of commands
    cmds = make_commands_rnaseq_analysis(
        sampleid=sampleid,
        read=read,
        trim_read=trim_read,
        trim_info=trim_info,
        trim_report=trim_report,
        trim_log=trim_log,
        aligned_read=aligned_read,
        aligned_log=aligned_log_final,
        rRNA_aligned_read=rRNA_aligned_read,
        rRNA_aligned_log=rRNA_aligned_log,
        scratch_dir=scratch_dir,
        adapter1=args.adapter1,
        trim_algorithm="adapter",
        threads=args.threads,
    )

    for cmd in cmds:
        if args.dryrun:
            print(cmd)
        else:
            run_cmd(cmd)

    if not args.dryrun:
        summarize_stats(
            pd.concat(
                [
                    get_trim_stats(trim_info, sampleid),
                    get_aln_stats(aligned_log_final, sampleid),
                    get_gene_count_stats(reads_per_gene, sampleid),
                    get_rRNA_stats(rRNA_aligned_log, sampleid),
                ],
                axis=1,
            ),
        ).to_csv(metrics)


def make_commands_rnaseq_analysis(
    sampleid: str,
    *,
    genome_index: PathLike = STAR_GENOME_INDEX,
    rRNA_index: PathLike = rRNA_FA,
    panel: str = "wts",
    read: PathLike,
    trim_read: PathLike,
    trim_info: PathLike,
    trim_report: PathLike,
    trim_log: PathLike,
    aligned_read: PathLike,
    aligned_log: PathLike,
    rRNA_aligned_read: PathLike,
    rRNA_aligned_log: PathLike,
    scratch_dir: PathLike,
    adapter1: str,
    trim_algorithm: str,
    threads: int = 1,
) -> list[Command]:
    """Make commands to process a single sample for rnseq analysis."""
    readgroup = f"@RG\\tID:{sampleid}\\tSM:{sampleid}\\tLB:{panel}\\tPL:mn\\tPU:run"
    scratch_dir = parse_path(scratch_dir)

    return [
        make_command_atropos_trim(
            algorithm=trim_algorithm,
            adapter1=[P5, *AG_TRIM],
            read1=read,
            trimmed_read1=trim_read,
            trim_info=trim_info,
            trim_report=trim_report,
            logfile=trim_log,
            nextseq_trim_qscore=20,
            sampleid=sampleid,
            threads=threads,
        ),
        make_command_star_alignreads(
            reference=genome_index,
            read1=trim_read,
            aligned_read=aligned_read,
            threads=threads,
        ),
        make_command_bwa_mem(
            reference=rRNA_index,
            read1=trim_read,
            aligned_read=rRNA_aligned_read,
            readgroup=readgroup,
            logfile=rRNA_aligned_log,
            threads=threads,
        ),
        make_command_samtools_view(
            aligned_read=rRNA_aligned_read,
            exclude_flags=256,
            output_read=scratch_dir / "rRNA_filtered_temp.bam",
        ),
        make_command_samtools_flagstat(
            aligned_read=scratch_dir / "rRNA_filtered_temp.bam",
            stat=rRNA_aligned_log,
        ),
    ]


def get_rRNA_stats(
    filename: PathLike,
    sampleid: str,
) -> pd.DataFrame:
    """Get rRNA alignment stat.

    Args:
        filename: rRNA alignment flagstat filename
        sampleid: sample_identifier

    Returns:
        rRNA mapped read, secondary & mapped

    """
    return (
        pd.read_json(
            filename,
            orient="index",
            dtype={"secondary": np.uint64, "mapped": np.uint64},
        )
        .loc[["QC-passed reads"], ["secondary", "mapped"]]
        .set_index(pd.Index([sampleid]))
    )


def get_aln_stats(
    filename: PathLike,
    sampleid: str,
) -> pd.DataFrame:
    r"""Get alignment stats information, from **Log.final.out**.

    Strips off trailing " \|" from first column.

    Example output from Log.final.out)
    UNIQUE READS:
    Uniquely mapped reads number |       558411
    Uniquely mapped reads % |       81.53%

    Args:
        filename: alignment final log file
        sampleid: sample identifier

    Returns:
        alignment statistics

    """
    return pd.read_csv(
        filename,
        sep="\t",
        index_col=0,
        names=["Name", sampleid],
        converters={"Name": lambda x: x.strip(" |")},
    ).T.astype(DTYPES_ALIGNMENT_STATS)


def get_gene_count_stats(
    filename: PathLike,
    sampleid: str,
) -> pd.DataFrame:
    """Get gene count stats information.

    Args:
        filename: gene counts file
        sampleid: sample identifier

    Returns:
        gene count statistics, total count from foward and reverse transcriptome data

    """
    d = pd.read_csv(
        filename,
        sep="\t",
        header=None,
        names=["name", "N_FR", "N_FWD", "N_REV"],
        dtype=DTYPES_GENE_COUNTS,
    )
    # Name column starting with letter N is not the expression count to be used
    d = pd.DataFrame(d.loc[~d["name"].str.startswith("N"), ["N_FWD", "N_REV"]].sum(), dtype=np.uint64).T
    return d.set_index(pd.Index([sampleid]))


def summarize_stats(data: pd.DataFrame) -> pd.DataFrame:
    """Summarize rnaseq stats.

    # TODO: Define reportable metrics, return with proper dtype with dtype conversion in future PR

    Args:
        data: stats output

    Returns:
        summary metrics statistics for the report

    """
    metrics = data[
        [
            "Number of input reads",
            "Mismatch rate per base, %",
            "MeanFragLen",
            "MeanAdapterLen",
            "PCTAdapter",
            "PCT_A",
            "PCT_G",
        ]
    ].rename(
        columns={
            "Number of input reads": "NumFrags",
            "Mismatch rate per base, %": "MismatchRate",
        },
    )

    n_unmap = data[
        [
            "Number of reads unmapped: too short",
            "Number of reads unmapped: other",
            "Number of reads mapped to too many loci",
        ]
    ]

    # FIXME: Add "PCT_mtRNA"
    metrics["PCT_UniqMapped"] = 100 * data["Uniquely mapped reads number"] / data["Number of input reads"]
    metrics["PCT_MultiMapped"] = 100 * data["Number of reads mapped to multiple loci"] / data["Number of input reads"]
    metrics["PCT_Unmapped"] = 100 * n_unmap.sum(axis=1) / data["Number of input reads"]
    metrics["PCT_FwdStranded"] = 100 * data["N_FWD"] / data["Number of input reads"]
    metrics["PCT_RevStranded"] = 100 * data["N_REV"] / data["Number of input reads"]
    metrics["PCT_rRNA"] = 100 * (data["mapped"] - data["secondary"]) / data["Number of input reads"]

    return metrics.rename_axis("SampleName", axis="index")


def get_trim_stats(
    filename: PathLike,
    sampleid: str,
) -> pd.DataFrame:
    """Get trim stats information.

    Args:
        filename: trim info output
        sample_id: sample identifier

    Returns:
        trim statistics

    """
    trimmed_rows = []
    untrimmed_rows = []
    with open(filename) as fin:
        for row in fin:
            r = row.rstrip("\n").split("\t")
            # Uneven number of columns, sometimes 4 columns if not trimmed, else 14 if trimmed
            if len(r) == 4:
                untrimmed_rows.append(r)
            else:
                # empty column is reported as truncated column
                trimmed_rows.append(r + [""] * (len(DTYPES_TRIMINFO) - len(r)))

    # FIXME: do proper typing during input, not just astype
    adapter_reads = pd.DataFrame(
        trimmed_rows,
        columns=DTYPES_TRIMINFO,
    ).astype(DTYPES_TRIMINFO)

    noadapter_reads = pd.DataFrame(
        untrimmed_rows,
        columns=DTYPES_NOTRIMINFO,
    ).astype(DTYPES_NOTRIMINFO)

    adapter_reads_cnt = len(adapter_reads.index)
    total_reads_cnt = adapter_reads_cnt + len(noadapter_reads.index)

    adapter_reads["adapter+"] = adapter_reads["adapter_seq"] + adapter_reads["right_seq"]
    len_adapter = adapter_reads["adapter+"].str.len().sum()
    len_remain = noadapter_reads["read"].str.len().sum() + adapter_reads["left_seq"].str.len().sum()

    # lt50_pct = 100 * len(d.loc[lt50_filter].index) / total_reads_cnt
    return pd.DataFrame(
        {
            "MeanFragLen": pd.Series([len_remain / total_reads_cnt], dtype=np.float64),
            "MeanAdapterLen": pd.Series([len_adapter / adapter_reads_cnt], dtype=np.float64),
            "PCTAdapter": pd.Series([100 * adapter_reads_cnt / total_reads_cnt], dtype=np.float64),
            "PCT_A": pd.Series(
                [100 * adapter_reads["adapter+"].str.count("A").sum() / len_adapter],
                dtype=np.float64,
            ),
            "PCT_G": pd.Series(
                [100 * adapter_reads["adapter+"].str.count("G").sum() / len_adapter],
                dtype=np.float64,
            ),
        },
    ).set_index(pd.Index([sampleid]))
