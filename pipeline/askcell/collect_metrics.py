"""collect_metrics.py to collect metrics for dna and rna analysis."""

import argparse
import itertools

from typing import Any, cast

import numpy as np
import pandas as pd
import pyfaidx
import pyranges as pr
import pysam

from scipy.stats import median_abs_deviation

from .utils.pathutils import PathLike, parse_path


def parse_cigar(cigar_str: str) -> list[tuple[int, str]]:
    """Parse cigar string."""
    cigar_iter = itertools.groupby(cigar_str, lambda k: k.isdigit())
    cigar_list = []
    for _, n in cigar_iter:
        op = int("".join(n)), "".join(next(cigar_iter)[1])
        cigar_list.append(op)

    return cigar_list


def compute_aln_length(cigartuples: list[tuple[int, str]]) -> int:
    """Compute alignment length."""
    return sum(op_len for op_len, op in cigartuples if op in "MD")


def get_gc(row: pd.DataFrame, genome: pyfaidx.Fasta, flank: int = 200) -> pd.DataFrame:
    """Get GC content %."""
    chrom, start, end = row["Chromosome"], row["Start"], row["End"]
    seq = genome[chrom][start - 200 : end + 200].seq
    row["amplicon_length"] = end - start
    row["gc"] = len([c for c in seq if c.upper() in ["G", "C"]]) / len(seq)
    return row


def collect_frags(bam_file: PathLike) -> tuple[pd.DataFrame, int]:
    """Collect flags."""
    seen = set()
    n_aln = 1
    n_frags = 0
    frags = []
    alignments = pysam.AlignmentFile(str(bam_file))
    mdup_or_not = 0

    for k, v in alignments.header.items():  # type: ignore  # incomplete pysam typing
        if k == "PG":
            for i in range(len(v)):
                if v[i]["ID"] == "MarkDuplicates":
                    mdup_or_not = 1

    for aln in alignments.fetch():
        qname = aln.query_name
        if qname in seen:
            continue
        if aln.is_qcfail or aln.is_secondary or aln.is_supplementary:
            n_frags += 1
            seen.add(qname)
            continue
        if aln.is_unmapped or aln.mate_is_unmapped:
            n_frags += 1
            seen.add(qname)
            continue

        # FIXME: isn't itself_aln_len = aln.reference_stop - aln.reference_start?
        itself_aln_len = compute_aln_length(parse_cigar(aln.cigarstring or ""))
        itself_start = aln.reference_start
        mate_start = aln.next_reference_start
        # for all the rest of alignments, the MC tag should be present
        mate_cigar = cast(str, aln.get_tag("MC"))
        mate_aln_len = compute_aln_length(parse_cigar(mate_cigar))
        frag_aln_start = itself_start if itself_start <= mate_start else mate_start
        frag_aln_end = mate_start + mate_aln_len if itself_start <= mate_start else itself_start + itself_aln_len
        frag_len = frag_aln_end - frag_aln_start
        # fwd_end = itself_start + itself_aln_len if itself_start <= mate_start else mate_start + mate_aln_len
        # rev_start = mate_start if itself_start <= mate_start else itself_start
        # proper = 1, non_proper_same_chr = 0, non_proper_diff_chr = -1
        aln_type = 1 if aln.is_proper_pair else 0
        aln_type = -1 if aln.reference_name != aln.next_reference_name else aln_type
        # print(aln.cigarstring)
        # print(mate_cigar)
        # print(mate_aln_len)
        # print(aln.reference_name, frag_aln_start, frag_aln_end, frag_len)
        dup = 1 if aln.is_duplicate else 0
        seen.add(qname)
        row = [
            aln.reference_name,
            frag_aln_start,
            frag_aln_end,
            frag_len,
            aln_type,
            dup,
        ]
        frags.append(row)
        # fOUT.write("{}\n".format(
        # ",".join(list(map(str, [aln.reference_name, frag_aln_start, frag_aln_end, frag_len, sample])))))
        n_frags += 1
        if n_aln % 10000 == 0:
            print(f"Processed {n_aln} alignments")
            n_aln += 1

    colnames = [
        "seqnames",
        "start",
        "end",
        "length",
        "aln_type",
        "dup",
    ]
    frags_df = pd.DataFrame(frags, columns=colnames)
    if not mdup_or_not:
        frags_df.loc[:, "dup"] = -1

    return frags_df, n_frags


def collect_dna_metrics(args: argparse.Namespace) -> int:
    """Collect metrics for dna analysis."""
    # FIXME: Focus on metrics on autosomes, chrX, and chrY, not scaffold chromosomes
    # Currently, some reads to scaffold chromosomes are included (i.e. frags.csv output)
    out_dir = args.out_dir / args.sample
    out_dir.mkdir(exist_ok=True)

    frags_df, n_frags = collect_frags(args.bam_file)
    genome = pyfaidx.Fasta(str(args.fa))
    target_df = pr.read_bed(str(args.target), as_df=True)
    if "Name" not in target_df.columns:
        target_df.loc[:, "Name"] = (
            target_df.Chromosome + ":" + target_df.Start.astype(str) + "-" + target_df.End.astype(str)
        )
    target_df = target_df.apply(get_gc, axis=1, args=(genome,))
    target_pr = pr.PyRanges(target_df, int64=True)

    frags_df.columns = [
        "Chromosome",
        "Start",
        "End",
        "length",
        "aln_type",
        "dup",
    ]
    frags_df.loc[:, "frag_id"] = np.arange(0, frags_df.shape[0], 1)
    frags_pr = pr.PyRanges(frags_df, int64=True)  # , chromosomes=frags_df.seqnames,
    mdup_or_not = 0 if frags_df[frags_df.dup == -1].shape[0] == frags_df.shape[0] else 1

    ovl_pr = frags_pr.join(target_pr, how="left", report_overlap=True)
    ovl_df = ovl_pr.df
    ovl_df.loc[:, "on_target_bases"] = ovl_df.Overlap
    ovl_df.drop(["Overlap"], axis=1, inplace=True)
    ovl_df.loc[ovl_df.Start_b == -1, "Name"] = "."
    ovl_df.loc[:, "on_or_off"] = np.where((ovl_df.aln_type == 1) & (ovl_df.Start_b != -1), 1, 0)
    ovl_df.loc[ovl_df.on_or_off != 1, "on_target_bases"] = 0
    # if a fragment mapped to multiple targets,
    # the fragment is a non-proper fragment
    ovl_df.drop_duplicates(subset=["frag_id"], inplace=True)
    assert ovl_df.shape[0] == frags_df.shape[0]
    assert ovl_df[ovl_df.on_or_off == 0].shape[0] + ovl_df[ovl_df.on_or_off == 1].shape[0] == frags_df.shape[0]

    ovl_df.loc[:, "usable"] = ovl_df.on_or_off.map(lambda x: 1 if x == 1 else 0)
    cov_df = ovl_df.groupby(["Name"]).agg({"usable": "sum"}).reset_index()
    cov_df = pd.merge(cov_df, target_df, how="inner", on=["Name"], copy=False)
    cov_df.columns = [
        "amplicon",
        "coverage",
        "seqnames",
        "start",
        "end",
        "length",
        "gc",
    ]
    cov_df = cov_df[["seqnames", "start", "end", "amplicon", "length", "gc", "coverage"]]
    median_target_cov = np.median(cov_df.coverage)  # type: ignore  # incomplete numpy typing

    cov_df.loc[:, "norm_coverage"] = cov_df.coverage / median_target_cov if median_target_cov > 0 else 0.0
    chrs = [f"chr{i}" for i in range(23)] + ["chrX", "chrY"]
    cov_df.loc[:, "seqnames"] = pd.Categorical(cov_df.seqnames, chrs)
    cov_df.sort_values(by=["seqnames", "start"], inplace=True)
    # FIXME: Add amplicon dropout metrics
    # n_amplicon_dropout = cov_df[cov_df.coverage == 0].shape[0]
    pct_ave_ovl = round(np.mean(100 * ovl_df.on_target_bases / ovl_df.amplicon_length), 2)

    frag_cols = [
        "Chromosome",
        "Start",
        "End",
        "length",
        "aln_type",
        "Name",
        "on_or_off",
        "on_target_bases",
        "dup",
    ]
    ovl_df = ovl_df[frag_cols]
    ovl_df.columns = [
        "seqnames",
        "start",
        "end",
        "frag_len",
        "aln_type",
        "amplicon",
        "on_or_off",
        "on_target_bases",
        "dup_or_not",
    ]
    ovl_df.loc[:, "sample"] = args.sample

    # output fragments table
    out_frags = out_dir / f"{args.sample}.frags.csv"
    ovl_df.to_csv(out_frags, sep=",", index=False)

    # output coverage table
    out_cov = out_dir / f"{args.sample}.targets.tsv"
    cov_df.to_csv(out_cov, sep="\t", index=False)

    # output metrics table
    median_frag = np.median(ovl_df.frag_len)  # type: ignore  # incomplete numpy typing
    mad_frag = median_abs_deviation(ovl_df.frag_len)
    median_proper_frag = np.median(ovl_df[ovl_df.aln_type == 1].frag_len)  # type: ignore  # incomplete numpy typing
    mad_proper_frag = median_abs_deviation(ovl_df[ovl_df.aln_type == 1].frag_len)
    proper_pct = round(100 * ovl_df[ovl_df.aln_type == 1].shape[0] / n_frags, 2)
    nonproper_pct_same = round(100 * ovl_df[ovl_df.aln_type == 0].shape[0] / n_frags, 2)
    nonproper_pct_diff = round(100 * ovl_df[ovl_df.aln_type == -1].shape[0] / n_frags, 2)
    n_on_target_frags = ovl_df[ovl_df.on_or_off == 1].shape[0]
    pct_on_target_frags = round(100 * n_on_target_frags / ovl_df.shape[0], 2)
    tot_on_target_bases = sum(ovl_df[ovl_df.on_or_off == 1].on_target_bases)
    tot_bases = sum(ovl_df[ovl_df.aln_type == 1].frag_len)
    pct_on_target_bases = round(
        100 * tot_on_target_bases / tot_bases,
        2,
    )
    dup_rate = ovl_df[ovl_df.dup_or_not == 0].shape[0] / ovl_df.shape[0] if mdup_or_not else -1.0
    cov_quantile: Any = np.quantile(cov_df.coverage, q=[0.1, 0.9])  # type: ignore  # incomplete numpy typing
    q90_q10_ratio = cov_quantile[1] / cov_quantile[0] if median_target_cov > 0 else 99999
    pct_lt_5 = round(100 * cov_df[cov_df.norm_coverage < 0.5].shape[0] / cov_df.shape[0], 2)
    # norm_cov_quantile = np.quantile(cov_df.norm_coverage, q=[0.1, 0.9])
    # norm_q90_q10_ratio = norm_cov_quantile[1]/norm_cov_quantile[0]
    # print(n_frags, ovl_df.shape[0])
    # print(proper_pct, nonproper_pct_same, nonproper_pct_diff)
    # print(n_on_target_frags, pct_on_target_frags, pct_on_target_bases, dup_rate)
    # print(median_amplicon_cov, np.mean(cov_df.coverage), q90_q10_ratio)
    metrics_dict = {
        "sample": [args.sample],
        "n_frags": [n_frags],
        "median_frag": [median_frag],
        "mad_frag": [mad_frag],
        "median_proper_frag": [median_proper_frag],
        "mad_proper_frag": [mad_proper_frag],
        "pct_proper": [proper_pct],
        "pct_nonproper_same_chr": [nonproper_pct_same],
        "pct_nonproper_diff_chr": [nonproper_pct_diff],
        "n_on_target_frags": [n_on_target_frags],
        "pct_on_target_frags": [pct_on_target_frags],
        "pct_on_target_bases": [pct_on_target_bases],
        "pct_ave_overlap": [pct_ave_ovl],
        "dup_rate": [dup_rate],
        "median_target_cov": [median_target_cov],
        "q90_q10_ratio": [q90_q10_ratio],
        "pct_lt_5": [pct_lt_5],
    }
    metrics_df = pd.DataFrame(metrics_dict)
    out_metric = out_dir / f"{args.sample}.metrics.csv"
    metrics_df.to_csv(out_metric, sep=",", index=False)
    return 0


def collect_rna_metrics(args: argparse.Namespace) -> None:
    """Collect metrics for rna analysis."""
    out_dir = args.out_dir / args.sample
    out_dir.mkdir(parents=True, exist_ok=True)


def parse_cmd() -> argparse.Namespace:
    """Parse command arguments."""
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(title="Commands", metavar="", dest="command")
    dna = subparsers.add_parser("dna", description="Collect DNA metrics")
    dna.add_argument(
        "--bam",
        metavar="FILE",
        dest="bam_file",
        type=parse_path,
        required=True,
        help="specify the bam file from which metrics will be extracted [Required]",
    )
    dna.add_argument(
        "--sample",
        metavar="STR",
        dest="sample",
        required=True,
        help="specify the sample name [Required]",
    )
    dna.add_argument(
        "--out",
        metavar="DIR",
        dest="out_dir",
        type=parse_path,
        required=True,
        help="specify the output directory [Required]",
    )
    dna.add_argument(
        "--panel",
        metavar="FILE",
        dest="target",
        type=parse_path,
        required=True,
        help="specify the panel BED file [Required]",
    )
    dna.add_argument(
        "--genome",
        metavar="FILE",
        dest="fa",
        type=parse_path,
        required=True,
        help="specify the genome fasta [Required]",
    )
    dna.set_defaults(func=collect_dna_metrics)

    rna = subparsers.add_parser("rna", description="Collect RNA metrics")
    rna.add_argument(
        "--bam",
        metavar="FILE",
        dest="bam_file",
        type=parse_path,
        required=True,
        help="specify the bam file from which metrics will be extracted [Required]",
    )
    rna.add_argument(
        "--gtf",
        metavar="FILE",
        dest="gtf_file",
        type=parse_path,
        required=True,
        help="specify the annotation file in GTF [Required]",
    )
    rna.add_argument(
        "--sample",
        metavar="STR",
        dest="sample",
        required=True,
        help="specify the sample name [Required]",
    )
    rna.add_argument(
        "--out",
        metavar="DIR",
        dest="out_dir",
        type=parse_path,
        required=True,
        help="specify the output directory [Required]",
    )
    rna.add_argument(
        "--genome",
        metavar="FILE",
        dest="fa",
        type=parse_path,
        required=True,
        help="specify the genome fasta [Required]",
    )
    rna.set_defaults(func=collect_rna_metrics)

    return parser.parse_args()


def main() -> None:
    """Run collect metrics."""
    args = parse_cmd()
    args.func(args)


if __name__ == "__main__":
    main()
