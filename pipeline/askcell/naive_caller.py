"""naive_caller.py."""

import argparse

import pandas as pd

from .pileup2count import run_pileup2counts
from .utils.pathutils import parse_path


def identify_alt_from_cnt(row: pd.DataFrame) -> pd.DataFrame:
    """Identify alt calls based on observed counts."""
    ref = row["REF"]
    dp = row["DEPTH"]
    rd = row["R"] + row["r"]
    bases = ["A", "C", "G", "T"]
    ads, afs = [], []
    alts = []
    for alt in bases:
        if alt == ref:
            continue
        else:
            ad = row[alt] + row[alt.lower()]
            af = ad / (dp * 1.0) if (dp > 0) else 0.0
            if ad > 1:
                alts.append(alt)
                ads.append(ad)
                afs.append(af)
    if alts:
        row["ALT"] = ",".join(alts)
        row["AD"] = ",".join(map(str, ads))
        row["AF"] = ",".join(map(str, afs))
    else:
        row["ALT"] = "."
        row["AD"] = 0
        row["AF"] = 0.0
    row["RD"] = rd
    return row


def call_variants_naively(args: argparse.Namespace) -> pd.DataFrame:
    """Count each base at each site from a pileup."""
    print("[Info] Start counting [START]")
    # set up output dir
    outdir = args.out.parent
    outdir.mkdir(parents=True, exist_ok=True)

    cnt_df = pd.DataFrame()
    if not args.out.exists():
        cnt_df = run_pileup2counts(args)
    else:
        cnt_df = pd.read_csv(args.out, sep=",")
        print(f"[Info] previous count table exists: {args.out}")
        print("[Info] Skip counting")
    print("[Info] Finished counting [DONE]")

    cnt_df = cnt_df.apply(identify_alt_from_cnt, axis=1)
    cnt_df = cnt_df[["CHR", "POSITION", "REF", "ALT", "DEPTH", "RD", "AD", "AF"]]
    cnt_df.loc[:, "Sample"] = args.sample
    out_snp = parse_path(str(args.out).split(".cnt.csv")[0]).with_suffix(".snp.tsv")
    cnt_df.to_csv(out_snp, sep="\t", index=False)


def parse_cmd() -> argparse.Namespace:
    """Command line interface."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bam",
        metavar="FILE",
        dest="bam",
        type=parse_path,
        required=True,
        help="specify the bam file used to generated pileups",
    )
    parser.add_argument(
        "--ref",
        metavar="FILE",
        dest="ref",
        type=parse_path,
        required=True,
        help="specify the human reference sequence",
    )
    parser.add_argument(
        "--target",
        metavar="FILE",
        dest="targets",
        type=parse_path,
        required=True,
        help="specify the targeted regions (BED) or positions (VCF)",
    )
    parser.add_argument(
        "--out",
        metavar="FILE",
        dest="out",
        type=parse_path,
        required=True,
        help="specify the output file, where prefix of this file is used for additional outputs",
    )
    parser.add_argument(
        "--sample",
        metavar="STR",
        dest="sample",
        required=True,
        help="specify the sample name",
    )
    parser.add_argument(
        "--pileup",
        dest="pileup_or_not",
        action="store_true",
        help="specify whether to output the final pileup file",
    )
    parser.add_argument(
        "--params",
        metavar="STR",
        dest="params",
        default="-aB -d99999 -q30 -Q25",
        help='specify the paramters to run samtools mpileup ["-a"]',
    )
    parser.add_argument(
        "--nproc",
        metavar="INT",
        dest="nproc",
        type=int,
        default=10,
        help="specify the number of processes to count alleles in parallel [10]",
    )
    parser.add_argument(
        "--batch_size",
        metavar="INT",
        dest="ntargets_per_proc",
        type=int,
        default=500,
        help="specify the number of regions/positions to process in parallel [500]",
    )
    return parser.parse_args()


def main() -> None:
    """Run naive caller."""
    args = parse_cmd()

    call_variants_naively(args)


if __name__ == "__main__":
    main()
