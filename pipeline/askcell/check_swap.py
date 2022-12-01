"""check_swap.py to check sample swap, especially between tumor and matched normal samples."""

import argparse
import sys

import pandas as pd
import pysam


def check_swap(args: argparse.Namespace) -> None:
    """Check swap."""
    variants = pysam.VariantFile(args.vcf)

    tot = 0
    n_low_dp = 0
    ibs0, ibs2 = 0, 0
    hets1, hets2 = 0, 0
    hom_alts1, hom_alts2 = 0, 0
    shared_hets, shared_hom_alts = 0, 0
    samples = []

    for variant in variants:
        if len(variant.alts or []) > 1:
            continue

        tot += 1
        format_fields = [k for k in variant.format if k in ["GT", "DP", "AD"]]
        if len(format_fields) != 3:
            print("[ERROR]: GT, DP, and AD are required fields for the FORMAT")
            sys.exit(1)

        sample1, sample2 = variant.samples.keys()

        if sample1 not in samples:
            samples.append(sample1)
        if sample2 not in samples:
            samples.append(sample2)

        dp1 = variant.samples[sample1]["DP"]
        dp2 = variant.samples[sample2]["DP"]
        if dp1 > args.min_dp and dp2 > args.min_dp:
            ad1 = variant.samples[sample1]["AD"][1]
            # ad2 = variant.samples[sample2]["AD"][1]
            af1 = round(ad1 / dp1, 2)
            # af2 = round(ad2 / dp2, 2)
            gt1 = variant.samples[sample1]["GT"]
            gt2 = variant.samples[sample2]["GT"]

            # ugly looking code but it works for now
            gt1 = 0 if gt1 == (0, 0) else (1 if gt1 == (0, 1) else 2)
            gt2 = 0 if gt2 == (0, 0) else (1 if gt2 == (0, 1) else 2)

            if gt1 == 1:
                # should use common SNPs, but the overlaps too few
                # perhaps check gnomad
                if args.min_het_af <= af1 <= args.max_het_af:
                    gt1 = 1
                else:
                    gt1 = -1

            if gt2 == 1:
                if args.min_het_af <= af1 <= args.max_het_af:
                    gt2 = 1
                else:
                    gt2 = -1

            if gt1 == gt2:
                ibs2 += 1
                if gt1 == 1:
                    shared_hets += 1
                    hets1 += 1
                    hets2 += 1
                elif gt1 == 2:
                    hom_alts1 += 1
                    hom_alts2 += 1
                    shared_hom_alts += 1
                else:
                    pass
            else:
                if (gt1 == 0 and gt2 == 2) or (gt1 == 2 and gt2 == 0):
                    ibs0 += 1
                if gt1 == 1:
                    hets1 += 1
                if gt2 == 1:
                    hets2 += 1
                if gt1 == 2:
                    hom_alts1 += 1
                if gt2 == 2:
                    hom_alts2 += 1
        else:
            n_low_dp += 1

    if min(hets1, hets2) > 0:
        relateness = (shared_hets - 2 * ibs0) / min(hets1, hets2)
    else:
        relateness = -1

    df = pd.DataFrame(
        {
            "Sample1": [samples[0]],
            "Sample2": [samples[1]],
            "TotalSitesChecked": [tot],
            "NumLowDP": [n_low_dp],
            "Hets1": [hets1],
            "HomAlts1": [hom_alts1],
            "Hets2": [hets2],
            "HomAlts2": [hom_alts2],
            "IBS2": [ibs2],
            "IBS0": [ibs0],
            "NumSharedHets": [shared_hets],
            "NumSharedHomAlts": [shared_hom_alts],
            "Relatedness": [relateness],
        },
    )
    df.to_csv(args.out, sep="\t", index=False)


def parse_cmd() -> argparse.Namespace:
    """Parse command arguments."""
    parser = argparse.ArgumentParser(description="Check sample swaps")
    parser.add_argument(
        "--vcf",
        metavar="FILE",
        dest="vcf",
        required=True,
        help="Specify the input vcf file that includes the two samples to check",
    )
    parser.add_argument(
        "--out",
        metavar="FILE",
        dest="out",
        required=True,
        help="Specify the output file",
    )
    parser.add_argument(
        "--min_dp",
        metavar="INT",
        dest="min_dp",
        type=int,
        default=100,
        help="Specify the minimum DP at a site to consider [100]",
    )
    parser.add_argument(
        "--min_het_af",
        metavar="FLOAT",
        dest="min_het_af",
        type=float,
        default=0.35,
        help="Specify the minimum AF to consider HET genotype [0.35]",
    )
    parser.add_argument(
        "--max_het_af",
        metavar="FLOAT",
        dest="max_het_af",
        type=float,
        default=0.65,
        help="Specify the maximum AF to consider HET genotype [0.65]",
    )
    return parser.parse_args()


def main() -> None:
    """Check sample swap."""
    args: argparse.Namespace = parse_cmd()

    check_swap(args)


if __name__ == "__main__":
    main()
