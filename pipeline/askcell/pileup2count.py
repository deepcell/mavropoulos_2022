"""pileup2count.py for generating the base counts from pileup."""

import argparse
import glob
import multiprocessing as mp
import os
import shlex
import signal
import subprocess as sp
import sys
import traceback

from collections import defaultdict

import pandas as pd
import pybedtools as pbt

from .utils.pathutils import PathLike, parse_path


NUCBASES = ["R", "r", "A", "a", "C", "c", "T", "t", "G", "g"]
CHRS = [f"chr{k + 1}" for k in range(22)] + ["chrX"]
PILECOUNTER_HEADER = ["CHR", "POSITION", "REF", "DEPTH"] + NUCBASES


class ReadPile:
    """ReadPile object."""

    def __init__(self) -> None:
        """Initialize the ReadPile object."""
        self.chrom = ""
        self.pos = ""
        self.ref = ""
        self.piles: list[tuple[int, str]] = []

    def parse(self, read_pile: str) -> None:
        """Parse the chromosome, position, and reference base information."""
        fields = read_pile.strip().split()
        self.chrom, self.pos, self.ref = fields[0:3]
        self.piles = [(int(fields[i]), fields[i + 1]) for i in range(3, len(fields), 3)]


class PileCounter:
    """PileCounter object."""

    def __init__(self, nsamples: int) -> None:
        """Initialize the Pilecounter object."""
        self.positions: list[str] = []
        self.nsamples = nsamples
        self.alleles = NUCBASES
        self.df = pd.DataFrame([], columns=PILECOUNTER_HEADER)
        self.counts: dict[int, dict[str, list[int]]] = defaultdict(lambda: defaultdict(list))
        self.indels: dict[int, dict[str, list[tuple[str, int]]]] = defaultdict(lambda: defaultdict(list))
        self.dp: list[int] = []
        self.zero_dp_sites: int = 0

    def count(self, pileup: str, skip_indels: bool = True) -> None:
        """Count the base coverage from the pileup data."""
        total_depth = 0  # total counts at a site across all samples
        readpile = ReadPile()
        readpile.parse(pileup)
        self.positions += ["-".join([readpile.chrom, readpile.pos, readpile.ref])]
        base_dict: dict[str, int] = defaultdict(int)
        for i in range(len(readpile.piles)):
            count, covering_bases = readpile.piles[i]
            total_depth += count
            j = 0
            while j < len(covering_bases):
                covering_base = covering_bases[j]
                if covering_base == ".":
                    base_dict["R"] += 1
                elif covering_base == ",":
                    base_dict["r"] += 1
                elif covering_base in "AaCcGgTt":
                    base_dict[covering_base] += 1
                # mark the base is at the beginning of a read # move one more index downstream
                elif covering_base == "^":
                    j += 1
                # mark the base is at the end of a read, skip
                elif covering_base in "$*":
                    pass
                # handle deletions
                elif covering_base == "-":
                    shift, del_seq = self.parse_indels(covering_bases[j + 1 :])
                    base_dict["-" + del_seq] += 1
                    j += shift
                # handle insertions
                elif covering_base == "+":
                    shift, ins_seq = self.parse_indels(covering_bases[j + 1 :])
                    base_dict["+" + ins_seq] += 1
                    j += shift
                else:
                    # should be all the cases
                    pass
                j += 1

            self.dp += [total_depth]
            for allele in self.alleles:
                self.counts[i + 1][allele] += [base_dict[allele]]
            if not skip_indels:
                position = f"{readpile.chrom}-{readpile.pos}"
                for allele, count in base_dict.items():
                    if allele.startswith("+") or allele.startswith("-"):
                        self.indels[i + 1][position] += [(allele, count)]

    def parse_indels(self, sub_read_pile: str) -> tuple[int, str]:
        """Parse indels.

        sub_read_pile start with the indel size number
        Thus, first need to get the size of the indel events

        """
        indel = ""
        for i in range(len(sub_read_pile)):
            if ord(sub_read_pile[i]) > 57:
                break
            # 0-9 ASCII value from 48 to 57 inclusive
            # avoid REGEX
            if ord(sub_read_pile[i]) >= 48 and ord(sub_read_pile[i]) <= 57:
                indel += sub_read_pile[i]
        shift = len(indel)
        indel_size = int(indel)
        indel_seq = sub_read_pile[i : i + indel_size]
        return (shift + indel_size, indel_seq)

    def check(self) -> int:
        """Check lengths."""
        for i in range(self.nsamples):
            all_allele_counts = [self.counts[i + 1][allele] for allele in self.alleles]
            if any(len(allele_counts) != len(self.positions) for allele_counts in all_allele_counts):
                print("[Error] Length of positions are not equal to length of allele counts")
                return 1
        return 0

    def tablize(self, outfile: PathLike) -> int:
        """Convert count dictionary to tables."""
        df = pd.DataFrame()
        for _, counts in self.counts.items():
            for allele in self.alleles:
                count_vec = pd.Series(counts[allele])
                df = pd.concat([df, count_vec], axis=1)
        if df.shape[0] == 0:
            print("Error: empty counts table")
            return 1
        first_3_cols = pd.DataFrame(self.positions, columns=["positions"])
        first_3_cols = first_3_cols["positions"].str.split("-", expand=True)
        self.df = pd.concat([first_3_cols, pd.Series(self.dp), df], axis=1)
        self.df.columns = ["CHR", "POSITION", "REF", "DEPTH"] + self.alleles
        self.zero_dp_sites = self.df[self.df.DEPTH == 0].shape[0]

        return 0


def assign_batch_jobs(
    targets: pbt.bedtool.BedTool,
    ntargets_per_proc: int,
    outp: PathLike,
) -> list[tuple[str, int, int, PathLike]]:
    """Assign batch jobs."""
    print("[Info] Assigning jobs")
    assignments: list[tuple[str, int, int, PathLike]] = []
    for chrom in CHRS:
        targets_per_chr = [k for k in targets if k.chrom == str(chrom)]
        assignment: list[str] = []
        chunk_start, chunk_stop = -1, -1
        if len(targets_per_chr) < ntargets_per_proc and targets_per_chr:
            chunk_start = targets_per_chr[0].start - 5000 if (targets_per_chr[0].start - 5000 >= 0) else 0
            chunk_stop = targets_per_chr[-1].stop + 5000
            assignment = [
                "\t".join(map(str, [target.chrom, target.start - 1, target.stop])) for target in targets_per_chr
            ]
            assigned_target_file = f"{outp}_{len(assignments)}.bed"
            with open(assigned_target_file, "w") as out:
                out.write("\n".join(assignment) + "\n")
            assignments.append((chrom, chunk_start, chunk_stop, assigned_target_file))
            continue
        for target in targets_per_chr:
            if chunk_start == -1:
                chunk_start = target.start - 5000 if (target.start - 5000 >= 0) else 0
            assignment.append("\t".join(map(str, [target.chrom, target.start - 1, target.stop])))
            if len(assignment) == ntargets_per_proc:
                # write to file and clear
                chunk_stop = target.stop + 5000
                assigned_target_file = f"{outp}_{len(assignments)}.bed"
                with open(assigned_target_file, "w") as out:
                    out.write("\n".join(assignment) + "\n")
                assignments.append((chrom, chunk_start, chunk_stop, assigned_target_file))
                assignment = []
                chunk_start, chunk_stop = -1, -1
        if assignment:
            chunk_stop = int(assignment[-1].split("\t")[-1]) + 5000
            assigned_target_file = f"{outp}_{len(assignments)}.bed"
            with open(assigned_target_file, "w") as out:
                out.write("\n".join(assignment) + "\n")
            assignments.append((chrom, chunk_start, chunk_stop, assigned_target_file))
    print(f"[Info] Assigned {len(assignments)} jobs")
    return assignments


def assign_batch_var_jobs(targets: pbt.bedtool.BedTool, outp: PathLike) -> list[tuple[str, int, int, PathLike]]:
    """Assign batch variant jobs."""
    assignments: list[tuple[str, int, int, PathLike]] = []
    print("pileup: targets")
    print(type(targets))
    for chrom in CHRS:
        vars_per_chr = [k for k in targets if k.chrom == str(chrom)]
        for var in vars_per_chr:
            chunk_start = var.start - 1000 if (var.start - 1000 >= 0) else 0
            chunk_stop = var.stop + 1000
            assignment = "\t".join(map(str, [var.chrom, var.start - 1, var.stop]))
            assigned_var_file = f"{outp}_{len(assignments)}.bed"
            with open(assigned_var_file, "w") as out:
                out.write(assignment + "\n")
            assignments.append((chrom, chunk_start, chunk_stop, assigned_var_file))
    return assignments


def count_worker(
    nth_job: int,
    assignment: tuple[str, int, int, PathLike],
    ref: PathLike,
    params: str,
    bams: PathLike,
    outp: PathLike,
    pileup_or_not: bool,
) -> tuple[int, str, pd.DataFrame, int]:
    """Work on Counting.

    The pileup2count worker performs the following tasks:

        1. generate pileup at each site
        2. count bases at each site
        3. tablize counts into table

    """
    nsamples = len(str(bams).split())
    out_pileup_file = f"{outp}_{nth_job}.pileup"
    out_tab_file = f"{outp}_{nth_job}.counts.h5"
    chrom, chunk_start, chunk_stop, assigned_target = assignment
    pilecounter = PileCounter(nsamples)
    assigned_target_file = parse_path(assigned_target)

    try:
        # if(pileup_or_not):
        #    fPILEUP = open(out_pileup_file, 'w', 1)
        mpileup_cmd = (
            f'samtools mpileup {params} -r "{chrom}:{chunk_start}-{chunk_stop}" '
            f"-l {assigned_target_file} -f {ref} {bams}"
        )
        print(f"[Info] job={nth_job}\tchrom={chrom}\tchunk_start={chunk_start}\tchunk_stop={chunk_stop} [START]")

        p = sp.Popen(shlex.split(mpileup_cmd), stdout=sp.PIPE, stderr=sp.PIPE, encoding="utf-8")

        assert p.stdout is not None and p.stderr is not None  # makes mypy happy

        for line in p.stdout:
            pilecounter.count(line)

        rc = p.returncode
        if rc:
            print(f"[Error] worker {nth_job} FAIL at count_worker (mpileup + pilecounter.count")
            traceback.print_exc()
            raise Exception(p.stderr.readline())

        if pilecounter.positions:
            rc = pilecounter.check()
            if rc:
                raise Exception(
                    f"[Error] worker {nth_job} reports an error when length of positions "
                    "are not equal to length of allele counts",
                )

            rc = pilecounter.tablize(out_tab_file)
            if rc:
                raise Exception(
                    f"[Error] worker {nth_job} reports an error " "when trying to tablize all counts",
                )

        assigned_target_file.unlink(missing_ok=True)

    except Exception as e:  # noqa: B902
        traceback.print_exc()
        print(
            f"[Info] job={nth_job}\tchrom={chrom}\tchunk_start={chunk_start}\t"
            f"chunk_stop={chunk_stop}\trcode={rc} [ERROR]",
        )
        raise Exception(
            f"[Error] worker {nth_job} reports an error when trying to count alleles from samtools mpileup",
        ) from e

    else:
        print(
            f"[Info] job={nth_job}\tchrom={chrom}\tchunk_start={chunk_start}\t"
            f"chunk_stop={chunk_stop}\trcode={rc} [DONE]",
        )
        return nth_job, out_pileup_file, pilecounter.df, pilecounter.zero_dp_sites


def run_pileup2counts(args: argparse.Namespace, skip_indels: bool = True) -> pd.DataFrame:
    """Run main gatekeeper for calling mpileup and generating counts from pileups."""
    print("[Info] pileupcount [START]")
    # set up parallel processes
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    pool = mp.Pool(processes=args.nproc)
    signal.signal(signal.SIGINT, original_sigint_handler)

    # if there are undeleted files from a previously interrupted run
    # then clean them first
    intermediates = glob.glob(os.path.join(os.path.dirname(args.out), os.path.basename(args.out) + "_*bed"))
    pool.map(os.remove, intermediates)

    # read targets and dispatch jobs in batches
    if args.targets.name.endswith("vcf"):
        target_beds = vcf2bed(args.targets)
        assignments = assign_batch_var_jobs(target_beds, args.out)
    elif args.targets.name.endswith("bed"):
        print("[Info] Reading target bed file")
        target_beds = pbt.BedTool(args.targets)
        target_beds = [(k.chrom, k.start + 1, k.stop) for k in target_beds.merge()]
        target_beds = pbt.BedTool(target_beds)
        print("[Info] Reading target bed file [FINISHED]")
        assignments = assign_batch_jobs(target_beds, args.ntargets_per_proc, args.out)
    else:
        print("not supported target file suffix")
        print("supported: bed, vcf")

    # set up parallel processes
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    pool = mp.Pool(processes=args.nproc)
    signal.signal(signal.SIGINT, original_sigint_handler)
    output = None
    output = pool.starmap(
        count_worker,
        [
            (
                i,
                assignments[i],
                args.ref,
                args.params,
                args.bam,
                args.out,
                args.pileup_or_not,
            )
            for i in range(len(assignments))
        ],
        chunksize=1,
    )

    # aggregate individual results
    dfs = []

    n_tot_zero_dp_sites = 0
    # for (nth_job, out_pileup_file, out_tab_file, n_zero_dp_sites) in output:
    for (_nth_job, _out_pileup_file, out_df, n_zero_dp_sites) in output:
        n_tot_zero_dp_sites += n_zero_dp_sites
        dfs.append(out_df)
    print(f"[Info] {len(assignments)} jobs returned")
    print(f"[Info] {len(dfs)} in-memory dataframes concated")
    if len(dfs) != len(assignments):
        sys.exit()
    pool.close()
    pool.join()

    merged_df = pd.concat(dfs, axis=0)
    merged_df.loc[:, "Sample"] = args.sample
    print(f"[Info] Final count table includes {merged_df.shape[0]} sites")
    print(f"[Info] Final count table includes {n_tot_zero_dp_sites} sites with zero DP")
    print("[Info] save count tables")
    merged_df.to_csv(args.out, sep=",", index=False)
    print("[Info] pileupcount [DONE]")
    return merged_df


def vcf2bed(vcf_file: PathLike) -> pbt.bedtool.BedTool:
    """Extract variant position in a VCF and convert to a BED file (only the positions)."""
    targets = []
    with open(vcf_file) as fIN:
        for line in fIN:
            if line.startswith("#"):
                continue
            tmpLine = line.strip().split()
            chrom = tmpLine[0]
            ref = tmpLine[3]
            alts = tmpLine[4].split(",")
            if len(alts) == 1 and len(alts[0]) != len(ref):
                # indel
                continue
            if len(alts) == 1 and len(ref) > 1 and len(ref) == len(alts[0]):
                # mnp
                pos = int(tmpLine[1])
                alts = list(alts[0])
                targets += [(chrom, p, p) for p in range(pos, pos + len(alts))]
                continue
            pos = int(tmpLine[1])
            targets.append((chrom, pos, pos))

    return pbt.BedTool(targets)
