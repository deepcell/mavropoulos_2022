"""askcell pipeline for DNA and RNA analysis."""

from __future__ import annotations

import argparse

from .dna.amplicon import (
    run_sid_sample,
    run_sid_samples,
    run_tgs_sample,
    run_tgs_samples,
)
from .dna.exome import run_exome_sample, run_exome_samples
from .dna.wgs import run_wgs_sample, run_wgs_samples
from .rna.rnaseq import run_rnaseq_sample, run_rnaseq_samples
from .static_data import P5, P7
from .utils.config import add_config_args
from .utils.pathutils import parse_path


def arg_parser() -> argparse.ArgumentParser:
    """Create an argument parser for AskCell pipeline."""
    # Top-level parser
    parser = argparse.ArgumentParser()
    # Options to override
    add_config_args(parser)

    commands = parser.add_subparsers(title="Commands", dest="command")

    parser.add_argument("--dryrun", action="store_true", help="print commands without execution")
    parser.add_argument("--threads", type=int, default=8, help="number of threads used for the run")
    parser.add_argument(
        "--logs-dir",
        type=parse_path,
        help="defines logs directory, default is <output_dir> / logs",
    )

    # TODO: Will refactor into DNA and RND workflows in future PR (DNA and RNA subpackages)
    # TODO: Add demultiplex entry points
    cmd = commands.add_parser("sid-sample", help="a single sample trim/alignment")
    cmd.add_argument("input_dir", metavar="DIR", type=parse_path, help="sequencing run folder w/ fastqs")
    cmd.add_argument("--output-dir", "-o", metavar="DIR", required=True, type=parse_path, help="output directory")
    cmd.add_argument("--scratch-dir", metavar="DIR", type=parse_path, help="scratch directory")
    cmd.add_argument("--sampleid", metavar="STR", required=True, help="sample identifier")
    cmd.add_argument("--panel", metavar="STR", help="swift panel name")
    cmd.add_argument(
        "--adapter1",
        metavar="STR",
        type=str,
        default=P5,
        help="sequence used to trim adapter sequences",
    )
    cmd.add_argument(
        "--adapter2",
        metavar="STR",
        type=str,
        default=P7,
        help="sequence used to trim adapter sequences",
    )
    cmd.set_defaults(func=run_sid_sample)

    cmd = commands.add_parser("sid", help="amplicon-based genotype panel workflow")
    cmd.add_argument("samplesheet", metavar="FILE", type=parse_path, help="samplesheet")
    cmd.add_argument("input_dir", metavar="DIR", type=parse_path, help="sequencing run folder w/ fastqs")
    cmd.add_argument("--output-dir", "-o", metavar="DIR", required=True, type=parse_path, help="output directory")
    cmd.add_argument("--scratch-dir", metavar="DIR", type=parse_path, help="scratch directory")
    cmd.add_argument(
        "--adapter1",
        metavar="STR",
        type=str,
        default=P5,
        help="sequence used to trim adapter sequences",
    )
    cmd.add_argument(
        "--adapter2",
        metavar="STR",
        type=str,
        default=P7,
        help="sequence used to trim adapter sequences",
    )
    cmd.set_defaults(func=run_sid_samples)

    cmd = commands.add_parser("tgs-sample", help="a single sample trim/alignment")
    cmd.add_argument("input_dir", metavar="DIR", type=parse_path, help="sequencing run folder w/ fastqs")
    cmd.add_argument("matched_dir", metavar="DIR", type=parse_path, help="sequencing run folder w/ fastqs")
    cmd.add_argument("--output-dir", "-o", metavar="DIR", required=True, type=parse_path, help="output directory")
    cmd.add_argument("--scratch-dir", metavar="DIR", type=parse_path, help="scratch directory")
    cmd.add_argument("--sampleid", metavar="STR", required=True, help="sample identifier")
    cmd.add_argument(
        "--matchedid",
        metavar="STR",
        required=True,
        type=str,
        help="sample identifier for matched normal",
    )
    cmd.add_argument("--panel", metavar="STR", help="swift panel name")
    cmd.add_argument(
        "--adapter1",
        metavar="STR",
        type=str,
        default=P5,
        help="sequence used to trim adapter sequences",
    )
    cmd.add_argument(
        "--adapter2",
        metavar="STR",
        type=str,
        default=P7,
        help="sequence used to trim adapter sequences",
    )
    cmd.set_defaults(func=run_tgs_sample)

    cmd = commands.add_parser("tgs", help="amplicon-based tumor panel workflow")
    cmd.add_argument("samplesheet", metavar="FILE", type=parse_path, help="samplesheet")
    cmd.add_argument("input_dir", metavar="DIR", type=parse_path, help="sequencing run folder w/ fastqs")
    cmd.add_argument("--output-dir", "-o", metavar="DIR", required=True, type=parse_path, help="output directory")
    cmd.add_argument("--scratch-dir", metavar="DIR", type=parse_path, help="scratch directory")
    cmd.add_argument(
        "--adapter1",
        metavar="STR",
        type=str,
        default=P5,
        help="sequence used to trim adapter sequences",
    )
    cmd.add_argument(
        "--adapter2",
        metavar="STR",
        type=str,
        default=P7,
        help="sequence used to trim adapter sequences",
    )
    cmd.set_defaults(func=run_tgs_samples)

    cmd = commands.add_parser("rnaseq-sample", help="a single sample rnaseq data processing")
    cmd.add_argument("metadata", metavar="FILE", type=parse_path, help="metadata")
    cmd.add_argument("input_dir", metavar="DIR", type=parse_path, help="sequencing run folder w/ fastqs")
    cmd.add_argument("--output-dir", "-o", metavar="DIR", required=True, type=parse_path, help="output directory")
    cmd.add_argument("--scratch-dir", metavar="DIR", type=parse_path, help="scratch directory")
    cmd.add_argument("--sampleid", metavar="STR", required=True, help="sample identifier")
    cmd.add_argument("--species", metavar="STR", default="human", help="human or mouse")
    cmd.add_argument(
        "--adapter1",
        metavar="STR",
        type=str,
        default=P5,
        help="sequence used to trim adapter sequences",
    )
    cmd.add_argument(
        "--adapter2",
        metavar="STR",
        type=str,
        help="sequence used to trim adapter sequences",
    )
    cmd.set_defaults(func=run_rnaseq_sample)

    cmd = commands.add_parser("rnaseq", help="RNA-seq workflow")
    cmd.add_argument("samplesheet", metavar="FILE", type=parse_path, help="samplesheet")
    cmd.add_argument("metadata", metavar="FILE", type=parse_path, help="metadata")
    cmd.add_argument("input_dir", metavar="DIR", type=parse_path, help="sequencing run folder w/ fastqs")
    cmd.add_argument("--output-dir", "-o", metavar="DIR", required=True, type=parse_path, help="output directory")
    cmd.add_argument("--scratch-dir", metavar="DIR", type=parse_path, help="scratch directory")
    cmd.add_argument("--species", metavar="STR", default="human", help="human or mouse")
    cmd.add_argument(
        "--adapter1",
        metavar="STR",
        type=str,
        default=P5,
        help="sequence used to trim adapter sequences",
    )
    cmd.add_argument(
        "--adapter2",
        metavar="STR",
        type=str,
        help="sequence used to trim adapter sequences",
    )
    cmd.set_defaults(func=run_rnaseq_samples)

    # TODO: Add CNV workflow

    cmd = commands.add_parser("wes-sample", help="a single sample whole exome sequencing data processing")
    cmd.add_argument("input_dir", metavar="DIR", type=parse_path, help="sequencing run folder w/ fastqs")
    cmd.add_argument("--output-dir", "-o", metavar="DIR", required=True, type=parse_path, help="output directory")
    cmd.add_argument("--scratch-dir", metavar="DIR", type=parse_path, help="scratch directory")
    cmd.add_argument("--sampleid", metavar="STR", required=True, help="sample identifier")
    cmd.add_argument(
        "--control-name",
        metavar="STR",
        help="if control, specify reference name such as NA12878 or NA12877",
    )
    cmd.set_defaults(func=run_exome_sample)

    cmd = commands.add_parser("wes", help="whole exome sequencing workflow")
    cmd.add_argument("samplesheet", metavar="FILE", type=parse_path, help="samplesheet")
    cmd.add_argument("input_dir", metavar="DIR", type=parse_path, help="sequencing run folder w/ fastqs")
    cmd.add_argument("--output-dir", "-o", metavar="DIR", required=True, type=parse_path, help="output directory")
    cmd.add_argument("--scratch-dir", metavar="DIR", type=parse_path, help="scratch directory")
    cmd.set_defaults(func=run_exome_samples)

    cmd = commands.add_parser("wgs-sample", help="run single whole genome sequencing sample")
    cmd.add_argument("input_dir", metavar="DIR", type=parse_path, help="sequencing run folder w/ fastqs")
    cmd.add_argument("--output-dir", "-o", metavar="DIR", required=True, type=parse_path, help="output directory")
    cmd.add_argument("--scratch-dir", metavar="DIR", type=parse_path, help="scratch directory")
    cmd.add_argument("--sampleid", metavar="STR", required=True, help="sample identifier")
    cmd.add_argument("--runid", metavar="STR", required=True, help="run identifier")
    # FIXME: limited choices for bsize once verification / validation done to identify reasonable values
    # FIXME: Also needs proper background files
    cmd.set_defaults(func=run_wgs_sample)

    cmd = commands.add_parser("wgs", help="whole genome sequencing workflow")
    cmd.add_argument("samplesheet", metavar="metasheet", type=parse_path, help="samplesheet")
    cmd.add_argument("input_dir", metavar="input", type=parse_path, help="sequencing run folder w/ fastqs")
    cmd.add_argument("--output-dir", "-o", metavar="DIR", required=True, type=parse_path, help="output directory")
    cmd.add_argument("--scratch-dir", metavar="DIR", type=parse_path, help="scratch directory")
    # FIXME: move some arguments to config file once the config functionality is available
    cmd.add_argument(
        "--run-samples",
        metavar="samples in metasheet to run",
        nargs="*",
        type=str,
        help="a subsert of samples that should be run",
    )
    cmd.set_defaults(func=run_wgs_samples)

    return parser


def main(
    *,
    argv: list[str] | None = None,
    args: argparse.Namespace | None = None,
) -> None:
    """Askcell CLI entry point.

    Args:
        argv: literal command line arguments
        args: parsed arguments

    """
    parser = arg_parser()

    if argv is not None and args is not None:
        raise ValueError("argv and args are mutually exclusive")
    elif args is None:
        args = parser.parse_args(argv)

    if args.command and args.func:
        args.func(args)
    else:
        parser.print_help()
