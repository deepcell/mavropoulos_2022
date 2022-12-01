"""scRNA-seq analysis workflow."""

from __future__ import annotations

import argparse
import logging
import shutil
import tempfile

from typing import Any, Callable

import pandas as pd

from .. import Path, only_one
from ..logging import add_logging_args, init_logging  # FIXME: use dcutil logger
from ..utils.config import PackageConfig, add_config_args
from ..utils.pathutils import parse_path  # FIXME: use dcutil pathutils
from ..utils.shellutils import Command
from ..utils.types import ImmutableStrMapping
from ..workflow import (
    StepHelper,
    bin_files_directories_from_config,
    check_exist,
    get_scrna_settings,
    load_config,
    make_dirs,
    stage_data,
    submit_slurm,
    write_data,
)
from .cell_analysis import run_scanpy_analysis
from .demultiplex import run_cellranger_mkfastq, run_sequencing_qc
from .request import AnalysisRequest
from .sample_analysis import make_sample_metadata, run_cellranger_count, run_qc


STEP_FUNC_MAP: dict[str, Callable[..., Any]] = {
    "demultiplex": run_cellranger_mkfastq,
    "demultiplex_qc": run_sequencing_qc,
    "10x": run_cellranger_count,
    "qc": run_qc,
    "scanpy": run_scanpy_analysis,
}


def run_workflow(args: argparse.Namespace) -> None:
    """Run scRNA-seq analysis workflow."""
    if args.command in ("demux", "count", "analyze"):
        run_task(args, task_name=args.command)

    else:
        request = AnalysisRequest(args.request)
        request.load_sample_metadata(args.metadata)

        demux_output_dir = args.output_dir / "reads"
        results_output_dir = args.output_dir / "counts"
        analysis_output_dir = args.output_dir / "analysis"
        args.slurm_logs_dir = args.slurm_logs_dir or args.output_dir / "logs"

        make_dirs(dirs=[args.output_dir, args.slurm_logs_dir])

        cmd = Command(["scrna"])
        if args.keep_temp:
            cmd += ["--keep-temp", args.keep_temp]

        if args.local_temp_dir is not None:
            cmd += ["--local-temp-dir", args.local_temp_dir]

        if args.shared_temp_dir is not None:
            cmd += ["--shared-temp-dir", args.shared_temp_dir]

        if args.slurm_logs_dir is not None:
            cmd += ["--slurm-logs-dir", args.slurm_logs_dir]

        if args.set_config is not None:
            for s in args.set_config:
                cmd += ["-s", s]

        if args.config is not None:
            for s in args.config:
                cmd += ["-c", s]

        demux_cmd = cmd + [
            "demux",
            "--run-dir",
            args.input_dir,
            "--output-dir",
            demux_output_dir,
        ]

        if args.scratch_dir is not None:
            demux_cmd += ["--scratch-dir", args.scratch_dir]

        demux_cmd += [args.request, args.metadata]
        demux_jobid = submit_slurm(demux_cmd, threads=args.threads, slurm_log=args.slurm_logs_dir / "demux.log")

        sample_cmds = []
        sample_jobids = []
        for sid in request.samples:
            sample_cmd = cmd + [
                "count",
                "--reads-dir",
                demux_output_dir / f"{request.runid}/fastq",
                "--output-dir",
                results_output_dir / sid.sid,
                "--sid",
                sid.sid,
            ]

            if args.scratch_dir is not None:
                sample_cmd += ["--scratch-dir", args.scratch_dir]

            sample_cmd += [args.request, args.metadata]
            sample_cmds.append(sample_cmd)

            sample_jobids.append(
                submit_slurm(
                    sample_cmds[-1],
                    threads=args.threads,
                    slurm_log=args.slurm_logs_dir / f"{sid.sid}.log",
                    dependencies=[demux_jobid],
                )
            )

        meta_df = []
        for s, row in request.samples.build_index("sid").items():
            meta_df.append(
                make_sample_metadata(
                    s,
                    only_one(row),
                    results_output_dir / s / "counts/raw_feature_bc_matrix.h5",
                )
            )
        write_data(args.output_dir / "metadata.tsv", pd.concat(meta_df))

        analysis_cmd = cmd + [
            "analyze",
            "--input-dir",
            results_output_dir,
            "--output-dir",
            analysis_output_dir,
            args.request,
            args.output_dir / "metadata.tsv",
        ]

        submit_slurm(
            analysis_cmd,
            threads=args.threads,
            slurm_log=args.slurm_logs_dir / "analysis.log",
            dependencies=sample_jobids,
        )


def run_task(
    args: argparse.Namespace,
    task_name: str | None = None,
) -> None:
    """Run task and its related task helper functions.

    * Update config, using parsed arguments, request, and metadata
    * Check input requirements
    * Make the local temp directories for the staged data and output directories
    * Stage the data, copy `<task>.input` from global config to `<task>.inputs` from local config
        - All Pathy objects are in the local config as files in local temporary directory paths
        - Stage inputs from config, using global and local configs (make s3 and gs inputs available locally)
    * Run task with 1 or more step functions, where each step runs with step analysis with step helper functions
    * Check all paths in `<task>.outputs` from local config exist
    * Copy files from `<task>.outputs` from local config to `<task>.outputs` from global config
    * Remove temporary files.

    """
    task_name = args.command if task_name is None else task_name
    logging.info(f"[Task] {task_name}")

    args.local_temp_dir = parse_path(tempfile.mkdtemp(dir=args.local_temp_dir))
    args.shared_temp_dir = parse_path(tempfile.mkdtemp(dir=args.shared_temp_dir or None))
    delete_shared_temp_dir = False

    request = AnalysisRequest(args.request)
    request.load_sample_metadata(args.metadata)

    global_config, local_config, sample_configs = load_config(
        request,
        args,
        settings=get_scrna_settings(args, request),
    )
    global_config_data = global_config.data
    local_config_data = local_config.data

    inputs = global_config_data[task_name]["input"]
    local_inputs = local_config_data[task_name]["input"]

    outputs = global_config_data[task_name]["output"]
    local_outputs = local_config_data[task_name]["output"]

    # scratch = global_config_data[task_name]["scratch"]
    local_scratch = local_config_data[task_name]["scratch"]

    global_input_files: list[Path] = []
    global_input_dirs: list[Path] = []
    global_output_files: list[Path] = []
    global_output_dirs: list[Path] = []

    local_input_files: list[Path] = []
    local_input_dirs: list[Path] = []
    local_output_files: list[Path] = []
    local_output_dirs: list[Path] = []

    bin_files_directories_from_config(inputs or {}, global_input_files, global_input_dirs)
    bin_files_directories_from_config(outputs or {}, global_output_files, global_output_dirs)
    bin_files_directories_from_config(local_inputs or {}, local_input_files, local_input_dirs)
    bin_files_directories_from_config(local_outputs or {}, local_output_files, local_output_dirs)

    for _, (global_s, local_s) in sample_configs.items():
        bin_files_directories_from_config(global_s.data["input"], global_input_files, global_input_dirs)
        bin_files_directories_from_config(local_s.data["input"], local_input_files, local_input_dirs)

    check_exist(dirs=global_input_dirs, files=global_input_files)
    make_dirs(dirs=local_input_dirs, files=local_input_files)
    make_dirs(dirs=local_output_dirs, files=local_output_files)
    make_dirs(dirs=global_output_dirs, files=global_output_files)

    stage_data(inputs, local_inputs)
    for _, (global_s, local_s) in sample_configs.items():
        stage_data(global_s.data["input"], local_s.data["input"])

    global_config.save_data(local_outputs["config_global"])
    local_config.save_data(local_outputs["config_local"])

    for step in inputs["steps"]:
        # not all steps require input/output/scratch/params sections in config
        logging.info(f"[Step] {step}")
        step_config = local_config_data[task_name][step]
        step_outputs = {}
        for key in step_config.get("output").keys() & local_outputs.keys():
            src, dst = step_config.get("output")[key], local_outputs[key]
            if src and isinstance(src, (str, Path)) and dst and isinstance(dst, (str, Path)):
                step_outputs[key] = dst

        run_step(
            inputs=step_config.get("input"),
            outputs=step_outputs,
            scratch=step_config.get("output"),
            params=step_config.get("params"),
            request=request,
            step_name=step,
            sample_configs=sample_configs,
        )

    # set_path_readonly(src)
    stage_data(local_outputs, outputs)
    check_exist(dirs=global_output_dirs, files=global_output_files)

    if not args.keep_temp:
        for _, temp_dir in local_scratch.items():
            try:
                logging.info(f"[Task] deleting temp scratch dir from {temp_dir}")
                shutil.rmtree(temp_dir)
            except OSError:
                pass

    if not args.keep_temp and args.local_temp_dir:
        try:
            logging.info(f"[Task] deleting temporary files from {args.local_temp_dir}")
            shutil.rmtree(args.local_temp_dir)
        except OSError:
            pass

    if delete_shared_temp_dir and args.shared_temp_dir:
        try:
            logging.info(f"[Task] deleting temporary shared files from {args.shared_temp_dir}")
            shutil.rmtree(args.shared_temp_dir)
        except OSError:
            pass


def run_step(
    inputs: ImmutableStrMapping,
    outputs: ImmutableStrMapping,
    scratch: ImmutableStrMapping,
    params: ImmutableStrMapping,
    request: AnalysisRequest,
    step_name: str,
    sample_configs: dict[str, tuple[PackageConfig, PackageConfig]],
) -> None:
    """Run step and its related step helper functions, using StepHelper."""
    logging.info(f"[Step] Executing {step_name}")

    with StepHelper(
        input_config=[inputs],
        output_config=[outputs],
        scratch_config=[scratch],
        copy_across_config=[(scratch, outputs)],
    ):
        # TODO: add a boilerplate function where function name can be interpreted as function
        # (something that is similar to eval(function_name), but safer to use)
        # for now, use name-function mapping here
        STEP_FUNC_MAP[step_name](
            inputs=inputs,
            outputs=outputs,
            scratch=scratch,
            params=params,
            request=request,
            sample_configs=sample_configs,
        )


def arg_parser() -> argparse.ArgumentParser:
    """Create an argument parser for AskCell pipeline."""
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(title="Commands", dest="command")

    add_logging_args(parser.add_argument_group("Logging"))

    parser.add_argument("--keep-temp", action="store_true", help="if set, keep intermediate files")
    parser.add_argument(
        "--slurm-logs-dir",
        metavar="DIR",
        type=parse_path,
        help="defines slurm logs directory, default is <output_dir_OR_url> / logs",
    )
    parser.add_argument(
        "--shared-temp-dir",
        metavar="DIR",
        type=parse_path,
        help="path to a shared storage location for temporary files",
    )
    parser.add_argument(
        "--local-temp-dir",
        metavar="DIR",
        type=parse_path,
        help="path to a local storage location for temporary files",
    )
    cmd = commands.add_parser("run", help="run sc-RNA seq analysis workflow")
    cmd.add_argument("request", metavar="REQUEST", type=parse_path, help="analysis request")
    cmd.add_argument("metadata", metavar="METADATA", type=parse_path, help="sample metadata")
    cmd.add_argument(
        "--input-dir",
        metavar="DIR_OR_URL",
        type=parse_path,
        required=True,
        help="path to sequencing run data",
    )
    cmd.add_argument(
        "--output-dir",
        metavar="DIR_OR_URL",
        type=parse_path,
        required=True,
        help="output directory or url",
    )
    cmd.add_argument(
        "--scratch-dir",
        metavar="DIR",
        type=parse_path,
        help="non-shared storage location for temporary files",
    )
    cmd.add_argument("--threads", type=int, default=8, help="number of threads to use")
    cmd.set_defaults(func=run_workflow)

    cmd = commands.add_parser("demux", help="a demultiplex scRNA-seq 10x data")
    cmd.add_argument("request", metavar="REQUEST", type=parse_path, help="analysis request")
    cmd.add_argument("metadata", metavar="METADATA", type=parse_path, help="sample metadata")
    cmd.add_argument(
        "--run-dir",
        metavar="DIR_OR_URL",
        type=parse_path,
        required=True,
        help="path to sequencing run data",
    )
    cmd.add_argument(
        "--output-dir",
        metavar="DIR_OR_URL",
        type=parse_path,
        required=True,
        help="output directory or url",
    )
    cmd.add_argument(
        "--scratch-dir",
        metavar="DIR",
        type=parse_path,
        help="non-shared storage location for temporary files",
    )
    cmd.add_argument("--threads", type=int, default=8, help="number of threads to use")
    cmd.set_defaults(func=run_task)

    cmd = commands.add_parser("count", help="align and create count data for a single sample")
    cmd.add_argument("request", metavar="FILE", type=parse_path, help="analysis request")
    cmd.add_argument("metadata", metavar="FILE", type=parse_path, help="sample metadata")
    cmd.add_argument(
        "--sid",
        metavar="STR",
        type=str,
        required=True,
        help="sample identifier, must be one in metadata",
    )
    cmd.add_argument(
        "--reads-dir",
        metavar="DIR_OR_URL",
        type=parse_path,
        required=True,
        help="path containing reads in fq format",
    )
    cmd.add_argument(
        "--output-dir",
        metavar="DIR_OR_URL",
        type=parse_path,
        required=True,
        help="output directory or url",
    )
    cmd.add_argument(
        "--scratch-dir",
        metavar="DIR",
        type=parse_path,
        help="non-shared storage location for temporary files",
    )
    cmd.add_argument("--threads", type=int, default=8, help="number of threads to use")
    cmd.set_defaults(func=run_task)

    cmd = commands.add_parser("analyze", help="perform post-count secondary analysis for scRNA-seq")
    cmd.add_argument("request", metavar="FILE", type=parse_path, help="analysis request")
    cmd.add_argument("metadata", metavar="FILE", type=parse_path, help="sample metadata")
    cmd.add_argument(
        "--input-dir",
        metavar="DIR_OR_URL",
        type=parse_path,
        required=True,
        help="path to directory with samples raw counts matrices results (in mtx or tsv or h5 format)",
    )
    cmd.add_argument(
        "--output-dir",
        metavar="DIR_OR_URL",
        type=parse_path,
        required=True,
        help="output directory or url",
    )
    cmd.add_argument(
        "--scratch-dir",
        metavar="DIR",
        type=parse_path,
        help="scratch directory",
    )
    cmd.add_argument("--threads", type=int, default=8, help="number of threads to use")
    cmd.set_defaults(func=run_task)

    add_config_args(parser)

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

    if not args.command or not args.func:
        parser.print_help()
        return

    init_logging(
        level=args.log_level,
        log_filename=(args.slurm_logs_dir or args.output_dir) / "logs/.log",
        testing=args.testing,
    )
    run_workflow(args)
