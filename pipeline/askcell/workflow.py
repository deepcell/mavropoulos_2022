"""scRNA-seq workflow-related functions."""

from __future__ import annotations

import argparse
import logging

from dataclasses import dataclass
from subprocess import PIPE, CalledProcessError, Popen
from types import TracebackType
from typing import Any, Iterable, Iterator, Sequence

import pandas as pd

from . import (  # set_path_readonly,
    Path,
    PathLike,
    Pathy,
    PosixPath,
    SomePathLikes,
    WindowsPath,
    only_one,
    parse_path,
    parse_paths,
)
from .scrna.request import AnalysisRequest
from .utils.config import (
    PackageConfig,
    PackageResource,
    load_package_data,
    load_package_schema,
)
from .utils.shellutils import Command, run_cmd
from .utils.stringutils import strip_prefix
from .utils.types import ImmutableStrMapping


@dataclass
class StepHelper:
    """This is a class to perform various checks during each step initialization and exit.

    It performs the following functions during step initialization:

       1. Check the input requirements
       2. Make the directories in `<task>.<step>.output` if it doesn't exist

    It performs the following functions during step exit:

       1. Check the paths in `<task>.<step>.output` exist, else throws exception
       2. Copy files from `<task>.<step>.output` to `<task>.output`

    Args:
        input_config:       TypedConfigSection for input
        output_config:      TypedConfigSection for output
        scratch_config:     TypedConfigSection for scratch
        input_files:        input files
        input_dirs:         input directories
        output_files:       output files
        output_dirs:        output directories
        scratch_files:      scratch files
        scratch_dirs:       scratch directories
        copy_across_config: List of paired TypedConfigSections

    """

    input_files: list[Path]
    input_dirs: list[Path]
    output_file: list[Path]
    output_dirs: list[Path]
    scratch_files: list[Path]
    scratch_dirs: list[Path]
    copy_across_config: list[tuple[ImmutableStrMapping, ImmutableStrMapping]]

    def __init__(
        self,
        *,
        input_config: Iterable[ImmutableStrMapping] | None = None,
        output_config: Iterable[ImmutableStrMapping] | None = None,
        scratch_config: Iterable[ImmutableStrMapping] | None = None,
        input_files: SomePathLikes | None = None,
        input_dirs: SomePathLikes | None = None,
        output_files: SomePathLikes | None = None,
        output_dirs: SomePathLikes | None = None,
        scratch_files: SomePathLikes | None = None,
        scratch_dirs: SomePathLikes | None = None,
        copy_across_config: Iterable[tuple[ImmutableStrMapping, ImmutableStrMapping]] | None = None,
    ) -> None:
        """Initialize files and directories associated with a Step."""
        self.input_files = list(parse_paths(input_files or []))
        self.input_dirs = list(parse_paths(input_dirs or []))
        self.output_files = list(parse_paths(output_files or []))
        self.output_dirs = list(parse_paths(output_dirs or []))
        self.scratch_files = list(parse_paths(scratch_files or []))
        self.scratch_dirs = list(parse_paths(scratch_dirs or []))

        for config in input_config or []:
            bin_files_directories_from_config(config or {}, self.input_files, self.input_dirs)

        for config in output_config or []:
            bin_files_directories_from_config(config or {}, self.output_files, self.output_dirs)

        for config in scratch_config or []:
            bin_files_directories_from_config(config or {}, self.scratch_files, self.scratch_dirs)

        self.copy_across_config = [cc for cc in copy_across_config or [] if cc]

    def __enter__(self) -> None:
        """Ensure input requirements are satisfied and create directories for outputs."""
        logging.info("[Step] check all step inputs exist")
        check_exist(dirs=self.input_dirs, files=self.input_files)

        make_dirs(dirs=self.scratch_dirs, files=self.scratch_files)
        make_dirs(dirs=self.output_dirs, files=self.output_files)

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        """Call teardown functions prior to Step exit."""
        for src_config, dst_config in self.copy_across_config or []:
            for src, dst in pair_paths_by_key(src_config, dst_config):
                if src and dst and isinstance(src, (str, Path)) and isinstance(dst, (str, Path)):
                    logging.info(f"[Step] copying output files from {src} to {dst}")
                    move_files(src, dst)

        logging.info("[Step] check all step outputs exist")
        check_exist(files=self.output_files)


def get_scrna_settings(
    args: argparse.Namespace,
    request: AnalysisRequest,
) -> list[str]:
    """Get data fragments to be used to update the config for sc-RNAseq analysis.

    settings will override the pre-existing values. It will include the
    parameters passed in as arguments as well as the information provided in the request.

    Args:
        args: arguments
        request: AnalysisRequest object with run-specific information

    Returns:
        setting options

    """
    settings = args.set_config or []
    settings += [
        f"input.request={args.request}",
        f"input.metadata={args.metadata}",
        f"workflow.runid={request.runid}",
        f'workflow.fcid={request.runid.rsplit("_", 1)[1]}',
        f"workflow.requestid={request.issueid}",
        f"output.output_dir={args.output_dir}",
        f"scratch.scratch_dir={args.local_temp_dir}",
        f"params.max_mito_pct_allowed={request.max_mito_pct_allowed}",
        f"params.max_ribo_pct_allowed={request.max_ribo_pct_allowed}",
        f"params.regress_out_cellcycle={request.regress_out_cellcycle}",
        f"params.n_pcs={request.n_pcs}",
        f"params.n_neighbors={request.n_neighbors}",
    ]

    if args.command == "analyze":
        settings += [f"input.input_dir={args.input_dir}"]

    if args.command == "count":
        settings += [f"input.reads_dir={args.reads_dir}"]
        settings += [f"workflow.sampleid={args.sid}"]

    if args.command == "demux":
        settings += [f"input.run_dir={args.run_dir}"]

    if args.keep_temp is not None:
        settings += [f"params.keep_temp={args.keep_temp}"]

    return settings


def load_config(
    request: AnalysisRequest,
    args: argparse.Namespace,
    data_files: SomePathLikes | None = None,
    settings: Iterable[str] | None = None,
) -> tuple[PackageConfig, PackageConfig, dict[str, tuple[PackageConfig, PackageConfig]]]:
    """Load a typed configuration file.

    Update the config to support temporary local files, originated from the cloud storage.
    FIXME: Update with config features from BIOIN-297_config.

    Args:
        request: analysis request
        args: arguments
        data_files: data paths
        settings: sequence of configuration settings data of the form section.option=value

    Returns:
        Typed configuration of the based on the `config/config.toml` schema and sample-level configs

    """
    data = PackageResource("askcell.config.data.scrna", "*.toml")
    # TODO: add data_files support

    settings = list(settings) if settings is not None else []
    schema = load_package_schema(PackageResource("askcell.config.schema.scrna", "*.toml"))

    global_data = load_package_data(schema, data=data, settings=settings)
    global_config = PackageConfig(schema, global_data)

    local_settings = settings.copy()
    local_settings += list(get_local_paths_settings(global_config.data, args.local_temp_dir))

    local_data = load_package_data(schema, data=data, settings=local_settings)
    local_config = PackageConfig(schema, local_data)

    sample_configs: dict[str, tuple[PackageConfig, PackageConfig]] = {}
    s_schema = schema["sample"]
    for sid, sample in request.samples.build_index(*["sid"]).items():
        s = only_one(sample)

        settings = [f"sampleid={sid}", f"label={s.label}", f"runid={s.runid}"]
        if s.raw_count is not None:
            settings += [f"input.h5_count={s.raw_count}"]

        s_data = load_package_data(
            s_schema,
            data=None,
            settings=settings,
        )

        sample_global_config = PackageConfig(s_schema, s_data)
        sample_local_settings = list(get_local_paths_settings(sample_global_config.data, args.local_temp_dir))

        sample_local_data = load_package_data(s_schema, data=s_data, settings=sample_local_settings)
        sample_local_config = PackageConfig(s_schema, sample_local_data)

        sample_configs[sid] = (sample_global_config, sample_local_config)

    return global_config, local_config, sample_configs


def get_local_paths_settings(
    data: ImmutableStrMapping,
    local_temp_dir: PathLike,
) -> Iterator[str]:
    """Get local paths settings for the paths objects in s3 or gcs.

    This function converts to s3 and gs path to posix path in temp directory,
    where bucket name is replaced with the local path.

    Args:
        data: data package resources, paths, or fragments
        local_temp_dir: local temporary directory

    Yields:
        input or output local path setting, in form of x.y=z

    """
    for section in {"input", "output"}.intersection(data.keys()):
        for k, v in data[section].items():
            value = parse_path(v) if isinstance(v, (str, Path)) else v
            if isinstance(value, Pathy):
                yield (f"{section}.{k}={parse_path(local_temp_dir) / value.prefix}")


def write_data(
    filename: PathLike,
    data: pd.DataFrame,
    use_index: bool = False,
) -> None:
    """Write data to a file.

    Args:
        filename: filename
        data: data in tabular format
        use_index: if set True, then include index in the output as column

    """
    filename = parse_path(filename)
    if filename.name.endswith(".tsv"):
        data.to_csv(filename, sep="\t", index=use_index)
    elif filename.name.endswith(".csv"):
        data.to_csv(filename, sep=",", index=use_index)
    elif filename.name.endswith(".pq"):
        data.to_parquet(filename, index=use_index)
    else:
        raise ValueError(f"Unknown file type: {filename}")


def stage_data(
    source: ImmutableStrMapping,
    destination: ImmutableStrMapping,
) -> None:
    """Stage files and directories from source to destination.

    Args:
        source: source files and directories
        destination: destination files and directories

    """
    for src, dst in pair_paths_by_key(source, destination):
        move_files(src, dst)


def move_files(source: PathLike, destination: PathLike) -> None:
    """Move files from source to destination within local paths or sync files.

    If the source is a directory, then the contents within a directory is copied to new destination.
    Creating the source directory-named directory within a new destination directory is not supported.

    Args:
        source: source data in PosixPath or WindowsPath or S3Path # FIXME: + GCS path
        destination: destination data in PosixPath or WindowsPath or S3Path # FIXME: + GCS path

    """
    source = parse_path(source)
    destination = parse_path(destination)

    if source == destination:
        return

    if destination.is_dir():
        make_directory(destination)

    # FIXME: any use case for gs <-> s3 moving files?
    if (isinstance(source, Pathy) and source.scheme == "gs") or (
        isinstance(destination, Pathy) and destination.scheme == "gs"
    ):
        if source.is_dir():
            run_cmd(Command(["gsutil", "-m", "cp", "-R", str(source), str(destination.parent)]))
        else:
            run_cmd(Command(["gsutil", "cp", str(source), str(destination)]))

    elif (isinstance(source, Pathy) and source.scheme == "s3") or (
        isinstance(destination, Pathy) and destination.scheme == "s3"
    ):
        run_cmd(Command(["aws", "s3", "sync" if source.is_dir() else "cp", str(source), str(destination)]))

    else:
        # FIXME: switch to 'mv' for local or 'rsync' for non-local (across the network) option
        cmd = Command(["cp", source, destination])
        if source.is_dir():
            cmd.append("-rT")
        run_cmd(cmd)


def bin_files_directories_from_config(
    config: ImmutableStrMapping,
    files: list[Path],
    dirs: list[Path],
) -> None:
    """Put files and directories into different bins from a single section.

    Args:
        config: Typed config
        files: zero or more path-like file instances
        dirs: zero or more path-like directory instances

    """
    for key, value in config.items():
        if not value or not isinstance(value, (Path, PosixPath, WindowsPath, Pathy)):
            continue
        if key.endswith("_dir"):
            dirs.append(value)
        else:
            files.append(value)


def pair_paths_by_key(
    source: ImmutableStrMapping,
    destination: ImmutableStrMapping,
) -> Iterator[tuple[PathLike, PathLike]]:
    """Pair paths by shared keys from source and destination.

    Args:
        source: source files
        destination: destination files

    Yields:
        paired source and destination paths

    Examples:
        >>> scratch = {'alignment': '/scratch/deduped.bam', 'score': '/scratch/score.txt'}
        >>> output = {'alignment': '/output/aligned.bam'}
        >>> list(pair_paths_by_key(scratch, output))
        [('/scratch/deduped.bam', '/output/aligned.bam')]

    """
    for key in source.keys() & destination.keys():
        src, dst = source[key], destination[key]
        if src and isinstance(src, (str, Path)) and dst and isinstance(dst, (str, Path)):
            yield src, dst


def check_exist(
    dirs: SomePathLikes | None = None,
    files: SomePathLikes | None = None,
) -> None:
    """Check whether file/directories exist.

    Args:
        dirs: list of directories
        files: list of files

    Raises:
        ValueError: Input/Output <filename> requirement is not met.

    """
    for p in parse_paths(dirs or []):
        if not p.exists():
            raise ValueError(f"Directory {p} does not exist.")
        elif not p.is_dir():
            raise ValueError(f"{p} is not a directory, where {p} is expected to be a directory.")

    for p in parse_paths(files or []):
        if not p.exists():
            raise ValueError(f"File {p} does not exist.")
        elif p.is_dir():
            raise ValueError(f"{p} is a directory, where {p} is expected to be a file.")


def make_dirs(
    dirs: SomePathLikes | None = None,
    files: SomePathLikes | None = None,
) -> None:
    """Make directories or parent directories for files.

    Args:
        dirs: list of directories
        files: list of files

    Raises:
        ValueError: Input/Output <filename> requirement is not met.

    """
    for p in parse_paths(dirs or []):
        make_directory(p, parents=True, exist_ok=True)

    for p in parse_paths(files or []):
        if not p.parent.exists():
            make_directory(p.parent, parents=True, exist_ok=True)


def make_directory(
    path: PathLike,
    parents: bool = True,
    exist_ok: bool = True,
) -> None:
    """Make directory or create key in cloud storage.

    Args:
        path: directory name
        parents: if set True, missing parent directory is allowed and will make one
        exist_ok: if set True, pre-existing directory with same name does not throw error

    """
    path = parse_path(path)

    if isinstance(path, Pathy):
        (path / ".log").touch(exist_ok=exist_ok)

    else:
        path.mkdir(parents=parents, exist_ok=exist_ok)


def submit_slurm(
    command: Command,
    *,
    threads: int,
    slurm_log: PathLike,
    dependencies: Sequence[str | None] | None = None,
    wait: bool = False,
) -> str:
    """Command to submit slurm job."""
    # FIXME
    sbatch = Command(
        [
            "sbatch",
            "-p",
            "c2-standard-16",
            "-n",
            1,
            "-c",
            threads,
            "-o",
            slurm_log,
            "--parsable",
        ]
    )

    if wait:
        sbatch += ["--wait"]

    if dependencies is not None:
        jobids = ":".join(filter(None, dependencies or []))
        sbatch += ["-d", f"afterok:{jobids}"]

    command = command.convert_to_shell()

    proc = Popen(sbatch.command_list, stdin=PIPE, stdout=PIPE, encoding="utf-8")
    out, _ = proc.communicate(f"#!/bin/sh\n\n{command.shell_string}\n")

    if proc.returncode != 0:
        raise CalledProcessError(proc.returncode, f"sbatch submission failed with exit code {proc.returncode}")

    return parse_jobid(out)


def parse_jobid(jobid: str) -> str:
    """Parse job id from batch job message."""
    return strip_prefix(jobid.rstrip(), "Submitted batch job")


def convert_pct_to_fraction(value: Any) -> float:
    """Convert the percentage to fraction.

    Args:
        value: value with % or without % at the end

    Returns:
        converted value or original value, if no conversion needed

    """
    if str(value).endswith("%"):
        value = value.rsplit("%")[0].strip()

    if value.replace(".", "", 1).isdigit():
        return float(value) * 0.01

    raise ValueError(f"{value} is not a numerical value")


def remove_thousand_separator(value: Any) -> int:
    """Remove the thousandth comma separator for numerical value.

    Args:
        value: value with thousandth comma separator to int

    Returns:
        converted value or original value, if no conversion needed

    """
    value = str(value).replace(",", "")

    if value.isnumeric():
        return int(value)

    raise ValueError(f"{value} is not a numerical value")


def check_qc() -> None:
    """Check the equality.

    TODO: Check other repo for similar functions.

    Check whether <left_value> <operator_type> <right_value> holds or not.
    where the operator types are GT = ">", GE = ">=", LT = "<", LE = "<=", NE = "!=", EQ = "=="

    Args:
        operator_type: GT, GE, LT, LE, NE, and EQ
        left_value: numerical value or a string or a boolean
        right_value: numerical value or a string or a boolean

    Returns:
        boolean, True if the inequality or equation holds

    """
    pass
