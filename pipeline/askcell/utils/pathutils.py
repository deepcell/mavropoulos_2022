"""Path related utility functions. FIXME: Remove once dcutil becomes dependency."""

from __future__ import annotations

import stat

from fnmatch import fnmatch
from itertools import chain
from pathlib import Path, PosixPath, WindowsPath
from typing import Iterable, Iterator, Union

from pathy import Pathy


__all__ = [
    "validate_path",  # validate path instance based on flags
    "validate_paths",  # validate path instances based on flags
    "parse_path",  # convert a string or other path-like instance into a concrete Path instance
    "parse_paths",  # convert zero or more string or other path-like instances and yield concrete Path instance
    "Path",  # re-export for convenience
    "PosixPath",  # re-export for convenience
    "WindowsPath",  # re-export for convenience
    "Pathy",  # re-export for convenience
    "PathLike",  # Type for a single path-like instance (str or Path)
    "SomePathLikes",  # Type for zero or more path-like instances
    "set_path_readonly",  # set read-only permissions
]


# Universal type definitions
PathLike = Union[str, Path]
SomePathLikes = Union[PathLike, Iterable[PathLike]]


def validate_path(
    path: Path,
    *,
    check_file: bool = False,
    check_dir: bool = False,
    check_exists: bool = False,
    allow_posix: bool = True,
    allow_s3: bool = True,
    allow_gcs: bool = True,
    allow_windows: bool = False,
) -> None:
    """Validate path is an allowed type or raise an exception.

    Args:
        check_file: If True then reject paths that do not refer to an existing file (default is False)
        check_dir: If True then reject paths that do not refer to an existing directory (default is False)
        check_exists: If True then reject paths that do not refer to an existing file or directory (default is False)
        allow_posix: If False then reject POSIX paths (default is True)
        allow_s3: If False then reject S3 paths (default is True)
        allow_gcs: If False then reject GCS paths (default is True)
        allow_windows: If False then reject Windows paths (default is False)

    """
    if isinstance(path, PosixPath):
        if not allow_posix:
            raise ValueError(f"POSIX path is not allowed in this context: {path}")

    elif isinstance(path, Pathy):
        if path.scheme == "s3" and not allow_s3:
            raise ValueError(f"S3 path is not allowed in this context: {path}")
        elif path.scheme == "gs" and not allow_gcs:
            raise ValueError(f"GCS path is not allowed in this context: {path}")

    elif isinstance(path, WindowsPath):
        if not allow_windows:
            raise ValueError(f"Windows path is not allowed in this context: {path}")

    else:
        raise ValueError(f"Unknown path type: {path}")

    if check_file and not path.is_file():
        raise ValueError(f"Path does not refer to an existing file: {path}")

    if check_dir and not path.is_dir():
        raise ValueError(f"Path does not refer to an existing directory: {path}")

    if check_exists and not path.exists():
        raise ValueError(f"Path does not exist: {path}")


def validate_paths(
    paths: Path | Iterable[Path],
    *,
    check_file: bool = False,
    check_dir: bool = False,
    check_exists: bool = False,
    allow_posix: bool = True,
    allow_s3: bool = True,
    allow_gcs: bool = True,
    allow_windows: bool = False,
) -> Iterator[Path]:
    """Validate each paths is an allowed type and yield it; otherwise raise an exception.

    Args:
        paths: a single path or an iterable of paths
        check_file: If True then reject paths that do not refer to an existing file (default is False)
        check_dir: If True then reject paths that do not refer to an existing directory (default is False)
        check_exists: If True then reject paths that do not refer to an existing file or directory (default is False)
        allow_posix: If False then reject POSIX paths (default is True)
        allow_s3: If False then reject S3 paths (default is True)
        allow_windows: If False then reject Windows paths (default is False)

    Yields:
        :class:`pathlib.Path` for every valid path

    """
    if isinstance(paths, Path):
        paths = [paths]

    for path in paths:
        validate_path(
            path,
            check_file=check_file,
            check_dir=check_dir,
            check_exists=check_exists,
            allow_posix=allow_posix,
            allow_s3=allow_s3,
            allow_gcs=allow_gcs,
            allow_windows=allow_windows,
        )
        yield path


def parse_path(
    path: PathLike,
    *pathsegments: str,
    expanduser: bool = True,
    resolve: bool = False,
    check_file: bool = False,
    check_dir: bool = False,
    check_exists: bool = False,
    allow_posix: bool = True,
    allow_s3: bool = True,
    allow_gcs: bool = True,
    allow_windows: bool = False,
) -> Path:
    r"""Create a new PosixPath, WindowsPath, or S3Path.

    Args:
        base: base path
        \*pathsegments: additional path segments
        expanduser: if True, expanded `~` and `~user` path prefixes (default is True)
        resolve: if True, make paths absolute, resolving any symlinks (default is False)
        check_file: If True then reject paths that do not refer to an existing file (default is False)
        check_dir: If True then reject paths that do not refer to an existing directory (default is False)
        check_exists: If True then reject paths that do not refer to an existing file or directory (default is False)
        allow_posix: If False then reject POSIX paths (default is True)
        allow_s3: If False then reject S3 paths (default is True)
        allow_windows: If False then reject Windows paths (default is False)

    Returns:
        New path instance (PosixPath, WindowsPath or S3Path depending on input)

    Examples:
        >>> parse_path('/absolute/posix/path/to/data')
        PosixPath('/absolute/posix/path/to/data')
        >>> parse_path('relative/posix/path/to/data')
        PosixPath('relative/posix/path/to/data')
        >>> parse_path('s3://bucket/path/to/data')
        Pathy('s3://bucket/path/to/data')

    """
    p: Path
    if isinstance(path, Path):
        p = path
    elif isinstance(path, str) and (path.startswith("s3:") or path.startswith("gs:")):
        p = Pathy(path)
    else:
        p = Path(path)

    if pathsegments:
        p = p.joinpath(*pathsegments)

    if isinstance(p, PosixPath):
        if expanduser:
            p = p.expanduser()

        if resolve:
            p = p.resolve()

    validate_path(
        p,
        check_file=check_file,
        check_dir=check_dir,
        check_exists=check_exists,
        allow_posix=allow_posix,
        allow_s3=allow_s3,
        allow_gcs=allow_gcs,
        allow_windows=allow_windows,
    )

    return p


def parse_paths(
    paths: SomePathLikes,
    *,
    expanduser: bool = True,
    resolve: bool = False,
    check_file: bool = False,
    check_dir: bool = False,
    check_exists: bool = False,
    allow_posix: bool = True,
    allow_s3: bool = True,
    allow_gcs: bool = True,
    allow_windows: bool = False,
) -> Iterator[Path]:
    """Yield a single PathLike or a sequence of PathLike objects as Path objects.

    Args:
        paths: a single path or an iterable of paths
        expanduser: if True, expanded `~` and `~user` path prefixes (default is True)
        resolve: if True, make paths absolute, resolving any symlinks (default is False)
        check_file: If True then reject paths that do not refer to an existing file (default is False)
        check_dir: If True then reject paths that do not refer to an existing directory (default is False)
        check_exists: If True then reject paths that do not refer to an existing file or directory (default is False)
        allow_posix: If False then reject POSIX paths (default is True)
        allow_s3: If False then reject S3 paths (default is True)
        allow_windows: If False then reject Windows paths (default is False)

    Yields:
        :class:`pathlib.Path` for every path

    """
    if isinstance(paths, (str, Path)):
        paths = [paths]

    for path in paths:
        yield parse_path(
            path,
            expanduser=expanduser,
            resolve=resolve,
            check_file=check_file,
            check_dir=check_dir,
            check_exists=check_exists,
            allow_posix=allow_posix,
            allow_s3=allow_s3,
            allow_gcs=allow_gcs,
            allow_windows=allow_windows,
        )


def set_path_readonly(path: PathLike, recursive: bool = True) -> None:
    """Set path readonly.

    Args:
        filepath:  file or directory, that should be read only
        recursive: If true and filepath is a directory path, \
            then dir itself plus all subdirectories/files change to read only

    """
    path = parse_path(path)

    if not isinstance(path, PosixPath):
        raise ValueError("requires POSIX path")

    paths: list[Iterable[Path]] = [[path]]

    if recursive and path.is_dir():
        paths.append(path.glob("**/*"))

    mask = ~(stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH)
    for filename in chain.from_iterable(paths):
        filename.chmod(filename.stat().st_mode & mask)


def match_path(path: PathLike, pattern: str) -> bool:
    """Match path's final filename against pattern.

    Args:
        path: path whose terminal filename is to be matched
        pattern: Glob-style filename pattern

    Returns:
        True if a match is found, False otherwise

    Examples:
        >>> match_path('~/foo/bar/baz.txt', '*.txt')
        True
        >>> match_path('~/foo/bar.csv/baz.txt', '*.csv')
        False
        >>> match_path('s3://bucket/foo/bar/baz.txt', '*.txt')
        True
        >>> match_path('s3://bucket/foo/bar.csv/baz.txt', '*.csv')
        False

    """
    return fnmatch(parse_path(path).name, pattern)
