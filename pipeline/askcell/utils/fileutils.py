"""fileutils.py to provide file related utility functions."""

from __future__ import annotations

import json

from typing import IO, Any

import toml

from askcell.utils.pathutils import Path, PathLike, parse_path
from askcell.utils.types import MutableStrMapping, StrMapping


def load_toml(
    filename_or_file: PathLike | IO[str],
    *,
    dictclass: type[MutableStrMapping] = dict,
    decoder: toml.TomlDecoder[Any] | None = None,
) -> MutableStrMapping:
    """Load data from a TOML file.

    Args:
        filename_or_file: input filename or file object
        dictclass: class to use for mappings
        decoder: optional TomlDecoder

    Returns:
        file contents as a StrMapping

    """
    data: MutableStrMapping

    if isinstance(filename_or_file, (str, Path)):
        with open(parse_path(filename_or_file)) as infile:
            data = toml.load(infile, _dict=dictclass, decoder=decoder)
    else:
        data = toml.load(filename_or_file, _dict=dictclass, decoder=decoder)

    return data


def save_toml(
    filename_or_file: PathLike | IO[str],
    data: StrMapping,
    *,
    mode: str = "w",
    encoder: toml.TomlEncoder | None = None,  # type: ignore
) -> None:
    """Write data to a TOML file.

    Args:
        filename_or_file: input filename or file object
        data: dictionary of arbitrary data to be written; values must be mappable to TOML data types
        mode: file mode (default=w)
        encoder: optional TOML data type encoder

    """
    if isinstance(filename_or_file, (str, Path)):
        with parse_path(filename_or_file).open(mode) as outfile:
            toml.dump(data, outfile, encoder=encoder)
    else:
        toml.dump(data, filename_or_file, encoder=encoder)


def load_json(
    filename_or_file: PathLike | IO[str],
    *,
    dictclass: type[MutableStrMapping] = dict,
) -> Any:
    """Load data from a JSON file.

    Args:
        filename_or_file: input filename or file object
        dictclass: class to use for mappings

    Returns:
        JSON file contents

    """
    if isinstance(filename_or_file, (str, Path)):
        with Path(filename_or_file).open() as infile:
            data: StrMapping = json.load(infile, object_hook=dictclass)
    else:
        data = json.load(filename_or_file, object_hook=dictclass)

    return data
