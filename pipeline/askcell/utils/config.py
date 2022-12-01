"""Common workflow configuration handling."""

from __future__ import annotations

import argparse
import importlib.resources as resources
import logging

from dataclasses import dataclass
from fnmatch import fnmatch
from types import MappingProxyType
from typing import Iterable, Iterator, Union

from .fileutils import load_toml, save_toml
from .freeze import freeze
from .interpolate import interpolate
from .merge import merge
from .pathutils import Path, PathLike, SomePathLikes, parse_path, parse_paths
from .types import (
    ImmutableStrMapping,
    MutableStrMapping,
    StrMapping,
    apply_schema,
    build_data_from_schema,
    parse_schema,
    unparse_data,
    validate_data,
)


logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


@dataclass(eq=True, frozen=True, init=False)
class PackageConfig:
    """Store a package's schema and data configuration.

    Attributes:
        schema: the package schema
        data:   the package data

    """

    schema: ImmutableStrMapping
    data: ImmutableStrMapping

    def __init__(self, schema: StrMapping, data: MutableStrMapping) -> None:
        """Validate data using schema and freeze them both.

        Args:
            schema: package schema
            data:   mutable package data to be validated and normalized

        """
        validate_data(schema, data)

        object.__setattr__(self, "schema", freeze(schema))
        object.__setattr__(self, "data", freeze(data))

    def save_data(self, path: PathLike) -> None:
        """Save package data to path.

        schema is used to unparse data. The unparsed data is then exported as a TOML file to path.

        Args:
            path: file to save package data to

        """
        save_toml(path, unparse_data(self.schema, self.data))


@dataclass(eq=True, frozen=True)
class PackageResource:
    """Define a package and a pattern for returning package resources.

    Attributes:
        package: a dotted path to a package
        pattern: Unix filename pattern

    """

    package: resources.Package
    pattern: str

    def get_resource(self) -> Iterator[Path]:
        """Yield a resource in package matching pattern.

        Yielded resources may correspond to a temporary file path and need to be processed or saved before moving on to
        other items. They are yielded in lexicographical order based on the name of the resource. This helps ensure a
        "deterministic" yield order as otherwise the yield order is subject to filesystem implementation.

        Yields:
            Path to package resource

        """
        for name in sorted(resources.contents(self.package)):
            if resources.is_resource(self.package, name) and fnmatch(parse_path(name).name, self.pattern):
                with resources.path(self.package, name) as path:
                    yield path


ConfigItem = Union[StrMapping, PackageResource, PathLike]
SomeConfigItems = Union[ConfigItem, Iterable[ConfigItem]]


def collect_config_files(paths: SomePathLikes, pattern: str = "*.toml") -> Iterator[Path]:
    r"""Collect from paths config files ending with the specified extension.

    Paths to files are yielded as-is and are not matched against the pattern. For paths that are directories, all files
    matching the pattern are yielded in lexicographical order. This search is not recursive, i.e. subdirectory files
    matching the pattern are not returned.

    The directory structure below will be used to illustrate the order in which paths are yielded::

         config/
            b.toml
            d.json

    collect_config_files(['config/d.json', 'config'], pattern='\*.toml') will yield paths from left to
    right from ['config/d.json', 'config/b.toml', 'config/c.toml']. As seen in the example, the overall
    left to right order of paths is respected, and the directory entry 'config' is replaced by matching configuration
    files in lexicographical order.

    Args:
        paths: one or a sequence of paths to configuration directories and files
        pattern: Unix filename pattern that defaults to \*.toml

    Yields:
        paths to specified configuration files

    """
    for path in parse_paths(paths):
        if path.is_dir():
            yield from sorted(path.glob(pattern))
        elif path.is_file():
            yield path
        else:
            logger.warning(f"Missing configuration file or directory: {path}")


def add_config_args(parser: argparse.ArgumentParser | argparse._ArgumentGroup) -> None:
    r"""Add configuration related parameters to an argument parser.

    Adds:
        \* -s / --set-config, which allows multiple arguments of the form key[.key]*=value to be specified
        \* -c / --config, which allows multiple configuration files or directories to be specified from the command line

    Args:
        parser: argument parser instance.

    """
    parser.add_argument(
        "-s",
        "--set-config",
        action="append",
        help="Set configuration options. Overrides contents of config files as key[.key]*=value",
    )
    parser.add_argument(
        "-c",
        "--config",
        type=Path,
        action="append",
        help="Config file or directory containing config files (option can be used multiple times)",
    )


def load_config_data(paths: SomePathLikes) -> Iterator[MutableStrMapping]:
    """Load data based on specified paths and settings.

    Each path in ``paths`` is deserialized into a string mapping or dictionary with string keys and yielded back to the
    caller.

    The current implementation only supports deserializing toml files.

    Args:
        paths: paths to data directories and files

    Yields:
        mutable string mapping

    """
    for path in parse_paths(paths):
        with path.open() as data_file:
            yield load_toml(data_file)


def load_package_schema(items: SomeConfigItems) -> MutableStrMapping:
    """Load a package's schema files.

    Args:
        items: unparsed schema package resources, paths, and fragments

    Returns:
        parsed package schema

    """
    return parse_schema(merge_config_data(items))


def merge_config_data(items: SomeConfigItems) -> MutableStrMapping:
    """Merge the package resources, data paths, and string mappings represented by items.

    All package resources will be assumed to be TOML files and will be treated as such.

    Args:
        items: package resources, files, or fragments to merge

    Returns:
        a mutable mapping

    """
    if isinstance(items, (dict, MappingProxyType, PackageResource, str, Path)):
        items = [items]

    configs: list[MutableStrMapping] = []
    for item in items:
        if isinstance(item, PackageResource):
            mappings = load_config_data(item.get_resource())
        elif isinstance(item, (str, Path)):
            mappings = load_config_data(collect_config_files(item, "*.toml"))
        else:
            mappings = [item]  # type: ignore
        configs += mappings
    return merge(configs)


def load_package_data(
    schema: StrMapping,
    data: SomeConfigItems | None = None,
    *,
    settings: Iterable[str] | None = None,
    interpolated: bool = True,
) -> MutableStrMapping:
    """Load data based on specified package resources, fragments, and settings.

    This function facilitates app usage of data mappings. Package resources, paths, and fragments (nested string
    mappings that may have missing values) are processed from left to right. The function outline is presented below:

    1. Build config data by merging parts from
       - default values from the schema
       - files
       - package resources
       - data fragments
       - settings
    2. Optionally interpolate the config data

    Args:
        schema: parsed schema mapping
        data: data package resources, files, and fragments
        settings: nested key assignment expressions of the form key[.key]*=value used to override data values
        interpolated: If True interpolate the merged data mapping

    Returns:
        (optionally interpolated) data that needs to validated, normalized, and (optionally) frozen

    """
    data = merge_config_data(data or ())

    data = merge(
        [
            build_data_from_schema(schema),
            apply_schema(schema, data),
        ]
    )

    for setting in settings or ():
        set_expr(schema, data, setting)

    if interpolated:
        data = interpolate(data)

    return data


def set_expr(schema: StrMapping, data: MutableStrMapping, expr: str) -> None:
    """Update data using an assignment expression.

    set_expr does the following:
        1. Splits expr (e.g. key1.key2.key2=value) into a nested key (e.g. key1.key2.key) and a value
        2. Lookups the schema attribute using the nested key (e.g. schema[key1][key2][key3]))
        3. Uses the schema attribute to convert value into the proper type
        4. Overwrites the old value in data to the new converted value

    Note:
        Any missing keys or levels in ``data`` are automatically added.

    Args:
        schema: nested dictionary with string keys of schema attributes
        data: nested dictionary with string keys
        expr: assignment expression of the form key[.key]*=value

    Raises:
        ValueError: an invalid assignment expression was passed

    Examples:
        >>> from .types import IntSchemaType
        >>> data = {}
        >>> schema = {'a': {'b': IntSchemaType()}}
        >>> set_expr(schema, data, 'a.b=-1')
        >>> data
        {'a': {'b': -1}}
        >>> set_expr(schema, data, 'a.b=')
        >>> data
        {'a': {'b': None}}

    """
    if "=" not in expr:
        raise ValueError(f"invalid assignment expression: {expr}")

    path, str_value = expr.split("=", 1)
    keys = path.split(".")

    for key in keys[:-1]:
        schema = schema[key]

    stype = schema[keys[-1]]
    value = stype.parse_string(str_value)

    if value is None and stype.required:
        raise ValueError(f"{path} is required")

    for key in keys[:-1]:
        data.setdefault(key, {})
        data = data[key]

    data[keys[-1]] = value
