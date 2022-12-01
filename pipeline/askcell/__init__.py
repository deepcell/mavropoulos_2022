"""Universally useful definitions and functions for askCell."""

from __future__ import annotations

from typing import Iterable, TypeVar, overload

from .utils.pathutils import (
    Path,
    PathLike,
    Pathy,
    PosixPath,
    SomePathLikes,
    WindowsPath,
    parse_path,
    parse_paths,
    set_path_readonly,
    validate_path,
    validate_paths,
)


__all__ = [
    # Export symbols defined in this module
    "Nothing",  # Singleton instance to represent Nothing (like None, but even less so)
    "only_one",  # Return the one and only item in iterator.
    "first",  # Return the one and only item in iterator.
    # path-related: re-export for convenience
    "validate_path",  # validate path instance based on flags
    "validate_paths",  # validate path instances based on flags
    "parse_path",  # convert a string or other path-like instance into a concrete Path instance
    "parse_paths",  # convert zero or more string or other path-like instances and yield concrete Path instance
    "set_path_readonly",  # sets permission of a file to read-only
    "Path",  # re-export for convenience
    "PosixPath",  # re-export for convenience
    "WindowsPath",  # re-export for convenience
    "Pathy",  # re-export for convenience
    "PathLike",  # Type for a single path-like instance (str or Path)
    "SomePathLikes",  # Type for zero or more path-like instances
]


class NothingType:
    """Class to represent true Nothing.

    To be used when None is a valid value and there is a need to represent Nothing.

    """


# Singleton instance to represent Nothing (like None, but even less so)
Nothing = NothingType()


T = TypeVar("T")
D = TypeVar("D")


@overload
def only_one(items: Iterable[T], *, default: NothingType = Nothing) -> T:
    """Explain to the type checker that calling with no default will return an item."""


@overload
def only_one(items: Iterable[T], *, default: D) -> T | D:
    """Explain to the type checker that calling with a default will return either an item or the default."""


def only_one(items: Iterable[T], *, default: NothingType | D = Nothing) -> T | D:
    """Return the one and only item in iterator.

    Args:
        items:   iterable
        default: default value to return if items is empty

    Returns:
        value in items or default value if specified and items is empty


    Raises:
        ValueError: items was empty and default was not specified or items contained more than one item

    Examples:
        >>> only_one([1])
        1
        >>> only_one((1,))
        1
        >>> only_one({1: 'a'})
        1
        >>> only_one(range(10, 11))
        10
        >>> only_one([], default=None)
        >>> only_one((), default=None)
        >>> only_one({}, default=None)
        >>> only_one(range(10, 10), default=None)
        >>> only_one([])
        Traceback (most recent call last):
            ...
        ValueError: items was empty
        >>> only_one(())
        Traceback (most recent call last):
            ...
        ValueError: items was empty
        >>> only_one({})
        Traceback (most recent call last):
            ...
        ValueError: items was empty
        >>> only_one(range(10, 10))
        Traceback (most recent call last):
            ...
        ValueError: items was empty
        >>> only_one([1, 2])
        Traceback (most recent call last):
            ...
        ValueError: items contained more than one value
        >>> only_one((1, 2))
        Traceback (most recent call last):
            ...
        ValueError: items contained more than one value
        >>> only_one({1: 'a', 2: 'b'})
        Traceback (most recent call last):
            ...
        ValueError: items contained more than one value
        >>> only_one(range(10, 12))
        Traceback (most recent call last):
            ...
        ValueError: items contained more than one value

    """
    items = iter(items)

    try:
        item = next(items)
    except StopIteration:
        if not isinstance(default, NothingType):
            return default
        raise ValueError("items was empty")

    # Try to get second result (supposed to fail)
    try:
        next(items)
        raise ValueError("items contained more than one value")
    except StopIteration:
        pass

    return item


@overload
def first(items: Iterable[T], *, default: NothingType = Nothing) -> T:
    """Explain to the type checker that calling with no default will return an item."""


@overload
def first(items: Iterable[T], *, default: D) -> T | D:
    """Explain to the type checker that calling with a default will return either an item or the default."""


def first(items: Iterable[T], *, default: NothingType | D = Nothing) -> T | D:
    """Return the one and only item in iterator.

    Args:
        items:   iterable
        default: default value to return if items is empty

    Returns:
        value in items or default value if specified and items is empty

    Raises:
        ValueError: items was empty and default was not specified or items contained more than one item

    Examples:
        >>> first([1])
        1
        >>> first([1, 2])
        1
        >>> first((1,))
        1
        >>> first({1: 'a', 2: 'b'})
        1
        >>> first(range(10, 11))
        10
        >>> first(range(10, 15))
        10
        >>> first([], default=None)
        >>> first((), default=None)
        >>> first({}, default=None)
        >>> first(range(10, 10), default=None)
        >>> first([])
        Traceback (most recent call last):
            ...
        ValueError: items was empty

    """
    items = iter(items)

    try:
        return next(items)
    except StopIteration:
        if not isinstance(default, NothingType):
            return default
        raise ValueError("items was empty")
