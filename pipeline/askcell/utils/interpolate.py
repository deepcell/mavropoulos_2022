"""Interpolation."""

from __future__ import annotations

import re

from pathlib import Path
from types import MappingProxyType
from typing import Any

from .merge import copy
from .types import AtomTypeTuple, ImmutableStrMapping, MutableStrMapping


__all__ = ["interpolate", "interpolate_expr"]


# Maxumim number of levels of recursion when resolving interpolated variables
# This is to protect against circular dependencies.
MAX_RECURSION = 20


# Regular expression for interpolation variables
interp_re = re.compile(r"{([A-Za-z0-9_.-]+)}")


def interpolate(data: MutableStrMapping) -> MutableStrMapping:
    """Process a nested-mapping and interpolate all strings.

    Absolute interpolation uses ``{section.section.value}`` to refer a value in a foreign section starting from the
    root of the nested-mapping.  Zero or more sections may be specified, but all lookups are relative to the root of
    the mapping.  Relative interpolation is similar, except that expressions begin with one or more ``.`` characters.
    The first ``.`` indicates that references are relative to the current section, while each additional ``.`` traverses
    to the parent of the current section.  i.e.  ``.value`` refers to a value in the current section, ``..value`` one
    in the parent of the current section, ``...section.value`` a value within a section that is the grandparent of the
    current, and so on.  Sections and values are composed of ASCII letters, ASCII digits, underscores, and dashes
    ``(A-Za-z0-9_-)``.

    The-level algorithm requires two or three passes through the data:

      1. a mutable copy of the data is made
      2. then interpolation is performed  in-place

    Copy operation has to take precautions against trying to copy data structures with reference
    cycles, which uses caching to store the mutable copy of each
    non-atomic value.  This is because recursive mappings may have reference cycles and a naive copying algorithm could
    result in infinite recursion. For example::

        a = "{a}"

        [section1]
        foo="{section2.bar}"

        [section2]
        bar="{section1.foo}"

    Interpolation needs to be in-place in an already mutable data structure, since interpolated values can create new
    reference cycles within the data structure.  This is because interpolated values can be anything -- not just strings
    -- which can introduce additional reference cycles.

    Notes:
      * Expressions are not recursive and cannot contain other nested expressions.  This does not preclude recursive
        expression resolution, just that one cannot have a recursive expression ``{a.{b}.c}`` (an expression within an
        expression).
      * Section and value names in expressions may contain only ASCII numbers,
        letters, underscore ``(_)`` and hyphen ``(-)``.
      * Expressions may resolve to non-string values, including boolean, numeric, datetime, date, tuples, and lists,
        mapping types.
      * Incorrectly formatted expressions do not raise errors (e.g. ``'{a'``).  They are left as-is.
      * Path objects are always converted to strings.

    Args:
        data: recursive string mapping

    Returns:
        interpolated copy of data

    Examples:
        >>> data = {
        ...     'foo': 1,
        ...     'bar': 'a',
        ...     'baz': {
        ...         'a': '{foo}',
        ...         'b': 'Absolute reference to {bar}',
        ...         'c': 'Relative reference to {..foo}',
        ...     }
        ... }

        >>> interpolate(data)
        {'foo': 1, 'bar': 'a', 'baz': {'a': 1, 'b': 'Absolute reference to a', 'c': 'Relative reference to 1'}}

    """
    # Must copy data or else references will not work correctly
    data = copy(data)

    _interpolate_mapping(data)

    return data


def _interpolate_mapping(
    data: MutableStrMapping,
    path: list[ImmutableStrMapping] | None = None,
    memo: dict[int, Any] | None = None,
) -> None:
    """Process a nested string-mapping in-place and interpolate all strings.

    Args:
        data: recursive string mapping
        path: list of parent mappings
        memo: memoization dictionary

    Examples:
        >>> data = {
        ...     'a':  '{foo}',
        ...     'foo': 1,
        ... }

        >>> _interpolate_mapping(data)

        >>> data == {'a': 1, 'foo': 1}
        True

    """
    if memo is None:
        memo = {}

    mkey = id(data)
    if mkey in memo:
        return

    memo[mkey] = data

    if path is None:
        path = [data]

    for key, value in data.items():
        data[key] = _interpolate_value(value, path, memo)


def _interpolate_value(value: Any, path: list[ImmutableStrMapping], memo: dict[int, Any]) -> Any:
    """Interpolate a value.

    Args:
        value: value to interpolate
        path: list of nested dictionary with string keys
        memo: memoization dictionary

    Returns:
        interpolated value

    Examples:
        >>> data = {
        ...     'a':  '{foo}',
        ...     'foo': 1,
        ... }

        >>> _interpolate_value('{foo}', [data], {})
        1

    """
    if isinstance(value, dict):
        _interpolate_mapping(value, path + [value], memo)
    elif isinstance(value, (str, Path)):
        value = _interpolate_string(str(value), path, memo)
    elif isinstance(value, list):
        mkey = id(value)
        if mkey not in memo:
            memo[mkey] = value
            value[:] = (_interpolate_value(item, path, memo) for item in value)
    elif not isinstance(value, AtomTypeTuple):
        raise ValueError(f"Unsupported value type: {type(value)}")

    return value


def _interpolate_string(s: str, path: list[ImmutableStrMapping], memo: dict[int, Any], level: int = 0) -> Any:
    """Interpolate a string value.

    Args:
        s: string to resolve
        path: list of nested dictionary with string keys
        memo: memoization dictionary
        level: recursion level

    Returns:
        interpolated string or value

    Examples:
          >>> _interpolate_string('no interpolation', [], {})
          'no interpolation'

          >>> _interpolate_string('a is {a}', [{'a': 1}], {})
          'a is 1'

    """
    if level >= MAX_RECURSION:
        raise ValueError(f"Maximum recursion limit exceeded: {MAX_RECURSION}")

    if not isinstance(s, str):
        raise TypeError("cannot interpolate non-string values")

    # Replace escaped braces with innocuous characters
    s = s.replace("{{", "\1").replace("}}", "\2")

    for match in reversed(list(interp_re.finditer(s))):
        start = match.start()
        stop = match.end()
        expr = match[1]
        value = _resolve_path(expr, path, memo, level)

        # Matches that span the entire string are returned as-is without conversion to string
        if not start and stop == len(s):
            if isinstance(value, (list, tuple)):
                raise ValueError("Cannot interpolate sequence types as values")
            return value

        s = f"{s[:start]}{value}{s[stop:]}"

    return s.replace("\1", "{").replace("\2", "}")


def _resolve_path(expr: str, path: list[ImmutableStrMapping], memo: dict[int, Any], level: int = 0) -> Any:
    """Lookup an interpolation expression relative to the current path.

    Args:
        expr: expression to resolve
        path: list of nested dictionary with string keys
        memo: memoization dictionary
        level: recursion level

    Returns:
        value of resolved expression

    Examples:
        >>> _resolve_path('a', [{'a': 1}], {})
        1

    """
    if level >= MAX_RECURSION:
        raise ValueError(f"Maximum recursion limit exceeded: {MAX_RECURSION}")

    if not expr:
        raise ValueError("no value to resolve")

    parts = expr.split(".")

    if parts[0]:
        # Absolute paths start at the root
        path = path[:1]
    else:
        # Relative paths start at the end of the current path
        path = path.copy()
        parts.pop(0)

    for key in parts[:-1]:
        if not key:
            if len(path) > 1:
                path.pop()
        else:
            pos = path[-1][key]
            if not isinstance(pos, (dict, MappingProxyType)):
                raise ValueError(f"Interpolation into non-dictionary: {key}")
            path += [pos]

    value = path[-1][parts[-1]]

    if isinstance(value, (str, Path)):
        value = _interpolate_string(str(value), path, memo, level + 1)

    return value


def interpolate_expr(data: MutableStrMapping, expr: str) -> Any:
    """Resolve an interpolation expression in data.

    Args:
        data: nested dictionary with string keys
        expr: expression to resolve

    Returns:
        value of resolved expression

    Examples:
        >>> interpolate_expr({'a': 1}, 'a')
        1

    """
    return _resolve_path(expr, [data], {})
