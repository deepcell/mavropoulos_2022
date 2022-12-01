"""Interpolation."""

from __future__ import annotations

import re

from typing import Iterable


# Regular expression for interpolation variables
interp_re = re.compile(r"{([A-Za-z0-9_.-]+)}")


def _parse_expression(expr: str) -> list[str]:
    """Parse a nested key expression into a list of nested keys.

    Args:
        expr: nested key expression of the form key[.key]*

    Returns:
        list of nested keys

    Examples:
        >>> _parse_expression('a')
        ['a']

        >>> _parse_expression('a.b')
        ['a', 'b']

        >>> _parse_expression('..a.')
        ['', '', 'a', '']

    """
    return expr.split(".")


def _join_keys(keys: Iterable[str]) -> str:
    """Return a nested key expression from a list of nested keys.

    Args:
        keys: list of nested keys

    Returns:
        nested key expression of the form key[.key]*

    """
    return ".".join(keys)
