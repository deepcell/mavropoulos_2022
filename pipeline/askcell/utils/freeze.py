"""Utilities to facilitate managing string-keyed nested-mapping data structures."""

from __future__ import annotations

from types import MappingProxyType

from .types import AtomTypeTuple, ImmutableStrMapping, MemoCache, StrMapping, Value


__all__ = ["freeze"]


def freeze(data: StrMapping, *, memo: MemoCache | None = None) -> ImmutableStrMapping:
    """Copy and make immutable a recursive string mapping.

    Args:
        data: nested dictionary with string keys
        memo: memoization dictionary.  This is an implementation detail that users should not specify.

    Returns:
        immutable deep copy of data

    Examples:
        >>> freeze({'a': [1, 2]})
        mappingproxy({'a': (1, 2)})

    """
    if not isinstance(data, (dict, MappingProxyType)):
        raise TypeError("data must be a mapping type")

    if memo is None:
        memo = {}

    mkey = id(data)
    result = memo.get(mkey)

    if result is None:
        # Use True as a placeholder in the memo cache for a non-atom that is being computed
        memo[mkey] = True
        result = memo[mkey] = MappingProxyType({key: _freeze_value(value, memo) for key, value in data.items()})
    elif result is True:
        raise ValueError("cannot freeze cyclic data structures")

    return result


def _freeze_value(value: Value, memo: MemoCache) -> Value:
    """Freeze a value (make an immutable copy).

    This function is not intended as a public API. Rather, use :func:`freeze`.

    Args:
        value: an atomic value (bool, int, str, float, date, datetime) or a non-atomic value (list, tuple).
        memo: memoization dictionary.

    Returns:
        an immutable copy of value

    """
    # Return atomic value (already immutable)
    if isinstance(value, AtomTypeTuple):
        return value

    mkey = id(value)
    result = memo.get(mkey)

    # Return already frozen memoized value
    if result is None:
        if isinstance(value, (dict, MappingProxyType)):
            # freeze takes care of memoizing mappings
            result = freeze(value, memo=memo)
        elif isinstance(value, (list, tuple)):
            # Use True as a placeholder for a non-atom that is being computed
            memo[mkey] = True
            result = memo[mkey] = tuple(_freeze_value(v, memo) for v in value)
        else:
            raise ValueError(f"Unsupported value type: {type(value)}")

    elif result is True:
        raise ValueError("cannot freeze cyclic data structures")

    return result
