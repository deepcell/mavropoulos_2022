"""Utilities to facilitate managing string-keyed nested-mapping data structures."""

from __future__ import annotations

from types import MappingProxyType
from typing import Any, Iterable

from askcell.utils.types import AtomTypeTuple

from .types import MutableStrMapping, StrMapping


__all__ = ["merge", "copy"]


def _merge2(a: MutableStrMapping, b: StrMapping, memo: dict[int, Any]) -> None:
    """Merge two nested mappings.

    When merging values from b into a, there are a number of situations to consider. First, given that a must remain
    a mutable nested mapping, any immutable collections in b such as MappingProxyType's and tuples must be made mutable.
    This is achieved by calling :func:`_copy_value` on all values of b before they are merged back to a.

    With the above in mind, the overlap between the keys of a and b are examined. For keys in b but not in a (i.e.
    b.keys() - a.keys()), the corresponding key-value pairs from b are merged to a. For keys common to both a and b
    (i.e. a.keys() intersect b.keys()), further examination is needed for the values of a and b at each key. If both a
    and b values correspond to mapping types, _merge2 is recursively called on those values.
    Otherwise, if a and b values differ, b's value overrides a's value, i.e. a[key]=b[key].

    This function is not intended as a public API. Rather, use :func:`merge`.

    Args:
        a: mutable nested dictionary with string keys
        b: nested dictionary with string keys. Overrides values set in a.
        memo: memoization dictionary

    Examples:
        >>> data1 = {'a': {'b': 1, 'd': {'will be': 'overwritten'}}}
        >>> data2 = {'a': MappingProxyType({'c': 2, 'd': 0})}
        >>> _merge2(data1, data2, {})
        >>> data1 == {'a': {'b': 1, 'd': 0, 'c': 2}}
        True

    """
    for key, value_b in b.items():
        if key in a:
            value_a = a[key]
            if isinstance(value_a, dict) and isinstance(value_b, (dict, MappingProxyType)):
                _merge2(value_a, value_b, memo)
            elif value_a != value_b:
                # Override current value
                a[key] = _copy_value(value_b, memo)
        else:
            # Merge in new value
            a[key] = _copy_value(value_b, memo)


def merge(data: Iterable[StrMapping]) -> MutableStrMapping:
    """Merge sequence of mappings.

    If data is empty, merge returns an empty dictionary. Else, _merge2 is cumulatively applied to items of data from
    left to right starting with an empty dictionary. This reduces data to a single merged mapping that is then returned.
    For example, if data equals [a, b, c, d] where a, b, c, and d are mappings; merge(data) will call
    _merge2(_merge2(_merge2(_merge2({}, a), b), c), d).

    It should be noted that merge will make immutable collections mutable where tuples are converted to lists and
    MappingProxyType's become regular dictionaries.

    Args:
        data: iterable of nested dictionary with string keys

    Returns:
        nested mapping that combines all items in data

    Examples:
        >>> data1 = {'a': {'b': 1, 'd': {'will be': 'overwritten'}}}
        >>> data2 = {'a': {'c': 2, 'd': 0}}
        >>> merge([data1, data2])
        {'a': {'b': 1, 'd': 0, 'c': 2}}

    """
    result: MutableStrMapping = {}
    memo: dict[int, Any] = {}
    for item in data:
        _merge2(result, item, memo)

    return result


def _copy_value(value: Any, memo: dict[int, Any]) -> Any:
    """Copy a value (making a mutable copy of non-atomic values).

    This function is not intended as a public API. Rather, use :func:`copy`.

    Args:
        value: an atomic value (bool, int, str, float, date, datetime) or a non-atomic value (list, tuple).
        memo: memoization dictionary.

    Returns:
        value unchanged if an atom, a mutable copy of value if not

    >>> _copy_value((1, 2), {})
    [1, 2]

    """
    # Return atomic values (since they are immutable)
    if isinstance(value, AtomTypeTuple):
        return value

    mkey = id(value)
    result: Any = memo.get(mkey)

    # Return already copied memoized value
    if result is not None:
        return result

    # Copy new non-atomic value
    if isinstance(value, (dict, MappingProxyType)):
        # copy takes care of memoizing mappings
        result = copy(value, memo=memo)
    elif isinstance(value, (list, tuple)):
        # Must memoize empty result before filling in values (or else we'll keep trying to copy it)
        result = memo[mkey] = []
        result[:] = (_copy_value(v, memo) for v in value)
    else:
        raise ValueError(f"Unsupported value type: {type(value)}")

    return result


def copy(data: StrMapping, *, memo: dict[int, Any] | None = None) -> MutableStrMapping:
    """Deeply copy a nested string-mapping, making non-atomic values mutable when possible.

    Args:
        data: nested dictionary with string keys
        memo: memoization dictionary.  This is an implementation detail that users should not specify.

    Returns:
        deep copy of data

    Examples:
        >>> data = {
        ...     'a': (1, 2),
        ...     'b': [3, 4],
        ... }

        >>> copy(data)
        {'a': [1, 2], 'b': [3, 4]}

    """
    if not isinstance(data, (dict, MappingProxyType)):
        raise TypeError("data must be a mapping type")

    if memo is None:
        memo = {}

    mkey = id(data)
    result = memo.get(mkey)

    if result is None:
        # Must memoize empty result before filling in values (or else we'll keep trying to copy it)
        result = memo[mkey] = {}
        result.update({key: _copy_value(value, memo) for key, value in data.items()})

    return result
