"""List-backed container storing dataclasses."""

from __future__ import annotations

from operator import attrgetter
from typing import Any, Generic, Iterable, Iterator, TypeVar

from .. import only_one


__all__ = ["RowContainer"]


T = TypeVar("T")
GenericRowContainer = TypeVar("GenericRowContainer", bound="RowContainer[Any]")  # TypeVar can't refer to Type variables


class RowContainer(Generic[T]):
    """Dataclass container of rows with list-like and query methods."""

    def __init__(self, row_type: type[T], rows: Iterable[T] | None = None) -> None:
        """Save dataclass type, validate passed rows, and create backing store using passed rows.

        Args:
            row_type: dataclass object type
            rows: iterable of dataclass objects

        """
        self.row_type = row_type

        # Fastpath: no rows were passed
        if not rows:
            self.rows: list[T] = []
            return

        # Slowpath: ensure list and also validate elements of rows are instances of self.row_type
        if not isinstance(rows, list):
            rows = list(rows)

        for row in rows:
            if not isinstance(row, row_type):
                raise TypeError(f"Row is of a type {type(row)} and is not an instance of {self.row_type}.")

        self.rows = rows

    def __eq__(self, other: object) -> bool:
        """Return True if self and other are equal."""
        if not isinstance(other, RowContainer):
            return NotImplemented
        return self.row_type == other.row_type and self.rows == other.rows

    def __iter__(self) -> Iterator[T]:
        """Iterate over rows.

        Returns:
            iterator over rows

        """
        return iter(self.rows)

    def __getitem__(self, i: int) -> T:
        """Return i'th row.

        Args:
            i: index of row object to return

        Returns:
            i'th row

        """
        return self.rows[i]

    def __len__(self) -> int:
        """Return the number of rows.

        Returns:
            number of rows

        """
        return len(self.rows)

    def append(self, row: T) -> None:
        """Append a new row.

        Args:
            row: row to append

        Raises:
            TypeError: Raised if not isinstance(row, self.row_type)

        """
        if not isinstance(row, self.row_type):
            raise TypeError(f"Row is of a type {type(row)} and is not an instance of {self.row_type}.")

        self.rows.append(row)

    def extend(self, rows: Iterable[T]) -> None:
        """Extend contents with additional rows.

        Args:
            row: row to append

        Raises:
            TypeError: Raised if not isinstance(row, self.row_type)

        """
        for row in rows:
            self.append(row)

    def count(self, **kwargs: Any) -> int:
        r"""Count all rows that match keyword arguments.

        Args:
            \*\*kwargs: keyword arguments as selection criteria for rows.

        Returns:
            count of matched rows

        """
        return sum(1 for _ in self.find(**kwargs))

    def find(self, **kwargs: Any) -> Iterator[T]:
        r"""Find all rows that match keyword arguments.

        Args:
            \*\*kwargs: keyword arguments as selection criteria for rows.

        Returns:
            iterator over matching requests

        """
        # Fastpath: no keyword arguments to filter
        if not kwargs:
            return iter(self)

        # Slowpath: one or more keyword arguments to filter
        attr_get = attrgetter(*kwargs)
        expected_values = tuple(kwargs.values())

        # Match cardinality adjustment done by attrgetter
        if len(expected_values) == 1:
            expected_values = expected_values[0]

        return (row for row in self if expected_values == attr_get(row))

    def find_one(self, **kwargs: Any) -> T:
        r"""Find the unique row object with attributes that match the keyword arguments. Raise exception otherwise.

        Args:
            \*\*kwargs: keyword arguments as selection criteria for rows.

        Returns:
            matching row

        Raises:
            ValueError: exception if found rows is not == 1

        """
        return only_one(self.find(**kwargs))

    def subset(self: GenericRowContainer, **kwargs: Any) -> GenericRowContainer:
        r"""Return a new RowContainer with a subset of rows with attributes that match the keyword arguments.

        Args:
            \*\*kwargs: keyword arguments as selection criteria for rows.

        Returns:
            iterator over matching requests

        """
        return type(self)(self.row_type, list(self.find(**kwargs)))

    def sort(self, *args: str) -> None:
        r"""Sort rows based on attributes in \*args.

        Args:
            args: column names to use as keys to sort

        """
        self.rows.sort(key=attrgetter(*args))

    def build_index(self, *args: str, unique: bool = False) -> dict[str, Any]:
        r"""Generate a multi-level nested index based on attributes in \*args.

        Args:
            args: column names to use as index dictionary keys, ordered by nesting level
            unique: When unique is False, the default, then the bottom level of the mapping is a list of rows. When
                unique is True, at most one row may have the same values as specified in the keys. If so, the bottom
                level mapping is a single row. If two or more rows share the same key values,
                then a KeyError is raised.

        Returns:
            nested dictionary(s) with unique row value keys from input arguments. If unique is False, then the bottom
                level is a list of rows, otherwise the bottom level is a single row.

        Raises:
            KeyError: If unique is True and two or more rows share the same key values.

        """  # noqa: RST301
        singleton = len(args) == 1
        get_keys = attrgetter(*args)
        results: dict[str, Any] = {}

        for row in self.rows:
            key_attrs = get_keys(row)
            keys = [key_attrs] if singleton else key_attrs  # Make mypy happy

            # For all but the last key, create nested dictionaries
            r = results
            for key in keys[:-1]:
                r = r.setdefault(key, {})

            # The last level is always a list of rows
            last_key = keys[-1]
            if not unique:
                r.setdefault(last_key, []).append(row)
            elif last_key not in r:
                r[last_key] = row
            else:
                kvs = ", ".join("=".join(kv) for kv in zip(args, keys))
                raise KeyError(f"Non-unique result for index on key(s): {kvs}")

        return results

    def attr_values(self, *args: str, ignore_missing: bool = True) -> set[Any]:
        """Find all unique row values for args.

        If one attribute is requested, then a set of attribute values is returned.  If more than one attribute is
        requested, then a tuple of attribute values is returned.

        Args:
            args: row attribute names
            ignore_missing: if true, ignore missing values (i.e. Nones).  Default is True.

        Returns:
            unique values found in args for all data rows

        """
        get_values = attrgetter(*args)
        values = {get_values(row) for row in self.rows}

        if ignore_missing:
            values.discard(None if len(args) == 1 else (None,) * len(args))  # type: ignore

        return values
