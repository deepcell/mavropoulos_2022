"""Test rowcontainer.py."""

from __future__ import annotations

from contextlib import ExitStack as does_not_raise
from dataclasses import dataclass
from typing import Any, Iterable

from pytest import mark, raises

from askcell.utils.rowcontainer import RowContainer


@dataclass
class Parent:
    """Simple test dataclass."""

    a: int
    b: str | None


@dataclass
class Child(Parent):
    """Test dataclass for isinstance checking."""

    c: str


@mark.parametrize(
    "row_type, rows, expected, exception",
    [
        # No error is raised, no rows are given
        (Parent, None, [], does_not_raise()),
        # No error is raised, correct rows are returned
        (Parent, [Parent(1, "x"), Child(2, "y", "z")], [Parent(1, "x"), Child(2, "y", "z")], does_not_raise()),
        (Parent, iter([Parent(1, "x")]), [Parent(1, "x")], does_not_raise()),
        # Raises a TypeError because the row type Child is not an instance of Parent
        (Child, [Parent(1, "x")], None, raises(TypeError)),
    ],
)
def test_init(
    row_type: Parent | Child,
    rows: Iterable[Parent | Child] | None,
    expected: Iterable[Parent | Child] | None,
    exception: does_not_raise,
) -> None:
    """Test constructor."""
    with exception:
        # Ignoring error as we expect mypy to complain
        assert RowContainer(row_type, rows).rows == expected  # type: ignore


@mark.parametrize(
    "object1, object2, expected",
    [
        (RowContainer(Parent), object(), False),
        (object(), RowContainer(Parent), False),
        (RowContainer(Parent), RowContainer(Child), False),
        (RowContainer(Parent), RowContainer(Parent), True),
        (RowContainer(Parent, [Parent(1, "x")]), RowContainer(Parent), False),
        (RowContainer(Parent, [Parent(1, "x")]), RowContainer(Parent, [Parent(2, "x")]), False),
        (RowContainer(Parent, [Parent(1, "x")]), RowContainer(Parent, [Parent(1, "x")]), True),
        (
            RowContainer(Parent, [Parent(1, "x"), Parent(2, "x")]),
            RowContainer(Parent, [Parent(2, "x"), Parent(1, "x")]),
            False,
        ),
        (
            RowContainer(Parent, [Parent(1, "x"), Parent(2, "x")]),
            RowContainer(Parent, [Parent(1, "x"), Parent(2, "x")]),
            True,
        ),
    ],
)
def test_eq(object1: Any, object2: Any, expected: bool) -> None:
    """Test equality."""
    assert (object1 == object2) is expected


def test_iter() -> None:
    """Test row iterator."""
    parent = Parent(1, "x")
    row_iter = iter(RowContainer(Parent, [parent]))
    assert parent is next(row_iter)

    with raises(StopIteration):
        next(row_iter)


def test_getitem() -> None:
    """Test getting a row."""
    parent = Parent(1, "x")
    rows = RowContainer(Parent, [parent])
    assert rows[0] is parent


def test_len() -> None:
    """Test getting number of rows."""
    rows = RowContainer(Parent, [Parent(1, "x")])
    assert len(rows) == 1


@mark.parametrize(
    "row_type, exception",
    [
        (Parent, does_not_raise()),
        (Child, raises(TypeError)),
    ],
)
def test_append(row_type: Any, exception: Any) -> None:
    """Test appending to internal backing store."""
    rows: RowContainer[Any] = RowContainer(row_type)
    parent = Parent(1, "x")
    with exception:
        rows.append(parent)
        assert rows[-1] is parent


@mark.parametrize(
    "row_type, exception",
    [
        (Parent, does_not_raise()),
        (Child, raises(TypeError)),
    ],
)
def test_extend(row_type: Any, exception: Any) -> None:
    """Test extending to internal backing store."""
    rows: RowContainer[Any] = RowContainer(row_type)
    parent = Parent(1, "x")
    with exception:
        rows.extend([parent])
        assert rows[-1] is parent


@mark.parametrize(
    "kwargs, expected",
    [
        ({}, [Parent(1, "x"), Child(1, "x", "z"), Child(2, "y", "z")]),
        ({"a": 2}, [Child(2, "y", "z")]),
        ({"a": 1, "b": "x"}, [Parent(1, "x"), Child(1, "x", "z")]),
    ],
)
def test_subset(kwargs: dict[str, int], expected: list[Parent | Child]) -> None:
    """Test subset."""
    rows = RowContainer(Parent, [Parent(1, "x"), Child(1, "x", "z"), Child(2, "y", "z")])
    assert list(rows.subset(**kwargs)) == expected


@mark.parametrize(
    "args, expected",
    [
        # Sort by a
        ("a", [Parent(1, "x"), Child(1, "x", "z"), Child(2, "y", "z"), Child(3, "y", "z"), Child(3, "a", "z")]),
    ],
)
def test_sorting(args: str, expected: list[Parent | Child]) -> None:
    """Test subset."""
    rows = RowContainer(
        Parent,
        [
            Parent(1, "x"),
            Child(1, "x", "z"),
            Child(3, "y", "z"),
            Child(2, "y", "z"),
            Child(3, "a", "z"),
        ],
    )
    rows.sort(args)
    assert list(rows) == expected


@mark.parametrize(
    "kwargs, expected",
    [
        ({}, 3),
        ({"a": 2}, 1),
        ({"a": 1, "b": "x"}, 2),
    ],
)
def test_count(kwargs: dict[str, int], expected: int) -> None:
    """Test count."""
    rows = RowContainer(Parent, [Parent(1, "x"), Child(1, "x", "z"), Child(2, "y", "z")])
    assert rows.count(**kwargs) == expected


@mark.parametrize(
    "kwargs, expected",
    [
        ({}, [Parent(1, "x"), Child(1, "x", "z"), Child(2, "y", "z")]),
        ({"a": 2}, [Child(2, "y", "z")]),
        ({"a": 1, "b": "x"}, [Parent(1, "x"), Child(1, "x", "z")]),
    ],
)
def test_find(kwargs: dict[str, int], expected: list[Parent | Child]) -> None:
    """Test querying."""
    rows = RowContainer(Parent, [Parent(1, "x"), Child(1, "x", "z"), Child(2, "y", "z")])
    assert list(rows.find(**kwargs)) == expected


def test_find_one() -> None:
    """Test exact querying."""
    child = Child(2, "y", "z")
    rows = RowContainer(Parent, [Parent(1, "x"), Child(1, "x", "z"), child])
    assert rows.find_one(a=2) is child


@mark.parametrize(
    "args, unique, expected, exception",
    [
        (["b"], False, {"x": [Parent(1, "x"), Parent(2, "x")]}, does_not_raise()),
        (["b"], True, None, raises(KeyError)),
        (["a", "b"], True, {1: {"x": Parent(1, "x")}, 2: {"x": Parent(2, "x")}}, does_not_raise()),
    ],
)
def test_build_index(args: list[str], unique: bool, expected: Any, exception: Any) -> None:
    """Test building indices."""
    rows = RowContainer(Parent, [Parent(1, "x"), Parent(2, "x")])
    with exception:
        assert rows.build_index(*args, unique=unique) == expected


@mark.parametrize(
    "attrs, ignore_missing, expected",
    [
        (("a",), False, {1, 2}),
        (("b",), False, {None, "x"}),
        (("b",), True, {"x"}),
        (("a", "b"), False, {(1, "x"), (2, None)}),
        (("a", "b"), True, {(1, "x"), (2, None)}),
        (("b", "a"), False, {("x", 1), (None, 2)}),
        (("b", "a"), True, {("x", 1), (None, 2)}),
        (("a", "a"), False, {(1, 1), (2, 2)}),
        (("a", "a"), True, {(1, 1), (2, 2)}),
        (("b", "b"), False, {("x", "x"), (None, None)}),
        (("b", "b"), True, {("x", "x")}),
    ],
)
def test_attr_values(
    attrs: tuple[str, ...],
    ignore_missing: bool,
    expected: set[int | str | tuple[int | str | None] | None],
) -> None:
    """Test getting attribute values across all rows."""
    rows = RowContainer(Parent, [Parent(1, "x"), Parent(2, None)])
    assert rows.attr_values(*attrs, ignore_missing=ignore_missing) == expected
