"""Test askcell.utils.merge."""
from __future__ import annotations

from datetime import date, datetime
from typing import Any

from pytest import mark, raises

from askcell.utils.merge import copy, merge
from askcell.utils.pathutils import Path


DATA = {
    "bool": True,
    "int": 0,
    "float": 1.2,
    "date": date(2022, 6, 2),
    "datetime": datetime(2022, 6, 3),
    "string": "{{escaped braces}}",
    "path": Path("/tmp"),
    "none": None,
    "dict": {
        "key1": "int equals {..int}",
        "key2": "{int}",
        "key3": "{.key2}",
    },
    "list": [1, 2, 3],
    "tuple": (4, 5, 6),
    "convoluted_list": [
        ("tuple element interpolation", "{float}"),
        ["list element interpolation", "Interpolating {bool} and returning a string"],
        "{dict}",
        "{path}",
    ],
    "convoluted_tuple": (
        ("tuple element interpolation", "Interpolating {date} and returning a string"),
        ["list element interpolation", "{datetime}"],
        "{dict}",
        Path("{path}"),
    ),
}


def test_copy() -> None:
    """Test copying a nested structure."""
    output = copy(DATA)
    assert output is not DATA
    assert output == {
        "bool": True,
        "int": 0,
        "float": 1.2,
        "date": date(2022, 6, 2),
        "datetime": datetime(2022, 6, 3),
        "string": "{{escaped braces}}",
        "none": None,
        "dict": {
            "key1": "int equals {..int}",
            "key2": "{int}",
            "key3": "{.key2}",
        },
        "path": Path("/tmp"),
        "list": [1, 2, 3],
        "tuple": [4, 5, 6],
        "convoluted_list": [
            ["tuple element interpolation", "{float}"],
            ["list element interpolation", "Interpolating {bool} and returning a string"],
            "{dict}",
            "{path}",
        ],
        "convoluted_tuple": [
            ["tuple element interpolation", "Interpolating {date} and returning a string"],
            ["list element interpolation", "{datetime}"],
            "{dict}",
            Path("{path}"),
        ],
    }


@mark.parametrize(
    "data, exception, message",
    [
        (object(), raises(TypeError), "data must be a mapping type"),
        ({"key": object()}, raises(ValueError), f"Unsupported value type: {type(object())}"),
    ],
)
def test_copy_unsupported_type(data: Any, exception: Any, message: str) -> None:
    """Test that an error is raised when copying an unsupported value type."""
    with exception as excinfo:
        copy(data)
    assert message in str(excinfo)


def test_copy_duplicate() -> None:
    """Test copy memoization."""
    list1 = [1, 2, "{int}"]
    dict1 = {"a": "{..int}"}

    data = {
        "int": 0,
        "list1": list1,
        "list2": list1,
        "dict1": dict1,
        "dict2": dict1,
        "dict3": {
            "list3": list1,
            "dict4": dict1,
        },
    }

    output = copy(data)

    assert output == {
        "int": 0,
        "list1": [1, 2, "{int}"],
        "list2": [1, 2, "{int}"],
        "dict1": {"a": "{..int}"},
        "dict2": {"a": "{..int}"},
        "dict3": {
            "list3": [1, 2, "{int}"],
            "dict4": {"a": "{..int}"},
        },
    }
    assert output["list2"] is output["list1"]
    assert output["dict2"] is output["dict1"]
    assert output["dict3"]["list3"] is output["list1"]
    assert output["dict3"]["dict4"] is output["dict1"]


def test_copy_recursive1() -> None:
    """Test that self-recursive data structures can be copied."""
    orig: dict[str, Any] = {}
    orig["key"] = orig

    new = copy(orig)

    assert new.keys() == orig.keys()
    assert new["key"] is new


def test_copy_recursive2() -> None:
    """Test that self-recursive data structures can be copied."""
    orig: dict[str, Any] = {}
    orig["key"] = [orig]

    new = copy(orig)

    assert new.keys() == orig.keys()
    assert len(new["key"]) == 1 and new["key"][0] is new


def test_merge1() -> None:
    """Test merge with overlapping and overriding values with some immutable elements returns a mutable mapping."""
    data1 = {
        "bool": True,
        "dict": {
            "a": 1,
            "b": 2,
            "z": 10,
        },
        "list": [1, 2, 3],
    }
    data2 = {
        "dict": {
            "x": 9,
            "z": 10,
            "c": {"d": 4},
        },
        "list": (7, 8, 9),
    }
    data12 = merge([data1, data2])

    assert data12 == {
        "bool": True,
        "dict": {
            "a": 1,
            "b": 2,
            "x": 9,
            "z": 10,
            "c": {"d": 4},
        },
        "list": [7, 8, 9],
    }

    data3 = {
        "dict": {
            "b": (1, 2, 3),
            "y": 8,
            "c": {"e": 5},
        },
    }

    assert merge([data12, data3]) == {
        "bool": True,
        "dict": {
            "a": 1,
            "b": [1, 2, 3],
            "x": 9,
            "y": 8,
            "z": 10,
            "c": {
                "d": 4,
                "e": 5,
            },
        },
        "list": [7, 8, 9],
    }


def test_merge2() -> None:
    """Test that merge returns an empty dict given an empty iterable."""
    assert merge([]) == {}


def test_merge3() -> None:
    """Test merge on an iterable with a single element with some immutable elements returns a mutable mapping."""
    assert merge([DATA]) == copy(DATA)


def test_merge_recursive1() -> None:
    """Test that self-recursive data structures can be merged."""
    orig: dict[str, Any] = {}
    orig["key"] = orig

    new = merge([orig])

    assert new.keys() == orig.keys()
    assert new["key"] == new


def test_merge_recursive2() -> None:
    """Test that self-recursive data structures can be merged."""
    orig: dict[str, Any] = {}
    orig["key"] = [orig]

    new = merge([orig])

    assert new.keys() == orig.keys()
    assert len(new["key"]) == 1 and new["key"][0] == new
