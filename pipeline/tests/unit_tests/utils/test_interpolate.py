"""Test askcell.utils.interpolate."""

from __future__ import annotations

from contextlib import nullcontext
from typing import Any

from pytest import mark, raises

from askcell.utils.fileutils import load_toml
from askcell.utils.interpolate import interpolate, interpolate_expr


@mark.parametrize(
    "expr, expected, exception, message",
    [
        ("a", None, raises(ValueError), "Maximum recursion limit exceeded: 20"),
        ("b", None, raises(ValueError), "Maximum recursion limit exceeded: 20"),
        ("", None, raises(ValueError), "no value to resolve"),
        ("d.g", None, raises(KeyError), "g"),
        ("g.h", None, raises(ValueError), "Interpolation into non-dictionary: g"),
        ("d.e", [1, 2], nullcontext(), None),
        ("f", {"e": [1, 2]}, nullcontext(), None),
    ],
)
def test_interpolate_expr(expr: str, expected: Any, exception: Any, message: str | None) -> None:
    """Test looking up expressions."""
    data = {
        "a": "{a}",
        "b": "{c}",
        "c": "{b}",
        "d": {"e": [1, 2]},
        "f": "{.d}",
        "g": 0,
    }
    with exception as excinfo:
        assert interpolate_expr(data, expr) == expected

    if message:
        assert message in str(excinfo)


def test_interpolate_toml() -> None:
    """Test interpolation of a TOML file."""
    config = load_toml("tests/data/test_utils/test_interpolate/static_data.toml")

    config["human"]["base"] = "/tmp_human"
    config["mouse"]["base"] = "/tmp_mouse"

    expected = load_toml("tests/data/test_utils/test_interpolate/static_data_interpolated.toml")

    assert interpolate(config) == expected
