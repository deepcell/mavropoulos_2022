"""Test askcell.utils.types module."""

from __future__ import annotations

from typing import Any

from pytest import mark, raises

from askcell.utils.fileutils import load_toml
from askcell.utils.pathutils import Path
from askcell.utils.types import (
    BoolSchemaType,
    DateSchemaType,
    DateTimeSchemaType,
    FloatSchemaType,
    IntSchemaType,
    ListSchemaType,
    MutableStrMapping,
    PathSchemaType,
    SchemaType,
    StrSchemaType,
    apply_schema,
    build_data_from_schema,
    parse_schema,
    parse_schema_type,
    unparse_data,
)


@mark.parametrize(
    "dtype, data, expected",
    [
        ("str", {"required": True, "default": "Foo"}, StrSchemaType),
        ("int", {"default": 0}, IntSchemaType),
        ("float", {"required": True}, FloatSchemaType),
        ("bool", {"required": False}, BoolSchemaType),
        ("date", {}, DateSchemaType),
        ("datetime", {}, DateTimeSchemaType),
        ("path", {}, PathSchemaType),
    ],
)
def test_parse_schema_type1(dtype: str, data: MutableStrMapping, expected: type[SchemaType[Any, Any]]) -> None:
    """Test parse_schema_type."""
    stype = expected(**data)
    assert parse_schema_type(type=dtype, **data) == stype
    assert data.get("required", False) == stype.required
    assert data.get("default") == stype.default


@mark.parametrize(
    "dtype, etype",
    [
        ("list[str]", StrSchemaType),
        ("list[int]", IntSchemaType),
        ("list[float]", FloatSchemaType),
        ("list[bool]", BoolSchemaType),
        ("list[date]", DateSchemaType),
        ("list[datetime]", DateTimeSchemaType),
        ("list[path]", PathSchemaType),
    ],
)
def test_parse_schema_type2(dtype: str, etype: Any) -> None:
    """Test parse_schema_type using lists."""
    stype = parse_schema_type(type=dtype)
    assert isinstance(stype, ListSchemaType)
    assert isinstance(stype.element_type, etype)


@mark.parametrize(
    "data",
    [
        {"type": "list[str]", "length": 2, "min_length": 1},
        {"type": "list[str]", "length": 2, "max_length": 3},
        {"type": "list[str]", "length": 2, "min_length": 1, "max_length": 2},
        {"type": "list[str]", "length": 2, "min_length": 2, "max_length": 3},
        {"type": "list[str]", "length": 2, "min_length": 1, "max_length": 3},
    ],
)
def test_parse_schema_type_errors(data: MutableStrMapping) -> None:
    """Test parse_schema_type raising an error for lists if given contradictory lengths."""
    with raises(ValueError) as excinfo:
        parse_schema_type(**data)
    assert "schema value definition has contradictory lengths" in str(excinfo)


def test_parse_schema() -> None:
    """Test parse_schema using schema.toml."""
    raw_config = load_toml(Path("tests/data/test_utils/test_types/schema.toml"))
    schema = parse_schema(raw_config)
    assert schema == {
        "type_string": {
            "string": StrSchemaType(default="a"),
        },
        "type_list": {
            "list_int": ListSchemaType(type="int", default=[1, 2, 3, 4, 5]),
            "list_float": ListSchemaType(type="float", default=[0.1, 0.2, 0.3, 0.4, 0.5]),
            "list_string": ListSchemaType(type="str", default=["a", "b", "c"]),
        },
    }


def test_build_data_from_schema() -> None:
    """Test build_data_from_schema using schema.toml."""
    raw_config = load_toml(Path("tests/data/test_utils/test_types/schema.toml"))
    assert build_data_from_schema(parse_schema(raw_config)) == {
        "type_string": {
            "string": "a",
        },
        "type_list": {
            "list_int": [1, 2, 3, 4, 5],
            "list_float": [0.1, 0.2, 0.3, 0.4, 0.5],
            "list_string": ["a", "b", "c"],
        },
    }


def test_apply_schema() -> None:
    """Test apply_schema using schema.toml."""
    raw_config = load_toml(Path("tests/data/test_utils/test_types/schema.toml"))
    data = {
        "type_string": {
            "attr0": "a",
            "attr1": "b",
        },
        "type_list": {
            "add_list_str": ["a", "b", "c", "d"],
        },
    }
    assert apply_schema(parse_schema(raw_config), data) == {
        "type_string": {
            "string": "a",
        },
        "type_list": {
            "list_int": [1, 2, 3, 4, 5],
            "list_float": [0.1, 0.2, 0.3, 0.4, 0.5],
            "list_string": ["a", "b", "c"],
        },
    }


def test_unparse_data() -> None:
    """Test unparse_data."""
    raw_config = load_toml(Path("tests/data/test_utils/test_types/unparse_data.toml"))
    schema = parse_schema(raw_config)
    schema_slice = {
        "a": {
            "b": {
                "params": schema["a"]["b"]["params"],
            },
        },
        "e": {
            "params": schema["e"]["params"],
        },
    }
    data = build_data_from_schema(schema_slice)
    assert unparse_data(schema, data) == {
        "a": {
            "b": {
                "params": {
                    "a": 1,
                    "b": "-p {sample.id} --s {.a}",
                    "c": 20,
                    "d": "-p {.c}",
                },
            },
        },
        "e": {
            "params": {
                "f": "{e.input.x}",
                "g": "{e.scratch.y}",
                "h": "{e.input.z}",
                "i": "{e.scratch.ss}",
            },
        },
    }
