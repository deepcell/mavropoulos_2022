"""Test askcell.utils.config."""
from __future__ import annotations

import argparse
import itertools
import sys

from datetime import date, datetime
from importlib import resources
from typing import Any, Sequence

from pytest import mark

from askcell import only_one
from askcell.utils.config import (
    PackageResource,
    SomeConfigItems,
    add_config_args,
    collect_config_files,
    load_config_data,
    load_package_data,
    load_package_schema,
    merge_config_data,
    set_expr,
)
from askcell.utils.fileutils import load_toml
from askcell.utils.merge import merge
from askcell.utils.pathutils import Path
from askcell.utils.types import MutableStrMapping, parse_schema


CONFIG_TEST_DIR = Path("tests/data/test_utils/test_config")
if CONFIG_TEST_DIR / "a-package.zip" not in sys.path:
    sys.path.append(str(CONFIG_TEST_DIR / "a-package.zip"))


def test_add_config_args() -> None:
    """Test add_config_args."""
    parser = argparse.ArgumentParser(prog="test_add_config_args")
    add_config_args(parser)
    argv = [
        "-c",
        str(CONFIG_TEST_DIR),
        "-s",
        "key=value",
        "-c",
        "/tmp/testing.toml",
        "-s",
        "key0.key1=value",
    ]
    args = parser.parse_args(argv)

    assert args.config == [CONFIG_TEST_DIR, Path("/tmp/testing.toml")]
    assert args.set_config == ["key=value", "key0.key1=value"]


def test_load_config_data() -> None:
    """Test load_config_data."""
    expected = {
        "gc": {
            "id": "",
            "region": "us-central-1",
            "gcs": "gs://us.artifacts.dc-bioinformatics.appspot.com/",
            "profile": "default",
        },
    }
    assert merge(load_config_data(CONFIG_TEST_DIR / "test_merge_config.toml")) == expected

    expected.update(
        {
            "ENV": {
                "TEST_DIR1": "/foo/bar",
                "TEST_DIR2": "/temp/test",
            },
        }
    )
    assert merge(load_config_data(collect_config_files([CONFIG_TEST_DIR], pattern="*config.toml"))) == expected


def test_collect_config_files() -> None:
    """Test collect_config_files."""
    expected = [
        CONFIG_TEST_DIR / "test_merge_config.toml",
        CONFIG_TEST_DIR / "test_update_environment_config.toml",
    ]
    assert list(collect_config_files(CONFIG_TEST_DIR, pattern="*config.toml")) == expected


@mark.parametrize(
    "item",
    [
        [],
        {},
    ],
)
def test_merge_config_data_empty(item: Any) -> None:
    """Test merge_config_data returning an empty dictionary."""
    assert merge_config_data(item) == {}


@mark.parametrize(
    "package, pattern, names",
    [
        ("a.b", "*config.toml", ["test_merge_config.toml", "test_update_environment_config.toml"]),
        ("a.b", "test_merge_config.toml", ["test_merge_config.toml"]),
    ],
)
def test_package_resource(package: resources.Package, pattern: str, names: Sequence[str]) -> None:
    """Test PackageResource."""
    expected = []
    for name in names:
        expected.append(only_one(load_config_data(CONFIG_TEST_DIR / name)))

    items = PackageResource(package, pattern)
    for i, path in enumerate(items.get_resource()):
        assert only_one(load_config_data(path)) == expected[i]


@mark.parametrize(
    "package, pattern",
    [
        ("a", "b"),
    ],
)
def test_package_resource_empty(package: str, pattern: str) -> None:
    """Test PackageResource returning no resources."""
    assert list(PackageResource(package, pattern).get_resource()) == []


def test_merge_config_data() -> None:
    """Test merge_config_data using a mixture of items."""
    fragment = {
        "foo": "bar",
    }
    package_resource = PackageResource("a.b", "test_merge_config.toml")
    path = CONFIG_TEST_DIR / "test_merge_config.toml"
    expected = {
        "gc": {
            "id": "",
            "region": "us-central-1",
            "gcs": "gs://us.artifacts.dc-bioinformatics.appspot.com/",
            "profile": "default",
        },
        "foo": "bar",
    }

    for permutation in itertools.permutations([fragment, path, package_resource]):
        assert merge_config_data(permutation) == expected  # type: ignore


def test_load_package_schema_empty() -> None:
    """Test load_package_schema using an empty list."""
    assert load_package_schema([]) == {}


def test_load_package_schema1() -> None:
    """Test load_package_schema using a package resource."""
    expected = parse_schema(load_toml(CONFIG_TEST_DIR / "test_merge_config_schema.toml"))
    schema = load_package_schema(PackageResource("a.b", "test_merge_config_schema.toml"))
    assert schema == expected
    assert isinstance(schema, dict)


def test_load_package_schema2() -> None:
    """Test load_package_schema using package resources and mappings."""
    fragment = {
        "gc": {
            "id": {"type": "str"},
            "region": {"type": "str", "default": "us-central-1"},
            "gcs": {"type": "str", "default": "gs://us.artifacts.dc-bioinformatics.appspot.com/"},
            "profile": {"type": "str", "default": "default"},
        },
    }
    expected = parse_schema(merge([load_toml(CONFIG_TEST_DIR / "test_merge_config_schema.toml"), fragment]))
    items: SomeConfigItems = [
        PackageResource("a.b", "test_merge_config.toml"),
        fragment,
    ]
    assert load_package_schema(items) == expected


SCHEMA = parse_schema(
    {
        "outer": {
            "bool": {"type": "bool"},
            "date": {"type": "date"},
            "datetime": {"type": "datetime"},
            "float": {"type": "float"},
            "int": {"type": "int"},
            "str": {"type": "str"},
            "path": {"type": "path"},
            "inner": {
                "int": {"type": "int"},
            },
        },
    }
)


def test_set_expr1() -> None:
    """Test set_expr."""
    data = {
        "outer": {
            "bool": True,
            "date": date(2020, 2, 10),
            "datetime": datetime(2020, 2, 10),
        },
    }
    set_expr(SCHEMA, data, "outer.bool=False")
    set_expr(SCHEMA, data, "outer.date=2020-02-01")
    set_expr(SCHEMA, data, "outer.datetime=2011-11-04 00:05:23.283")
    set_expr(SCHEMA, data, "outer.float=-1")
    set_expr(SCHEMA, data, "outer.int=-1")
    set_expr(SCHEMA, data, "outer.str=bar")
    assert data == {
        "outer": {
            "bool": False,
            "date": date(2020, 2, 1),
            "datetime": datetime(2011, 11, 4, 0, 5, 23, 283000),
            "float": -1.0,
            "int": -1,
            "str": "bar",
        },
    }


def test_set_expr2() -> None:
    """Test set_expr where values are empty strings."""
    data = {
        "outer": {
            "float": 0.0,
            "int": 0,
            "str": "foo",
        },
    }
    set_expr(SCHEMA, data, "outer.bool=")
    set_expr(SCHEMA, data, "outer.date=")
    set_expr(SCHEMA, data, "outer.datetime=")
    set_expr(SCHEMA, data, "outer.float=")
    set_expr(SCHEMA, data, "outer.int=")
    set_expr(SCHEMA, data, "outer.str=")
    set_expr(SCHEMA, data, "outer.path=")
    assert data == {
        "outer": {
            "bool": None,
            "date": None,
            "datetime": None,
            "float": None,
            "int": None,
            "str": None,
            "path": None,
        },
    }


def test_set_expr3() -> None:
    """Test set_expr where bad values are passed."""
    data: dict[str, Any] = {
        "outer": {},
    }
    set_expr(SCHEMA, data, "outer.bool=Foo")
    set_expr(SCHEMA, data, "outer.date=Foo")
    set_expr(SCHEMA, data, "outer.datetime=Foo")
    set_expr(SCHEMA, data, "outer.float=Foo")
    set_expr(SCHEMA, data, "outer.int=Foo")
    assert data == {
        "outer": {
            "bool": "Foo",
            "date": "Foo",
            "datetime": "Foo",
            "float": "Foo",
            "int": "Foo",
        },
    }


def test_set_expr4() -> None:
    """Test set_expr fills in missing levels."""
    data: MutableStrMapping = {}
    set_expr(SCHEMA, data, "outer.inner.int=0")
    assert data == {
        "outer": {
            "inner": {
                "int": 0,
            },
        },
    }


def test_load_package_data1() -> None:
    """Test load_package_data."""
    schema = load_package_schema(PackageResource("a.b", "test_merge_config_schema.toml"))
    assert load_package_data(schema) == {
        "gc": {
            "id": None,
            "region": "us-central-1",
            "gcs": "gs://us.artifacts.dc-bioinformatics.appspot.com/",
            "profile": "default",
        },
    }


def test_load_package_data2() -> None:
    """Test load_package_data with settings."""
    schema = load_package_schema(
        {
            "gc": {
                "id": {"type": "str"},
                "region": {"type": "str", "default": "us-central-1"},
                "gcs": {"type": "str", "default": "gs://us.artifacts.dc-bioinformatics.appspot.com/"},
                "profile": {"type": "str", "default": "default"},
            },
            "network": {"port": {"type": "int", "default": "8080"}, "address": {"type": "str", "default": "127.0.0.1"}},
        }
    )
    settings = ["gc.id=test"]
    data = {
        "network": {
            "port": None,
            "address": "192.168.0.1",
        },
        "not_in_schema": "something",
    }
    assert load_package_data(schema, data, settings=settings) == {
        "gc": {
            "id": "test",
            "region": "us-central-1",
            "gcs": "gs://us.artifacts.dc-bioinformatics.appspot.com/",
            "profile": "default",
        },
        "network": {
            "port": None,
            "address": "192.168.0.1",
        },
    }
