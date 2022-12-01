"""Test pathutils.py."""

from __future__ import annotations

from typing import Sequence

from pytest import mark, raises

from askcell import Path, PathLike, SomePathLikes
from askcell.utils.pathutils import (
    match_path,
    parse_path,
    parse_paths,
    set_path_readonly,
    validate_path,
)


@mark.parametrize(
    "path, expected",
    [
        ("~root/temp", Path("/root/temp")),
        ("askcell", Path.cwd() / "askcell"),
        (Path("/home/user"), Path("/home/user")),
    ],
)
def test_parse_path(path: PathLike, expected: Path) -> None:
    """Test parse_path."""
    assert expected == parse_path(path, resolve=True)


@mark.parametrize(
    "paths, expected",
    [
        (".", [Path(".")]),
        (Path("."), [Path(".")]),
        ((".", Path("/tmp")), [Path("."), Path("/tmp")]),
    ],
)
def test_parse_paths(paths: SomePathLikes, expected: Sequence[Path]) -> None:
    """Test parse_paths."""
    assert list(parse_paths(paths)) == expected


@mark.parametrize(
    "permission_octal, expected_octal",
    [
        (0o100466, 0o100444),
        (0o100664, 0o100444),
        (0o100646, 0o100444),
        (0o100666, 0o100444),
        (0o100744, 0o100544),
        (0o100674, 0o100454),
        (0o100647, 0o100445),
        (0o100777, 0o100555),
    ],
)
def test_set_path_readonly(tmp_dirname: Path, permission_octal: int, expected_octal: int) -> None:
    """Test set path readonly permission."""
    p = Path(tmp_dirname / "temp.txt")
    p.touch()
    p.chmod(permission_octal)
    assert p.stat().st_mode == permission_octal

    set_path_readonly(p)
    assert p.stat().st_mode == expected_octal


@mark.parametrize(
    "path, pattern, expected",
    [
        ("", "a", False),
        ("", "", True),
        (".", "a", False),
        (".", "", True),
        ("a.b", "*.b", True),
        ("a.b", "a.b", True),
        ("a.b", ".b", False),
        ("a.b", "", False),
        ("/tmp/a.b", "*.b", True),
        ("/tmp/a.b", "a.b", True),
        ("/tmp/a.b", ".b", False),
        ("/tmp/a.b", "/tmp/a.b", False),
        (Path(""), "a", False),
        (Path(""), "", True),
        (Path("."), "a", False),
        (Path("."), "", True),
        (Path("a.b"), "*.b", True),
        (Path("a.b"), "a.b", True),
        (Path("a.b"), ".b", False),
        (Path("a.b"), "", False),
        (Path("/tmp/a.b"), "*.b", True),
        (Path("/tmp/a.b"), "a.b", True),
        (Path("/tmp/a.b"), ".b", False),
        (Path("/tmp/a.b"), "/tmp/a.b", False),
    ],
)
def test_match_path(path: PathLike, pattern: str, expected: bool) -> None:
    """Test match_path."""
    assert match_path(path, pattern) == expected


@mark.parametrize(
    "path, check_exists, allow_posix, allow_s3, allow_gcs, allow_windows, msg",
    [
        ("/dir_does_not_exist", True, True, False, False, False, "Path does not exist: /dir_does_not_exist"),
        ("s3://no/nonono", True, True, False, False, False, "S3 path is not allowed in this context: s3://no/nonono"),
        ("s3://no/nonono", False, True, False, False, False, "S3 path is not allowed in this context: s3://no/nonono"),
        ("gs://no/nonono", True, True, False, False, False, "GCS path is not allowed in this context: gs://no/nonono"),
        ("gs://no/nonono", False, True, False, False, False, "GCS path is not allowed in this context: gs://no/nonono"),
    ],
)
def test_validate_path_error(
    path: PathLike,
    check_exists: bool,
    allow_posix: bool,
    allow_s3: bool,
    allow_gcs: bool,
    allow_windows: bool,
    msg: str,
) -> None:
    """Test that an error is raised from validate_path."""
    with raises(ValueError) as excinfo:
        validate_path(
            parse_path(path),
            check_exists=check_exists,
            allow_posix=allow_posix,
            allow_s3=allow_s3,
            allow_gcs=allow_gcs,
            allow_windows=allow_windows,
        )
    assert msg == str(excinfo.value)
