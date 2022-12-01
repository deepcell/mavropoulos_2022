"""Set up temp directory factory function."""

import shutil

from pathlib import Path
from typing import Iterator

from pytest import TempdirFactory, fixture


@fixture(scope="session")
def tmp_dirname(tmpdir_factory: TempdirFactory) -> Iterator[Path]:
    """Temp directory factory function."""
    p = Path(tmpdir_factory.mktemp("tmp_test_askcell"))
    yield p
    shutil.rmtree(p)
