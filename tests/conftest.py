"""Shared test setup."""

import os
import shutil
import tempfile
import time
from pathlib import Path

import pytest

from generoutes import fda

# How long an abandoned run directory is kept before a later run clears
# it. Long enough to survive the slowest full suite several times over,
# short enough that a crashed run does not leave gigabytes behind.
_RUN_DIR_MAX_AGE = 12 * 3600


def pytest_configure(config):
    """Put this run's temporary files somewhere no other run will touch.

    Left to itself pytest puts every `tmp_path` under one shared root,
    `<temp>/pytest-of-<user>/pytest-N`, keeps the last three of those and
    **prunes the older ones at startup**. One checkout with several
    sessions running the suite at once is all it takes: a run that starts
    while another is working can delete the directory the first run's
    fixtures are living in, halfway through.

    What that looks like is not an error about missing files. It is a
    test failing on its own subject matter - the memoisation tests fail
    first, because they are the ones that need a file written on one line
    to still be there on the next - so the failure reads as a bug in
    whatever those tests cover, and bisecting it finds nothing.

    The run's own process id is enough to separate them: two live
    processes never share one, and a directory left by a dead process is
    safe for a later run to clear.

    An explicit `--basetemp` always wins, so a caller who wants their
    files somewhere specific still gets them there.
    """
    if config.option.basetemp:
        return
    root = Path(tempfile.gettempdir()) / "generoutes-pytest"
    # pytest makes the basetemp itself and not its parents, so the shared
    # root has to exist before it is handed one underneath it.
    root.mkdir(parents=True, exist_ok=True)
    _clear_abandoned_runs(root)
    config.option.basetemp = str(root / f"run-{os.getpid()}")


def _clear_abandoned_runs(root: Path) -> None:
    """Remove run directories old enough that nothing can still be using one.

    Never fatal: a directory that cannot be read or removed belongs to
    somebody else or to a locked file, and leaving it costs only disk.
    """
    cutoff = time.time() - _RUN_DIR_MAX_AGE
    try:
        stale = [d for d in root.glob("run-*")
                 if d.is_dir() and d.stat().st_mtime < cutoff]
    except OSError:
        return
    for d in stale:
        shutil.rmtree(d, ignore_errors=True)


@pytest.fixture(autouse=True)
def _no_fda_filter_unless_a_test_asks(monkeypatch):
    """Tests build results from invented drug names, and a machine with an
    Orange Book in data/ would hide every one of them. The filter is off here,
    as it is on a machine with no Orange Book, and the tests of the filter
    switch it on with their own table."""
    monkeypatch.setattr(fda, "available", lambda: False)
    fda.reset()


@pytest.fixture(autouse=True)
def _no_local_cmax_table(monkeypatch):
    """The dose check reads data/fda_label_cmax.tsv when it exists. A test
    must not depend on whether the person running it has built that file."""
    from generoutes import exposure
    monkeypatch.setattr(exposure, "_default_tables", lambda: {})
