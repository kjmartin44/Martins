"""Shared test setup."""

import pytest

from promoter_tf_drugs import fda


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
    from promoter_tf_drugs import exposure
    monkeypatch.setattr(exposure, "_default_tables", lambda: {})
