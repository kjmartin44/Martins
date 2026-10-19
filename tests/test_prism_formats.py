"""The PRISM loader reads the layouts the releases actually ship.

The first version assumed cell lines in rows. The 24Q2 matrix has a compound
per row and a cell line per column, and the file beside it is a long table
with one row per well. With the file the package named as its default, the
loader returned 1.5 million "cell lines" and every phenocopy search failed
for want of overlap with DepMap.
"""

from __future__ import annotations

import logging

import numpy as np

from promoter_tf_drugs import sensitivity

NL = chr(10)


def _write(tmp_path, name, rows):
    p = tmp_path / name
    p.write_text(NL.join(",".join(r) for r in rows) + NL, encoding="utf-8")
    return p


def test_a_matrix_with_cell_lines_in_rows_is_read_as_is(tmp_path):
    p = _write(tmp_path, "old.csv", [
        ["", "BRD:BRD-K1", "BRD:BRD-K2"],
        ["ACH-000001", "-1.5", ""],
        ["ACH-000002", "0.25", "-0.5"]])
    got = sensitivity.load_prism(str(p))
    assert got["models"] == ["ACH-000001", "ACH-000002"]
    assert got["columns"] == ["BRD:BRD-K1", "BRD:BRD-K2"]
    assert got["values"].shape == (2, 2)
    assert got["values"][0, 0] == -1.5 and np.isnan(got["values"][0, 1])


def test_a_matrix_with_cell_lines_in_columns_is_turned_round(tmp_path):
    """The 24Q2 layout: a row per compound, a column per cell line. The
    result must equal the other layout's, so every caller sees one shape."""
    p = _write(tmp_path, "new.csv", [
        ["", "ACH-000001", "ACH-000002"],
        ["BRD:BRD-K1", "-1.5", "0.25"],
        ["BRD:BRD-K2", "", "-0.5"]])
    got = sensitivity.load_prism(str(p))
    assert got["models"] == ["ACH-000001", "ACH-000002"]
    assert got["columns"] == ["BRD:BRD-K1", "BRD:BRD-K2"]
    assert got["values"].shape == (2, 2)
    assert got["values"][0, 0] == -1.5 and got["values"][1, 0] == 0.25
    assert got["values"][1, 1] == -0.5 and np.isnan(got["values"][0, 1])


def test_both_layouts_give_the_same_matrix(tmp_path):
    rows = [["ACH-000001", "-1.5", ""], ["ACH-000002", "0.25", "-0.5"]]
    a = _write(tmp_path, "a.csv", [["", "BRD:K1", "BRD:K2"]] + rows)
    b = _write(tmp_path, "b.csv", [
        ["", "ACH-000001", "ACH-000002"],
        ["BRD:K1", "-1.5", "0.25"], ["BRD:K2", "", "-0.5"]])
    x, y = sensitivity.load_prism(str(a)), sensitivity.load_prism(str(b))
    assert x["models"] == y["models"] and x["columns"] == y["columns"]
    assert np.array_equal(x["values"], y["values"], equal_nan=True)


def test_the_long_one_row_per_well_table_is_refused_and_named(tmp_path, caplog):
    p = _write(tmp_path, "long.csv", [
        ["row_id", "broad_id", "dose", "compound_plate", "LFC"],
        ["ACH-000001::P1", "BRD-A1", "2.5", "PREP053", "-0.3"]])
    with caplog.at_level(logging.WARNING):
        assert sensitivity.load_prism(str(p)) == {}
    assert any("long format" in r.message and "Data_Matrix" in r.message
               for r in caplog.records)


def test_the_default_search_order_puts_the_matrix_first():
    first = sensitivity.DEFAULT_PRISM_PATHS[0]
    assert first.endswith("Extended_Primary_Data_Matrix.csv")
    assert "data/Repurposing_Public_24Q2_LFC_COLLAPSED.csv" in (
        sensitivity.DEFAULT_PRISM_PATHS)


def test_the_24q2_compound_list_columns_give_names(tmp_path):
    p = _write(tmp_path, "compounds.csv", [
        ["screen", "dose", "repurposing_target", "MOA", "IDs", "Drug.Name",
         "Synonyms"],
        ["REP.300", "2.5", "", "", "BRD:BRD-K00104122-001-01-9", "MRTX1133",
         "MRTX1133"],
        ["REP.300", "2.5", "", "", "BRD:BRD-K00104124-001-01-9", "PEMIGATINIB",
         ""]])
    names = sensitivity.load_compound_names(str(p))
    assert names == {"BRD:BRD-K00104122-001-01-9": "MRTX1133",
                     "BRD:BRD-K00104124-001-01-9": "PEMIGATINIB"}
