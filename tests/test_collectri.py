"""CollecTRI as the signed table behind --tf-activity.

The activity test needs edges that say whether a factor raises or lowers its
targets. CollecTRI signs every one; ReMap occupancy signs none. These tests pin
that the right table reaches the test, that an explicit choice is respected,
and that a CollecTRI file is recognised for what it is.
"""

from __future__ import annotations

import logging

import pytest

from generoutes import pipeline, regulon, tfactivity

TAB = chr(9)
NL = chr(10)


def _write(path, header, rows):
    path.write_text(NL.join([TAB.join(header)] + [TAB.join(map(str, r))
                                                  for r in rows]) + NL,
                    encoding="utf-8")
    return path


def _collectri(path, n_factors=3, n_targets=6):
    rows = []
    for f in range(n_factors):
        for t in range(n_targets):
            rows.append([f"TF{f}", f"G{t}", 1.0 if t % 3 else -1.0,
                         "ExTRI;TRRUST", "12345", "PMID"])
    return _write(path, ["source", "target", "weight", "resources",
                         "references", "sign_decision"], rows)


def _chip(path):
    rows = [[f"CF{f}", f"G{t}", "x", "ReMap:cells"]
            for f in range(3) for t in range(6)]
    return _write(path, ["tf", "target", "mode", "evidence"], rows)


def _trrust(path):
    rows = [[f"TR{f}", f"G{t}", "Activation" if t % 2 else "Unknown", "999"]
            for f in range(3) for t in range(6)]
    path.write_text(NL.join(TAB.join(map(str, r)) for r in rows) + NL,
                    encoding="utf-8")
    return path


# ---------------------------------------------------------------- the file
def test_a_collectri_file_is_recognised_as_collectri(tmp_path):
    kind, phrase = regulon.describe_source(_collectri(tmp_path / "c.tsv"))
    assert kind == "collectri" and "CollecTRI" in phrase
    assert "signed" in phrase


def test_the_other_tables_are_still_told_apart(tmp_path):
    assert regulon.describe_source(_chip(tmp_path / "chip.tsv"))[0] == "chip"
    assert regulon.describe_source(_trrust(tmp_path / "t.tsv"))[0] == "curated"


def test_every_collectri_edge_is_signed_from_its_weight(tmp_path):
    got = tfactivity.load_signed_regulons(_collectri(tmp_path / "c.tsv"))
    assert set(got) == {"TF0", "TF1", "TF2"}
    assert got["TF0"]["G0"] == -1 and got["TF0"]["G1"] == 1
    assert all(s in (1, -1) for tg in got.values() for s in tg.values()), \
        "no edge is left unsigned"


def test_trrust_leaves_edges_unsigned_and_chip_signs_none(tmp_path):
    tr = tfactivity.load_signed_regulons(_trrust(tmp_path / "t.tsv"))
    assert any(s == 0 for tg in tr.values() for s in tg.values())
    ch = tfactivity.load_signed_regulons(_chip(tmp_path / "chip.tsv"))
    assert all(s == 0 for tg in ch.values() for s in tg.values())


def test_a_comma_separated_collectri_loads_too(tmp_path):
    p = _collectri(tmp_path / "c.csv")
    p.write_text(p.read_text().replace(TAB, ","), encoding="utf-8")
    assert len(tfactivity.load_signed_regulons(p)) == 3


def test_autoload_prefers_collectri_over_trrust(tmp_path, monkeypatch):
    (tmp_path / "data").mkdir()
    _collectri(tmp_path / "data" / "collectri.tsv")
    _trrust(tmp_path / "data" / "trrust.tsv")
    monkeypatch.chdir(tmp_path)
    regs, path = tfactivity.autoload_regulons()
    assert path.endswith("collectri.tsv") and "TF0" in regs
    (tmp_path / "data" / "collectri.tsv").unlink()
    assert tfactivity.autoload_regulons()[1].endswith("trrust.tsv")


# ------------------------------------------------- which table the test gets
@pytest.fixture
def data(tmp_path, monkeypatch):
    (tmp_path / "data").mkdir()
    monkeypatch.chdir(tmp_path)
    return tmp_path / "data"


def test_an_occupancy_table_gives_way_to_a_signed_one(data, caplog):
    _collectri(data / "collectri.tsv")
    chip = _chip(data / "chip.tsv")
    with caplog.at_level(logging.INFO):
        regs, source = pipeline._activity_regulons(str(chip))
    assert source.endswith("collectri.tsv") and "TF0" in regs
    assert "unsigned occupancy" in caplog.text
    assert "still used for the regulon breadth count" in caplog.text


def test_an_occupancy_table_alone_is_used_with_a_warning(data, caplog):
    chip = _chip(data / "chip.tsv")
    with caplog.at_level(logging.WARNING):
        regs, source = pipeline._activity_regulons(str(chip))
    assert source == str(chip) and "CF0" in regs
    assert "which way" in caplog.text and "CollecTRI" in caplog.text


def test_a_curated_or_collectri_file_given_on_purpose_wins(data):
    _collectri(data / "collectri.tsv")
    trrust = _trrust(data / "trrust.tsv")
    regs, source = pipeline._activity_regulons(str(trrust))
    assert source == str(trrust) and "TR0" in regs and "TF0" not in regs
    other = _collectri(data / "mine.tsv", n_factors=2)
    regs2, source2 = pipeline._activity_regulons(str(other))
    assert source2 == str(other) and set(regs2) == {"TF0", "TF1"}


def test_with_no_file_given_the_one_in_data_is_used(data):
    _collectri(data / "collectri.tsv")
    regs, source = pipeline._activity_regulons(None)
    assert source.endswith("collectri.tsv") and len(regs) == 3


def test_with_nothing_anywhere_it_says_where_to_get_one(data, caplog):
    with caplog.at_level(logging.WARNING):
        regs, source = pipeline._activity_regulons(None)
    assert not regs and source is None
    assert "found no regulon table" in caplog.text
    assert "decoupler" in caplog.text


def test_an_unreadable_given_file_does_not_stop_the_run(data):
    bad = data / "bad.tsv"
    bad.write_bytes(b"\xff\xfe\x00")
    regs, source = pipeline._activity_regulons(str(bad))
    assert source == str(bad)        # parsed or failed softly, never raised
