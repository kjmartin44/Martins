"""Resolving a ChEMBL id for a drug that arrived with an RxNorm or NCIt id.

A drug given up on used to be scored as average and have its safety lookups
skipped without a note. These tests pin that a resolved id is real, an
unresolved one is reported as unresolved, and an id ChEMBL supplied is never
replaced.
"""

from __future__ import annotations

import logging

import pytest

from generoutes import chemblid as cid
from generoutes.http import HttpError

ROWS = [{"name": "RIFAMPIN", "chembl_id": "CHEMBL374478"},
        {"name": "IMATINIB MESYLATE", "chembl_id": "CHEMBL1642"},
        {"name": "TRETINOIN", "chembl_id": "CHEMBL38"},
        {"name": "ODDITY", "chembl_id": "rxcui:77"},          # not a ChEMBL id
        {"name": "RIFAMPIN", "chembl_id": "CHEMBLDUP"}]       # a duplicate


class _Client:
    def __init__(self, parents=None, search=None, fail=()):
        self.parents, self.search_results = parents or {}, search or []
        self.fail, self.calls = set(fail), []

    def get(self, url, params=None, **kw):
        self.calls.append((url, params))
        for marker in self.fail:
            if marker in url:
                raise HttpError("503")
        if "/molecule/search" in url:
            return {"molecules": self.search_results}
        cid_ = url.rsplit("/", 1)[-1].replace(".json", "")
        parent = self.parents.get(cid_)
        return {"molecule_hierarchy": {"parent_chembl_id": parent}
                if parent else None}


@pytest.fixture(autouse=True)
def _table():
    cid.library_index(ROWS)
    yield
    cid.reset()


def test_an_id_chembl_supplied_is_returned_untouched_with_no_request():
    c = _Client()
    assert cid.resolve(c, "Anything", "CHEMBL25") == ("CHEMBL25", "given")
    assert c.calls == []


def test_an_rxnorm_drug_is_found_by_name_in_the_approved_table():
    c = _Client()
    got = cid.resolve(c, "Rifampin", "rxcui:9384")
    assert got == ("CHEMBL374478", "approved table")


def test_case_hyphens_and_salt_words_do_not_stop_a_match():
    c = _Client()
    assert cid.resolve(c, "imatinib", "ncit:C1")[0] == "CHEMBL1642"
    assert cid.resolve(c, "Imatinib Mesylate", None)[0] == "CHEMBL1642"
    assert cid.resolve(c, "TRETINOIN", None)[0] == "CHEMBL38"


def test_the_first_table_row_for_a_name_wins_and_non_chembl_ids_are_ignored():
    idx = cid.library_index(ROWS)
    assert idx["rifampin"] == "CHEMBL374478"
    assert "oddity" not in idx


def test_the_id_moves_to_the_parent_molecule():
    c = _Client(parents={"CHEMBL1642": "CHEMBL941"})
    assert cid.resolve(c, "imatinib", "rxcui:1") == (
        "CHEMBL941", "approved table")


def test_a_failed_parent_lookup_keeps_the_found_id_and_says_so(caplog):
    c = _Client(fail=("/molecule/CHEMBL38",))
    with caplog.at_level(logging.WARNING):
        got = cid.resolve(c, "tretinoin", "rxcui:2")
    assert got == ("CHEMBL38", "approved table")
    assert "could not give the parent of CHEMBL38" in caplog.text


def test_a_name_missing_from_the_table_falls_back_to_an_exact_search_hit():
    c = _Client(search=[
        {"pref_name": "NEWDRUG ACID", "molecule_chembl_id": "CHEMBL9"},
        {"pref_name": "NEWDRUG", "molecule_chembl_id": "CHEMBL7"}],
        parents={"CHEMBL7": "CHEMBL70"})
    assert cid.resolve(c, "Newdrug", "rxcui:3") == ("CHEMBL70", "name search")


def test_a_near_match_is_not_the_drug():
    c = _Client(search=[{"pref_name": "RIFAMPICIN SV",
                         "molecule_chembl_id": "CHEMBL5"}])
    assert cid.resolve(c, "Rifapentine", "rxcui:4") == (None, "unresolved")


def test_a_failed_search_is_unresolved_not_an_exception(caplog):
    c = _Client(fail=("/molecule/search",))
    with caplog.at_level(logging.WARNING):
        assert cid.resolve(c, "Newdrug", None) == (None, "unresolved")
    assert "name search failed" in caplog.text


def test_an_empty_name_is_unresolved_without_a_request():
    c = _Client()
    assert cid.resolve(c, "  ", "rxcui:5") == (None, "unresolved")
    assert cid.resolve(c, None, None) == (None, "unresolved")
    assert c.calls == []


def test_with_no_approved_table_only_the_search_can_resolve(monkeypatch):
    from generoutes import localdrugs
    cid.reset()
    monkeypatch.setattr(localdrugs, "library", lambda path=None: [])
    c = _Client(search=[{"pref_name": "RIFAMPIN",
                         "molecule_chembl_id": "CHEMBL374478"}])
    assert cid.resolve(c, "Rifampin", None) == ("CHEMBL374478", "name search")
    assert cid.library_index() == {}


# ----------------------------------------------------------- target ids
class _TargetClient:
    def __init__(self, targets=None, fail=False):
        self.targets, self.fail, self.calls = targets or [], fail, []

    def get(self, url, params=None, **kw):
        self.calls.append((url, params))
        if self.fail:
            raise HttpError("503")
        return {"targets": self.targets}


@pytest.fixture
def uniprot():
    cid.reset()
    yield
    cid.reset()


def test_a_symbol_becomes_the_single_protein_target_id(uniprot):
    c = _TargetClient([
        {"target_chembl_id": "CHEMBL6193812",
         "target_type": "PROTEIN-PROTEIN INTERACTION"},
        {"target_chembl_id": "CHEMBL3401", "target_type": "SINGLE PROTEIN"}])
    assert cid.target_id(c, "nr1i2") == "CHEMBL3401"
    params = c.calls[0][1]
    assert params["target_type"] == "SINGLE PROTEIN"
    assert params["organism"] == "Homo sapiens"
    assert params[
        "target_components__target_component_synonyms__component_synonym"
        "__iexact"] == "NR1I2"


def test_a_symbol_that_matches_two_targets_is_not_guessed(uniprot):
    c = _TargetClient([
        {"target_chembl_id": "CHEMBL1", "target_type": "SINGLE PROTEIN"},
        {"target_chembl_id": "CHEMBL2", "target_type": "SINGLE PROTEIN"}])
    assert cid.target_id(c, "AMBIG") is None


def test_complexes_alone_are_not_the_factor(uniprot):
    c = _TargetClient([{"target_chembl_id": "CHEMBL6193812",
                        "target_type": "PROTEIN-PROTEIN INTERACTION"}])
    assert cid.target_id(c, "NR1I2") is None


def test_an_empty_symbol_makes_no_request_and_an_unknown_one_is_none(uniprot):
    c = _TargetClient()
    assert cid.target_id(c, "") is None and cid.target_id(c, None) is None
    assert c.calls == []
    assert cid.target_id(c, "NOSUCH") is None and len(c.calls) == 1


def test_a_found_target_is_remembered_and_a_failure_is_not(uniprot, caplog):
    ok = _TargetClient([{"target_chembl_id": "CHEMBL3401",
                         "target_type": "SINGLE PROTEIN"}])
    cid.target_id(ok, "NR1I2")
    cid.target_id(ok, "NR1I2")
    assert len(ok.calls) == 1
    cid.reset()
    bad = _TargetClient(fail=True)
    with caplog.at_level(logging.WARNING):
        assert cid.target_id(bad, "NR1I2") is None
    assert "target lookup failed for NR1I2" in caplog.text
    cid.target_id(bad, "NR1I2")
    assert len(bad.calls) == 2, "a failure is retried, not remembered as none"
