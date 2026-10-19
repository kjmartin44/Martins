"""Folding near-identical drugs into one row, and the sentence on each row.

A family of close analogues fills the top of a ranking with one idea. These
tests pin when drugs are folded and when they are not: a better-ranked drug
leads, a chain does not merge into one blob, a different target keeps two
routes apart, and a drug with no structure is never folded on a guess.
"""

from __future__ import annotations

import pytest

from promoter_tf_drugs import clusters as cl

BENZOATES = {   # methyl to hexyl benzoate: similarity falls off along the series
    "methyl": "COC(=O)c1ccccc1", "ethyl": "CCOC(=O)c1ccccc1",
    "propyl": "CCCOC(=O)c1ccccc1", "butyl": "CCCCOC(=O)c1ccccc1",
    "pentyl": "CCCCCOC(=O)c1ccccc1", "hexyl": "CCCCCCOC(=O)c1ccccc1"}


def _items(*names):
    return [{"name": n} for n in names]


def _smiles(table):
    return lambda item: table.get(item["name"])


@pytest.fixture(autouse=True)
def _fresh():
    cl.reset()
    yield
    cl.reset()


# -------------------------------------------------------------- structures
def test_a_salt_and_its_parent_have_the_same_fingerprint():
    a, b = cl.fingerprint("CCO"), cl.fingerprint("CCO.Cl")
    assert cl.tanimoto(a, b) == 1.0, "counter-ions are removed first"


def test_unrelated_drugs_are_not_similar():
    a = cl.fingerprint("CC(=O)Oc1ccccc1C(=O)O")
    b = cl.fingerprint("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
    assert cl.tanimoto(a, b) < 0.3


def test_no_structure_or_an_unreadable_one_gives_no_fingerprint():
    assert cl.fingerprint(None) is None and cl.fingerprint("") is None
    assert cl.fingerprint("not a structure ((") is None


# ---------------------------------------------------------------- grouping
def test_close_analogues_fold_into_the_best_ranked_one():
    got = cl.group(_items("pentyl", "hexyl", "methyl"), smiles_of=_smiles(
        BENZOATES))
    assert [c["leader"]["name"] for c in got] == ["pentyl", "methyl"]
    assert [m["name"] for m in got[0]["members"]] == ["hexyl"]
    assert got[0]["similarity"]["hexyl"] > 0.9


def test_the_leader_is_the_better_ranked_drug_not_the_alphabetical_one():
    got = cl.group(_items("hexyl", "pentyl"), smiles_of=_smiles(BENZOATES))
    assert got[0]["leader"]["name"] == "hexyl"
    assert [m["name"] for m in got[0]["members"]] == ["pentyl"]


def test_a_chain_of_analogues_does_not_merge_into_one_cluster():
    """Each is close to the next, but only the leader counts. Propyl joins
    ethyl at 0.70; butyl is 0.79 from propyl and only 0.66 from ethyl, and a
    member is not a leader, so butyl starts a cluster of its own."""
    got = cl.group(_items("ethyl", "propyl", "butyl"),
                   smiles_of=_smiles(BENZOATES), threshold=0.69)
    assert [c["leader"]["name"] for c in got] == ["ethyl", "butyl"]
    assert [m["name"] for m in got[0]["members"]] == ["propyl"]


def test_the_threshold_decides_what_is_close():
    items = _items("ethyl", "propyl")
    strict = cl.group(items, smiles_of=_smiles(BENZOATES), threshold=0.9)
    loose = cl.group(items, smiles_of=_smiles(BENZOATES), threshold=0.6)
    assert len(strict) == 2 and len(loose) == 1


def test_the_default_threshold_is_the_documented_conventional_one():
    assert cl.SIMILARITY == 0.75
    got = cl.group(_items("methyl", "ethyl"), smiles_of=_smiles(BENZOATES))
    assert len(got) == 2, "0.58 is related, not nearly identical"


def test_a_condition_can_keep_close_drugs_apart():
    items = [{"name": "pentyl", "target": "A"}, {"name": "hexyl",
                                                 "target": "B"}]
    apart = cl.group(items, smiles_of=_smiles(BENZOATES),
                     same_route=lambda a, b: a["target"] == b["target"])
    together = cl.group(items, smiles_of=_smiles(BENZOATES))
    assert len(apart) == 2 and len(together) == 1


def test_a_drug_with_no_structure_is_never_folded_on_a_guess():
    got = cl.group(_items("pentyl", "mystery"), smiles_of=_smiles(BENZOATES))
    assert len(got) == 2


def test_the_same_drug_under_an_ester_name_folds_without_a_structure():
    got = cl.group(_items("Hydrocortisone", "Hydrocortisone Valerate",
                          "Other"), smiles_of=lambda i: None)
    assert [c["leader"]["name"] for c in got] == ["Hydrocortisone", "Other"]
    assert got[0]["similarity"] == {"Hydrocortisone Valerate": 1.0}


def test_without_rdkit_only_names_fold_and_it_says_so(monkeypatch, caplog):
    import logging
    monkeypatch.setitem(__import__("sys").modules, "rdkit", None)
    cl._FP.clear()
    cl._WARNED["rdkit"] = False
    with caplog.at_level(logging.WARNING):
        got = cl.group(_items("pentyl", "hexyl"), smiles_of=_smiles(BENZOATES))
    assert len(got) == 2
    assert "folded by name only" in caplog.text


def test_names_are_listed_with_a_count_when_there_are_many():
    c = {"members": [{"name": n} for n in "ABCDEFGH"]}
    assert cl.names(c, limit=3) == "A, B, C and 5 more"
    assert cl.names({"members": [{"name": "A"}]}) == "A"


def test_the_approved_table_supplies_structures_by_normalised_name():
    cl.smiles_index([{"name": "RIFAMPIN", "smiles": "CC"},
                     {"name": "RIFAMPIN", "smiles": "CCC"},
                     {"name": "NOSMILES", "smiles": ""}])
    idx = cl.smiles_index()
    assert idx == {"rifampin": "CC"}, "first wins; empty structures skipped"
    got = cl.group(_items("Rifampin Sodium"))
    assert len(got) == 1


def test_no_table_gives_no_structures_and_no_error(monkeypatch):
    from promoter_tf_drugs import localdrugs
    monkeypatch.setattr(localdrugs, "library", lambda path=None: [])
    cl.reset()
    assert cl.smiles_index() == {}
    assert len(cl.group(_items("a", "b"))) == 2


def test_a_drug_exactly_at_the_threshold_is_folded():
    s = cl.tanimoto(cl.fingerprint(BENZOATES["pentyl"]),
                    cl.fingerprint(BENZOATES["hexyl"]))
    items = _items("pentyl", "hexyl")
    assert len(cl.group(items, smiles_of=_smiles(BENZOATES),
                        threshold=s)) == 1
    assert len(cl.group(items, smiles_of=_smiles(BENZOATES),
                        threshold=s + 1e-6)) == 2


def test_a_drug_close_to_two_leaders_joins_only_the_better_ranked_one():
    """Methyl and pentyl benzoate are 0.48 apart, so both lead. Ethyl is
    0.58 from methyl and 0.63 from pentyl: it goes to methyl, which ranks
    first, and is not listed under both."""
    got = cl.group(_items("methyl", "pentyl", "ethyl"),
                   smiles_of=_smiles(BENZOATES), threshold=0.55)
    assert [c["leader"]["name"] for c in got] == ["methyl", "pentyl"]
    assert [m["name"] for m in got[0]["members"]] == ["ethyl"]
    assert got[1]["members"] == []
