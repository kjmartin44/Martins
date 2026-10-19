"""Ranking consistency: one scoring for one drug-factor pair, and ids that resolve.

A review of the ranking found the top of every list was a tie (26 to 30 of the
top 30 drugs shared a score with another), because half the factor drugs never
got a potency lookup and every ligand found by a receptor screen had the same
constant strength. These tests pin the fixes. They do not claim the ranking is
more accurate, which the available answer keys cannot show.
"""

from __future__ import annotations

import pytest

from promoter_tf_drugs import catalogue as cat
from promoter_tf_drugs import chemblid, drug_scoring as ds


# ----------------------------------------------------------- the helpers
def test_potency_maps_pchembl_5_to_0_and_9_to_1_and_unknown_to_neutral():
    assert cat._potency01(5.0) == 0.0 and cat._potency01(9.0) == 1.0
    assert cat._potency01(7.0) == 0.5
    assert cat._potency01(3.0) == 0.0 and cat._potency01(11.0) == 1.0
    assert cat._potency01(None) == 0.5


def test_tf_info_reads_confidence_role_and_each_drugs_direction():
    tfs = [{"symbol": "NR1I2", "score": 55.6, "role": "activator",
            "drugs": [{"name": "Rifampin", "direction": "up"},
                      {"name": "Mifepristone Hydrochloride", "direction": "down"}]},
           {"symbol": "X", "score": None, "role": None, "drugs": []}]
    info = cat.tf_info_from_tfs(tfs)
    assert info["NR1I2"]["confidence"] == pytest.approx(0.556)
    assert info["NR1I2"]["role"] == "activator"
    assert info["NR1I2"]["directions"] == {"rifampin": "up",
                                           "mifepristone": "down"}
    assert info["X"]["confidence"] == 0.0
    assert cat.tf_info_from_tfs(None) == {}


# ----------------------------------------------- screen-found ligand routes
EST = "ESTABLISHED LIGAND, untested on this gene"


def _screen(*hits, sym="NR1I2"):
    return {sym: {"enrichment": {"usable": True},
                  "hits": [{"quadrant": EST, **h} for h in hits]}}


def _info(role="activator", conf=0.6, **directions):
    return {"NR1I2": {"confidence": conf, "role": role,
                      "directions": directions}}


def _routes(hits, info, goal="up"):
    return cat.from_screen(_screen(*hits), tf_info=info, goal=goal)


def _score(r):
    return cat._score(r)


def test_more_potent_ligands_score_higher_where_all_used_to_tie():
    got = _routes([{"drug": "Strong", "measured_pchembl": 9.4},
                   {"drug": "Weak", "measured_pchembl": 6.0},
                   {"drug": "Unmeasured"}],
                  _info(strong="up", weak="up", unmeasured="up"))
    s = {r.drug: _score(r) for r in got}
    assert s["Strong"] > s["Unmeasured"] > s["Weak"]
    assert len(set(s.values())) == 3


def test_a_known_matching_action_beats_an_unrecorded_one():
    got = _routes([{"drug": "Known", "measured_pchembl": 7.0},
                   {"drug": "Unknown", "measured_pchembl": 7.0}],
                  _info(known="up"))
    by = {r.drug: r for r in got}
    assert _score(by["Known"]) > _score(by["Unknown"])
    assert any("not recorded" in n for n in by["Unknown"].notes)
    assert not any("not recorded" in n for n in by["Known"].notes)
    assert by["Known"].tier == by["Unknown"].tier == 2


def test_a_drug_known_to_move_the_factor_the_wrong_way_is_not_a_route():
    got = _routes([{"drug": "Blocker", "measured_pchembl": 8.0},
                   {"drug": "Booster", "measured_pchembl": 8.0}],
                  _info(role="activator", blocker="down", booster="up"),
                  goal="up")
    assert [r.drug for r in got] == ["Booster"]
    flipped = _routes([{"drug": "Blocker", "measured_pchembl": 8.0}],
                      _info(role="activator", blocker="down"), goal="down")
    assert [r.drug for r in flipped] == ["Blocker"]


def test_a_repressor_reverses_which_action_is_wanted():
    got = _routes([{"drug": "Blocker", "measured_pchembl": 8.0},
                   {"drug": "Booster", "measured_pchembl": 8.0}],
                  _info(role="repressor", blocker="down", booster="up"),
                  goal="up")
    assert [r.drug for r in got] == ["Blocker"]


def test_a_factor_of_unknown_role_cannot_rule_a_direction_out():
    got = _routes([{"drug": "Blocker", "measured_pchembl": 8.0}],
                  _info(role="unknown", blocker="down"), goal="up")
    assert [r.drug for r in got] == ["Blocker"]


def test_the_factors_own_confidence_moves_the_strength():
    hits = [{"drug": "D", "measured_pchembl": 8.0}]
    hi = _routes(hits, _info(conf=0.9, d="up"))[0]
    lo = _routes(hits, _info(conf=0.2, d="up"))[0]
    assert hi.axes.strength > lo.axes.strength
    none = _routes(hits, {"NR1I2": {"confidence": None, "role": "activator",
                                    "directions": {"d": "up"}}})[0]
    assert lo.axes.strength < none.axes.strength < hi.axes.strength, \
        "an unknown confidence is neutral (0.5), not zero"


def test_the_strength_is_the_graded_formula_so_the_two_collectors_agree():
    hit = {"drug": "D", "measured_pchembl": 8.0}
    r = _routes([hit], _info(conf=0.6, d="up"))[0]
    expected = round(0.45 * 0.75 + 0.30 * 1.0 + 0.25 * 0.6, 3)
    assert r.axes.strength == expected
    g = cat.from_grades([{
        "drug_name": "D", "tf_symbol": "NR1I2", "tf_confidence": 0.6,
        "components": {"potency": 0.75, "action": 1.0},
        "evidence": {"pchembl_median": 8.0, "pchembl_n": 4}, "notes": []}])[0]
    assert g.axes.strength == r.axes.strength


def test_without_tf_info_a_screen_ligand_keeps_its_old_constant_axes():
    got = cat.from_screen(_screen({"drug": "D", "measured_pchembl": 9.0}))
    assert got[0].axes.strength == 0.7 and got[0].tier == 2


def test_novel_hypotheses_are_untouched_by_the_new_scoring():
    sc = {"NR1I2": {"enrichment": {"usable": True}, "hits": [
        {"quadrant": "NOVEL HYPOTHESIS", "drug": "N", "probability": 0.8}]}}
    got = cat.from_screen(sc, tf_info=_info(), goal="up")
    assert got[0].tier == 5 and got[0].category == "virtual screen hypothesis"


def test_build_hands_the_blocks_tf_info_to_the_screen_collector():
    block = {"screens": _screen({"drug": "D", "measured_pchembl": 8.0}),
             "tf_info": _info(conf=0.9, d="up")}
    built = cat.build(block, goal="up", include_chromatin=False)
    assert built["routes"][0]["axes"]["strength"] > 0.7


# ---------------------------------------------- graded routes: measured tier
def _grade(p, action=1.0, name="D"):
    return {"drug_name": name, "tf_symbol": "T", "tf_confidence": 0.6,
            "components": {"potency": 0.5, "action": action},
            "evidence": {"pchembl_median": p, "pchembl_n": 3}, "notes": []}


def test_measured_binding_at_one_micromolar_with_a_recorded_action_is_tier_2():
    r = cat.from_grades([_grade(6.0)])[0]
    assert r.tier == 2 and r.axes.evidence == "measured binding"
    assert any("pChEMBL 6.0 against T" in n for n in r.notes)


def test_weak_unmeasured_or_actionless_drugs_stay_tier_3():
    for g in (_grade(5.9), _grade(None), _grade(8.0, action=0.4),
              _grade(8.0, action=0.69)):
        r = cat.from_grades([g])[0]
        assert r.tier == 3 and r.axes.evidence == "predicted"
        assert not any(n.startswith("measured:") for n in r.notes)


def test_a_graded_route_with_a_good_action_but_no_evidence_block_is_tier_3():
    g = _grade(7.0)
    del g["evidence"]
    assert cat.from_grades([g])[0].tier == 3


# ------------------------------------------------ evidence lookups by name
class _Chembl:
    def __init__(self):
        self.calls = []

    def get(self, url, params=None, **kw):
        self.calls.append((url, params))
        if url.endswith("/activity"):
            return {"activities": [{"pchembl_value": "6.5"},
                                   {"pchembl_value": "5.5"}]}
        if url.endswith("/mechanism"):
            return {"mechanisms": [{"target_chembl_id": "CHEMBL1"}]}
        if url.endswith("/drug_warning"):
            return {"drug_warnings": []}
        raise AssertionError(url)


class _Drug:
    def __init__(self, chembl_id, name="Rifampin", target_chembl_id=None):
        self.chembl_id, self.name = chembl_id, name
        self.target_chembl_id, self.target_symbol = target_chembl_id, "NR1I2"
        self.withdrawn, self.properties = False, {}
        self.action_type, self.max_phase = "AGONIST", 4.0


@pytest.fixture(autouse=True)
def _offline_helpers(monkeypatch):
    from promoter_tf_drugs import localannot
    monkeypatch.setattr(localannot, "uniprot_for", lambda s: None)
    chemblid.reset()
    yield
    chemblid.reset()


def test_a_drug_with_an_rxnorm_id_gets_its_evidence_through_a_found_id(
        monkeypatch):
    monkeypatch.setattr(chemblid, "resolve",
                        lambda c, n, g: ("CHEMBL374478", "approved table"))
    monkeypatch.setattr(chemblid, "target_id", lambda c, s: "CHEMBL3401")
    c = _Chembl()
    ev = ds.fetch_evidence(c, _Drug("rxcui:9384"))
    assert ev.pchembl_median == 6.0 and ev.pchembl_n == 2
    assert ev.id_note == "ChEMBL id CHEMBL374478 found by approved table"
    assert not ev.id_unresolved and ev.source == "ChEMBL"
    act = next(p for u, p in c.calls if u.endswith("/activity"))
    assert act["molecule_chembl_id"] == "CHEMBL374478"
    assert act["target_chembl_id"] == "CHEMBL3401"


def test_a_drug_whose_id_cannot_be_found_is_marked_unresolved_not_average(
        monkeypatch):
    monkeypatch.setattr(chemblid, "resolve", lambda c, n, g: (None,
                                                              "unresolved"))
    c = _Chembl()
    ev = ds.fetch_evidence(c, _Drug("ncit:C1"))
    assert ev.id_unresolved and ev.source == "unavailable"
    assert c.calls == [], "nothing is looked up under an id that is not one"
    grade = ds.grade_drug(_Drug("ncit:C1"), 0.6, ev)
    assert any("no ChEMBL id could be found" in n and "safety flags" in n
               for n in grade.notes)


def test_an_id_chembl_supplied_is_used_as_it_is(monkeypatch):
    monkeypatch.setattr(chemblid, "resolve",
                        lambda *a: pytest.fail("an id was supplied"))
    monkeypatch.setattr(chemblid, "target_id", lambda c, s: pytest.fail("x"))
    c = _Chembl()
    ev = ds.fetch_evidence(c, _Drug("CHEMBL25", target_chembl_id="CHEMBL2"))
    assert ev.id_note is None and not ev.id_unresolved
    act = next(p for u, p in c.calls if u.endswith("/activity"))
    assert act["molecule_chembl_id"] == "CHEMBL25"
    assert act["target_chembl_id"] == "CHEMBL2"


def test_a_missing_target_id_is_found_from_the_factor_symbol(monkeypatch):
    seen = {}

    def tid(client, symbol):
        seen["symbol"] = symbol
        return "CHEMBL3401"
    monkeypatch.setattr(chemblid, "target_id", tid)
    c = _Chembl()
    ev = ds.fetch_evidence(c, _Drug("CHEMBL25"))
    assert seen["symbol"] == "NR1I2" and ev.pchembl_median == 6.0


def test_the_found_id_is_reported_in_the_grade():
    ev = ds.DrugEvidence(id_note="ChEMBL id CHEMBL9 found by name search")
    g = ds.grade_drug(_Drug("rxcui:1"), 0.6, ev)
    assert "ChEMBL id CHEMBL9 found by name search" in g.notes


# -------------------------------------------------------- the pipeline wire
def test_the_pipeline_gives_the_catalogue_each_factors_confidence(monkeypatch):
    from promoter_tf_drugs import pipeline
    seen = {}
    real = cat.build

    def spy(block, **kw):
        seen.setdefault("info", block.get("tf_info"))
        return real(block, **kw)
    monkeypatch.setattr(cat, "build", spy)
    pipeline.run("SMN2", demo=True, offline=True)
    info = seen["info"]
    assert info, "demo mode has factors"
    for sym, v in info.items():
        assert 0.0 <= v["confidence"] <= 1.0
        assert set(v) == {"confidence", "role", "directions"}
