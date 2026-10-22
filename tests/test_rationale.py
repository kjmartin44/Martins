"""The sentence on how a candidate would act, and the folded summary rows.

The sentence states the mechanism a route assumes. It must not read as
established where the route says it is not, and must not claim to treat the
condition.
"""

from __future__ import annotations

import pytest

from generoutes import combined as cb
from generoutes import rationale as rat


def _route(cat, tf="NR1I2", mech="a mechanism", evidence="predicted", **kw):
    r = {"category": cat, "acts_on": tf, "mechanism": mech,
         "axes": {"evidence": evidence}, "drug": "D", "tier": 2,
         "score": 50.0, "source": "x"}
    r.update(kw)
    return r


def S(route, goal="up", roles=None, gene="G1"):
    return rat.route_sentence(route, gene, goal, roles or {})


# ------------------------------------------------------------ factor routes
def test_an_activator_is_activated_to_raise_and_blocked_to_lower():
    r = _route("TF: orthosteric ligand", evidence="measured binding")
    up = S(r, "up", {"NR1I2": "activator"})
    assert "switch on G1" in up and "activating NR1I2 should raise G1" in up
    down = S(r, "down", {"NR1I2": "activator"})
    assert "blocking NR1I2 should lower G1" in down


def test_a_repressor_reverses_the_action():
    r = _route("TF: orthosteric ligand")
    assert "hold down G1" in S(r, "up", {"NR1I2": "repressor"})
    assert "blocking NR1I2 should raise G1" in S(r, "up",
                                                  {"NR1I2": "repressor"})


def test_a_factor_of_unknown_role_does_not_get_an_invented_action():
    for roles in ({}, {"NR1I2": None}, {"NR1I2": "unknown"}):
        s = S(_route("TF: orthosteric ligand"), "up", roles)
        assert "role there is not known" in s
        assert "whether to activate or block it" in s
        assert "should raise" not in s


def test_measured_binding_is_said_apart_from_the_untested_step():
    m = S(_route("TF: orthosteric ligand", evidence="measured binding"),
          roles={"NR1I2": "activator"})
    assert "Its binding is measured; the step from NR1I2 to G1 is not." in m
    p = S(_route("TF: orthosteric ligand"), roles={"NR1I2": "activator"})
    assert "Neither its effect on NR1I2 nor the step" in p


def test_other_factor_routes_say_what_they_change():
    assert "how much NR1I2 the cell makes" in S(
        _route("TF: abundance change"))
    assert "activity of NR1I2 in cells, inferred from signatures" in S(
        _route("TF: inferred activity change"))
    assert "how long NR1I2 lasts" in S(
        _route("TF: degradation blocker", mech="blocks degradation"))
    assert "natural ligand" in S(_route("TF: endogenous ligand"))


# ----------------------------------------------------------- the other kinds
def test_a_drug_on_the_gene_product_changes_activity_not_amount():
    s = S(_route("direct: gene product", tf="G1", mech="antagonist"),
          goal="down")
    assert "changing its activity and not its amount" in s
    assert "(antagonist)" in s


def test_signature_upstream_chromatin_and_dna_routes():
    assert "cell-line signatures at lab doses" in S(
        _route("empirical signature"))
    assert "only by inference through the network" in S(
        _route("upstream signalling", tf="S1PR1"))
    assert "thousands of other genes" in S(_route("chromatin-level"))
    assert "unproven" in S(_route("promoter DNA binder"))
    assert "never tested at that receptor" in S(
        _route("virtual screen hypothesis"))
    assert "whether it blocks or activates it is not known" in S(
        _route("direct: predicted product binder"))


def test_an_unknown_category_falls_back_to_its_mechanism_text():
    assert S(_route("something new", mech="does a thing")) == "does a thing."
    assert "something new" in S(_route("something new", mech=""))


def test_a_docking_only_drug_says_binder_not_direction():
    gene = rat.docking_sentence("g1", "G1")
    assert "G1 protein pocket" in gene and "not known" in gene
    factor = rat.docking_sentence("JUN", "G1")
    assert "pocket of JUN" in factor and "step from JUN to G1" in factor


# ----------------------------------------------------------------- the aim
def test_the_aim_states_the_dosage_requirement_when_there_is_one():
    req = {"mechanism": "one allele lost, raising the other",
           "requirement": "about 2x from the working allele"}
    a = rat.aim("BUB1B", "up", req)
    assert "one allele lost, raising the other" in a
    assert "about 2x from the working allele" in a


def test_the_declared_scenario_is_used_when_there_is_no_requirement():
    a = rat.aim("G1", "down", None, None, "toxic_gain")
    assert a.startswith("For G1, the harmful product comes from the gene")
    assert rat.aim("G1", "up", None, None, "no such scenario").startswith(
        "The aim is to raise G1.")


def test_the_aim_without_a_scenario_names_the_goal_and_condition():
    assert rat.aim("G1", "down").startswith("The aim is to lower G1.")
    assert "to help with asthma" in rat.aim("G1", "up", None, "asthma")


def test_the_aim_never_claims_the_drug_treats_the_condition():
    for a in (rat.aim("G1", "up"), rat.aim("G1", "down", None, "x")):
        assert "Nothing here shows that moving the gene treats the " \
               "condition" in a
        assert "has not been measured" in a


def test_there_is_a_sentence_for_every_scenario_the_launcher_offers():
    assert set(rat._SCENARIO_HELP) == {
        "haploinsufficiency", "hypomorph", "toxic_gain", "null",
        "silenced_allele", "modifier"}


# --------------------------------------------- folding in the summary report
def _demo_with_a_family():
    from generoutes import pipeline
    res = pipeline.run("SMN2", demo=True, offline=True)
    blk = res["results"][0]
    routes = blk["all_routes"]["routes"]
    top = next(r for r in routes if r.get("acts_on"))
    twin = dict(top, drug=top["drug"] + " Valerate",
                score=top["score"] - 0.1)
    routes.append(twin)
    routes.sort(key=lambda r: (r["tier"], -r["score"]))
    return res, top, twin


def test_the_summary_folds_a_family_and_names_the_folded_drug():
    from generoutes import report
    res, top, twin = _demo_with_a_family()
    text = report.summary_report(res)
    assert "| How it would work |" in text
    assert "(+1 similar)" in text
    assert f"{top['drug']}: {twin['drug']}" in text, "the folded drug is named"
    shortlist = text.split("## What to test first")[1].split("**Every row is off-label")[0]
    rows = [ln for ln in shortlist.splitlines()
            if ln.startswith("| ") and f"**{twin['drug']}**" in ln]
    assert not rows, "the folded drug has no row of its own in the shortlist"


def test_the_summary_states_the_aim_and_a_sentence_per_row():
    from generoutes import report
    res, *_ = _demo_with_a_family()
    text = report.summary_report(res)
    assert "SMN2" in text.split("## What to test first")[1].split("### ")[0]
    assert "Nothing here shows that moving the gene treats the condition" \
        in text


def test_a_drug_reached_through_two_factors_is_one_row_with_both_targets():
    from generoutes import report
    res, top, _ = _demo_with_a_family()
    other = dict(top, acts_on="OTHERTF", score=top["score"] - 5)
    res["results"][0]["all_routes"]["routes"].append(other)
    text = report.summary_report(res)
    shortlist = text.split("## What to test first")[1].split("**Every row is off-label")[0]
    rows = [ln for ln in shortlist.splitlines()
            if ln.startswith("| ") and f"**{top['drug']}**" in ln]
    assert len(rows) == 1
    assert "OTHERTF" in rows[0] and top["acts_on"] in rows[0]


# ----------------------------------------------- folding in the combined list
def _crow(drug, tier=2, score=50.0, source="screening"):
    return {"drug": drug, "tier": tier, "score": score, "source": source,
            "category": "TF: orthosteric ligand", "acts_on": "T",
            "mechanism": "m", "axes": {"evidence": "predicted"}}


@pytest.fixture
def structures():
    from generoutes import clusters
    clusters.smiles_index([
        {"name": "PENTYL", "smiles": "CCCCCOC(=O)c1ccccc1"},
        {"name": "HEXYL", "smiles": "CCCCCCOC(=O)c1ccccc1"},
        {"name": "ASPIRIN", "smiles": "CC(=O)Oc1ccccc1C(=O)O"}])
    yield
    clusters.reset()


def _results(*rows):
    return {"query": {"gene": "G1", "goal": "up"},
            "results": [{"all_routes": {"routes": list(rows)}}]}


def test_the_combined_list_folds_by_structure_and_keeps_every_rank(
        structures):
    res = cb.build(results=_results(_crow("Pentyl", 2, 60),
                                    _crow("Aspirin", 2, 55),
                                    _crow("Hexyl", 2, 50)))
    assert [t["drug"] for t in res["table"]] == ["Pentyl", "Aspirin"]
    assert res["table"][0]["n_similar"] == 1
    assert res["table"][0]["similar"] == "Hexyl"
    f = res["folded"][0]
    assert f["leader"] == "Pentyl" and f["members"][0]["name"] == "Hexyl"
    assert f["members"][0]["rank"] > f["leader_rank"]
    assert res["n_folded"] == 1 and res["n_clusters_folding"] == 1
    md = cb.markdown(res)
    assert "(+1 similar)" in md and "Folded into the rows above" in md
    assert "How it would work" in md


def test_the_combined_list_can_list_every_drug(structures):
    results = _results(_crow("Pentyl", 2, 60), _crow("Hexyl", 2, 50))
    assert len(cb.build(results=results, fold=False)["table"]) == 2
    assert len(cb.build(results=results)["table"]) == 1
    loose = cb.build(results=results, similarity=0.99)
    assert len(loose["table"]) == 2, "a stricter threshold folds fewer"


def test_every_combined_row_has_a_sentence_and_docking_only_rows_too():
    res = cb.build(
        results=_results(_crow("A")),
        dockings=[("G1", {"usable": True, "auc": 0.9, "ef1": 5,
                          "n_docked": 10, "hits": [
                              {"drug": "Z", "score": -9.0,
                               "cnn_score": 0.5}]})])
    how = {t["drug"]: t["how"] for t in res["table"]}
    assert "Binds T" in how["A"] and "G1 protein pocket" in how["Z"]


def test_the_command_line_checks_the_similarity_option(capsys):
    with pytest.raises(SystemExit):
        cb.main(["--results", "x", "--similarity", "1.5"])
    assert "above 0 and at most 1" in capsys.readouterr().err


def test_close_analogues_acting_through_different_targets_stay_two_rows():
    from generoutes import report
    res, top, twin = _demo_with_a_family()
    twin["acts_on"] = "SOMEOTHERTF"
    text = report.summary_report(res)
    shortlist = text.split("## What to test first")[1].split("**Every row is off-label")[0]
    rows = {ln for ln in shortlist.splitlines() if ln.startswith("| ")}
    mine = [ln for ln in rows if f"**{twin['drug']}**" in ln
            or f"**{top['drug']}**" in ln]
    assert len(mine) == 2, "a different route is kept as its own row"
    assert not any("similar" in ln for ln in mine)


def test_the_summary_uses_the_declared_scenario_when_no_requirement_was_made():
    from generoutes import report
    res, *_ = _demo_with_a_family()
    res["requirement"] = None
    res["query"]["scenario"] = "toxic_gain"
    shortlist = report.summary_report(res).split("## What to test first")[1]
    assert "the harmful product comes from the gene" in shortlist


def test_the_combined_list_takes_its_goal_from_the_panel_when_given_only_that():
    panel = {"query": {"genes": ["G1"], "goal": "down"},
             "surviving": [dict(_crow("A"), gene="G1")]}
    res = cb.build(panel=panel)
    assert res["goal"] == "down"
    assert "lower G1" in res["table"][0]["how"] or "block" in \
        res["table"][0]["how"]


# ------------------------------------------- one row per mechanism family
def _fam_route(drug, tier, score, acts_on="VDR", category="TF: orthosteric "
               "ligand"):
    return {"drug": drug, "tier": tier, "score": score, "source": "t",
            "category": category, "acts_on": acts_on, "mechanism": "m",
            "axes": {"evidence": "predicted"}}


def _shortlist_rows(text):
    part = text.split("## What to test first")[1].split("**Every row is off-label")[0]
    return [ln for ln in part.splitlines()
            if ln.startswith("| ") and "**" in ln]


def test_unrelated_drugs_through_one_factor_take_one_row(structures):
    from generoutes import report
    text = report.summary_report(_results(
        _fam_route("PENTYL", 3, 60), _fam_route("ASPIRIN", 3, 50),
        _fam_route("OTHER", 3, 40, acts_on="RXRA")))
    rows = _shortlist_rows(text)
    assert len(rows) == 2
    # Position carries nothing now: rows inside a group are ordered by
    # independent method count, then alphabetically, never by score.
    pentyl = next(r for r in rows if "**PENTYL**" in r)
    assert "(+1 similar)" in pentyl
    assert "PENTYL: ASPIRIN" in text


def test_the_same_factor_by_a_different_kind_of_route_is_another_row(
        structures):
    from generoutes import report
    text = report.summary_report(_results(
        _fam_route("PENTYL", 3, 60),
        _fam_route("ASPIRIN", 3, 50, category="TF: known regulator")))
    assert len(_shortlist_rows(text)) == 2


def test_drugs_with_no_stated_factor_are_never_folded_as_a_family(structures):
    from generoutes import report
    lit = "literature: reported effect"
    text = report.summary_report(_results(
        _fam_route("PENTYL", 4, 60, acts_on="not stated", category=lit),
        _fam_route("ASPIRIN", 4, 50, acts_on="not stated", category=lit)))
    assert len(_shortlist_rows(text)) == 2


def test_folding_never_drops_a_drug_from_the_summary(structures):
    """Every drug that reaches the pool is a row or is named under the table."""
    from generoutes import report
    drugs = ["PENTYL", "ASPIRIN", "HEXYL", "OTHER"]
    text = report.summary_report(_results(
        _fam_route("PENTYL", 3, 60), _fam_route("ASPIRIN", 3, 50),
        _fam_route("HEXYL", 3, 45), _fam_route("OTHER", 3, 40,
                                              acts_on="RXRA")))
    rows = " ".join(_shortlist_rows(text))
    folded = text.split("**Folded into the rows above**")[1]
    for d in drugs:
        assert f"**{d}**" in rows or d in folded, f"{d} vanished"


def test_a_folded_drug_stays_in_the_full_results(structures):
    from generoutes import report
    res = _results(_fam_route("PENTYL", 3, 60), _fam_route("ASPIRIN", 3, 50))
    report.summary_report(res)
    kept = [r["drug"] for r in res["results"][0]["all_routes"]["routes"]]
    assert kept == ["PENTYL", "ASPIRIN"], "folding only changes the summary"
