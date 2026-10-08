"""The search for druggable nodes upstream of a transcription factor.

It used to ask about drugs for only the 25 best-scoring nodes. Those are
mostly chromatin and DNA-repair proteins with no approved drug, so a factor
whose first druggable node ranked 103rd of 1,189 was reported as having no
druggable route at all. These tests pin the walk past the top, the batched
index that makes it cheap, and what the ledger says about how far it got.
"""

from __future__ import annotations

import gzip
import random

import pytest

from promoter_tf_drugs import drugs as drugs_mod
from promoter_tf_drugs import network as net
from promoter_tf_drugs.drugs import Drug
from promoter_tf_drugs.http import HttpError

NL = chr(10)
TAB = chr(9)


def _lf_drug(name, symbol, direction="up", action="AGONIST"):
    return Drug(name=name, chembl_id="CHEMBL" + name[:3].upper(),
                target_symbol=symbol, action_type=action, mechanism=None,
                approved=True, max_phase=4.0, first_approval=None,
                withdrawn=False, oral=None, direction=direction,
                source="test")


def _lf_star(n, tf="TFX"):
    """`n` nodes one signed hop upstream of the factor, all activating."""
    return [net.Edge(f"NODE{i:03d}", tf, 1, 2, "test") for i in range(n)]


# ------------------------------------------------- the walk past the top
def test_a_druggable_node_far_down_the_ranking_is_found(monkeypatch):
    calls = []

    def fake(client, symbol, uniprot, approved_only):
        calls.append(symbol)
        return [_lf_drug("Rescue", symbol)] if symbol == "NODE090" else []
    monkeypatch.setattr(drugs_mod, "chembl_drugs", fake)
    found = net.find_druggable_upstream(None, "TFX", goal="up",
                                        tf_role="activator",
                                        edges=_lf_star(120))
    assert [h.node for h in found.hits] == ["NODE090"]
    assert found.n_reachable == 120
    assert found.n_checked == net.UPSTREAM_MAX_CHECKED == 100
    assert found.truncated is True, (
        "it stopped at the limit with more reachable, and said so")
    assert found.hits[0].drugs[0]["name"] == "Rescue"


def test_the_walk_is_bounded_and_says_it_stopped_short(monkeypatch):
    monkeypatch.setattr(drugs_mod, "chembl_drugs",
                        lambda c, s, u, a: [_lf_drug("Far", s)]
                        if s == "NODE140" else [])
    found = net.find_druggable_upstream(None, "TFX", goal="up",
                                        tf_role="activator",
                                        edges=_lf_star(150), max_checked=100)
    assert found.hits == []
    assert found.n_checked == 100 and found.n_reachable == 150
    assert found.truncated is True, "a route past the limit is not ruled out"


def test_the_walk_stops_once_it_has_enough(monkeypatch):
    asked = []
    monkeypatch.setattr(
        drugs_mod, "chembl_drugs",
        lambda c, s, u, a: asked.append(s) or [_lf_drug("D" + s[-3:], s)])
    found = net.find_druggable_upstream(None, "TFX", goal="up",
                                        tf_role="activator",
                                        edges=_lf_star(60), want=3)
    assert len(found.hits) == 3 and found.n_checked == 3
    assert len(asked) == 3 and found.truncated is False


def test_a_node_asked_about_twice_is_looked_up_once(monkeypatch):
    asked = []
    monkeypatch.setattr(
        drugs_mod, "chembl_drugs",
        lambda c, s, u, a: asked.append(s) or [])
    cache = {}
    edges = _lf_star(5, "TFA") + _lf_star(5, "TFB")
    net.find_druggable_upstream(None, "TFA", goal="up", tf_role="activator",
                                edges=edges, drug_cache=cache)
    net.find_druggable_upstream(None, "TFB", goal="up", tf_role="activator",
                                edges=edges, drug_cache=cache)
    assert sorted(asked) == sorted(set(asked)) and len(asked) == 5


def test_an_outage_is_counted_as_failures_not_as_no_drug(monkeypatch):
    def down(client, symbol, uniprot, approved_only):
        raise HttpError("ChEMBL is down")
    monkeypatch.setattr(drugs_mod, "chembl_drugs", down)
    found = net.find_druggable_upstream(None, "TFX", goal="up",
                                        tf_role="activator",
                                        edges=_lf_star(40))
    assert found.hits == []
    assert found.n_lookup_failed == net.UPSTREAM_MAX_CONSECUTIVE_FAILURES
    assert found.n_checked == found.n_lookup_failed, "it stopped asking"
    assert "service looks unavailable" in found.stopped


def test_a_failed_lookup_is_not_remembered_as_an_answer(monkeypatch):
    state = {"up": False}

    def flaky(client, symbol, uniprot, approved_only):
        if not state["up"]:
            raise HttpError("503")
        return [_lf_drug("Later", symbol)]
    monkeypatch.setattr(drugs_mod, "chembl_drugs", flaky)
    cache = {}
    first = net.find_druggable_upstream(None, "TFX", goal="up",
                                        tf_role="activator",
                                        edges=_lf_star(3), drug_cache=cache)
    assert first.hits == [] and cache == {}
    state["up"] = True
    second = net.find_druggable_upstream(None, "TFX", goal="up",
                                         tf_role="activator",
                                         edges=_lf_star(3), drug_cache=cache)
    assert len(second.hits) == 3


def test_attach_drugs_still_works_and_says_when_lookups_failed(monkeypatch,
                                                              caplog):
    import logging
    monkeypatch.setattr(
        drugs_mod, "chembl_drugs",
        lambda c, s, u, a: (_ for _ in ()).throw(HttpError("x")))
    hits = net.find_upstream_targets(None, "TFX", goal="up",
                                     tf_role="activator", edges=_lf_star(3))
    with caplog.at_level(logging.WARNING):
        net.attach_drugs(None, hits)
    assert any("drug lookup failed for 3 of 3" in r.message
               for r in caplog.records)


# ------------------------------------------------------ the batched index
def _lf_library(n):
    return [{"chembl_id": f"CHEMBL{i}", "name": f"drug {i}",
             "withdrawn": "1" if i == 2 else ""} for i in range(n)]


class _LfChembl:
    """Answers the two batched ChEMBL endpoints the index reads."""

    def __init__(self, mechanisms, targets, fail_on=None):
        self.mechanisms, self.targets = mechanisms, targets
        self.calls, self.fail_on = [], fail_on

    def get(self, url, params=None, **kw):
        self.calls.append((url.rsplit("/", 1)[-1], params))
        if self.fail_on and self.fail_on in url:
            raise HttpError("500")
        if url.endswith("/mechanism"):
            ids = set(params["molecule_chembl_id__in"].split(","))
            return {"mechanisms": [m for m in self.mechanisms
                                   if m["molecule_chembl_id"] in ids]}
        ids = set(params["target_chembl_id__in"].split(","))
        return {"targets": [t for t in self.targets
                            if t["target_chembl_id"] in ids]}


def _lf_target(tid, symbol, organism="Homo sapiens", kind="SINGLE PROTEIN"):
    return {"target_chembl_id": tid, "organism": organism,
            "target_type": kind,
            "target_components": [{"target_component_synonyms": [
                {"component_synonym": symbol, "syn_type": "GENE_SYMBOL"},
                {"component_synonym": "OLD" + symbol,
                 "syn_type": "GENE_SYMBOL_OTHER"},
                {"component_synonym": "Some protein name",
                 "syn_type": "UNIPROT"}]}]}


def test_the_index_maps_gene_symbols_to_approved_drugs():
    mechs = [
        {"molecule_chembl_id": "CHEMBL0", "target_chembl_id": "T1",
         "action_type": "AGONIST", "mechanism_of_action": "agonises"},
        {"molecule_chembl_id": "CHEMBL1", "target_chembl_id": "T1",
         "action_type": "INHIBITOR", "mechanism_of_action": "blocks"},
        {"molecule_chembl_id": "CHEMBL1", "target_chembl_id": "T1",
         "action_type": "INHIBITOR", "mechanism_of_action": "blocks"},
        {"molecule_chembl_id": "CHEMBL2", "target_chembl_id": "T2",
         "action_type": "BINDING AGENT", "mechanism_of_action": "x"},
        {"molecule_chembl_id": "CHEMBL3", "target_chembl_id": "TMOUSE",
         "action_type": "AGONIST"},
        {"molecule_chembl_id": "CHEMBL3", "target_chembl_id": "TFAM",
         "action_type": "AGONIST"},
        {"molecule_chembl_id": "NOTINLIBRARY", "target_chembl_id": "T1",
         "action_type": "AGONIST"},
    ]
    targets = [_lf_target("T1", "esr1"), _lf_target("T2", "ABC1"),
               _lf_target("TMOUSE", "Esr1", organism="Mus musculus"),
               _lf_target("TFAM", "XYZ", kind="PROTEIN FAMILY")]
    client = _LfChembl(mechs, targets)
    idx = drugs_mod.approved_target_index(client, _lf_library(4))
    assert sorted(idx) == ["ABC1", "ESR1"], (
        "human single proteins only, by GENE_SYMBOL, upper-cased")
    assert len(idx["ESR1"]) == 2, "the repeated mechanism is listed once"
    esr1 = {d.name: d for d in idx["ESR1"]}
    assert sorted(esr1) == ["Drug 0", "Drug 1"]
    assert esr1["Drug 0"].direction == "up" and esr1["Drug 1"].direction == "down"
    assert idx["ABC1"][0].direction == "unclear"
    assert esr1["Drug 0"].approved and esr1["Drug 0"].source.startswith("ChEMBL")
    assert esr1["Drug 0"].target_symbol == "ESR1"


def test_the_index_reads_fifty_drugs_per_request():
    client = _LfChembl([], [])
    drugs_mod.approved_target_index(client, _lf_library(120))
    assert [c[0] for c in client.calls] == ["mechanism"] * 3
    sizes = [len(c[1]["molecule_chembl_id__in"].split(",")) for c in client.calls]
    assert sizes == [50, 50, 20]


def test_the_index_carries_the_withdrawn_flag_from_the_library():
    mechs = [{"molecule_chembl_id": "CHEMBL2", "target_chembl_id": "T1",
              "action_type": "AGONIST"}]
    idx = drugs_mod.approved_target_index(
        _LfChembl(mechs, [_lf_target("T1", "G1")]), _lf_library(4))
    assert idx["G1"][0].withdrawn is True


def test_a_failed_batch_raises_so_the_index_is_never_silently_partial():
    with pytest.raises(HttpError):
        drugs_mod.approved_target_index(
            _LfChembl([], [], fail_on="mechanism"), _lf_library(10))


def test_with_an_index_every_reachable_node_is_checked_for_free():
    index = {"NODE140": [_lf_drug("Far", "NODE140")]}
    found = net.find_druggable_upstream(
        None, "TFX", goal="up", tf_role="activator", edges=_lf_star(150),
        max_checked=10, drug_index=index)
    assert [h.node for h in found.hits] == ["NODE140"]
    assert found.n_checked == found.n_reachable == 150
    assert found.truncated is False, "the cap does not apply to an index"


def test_the_index_respects_the_direction_the_path_needs():
    """A repressor's upstream node has to be moved the other way."""
    edges = [net.Edge("KIN", "REP", -1, 2, "t")]
    both = {"KIN": [_lf_drug("Raiser", "KIN", "up"),
                    _lf_drug("Lowerer", "KIN", "down")]}
    found = net.find_druggable_upstream(None, "REP", goal="up",
                                        tf_role="repressor", edges=edges,
                                        drug_index=both)
    assert {d["direction"] for d in found.hits[0].drugs} == {"up"}


def test_lowering_the_gene_through_an_activator_is_reached_by_an_inhibitory_edge():
    """To lower a gene through an activator the activator has to go down,
    so only a node that INHIBITS it counts, and that node has to be raised.
    An activating edge would raise the gene, so it is not a path at all."""
    both = {"KIN": [_lf_drug("Raiser", "KIN", "up"),
                    _lf_drug("Lowerer", "KIN", "down")]}
    inhibits = [net.Edge("KIN", "ACT", -1, 2, "t")]
    found = net.find_druggable_upstream(None, "ACT", goal="down",
                                        tf_role="activator", edges=inhibits,
                                        drug_index=both)
    assert {d["direction"] for d in found.hits[0].drugs} == {"up"}
    activates = [net.Edge("KIN", "ACT", 1, 2, "t")]
    none = net.find_druggable_upstream(None, "ACT", goal="down",
                                       tf_role="activator", edges=activates,
                                       drug_index=both)
    assert none.hits == [] and none.n_reachable == 0


# ----------------------------------------- what a whole run says about it
def _lf_offline_run(tmp_path, monkeypatch, *, edges, library, index=None,
                    index_error=None, chembl=None):
    from promoter_tf_drugs import eqtl, genes, neighbours
    from promoter_tf_drugs import localdrugs, pipeline

    monkeypatch.setattr(genes, "_TX_BY_GENE", {})
    monkeypatch.setattr(genes, "_TX_SOURCE", {"path": None})
    monkeypatch.setattr(genes, "_GENOME", {"path": None, "index": {}})
    monkeypatch.setattr(neighbours, "_GTF_GENES", [])
    monkeypatch.setattr(neighbours, "_GTF_SOURCE", None)
    monkeypatch.setattr(eqtl, "_CACHE", {"rows": None, "targets": None})
    monkeypatch.setattr(net, "fetch_network", lambda client, **k: edges)
    monkeypatch.setattr(localdrugs, "library", lambda *a, **k: library)

    def build(client, lib):
        if index_error:
            raise HttpError(index_error)
        return index
    monkeypatch.setattr(drugs_mod, "approved_target_index", build)
    if chembl is not None:
        monkeypatch.setattr(drugs_mod, "chembl_drugs", chembl)

    rng = random.Random(1)
    seq = "".join(rng.choice("ACGT") for _ in range(6000))
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1" + NL + NL.join(seq[i:i + 60]
                                         for i in range(0, 6000, 60)) + NL,
                  encoding="utf-8")
    gtf = tmp_path / "g.gtf"
    gtf.write_text(NL.join([
        TAB.join(["chr1", "H", "gene", "2001", "4000", ".", "+", ".",
                  'gene_id "G1"; gene_name "TESTG"; '
                  'gene_type "protein_coding";']),
        TAB.join(["chr1", "H", "transcript", "2001", "4000", ".", "+", ".",
                  'gene_id "G1"; gene_name "TESTG"; '
                  'transcript_id "ENST1.1"; '
                  'transcript_type "protein_coding"; tag "MANE_Select"; '
                  'tag "Ensembl_canonical";'])]) + NL, encoding="utf-8")
    pfm = tmp_path / "m.txt"
    pfm.write_text(NL.join([
        ">MA0001.1" + TAB + "TEST", "A [ 20 0 0 0 ]", "C [ 0 20 0 0 ]",
        "G [ 0 0 20 0 ]", "T [ 0 0 0 20 ]"]) + NL, encoding="utf-8")
    res = pipeline.run("TESTG", gtf=str(gtf), genome=str(fa),
                       pfm_file=str(pfm), offline=True, goal="down",
                       upstream_network=True, include_direct=False,
                       upstream=500, downstream=100,
                       cache_dir=str(tmp_path / "cache"), grade_drugs=False,
                       estimate_outcome=False, find_improvements=False,
                       compare_mechanisms=False, include_stabilizers=False)
    src = res["evidence_sources"]["sources"]
    rows = {r["source"]: r for r in src} if isinstance(src, list) else src
    return res["results"][0], rows["upstream"]


def _lf_tf_edges(n=6):
    return [net.Edge(f"KIN{i}", "TEST", 1, 3, "t") for i in range(n)]


def test_a_run_reports_a_found_route_and_how_much_it_covered(tmp_path,
                                                            monkeypatch):
    index = {"KIN3": [_lf_drug("Hit", "KIN3", "up")]}
    block, row = _lf_offline_run(
        tmp_path, monkeypatch, edges=_lf_tf_edges(),
        library=[{"chembl_id": "C1", "name": "x"}], index=index)
    assert row["ok"] is True and row["yielded"] == 1
    assert "6 of 6 reachable upstream node(s) checked" in row["detail"]
    assert [h["node"] for h in block["upstream"]["TEST"]] == ["KIN3"]
    assert block["upstream_search"]["drug_index"] is True
    assert block["upstream_search"]["nodes_checked"] == 6


def test_a_run_with_no_druggable_node_says_among_how_many(tmp_path,
                                                          monkeypatch):
    block, row = _lf_offline_run(
        tmp_path, monkeypatch, edges=_lf_tf_edges(),
        library=[{"chembl_id": "C1", "name": "x"}], index={})
    assert row["ok"] is True and row["yielded"] == 0
    assert "no node checked has an approved drug" in row["detail"]
    assert "6 of 6" in row["detail"]
    assert block["upstream"] == {}


def test_a_run_without_the_network_says_so_not_that_no_path_exists(
        tmp_path, monkeypatch):
    block, row = _lf_offline_run(
        tmp_path, monkeypatch, edges=[],
        library=[{"chembl_id": "C1", "name": "x"}], index={})
    assert row["ok"] is False
    assert "could not be loaded" in row["detail"]


def test_if_the_index_cannot_be_built_the_search_falls_back_and_is_capped(
        tmp_path, monkeypatch, caplog):
    import logging
    asked = []

    def chembl(client, symbol, uniprot, approved_only):
        asked.append(symbol)
        return [_lf_drug("Hit", symbol, "up")] if symbol == "KIN2" else []
    with caplog.at_level(logging.WARNING):
        block, row = _lf_offline_run(
            tmp_path, monkeypatch, edges=_lf_tf_edges(),
            library=[{"chembl_id": "C1", "name": "x"}],
            index_error="ChEMBL 500", chembl=chembl)
    assert any("could not build the approved-drug index" in r.message
               for r in caplog.records)
    assert block["upstream_search"]["drug_index"] is False
    assert row["ok"] is True and "KIN2" in asked
    assert [h["node"] for h in block["upstream"]["TEST"]] == ["KIN2"]


def test_a_run_whose_lookups_all_fail_is_not_reported_as_a_clean_negative(
        tmp_path, monkeypatch):
    def down(client, symbol, uniprot, approved_only):
        raise HttpError("503")
    block, row = _lf_offline_run(
        tmp_path, monkeypatch, edges=_lf_tf_edges(),
        library=[], chembl=down)
    assert row["ok"] is False
    assert "says nothing about whether a druggable route exists" \
        in row["detail"]
    assert block["upstream_search"]["lookups_failed"] > 0
