"""The gates added after the BUB1B review.

Each one exists because the BUB1B run was confidently wrong in a specific
way, and each test names that way. They are offline: the one that needs a
LINCS signature file writes a small one.
"""

from __future__ import annotations

import json

from promoter_tf_drugs import catalogue as cat
from promoter_tf_drugs import crosscheck as xc
from promoter_tf_drugs import magnitude as mag
from promoter_tf_drugs import proliferation as pr
from promoter_tf_drugs import regulators as reg
from promoter_tf_drugs import report as rep
from promoter_tf_drugs import systemic as sy


def _grade(name="D", tf="T"):
    return {"drug_name": name, "tf_symbol": tf, "tf_confidence": 0.6,
            "components": {"potency": 0.5, "action": 1.0},
            "evidence": {"pchembl_median": 6.0, "pchembl_n": 3},
            "notes": []}


# ===================================================================== #
# C1. The regulators with no drug are the biology.                      #
# ===================================================================== #
def test_recorded_regulators_with_no_drug_are_reported_not_dropped():
    """CollecTRI names FOXM1, MYC, TP53 and ZNF143 for BUB1B. None has an
    approved drug, so the regulator section showed AHR and three HDACs
    instead - factors that are there because somebody made a ligand."""
    regs = [reg.Regulator("FOXM1", 1, ["CollecTRI"]),
            reg.Regulator("ZNF143", 1, ["CollecTRI"]),
            reg.Regulator("ESR2", 1, ["CollecTRI"])]

    class D:
        direction = "up"

        def to_dict(self):
            return {"name": "Bazedoxifene", "direction": "up"}

    index = {"ESR2": [D()]}
    drugged = reg.find_drugs(regs, "up", index)
    undrugged = reg.without_drugs(regs, "up", index)

    assert set(drugged) == {"ESR2"}, "unchanged: only the druggable one"
    assert [r["symbol"] for r in undrugged] == ["FOXM1", "ZNF143"]
    assert all("drugs" not in r for r in undrugged), "nothing to offer"
    assert undrugged[0]["sources"] == ["CollecTRI"]


def test_the_undrugged_regulators_reach_the_report():
    rows = [{"symbol": "FOXM1", "sign": 1, "need": "up",
             "sources": ["CollecTRI"], "via": None}]
    out = "\n".join(rep._undrugged_regulators(
        {"regulators_without_drugs": rows}))
    assert "FOXM1" in out and "CollecTRI" in out
    assert "no approved drug" in out
    assert rep._undrugged_regulators({}) == []


# ===================================================================== #
# C2. Did the compound move the gene, or move the cell cycle?           #
# ===================================================================== #
def _gmt(tmp_path, rows):
    """A small Enrichr LINCS .gmt. Names follow the real shape,
    `CPC001 HA1E 24H-<drug>-10.0`, because that is what the parser reads.
    """
    p = tmp_path / "sig.gmt"
    p.write_text("".join(
        f"CPC001 HA1E 24H-{drug}-{i}.0\t\t" + "\t".join(genes) + "\n"
        for i, (drug, genes) in enumerate(rows)), encoding="utf-8")
    return p


def test_a_compound_is_judged_against_the_background_not_a_flat_share(
        tmp_path):
    """The background is not zero: 1.9% of all LINCS signatures move three
    or more mitotic markers on their own. A flat threshold called that
    signal; the bar has to be what other compounds do."""
    programme = ["AURKA", "CDK1", "PLK1", "BUB1B"]
    quiet = ["ALB", "GAPDH", "ACTB", "BUB1B"]
    rows = [("quiet", quiet) for _ in range(40)]
    rows += [("tracker", programme) for _ in range(10)]
    path = _gmt(tmp_path, rows)

    got = pr.assess("BUB1B", ["tracker", "quiet"], direction="up",
                    up_path=path)
    assert got["tracker"]["tracks_proliferation"] is True
    assert got["quiet"]["tracks_proliferation"] is False
    # BUB1B itself must never count towards its own confound.
    assert "BUB1B" not in got["tracker"]["markers"]


def test_a_compound_with_too_few_signatures_gets_no_verdict(tmp_path):
    programme = ["AURKA", "CDK1", "PLK1"]
    rows = [("quiet", ["ALB", "GAPDH", "ACTB"]) for _ in range(40)]
    rows += [("rare", programme)]
    got = pr.assess("BUB1B", ["rare"], direction="up",
                    up_path=_gmt(tmp_path, rows))
    assert got["rare"]["tracks_proliferation"] is False, "1 signature is not evidence"


def test_a_measured_proliferation_effect_halves_the_score():
    rec = {"tracks_proliferation": True, "n_with_programme": 40,
           "n_signatures": 100, "share": 0.4, "background_share": 0.019,
           "markers": ["AURKA", "CDK1", "PLK1"], "percentile": 95}
    blk = {"drug_grades": [_grade()], "gene_symbol": "BUB1B"}
    plain = cat.build(dict(blk, drug_grades=[_grade()]), goal="up")
    cut = cat.build(dict(blk, drug_grades=[_grade()]), goal="up",
                    proliferation={"d": rec})
    assert cut["routes"][0]["score"] == round(
        plain["routes"][0]["score"] * cat.MITOTIC_SCORE_FACTOR, 1)
    assert cut["routes"][0]["proliferation_measured"] is True


# ===================================================================== #
# C3. What KIND of chemical moves this gene?                            #
# ===================================================================== #
def test_the_literature_profile_says_when_the_movers_are_not_medicines():
    """Of 59 chemicals CTD records as raising BUB1B, most are carcinogens,
    tumour promoters and industrial chemicals. The run held that list and
    used it only to decide which candidates were novel."""
    names = ["1,2-Dimethylhydrazine", "Methylcholanthrene", "Acrylamide",
             "Tetradecanoylphorbol Acetate", "Croton Oil", "Nickel",
             "Lipopolysaccharides", "Titanium Dioxide", "Chlorpyrifos",
             "Estradiol"]
    prof = xc.mover_profile(names, "up")
    assert prof["n"] == 10
    assert prof["n_hazard"] == 10
    assert "not a medicine" in prof["verdict"]
    assert prof["counts"]["carcinogen or genotoxin"] >= 3
    assert prof["counts"]["tumour promoter"] >= 2


def test_an_unmatched_chemical_is_unclassified_and_not_called_safe():
    prof = xc.mover_profile(["Zzzyxine", "Acrylamide"], "up")
    assert prof["counts"]["unclassified"] == 1
    assert prof["n_hazard"] == 1


def test_no_movers_gives_no_profile():
    assert xc.mover_profile([], "up") == {}


# ===================================================================== #
# C4. A route cannot be "measured" over a site list that is noise.      #
# ===================================================================== #
def test_a_noisy_scan_caps_the_tier_of_routes_that_rest_on_it():
    """46% of the BUB1B sites are expected on shuffled sequence, and the
    run still produced tier-2 routes through factors found by that scan."""
    clean = cat.build({"drug_grades": [_grade()]}, goal="up", scan_fdp=0.10)
    noisy = cat.build({"drug_grades": [_grade()]}, goal="up", scan_fdp=0.46)
    assert clean["routes"][0]["tier"] == 2
    assert noisy["routes"][0]["tier"] == cat.NOISY_SCAN_TIER
    assert noisy["routes"][0]["tier_before_scan_penalty"] == 2
    assert any("46%" in n for n in noisy["routes"][0]["notes"])


def test_an_unknown_scan_quality_changes_nothing():
    a = cat.build({"drug_grades": [_grade()]}, goal="up")
    b = cat.build({"drug_grades": [_grade()]}, goal="up", scan_fdp=None)
    assert a["routes"][0]["tier"] == b["routes"][0]["tier"] == 2


def test_a_route_that_does_not_rest_on_the_scan_keeps_its_tier():
    """CTD and LINCS routes do not come from the motif scan, so a noisy
    scan says nothing about them."""
    lit = {"gene": "G", "reported_links": {
        "wanted_direction": [{"chemical": "Phenobarbital"}]}}
    noisy = cat.build({"literature": lit}, goal="up", scan_fdp=0.9)
    assert noisy["routes"], "the CTD route should still be built"
    assert all(r["tier_before_scan_penalty"] is None
               for r in noisy["routes"] if r["source"] == "ctd")


# ===================================================================== #
# C5. A scenario for compound heterozygotes, checked against the data.  #
# ===================================================================== #
def test_the_biallelic_scenario_does_not_promise_a_wild_type_allele():
    a = mag.ask_for("biallelic_hypomorph", "BUB1B",
                    broken_allele_makes="truncated")
    assert a.goal == "up"
    text = " ".join(a.notes).lower()
    assert "impaired" in text
    assert "truncated" in text, "the amplified broken product must be said"
    assert "stabilis" in text or "readthrough" in text


def test_the_biallelic_scenario_still_agrees_with_an_up_direction():
    assert mag.goal_conflict("biallelic_hypomorph", "up") is None
    assert mag.goal_conflict("biallelic_hypomorph", "down") is not None


def test_a_gene_that_tolerates_losing_one_copy_contradicts_haploinsufficiency():
    """BUB1B's pLI is about zero because MVA is recessive and carriers are
    healthy. The run declared haploinsufficiency and nothing noticed."""
    said = mag.constraint_conflict("haploinsufficiency", 2.6e-14, 0.707)
    assert said and "biallelic_hypomorph" in said
    assert mag.constraint_conflict("haploinsufficiency", 0.99) is None
    assert mag.constraint_conflict("hypomorph", 2.6e-14) is None
    assert mag.constraint_conflict("haploinsufficiency", None) is None


# ===================================================================== #
# C6. Not an agent that causes the defect being treated.                #
# ===================================================================== #
def test_dna_damaging_agents_are_recognised():
    for d in ("Azacitidine", "Decitabine", "Doxorubicin", "Methotrexate",
              "Cyclophosphamide", "Paclitaxel"):
        assert sy.is_genotoxic(d), d
    for d in ("Bexarotene", "Ursodiol", "Levothyroxine"):
        assert not sy.is_genotoxic(d), d


def test_a_genotoxic_route_sinks_only_when_the_patient_cannot_afford_one():
    blk = {"drug_grades": [_grade("Azacitidine")]}
    def graded(out):
        # By source, because the chromatin collector also offers
        # azacitidine and the sort moves a demoted row to the bottom.
        return next(r for r in out["routes"] if r["source"] == "drug_scoring")

    off = graded(cat.build(dict(blk, drug_grades=[_grade("Azacitidine")]),
                           goal="up"))
    out_on = cat.build(dict(blk, drug_grades=[_grade("Azacitidine")]),
                       goal="up", avoid_genotoxic=True)
    on = graded(out_on)
    assert off["tier"] == 2 and off["genotoxic"] is False
    assert on["tier"] == cat.GENOTOXIC_TIER
    assert on["genotoxic"] is True
    assert on["score"] < off["score"]
    assert on["tier"] == max(r["tier"] for r in out_on["routes"]), \
        "a genotoxic route must not sit above anything"
    # Kept, not deleted: a silent removal is not an answer.
    assert len(out_on["routes"]) == len(
        cat.build(dict(blk, drug_grades=[_grade("Azacitidine")]),
                  goal="up")["routes"])


# ===================================================================== #
# C8. When the problem is the protein, say so above the drug tables.    #
# ===================================================================== #
def test_the_protein_banner_fires_on_what_was_declared():
    base = {"query": {"gene": "BUB1B"}, "results": [{}]}
    assert rep._protein_level_section(base) == [], "nothing declared"

    for q in ({"scenario": "biallelic_hypomorph"},
              {"scenario": "hypomorph"},
              {"has_ptc": True},
              {"broken_allele_makes": "truncated"}):
        out = "\n".join(rep._protein_level_section(
            {"query": {"gene": "BUB1B", **q}, "results": [{}]}))
        assert "not the amount" in out, q
        assert "more of the same impairment" in out, q


def test_the_protein_banner_points_at_the_sections_that_ran():
    res = {"query": {"gene": "BUB1B", "has_ptc": True},
           "results": [{"stabilisation": {"FOXM1": {}},
                        "beyond_transcription": {
                            "post_transcriptional": {"routes": [1, 2]}}}]}
    out = "\n".join(rep._protein_level_section(res))
    assert "Stabilisation" in out and "2 route(s)" in out
    assert "Readthrough" in out, "a premature stop names readthrough"


def test_the_banner_says_so_when_no_protein_route_ran():
    out = "\n".join(rep._protein_level_section(
        {"query": {"gene": "G", "has_ptc": True}, "results": [{}]}))
    assert "None ran" in out and "--design-modality" in out


# ===================================================================== #
# Everything above has to survive a real run.                           #
# ===================================================================== #
def test_a_demo_run_carries_the_new_blocks(tmp_path):
    from promoter_tf_drugs import pipeline as pl

    res = pl.run("DEMOGENE", demo=True, scenario="biallelic_hypomorph",
                 broken_allele_makes="truncated")
    assert res["query"]["scenario"] == "biallelic_hypomorph"
    assert "2x from the weak allele" in (res["requirement"] or {})["requirement"]
    blk = res["results"][0]
    assert "proliferation" in blk
    assert "regulators_without_drugs" in blk
    report = pl.markdown_report(res)
    assert "The problem may be the protein, not the amount" in report
    # And it still round-trips to disk.
    written = pl.write_outputs(res, tmp_path)
    assert any(p.name == "results.json" for p in written)
    json.loads((tmp_path / "results.json").read_text(encoding="utf-8"))
