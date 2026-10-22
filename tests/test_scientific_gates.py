"""The gates added after the BUB1B review.

Each one exists because the BUB1B run was confidently wrong in a specific
way, and each test names that way. They are offline: the one that needs a
LINCS signature file writes a small one.
"""

from __future__ import annotations

import json

from generoutes import catalogue as cat
from generoutes import abundance as ab
from generoutes import crosscheck as xc
from generoutes import fame as fm
from generoutes import magnitude as mag
from generoutes import proliferation as pr
from generoutes import regulators as reg
from generoutes import report as rep
from generoutes import systemic as sy


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


def test_the_direct_regulators_are_not_buried_by_the_two_step_ones():
    """With --regulator-steps 2 the two-step factors outnumbered the direct
    ones twenty to one. Sorted by symbol and capped, the rows shown were
    ABL1 to CREB1 and FOXM1, MYC, TP53 and ZNF143 appeared nowhere - the
    four this section exists to name."""
    direct = [{"symbol": s, "sign": 1, "need": "up",
               "sources": ["CollecTRI"], "via": None}
              for s in ("FOXM1", "MYC", "TP53", "ZNF143")]
    # Alphabetically ahead of every direct one, as the real data was.
    two_step = [{"symbol": f"A{i:03d}", "sign": 1, "need": "up",
                 "sources": ["CollecTRI"], "via": "TP53"}
                for i in range(200)]
    out = "\n".join(rep._undrugged_regulators(
        {"regulators_without_drugs": two_step + direct}))

    for sym in ("FOXM1", "MYC", "TP53", "ZNF143"):
        assert f"| {sym} |" in out, sym
    assert "A000" not in out, "a two-step factor must not take a row"
    assert "**200**" in out, "they are counted instead"


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
    from generoutes import pipeline as pl

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
# ===================================================================== #
# C9. How much of the order is drug fame?                               #
# ===================================================================== #
def _fame_table(tmp_path, counts, n_genes):
    path = tmp_path / "ctd_fame.tsv"
    fm.write_table(path, counts, n_genes=n_genes)
    return fm.load(path)


def test_a_compound_recorded_for_every_gene_is_called_promiscuous(tmp_path):
    """Tetrachlorodibenzodioxin, benzo(a)pyrene and bisphenol A are in the
    CTD record of every gene in the saved runs here. A chemical recorded
    for all of them is telling you about the literature, not your gene."""
    table = _fame_table(tmp_path, {
        "Tetrachlorodibenzodioxin": 11, "Bisphenol A": 11,
        "Valproic Acid": 8, "Obeticholic Acid": 1, "Ursodiol": 1,
        **{f"Rare {i}": 1 for i in range(30)},
    }, 11)
    got = fm.assess(["Tetrachlorodibenzodioxin", "Valproic acid",
                     "Obeticholic Acid"], table)
    assert got["tetrachlorodibenzodioxin"]["promiscuous"] is True
    assert got["valproicacid"]["promiscuous"] is True
    assert got["obeticholicacid"]["promiscuous"] is False
    # Salt and spelling differences must meet.
    assert got["valproicacid"]["n_genes"] == 8
    assert got["valproicacid"]["of"] == 11
    assert got["valproicacid"]["name"] == "Valproic acid", "the caller's spelling"


def test_a_compound_the_table_does_not_hold_gets_no_verdict(tmp_path):
    table = _fame_table(tmp_path, {"A": 9, **{f"R{i}": 1 for i in range(30)}},
                        11)
    assert fm.assess(["Nobody"], table) == {}, "absent is not a pass"


def test_a_table_over_too_few_genes_is_refused(tmp_path):
    """A count over three genes says nothing about how well known a drug is."""
    path = tmp_path / "thin.tsv"
    fm.write_table(path, {"A": 3, "B": 1}, n_genes=3)
    assert fm.load(path) is None


def test_no_table_at_all_is_reported_rather_than_assumed(tmp_path):
    assert fm.load(tmp_path / "missing.tsv") is None
    assert fm.assess(["A"], None) == {}
    said = fm.summarise({}, None)
    assert "No fame table" in said and "build_fame_table" in said


def test_a_promiscuous_compound_loses_a_little_score_and_keeps_its_tier():
    """Slight on purpose: being well studied is not evidence of being
    wrong. The claim is that such a row says little about THIS gene."""
    rec = {"name": "Valproic acid", "n_genes": 8, "of": 11, "share": 0.73,
           "cut": 7, "percentile": 90, "promiscuous": True}
    blk = {"drug_grades": [_grade("Valproic acid")]}
    plain = cat.build(dict(blk, drug_grades=[_grade("Valproic acid")]),
                      goal="up")
    cut = cat.build(dict(blk, drug_grades=[_grade("Valproic acid")]),
                    goal="up", fame={"valproicacid": rec})

    def graded(out):
        return next(r for r in out["routes"] if r["source"] == "drug_scoring")

    a, b = graded(plain), graded(cut)
    assert b["tier"] == a["tier"], "fame never changes a tier"
    assert b["score"] == round(a["score"] * cat.FAME_SCORE_FACTOR, 1)
    assert b["fame_penalty"] == cat.FAME_SCORE_FACTOR
    assert a["fame_penalty"] is None
    assert any("8 of the 11 genes" in n for n in b["notes"])


def test_the_penalty_is_slight_rather_than_a_demotion():
    """A heavier hand would be a verdict on the drug. It is a caution
    about what the row can tell you."""
    assert 0.75 < cat.FAME_SCORE_FACTOR < 1.0
    assert cat.FAME_SCORE_FACTOR > cat.BROAD_SCORE_FACTOR


def test_the_report_names_the_compound_the_way_the_caller_did(tmp_path):
    table = _fame_table(tmp_path, {
        "Valproic Acid": 9, **{f"R{i}": 1 for i in range(30)}}, 11)
    got = fm.assess(["Valproic acid"], table)
    out = "\n".join(rep._fame_section({"results": [{"fame": got}]}))
    assert "Valproic acid" in out and "valproicacid |" not in out
    assert "9 of 11" in out


def test_the_builder_counts_a_gene_once_however_many_runs_it_has(tmp_path):
    """Eight reruns of BUB1B must not make a chemical eight times famous."""
    import importlib.util
    import json
    from pathlib import Path

    src = Path(__file__).resolve().parents[1] / "scripts" / "build_fame_table.py"
    spec = importlib.util.spec_from_file_location("_bft", src)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    chemicals_by_gene = mod.chemicals_by_gene

    def run(folder, gene, chems):
        d = tmp_path / folder
        d.mkdir()
        (d / "results.json").write_text(json.dumps({
            "query": {"gene": gene},
            "results": [{"improvements": {"reported_links": {
                "wanted_direction": [{"chemical": c} for c in chems]}}}],
        }), encoding="utf-8")

    run("BUB1B_up", "BUB1B", ["Estradiol", "Valproic Acid"])
    run("BUB1B_up_v2", "BUB1B", ["Estradiol"])
    run("CCR5_down", "CCR5", ["Estradiol"])
    by_gene = chemicals_by_gene(tmp_path)
    assert set(by_gene) == {"BUB1B", "CCR5"}
    counts = {}
    for names in by_gene.values():
        for n in names:
            counts[n] = counts.get(n, 0) + 1
    assert counts["Estradiol"] == 2, "two genes, not three runs"
    assert counts["Valproic Acid"] == 1
# ===================================================================== #
# C10. The abundance scan asks about the regulators nothing binds.      #
# ===================================================================== #
def test_a_regulators_sign_becomes_the_role_the_abundance_scan_needs():
    """The pipeline turns a recorded regulator into a factor the abundance
    scan can use, and the whole value of that is the sign. A factor that
    RAISES the gene has to go up to raise it; one that lowers it has to go
    down. Inverting this would recommend the exact opposite compound, and
    nothing downstream would notice."""
    raises = "activator"      # sign +1
    lowers = "repressor"      # sign -1
    assert ab.wanted_direction("up", raises) == "up"
    assert ab.wanted_direction("up", lowers) == "down"
    assert ab.wanted_direction("down", raises) == "down"
    assert ab.wanted_direction("down", lowers) == "up"


def test_the_abundance_scan_reads_a_regulator_shaped_row(tmp_path):
    """The rows the pipeline appends carry only `symbol` and `role`, not
    the `TFResult.to_dict` shape the scan was written against."""
    up = tmp_path / "up.gmt"
    down = tmp_path / "down.gmt"
    # Three signatures each, because one profiling is the weakest thing a
    # corpus can say and the scan asks for more than that.
    up.write_text("".join(
        f"CPC001 HA1E 24H-raiser-{i}.0\t\tFOXM1\tALB\n" for i in range(4)),
        encoding="utf-8")
    down.write_text("".join(
        f"CPC001 HA1E 24H-lowerer-{i}.0\t\tSRF\tALB\n" for i in range(4)),
        encoding="utf-8")

    rows = [{"symbol": "FOXM1", "role": "activator"},
            {"symbol": "SRF", "role": "repressor"}]
    hits = ab.find(rows, goal="up", up_path=up, down_path=down)
    got = {(h.tf, h.drug) for h in hits}
    assert ("FOXM1", "raiser") in got, "an activator wanted up"
    assert ("SRF", "lowerer") in got, "a repressor wanted down"


def test_the_pipeline_asks_the_undrugged_regulators_too():
    """They are the factors with no ligand, which is what this step is
    for - and they come from CollecTRI, not from the motif scan, so the
    `ranked[:N]` slice never reached them."""
    import inspect
    from generoutes import pipeline as pl

    src = inspect.getsource(pl.run)
    assert "for _reg in regulators_undrugged:" in src
    assert '"activator" if _reg.get("sign") == 1' in src
    assert '_reg.get("via")' in src, "two-step regulators stay out"
# ===================================================================== #
# C11. Not a route that raises an oncogene in a cancer-prone patient.   #
# ===================================================================== #
def _ab_hit(drug, tf, n=20, space="landmark"):
    return {"drug": drug, "tf": tf, "tf_role": "activator",
            "moved_factor": "up", "n_signatures": n, "feature_space": space}


def test_raising_an_oncogene_is_a_risk_and_raising_a_suppressor_is_not():
    """The abundance scan offered to raise MYC, FOXM1 and ETS2 for a child
    who is cancer-predisposed, and nothing said so. Raising TP53 is the
    control: it is a tumour suppressor, so more of it is not the risk."""
    assert sy.cancer_risk_direction("MYC", "up")
    assert sy.cancer_risk_direction("FOXM1", "up")
    assert sy.cancer_risk_direction("TP53", "down")
    assert sy.cancer_risk_direction("TP53", "up") is None
    assert sy.cancer_risk_direction("MYC", "down") is None, "lowering it is not the risk"
    assert sy.cancer_risk_direction("OSR2", "up") is None, "not on either list"


def test_an_unrecorded_direction_makes_no_claim():
    """Guessing the direction would invent the thing the caution rests on."""
    assert sy.cancer_risk_direction("MYC", None) is None
    assert sy.cancer_risk_direction("MYC", "sideways") is None
    assert sy.cancer_risk_direction("", "up") is None


def test_the_gate_sinks_the_oncogene_route_and_leaves_the_suppressor():
    blk = {"factor_abundance": [_ab_hit("torin-1", "MYC"),
                                _ab_hit("wortmannin", "TP53")],
           "gene_symbol": "BUB1B",
           "tf_roles": {"MYC": "activator", "TP53": "activator"}}
    off = cat.build(dict(blk), goal="up")
    on = cat.build(dict(blk), goal="up", avoid_genotoxic=True)

    def by_drug(out, name):
        return next(r for r in out["routes"] if r["drug"] == name)

    assert by_drug(off, "torin-1")["cancer_risk"] is None
    sunk = by_drug(on, "torin-1")
    assert sunk["tier"] == cat.GENOTOXIC_TIER
    assert "proto-oncogene" in sunk["cancer_risk"]
    assert sunk["score"] < by_drug(off, "torin-1")["score"]
    # The suppressor route is untouched: more TP53 is not the hazard.
    assert by_drug(on, "wortmannin")["tier"] == by_drug(
        off, "wortmannin")["tier"]
    assert by_drug(on, "wortmannin")["cancer_risk"] is None
    # Kept, not deleted, like every other demotion here.
    assert len(on["routes"]) == len(off["routes"])


def test_the_gate_is_silent_without_the_flag():
    blk = {"factor_abundance": [_ab_hit("torin-1", "MYC")],
           "gene_symbol": "BUB1B", "tf_roles": {"MYC": "activator"}}
    out = cat.build(blk, goal="up")
    assert all(r["cancer_risk"] is None for r in out["routes"])


def test_a_regulators_role_reaches_the_catalogue():
    """MYC and FOXM1 come from CollecTRI, not the motif scan, so without
    their sign in tf_roles the check can make no claim about them."""
    import inspect
    from generoutes import pipeline as pl

    src = inspect.getsource(pl.run)
    assert '"tf_roles": {' in src
    assert 'for r in regulators_undrugged if not r.get("via")' in src
# ===================================================================== #
# C12. Screen the endpoint the patient needs, not a reporter.          #
# ===================================================================== #
class _Prom:
    gene_symbol = "BUB1B"
    sequence = "ACGT" * 200
    chrom, start, end = "15", 1000, 2000
    upstream, downstream = 1000, 200
    transcript_id = "ENST1"

    @staticmethod
    def offset_to_tss(i):
        return i - 1000


def test_a_declared_protein_problem_screens_the_protein():
    """A promoter reporter cannot see a stabiliser, a splice modifier or
    a readthrough agent, and those are the mechanisms in play when the
    protein is the problem."""
    from generoutes import screens as sc

    plan = sc.design(_Prom(), ["FOXM1"], [], [],
                     protein_endpoint=True,
                     why="the only working copy is impaired").to_dict()
    assert plan["reporter"] is None, "no reporter when it is the wrong assay"
    assert plan["protein"], "the protein plan is there instead"
    assert "protein" in plan["endpoint"]
    assert "impaired" in plan["protein"]["why"]
    # The CRISPRa arm stays: it raises the endogenous gene either way.
    assert plan["crispra"]


def test_the_default_is_still_the_reporter():
    from generoutes import screens as sc

    plan = sc.design(_Prom(), ["FOXM1"], [], []).to_dict()
    assert plan["reporter"] and plan["protein"] is None
    assert plan["endpoint"] == "promoter reporter"


def test_the_screen_controls_for_the_cell_cycle_confound():
    """For a gene whose protein tracks the cycle, a compound that arrests
    cells raises the signal without touching the gene. Every other
    finding in this package says that is the error to design against."""
    from generoutes import screens as sc

    pp = sc.protein_abundance_plan("BUB1B")
    text = " ".join(pp["normalisation"] + pp["controls"] + pp["cells"]).lower()
    assert "phospho-histone h3" in text or "dna content" in text
    assert "cell-cycle" in text or "cell cycle" in text
    assert "not pcna" in text, "the housekeeping trap is named"
    assert "cytotoxicity" in text, "most primary hits are simply toxic"
    assert any("sub-confluent" in c for c in pp["cells"])


def test_the_screen_is_against_approved_drugs_only():
    from generoutes import screens as sc

    lib = sc.protein_abundance_plan("BUB1B")["library"]
    assert "approved" in lib["what"]
    assert any("Prestwick" in s or "Repurposing" in s for s in lib["sources"])


def test_a_negative_screen_is_reported_as_informative():
    from generoutes import screens as sc

    said = sc.protein_expectations()["what_a_negative_result_means"]
    assert "informative" in said
    assert "repurposing is unlikely" in said


def test_the_endpoint_follows_what_was_declared(tmp_path):
    from generoutes import pipeline as pl

    protein = pl.run("DEMOGENE", demo=True, design_screen=True,
                     scenario="biallelic_hypomorph",
                     broken_allele_makes="truncated")
    reporter = pl.run("DEMOGENE", demo=True, design_screen=True)
    assert protein["results"][0]["screen_plan"]["protein"]
    assert protein["results"][0]["screen_plan"]["reporter"] is None
    assert reporter["results"][0]["screen_plan"]["reporter"]
    assert reporter["results"][0]["screen_plan"]["protein"] is None
    out = pl.markdown_report(protein)
    assert "Primary screen:" in out and "Reporter construct" not in out
# ===================================================================== #
# C13. Drug to PROTEIN abundance, where the rest reads transcripts.     #
# ===================================================================== #
def _prot_table(tmp_path, rows):
    p = tmp_path / "prot.tsv"
    p.write_text("drug\tprotein\tlog2fc\tcell_line\n" + "".join(
        f"{d}\t{g}\t{fc}\t{c}\n" for d, g, fc, c in rows),
        encoding="utf-8")
    return p


def test_one_observation_is_not_a_measurement(tmp_path):
    """Proteomics is noisier per run than RNA-seq, and a single 1.5x
    shift in one cell line is a lead to repeat, not a finding."""
    from generoutes import proteomics as pm

    path = _prot_table(tmp_path, [
        ("MG132", "BUB1B", 1.3, "HEK293"),
        ("MG132", "BUB1B", 1.1, "Jurkat"),
        ("onceonly", "BUB1B", 2.0, "Jurkat"),
        ("tiny", "BUB1B", 0.1, "Jurkat"),
        ("tiny", "BUB1B", 0.2, "A549"),
    ])
    got = {h.drug for h in pm.load("BUB1B", path=path, want="up")}
    assert got == {"MG132"}, "one observation and a small shift both go"


def test_the_direction_asked_for_is_the_direction_returned(tmp_path):
    from generoutes import proteomics as pm

    path = _prot_table(tmp_path, [
        ("raiser", "BUB1B", 1.2, "A"), ("raiser", "BUB1B", 1.4, "B"),
        ("lowerer", "BUB1B", -1.2, "A"), ("lowerer", "BUB1B", -1.4, "B"),
    ])
    assert [h.drug for h in pm.load("BUB1B", path=path, want="up")] == ["raiser"]
    assert [h.drug for h in pm.load("BUB1B", path=path, want="down")] == ["lowerer"]


def test_no_table_is_reported_and_names_the_right_datasets(tmp_path):
    """LINCS P100 is ~90 phosphopeptides and CPTAC has no drug arm, so
    naming them here would send somebody after the wrong file."""
    from generoutes import proteomics as pm

    assert pm.load("BUB1B", path=tmp_path / "nope.tsv") == []
    said = pm.summarise([], "BUB1B", path=tmp_path / "nope.tsv")
    assert "decryptE" in pm.SOURCE and "ProTargetMiner" in pm.SOURCE
    assert "do not" in pm.SOURCE and "P100" in pm.SOURCE
    assert "transcript" in said


def test_a_protein_measurement_outranks_a_transcript_one(tmp_path):
    """It is measured on the thing the patient needs, with nothing
    inferred between compound and product."""
    hits = [{"drug": "MG132", "protein": "BUB1B", "direction": "up",
             "median_log2fc": 1.2, "n_observations": 3,
             "cell_lines": ["Jurkat"], "sources": [], "notes": []}]
    out = cat.build({"protein_abundance": hits, "gene_symbol": "BUB1B"},
                    goal="up")
    r = next(x for x in out["routes"] if x["source"] == "proteomics")
    assert r["tier"] == 2 and r["directness"] == 0
    assert r["category"] == "protein abundance change"
    assert any("says nothing about HOW" in n for n in r["notes"])


# ===================================================================== #
# C14. Six steps to a protein, and which were examined.                 #
# ===================================================================== #
def test_a_step_that_is_only_named_is_not_counted_as_covered():
    """A full table of transcriptional routes is not evidence that the
    other five steps were considered, and nothing else would say so."""
    from generoutes import coverage as cov

    res = {"results": [{"all_routes": {"routes": [
        {"category": "TF: orthosteric ligand"},
        {"category": "chromatin-level"}]},
        "beyond_transcription": {
            "post_transcriptional": {"routes": [{"mechanism": "x"}]},
            "product_stability": {"routes": [{"mechanism": "y"}]}}}],
        "brake_exons": {"candidates": []}}
    rows = {r["step"]: r for r in cov.audit(res)}
    assert rows["transcription"]["state"] == cov.SEARCHED
    assert rows["transcription"]["found"] == 2
    assert rows["mRNA stability"]["state"] == cov.LISTED
    assert rows["translation"]["state"] == cov.LISTED
    for step in ("mRNA stability", "translation"):
        assert rows[step]["enable"], "a gap has to say what would close it"


def test_a_gene_specific_lookup_counts_as_searched():
    from generoutes import coverage as cov

    res = {"results": [{"all_routes": {"routes": []},
                        "product_screen": {"hits": [{"drug": "a"}]},
                        "protein_abundance": [{"drug": "MG132"}],
                        "beyond_transcription": {}}],
           "brake_exons": {"candidates": [1, 2, 3]}}
    rows = {r["step"]: r for r in cov.audit(res)}
    assert rows["splicing and isoform choice"]["state"] == cov.SEARCHED
    assert rows["splicing and isoform choice"]["found"] == 3
    assert rows["folding"]["state"] == cov.SEARCHED
    assert rows["degradation"]["state"] == cov.SEARCHED
    assert rows["degradation"]["found"] == 1


def test_the_summary_states_the_ratio_plainly():
    from generoutes import coverage as cov

    res = {"results": [{"all_routes": {"routes": [
        {"category": "TF: orthosteric ligand"}] * 300},
        "beyond_transcription": {}}],
        "brake_exons": {"candidates": [1]}}
    said = cov.summarise(cov.audit(res))
    assert "300 transcriptional route(s)" in said
    assert "property of this package, not of the biology" in said


def test_the_audit_reaches_the_report():
    from generoutes import pipeline as pl

    res = pl.run("DEMOGENE", demo=True)
    out = pl.markdown_report(res)
    assert "How much of the problem this run examined" in out
    assert "listed only" in out
    assert "mRNA stability" in out and "translation" in out
# ===================================================================== #
# C15. Translation and decay, read off the transcript's own UTRs.       #
# ===================================================================== #
_UORF = "ATGGCGTTGCCGTAA"          # ATG, three codons, stop
_FILLER = "GCCACCGCGAGCGCGGCACCG"


def test_an_upstream_orf_is_found_with_its_kozak():
    from generoutes import utr

    five = "GGCACGCG" + _UORF + _FILLER
    got = utr.find_uorfs(five)
    assert len(got) == 1
    assert got[0]["codons"] == 4 and got[0]["overlaps_cds"] is False
    assert got[0]["kozak"] in ("strong", "adequate", "weak")


def test_an_orf_with_no_stop_is_marked_as_overlapping_the_cds():
    """It is out-of-frame overlap with the main ORF, which represses
    differently from a discrete uORF, so it is not counted as one."""
    from generoutes import utr

    got = utr.find_uorfs("GGCACG" + "ATG" + "GCGGCGGCGGCGGCG")
    assert got and got[0]["overlaps_cds"] is True


def test_a_feature_is_judged_against_a_shuffle_of_its_own_sequence():
    """AUUUA occurs constantly in an AU-rich 3'UTR by composition alone.
    Reporting that as a finding would be reporting the base content."""
    from generoutes import utr

    # AU-rich but with no clustered element beyond what composition gives.
    plain = "ATATATATTAATATTAATATATTAAT" * 8
    s = utr.scan("G" * 60, plain, shuffles=30)
    pent = s["stability"]["pentamer"]
    assert pent["n_observed"] >= 0
    assert not pent["enriched"], "composition alone must not count"
    assert "expected" in pent and "p_value" in pent


def test_a_real_cluster_does_beat_the_shuffle():
    from generoutes import utr

    three = ("TTATTTATT" + "GCGCGCGCGCGC") * 10
    s = utr.scan("G" * 60, three, shuffles=50)
    assert s["stability"]["nonamer"]["enriched"] is True


def test_a_utr_too_short_to_scan_says_so_rather_than_returning_zero():
    from generoutes import utr

    s = utr.scan("ACGT", "ACGT", shuffles=5)
    assert "error" in s["translation"] and "error" in s["stability"]
    said = utr.summarise(s, "G")
    assert "too short" in said


def test_the_scan_makes_the_two_listed_steps_searched():
    """These were the rows the coverage audit marked 'listed only'."""
    from generoutes import coverage as cov

    base = {"results": [{"all_routes": {"routes": []},
                         "beyond_transcription": {
                             "post_transcriptional": {"routes": [{"m": 1}]}}}],
            "brake_exons": {}}
    before = {r["step"]: r["state"] for r in cov.audit(base)}
    assert before["translation"] == cov.LISTED
    assert before["mRNA stability"] == cov.LISTED

    scanned = dict(base)
    scanned["results"] = [dict(base["results"][0], utr={
        "five_prime_length": 70, "three_prime_length": 206,
        "translation": {"n_observed": 2, "enriched": True},
        "stability": {"pentamer": {"enriched": False},
                      "nonamer": {"enriched": True}}})]
    after = {r["step"]: r for r in cov.audit(scanned)}
    assert after["translation"]["state"] == cov.SEARCHED
    assert after["mRNA stability"]["state"] == cov.SEARCHED
    assert after["mRNA stability"]["found"] == 1
    # And it still says no drug was matched, because none was.
    assert "no drug was matched" in after["mRNA stability"]["gap"]


def test_a_scan_that_found_nothing_is_still_a_search():
    """'No uORF in a 200 bp 5\'UTR' is a result about this transcript."""
    from generoutes import coverage as cov

    res = {"results": [{"all_routes": {"routes": []},
                        "beyond_transcription": {},
                        "utr": {"five_prime_length": 200,
                                "three_prime_length": 300,
                                "translation": {"n_observed": 0,
                                                "enriched": False},
                                "stability": {"pentamer": {"enriched": False},
                                              "nonamer": {"enriched": False}}}}],
           "brake_exons": {}}
    rows = {r["step"]: r for r in cov.audit(res)}
    assert rows["translation"]["state"] == cov.SEARCHED
    assert rows["translation"]["found"] == 0
# ===================================================================== #
# C16. The other half of the degradation balance.                       #
# ===================================================================== #
def test_the_dub_axis_states_the_direction_before_the_drug():
    """To RAISE a protein you need MORE deubiquitination, so you have to
    enhance the DUB. Essentially all DUB chemistry is inhibitors, which
    lower the substrate - a reader who sees "a USP7 inhibitor exists" and
    reaches for it would be moving the protein the wrong way."""
    from generoutes import stabilizers as st

    said = st.summarise_rescue(st.rescue_route("TP53", drug_index={}))
    assert "USP7" in said
    assert said.index("ENHANCING") < said.index("inhibitors"), (
        "the direction has to come before the chemistry")
    assert "LOWERS that protein" in said


def test_an_unrecorded_dub_is_reported_as_a_gap_not_an_absence():
    """Not knowing which DUB rescues a protein says the counter-
    regulation is unexplored, which is where a selective lever would be."""
    from generoutes import stabilizers as st

    rec = st.rescue_route("BUB1B", drug_index={})
    assert rec["established"] is False and rec["dubs"] == []
    said = st.summarise_rescue(rec)
    assert "No deubiquitinase is recorded" in said
    assert "gap in the biology" in said


def test_both_halves_of_the_balance_reach_the_product_block():
    from generoutes import stabilizers as st

    out = st.product_stabilisation("TP53", goal="up", drug_index={})
    assert out["degradation"]["ligase"] == "MDM2"
    assert out["rescue"]["dubs"] == ["USP7", "USP10"]
    # Still only meaningful for raising a protein.
    assert st.product_stabilisation("TP53", goal="down") == {}


def test_the_report_shows_the_ligase_and_the_enzyme_that_reverses_it():
    from generoutes import pipeline as pl

    res = pl.run("DEMOGENE", demo=True)
    out = pl.markdown_report(res)
    # DEMOGENE is not in either table, so the report must say so rather
    # than omit the question.
    assert "deubiquitinase" in out.lower() or "Raising the protein" in out
