"""Starting from the factors known to regulate the gene.

The route reverses the usual search: take the factors the literature and the
knockout data record as regulating the gene, then ask for their drugs. What
matters is the sign arithmetic (which way a drug must move a factor), which
factors are skipped, and that each source is labelled for what it is.
"""

from __future__ import annotations

from promoter_tf_drugs import catalogue as cat
from promoter_tf_drugs import knockout as ko
from promoter_tf_drugs import regulators as reg
from promoter_tf_drugs.drugs import Drug


def gmt(tmp_path, lines, name="ko.gmt"):
    p = tmp_path / name
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return [str(p)]


def row(term, *genes):
    return term + "\t\t" + "\t".join(genes)


def drug(name, symbol, direction, **kw):
    return Drug(name=name, chembl_id=None, target_symbol=symbol,
                action_type="X", mechanism=None, approved=True, max_phase=4,
                first_approval=2000, withdrawn=False, oral=None,
                direction=direction, source="ChEMBL (approved-drug index)",
                **kw)


# ------------------------------------------------- knockout reverse lookup
def test_a_factor_whose_loss_lowers_the_gene_activates_it(tmp_path):
    p = gmt(tmp_path, [
        row("AAA Down", "G1", "G2"),            # loss, gene falls: activates
        row("BBB Up", "G1"),                    # loss, gene rises: represses
        row("CCC OE HUMAN GSE1 X GENE 1 UP", "G1"),     # gain, rises: activates
        row("DDD KO MOUSE GSE2 X GENE 2 DOWN", "G9"),   # gene not listed
    ])
    out = ko.regulators_of("g1", p)
    assert out["AAA"] == {"activates": 1, "represses": 0, "sign": 1}
    assert out["BBB"] == {"activates": 0, "represses": 1, "sign": -1}
    assert out["CCC"]["sign"] == 1
    assert "DDD" not in out


def test_experiments_are_counted_not_list_memberships(tmp_path):
    p = gmt(tmp_path, [
        row("AAA KO HUMAN GSE1 X GENE 1 DOWN", "G1"),
        row("AAA KO HUMAN GSE1 X GENE 1 UP", "G2"),      # same experiment
        row("AAA KO HUMAN GSE2 X GENE 2 DOWN", "G1"),
        row("AAA KO HUMAN GSE3 X GENE 3 UP", "G1"),
    ])
    assert ko.regulators_of("G1", p)["AAA"] == {
        "activates": 2, "represses": 1, "sign": 1}


def test_a_tie_has_no_sign_and_a_mutation_is_not_read(tmp_path):
    p = gmt(tmp_path, [
        row("AAA KO HUMAN GSE1 X GENE 1 DOWN", "G1"),
        row("AAA KO HUMAN GSE2 X GENE 2 UP", "G1"),
        row("BBB MUTATION HUMAN GSE3 X GENE 3 DOWN", "G1"),
    ])
    out = ko.regulators_of("G1", p)
    assert out["AAA"]["sign"] == 0
    assert "BBB" not in out


def test_a_file_with_windows_line_endings_still_finds_the_last_gene(tmp_path):
    path = tmp_path / "crlf.gmt"
    path.write_bytes((row("AAA Down", "G2", "G1") + "\r\n").encode())
    assert ko.regulators_of("G1", [str(path)])["AAA"]["sign"] == 1


# --------------------------------------------------------------- candidates
SIGNED = {"ACT": {"G1": 1, "G3": 1}, "REP": {"G1": -1}, "NOSIGN": {"G1": 0},
          "OTHER": {"G2": 1}, "KOONLY": {"G7": 1}, "ALSO": {"G8": 1}}


def test_collectri_gives_signed_regulators_and_leaves_out_unsigned_edges():
    assert reg.collectri_regulators("g1", SIGNED) == {"ACT": 1, "REP": -1}
    assert reg.collectri_regulators("G1", None) == {}


def test_knockout_alone_needs_two_agreeing_experiments_and_a_known_tf(tmp_path):
    p = gmt(tmp_path, [
        row("KOONLY KO HUMAN GSE1 X GENE 1 DOWN", "G1"),
        row("KOONLY KO HUMAN GSE2 X GENE 2 DOWN", "G1"),     # two: kept
        row("ALSO KO HUMAN GSE3 X GENE 3 DOWN", "G1"),       # one: dropped
        row("NOTATF KO HUMAN GSE4 X GENE 4 DOWN", "G1"),
        row("NOTATF KO HUMAN GSE5 X GENE 5 DOWN", "G1"),     # two, not a TF
    ])
    got = {r.symbol: r for r in reg.candidates(
        "G1", SIGNED, knockout_paths=p, use_collectri=False)}
    assert set(got) == {"KOONLY"}
    assert got["KOONLY"].sources == [reg.KNOCKOUT] and got["KOONLY"].sign == 1
    anyf = {r.symbol for r in reg.candidates(
        "G1", SIGNED, knockout_paths=p, use_collectri=False, tf_only=False)}
    assert anyf == {"KOONLY", "NOTATF"}


def test_the_two_sources_are_joined_and_a_disagreement_is_recorded(tmp_path):
    p = gmt(tmp_path, [
        row("ACT KO HUMAN GSE1 X GENE 1 DOWN", "G1"),        # agrees: +1
        row("REP KO HUMAN GSE2 X GENE 2 DOWN", "G1"),        # says +1, table -1
    ])
    got = {r.symbol: r for r in reg.candidates("G1", SIGNED,
                                               knockout_paths=p)}
    assert got["ACT"].sources == [reg.COLLECTRI, reg.KNOCKOUT]
    assert got["ACT"].n_responses == 1 and not got["ACT"].disagree
    assert got["REP"].disagree and got["REP"].sign == -1, (
        "the curated sign is kept and the disagreement is shown")


def test_factors_that_already_have_a_route_are_skipped(tmp_path):
    got = reg.candidates("G1", SIGNED, knockout_paths=[], skip=["act"],
                         use_knockout=False)
    assert [r.symbol for r in got] == ["REP"]


# ---------------------------------------------------------------- find_drugs
def test_the_drug_must_move_the_factor_the_way_the_sign_and_goal_need():
    index = {"ACT": [drug("Agonist", "ACT", "up"),
                     drug("Blocker", "ACT", "down"),
                     drug("Vague", "ACT", "unclear")],
             "REP": [drug("Agonist2", "REP", "up"),
                     drug("Blocker2", "REP", "down")]}
    regs = [reg.Regulator("ACT", 1, ["CollecTRI"]),
            reg.Regulator("REP", -1, ["CollecTRI"])]
    up = reg.find_drugs(regs, "up", index)
    # raise the gene: switch an activator on, switch a repressor off
    assert [d["name"] for d in up["ACT"]["drugs"]] == ["Agonist", "Vague"]
    assert [d["name"] for d in up["REP"]["drugs"]] == ["Blocker2"]
    down = reg.find_drugs(regs, "down", index)
    assert [d["name"] for d in down["ACT"]["drugs"]] == ["Blocker", "Vague"]
    assert [d["name"] for d in down["REP"]["drugs"]] == ["Agonist2"]
    assert up["ACT"]["need"] == "up" and down["ACT"]["need"] == "down"


def test_a_regulator_with_no_matching_drug_is_absent_and_the_list_is_capped():
    index = {"ACT": [drug(f"D{i}", "ACT", "up")
                     for i in range(reg.MAX_DRUGS_PER_FACTOR + 5)]}
    regs = [reg.Regulator("ACT", 1, ["CollecTRI"]),
            reg.Regulator("NODRUG", 1, ["CollecTRI"])]
    got = reg.find_drugs(regs, "up", index)
    assert list(got) == ["ACT"]
    assert len(got["ACT"]["drugs"]) == reg.MAX_DRUGS_PER_FACTOR
    assert got["ACT"]["n_drugs"] == reg.MAX_DRUGS_PER_FACTOR + 5


# ------------------------------------------------------------------- routes
def _record(**kw):
    base = {"sign": 1, "need": "up", "sources": ["CollecTRI"],
            "n_responses": 0, "n_against": 0, "disagree": False,
            "drugs": [drug("Agonist", "ACT", "up").to_dict()]}
    base.update(kw)
    return {"ACT": base}


def test_routes_are_tier_four_and_say_what_they_do():
    (r,) = cat.from_regulators(_record())
    assert r.category == "TF: known regulator" == reg.CATEGORY
    assert r.tier == 4 and r.directness == 1 and r.source == reg.SOURCE
    assert r.mechanism.startswith("activate ACT")
    assert any("raises the gene according to CollecTRI" in n for n in r.notes)
    (r2,) = cat.from_regulators(_record(sign=-1, need="down"))
    assert r2.mechanism.startswith("inhibit ACT")
    assert any("lowers the gene" in n for n in r2.notes)


def test_a_curated_edge_with_knockout_support_scores_above_either_alone():
    def strength(**kw):
        (r,) = cat.from_regulators(_record(**kw))
        return r.axes.strength
    both = strength(sources=["CollecTRI", "knockout data"], n_responses=3)
    assert both > strength() > strength(sources=["knockout data"])
    assert strength(sources=["CollecTRI", "knockout data"],
                    disagree=True) == strength()


def test_unclear_family_and_disagreement_are_stated_in_the_notes():
    fam = drug("Famdrug", "ACT", "unclear", target_type="PROTEIN FAMILY")
    (r,) = cat.from_regulators(_record(
        drugs=[fam.to_dict()], sources=["CollecTRI", "knockout data"],
        disagree=True, n_responses=1, n_against=2))
    text = " ".join(r.notes)
    assert "not recorded as raising or lowering" in text
    assert "protein family" in text
    assert "opposite sign" in text
    assert "1 knockout experiment(s) agree and 2 disagree" in text


def test_the_catalogue_builds_these_routes_from_the_result_block():
    out = cat.build({"regulators": _record(), "gene_symbol": "G1"})
    cats = {r["category"] for r in out["routes"]}
    assert "TF: known regulator" in cats
    assert "TF: known regulator" not in out["categories_absent"]
    assert "TF: known regulator" in cat.CATEGORY_GUIDE


def test_a_run_without_the_flag_has_no_regulator_routes():
    out = cat.build({"gene_symbol": "G1"})
    assert "TF: known regulator" in out["categories_absent"]


def test_the_flag_reaches_run_from_both_command_lines():
    from promoter_tf_drugs import cli, panel
    for mod, base in ((cli, ["--gene", "CCR5"]), (panel, ["--genes", "CCR5"])):
        parser = mod.build_parser()
        assert cli.run_kwargs_from_args(
            parser.parse_args(base))["regulator_route"] is False
        assert cli.run_kwargs_from_args(parser.parse_args(
            base + ["--regulator-route"]))["regulator_route"] is True


# ------------------------------------------------- the evaluation script
def _gain_script():
    import importlib.util
    import sys
    from pathlib import Path
    path = Path(__file__).resolve().parents[1] / "scripts" / "regulator_gain.py"
    spec = importlib.util.spec_from_file_location("_regulator_gain", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_the_bootstrap_resamples_genes_and_gives_an_ordered_interval():
    rg = _gain_script()
    per_run = [(3, 1), (0, 2), (2, 2), (4, 0), (1, 3), (2, 1)]
    chance = [(2, 2)] * 6
    full, lo, hi, glo, ghi = rg.bootstrap_share(per_run, chance)
    assert abs(full - 12 / 21) < 1e-9
    assert 0 <= lo <= full <= hi <= 1
    assert glo <= ghi
    again = rg.bootstrap_share(per_run, chance)
    assert again == (full, lo, hi, glo, ghi), "fixed seed, same answer"


def test_a_gene_with_no_verdicts_cannot_make_the_interval_up():
    rg = _gain_script()
    full, lo, hi, glo, ghi = rg.bootstrap_share([(0, 0)] * 4, [(1, 1)] * 4)
    assert full is None and lo is None and ghi is None


def test_the_evaluation_variants_switch_sources_on_and_off():
    rg = _gain_script()
    signed = {"ACT": {"G1": 1}}
    index = {"ACT": [drug("Agonist", "ACT", "up")]}
    # no knockout files here, so only CollecTRI can supply a regulator
    kw = dict(gene="G1", goal="up", skip=set(), signed=signed, index=index,
              ko_paths=[])
    assert [r["drug"] for r in rg.regulator_routes("collectri", **kw)] == [
        "Agonist"]
    assert rg.regulator_routes("baseline", **kw) == []
    assert rg.regulator_routes("knockout, TFs", **kw) == []
