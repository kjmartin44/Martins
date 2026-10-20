"""How many motif sites chance gives, and the scan that has to honour its p-value.

Two faults were found by shuffling a real masked promoter. A motif whose best
possible match is more common than the requested p-value was scanned at its
consensus word, which matched at several times the stated rate (254 of 879
matrices in one window). And both strands shared one threshold although the
reverse strand of a matrix has a different score distribution when the bases
are not balanced. The estimate of chance sites then used a formula that
assumed neither happened, and undercounted by half.
"""

from __future__ import annotations

import random
from collections import Counter

import numpy as np
import pytest

from promoter_tf_drugs import motifs as m
from promoter_tf_drugs import pwm


def _pairs(seq):
    return Counter(zip(seq, seq[1:]))


def _motif(rows, mid="MA0001.1", name="TEST"):
    return m.Motif(matrix_id=mid, name=name,
                   counts=np.array(rows, dtype=float))


def _sharp(consensus, n=20, mid="MA0001.1"):
    """A matrix that is very sure of every base of `consensus`."""
    rows = []
    for b in consensus:
        r = [0.2, 0.2, 0.2, 0.2]
        r["ACGT".index(b)] = n
        rows.append(r)
    return _motif(rows, mid)


# ------------------------------------------------------------ the shuffle
def test_a_dinucleotide_shuffle_keeps_every_two_base_count_and_both_ends():
    rng = random.Random(3)
    src = random.Random(1)
    seq = "".join(src.choice("ACGT") for _ in range(400))
    for _ in range(20):
        out = m.dinucleotide_shuffle(seq, rng)
        assert len(out) == len(seq)
        assert _pairs(out) == _pairs(seq)
        assert out[0] == seq[0] and out[-1] == seq[-1]
    assert m.dinucleotide_shuffle(seq, rng) != seq


def test_short_runs_are_returned_as_they_are():
    assert m.dinucleotide_shuffle("ACG", random.Random(0)) == "ACG"


def test_the_shuffle_leaves_the_repeat_mask_and_gaps_where_they_were():
    seq = "ACGTACGTAC" + "ggggcccc" + "ACGTTGCAAC" + "NNNN" + "TTAACCGGTT"
    out = m.shuffle_scannable(seq, random.Random(2), skip_repeats=True)
    assert out[10:18] == "ggggcccc" and out[28:32] == "NNNN"
    for a, b in ((0, 10), (18, 28), (32, 42)):
        assert sorted(out[a:b]) == sorted(seq[a:b])


def test_a_null_run_is_repeatable_and_has_the_requested_size():
    src = random.Random(5)
    seq = "".join(src.choice("ACGT") for _ in range(300))
    motifs = [_sharp("ACGTAC"), _sharp("GGATCC", mid="MA0002.1")]
    a = m.empirical_null(seq, motifs, pvalue=1e-3, min_ic=2.0, n_shuffles=12)
    b = m.empirical_null(seq, motifs, pvalue=1e-3, min_ic=2.0, n_shuffles=12)
    assert a["counts"] == b["counts"] and len(a["counts"]) == 12
    assert a["n_shuffles"] == 12 and len(a["scores"]) == 12


# ------------------------------------------------------- the scan itself
def test_a_motif_that_cannot_reach_the_pvalue_is_not_scanned():
    """A 4-base matrix: its best word has probability about 1/256, so p = 1e-5
    is unreachable, and the old scan reported every match of that word."""
    seq = "ACGT" * 40
    info: dict = {}
    hits = m.scan_promoter(seq, [_sharp("ACGT")], pvalue=1e-5, min_ic=0.0,
                           info=info)
    assert hits == []
    assert info == {"n_scanned": 0, "n_low_information": 0,
                    "n_unreachable": 1}
    assert pwm.threshold_from_tail(
        0, np.array([0.01, 0.004]), 1e-5, strict=True) is None
    assert pwm.threshold_from_tail(
        0, np.array([0.01, 0.004]), 1e-5) == 0.01, "the old fallback"


def test_a_reachable_motif_is_scanned_and_every_hit_passes_its_pvalue():
    rng = random.Random(9)
    seq = "".join(rng.choice("ACGT") for _ in range(600))
    seq = seq[:200] + "ACGTACGTACGT" + seq[200:]
    info: dict = {}
    hits = m.scan_promoter(seq, [_sharp("ACGTACGTACGT")], pvalue=1e-5,
                           min_ic=0.0, info=info)
    assert info["n_scanned"] == 1 and hits
    assert all(h.pvalue is not None and h.pvalue <= 1e-5 + 1e-12 for h in hits)
    assert any(h.start == 200 for h in hits)


def test_each_strand_has_its_own_threshold_when_the_bases_are_unbalanced():
    enc = pwm.encode("CCAAAAAACCGGTTTTTTGG")
    pssm = pwm.pfm_to_pssm(_sharp("AAAAAA").counts,
                           np.array([0.4, 0.1, 0.1, 0.4]))
    best = float(pssm.max(axis=1).sum())
    plus = list(pwm.scan_strands(enc, pssm, best - 0.01, None))
    assert plus and {s for _, _, s, _ in plus} == {"+"}
    minus = list(pwm.scan_strands(enc, pssm, None, best - 0.01))
    assert minus and {s for _, _, s, _ in minus} == {"-"}
    both = {s for _, _, s, _ in pwm.scan_strands(enc, pssm, best - 0.01,
                                                 best - 0.01)}
    assert both == {"+", "-"}
    assert list(pwm.scan_strands(enc, pssm, None, None)) == []


def test_the_reverse_strand_pvalue_comes_from_the_reverse_distribution():
    bg = np.array([0.5, 0.1, 0.1, 0.3])
    pssm = pwm.pfm_to_pssm(_sharp("AAAACC", n=30).counts, bg)
    bf, tf = pwm.score_tail(pssm, bg)
    br, tr = pwm.score_tail(pssm[::-1, ::-1], bg)
    # The best forward word is AAAACC and the best reverse word is GGTTTT, and
    # with A and T unequal they are not equally likely.
    assert tf[-1] != pytest.approx(tr[-1], rel=0.05)
    assert pwm.tail_pvalue(1e9, bf, tf) == tf[-1]
    assert pwm.tail_pvalue(-1e9, bf, tf) == tf[0]


# --------------------------------------------------------- the judgement
def test_without_a_null_the_formula_is_used_and_says_so():
    q = m.scan_quality(200, 800, 1200, 1e-5)
    assert q["null"] is None and q["expected_by_chance"] == q[
        "expected_by_formula"]
    assert "formula" in q["fdp_note"]


def test_with_a_null_the_expectation_is_its_mean_and_the_range_follows():
    null = {"counts": [20, 22, 24, 26, 28, 30, 22, 24, 26, 28, 24, 26],
            "scores": [[1e-6, 5e-6, 9e-6] + [9e-6] * (c - 3)
                       for c in [20, 22, 24, 26, 28, 30, 22, 24, 26, 28, 24, 26]],
            "n_shuffles": 12}
    obs = [1e-9] * 5 + [1e-6] * 5 + [9e-6] * 30
    q = m.scan_quality(40, 800, 1200, 1e-5, null=null, observed_scores=obs)
    assert q["expected_by_chance"] == pytest.approx(25.0, abs=0.1)
    assert q["expected_by_formula"] != q["expected_by_chance"]
    lo, hi = q["false_discovery_range"]
    assert lo <= q["false_discovery_proportion"] <= hi
    assert q["null"]["n_shuffles"] == 12
    assert "shuffled" in q["fdp_note"] and "formula" in q["fdp_note"]
    tops = {t["top"]: t["fdp"] for t in q["top_site_fdp"]}
    assert set(tops) == {5, 10, 20}
    assert tops[5] < tops[20], "the best sites are less likely to be chance"


def test_fewer_than_ten_shuffles_is_not_trusted():
    q = m.scan_quality(200, 800, 1200, 1e-5,
                       null={"counts": [1, 2, 3], "scores": [[], [], []],
                             "n_shuffles": 3})
    assert q["null"] is None and q["expected_by_chance"] == q[
        "expected_by_formula"]


def test_a_noise_scan_is_still_called_noise_with_a_null():
    null = {"counts": [900, 950, 1000, 980, 920, 960, 940, 990, 970, 930],
            "scores": [[] for _ in range(10)], "n_shuffles": 10}
    q = m.scan_quality(953, 800, 6000, 1e-4, null=null)
    assert q["verdict"] == "no signal" and q["false_discovery_proportion"] > 0.9


# ------------------------------------------------------ the run and report
def test_the_run_takes_the_option_and_the_demo_records_the_null():
    import inspect
    from promoter_tf_drugs import pipeline
    assert inspect.signature(pipeline.run).parameters[
        "null_shuffles"].default == m.NULL_SHUFFLES == 30
    res = pipeline.run("DEMOGENE", demo=True, null_shuffles=10)
    q = res["results"][0]["scan_quality"]
    assert q["null"] and q["null"]["n_shuffles"] == 10
    assert q["motifs_scanned"] is not None
    off = pipeline.run("DEMOGENE", demo=True, null_shuffles=0)
    assert off["results"][0]["scan_quality"]["null"] is None
    with pytest.raises(Exception):
        pipeline.run("DEMOGENE", demo=True, null_shuffles=-1)


def test_the_summary_line_gives_the_count_the_range_and_the_limit_of_a_site():
    from promoter_tf_drugs import pipeline, report
    res = pipeline.run("DEMOGENE", demo=True, null_shuffles=0)
    res["evidence_verdict"] = {
        "level": "degraded", "headline": "h", "reasons": [],
        "false_discovery_proportion": 0.455,
        "false_discovery_range": [0.235, 0.647], "scan_sites": 17}
    text = report.summary_report(res)
    assert "46% of the 17 motif sites" in text
    assert "range 24% to 65%" in text
    assert "not evidence of binding" in text


def test_a_window_with_unscanned_motifs_says_so_in_the_report():
    from promoter_tf_drugs import pipeline, report
    res = pipeline.run("DEMOGENE", demo=True, null_shuffles=0)
    res["results"][0]["scan_quality"]["motifs_unreachable"] = 254
    text = report.markdown_report(res)
    assert "254 of those motifs were not scanned" in text


# ---------------------------------------- widening a window with no signal
def _run_with_windows(tmp_path, monkeypatch, *, planted_narrow: bool,
                      planted_wide: bool, **run_kw):
    """An offline run whose narrow window does or does not hold the site and
    whose wide window does or does not."""
    import dataclasses

    from test_upstream_search import _lf_offline_run
    from promoter_tf_drugs import genes, pipeline
    real = genes.get_promoters
    site = "ACGTACGTACGT"
    asked = []

    def fake(client, gene, **kw):
        asked.append((kw["upstream"], kw["downstream"]))
        wide = kw["upstream"] >= pipeline.WIDE_UPSTREAM
        base = real(client, gene, **{**kw, "upstream": 500, "downstream": 100})
        planted = planted_wide if wide else planted_narrow
        out = []
        for p in base:
            # A sequence with no C or G cannot hold the site by chance.
            seq = ("AT" * len(p.sequence))[:len(p.sequence)]
            if planted:
                seq = seq[:60] + site + seq[60 + len(site):]
            out.append(dataclasses.replace(
                p, upstream=kw["upstream"], downstream=kw["downstream"],
                sequence=seq))
        return out
    monkeypatch.setattr(pipeline.genes, "get_promoters", fake)
    block, _ = _lf_offline_run(
        tmp_path, monkeypatch, edges=[],
        library=[{"chembl_id": "C1", "name": "x"}], index={},
        ledger_row="upstream", **run_kw)
    return block, asked


def test_a_window_with_no_signal_is_rescanned_wider_and_kept_if_it_helps(
        tmp_path, monkeypatch):
    block, asked = _run_with_windows(tmp_path, monkeypatch,
                                     planted_narrow=False, planted_wide=True)
    ww = block["scan_quality"]["window_widened"]
    assert ww and ww["kept"] is True
    assert ww["from"] == [500, 100] and ww["to"] == [5000, 1000]
    assert ww["sites_before"] == 0 and ww["sites_after"] >= 1
    assert block["scan_quality"]["verdict"] != "no signal"
    assert [a for a in asked if a[0] >= 5000], "the wide window was asked for"


def test_a_wider_window_that_does_not_help_is_not_kept(tmp_path, monkeypatch):
    block, _ = _run_with_windows(tmp_path, monkeypatch, planted_narrow=False,
                                 planted_wide=False)
    ww = block["scan_quality"]["window_widened"]
    assert ww and ww["kept"] is False
    assert block["scan_quality"]["verdict"] == "no signal"
    assert block["promoter"]["upstream"] == 500, "the requested window stays"


def test_a_window_with_signal_is_left_alone_and_no_widen_turns_it_off(
        tmp_path, monkeypatch):
    block, asked = _run_with_windows(tmp_path, monkeypatch,
                                     planted_narrow=True, planted_wide=True)
    assert block["scan_quality"]["window_widened"] is None
    assert not [a for a in asked if a[0] >= 5000]


def test_no_widen_keeps_the_requested_window(tmp_path, monkeypatch):
    block, asked = _run_with_windows(
        tmp_path, monkeypatch, planted_narrow=False, planted_wide=True,
        widen_no_signal=False)
    assert block["scan_quality"]["window_widened"] is None
    assert block["scan_quality"]["verdict"] == "no signal"
    assert not [a for a in asked if a[0] >= 5000]


def test_the_option_reaches_the_run_and_the_report_says_what_happened():
    import inspect
    from promoter_tf_drugs import cli, pipeline, report
    assert inspect.signature(pipeline.run).parameters[
        "widen_no_signal"].default is True
    kw = cli.run_kwargs_from_args(cli.build_parser().parse_args(
        ["--gene", "X", "--no-widen"]))
    assert kw["widen_no_signal"] is False
    res = pipeline.run("DEMOGENE", demo=True, null_shuffles=0)
    q = res["results"][0]["scan_quality"]
    q["window_widened"] = {"from": [1000, 200], "to": [5000, 1000],
                           "sites_before": 3, "sites_after": 40,
                           "chance_after": 20.5, "verdict_after": "weak",
                           "kept": True}
    text = report.markdown_report(res)
    assert "The window was widened" in text and "-5000/+1000" in text
    q["window_widened"]["kept"] = False
    assert "did not help either" in report.markdown_report(res)


def test_a_scan_with_no_sites_is_never_called_enriched():
    null = {"counts": [0] * 12, "scores": [[] for _ in range(12)],
            "n_shuffles": 12}
    assert m.scan_quality(0, 500, 800, 1e-5, null=null)["verdict"] == "no signal"
    assert m.scan_quality(0, 0, 800, 1e-5)["verdict"] == "no signal"
    assert m.scan_quality(3, 500, 800, 1e-5, null=null)["verdict"] == "enriched"
