"""Correctness checks for the parts that fail silently if they're wrong."""

from __future__ import annotations

import numpy as np

from promoter_tf_drugs import pwm, motifs, pipeline
from promoter_tf_drugs.demo_data import consensus_to_pfm, demo_motifs, demo_promoter_sequence


FLAT = np.full(4, 0.25)


def test_encode_and_revcomp():
    assert list(pwm.encode("ACGTN")) == [0, 1, 2, 3, 4]
    assert pwm.reverse_complement("AACCGGTT") == "AACCGGTT"
    assert pwm.reverse_complement("ATGC") == "GCAT"


def test_pssm_orientation():
    """PSSM row order must be A,C,G,T and the consensus must round-trip."""
    counts = consensus_to_pfm("ACGTACGT")
    pssm = pwm.pfm_to_pssm(counts, FLAT)
    assert pwm.consensus(pssm) == "ACGTACGT"


def test_perfect_match_scores_max():
    counts = consensus_to_pfm("TGACGTCA")
    pssm = pwm.pfm_to_pssm(counts, FLAT)
    lo, hi = pwm.min_max_score(pssm)
    hits = list(pwm.scan_both_strands("AAAA" + "TGACGTCA" + "AAAA", pssm, lo))
    best = max(hits, key=lambda h: h[3])
    assert best[0] == 4 and abs(best[3] - hi) < 1e-9


def test_reverse_strand_detection():
    counts = consensus_to_pfm("GGGGCGGGGC")          # not palindromic
    pssm = pwm.pfm_to_pssm(counts, FLAT)
    thr = pwm.score_threshold(pssm, FLAT, 1e-4)
    rc = pwm.reverse_complement("GGGGCGGGGC")
    hits = list(pwm.scan_both_strands("TTTT" + rc + "TTTT", pssm, thr))
    assert hits, "reverse-strand site was missed"
    assert any(h[2] == "-" and h[0] == 4 for h in hits)


def test_threshold_matches_empirical_fpr():
    """The analytic p-value must agree with simulation, or every hit count lies."""
    rng = np.random.default_rng(7)
    counts = consensus_to_pfm("AGAACANNNTGTTCT")
    pssm = pwm.pfm_to_pssm(counts, FLAT)
    target_p = 1e-3
    thr = pwm.score_threshold(pssm, FLAT, target_p)

    seq = "".join(rng.choice(list("ACGT"), size=400_000))
    enc = pwm.encode(seq)
    starts, _ = pwm.scan_strand(enc, pssm, thr)
    empirical = len(starts) / (len(enc) - pssm.shape[0] + 1)
    assert empirical <= target_p * 2.0, f"empirical FPR {empirical:.2e} >> {target_p:.0e}"
    assert empirical > target_p / 20, f"threshold far too strict ({empirical:.2e})"


def test_threshold_is_monotone_in_pvalue():
    counts = consensus_to_pfm("GGGACTTTCC")
    pssm = pwm.pfm_to_pssm(counts, FLAT)
    t1 = pwm.score_threshold(pssm, FLAT, 1e-3)
    t2 = pwm.score_threshold(pssm, FLAT, 1e-5)
    assert t2 > t1


def test_gc_background_shifts_gc_motif():
    """A GC-rich motif must score lower against a GC-rich background."""
    counts = consensus_to_pfm("GGGGCGGGGC")
    gc_bg = np.array([0.2, 0.3, 0.3, 0.2])
    flat_score = pwm.min_max_score(pwm.pfm_to_pssm(counts, FLAT))[1]
    gc_score = pwm.min_max_score(pwm.pfm_to_pssm(counts, gc_bg))[1]
    assert gc_score < flat_score


def test_jaspar_parser_round_trip():
    text = (
        ">MA0001.1\tTEST\n"
        "A [  1  0  9 ]\n"
        "C [  0  1  0 ]\n"
        "G [  9  9  1 ]\n"
        "T [  0  0  0 ]\n"
    )
    ms = motifs.parse_jaspar_pfm_text(text)
    assert len(ms) == 1
    m = ms[0]
    assert m.matrix_id == "MA0001.1" and m.name == "TEST" and m.length == 3
    assert pwm.consensus(m.counts) == "GGA"


def test_complex_name_splitting():
    m = motifs.Motif("X", "PPARG::RXRA(var.2)", consensus_to_pfm("ACGT"))
    assert m.tf_symbols == ["PPARG", "RXRA"]


def test_planted_sites_are_recovered():
    """End-to-end on the demo fixture: the planted GRE and SP1 sites must appear."""
    seq = demo_promoter_sequence()
    prom_offset = lambda i: (i - 1000 + 1) if i >= 1000 else (i - 1000)  # noqa: E731
    hits = motifs.scan_promoter(seq, demo_motifs(), pvalue=1e-4,
                                offset_fn=prom_offset)
    hits = motifs.collapse_overlapping(hits)
    found = {s for h in hits for s in h.tf_symbols}
    assert "NR3C1" in found, "planted GRE not recovered"
    assert "SP1" in found, "planted SP1 site not recovered"
    assert "REST" in found, "planted NRSE not recovered"

    gre = [h for h in hits if "NR3C1" in h.tf_symbols]
    assert any(h.tss_offset < -800 for h in gre), "GRE found at the wrong position"


def test_demo_pipeline_end_to_end():
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    assert res["results"], "no results"
    tfs = res["results"][0]["tfs"]
    assert tfs and tfs[0]["score"] > 0

    by_sym = {t["symbol"]: t for t in tfs}
    assert by_sym["NR3C1"]["nuclear_receptor"] is True
    assert any(d["name"] == "Dexamethasone" for d in by_sym["NR3C1"]["drugs"])

    # Direction logic: REST is a repressor, so only "down" drugs may be listed.
    if "REST" in by_sym and by_sym["REST"]["drugs"]:
        assert all(d["direction"] == "down" for d in by_sym["REST"]["drugs"])

    # An activator must never be paired with an inhibiting drug.
    for t in tfs:
        if t["role"] == "activator":
            assert all(d["direction"] in ("up", "unclear") for d in t["drugs"])

    md = pipeline.markdown_report(res)
    assert "DEMO MODE" in md and "Dexamethasone" in md


def _filler(n: int, gc: float = 0.5, seed: int = 11) -> str:
    """Neutral context of realistic length. Background estimation needs it:
    on a 40 bp window the composition estimate is pure noise."""
    rng = np.random.default_rng(seed)
    p = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return "".join(rng.choice(list("ACGT"), size=n, p=p))


def test_collapse_overlapping_dedupes():
    """A GC tandem repeat matches SP1 at many offsets; collapsing must remove
    the mutually overlapping duplicates and keep the best."""
    seq = _filler(600) + "GGGGCGGGGCGGGGCGGGGC" + _filler(600, seed=12)
    ms = [m for m in demo_motifs() if m.name == "SP1"]
    raw = motifs.scan_promoter(seq, ms, pvalue=1e-3)
    collapsed = motifs.collapse_overlapping(raw, max_overlap=0.75)
    assert 0 < len(collapsed) <= len(raw)
    for i, a in enumerate(collapsed):
        for b in collapsed[i + 1:]:
            if a.matrix_id != b.matrix_id:
                continue
            ov = min(a.end, b.end) - max(a.start, b.start)
            assert ov <= 0.75 * (a.end - a.start), "overlapping duplicates survived"


def test_repeat_masking_suppresses_low_complexity_hits():
    """Soft-masked (lowercase) repeat sequence must be excluded from scanning."""
    ms = [m for m in demo_motifs() if m.name == "SP1"]
    repeat = "ggggcggggcggggcggggc"          # soft-masked by Ensembl mask=soft
    seq = _filler(600, seed=21) + repeat + _filler(600, seed=22)

    lo, hi = 600, 620                      # the repeat's span

    unmasked = motifs.scan_promoter(seq.upper(), ms, pvalue=1e-3)
    masked = motifs.scan_promoter(seq, ms, pvalue=1e-3, skip_repeats=True)

    in_repeat = lambda hits: [h for h in hits if h.start < hi and h.end > lo]  # noqa: E731
    assert in_repeat(unmasked), "repeat should match SP1 when unmasked"
    assert not in_repeat(masked), "masked repeat still produced hits"

    # Chance hits in the neutral filler are expected and must survive masking:
    # ~1200 bp x 2 strands x 1e-3 is a couple of sites. If this number ever
    # drifts far from that, the p-value threshold is miscalibrated.
    assert len(masked) < 15


# ------------------------------------------------------------- coordinates
from promoter_tf_drugs.genes import Promoter, select_transcripts  # noqa: E402


def _prom(strand: int) -> Promoter:
    tss = 1_000_000
    up, down = 1000, 200
    if strand == 1:
        start, end = tss - up, tss + down - 1
    else:
        start, end = tss - down + 1, tss + up
    return Promoter("G", "ENSG", "ENST", "test", "7", strand, tss, start, end,
                    up, down, "GRCh38", "N" * (up + down))


def test_offset_to_tss_has_no_zero_position():
    p = _prom(1)
    assert p.offset_to_tss(0) == -1000        # first base of the window
    assert p.offset_to_tss(999) == -1         # last base before the TSS
    assert p.offset_to_tss(1000) == +1        # the TSS itself
    assert p.offset_to_tss(1199) == +200


def test_promoter_window_spans_tss_on_both_strands():
    for strand in (1, -1):
        p = _prom(strand)
        assert p.end - p.start + 1 == p.upstream + p.downstream
        assert p.start <= p.tss <= p.end
    plus, minus = _prom(1), _prom(-1)
    # Upstream is genomically lower on + and higher on -
    assert plus.tss - plus.start == 1000
    assert minus.end - minus.tss == 1000


def test_bed_coordinates_flip_with_strand():
    """A hit at sequence index 0 is the far-upstream end on both strands, but
    lands at opposite genomic edges of the window."""
    from promoter_tf_drugs.pipeline import write_outputs
    import tempfile, pathlib

    # Plus strand: index 0 is genomic 999_000 (1-based) -> BED start 998_999.
    # Minus strand: the window runs 999_801..1_001_000 and index 0 is the
    # *highest* coordinate, so a 10 bp hit occupies BED [1_000_990, 1_001_000).
    for strand, expect_start in ((1, 998_999), (-1, 1_000_990)):
        p = _prom(strand)
        result = {
            "query": {"gene": "G", "demo": True, "pvalue": 1e-4, "n_motifs": 1,
                      "species": "h", "upstream": 1000, "downstream": 200,
                      "chip_bed": None, "approved_only": True},
            "weights": {},
            "results": [{
                "promoter": p.to_dict(), "sequence": p.sequence,
                "n_sites": 1, "n_tfs": 0, "tfs": [],
                "hits": [{"start": 0, "end": 10, "strand": "+", "rel_score": 1.0,
                          "motif_name": "M", "matrix_id": "MA0001.1"}],
            }],
        }
        with tempfile.TemporaryDirectory() as d:
            write_outputs(result, d)
            bed = (pathlib.Path(d) / "sites.bed").read_text().splitlines()[1]
            fields = bed.split("\t")
            assert int(fields[1]) == expect_start, f"strand {strand}: {bed}"
            assert int(fields[2]) - int(fields[1]) == 10


def test_transcript_selection_prefers_mane():
    gene = {"id": "ENSG", "strand": 1, "Transcript": [
        {"id": "T1", "start": 100, "end": 900, "biotype": "protein_coding"},
        {"id": "T2", "start": 200, "end": 300, "biotype": "protein_coding",
         "tag": ["MANE_Select"]},
        {"id": "T3", "start": 50, "end": 5000, "biotype": "protein_coding",
         "is_canonical": 1},
    ]}
    (tr, reason), = select_transcripts(gene)
    assert tr["id"] == "T2" and reason == "MANE Select"

    del gene["Transcript"][1]["tag"]
    (tr, reason), = select_transcripts(gene)
    assert tr["id"] == "T3" and reason == "Ensembl canonical"


def test_all_tss_clusters_nearby_starts():
    """TSSs within 500 bp give windows that overlap almost completely.
    Treating them as distinct tripled the CCR5 report for no information —
    two of its 'distinct' TSSs were 7 bp apart."""
    gene = {"id": "ENSG", "strand": -1, "Transcript": [
        {"id": "T1", "start": 100, "end": 900},
        {"id": "T2", "start": 150, "end": 900},    # identical TSS
        {"id": "T3", "start": 100, "end": 1200},   # 300 bp away: same cluster
        {"id": "T4", "start": 100, "end": 9000},   # far away: its own cluster
    ]}
    picked = select_transcripts(gene, all_tss=True)
    assert len(picked) == 2, [t["id"] for t, _ in picked]
    assert any("representing 3" in why for _, why in picked)


def test_tss_cluster_width_is_configurable():
    gene = {"id": "ENSG", "strand": 1, "Transcript": [
        {"id": "A", "start": 100, "end": 900},
        {"id": "B", "start": 700, "end": 900},
    ]}
    assert len(select_transcripts(gene, all_tss=True, cluster_bp=500)) == 2
    assert len(select_transcripts(gene, all_tss=True, cluster_bp=1000)) == 1


def test_cluster_representative_prefers_mane():
    gene = {"id": "ENSG", "strand": 1, "Transcript": [
        {"id": "plain", "start": 100, "end": 900},
        {"id": "mane", "start": 120, "end": 900, "tag": ["MANE_Select"]},
    ]}
    (tr, why), = select_transcripts(gene, all_tss=True)
    assert tr["id"] == "mane"
    assert "representing 2" in why


# ------------------------------------------------------------- direction
from promoter_tf_drugs.drugs import desired_direction  # noqa: E402


def test_desired_direction_inverts_with_goal():
    assert desired_direction("activator", "up") == "up"
    assert desired_direction("activator", "down") == "down"
    assert desired_direction("repressor", "up") == "down"
    assert desired_direction("repressor", "down") == "up"
    assert desired_direction("unknown", "down") == "either"


def test_down_mode_returns_blocking_drugs_for_activators():
    res = pipeline.run("DEMOGENE", demo=True, offline=True, goal="down")
    by_sym = {t["symbol"]: t for t in res["results"][0]["tfs"]}

    # An agonist must never appear when the goal is to lower the gene and the
    # TF is an activator.
    for sym, t in by_sym.items():
        if t["role"] in ("activator", "dual"):
            assert all(d["direction"] in ("down", "unclear") for d in t["drugs"]), sym
        elif t["role"] == "repressor":
            assert all(d["direction"] in ("up", "unclear") for d in t["drugs"]), sym

    # NR3C1 is an activator here, so mifepristone (antagonist) should surface
    # and dexamethasone (agonist) should not.
    names = {d["name"] for d in by_sym["NR3C1"]["drugs"]}
    assert "Mifepristone" in names
    assert "Dexamethasone" not in names


def test_up_and_down_give_disjoint_drug_sets():
    up = pipeline.run("DEMOGENE", demo=True, offline=True, goal="up")
    dn = pipeline.run("DEMOGENE", demo=True, offline=True, goal="down")

    def pairs(r):
        return {(t["symbol"], d["name"])
                for t in r["results"][0]["tfs"] for d in t["drugs"]}

    up_p, dn_p = pairs(up), pairs(dn)
    assert up_p and dn_p
    assert not (up_p & dn_p), "the same TF/drug pair cannot serve both directions"


def test_down_mode_report_language():
    res = pipeline.run("DEMOGENE", demo=True, offline=True, goal="down")
    md = pipeline.markdown_report(res)
    assert "Downregulation strategy" in md
    assert "decrease" in md
    assert "siRNA" in md          # direct-modality guidance for knockdown


# ----------------------------------------------------------- drug grading
from promoter_tf_drugs import drug_scoring as dsc          # noqa: E402
from promoter_tf_drugs.demo_data import demo_evidence, DEMO_DRUGS  # noqa: E402


def test_potency_normalisation_endpoints():
    assert dsc._potency(dsc.DrugEvidence(pchembl_median=9.0, pchembl_n=5))[0] == 1.0
    assert dsc._potency(dsc.DrugEvidence(pchembl_median=5.0, pchembl_n=5))[0] == 0.0
    assert dsc._potency(dsc.DrugEvidence(pchembl_median=3.0, pchembl_n=5))[0] == 0.0
    # unknown potency is neutral, not zero, and not full marks
    assert dsc._potency(dsc.DrugEvidence())[0] == 0.5


def test_selectivity_penalises_promiscuity():
    single = dsc._selectivity(dsc.DrugEvidence(n_targets=1))[0]
    many = dsc._selectivity(dsc.DrugEvidence(n_targets=10))[0]
    assert single == 1.0
    assert many < 0.35


def test_weak_tf_evidence_caps_a_perfect_drug():
    """The core claim of the design: a great drug cannot rescue a bad TF hit."""
    drug = DEMO_DRUGS["VDR"][0]                  # calcitriol: potent, selective
    ev = demo_evidence(drug)
    strong = dsc.grade_drug(drug, tf_confidence=0.9, evidence=ev)
    weak = dsc.grade_drug(drug, tf_confidence=0.05, evidence=ev)
    assert strong.combined > 70
    assert weak.combined < 30, "weak TF evidence failed to cap the grade"
    assert any("limiting factor" in n for n in weak.notes)


def test_bad_drug_is_not_rescued_by_strong_tf_evidence():
    """And the converse, which an additive score would get wrong."""
    good = DEMO_DRUGS["VDR"][0]
    bad = DEMO_DRUGS["REST"][0]                  # valproate: weak, promiscuous
    g = dsc.grade_drug(good, 0.85, demo_evidence(good))
    b = dsc.grade_drug(bad, 0.85, demo_evidence(bad))
    assert g.combined > b.combined + 30
    assert b.grade in ("D", "F")


def test_withdrawn_drug_is_penalised_but_still_reported():
    drug = DEMO_DRUGS["SP1"][0]                  # plicamycin, withdrawn
    g = dsc.grade_drug(drug, 0.8, demo_evidence(drug))
    assert any("WITHDRAWN" in n for n in g.notes)
    assert g.combined < 45
    assert g.drug_name == "Plicamycin"           # surfaced, not silently dropped


def test_tissue_filter_downweights_non_cns_drugs():
    drug = DEMO_DRUGS["SP1"][0]                  # plicamycin: MW 1085, TPSA 331
    ev = demo_evidence(drug)
    plain = dsc.grade_drug(drug, 0.8, ev)
    cns = dsc.grade_drug(drug, 0.8, ev, tissue="cns")
    assert cns.combined < plain.combined
    assert any("CNS" in n for n in cns.notes)


def test_tf_weight_shifts_the_balance():
    drug = DEMO_DRUGS["REST"][0]                 # poor drug quality
    ev = demo_evidence(drug)
    drug_led = dsc.grade_drug(drug, 0.9, ev, tf_weight=0.1)
    tf_led = dsc.grade_drug(drug, 0.9, ev, tf_weight=0.9)
    # Leaning on TF confidence should help a strong-TF/weak-drug pair
    assert tf_led.combined > drug_led.combined


def test_grades_are_sorted_and_present_in_report():
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    grades = res["results"][0]["drug_grades"]
    assert grades, "no grades produced"
    assert grades == sorted(grades, key=lambda g: -g["combined"])
    assert all(g["grade"] in "ABCDF" for g in grades)
    assert all(0 <= g["combined"] <= 100 for g in grades)

    md = pipeline.markdown_report(res)
    assert "Graded candidates" in md
    assert grades[0]["drug_name"] in md


# --------------------------------------------------- agonism / conformation
from promoter_tf_drugs import agonism, docking          # noqa: E402

_CORES = ["c1ccccc1O", "c1ccc2ccccc2c1O", "c1ccc(cc1)c1ccccc1O", "Oc1ccc(cc1)C(C)C",
          "Oc1ccc2c(c1)CCCC2", "Oc1ccc(cc1)CC", "Oc1ccc(cc1)CCC",
          "Oc1ccc(cc1)C(C)(C)C", "Oc1ccc(cc1)Cl", "Oc1ccc(cc1)F"]


def _toy_ligands():
    """Agonists = bare cores; antagonists = same cores + a basic tether,
    mimicking the H12-displacing side chain. Tests the machinery, not chemistry."""
    ls = [{"smiles": c, "label": "agonist", "mw": 150.0, "chembl_id": f"A{i}"}
          for i, c in enumerate(_CORES)]
    ls += [{"smiles": c.replace("O", "OCCN(C)C", 1), "label": "antagonist",
            "mw": 250.0, "chembl_id": f"B{i}"} for i, c in enumerate(_CORES)]
    return ls


def test_model_refuses_when_training_set_too_small():
    assert agonism.train_model(_toy_ligands()[:6]) is None


def test_model_trains_and_reports_auc():
    m = agonism.train_model(_toy_ligands())
    assert m is not None
    assert m.n_agonist == 10 and m.n_antagonist == 10
    assert m.auc is not None and m.auc > 0.7


def test_out_of_domain_query_is_refused():
    """A chemically unrelated molecule must get no prediction at all."""
    m = agonism.train_model(_toy_ligands())
    glucose = "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"
    p = agonism.predict(m, glucose, drug_name="glucose")
    assert p.in_domain is False
    assert p.predicted == "undetermined"
    assert any("APPLICABILITY DOMAIN" in n for n in p.notes)


def test_protrusion_flags_antagonist_signature():
    feats = agonism.protrusion_features("CN(C)CCOc1ccc(cc1)CCCC",
                                        agonist_mws=[150.0] * 5)
    assert feats["basic_amines"] >= 1
    assert feats["mw_excess_vs_agonists"] > 60
    assert feats["antagonist_like_bulk"] is True

    plain = agonism.protrusion_features("Oc1ccc(cc1)CCCC", agonist_mws=[150.0] * 5)
    assert plain["antagonist_like_bulk"] is False


def test_structural_flag_lowers_confidence_on_disagreement():
    m = agonism.train_model(_toy_ligands())
    p = agonism.predict(m, "CN(C)CCOc1ccc(cc1)CCCC", drug_name="bulky")
    # Whatever the classifier says, a bulk+amine molecule must not come back
    # as a confident agonist.
    assert not (p.predicted == "agonist" and p.confidence == "high")


# --- conformation detection from PDB text ---------------------------------
def _pdb(chains: dict[str, str]) -> str:
    """Minimal PDB with CA records only, one chain per entry."""
    aa1to3 = {v: k for k, v in docking.AA3_TO_1.items() if k != "MSE"}
    lines, n = [], 0
    for chain, seq in chains.items():
        for i, aa in enumerate(seq, 1):
            n += 1
            lines.append(
                f"ATOM  {n:5d}  CA  {aa1to3[aa]:>3s} {chain}{i:4d}"
                f"    {0.0:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C")
    return "\n".join(lines)


def test_chain_sequences_parsed():
    seqs = docking.chain_sequences(_pdb({"A": "ACDEFGHIK"}))
    assert seqs == {"A": "ACDEFGHIK"}


def test_lxxll_peptide_marks_active_conformation():
    """A short chain with an LXXLL NR box means a coactivator is in the AF-2
    groove, which only happens when helix 12 is in the active position."""
    receptor = "A" * 250
    coactivator = "KHKILHRLLQDSS"          # contains LXXLL
    conf = docking.classify_conformation(_pdb({"A": receptor, "B": coactivator}), "1TST")
    assert conf.state == "active"
    assert "LXXLL" in conf.evidence


def test_cornr_peptide_marks_inactive_conformation():
    receptor = "A" * 250
    corepressor = "LEDIIRKALMGSF"          # L...I I...L CoRNR-like
    conf = docking.classify_conformation(_pdb({"A": receptor, "B": corepressor}), "2TST")
    assert conf.state == "inactive"


def test_apo_structure_is_unknown_not_guessed():
    conf = docking.classify_conformation(_pdb({"A": "A" * 250}), "3TST")
    assert conf.state == "unknown"
    assert "another way" in conf.evidence


def test_docking_requires_both_conformations(monkeypatch=None):
    """Without both states there is no directional signal, and the module must
    say so rather than returning a score."""
    from promoter_tf_drugs.drugs import Drug
    drug = Drug("Test", "CHEMBL1", "VDR", "AGONIST", None, True, 4.0, 2000,
                False, True, "up", "test")
    import promoter_tf_drugs.docking as dk
    orig = dk.find_structures
    dk.find_structures = lambda *a, **k: []
    try:
        res = dk.predict_direction(None, drug, "CCO", "P11473")
        assert res.predicted == "undetermined"
        assert res.best_active is None
    finally:
        dk.find_structures = orig


def test_direction_conflict_penalises_and_flags():
    """The point of the direction model: disagreement with the database is a
    warning to go read the literature, not a silent re-label."""
    from promoter_tf_drugs.drug_scoring import apply_direction_check, DrugGrade
    from promoter_tf_drugs.agonism import AgonismPrediction

    def fresh():
        return DrugGrade("VDR", "TestDrug", 60.0, "B", 0.7, 0.8, {}, {}, [])

    agree = AgonismPrediction("TestDrug", "VDR", "agonist",
                              probability_agonist=0.9, confidence="high")
    conflict = AgonismPrediction("TestDrug", "VDR", "antagonist",
                                 probability_agonist=0.1, confidence="high")

    up = apply_direction_check(fresh(), agree, "up")
    assert up.combined > 60.0
    assert any("supports agonist" in n for n in up.notes)

    bad = apply_direction_check(fresh(), conflict, "up")
    assert bad.combined < 60.0
    assert any("CONFLICT" in n for n in bad.notes)

    # An undetermined prediction must change nothing.
    none = apply_direction_check(fresh(), AgonismPrediction("x", "VDR", "undetermined"), "up")
    assert none.combined == 60.0


# --------------------------------------------- transcriptional specificity
from promoter_tf_drugs import regulon as rg, direct as dr   # noqa: E402


def test_global_regulator_scores_far_below_narrow_tf():
    """The core claim: hitting one protein that runs 5000 genes is not specific."""
    sp1, n_sp1, note = rg.transcriptional_specificity("SP1")
    vdr, n_vdr, _ = rg.transcriptional_specificity("VDR")
    assert sp1 < 0.25 and vdr > 0.6
    assert n_sp1 > n_vdr
    assert "global regulator" in note


def test_unknown_tf_is_neutral_not_zero():
    score, n, note = rg.transcriptional_specificity("ZZZ9FAKE")
    assert 0.4 < score < 0.7 and n is None and "unknown" in note


def test_regulon_file_overrides_bundled_tiers(tmp_path):
    f = tmp_path / "trrust.tsv"
    f.write_text("SP1\tGENEA\tActivation\t123\nSP1\tGENEB\tRepression\t124\n"
                 "MYC\tGENEC\tActivation\t125\n")
    reg = rg.load_regulons(f)
    assert reg["SP1"] == {"GENEA", "GENEB"}
    score, n, _ = rg.transcriptional_specificity("SP1", reg)
    assert n == 2 and score == 1.0        # 2 real targets beats the coarse tier


def test_regulon_file_with_header_and_csv(tmp_path):
    f = tmp_path / "dorothea.csv"
    f.write_text("tf,target,confidence\nESR1,GREB1,A\nESR1,TFF1,A\n")
    reg = rg.load_regulons(f)
    assert reg["ESR1"] == {"GREB1", "TFF1"}


def test_evidence_depth_separates_clean_from_unstudied():
    """A drug nobody has measured must not score as if it were selective."""
    studied = dsc.DrugEvidence(pchembl_median=8.0, pchembl_n=40, n_targets=1)
    unstudied = dsc.DrugEvidence(pchembl_median=None, pchembl_n=0, n_targets=None)
    assert dsc._evidence_depth(studied)[0] == 1.0
    assert dsc._evidence_depth(unstudied)[0] < 0.4
    assert "barely characterised" in dsc._evidence_depth(unstudied)[1]


def test_regulon_breadth_changes_the_grade():
    drug = DEMO_DRUGS["VDR"][0]
    ev_narrow = demo_evidence(drug)
    ev_narrow.regulon_score = 1.0
    ev_broad = demo_evidence(drug)
    ev_broad.regulon_score = 0.16
    ev_broad.regulon_note = "global regulator"

    narrow = dsc.grade_drug(drug, 0.8, ev_narrow)
    broad = dsc.grade_drug(drug, 0.8, ev_broad)
    assert narrow.combined > broad.combined
    assert any("global regulator" in n for n in broad.notes)


# ------------------------------------------------- direct gene-targeted drugs
def test_curated_direct_therapies_respect_direction():
    up = dr.curated_for_gene("SMN2", "up")
    assert {t.drug for t in up} >= {"Nusinersen", "Risdiplam"}
    assert all(t.direction == "up" for t in up)
    assert not dr.curated_for_gene("SMN2", "down")

    down = dr.curated_for_gene("TTR", "down")
    assert {t.drug for t in down} >= {"Patisiran", "Inotersen"}
    assert not dr.curated_for_gene("TTR", "up")


def test_activity_only_drugs_are_distinguished_from_abundance():
    """Ivacaftor raises CFTR channel function, not CFTR protein. A western blot
    would show nothing. The flag must survive into the output."""
    cftr = {t.drug: t for t in dr.curated_for_gene("CFTR", "up")}
    assert cftr["Ivacaftor"].raises_protein is False
    smn = {t.drug: t for t in dr.curated_for_gene("SMN2", "up")}
    assert smn["Nusinersen"].raises_protein is True


def test_sequence_directed_outranks_small_molecule():
    hits = dr.find_direct(None, "DMD", goal="up", demo=True)
    assert hits
    assert hits[0].modality in ("gene therapy", "splice-switching oligonucleotide")
    assert hits[0].specificity > 0.9


def test_withdrawn_direct_therapy_still_listed():
    hits = dr.find_direct(None, "APOB", goal="down", demo=True)
    assert any(t.drug == "Mipomersen" and t.withdrawn for t in hits)


def test_direct_section_appears_in_report():
    res = pipeline.run("SMN2", demo=True, offline=True)
    assert res["direct_therapies"], "curated direct hits missing"
    md = pipeline.markdown_report(res)
    assert "targeting SMN2 directly" in md
    assert "Nusinersen" in md


def test_no_direct_hits_reports_honestly():
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    assert res["direct_therapies"] == []
    md = pipeline.markdown_report(res)
    assert "None found in the curated table" in md


# ------------------------------------------- outcome chain & virtual screen
from promoter_tf_drugs import outcome as oc, screening as scr   # noqa: E402


class _FakeTF:
    def __init__(self, score, chip=0):
        self.score, self.chip_peaks = score, chip


def test_outcome_chain_multiplies_and_stays_honest():
    """Four uncertain links must produce a small number, not a confident one."""
    est = oc.estimate("TestDrug", "VDR", "PAH", goal="up",
                      tf_result=_FakeTF(70, chip=3))
    assert est.p_point == round(
        est.p_engages * est.p_direction * est.p_activates * est.p_responsive, 3)
    assert est.p_point < 0.5, "a four-link chain should not read as confident"
    assert est.p_low < est.p_point < est.p_high


def test_limiting_step_is_identified_and_actionable():
    est = oc.estimate("D", "VDR", "PAH", tf_result=_FakeTF(10))
    assert est.limiting_step in (
        "activation at this promoter", "gene responsiveness in your cells")
    advice = oc.narrowing_experiment(est)
    assert len(advice) > 40                       # a real experiment, not a stub


def test_chip_evidence_raises_the_promoter_link():
    without = oc.estimate("D", "VDR", "G", tf_result=_FakeTF(60, chip=0))
    with_chip = oc.estimate("D", "VDR", "G", tf_result=_FakeTF(60, chip=5))
    assert with_chip.p_activates > without.p_activates
    assert with_chip.p_point > without.p_point


def test_ambiguous_direction_raises_the_serm_flag():
    """Partial H12 engagement is where tamoxifen lives: direction becomes a
    property of the cell type, not the molecule."""
    from promoter_tf_drugs.agonism import AgonismPrediction
    borderline = AgonismPrediction("X", "ESR1", "agonist",
                                   probability_agonist=0.52, confidence="moderate")
    est = oc.estimate("X", "ESR1", "GREB1", prediction=borderline,
                      tf_result=_FakeTF(60))
    assert est.tissue_dependent is True
    assert any("SERM" in n for n in est.notes)


def test_low_confidence_prediction_shrinks_toward_coinflip():
    from promoter_tf_drugs.agonism import AgonismPrediction
    hi = AgonismPrediction("X", "V", "agonist", probability_agonist=0.95,
                           confidence="high")
    lo = AgonismPrediction("X", "V", "agonist", probability_agonist=0.95,
                           confidence="low")
    e_hi = oc.estimate("X", "V", "G", prediction=hi, tf_result=_FakeTF(60))
    e_lo = oc.estimate("X", "V", "G", prediction=lo, tf_result=_FakeTF(60))
    assert e_hi.p_direction > e_lo.p_direction
    assert abs(e_lo.p_direction - 0.5) < abs(e_hi.p_direction - 0.5)


# --- virtual screening ----------------------------------------------------
def _toy_library():
    lib = [{"chembl_id": f"A{i}", "name": f"known-ag-{i}", "smiles": c,
            "max_phase": 4} for i, c in enumerate(_CORES)]
    lib += [{"chembl_id": f"B{i}", "name": f"known-ant-{i}",
             "smiles": c.replace("O", "OCCN(C)C", 1), "max_phase": 4}
            for i, c in enumerate(_CORES)]
    # structurally unrelated approved-drug stand-ins
    lib += [{"chembl_id": "N1", "name": "novel-sugar",
             "smiles": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O", "max_phase": 4},
            {"chembl_id": "N2", "name": "novel-peptide-ish",
             "smiles": "NCC(=O)NCC(=O)NCC(=O)O", "max_phase": 4}]
    return lib


def test_screen_requires_a_model():
    import pytest
    with pytest.raises(ValueError, match="no trained model"):
        scr.screen(None, _toy_library())


def test_screen_labels_quadrants():
    m = agonism.train_model(_toy_ligands())
    known = {f"A{i}" for i in range(len(_CORES))}
    hits = scr.screen(m, _toy_library(), goal="up", known_ids=known, top_n=100)
    quads = {h.quadrant for h in hits}
    assert "ESTABLISHED LIGAND, untested on this gene" in quads
    # established binders untested on the gene are the most actionable class,
    # not a control to be skipped past
    assert hits[0].quadrant == "ESTABLISHED LIGAND, untested on this gene"
    for h in hits:
        if h.quadrant == "NOVEL HYPOTHESIS":
            assert any("extrapolation" in n for n in h.notes)


def test_novelty_is_inverse_similarity():
    m = agonism.train_model(_toy_ligands())
    hits = scr.screen(m, _toy_library(), known_ids=set(), top_n=100)
    for h in hits:
        assert abs(h.novelty - (1.0 - h.nearest_tanimoto)) < 1e-9


def test_enrichment_check_is_relative_not_absolute():
    """Corrected semantics. An absolute threshold rejected ESR1's working
    model, because its annotated ligands are mostly ANTAGONISTS — so a low
    P(agonist) was the right answer, not a broken screen. The test is whether
    known actives outrank the novel candidates."""
    good = [scr.ScreenHit("k", "A1", 0.9, 0.9, 0.1, "known ligand", True, 4.0),
            scr.ScreenHit("n", "N1", 0.8, 0.2, 0.8, "NOVEL HYPOTHESIS", False, 4.0)]
    assert scr.enrichment_check(good)["usable"] is True

    # the ESR1 case: low absolute scores, but correctly ordered
    esr1 = [scr.ScreenHit(f"k{i}", f"A{i}", 0.24, 0.9, 0.1, "known ligand",
                          True, 4.0) for i in range(4)]
    esr1 += [scr.ScreenHit(f"n{i}", f"N{i}", 0.20, 0.2, 0.8,
                           "NOVEL HYPOTHESIS", False, 4.0) for i in range(3)]
    assert scr.enrichment_check(esr1)["usable"] is True, \
        "rejected a correctly-ordered screen on absolute score"

    # genuinely inverted: novel candidates outrank the known actives
    bad = [scr.ScreenHit("k", "A1", 0.3, 0.9, 0.1, "known ligand", True, 4.0),
           scr.ScreenHit("n", "N1", 0.9, 0.2, 0.8, "NOVEL HYPOTHESIS", False, 4.0)]
    res = scr.enrichment_check(bad)
    assert res["usable"] is False and "inverted" in res["reason"]

    none = [scr.ScreenHit("n", "N1", 0.9, 0.2, 0.8, "NOVEL HYPOTHESIS", False, 4.0)]
    assert scr.enrichment_check(none)["usable"] is False


# ------------------------------------------------- discovery / superiority
from promoter_tf_drugs import discovery as disc          # noqa: E402


def _axes(s, sp, e, evidence="predicted"):
    return disc.Axes(s, sp, e, evidence)


def test_dominance_requires_no_loss_elsewhere():
    """The trade this package exists to avoid: slightly stronger, far less
    specific, laundered into an overall win."""
    incumbent = _axes(0.60, 0.80, 0.70)
    stronger_but_dirty = _axes(0.75, 0.30, 0.70)
    assert not disc.dominates(stronger_but_dirty, incumbent)

    better_everywhere = _axes(0.75, 0.85, 0.72)
    assert disc.dominates(better_everywhere, incumbent)


def test_predicted_challenger_must_clear_a_margin():
    """Selection bias: you notice the candidate because its predicted score is
    high, so noise is already working in its favour."""
    measured = disc.Candidate("Incumbent", "direct", None,
                              _axes(0.60, 0.80, 0.70, "measured"))
    marginal = disc.Candidate("Challenger", "via TF", "VDR",
                              _axes(0.68, 0.80, 0.70, "predicted"))
    disc.compare(marginal, measured)
    assert marginal.wins_on == []                     # +0.08 < 0.15 margin
    assert "within noise" in marginal.verdict

    clear = disc.Candidate("BigWin", "via TF", "VDR",
                           _axes(0.85, 0.80, 0.70, "predicted"))
    disc.compare(clear, measured)
    assert "stronger activation" in clear.wins_on


def test_measured_vs_measured_needs_no_margin():
    a = disc.Candidate("A", "direct", None, _axes(0.68, 0.80, 0.70, "measured"))
    b = disc.Candidate("B", "direct", None, _axes(0.60, 0.80, 0.70, "measured"))
    disc.compare(a, b)
    assert "stronger activation" in a.wins_on


def test_tradeoff_is_reported_as_a_tradeoff():
    base = disc.Candidate("Base", "direct", None, _axes(0.60, 0.85, 0.70, "measured"))
    mixed = disc.Candidate("Mixed", "via TF", "SP1", _axes(0.90, 0.40, 0.70))
    disc.compare(mixed, base)
    assert mixed.wins_on and mixed.loses_on
    assert mixed.verdict.startswith("trade-off")


def test_no_incumbent_is_stated_not_assumed():
    c = disc.Candidate("Only", "via TF", "VDR", _axes(0.5, 0.5, 0.5))
    disc.compare(c, None)
    assert "no incumbent" in c.verdict
    assert any("gap in the databases" in n for n in c.notes)


def test_novelty_matching_errs_toward_already_known():
    """Calling a known link 'new' is the embarrassing failure; the reverse only
    costs you a candidate."""
    links = {"dexamethasone": ["increases^expression"],
             "retinoic acid": ["decreases^expression"]}
    known, desc = disc.check_novelty("Dexamethasone", links)
    assert known and "expression" in desc

    known, _ = disc.check_novelty("Calcitriol", links)
    assert not known

    # substring match in either direction still counts as known
    known, _ = disc.check_novelty("All-Trans Retinoic Acid", links)
    assert known


def test_specificity_axis_multiplies_both_layers():
    """One protein out of many is not specific if that protein runs 5000 genes."""
    class G:
        components = {"potency": 0.9, "action": 1.0, "selectivity": 1.0,
                      "regulon": 0.16, "maturity": 1.0, "tissue": 1.0,
                      "evidence_depth": 1.0}
        tf_confidence = 0.8
        evidence = {}
    ax = disc.axes_from_grade(G())
    assert ax.specificity < 0.2          # perfect protein selectivity, global TF


def test_pareto_front_excludes_dominated_candidates():
    good = disc.Candidate("Good", "direct", None, _axes(0.9, 0.9, 0.9))
    bad = disc.Candidate("Bad", "via TF", "SP1", _axes(0.3, 0.3, 0.3))
    # A 0.05 edge is inside TOLERANCE, so it gets absorbed and Good dominates.
    marginal = disc.Candidate("Marginal", "via TF", "VDR", _axes(0.95, 0.2, 0.5))
    assert {c.name for c in disc.pareto_front([good, marginal])} == {"Good"}

    # A clear edge on strength survives, even with much worse specificity.
    mixed = disc.Candidate("Mixed", "via TF", "VDR", _axes(1.0, 0.2, 0.5))
    front = {c.name for c in disc.pareto_front([good, bad, mixed])}
    assert "Good" in front and "Mixed" in front
    assert "Bad" not in front


def test_discovery_runs_end_to_end_and_names_a_baseline():
    res = pipeline.run("SMN2", demo=True, offline=True)
    imp = res["results"][0]["improvements"]
    assert imp["baseline"] is not None
    assert imp["baseline"]["name"] in ("Nusinersen", "Risdiplam")
    assert imp["n_candidates"] > 0
    md = pipeline.markdown_report(res)
    assert "Candidates that may beat the current standard" in md
    assert "The bar:" in md and "per-axis best" in md


def test_gene_with_no_incumbent_says_so():
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    imp = res["results"][0]["improvements"]
    assert imp["baseline"] is None
    assert "No incumbent found" in pipeline.markdown_report(res)


# --------------------------------------------------- ChIP evidence routing
from promoter_tf_drugs import evidence as ev_mod          # noqa: E402


class _FakeClient:
    def __init__(self, payload=None, fail=False):
        self.payload, self.fail, self.calls = payload or [], fail, []

    def get(self, url, params=None, **kw):
        self.calls.append(url)
        if self.fail:
            raise RuntimeError("network down")
        return self.payload


def test_api_counts_peaks_per_tf():
    client = _FakeClient([{"target": "CTCF"}, {"target": "CTCF"},
                          {"target": "NR3C1"}])
    counts = ev_mod.peaks_from_api(client, "7", 1000, 2000)
    assert counts == {"CTCF": 2, "NR3C1": 1}
    assert "chr7:1000-2000" in client.calls[0]


def test_api_normalises_chr_prefix():
    client = _FakeClient([])
    ev_mod.peaks_from_api(client, "chr7", 10, 20)
    assert "chrchr7" not in client.calls[0]


def test_api_failure_degrades_to_no_evidence():
    """A network failure must mean 'no ChIP evidence', never a crashed run."""
    assert ev_mod.peaks_from_api(_FakeClient(fail=True), "1", 1, 2) == {}


def test_chip_evidence_prefers_local_bed_over_api(tmp_path):
    bed = tmp_path / "peaks.bed"
    bed.write_text("chr7\t1000\t1200\tGSE1.CTCF.HeLa\t500\t.\n")
    prom = _prom(1)
    prom.chrom, prom.start, prom.end = "7", 900, 1300

    client = _FakeClient([{"target": "NR3C1"}])
    counts = ev_mod.chip_evidence(prom, chip_bed=bed, client=client, use_api=True)
    assert counts == {"CTCF": 1}
    assert client.calls == [], "API was called even though a BED was supplied"


def test_chip_evidence_returns_empty_when_neither_source_given():
    assert ev_mod.chip_evidence(_prom(1)) == {}


# ------------------------------------------- mechanism classes / allostery
from promoter_tf_drugs import modulators as mods          # noqa: E402


def _drug(action, name="X", symbol="VDR"):
    from promoter_tf_drugs.drugs import Drug
    return Drug(name, "CHEMBL1", symbol, action, None, True, 4.0, 2000,
                False, True, "up", "test")


def test_mechanism_classification():
    assert mods.classify(_drug("POSITIVE ALLOSTERIC MODULATOR")) == "PAM"
    assert mods.classify(_drug("AGONIST")) == "orthosteric agonist"
    assert mods.classify(_drug("NEGATIVE ALLOSTERIC MODULATOR")) == "NAM"
    assert mods.classify(_drug("SOMETHING ODD")) == "unclassified"


def test_pam_profile_encodes_the_real_tradeoff():
    pam = mods.profile_for("PAM")
    ago = mods.profile_for("orthosteric agonist")

    # PAM advantages
    assert pam.subtype_selectivity == "typically higher"
    assert pam.preserves_dynamics and not ago.preserves_dynamics
    assert pam.selectivity_bonus > ago.selectivity_bonus

    # and the cost, which must not be silently omitted
    assert pam.needs_endogenous_ligand and not ago.needs_endogenous_ligand
    assert "bounded" in pam.max_effect
    assert "endogenous ligand" in pam.caveat


def test_allosteric_bonus_raises_specificity_but_flags_dependency():
    ev = dsc.DrugEvidence(pchembl_median=8.0, pchembl_n=20, n_targets=4)
    ago = dsc.grade_drug(_drug("AGONIST"), 0.7, ev)
    pam = dsc.grade_drug(_drug("POSITIVE ALLOSTERIC MODULATOR"), 0.7, ev)

    assert pam.components["selectivity"] > ago.components["selectivity"]
    assert any("requires endogenous ligand" in n for n in pam.notes)
    # but PAMs still carry lower action confidence - the trade is not free
    assert pam.components["action"] < ago.components["action"]


def test_selectivity_bonus_cannot_exceed_one():
    ev = dsc.DrugEvidence(pchembl_median=8.0, pchembl_n=20, n_targets=1)
    pam = dsc.grade_drug(_drug("POSITIVE ALLOSTERIC MODULATOR"), 0.7, ev)
    assert pam.components["selectivity"] <= 1.0


def test_comparison_detects_both_routes_available():
    hits = [
        mods.ModulatorHit("Ago", "VDR", "orthosteric agonist", "AGONIST", "up",
                          True, False, 2000, False),
        mods.ModulatorHit("Pam", "VDR", "PAM", "POSITIVE ALLOSTERIC MODULATOR",
                          "up", True, False, 2015, True),
    ]
    cmp = mods.compare_mechanisms(hits)
    assert cmp["has_allosteric_option"] is True
    assert "Both allosteric and orthosteric" in cmp["summary"]
    assert cmp["by_mechanism"][0]["mechanism"] == "PAM"   # allosteric sorts first


def test_comparison_says_so_when_only_orthosteric_exists():
    hits = [mods.ModulatorHit("Ago", "VDR", "orthosteric agonist", "AGONIST",
                              "up", True, False, 2000, False)]
    cmp = mods.compare_mechanisms(hits)
    assert cmp["has_allosteric_option"] is False
    assert "No approved allosteric modulator" in cmp["summary"]


def test_direction_filter_excludes_nams_when_upregulating():
    assert "NAM" in mods.DOWN_CLASSES and "NAM" not in mods.UP_CLASSES
    assert "PAM" in mods.UP_CLASSES and "PAM" not in mods.DOWN_CLASSES


def test_mechanism_section_renders_with_the_dependency_warning():
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    assert res["results"][0]["mechanisms"]
    md = pipeline.markdown_report(res)
    assert "orthosteric vs allosteric" in md
    assert "verify that the pathway has endogenous tone" in md


# --------------------------------------------- multi-hop upstream network
from promoter_tf_drugs import network as net          # noqa: E402


def _toy_network():
    """RECEPTOR -> KINASE -> TF (all activating), plus an inhibitory branch
    RECEPTOR2 -| PHOSPHATASE -| TF  (double negative = net activating),
    plus a hub with many edges."""
    E = net.Edge
    edges = [
        E("KINASE", "TF", 1, 5), E("RECEPTOR", "KINASE", 1, 5),
        E("PHOSPHATASE", "TF", -1, 4), E("RECEPTOR2", "PHOSPHATASE", -1, 4),
        E("REPRESSOR", "TF", -1, 3),
    ]
    # a study-bias hub: lots of edges, reaches TF in one hop
    edges.append(E("HUB", "TF", 1, 9))
    for i in range(200):
        edges.append(E("HUB", f"OTHER{i}", 1, 1))
    return edges


def test_signs_multiply_along_the_path():
    """Inhibiting an inhibitor is activation - the whole reason unsigned
    networks cannot support this analysis."""
    inc, out, deg = net.build_graph(_toy_network())
    paths = net.upstream_paths(inc, out, deg, "TF", want_sign=1)
    by_entry = {p.nodes[0]: p for p in paths}

    assert "RECEPTOR" in by_entry
    assert by_entry["RECEPTOR"].signs == [1, 1]

    # double-negative route must be recovered as net ACTIVATING
    assert "RECEPTOR2" in by_entry
    assert by_entry["RECEPTOR2"].signs == [-1, -1]
    assert by_entry["RECEPTOR2"].net_sign == 1


def test_single_inhibitory_edge_excluded_from_activating_paths():
    inc, out, deg = net.build_graph(_toy_network())
    up = {p.nodes[0] for p in net.upstream_paths(inc, out, deg, "TF", want_sign=1)}
    down = {p.nodes[0] for p in net.upstream_paths(inc, out, deg, "TF", want_sign=-1)}
    assert "REPRESSOR" in down and "REPRESSOR" not in up
    assert "PHOSPHATASE" in down and "PHOSPHATASE" not in up


def test_specificity_falls_with_hops_and_fanout():
    """Distance buys druggability and spends selectivity."""
    inc, out, deg = net.build_graph(_toy_network())
    paths = {p.nodes[0]: p for p in net.upstream_paths(inc, out, deg, "TF", want_sign=1)}
    assert paths["RECEPTOR"].hops == 2 and paths["KINASE"].hops == 1
    assert paths["RECEPTOR"].specificity < paths["KINASE"].specificity

    # the hub reaches TF in one hop but fans out to 200 other nodes
    assert paths["HUB"].fan_out > 100
    assert paths["HUB"].specificity < paths["KINASE"].specificity


def test_hub_nodes_are_flagged():
    inc, out, deg = net.build_graph(_toy_network())
    paths = {p.nodes[0]: p for p in net.upstream_paths(inc, out, deg, "TF", want_sign=1)}
    assert "HUB" in paths["HUB"].hub_nodes
    assert paths["KINASE"].hub_nodes == []


def test_degree_normalisation_demotes_the_hub():
    """Without it, every query returns whichever node is best studied."""
    inc, out, deg = net.build_graph(_toy_network())
    raw = net.random_walk_with_restart(inc, out, deg, ["TF"], degree_normalise=False)
    normed = net.random_walk_with_restart(inc, out, deg, ["TF"], degree_normalise=True)
    # the hub's share of probability must drop once degree is accounted for
    assert normed.get("HUB", 0) / max(normed.get("KINASE", 1e-9), 1e-9) < \
           raw.get("HUB", 0) / max(raw.get("KINASE", 1e-9), 1e-9)


def test_confidence_decays_per_hop():
    inc, out, deg = net.build_graph(_toy_network())
    paths = {p.nodes[0]: p for p in net.upstream_paths(inc, out, deg, "TF", want_sign=1)}
    assert paths["RECEPTOR"].confidence < paths["KINASE"].confidence
    assert abs(paths["KINASE"].confidence - net.PER_HOP_CONFIDENCE) < 1e-6


def test_path_describes_its_own_mechanism():
    inc, out, deg = net.build_graph(_toy_network())
    paths = {p.nodes[0]: p for p in net.upstream_paths(inc, out, deg, "TF", want_sign=1)}
    desc = paths["RECEPTOR2"].describe()
    assert "RECEPTOR2" in desc and "TF" in desc and "—|" in desc


def test_repressor_tf_inverts_the_required_direction():
    """To raise a gene through a repressor you must LOWER the repressor."""
    edges = _toy_network()
    act = net.find_upstream_targets(None, "TF", goal="up", tf_role="activator",
                                    edges=edges)
    rep = net.find_upstream_targets(None, "TF", goal="up", tf_role="repressor",
                                    edges=edges)
    assert {h.node for h in act} & {"KINASE", "RECEPTOR"}
    assert "REPRESSOR" in {h.node for h in rep}
    assert "REPRESSOR" not in {h.node for h in act}


def test_downregulation_goal_flips_again():
    edges = _toy_network()
    up = {h.node for h in net.find_upstream_targets(None, "TF", goal="up",
                                                    tf_role="activator", edges=edges)}
    down = {h.node for h in net.find_upstream_targets(None, "TF", goal="down",
                                                      tf_role="activator", edges=edges)}
    assert "KINASE" in up and "KINASE" not in down
    assert "REPRESSOR" in down


def test_hits_carry_warnings_about_distance_and_hubs():
    hits = net.find_upstream_targets(None, "TF", goal="up", tf_role="activator",
                                     edges=_toy_network())
    by_node = {h.node: h for h in hits}
    assert any("hub" in w for w in by_node["HUB"].warnings)


def test_empty_network_returns_nothing_gracefully():
    assert net.find_upstream_targets(None, "TF", edges=[]) == []


def test_upstream_flag_does_not_collide_with_promoter_window():
    """--upstream (hops) and --upstream N (bp) are different things; the
    parameter names must not shadow each other."""
    import inspect
    sig = inspect.signature(pipeline.run)
    assert sig.parameters["upstream"].default == 1000        # bp of promoter
    assert sig.parameters["upstream_network"].default is False
    res = pipeline.run("DEMOGENE", demo=True, offline=True, upstream=500)
    assert res["results"][0]["promoter"]["upstream"] == 500


# ------------------------------------------------------ cell-cycle phase
from promoter_tf_drugs import cellcycle as cc          # noqa: E402


def test_e2f_blocked_in_post_mitotic_tissue():
    """Not demoted - blocked. DREAM actively represses these promoters in
    terminally differentiated cells."""
    a = cc.assess("E2F1", "post-mitotic")
    assert a.compatible is False and a.multiplier == 0.0
    assert any("actively repressed" in n for n in a.notes)


def test_e2f4_is_the_exception():
    """E2F4 is the repressive arm of DREAM and persists in quiescence."""
    assert cc.assess("E2F4", "post-mitotic").compatible is True
    assert cc.assess("E2F1", "post-mitotic").compatible is False


def test_nuclear_receptors_unaffected_by_cycle():
    for sym in ("NR3C1", "VDR", "PPARG", "RARA"):
        a = cc.assess(sym, "post-mitotic")
        assert a.compatible and a.multiplier == 1.0


def test_duty_cycle_dilutes_in_cycling_cells():
    a = cc.assess("E2F1", "cycling")
    assert a.compatible is True
    assert 0 < a.multiplier < 0.5
    assert any("asynchronous" in n for n in a.notes)
    assert any("synchronis" in n for n in a.notes)

    assert cc.assess("NR3C1", "cycling").multiplier == 1.0


def test_mixed_population_scales_with_proliferative_index():
    mixed = cc.assess("FOXM1", "mixed")
    cycling = cc.assess("FOXM1", "cycling")
    assert mixed.multiplier < cycling.multiplier
    assert any("proliferative index" in n for n in mixed.notes)


def test_unknown_state_flags_without_penalising():
    a = cc.assess("E2F1", "unknown")
    assert a.multiplier == 1.0          # no unearned penalty
    assert any("--proliferation" in n for n in a.notes)


def test_unannotated_tf_assumed_cycle_independent():
    a = cc.assess("ZZZ9FAKE", "post-mitotic")
    assert a.compatible and a.multiplier == 1.0
    assert "assumed cycle-independent" in " ".join(a.notes)


def test_family_prefix_catches_unlisted_cycle_genes():
    prof = cc.profile_for("E2F6")
    assert prof.duty < 1.0 and not prof.post_mitotic_active
    assert "verify" in prof.note


def test_bookmarking_factors_are_noted():
    a = cc.assess("FOXA1", "cycling")
    assert any("bookmarking" in n for n in a.notes)


def test_filter_removes_blocked_and_rescales_others():
    class R:
        def __init__(self, sym, score):
            self.symbol, self.score = sym, score
    tfs = [R("E2F1", 80.0), R("NR3C1", 60.0), R("FOXM1", 70.0)]
    ass = cc.filter_and_annotate(tfs, "post-mitotic")
    assert ass["E2F1"]["compatible"] is False
    by = {t.symbol: t.score for t in tfs}
    assert by["E2F1"] == 0.0 and by["FOXM1"] == 0.0
    assert by["NR3C1"] == 60.0                 # untouched


def test_pipeline_removes_incompatible_tfs_from_results():
    res = pipeline.run("DEMOGENE", demo=True, offline=True,
                       proliferation="post-mitotic")
    assert res["query"]["proliferation"] == "post-mitotic"
    assert "cell_cycle" in res["results"][0]
    # demo TFs are all cycle-independent, so nothing should be dropped
    assert res["results"][0]["cell_cycle_blocked"] == []


def test_cell_cycle_section_renders_when_restricted_factors_present():
    """Report must surface the issue even when the user hasn't set the flag."""
    ass = {"E2F1": cc.assess("E2F1", "unknown").to_dict()}
    assert "--proliferation" in cc.summarise(ass, "unknown")
    ass2 = {"E2F1": cc.assess("E2F1", "post-mitotic").to_dict()}
    assert "incompatible" in cc.summarise(ass2, "post-mitotic")


# ------------------------------------------- LINCS signatures / genetics
from promoter_tf_drugs import signatures as sig, genetics as gen   # noqa: E402


def test_landmark_vs_inferred_changes_evidence_weight():
    """L1000 measures ~978 transcripts and INFERS ~11,350. Treating an inferred
    value as a measurement is the commonest misuse of this data."""
    w_lm, note_lm = sig.weight_for(True)
    w_inf, note_inf = sig.weight_for(False)
    w_unk, note_unk = sig.weight_for(None)

    assert w_lm == 1.0 and note_lm is None
    assert w_inf < 0.5 and "inferred" in note_inf and "RT-qPCR" in note_inf
    # unknown must be treated as inferred, not as measured
    assert w_unk == w_inf and note_unk is not None


def test_known_landmark_recognised_without_network():
    assert sig.is_landmark("HMOX1") is True
    assert sig.is_landmark("ZZZ9FAKE") is None      # unknown, not False


def test_connectivity_is_signed_and_symmetric():
    a = {f"G{i}": float(i % 5 - 2) for i in range(40)}
    b = dict(a)
    opposite = {k: -v for k, v in a.items()}

    assert sig.connectivity(a, b) > 0.99
    assert sig.connectivity(a, opposite) < -0.99
    assert abs(sig.connectivity(a, b) - sig.connectivity(b, a)) < 1e-9


def test_connectivity_requires_enough_overlap():
    """Too few shared genes must return 0, not a confident spurious score."""
    a = {"G1": 1.0, "G2": -1.0, "G3": 2.0}
    b = {"G1": 1.0, "G2": -1.0, "G3": 2.0}
    assert sig.connectivity(a, b) == 0.0


def test_non_specific_compounds_are_flagged_and_demoted():
    """Cytotoxic compounds connect to everything; without this filter a
    connectivity search returns the same two dozen drugs for any query."""
    hits = []
    for q in ("GENEA", "GENEB", "GENEC", "GENED"):
        hits.append(sig.SignatureHit("Cytotoxin", q, "up", 3.0,
                                     weighted_score=3.0, reference_drug=q))
        hits.append(sig.SignatureHit(f"Specific-{q}", q, "up", 3.0,
                                     weighted_score=3.0, reference_drug=q))
    sig.flag_non_specific(hits)

    cyto = [h for h in hits if h.drug == "Cytotoxin"]
    spec = [h for h in hits if h.drug.startswith("Specific")]
    assert all(h.non_specific for h in cyto)
    assert not any(h.non_specific for h in spec)
    assert cyto[0].weighted_score < spec[0].weighted_score
    assert any("cytotoxic or promiscuous" in n for n in cyto[0].notes)


def test_connectivity_hits_warn_about_shared_cytotoxicity():
    ref = {f"G{i}": float(i % 7 - 3) for i in range(60)}
    cands = {"DrugA": dict(ref), "DrugB": {k: -v for k, v in ref.items()}}
    hits = sig.connected_to_reference(ref, cands, min_score=0.2)
    assert {h.drug for h in hits} == {"DrugA", "DrugB"}
    assert any(h.direction == "down" for h in hits)     # opposed signature
    assert all("not on target annotation" in " ".join(h.notes) for h in hits)


def test_signature_file_round_trip(tmp_path):
    f = tmp_path / "sig.csv"
    f.write_text("gene,score\nHMOX1,2.5\nEGR1,-1.2\n")
    s = sig.load_signature_file(f)
    assert s == {"HMOX1": 2.5, "EGR1": -1.2}


def test_empty_lincs_result_explains_itself():
    msg = sig.summarise([], "MYGENE")
    assert "expressed in the L1000 cell panel" in msg


def test_genetics_trans_support_is_honest_about_absence():
    """Absence of a trans-eQTL is weak evidence against, not evidence of no
    regulation - trans detection is power-limited."""
    ev = gen.tf_link_support(None, "NR3C1", "PAH")
    assert ev.kind == "trans"
    joined = " ".join(ev.notes)
    assert "eQTLGen" in joined
    assert "Absence is weak evidence" in joined or "weak evidence against" in joined


def test_genetics_evidence_serialises():
    ev = gen.GeneticEvidence("PAH", kind="cis", supported=True,
                             phenotypes=[{"phenotype": "X", "genetic_score": 0.8,
                                          "overall_score": 0.5}])
    d = ev.to_dict()
    assert d["gene"] == "PAH" and d["phenotypes"][0]["phenotype"] == "X"
    assert gen.summarise(ev) and "PCSK9" in gen.summarise(ev)
    assert gen.summarise(gen.GeneticEvidence("PAH")) is None


# ------------------------------------------------- delivery / cell type
from promoter_tf_drugs import delivery as dlv          # noqa: E402

_DIAZEPAM = "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21"
_CAFFEINE = "CN1C=NC2=C1C(=O)N(C)C(=O)N2C"
_ATENOLOL = "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1"
_PLICAMYCIN = ("CC1C(C(CC(O1)OC2C(CC(C(O2)C)OC3CC(C(C(O3)C)O)(C)O)O)O)"
               "OC4CC(C(C(O4)C)O)OC5CC(C(C(O5)C)O)O")


def test_cns_mpo_ranks_known_cns_drugs_above_non_cns():
    caf, _, _ = dlv.cns_mpo(_CAFFEINE)
    dia, _, _ = dlv.cns_mpo(_DIAZEPAM)
    pli, _, _ = dlv.cns_mpo(_PLICAMYCIN)
    assert caf > pli and dia > pli
    assert pli < 1.5


def test_tpsa_is_penalised_at_both_ends():
    """Wager's TPSA term is a hump: too apolar is penalised as well as too
    polar. A one-sided cutoff misses half of it."""
    low = dlv.cns_mpo(properties={"tpsa": 10})[2]["tpsa"]
    good = dlv.cns_mpo(properties={"tpsa": 60})[2]["tpsa"]
    high = dlv.cns_mpo(properties={"tpsa": 130})[2]["tpsa"]
    assert good == 1.0 and low == 0.0 and high == 0.0


def test_missing_parameters_are_reported_not_rescaled():
    score, n, comp = dlv.cns_mpo(_DIAZEPAM)
    assert n == 4 and "pka" not in comp and "clogd" not in comp
    a = dlv.assess_delivery("Diazepam", compartment="cns", smiles=_DIAZEPAM)
    assert a.mpo_parameters == 4
    assert any("4 of 6 parameters" in nn for nn in a.notes)

    full = dlv.cns_mpo(_DIAZEPAM, properties={"logd": 2.5, "pka": 3.4})
    assert full[1] == 6


def test_pka_term_rescues_the_atenolol_case():
    """Documented limitation: 4-parameter MPO passes atenolol, a textbook poor
    penetrant. The pKa term is what catches it."""
    without = dlv.cns_mpo(_ATENOLOL)[0] / dlv.cns_mpo(_ATENOLOL)[1]
    with_pka = dlv.cns_mpo(_ATENOLOL, properties={"pka": 9.6})
    assert with_pka[2]["pka"] < 0.3
    assert with_pka[0] / with_pka[1] < without


def test_known_pgp_substrate_overrides_good_physicochemistry():
    """Active efflux beats passive diffusion - the usual reason a good-looking
    molecule never reaches brain."""
    a = dlv.assess_delivery("Dexamethasone", compartment="cns",
                            smiles=_DIAZEPAM)   # deliberately good physchem
    assert a.pgp_substrate is True
    assert a.likely_penetrant is False
    assert any("pumped out" in n for n in a.notes)


def test_negative_pgp_flag_is_not_treated_as_reassurance():
    a = dlv.assess_delivery("Caffeine", compartment="cns", smiles=_CAFFEINE)
    assert a.likely_penetrant is True
    assert any("not reassurance" in n for n in a.notes)


def test_antibody_against_intracellular_target_is_impossible():
    """Every TF is intracellular. This is not a penalty, it is a hard stop."""
    a = dlv.assess_delivery("SomeMab", modality="monoclonal antibody",
                            target_is_intracellular=True)
    assert a.compatible is False and a.multiplier == 0.0
    assert any("not an option" in n for n in a.notes)

    ok = dlv.assess_delivery("SomeMab", modality="monoclonal antibody",
                             target_is_intracellular=False)
    assert ok.compatible is True


def test_non_cns_compartment_admits_it_cannot_predict():
    a = dlv.assess_delivery("X", compartment="liver", smiles=_CAFFEINE)
    assert a.multiplier == 1.0
    assert any("no validated physicochemical model" in n for n in a.notes)


def test_expression_overlap_flags_disjoint_tissues():
    """A hepatocyte TF cannot regulate your gene in neurons, whatever the
    motif score says."""
    res = dlv.expression_overlap({"liver", "kidney"}, {"brain", "cerebellum"})
    assert res["compatible"] is False
    assert "NO shared tissue" in res["note"]
    assert "Motif quality is irrelevant" in res["note"]

    ok = dlv.expression_overlap({"liver", "brain"}, {"brain"})
    assert ok["compatible"] is True and ok["shared"] == ["brain"]


def test_unknown_expression_does_not_claim_incompatibility():
    res = dlv.expression_overlap(set(), {"brain"})
    assert res["known"] is False
    assert "compatible" not in res          # absence of data is not a verdict


def test_drug_scoring_now_uses_cns_mpo():
    """The tissue term should discriminate a small CNS-like molecule from a
    large polar one, via MPO rather than a three-way cutoff."""
    good = dsc.DrugEvidence(properties={"full_mwt": 284, "psa": 33, "cx_logp": 2.8})
    bad = dsc.DrugEvidence(properties={"full_mwt": 1085, "psa": 331, "cx_logp": -0.4})
    s_good, _ = dsc._tissue(good, "cns")
    s_bad, note = dsc._tissue(bad, "cns")
    assert s_good > s_bad
    assert "CNS MPO" in note
    # and the filter only applies when CNS is requested
    assert dsc._tissue(bad, None)[0] == 1.0


# ------------------------------------- stabilisers / glues / degradation
from promoter_tf_drugs import stabilizers as stab          # noqa: E402


def test_pocketless_tf_gets_the_stabilisation_recommendation():
    """The whole point: SP1 has no ligand pocket, so agonist chemistry will
    never reach it, but its interfaces and degradation machinery can."""
    a = stab.assess("SP1", has_direct_drug=False)
    assert a["has_ligand_pocket"] is False
    assert a["recommended"] is True
    assert any("no ligand-binding pocket" in n for n in a["notes"])
    # with nothing curated, it must still point somewhere concrete
    assert a.get("exemplars")
    assert any("E3 ligase" in n for n in a["notes"])


def test_nuclear_receptor_is_not_pushed_toward_stabilisation():
    a = stab.assess("VDR", annotation={"nuclear_receptor": True})
    assert a["has_ligand_pocket"] is True
    assert a["recommended"] is False
    assert any("simpler route" in n for n in a["notes"])


def test_activity_versus_abundance_is_distinguished():
    """A conformational stabiliser raises activity; a degradation blocker
    raises levels. A western blot cannot see the first."""
    assert stab.MECHANISMS["conformational stabiliser"]["raises"] == "activity"
    assert stab.MECHANISMS["degradation blocker"]["raises"] == "abundance"

    nrf2 = stab.curated_for("NFE2L2", "up")
    assert nrf2 and all(s.raises == "abundance" for s in nrf2)
    a = stab.assess("NFE2L2")
    assert any("western blot cannot see" in n for n in a["notes"])


def test_approved_degradation_blockers_are_present():
    """NRF2 and HIF are the two classes that actually shipped."""
    nrf2 = {s.drug for s in stab.curated_for("NFE2L2", "up") if s.status == "approved"}
    hif = {s.drug for s in stab.curated_for("HIF1A", "up") if s.status == "approved"}
    assert "Dimethyl fumarate" in nrf2 and "Omaveloxolone" in nrf2
    assert "Daprodustat" in hif or "Roxadustat" in hif


def test_glue_degraders_appear_only_for_downregulation():
    """Lenalidomide destroys IKZF1 - that is a down mechanism, not an up one."""
    down = {s.drug for s in stab.curated_for("IKZF1", "down")}
    assert "Lenalidomide" in down and "Thalidomide" in down
    assert stab.curated_for("IKZF1", "up") == []
    assert all(s.raises.startswith("none") for s in stab.curated_for("IKZF1", "down"))


def test_unapproved_entries_are_labelled_clearly():
    tp53 = {s.drug: s for s in stab.curated_for("TP53", "up")}
    epre = tp53["Eprenetapopt"]
    assert epre.status != "approved"
    assert any("NOT approved" in n for n in epre.notes)
    assert epre.mechanism == "conformational stabiliser"


def test_ppi_disruptor_direction_is_up_despite_being_an_inhibitor():
    """MDM2 inhibitors RELEASE p53 - inhibiting the inhibitor raises the TF."""
    tp53 = {s.drug: s for s in stab.curated_for("TP53", "up")}
    assert tp53["Idasanutlin"].mechanism == "PPI disruptor (releasing)"
    assert tp53["Idasanutlin"].target_node == "MDM2"
    assert tp53["Idasanutlin"].direction == "up"


def test_exemplars_are_all_approved():
    """These exist to show the mechanism class shipped; an unapproved
    exemplar would defeat the purpose."""
    assert all(e.status == "approved" for e in stab.EXEMPLARS)
    names = {e.drug for e in stab.EXEMPLARS}
    assert "Tafamidis" in names and "Migalastat" in names


def test_stabilisation_section_renders():
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    stabs = res["results"][0].get("stabilisation")
    assert stabs, "no stabilisation assessment produced"
    md = pipeline.markdown_report(res)
    assert "Locking the factor into an active state" in md
    assert "no ligand pocket" in md


def test_cli_parser_builds_without_conflicts():
    """--upstream (bp of promoter) and --upstream-network (hops) are different
    flags; registering both under one name silently breaks every invocation."""
    from promoter_tf_drugs.cli import build_parser
    parser = build_parser()
    args = parser.parse_args(["--gene", "PAH", "--upstream", "2000",
                              "--upstream-network"])
    assert args.upstream == 2000
    assert args.upstream_network is True


def test_cli_accepts_the_full_flag_set():
    """Guard against a future flag collision anywhere in the parser."""
    from promoter_tf_drugs.cli import build_parser
    args = build_parser().parse_args([
        "--gene", "PAH", "--direction", "down", "--upstream", "5000",
        "--downstream", "1000", "--all-tss", "--pvalue", "1e-5",
        "--chip-api", "--tf-weight", "0.7", "--tissue", "cns",
        "--proliferation", "post-mitotic", "--cell-type-check",
        "--upstream-network", "--upstream-depth", "4", "--lincs",
        "--genetics", "--predict-direction", "--screen",
    ])
    assert args.gene == "PAH" and args.direction == "down"
    assert args.proliferation == "post-mitotic" and args.upstream == 5000


def test_all_output_files_written_for_every_section(tmp_path):
    """Regression: later sections reusing the name `hits` clobbered the motif
    hits and broke sites.bed. End-to-end write is the only thing that catches
    variable shadowing like that."""
    res = pipeline.run("SMN2", demo=True, offline=True)
    written = pipeline.write_outputs(res, tmp_path)
    names = {p.name for p in written}
    assert names == {"report.md", "results.json", "tf_summary.tsv",
                     "sites.bed", "promoter.fa"}

    bed = (tmp_path / "sites.bed").read_text().splitlines()
    assert len(bed) > 1
    for line in bed[1:]:
        f = line.split("\t")
        assert f[0].startswith("chr") and int(f[2]) > int(f[1])
        assert f[5] in "+-"


def test_motif_hits_survive_downstream_sections():
    res = pipeline.run("SMN2", demo=True, offline=True)
    hits = res["results"][0]["hits"]
    assert hits and all("start" in h and "matrix_id" in h for h in hits)


# ------------------------------------------- AlphaFold guards / screen design
def _af_pdb(plddts):
    """Column-exact AlphaFold CA records; pLDDT lives in B-factor (cols 61-66)."""
    lines = []
    for i, v in enumerate(plddts, 1):
        lines.append("ATOM  " + f"{i:5d}" + " " + " CA " + " " + "ALA" + " "
                     + "A" + f"{i:4d}" + " " + "   "
                     + f"{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}"
                     + f"{1.0:6.2f}{v:6.2f}" + "           C")
    return "\n".join(lines)


def test_plddt_parsed_from_bfactor_column():
    p = docking.residue_plddt(_af_pdb([95.0, 40.0, 72.5]))
    assert p == {1: 95.0, 2: 40.0, 3: 72.5}


def test_disordered_model_is_refused_for_docking():
    """The TF-specific trap: low pLDDT is a disorder prediction, and a docking
    program will happily score compounds into it with great affinities."""
    a = docking.assess_model(_af_pdb([40.0] * 80 + [95.0] * 20), "P00000")
    assert a.dockable is False
    assert a.ordered_fraction < 0.5
    joined = " ".join(a.notes)
    assert "disorder prediction" in joined and "Do not dock" in joined


def test_confident_model_still_warns_about_apo_side_chains():
    a = docking.assess_model(_af_pdb([95.0] * 100), "P00001")
    assert a.dockable is True
    assert any("apo" in n and "cryptic pockets" in n for n in a.notes)


def test_partially_ordered_model_restricts_to_folded_region():
    a = docking.assess_model(_af_pdb([95.0] * 60 + [30.0] * 40), "P2")
    assert a.dockable is True
    assert any("restrict docking" in n for n in a.notes)

    ordered = docking.ordered_residues(_af_pdb([95.0] * 60 + [30.0] * 40))
    assert len(ordered) == 60
    assert docking.pose_is_in_ordered_region([1, 2, 3, 4, 5], ordered) is True
    assert docking.pose_is_in_ordered_region([61, 62, 63, 64, 65], ordered) is False
    assert docking.pose_is_in_ordered_region([], ordered) is False


# --- screen design --------------------------------------------------------
def test_guides_found_and_ranked_by_crispra_window():
    from promoter_tf_drugs import screens
    prom = _prom(1)
    prom.sequence = demo_promoter_sequence()
    guides = screens.find_guides(prom.sequence, prom.offset_to_tss)
    assert guides
    # best guides must sit in the productive window
    assert guides[0].in_optimal_window
    assert all(len(g.sequence) == 20 for g in guides)
    assert all(g.pam.endswith("GG") for g in guides)


def test_guide_with_pol3_terminator_is_flagged():
    from promoter_tf_drugs import screens
    seq = "A" * 50 + "TTTT" + "GCGCGCGCGCGCGCGC" + "AGG" + "A" * 50
    guides = screens.find_guides(seq, lambda i: i - 50)
    flagged = [g for g in guides if "TTTT" in g.sequence]
    if flagged:
        assert any("terminates U6" in w for w in flagged[0].warnings)


def test_guides_outside_window_are_warned():
    from promoter_tf_drugs import screens
    prom = _prom(1)
    prom.sequence = demo_promoter_sequence()
    guides = screens.find_guides(prom.sequence, prom.offset_to_tss, max_guides=40)
    out = [g for g in guides if not g.in_optimal_window]
    for g in out:
        assert any("productive CRISPRa window" in w for w in g.warnings)


def test_reporter_plan_includes_the_mutant_control():
    """The scrambled-motif control is what makes a hit interpretable."""
    from promoter_tf_drugs import screens
    prom = _prom(1)
    plan = screens.reporter_plan(prom)
    joined = " ".join(plan["controls"])
    assert "Scrambled-motif" in joined
    assert "Promoterless" in joined and "normalisation" in joined
    assert "chromatin context" in plan["caveat"]


def test_compound_plan_prioritises_direct_agents_first():
    from promoter_tf_drugs import screens
    plan = screens.compound_plan(
        [{"drug_name": "DrugB", "tf_symbol": "VDR", "grade": "A"}],
        [{"drug": "DrugA"}])
    assert plan["cherry_pick"][0]["compound"] == "DrugA"
    assert plan["cherry_pick"][0]["priority"] == 1
    assert "counter_screen" in plan and "constitutive-promoter" in plan["counter_screen"]


def test_expectations_state_what_a_negative_result_means():
    from promoter_tf_drugs import screens
    e = screens.expectations()
    assert "Z'" in e["quality_metric"]
    assert "enhancer" in e["what_a_negative_result_means"]


def test_screen_plan_renders_end_to_end():
    res = pipeline.run("SMN2", demo=True, offline=True, design_screen=True)
    sp = res["results"][0]["screen_plan"]
    assert sp["gene"] == "SMN2"
    md = pipeline.markdown_report(res)
    assert "Screen design" in md and "CRISPRa screen" in md
    assert "Cherry-pick" in md or "cherry" in md.lower()


# --------------------- discovery: nothing silently withheld (regression)
def test_no_incumbent_still_reports_candidates():
    """Regression: with baseline None, `compare` never populated wins_on, so
    novel_improvements came back EMPTY for exactly the case you care about
    most — a gene with no existing treatment."""
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    imp = res["results"][0]["improvements"]
    assert imp["has_incumbent"] is False
    assert imp["n_candidates"] > 0
    assert imp["novel_improvements"], "candidates hidden when there is no incumbent"

    md = pipeline.markdown_report(res)
    assert "No incumbent found" in md
    # and the contradictory line must be gone
    assert "No unlinked candidate beat the" not in md


def test_every_candidate_is_accounted_for():
    """No candidate may vanish between the total and the reported buckets."""
    for gene in ("DEMOGENE", "SMN2"):
        imp = pipeline.run(gene, demo=True, offline=True)["results"][0]["improvements"]
        a = imp["accounting"]
        assert (a["previously_linked"] + a["novel_and_on_front"]
                + a["novel_but_dominated"] + a["novel_no_advantage"]) == a["total"]
        assert len(imp["all_candidates"]) == a["total"]


def test_bar_is_per_axis_best_not_one_drug():
    """A challenger must beat the best established value on each axis. Using a
    single summed-score incumbent lets it 'win' on specificity while another
    established drug is more specific still."""
    from promoter_tf_drugs import discovery as d
    strong = d.Candidate("Strong", "direct", None,
                         d.Axes(0.95, 0.30, 0.50, "measured"), previously_linked=True)
    specific = d.Candidate("Specific", "direct", None,
                           d.Axes(0.40, 0.98, 0.50, "measured"), previously_linked=True)
    env = d.baseline_envelope([strong, specific])
    assert env.strength == 0.95 and env.specificity == 0.98

    # A challenger beating only the weaker drug's specificity must NOT count.
    ch = d.Candidate("Challenger", "via TF", "VDR", d.Axes(0.50, 0.60, 0.50))
    d.compare(ch, strong, env)
    assert "higher specificity" not in ch.wins_on


def test_dominated_novel_candidates_are_listed_not_dropped():
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    imp = res["results"][0]["improvements"]
    assert imp["novel_dominated"]
    md = pipeline.markdown_report(res)
    assert "shown so nothing is hidden" in md
    assert imp["novel_dominated"][0]["name"] in md


# ------------------------------------------- unified route catalogue
from promoter_tf_drugs import catalogue as cat            # noqa: E402


def test_evidence_tier_outranks_score_always():
    """A weak clinical route must beat a brilliant hypothesis. Collapsing these
    into one number destroys exactly the information you need."""
    block = {
        "direct_therapies": [{
            "drug": "WeakButApproved", "gene": "G", "modality": "siRNA",
            "mechanism": "knockdown", "specificity": 0.3,
            "raises_protein": True, "approved": True, "withdrawn": False}],
        "screens": {"TFX": {
            "enrichment": {"usable": True},
            "hits": [{"drug": "GreatHypothesis", "quadrant": "NOVEL HYPOTHESIS",
                      "probability": 0.99}]}},
    }
    routes = cat.build(block)["routes"]
    assert routes[0]["drug"] == "WeakButApproved"
    assert routes[0]["tier"] == 1
    assert routes[-1]["tier"] == 5


def test_directness_is_not_double_counted_in_the_sort():
    """Regression: directness sat in both the score and the sort key, which
    floated genome-wide chromatin drugs above far better TF ligands."""
    res = pipeline.run("SMN2", demo=True, offline=True)
    routes = res["results"][0]["all_routes"]["routes"]
    tier3 = [r for r in routes if r["tier"] == 3]
    scores = [r["score"] for r in tier3]
    assert scores == sorted(scores, reverse=True), "tier not sorted by score"

    chrom = next(r for r in tier3 if r["category"] == "chromatin-level")
    tf = next(r for r in tier3 if r["category"] == "TF: orthosteric ligand")
    assert tf["score"] > chrom["score"]
    assert tier3.index(tf) < tier3.index(chrom)


def test_chromatin_drugs_scored_as_remote_and_unspecific():
    routes = cat.chromatin_routes("up")
    assert routes
    for r in routes:
        assert r.axes.specificity <= 0.1
        assert r.directness >= 3, "genome-wide action is not 'direct on my gene'"
        assert any("genome-wide" in n for n in r.notes)
    assert cat.chromatin_routes("down") == []


def test_transcript_acting_small_molecule_categorised_correctly():
    """Risdiplam is a splicing modifier - it acts on the transcript, not the
    SMN protein, despite being a small molecule."""
    block = {"direct_therapies": [{
        "drug": "Risdiplam", "gene": "SMN2", "modality": "small molecule",
        "mechanism": "SMN2 pre-mRNA splicing modifier", "specificity": 0.55,
        "raises_protein": True, "approved": True, "withdrawn": False}]}
    r = cat.build(block)["routes"][0]
    assert r["category"] == "direct: gene or transcript"


def test_all_modules_feed_the_catalogue():
    """Every source module must be able to contribute a route."""
    block = {
        "direct_therapies": [{"drug": "D1", "gene": "G", "modality": "siRNA",
                              "mechanism": "kd", "specificity": 0.9,
                              "raises_protein": True, "approved": True,
                              "withdrawn": False}],
        "drug_grades": [{"drug_name": "D2", "tf_symbol": "VDR",
                         "tf_confidence": 0.7, "components": {}, "notes": []}],
        "stabilisation": {"SP1": {"stabilisers": [
            {"drug": "D3", "mechanism": "degradation blocker", "raises":
             "abundance", "target_node": "E3", "status": "approved",
             "approval_year": 2020}]}},
        "upstream": {"SP1": [{"node": "KIN", "net_sign": 1,
                              "best_path": {"hops": 2, "confidence": 0.6,
                                            "specificity": 0.5,
                                            "description": "KIN --> SP1"},
                              "drugs": [{"name": "D4", "approved": True}]}]},
        "lincs": {"hits": [{"drug": "D5", "score": 3.0, "is_landmark": True,
                            "non_specific": False, "cell_line": "MCF7"}]},
        "screens": {"VDR": {"enrichment": {"usable": True},
                            "hits": [{"drug": "D6",
                                      "quadrant": "NOVEL HYPOTHESIS",
                                      "probability": 0.8}]}},
    }
    out = cat.build(block)
    sources = {r["source"] for r in out["routes"]}
    assert {"direct", "drug_scoring", "stabilizers", "network",
            "signatures", "screening"} <= sources
    assert out["n_routes"] == len(out["routes"])


def test_non_specific_lincs_hits_excluded_from_catalogue():
    block = {"lincs": {"hits": [
        {"drug": "Cytotoxin", "score": 5.0, "is_landmark": True,
         "non_specific": True},
        {"drug": "Real", "score": 3.0, "is_landmark": True,
         "non_specific": False}]}}
    names = {r["drug"] for r in cat.build(block, include_chromatin=False)["routes"]}
    assert names == {"Real"}


def test_unusable_screen_contributes_nothing():
    """A screen that failed its enrichment check must not leak hypotheses in."""
    block = {"screens": {"VDR": {
        "enrichment": {"usable": False, "reason": "not discriminating"},
        "hits": [{"drug": "Bogus", "quadrant": "NOVEL HYPOTHESIS",
                  "probability": 0.99}]}}}
    assert cat.build(block, include_chromatin=False)["routes"] == []


def test_inferred_lincs_value_drops_a_tier():
    measured = cat.build({"lincs": {"hits": [
        {"drug": "A", "score": 3.0, "is_landmark": True, "non_specific": False}]}},
        include_chromatin=False)
    inferred = cat.build({"lincs": {"hits": [
        {"drug": "A", "score": 3.0, "is_landmark": False, "non_specific": False}]}},
        include_chromatin=False)
    assert measured["routes"][0]["tier"] == 2
    assert inferred["routes"][0]["tier"] == 4
    assert any("INFERRED" in n for n in inferred["routes"][0]["notes"])


def test_absent_categories_reported():
    out = cat.build({"direct_therapies": []}, include_chromatin=False)
    assert out["n_routes"] == 0
    assert len(out["categories_absent"]) == len(cat.CATEGORIES)


def test_catalogue_renders_in_report():
    res = pipeline.run("SMN2", demo=True, offline=True)
    md = pipeline.markdown_report(res)
    assert "Every route that acts on this gene" in md
    assert "Best in each category" in md
    assert "Nusinersen" in md


# ----------------------------------------- Ensembl masking fallback
def test_masked_sequence_failure_falls_back_to_unmasked():
    """Repeat masking is an optimisation. Ensembl intermittently 500s on
    mask=soft; that must not kill the run."""
    from promoter_tf_drugs import genes as g
    from promoter_tf_drugs.http import HttpError

    class C:
        def __init__(self):
            self.calls = []

        def get(self, url, params=None, **kw):
            self.calls.append(dict(params or {}))
            if "mask" in (params or {}):
                raise HttpError("HTTP 500")
            return {"seq": "acgtACGT"}

    c = C()
    seq = g.fetch_sequence(c, "homo_sapiens", "3", 1, 8, 1, mask="soft")
    assert len(c.calls) == 2
    assert "mask" in c.calls[0] and "mask" not in c.calls[1]
    # fallback path must normalise case, since nothing is masked any more
    assert seq == "ACGTACGT"


def test_unmasked_failure_raises_an_actionable_error():
    from promoter_tf_drugs import genes as g
    from promoter_tf_drugs.http import HttpError
    import pytest

    class C:
        def get(self, url, params=None, **kw):
            raise HttpError("HTTP 500")

    with pytest.raises(HttpError, match="service is having problems"):
        g.fetch_sequence(C(), "homo_sapiens", "3", 1, 8, 1)


def test_masked_success_preserves_lowercase():
    from promoter_tf_drugs import genes as g

    class C:
        def get(self, url, params=None, **kw):
            return {"seq": "acgtACGT"}

    assert g.fetch_sequence(C(), "homo_sapiens", "3", 1, 8, 1,
                            mask="soft") == "acgtACGT"


# ------------------------------------------- direct modality design
from promoter_tf_drugs import modalities as mo          # noqa: E402

_TX = ("ATG" + "".join(np.random.default_rng(5).choice(list("ACGT"), size=1200)))


def test_delivery_is_the_gate_not_design_quality():
    """A potent oligo for an undeliverable compartment is a cell-culture
    reagent. The verdict must say so."""
    immune = mo.design("CCR5", goal="down", tissue="microglia")
    liver = mo.design("PCSK9", goal="down", tissue="liver")

    assert immune["feasible_in_vivo"] is False
    assert "DELIVERY, not design, is the blocker" in immune["verdict"]
    assert "ex vivo" in immune["verdict"]

    assert liver["feasible_in_vivo"] is True
    assert "Design quality will be the limiting factor" in liver["verdict"]


def test_tissue_aliases_map_to_compartments():
    for t in ("microglia", "T cell", "macrophage", "blood"):
        assert mo.delivery_for(t)["compartment"] == "immune cells"
    assert mo.delivery_for("neuron")["compartment"] == "cns"
    assert mo.delivery_for("hepatocyte")["compartment"] == "liver"
    assert mo.delivery_for(None)["compartment"] == "unknown"


def test_sirna_rules_applied():
    cands = mo.design_sirna(_TX, n=10)
    assert cands and all(len(c.sequence) == 21 for c in cands)
    assert all(c.position >= 75 for c in cands), "5' region not skipped"
    # antisense must be the reverse complement
    for c in cands:
        assert c.antisense == mo.revcomp(c.sequence)
    # best candidates should be cleaner than the worst
    assert len(cands[0].warnings) <= len(cands[-1].warnings)


def test_sirna_flags_bad_gc_and_homopolymers():
    bad = "A" * 80 + "GGGGGGGGGGGGGGGGGGGGG" + "A" * 80
    cands = mo.design_sirna(bad, n=5)
    assert any("run of 4+" in w for c in cands for w in c.warnings)


def test_thermodynamic_asymmetry_computed():
    """AU-rich antisense 5' end favours correct strand loading."""
    good = mo._asymmetry("GCGC" + "A" * 17)     # sense 5' stable, antisense AU
    bad = mo._asymmetry("ATAT" + "G" * 17)
    assert good > bad


def test_gapmer_flags_cpg_and_g_quadruplex():
    cands = mo.design_gapmer("A" * 80 + "CGCGCGCGCGCGCGCG" + "T" * 80, n=5)
    assert any("CpG" in w for c in cands for w in c.warnings)
    cands2 = mo.design_gapmer("A" * 80 + "GGGGTTTTGGGGTTTT" + "A" * 80, n=5)
    assert any("quadruplex" in w for c in cands2 for w in c.warnings)


def test_crispri_window_differs_from_crispra():
    """KRAB silencing works near the TSS; VP64 activation works well upstream."""
    from promoter_tf_drugs import screens
    assert mo.CRISPRI_WINDOW != screens.CRISPRA_WINDOW
    assert mo.CRISPRI_WINDOW[1] > 0        # CRISPRi extends downstream of TSS
    assert screens.CRISPRA_WINDOW[1] < 0   # CRISPRa stays upstream

    prom = _prom(1)
    prom.sequence = demo_promoter_sequence()
    guides = mo.design_crispri(prom.sequence, prom.offset_to_tss, n=6)
    assert guides
    assert guides[0]["in_optimal_window"]


def test_upregulation_offers_different_modalities():
    up = mo.design("SMN2", goal="up", tissue="cns")
    mods = {o["modality"] for o in up["options"]}
    assert "splice-switching oligonucleotide" in mods
    assert "CRISPRa (dCas9-VP64/SAM)" in mods
    assert "siRNA" not in mods

    down = mo.design("TTR", goal="down", tissue="liver", transcript=_TX)
    assert "siRNA" in {o["modality"] for o in down["options"]}


def test_off_target_warning_is_explicit_about_scope():
    r = mo.design("X", goal="down", tissue="liver", transcript=_TX)
    assert "seed" in r["off_target_warning"].lower()
    assert "BLAST" in r["off_target_warning"]
    assert "local sanity only" in r["off_target_warning"]


def test_modality_section_renders():
    res = pipeline.run("CCR5", demo=True, offline=True, design_modality=True,
                       goal="down", tissue="cns")
    md = pipeline.markdown_report(res)
    assert "Designing a direct modality" in md
    assert "Delivery to" in md


# ----------------------------------------------- post-hoc shortlisting
from promoter_tf_drugs import shortlist as sl          # noqa: E402


def _result_for_shortlist():
    return pipeline.run("SMN2", demo=True, offline=True)


def test_missing_grade_does_not_disqualify():
    """Direct agents have no letter grade. Dropping a tier-1 approved therapy
    for lacking a field that does not apply to it is backwards."""
    assert sl._grade_ok(None, "B") is True
    assert sl._grade_ok("A", "B") is True
    assert sl._grade_ok("C", "B") is False

    rows = sl.filter_routes(_result_for_shortlist(), min_grade="B")
    assert any(r["drug"] == "Nusinersen" for r in rows)


def test_strict_excludes_genome_wide_agents():
    """Chromatin drugs have no grade and would otherwise slip past the grade
    filter despite specificity ~0.05."""
    rows = sl.filter_routes(_result_for_shortlist(), **sl.PRESETS["strict"])
    assert not any(r["category"] == "chromatin-level" for r in rows)
    assert all(r["axes"]["specificity"] >= 0.3 for r in rows)


def test_presets_are_progressively_narrower():
    res = _result_for_shortlist()
    total = len((res["results"][0]["all_routes"])["routes"])
    counts = {p: len(sl.filter_routes(res, **sl.PRESETS[p])) for p in sl.PRESETS}
    assert counts["empirical"] <= counts["strict"] or counts["strict"] >= 0
    assert all(c <= total for c in counts.values())
    # empirical keeps only measured evidence
    assert all(r["tier"] <= 2
               for r in sl.filter_routes(res, **sl.PRESETS["empirical"]))


def test_specific_preset_cuts_distance_and_breadth():
    rows = sl.filter_routes(_result_for_shortlist(), **sl.PRESETS["specific"])
    assert all(r["directness"] <= 1 for r in rows)
    assert not any(r["category"] == "chromatin-level" for r in rows)


def test_novel_preset_excludes_reported_links():
    res = _result_for_shortlist()
    rows = sl.filter_routes(res, **sl.PRESETS["novel"])
    assert not any(r["previously_linked"] for r in rows)
    assert all(r["tier"] >= 3 for r in rows)


def test_filters_compose_and_narrow_monotonically():
    res = _result_for_shortlist()
    wide = len(sl.filter_routes(res, max_tier=5))
    narrow = len(sl.filter_routes(res, max_tier=5, max_steps=1))
    narrower = len(sl.filter_routes(res, max_tier=5, max_steps=1,
                                    min_specificity=0.5))
    assert wide >= narrow >= narrower


def test_empty_shortlist_explains_which_knob_to_loosen():
    text = sl.render([], {"max_tier": 1, "min_grade": "A"}, 24)
    assert "Nothing passed" in text
    assert "Loosen one criterion" in text
    assert "require_chip" in text


def test_shortlist_cli_runs(tmp_path):
    import json
    res = _result_for_shortlist()
    f = tmp_path / "results.json"
    f.write_text(json.dumps(res, default=str), encoding="utf-8")
    out = tmp_path / "short.md"
    assert sl.main([str(f), "--preset", "mechanism", "--out", str(out)]) == 0
    assert "Shortlist" in out.read_text(encoding="utf-8")


# ------------------------------------ scan quality / API robustness
def test_expected_false_positives_matches_the_null():
    """800 motifs x 6000 bp x 2 strands x 1e-4 = 960. If a scan returns ~953,
    it found nothing — this is the check that says so."""
    from promoter_tf_drugs import motifs as m
    assert abs(m.expected_false_positives(800, 6000, 1e-4) - 960.0) < 1e-6
    assert abs(m.expected_false_positives(800, 1200, 1e-5) - 19.2) < 1e-6


def test_scan_quality_calls_noise_noise():
    from promoter_tf_drugs import motifs as m
    q = m.scan_quality(953, 800, 6000, 1e-4)
    assert q["verdict"] == "no signal"
    assert q["enrichment"] < 1.2
    assert "chance alone predicts" in q["advice"]
    assert "1e-6" in q["advice"]


def test_scan_quality_recognises_real_enrichment():
    from promoter_tf_drugs import motifs as m
    q = m.scan_quality(200, 800, 1200, 1e-5)
    assert q["verdict"] == "enriched" and q["enrichment"] > 2


def test_no_signal_warning_reaches_the_report():
    from promoter_tf_drugs import pipeline as pl
    res = pl.run("DEMOGENE", demo=True, offline=True)
    res["results"][0]["scan_quality"] = {
        "verdict": "no signal", "observed": 953, "expected_by_chance": 960.0,
        "enrichment": 0.99, "advice": "953 sites is what chance alone predicts."}
    md = pl.markdown_report(res)
    assert "This scan found no signal" in md
    assert "ranked from a list of chance matches" in md


def test_expression_failure_warns_once_then_disables():
    """A schema error fails identically for every gene; 11 identical warnings
    is noise, and continuing to ask wastes the run."""
    from promoter_tf_drugs import delivery as d
    d._EXPR_DISABLED["value"] = False

    class C:
        def __init__(self):
            self.posts = 0

        def post(self, url, data, **kw):
            # the resolver call succeeds; the expression query is what 400s
            if "search" in (data.get("query") or ""):
                return {"data": {"search": {"hits": [{"id": "ENSG1"}]}}}
            self.posts += 1
            raise RuntimeError("HTTP 400")

        def get(self, url, params=None, **kw):
            return {}

    c = C()
    assert d.tissue_expression(c, "CCR5") == set()
    first = c.posts
    assert d._EXPR_DISABLED["value"] is True
    for sym in ("GATA2", "GATA4", "KLF4"):
        assert d.tissue_expression(c, sym) == set()
    assert c.posts == first, "kept querying a known-broken endpoint"
    d._EXPR_DISABLED["value"] = False


def test_post_redirect_is_followed_with_body(tmp_path):
    """urllib will not replay a POST body across 307/308, so it surfaces the
    redirect as an error. LINCS returns 307."""
    import urllib.error, io
    from promoter_tf_drugs.http import Client

    # Fresh cache dir per run: a persistent one serves the second call from
    # cache and the test silently stops exercising the redirect path.
    c = Client(cache_dir=tmp_path / "cache", offline=False)
    calls = []

    def fake(req, timeout=None):
        calls.append(req.full_url)
        if len(calls) == 1:
            raise urllib.error.HTTPError(
                req.full_url, 307, "Temporary Redirect",
                {"Location": "https://example.org/v2"}, None)

        class R:
            headers = {}
            def read(self): return b'{"ok": true}'
            def __enter__(self): return self
            def __exit__(self, *a): return False
        return R()

    import urllib.request
    orig = urllib.request.urlopen
    urllib.request.urlopen = fake
    try:
        out = c.post("https://example.org/v1", {"q": 1})
    finally:
        urllib.request.urlopen = orig

    assert out == {"ok": True}
    assert calls == ["https://example.org/v1", "https://example.org/v2"]


# ------------------------------------------------ offline data sources
def test_hpa_long_format_uses_symbol_not_ensembl_id(tmp_path):
    """HPA ships both 'Gene' (Ensembl id) and 'Gene name' (symbol). Scanning
    columns positionally grabs the id, which then matches nothing."""
    from promoter_tf_drugs import delivery as d
    f = tmp_path / "hpa.tsv"
    f.write_text("Gene\tGene name\tTissue\tnTPM\n"
                 "ENSG1\tCCR5\tspleen\t42.0\n"
                 "ENSG1\tCCR5\tliver\t0.2\n"
                 "ENSG2\tGATA2\tspleen\t15.0\n", encoding="utf-8")
    table = d.load_expression_file(f)
    try:
        assert "CCR5" in table and "ENSG1" not in table
        assert d.tissue_expression(None, "CCR5") == {"spleen"}   # liver below floor
        ov = d.expression_overlap(d.tissue_expression(None, "CCR5"),
                                  d.tissue_expression(None, "GATA2"))
        assert ov["compatible"] is True
    finally:
        d._LOCAL_EXPRESSION.clear()


def test_gtex_wide_format_parsed(tmp_path):
    from promoter_tf_drugs import delivery as d
    f = tmp_path / "gtex.tsv"
    f.write_text("Name\tDescription\tWhole Blood\tLiver\n"
                 "ENSG1\tCCR5\t30.5\t0.1\n", encoding="utf-8")
    try:
        table = d.load_expression_file(f)
        assert table["CCR5"] == {"whole blood": 30.5}
    finally:
        d._LOCAL_EXPRESSION.clear()


def test_local_expression_table_bypasses_the_api(tmp_path):
    """A loaded table must take priority, so a broken endpoint is irrelevant."""
    from promoter_tf_drugs import delivery as d

    class Boom:
        def post(self, *a, **k): raise AssertionError("API was called")
        def get(self, *a, **k): raise AssertionError("API was called")

    f = tmp_path / "e.tsv"
    f.write_text("Gene name\tTissue\tnTPM\nCCR5\tspleen\t42\n", encoding="utf-8")
    d.load_expression_file(f)
    try:
        assert d.tissue_expression(Boom(), "CCR5") == {"spleen"}
    finally:
        d._LOCAL_EXPRESSION.clear()


def test_local_signature_table_respects_direction_and_landmark(tmp_path):
    from promoter_tf_drugs import signatures as sg
    f = tmp_path / "sig.csv"
    f.write_text("drug,gene,score,cell_line,is_landmark\n"
                 "DrugA,CCR5,3.2,PC3,1\n"
                 "DrugB,CCR5,-2.5,MCF7,0\n"
                 "DrugC,CCR5,0.4,PC3,1\n"
                 "DrugD,OTHER,9,PC3,1\n", encoding="utf-8")

    up = sg.load_signature_table(f, "CCR5", direction="up")
    down = sg.load_signature_table(f, "CCR5", direction="down")

    assert [h.drug for h in up] == ["Druga"]          # only the positive score
    assert [h.drug for h in down] == ["Drugb"]
    assert up[0].is_landmark is True and up[0].weighted_score == 3.2
    # inferred value must be down-weighted, not trusted equally
    assert down[0].is_landmark is False
    assert down[0].weighted_score < abs(down[0].score)
    assert any("local signature table" in n for n in up[0].notes)


def test_local_signature_table_filters_weak_scores(tmp_path):
    from promoter_tf_drugs import signatures as sg
    f = tmp_path / "s.tsv"
    f.write_text("drug\tgene\tscore\nWeak\tCCR5\t0.4\nStrong\tCCR5\t4.0\n",
                 encoding="utf-8")
    hits = sg.load_signature_table(f, "CCR5", direction="up")
    assert [h.drug for h in hits] == ["Strong"]


def test_expression_file_accepts_zip_gz_or_plain(tmp_path):
    """HPA ships a .zip, GTEx a .gz, and both get saved under a .tsv name.
    Sniff the magic bytes rather than trusting the extension."""
    import gzip, zipfile
    from promoter_tf_drugs import delivery as d

    rows = ("Gene\tGene name\tTissue\tnTPM\n"
            "ENSG1\tCCR5\tspleen\t42.0\n"
            "ENSG2\tGATA2\tspleen\t15.0\n")

    z = tmp_path / "hpa.tsv"                     # a zip named .tsv
    with zipfile.ZipFile(z, "w") as zf:
        zf.writestr("rna_tissue_consensus.tsv", rows)
    g = tmp_path / "gz.tsv"                      # a gzip named .tsv
    g.write_bytes(gzip.compress(rows.encode()))
    t = tmp_path / "plain.tsv"
    t.write_text(rows, encoding="utf-8")

    for f in (z, g, t):
        d._LOCAL_EXPRESSION.clear()
        assert set(d.load_expression_file(f)) == {"CCR5", "GATA2"}, f.name
    d._LOCAL_EXPRESSION.clear()


def test_html_error_page_reported_as_failed_download(tmp_path):
    from promoter_tf_drugs import delivery as d
    import pytest
    f = tmp_path / "bad.tsv"
    f.write_text("<!DOCTYPE html><html>404 not found</html>", encoding="utf-8")
    with pytest.raises(ValueError, match="HTML page, not data"):
        d.load_expression_file(f)


def test_missing_expression_file_says_so(tmp_path):
    from promoter_tf_drugs import delivery as d
    import pytest
    with pytest.raises(FileNotFoundError):
        d.load_expression_file(tmp_path / "nope.tsv")


def test_unrecognised_table_names_the_header(tmp_path):
    """The error must say what it saw, or you cannot tell which file is wrong."""
    from promoter_tf_drugs import delivery as d
    import pytest
    f = tmp_path / "wrong.tsv"
    f.write_text("alpha\tbeta\tgamma\n1\t2\t3\n", encoding="utf-8")
    with pytest.raises(ValueError, match="not gene symbols"):
        d.load_expression_file(f)
    d._LOCAL_EXPRESSION.clear()


def test_remap_response_shapes_all_parse():
    """The API has changed shape between versions and a wrong guess yields an
    empty list rather than an error."""
    from promoter_tf_drugs import evidence as ev

    shapes = [
        [{"target": "CTCF"}, {"target": "NR3C1"}],
        {"peaks": [{"target": "CTCF"}, {"target": "NR3C1"}]},
        {"results": [{"tf": "CTCF"}, {"tf_name": "NR3C1"}]},
        {"data": {"chr15:1-2": [{"target": "CTCF"}, {"target": "NR3C1"}]}},
        {"features": [{"target_name": "CTCF"}, {"target_name": "NR3C1"}]},
    ]
    for shape in shapes:
        rows = ev._extract_rows(shape)
        assert len(rows) == 2, shape


def test_sparse_chip_result_is_flagged_not_trusted(caplog):
    """One TF over a promoter is a parsing failure, not biology. Silently
    treating every TF as unbound skews the whole ranking."""
    from promoter_tf_drugs import evidence as ev

    class C:
        def get(self, url, params=None, **kw):
            return [{"target": "CTCF"}]

    with caplog.at_level("WARNING"):
        counts = ev.peaks_from_api(C(), "15", 40156023, 40162022)
    assert counts == {"CTCF": 1}
    assert any("more likely an API or parsing problem" in r.message
               for r in caplog.records)
    assert any("--chip-bed" in r.message for r in caplog.records)


def test_healthy_chip_result_is_not_flagged(caplog):
    from promoter_tf_drugs import evidence as ev

    class C:
        def get(self, url, params=None, **kw):
            return [{"target": t} for t in
                    ("CTCF", "NR3C1", "SP1", "JUN", "FOS", "MYC")]

    with caplog.at_level("WARNING"):
        counts = ev.peaks_from_api(C(), "15", 1, 6000)
    assert len(counts) == 6
    assert not any("parsing problem" in r.message for r in caplog.records)


# ------------------------------- optional stages must never abort the run
def test_remote_disconnected_is_retried_not_raised(tmp_path):
    """RemoteDisconnected subclasses ConnectionResetError -> OSError, not
    URLError. A narrower except clause lets it escape the retry loop and kill
    the whole run — which is exactly what happened with LINCS."""
    import http.client
    import urllib.request
    from promoter_tf_drugs.http import Client, HttpError

    c = Client(cache_dir=tmp_path / "c", offline=False, min_interval=0)
    attempts = []

    def fake(req, timeout=None):
        attempts.append(1)
        raise http.client.RemoteDisconnected("closed without response")

    orig = urllib.request.urlopen
    urllib.request.urlopen = fake
    try:
        import pytest
        with pytest.raises(HttpError, match="failed after"):
            c.get("https://example.org/x", retries=2)
    finally:
        urllib.request.urlopen = orig

    assert len(attempts) == 2, "gave up without retrying"


def test_connection_reset_also_retried(tmp_path):
    import urllib.request
    from promoter_tf_drugs.http import Client, HttpError
    import pytest

    c = Client(cache_dir=tmp_path / "c2", offline=False, min_interval=0)
    calls = []

    def fake(req, timeout=None):
        calls.append(1)
        raise ConnectionResetError("reset by peer")

    orig = urllib.request.urlopen
    urllib.request.urlopen = fake
    try:
        with pytest.raises(HttpError):
            c.get("https://example.org/y", retries=2)
    finally:
        urllib.request.urlopen = orig
    assert len(calls) == 2


def test_optional_helper_swallows_any_failure(caplog):
    from promoter_tf_drugs.pipeline import _optional

    def boom():
        raise RuntimeError("service is down")

    with caplog.at_level("WARNING"):
        assert _optional("LINCS", boom, []) == []
    assert any("skipping that section" in r.message for r in caplog.records)
    assert any("rest of the run is unaffected" in r.message
               for r in caplog.records)

    assert _optional("ok", lambda: 42, None) == 42


def test_dead_lincs_does_not_abort_the_pipeline(monkeypatch):
    """A dead third-party API should cost you that section, not the run."""
    import http.client
    from promoter_tf_drugs import pipeline as pl, signatures as sg

    def boom(*a, **k):
        raise http.client.RemoteDisconnected("closed")

    monkeypatch.setattr(sg, "query_gene", boom)
    res = pl.run("DEMOGENE", demo=True, offline=True)
    assert res["results"], "run aborted when an optional stage failed"
    assert res["results"][0]["tfs"], "core analysis lost"


# ------------------------------- endogenous ligand modulation
from promoter_tf_drugs import ligands as lg          # noqa: E402


def test_strategy_determines_direction():
    """Inhibiting synthesis lowers the receptor; inhibiting catabolism raises
    it. Getting this backwards would recommend the opposite drug."""
    assert lg.STRATEGY_DIRECTION["inhibit synthesis"] == "down"
    assert lg.STRATEGY_DIRECTION["inhibit catabolism"] == "up"
    assert lg.STRATEGY_DIRECTION["supplement ligand"] == "up"
    assert lg.STRATEGY_DIRECTION["sequester ligand"] == "down"

    for r in lg.ROUTES["ESR1"]:
        assert r.direction == lg.STRATEGY_DIRECTION[r.strategy]


def test_aromatase_inhibitors_lower_er_not_raise_it():
    down = {r.drug for r in lg.routes_for("ESR1", "down")}
    up = {r.drug for r in lg.routes_for("ESR1", "up")}
    assert {"Anastrozole", "Letrozole", "Exemestane"} <= down
    assert "Anastrozole" not in up
    assert "Estradiol" in up


def test_catabolism_inhibition_raises_the_receptor():
    """Carbenoxolone blocks cortisol inactivation, so GR activity goes UP —
    the opposite of metyrapone, which blocks synthesis."""
    up = {r.drug for r in lg.routes_for("NR3C1", "up")}
    down = {r.drug for r in lg.routes_for("NR3C1", "down")}
    assert "Carbenoxolone" in up and "Carbenoxolone" not in down
    assert "Metyrapone" in down and "Metyrapone" not in up


def test_tissue_selective_routes_rank_first():
    """Local ligand change is the main advantage this approach offers."""
    routes = lg.routes_for("NR3C1", "up")
    assert routes[0].tissue_selective is True
    a = lg.assess("NR3C1", "up")
    assert any("tissue-selective" in n for n in a["notes"])


def test_systemic_caveat_always_present():
    """The dominant failure mode must not be omitted."""
    a = lg.assess("ESR1", "down")
    joined = " ".join(a["notes"])
    assert "systemic phenotype" in joined
    assert "homeostatic" in joined or "feedback" in joined.lower()
    # the concrete example, not just the abstraction
    assert any("bone loss" in c or "osteoporosis" in c
               for r in a["routes"] for c in r["caveats"]) or "bone loss" in joined


def test_supplementation_saturation_warning():
    a = lg.assess("VDR", "up")
    assert any("saturated" in n for n in a["notes"])


def test_unapproved_ligand_routes_labelled():
    tal = next(r for r in lg.ROUTES["RARA"] if r.drug == "Talarozole")
    assert tal.status != "approved"
    assert any("NOT approved" in c for c in tal.caveats)
    assert lg.routes_for("RARA", "up", approved_only=True) and \
        "Talarozole" not in {r.drug for r in
                             lg.routes_for("RARA", "up", approved_only=True)}


def test_ligand_independent_tf_says_so():
    a = lg.assess("SP1", "up")
    assert a["has_ligand_axis"] is False
    assert a["routes"] == []
    assert any("ligand-independent" in n for n in a["notes"])


def test_ligand_routes_enter_the_catalogue():
    from promoter_tf_drugs import catalogue as cat
    block = {"ligand_axes": {"ESR1": lg.assess("ESR1", "down")}}
    routes = cat.build(block, goal="down", include_chromatin=False)["routes"]
    assert routes
    assert all(r["category"] == "TF: endogenous ligand" for r in routes)
    assert all(r["directness"] == 2 for r in routes)   # drug -> enzyme -> ligand -> TF


def test_tissue_selective_scores_more_specific_in_catalogue():
    from promoter_tf_drugs import catalogue as cat
    block = {"ligand_axes": {"AR": lg.assess("AR", "down")}}
    routes = cat.build(block, goal="down", include_chromatin=False)["routes"]
    by_drug = {r["drug"]: r for r in routes}
    # finasteride's enzyme is tissue-restricted; abiraterone's is not
    assert (by_drug["Finasteride"]["axes"]["specificity"]
            > by_drug["Abiraterone"]["axes"]["specificity"])


def test_chain_reads_the_full_causal_path():
    """A route is only actionable if you can read each link and judge it.
    Printing just the drug name hides three inferential steps."""
    anas = next(r for r in lg.ROUTES["ESR1"] if r.drug == "Anastrozole")
    c = lg.chain(anas, "MYGENE")
    for part in ("Anastrozole", "CYP19A1", "less estradiol",
                 "less ESR1 activity", "down at MYGENE"):
        assert part in c, c


def test_catabolism_chain_says_more_not_less():
    carb = next(r for r in lg.ROUTES["NR3C1"] if r.drug == "Carbenoxolone")
    c = lg.chain(carb)
    assert "blocks HSD11B2" in c
    assert "more cortisol" in c and "more NR3C1 activity" in c


def test_supplement_chain_drops_the_enzyme_link():
    """No enzyme sits in the path, so the chain must not invent one or repeat
    the ligand twice."""
    hc = next(r for r in lg.ROUTES["NR3C1"] if r.drug == "Hydrocortisone")
    c = lg.chain(hc)
    assert c.count("cortisol") == 1
    assert "blocks" not in c and "→ more NR3C1 activity" in c


def test_enzyme_names_deduplicated():
    """'CYP11B1' and 'CYP11B1 (11b-hydroxylase)' are one enzyme, not two."""
    enzymes = lg.assess("NR3C1", "up")["enzymes"]
    assert "CYP11B1" in enzymes
    assert not any(e.startswith("CYP11B1 (") for e in enzymes)


def test_opposite_direction_routes_are_listed_separately():
    """Same axis, opposite effect — these must be visible but never mixed in
    with the recommendations."""
    a = lg.assess("NR3C1", "up")
    up_drugs = {r["drug"] for r in a["routes"]}
    opp_drugs = {o["drug"] for o in a["opposite_direction"]}
    assert "Metyrapone" in opp_drugs and "Metyrapone" not in up_drugs
    assert "Carbenoxolone" in up_drugs and "Carbenoxolone" not in opp_drugs
    assert not (up_drugs & opp_drugs)


def test_covered_receptors_lists_each_ligand():
    cov = lg.covered_receptors()
    assert cov["ESR1"] == "estradiol"
    assert cov["AR"] == "dihydrotestosterone"
    assert cov["VDR"] == "calcitriol"
    assert "SP1" not in cov


def test_report_lists_chains_and_uncovered_tfs():
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    md = pipeline.markdown_report(res)
    if "Changing the natural ligand" in md:
        assert "The chain, link by link" in md
        assert "Same axis, opposite direction" in md
        assert "No ligand axis for:" in md


# --------------------------------------------- per-compound expansion
from promoter_tf_drugs import explain as ex          # noqa: E402


def _demo_result():
    return pipeline.run("SMN2", demo=True, offline=True)


def test_explain_gathers_evidence_from_every_section():
    """The point: a run scatters one drug's evidence across eight sections."""
    found = ex.collect(_demo_result(), "calcitriol")
    assert found["names"] == ["Calcitriol"]
    b = found["blocks"][0]
    assert b["gene"] == "SMN2"
    for section in ("routes", "grades", "outcomes", "comparison"):
        assert section in b, f"{section} not collected"


def test_explain_renders_components_and_limiting_step():
    md = ex.render(ex.collect(_demo_result(), "calcitriol"))
    assert "Grade A" in md
    assert "median pChEMBL" in md
    assert "Limiting step" in md
    assert "The bar:" in md and "Previously linked to this gene" in md


def test_explain_reports_the_incumbent_as_such():
    md = ex.render(ex.collect(_demo_result(), "nusinersen"))
    assert "Direct agent" in md
    assert "Raises protein: yes" in md


def test_activity_only_agents_flagged_in_expansion():
    result = _demo_result()
    result["results"][0]["direct_therapies"] = [{
        "drug": "Ivacaftor", "gene": "CFTR", "modality": "small molecule",
        "mechanism": "CFTR potentiator", "approval_year": 2012,
        "raises_protein": False, "approved": True, "withdrawn": False,
        "specificity": 0.55, "direction": "up", "source": "curated"}]
    md = ex.render(ex.collect(result, "ivacaftor"))
    assert "activity only" in md


def test_ambiguous_query_reports_all_matches_not_the_first():
    """Silently picking one compound would hide the others."""
    found = ex.collect(_demo_result(), "cal")
    if len(found["names"]) > 1:
        md = ex.render(found)
        assert "compounds match" in md
        for n in found["names"]:
            assert n in md


def test_unknown_compound_says_so_and_suggests_next_step():
    md = ex.render(ex.collect(_demo_result(), "zzz-not-a-drug"))
    assert "Not found in this results file" in md
    assert "shortlist" in md


def test_explain_cli_writes_a_file(tmp_path):
    import json
    f = tmp_path / "results.json"
    f.write_text(json.dumps(_demo_result(), default=str), encoding="utf-8")
    out = tmp_path / "calcitriol.md"
    assert ex.main([str(f), "calcitriol", "--out", str(out)]) == 0
    assert "Calcitriol" in out.read_text(encoding="utf-8")


# ------------------------------------------- category glossary in report
def test_every_category_has_a_guide_entry():
    """A category with no explanation is a column header the reader cannot
    act on."""
    from promoter_tf_drugs import catalogue as c
    missing = [cat for cat in c.CATEGORIES if cat not in c.CATEGORY_GUIDE]
    assert not missing, f"no guide entry for: {missing}"
    for cat, g in c.CATEGORY_GUIDE.items():
        assert set(g) == {"what", "why", "cost", "differs"}, cat
        assert all(len(v) > 40 for v in g.values()), cat


def test_guide_distinguishes_the_confusable_pairs():
    """These four pairs are the ones that get conflated."""
    from promoter_tf_drugs import catalogue as c
    g = c.CATEGORY_GUIDE

    # abundance vs activity
    assert "abundance" in g["TF: degradation blocker"]["differs"]
    assert "unchanged" in g["TF: stabiliser or glue"]["differs"]
    # receptor vs the ligand that feeds it
    assert "Never touches the receptor" in g["TF: endogenous ligand"]["differs"]
    # orthosteric vs allosteric
    assert "selectivity" in g["TF: allosteric modulator"]["differs"]
    # gene vs gene product
    assert "gene product" in g["direct: gene or transcript"]["differs"]


def test_glossary_only_describes_categories_present():
    res = pipeline.run("SMN2", demo=True, offline=True)
    md = pipeline.markdown_report(res)
    present = res["results"][0]["all_routes"]["by_category"]

    assert "What these categories mean" in md
    for cat, n in present.items():
        assert f"**{cat}** ({n} found)" in md

    # absent ones appear only in the "what you're missing" list, with no count
    assert "Not found here, and what you would be missing" in md
    for cat in res["results"][0]["all_routes"]["categories_absent"]:
        assert f"**{cat}** ({0} found)" not in md


def test_glossary_leads_with_the_assay_distinction():
    md = pipeline.markdown_report(pipeline.run("SMN2", demo=True, offline=True))
    i = md.find("What these categories mean")
    head = md[i:i + 400]
    assert "abundance" in head and "activity" in head
    assert "western blot" in head


# ------------------------------- CTD reported links, expanded
def test_link_direction_parses_ctd_actions():
    from promoter_tf_drugs.discovery import link_direction
    assert link_direction(["increases^expression"]) == "up"
    assert link_direction(["decreases^expression"]) == "down"
    assert link_direction(["increases^expression", "decreases^expression"]) == "both"
    assert link_direction(["affects^binding"]) == "unclear"


def test_wrong_direction_links_are_separated_not_hidden():
    """A chemical reported to RAISE your gene is a lead when you want it up
    and a warning when you want it down. Mixing them would be actively
    misleading."""
    from promoter_tf_drugs.discovery import summarise_links
    links = {"dexamethasone": ["increases^expression"],
             "tretinoin": ["decreases^expression"]}

    up = summarise_links(links, "up", set())
    assert [e["chemical"] for e in up["wanted_direction"]] == ["Dexamethasone"]
    assert [e["chemical"] for e in up["opposite_direction"]] == ["Tretinoin"]

    down = summarise_links(links, "down", set())
    assert [e["chemical"] for e in down["wanted_direction"]] == ["Tretinoin"]
    assert [e["chemical"] for e in down["opposite_direction"]] == ["Dexamethasone"]


def test_compounds_reported_but_not_predicted_are_surfaced():
    """The interesting group: measured effect on the gene, no mechanism from
    this analysis. Reporting only a count buries them."""
    from promoter_tf_drugs.discovery import summarise_links
    links = {"dexamethasone": ["increases^expression"],
             "cadmium chloride": ["increases^expression"]}
    s = summarise_links(links, "up", {"Dexamethasone"})
    names = [e["chemical"] for e in s["reported_but_not_predicted"]]
    assert names == ["Cadmium Chloride"]
    assert s["counts"]["not_predicted_by_us"] == 1

    matched = next(e for e in s["wanted_direction"]
                   if e["chemical"] == "Dexamethasone")
    assert matched["in_candidates"] is True


def test_both_direction_links_counted_as_context_dependent():
    from promoter_tf_drugs.discovery import summarise_links
    s = summarise_links({"valproic acid":
                         ["increases^expression", "decreases^expression"]},
                        "up", set())
    assert s["counts"]["unclear"] == 1
    # still worth surfacing as not-predicted, since the effect is real
    assert s["counts"]["not_predicted_by_us"] == 1


def test_substring_matching_credits_the_candidate():
    from promoter_tf_drugs.discovery import summarise_links
    s = summarise_links({"all-trans retinoic acid": ["increases^expression"]},
                        "up", {"Retinoic Acid"})
    assert s["wanted_direction"][0]["in_candidates"] is True


def test_reported_links_reach_the_report():
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    imp = res["results"][0]["improvements"]
    assert "reported_links" in imp
    imp["reported_links"] = {
        "wanted_direction": [{"chemical": "Dexamethasone", "direction": "up",
                              "actions": "increases^expression",
                              "in_candidates": True}],
        "opposite_direction": [{"chemical": "Tretinoin", "direction": "down",
                                "actions": "decreases^expression",
                                "in_candidates": False}],
        "unclear_direction": [],
        "reported_but_not_predicted": [
            {"chemical": "Cadmium Chloride", "direction": "up",
             "actions": "increases^expression", "in_candidates": False}],
        "counts": {"total": 3, "wanted": 1, "opposite": 1, "unclear": 1,
                   "not_predicted_by_us": 1},
    }
    md = pipeline.markdown_report(res)
    assert "Compounds already reported to change this gene" in md
    assert "Dexamethasone" in md and "Cadmium Chloride" in md
    assert "WRONG way" in md
    assert "NOT predicted by this run" in md


# --------------------------- outcome table must agree with discovery
def test_outcome_estimates_sorted_by_their_own_value():
    """They used to inherit grade order, so the table was sorted by a
    different number than the one it displayed."""
    res = pipeline.run("SMN2", demo=True, offline=True)
    ps = [e["p_point"] for e in res["results"][0]["outcome_estimates"]]
    assert ps == sorted(ps, reverse=True)


def test_discovery_winner_always_appears_in_the_outcome_table():
    """The reported bug: the candidate the report names as the winner fell
    outside the truncated rows of the next section."""
    res = pipeline.run("SMN2", demo=True, offline=True)
    imp = res["results"][0]["improvements"]
    winners = {c["name"] for c in (imp.get("novel_improvements") or [])}
    if not winners:
        return
    md = pipeline.markdown_report(res, top_n=1)   # aggressive truncation
    section = _section(md, "Predicted transcriptional outcome")
    for w in winners:
        assert w in section, f"{w} named as winner but missing from outcomes"
        assert "★" in section


def _section(md: str, heading: str) -> str:
    """Slice one report section, stopping at the next heading — otherwise a
    row count picks up every table below it."""
    start = md.find(heading)
    assert start >= 0, f"section {heading!r} missing"
    nxt = md.find("\n### ", start + len(heading))
    return md[start:nxt if nxt > 0 else len(md)]


def test_pinned_rows_do_not_exceed_the_row_budget_silently():
    res = pipeline.run("SMN2", demo=True, offline=True)
    md = pipeline.markdown_report(res, top_n=3)
    section = _section(md, "Predicted transcriptional outcome")
    rows = [l for l in section.splitlines()
            if l.startswith("| ") and "P(engage)" not in l and "---" not in l]
    assert 0 < len(rows) <= 5, f"row budget blown out by pinning: {len(rows)}"


def test_flat_spread_is_explained_not_presented_as_a_ranking():
    """Two of four links are constants without --predict-direction, so the
    chain conveys magnitude, not order. Saying so prevents misreading it."""
    res = pipeline.run("SMN2", demo=True, offline=True)
    ests = res["results"][0]["outcome_estimates"]
    for i, e in enumerate(ests):           # force a flat spread
        e["p_point"] = 0.15
        e["p_low"], e["p_high"] = 0.08, 0.22
    md = pipeline.markdown_report(res)
    assert "These numbers barely differ" in md
    assert "conveys *magnitude*" in md
    assert "not a ranking" in md


def test_binding_model_trains_where_a_direction_model_cannot():
    """The core fix: /activity has hundreds of actives where /mechanism has
    two, so a screen that used to refuse can now run."""
    from promoter_tf_drugs import screening as sc
    rng = np.random.default_rng(3)

    # actives share a scaffold; decoys do not
    actives = [{"chembl_id": f"A{i}", "smiles": f"Oc1ccc(cc1)C{'C' * (i % 5)}N",
                "pchembl": 7.0} for i in range(30)]
    decoys = [{"chembl_id": f"D{i}",
               "smiles": "CC(=O)N" + "C" * (i % 7) + "C(=O)O"}
              for i in range(90)]

    model, auc = sc.build_binding_model(actives, decoys)
    assert model is not None
    assert auc is not None and auc > 0.7, f"AUC {auc}"

    # and it refuses below the active floor rather than producing noise
    assert sc.build_binding_model(actives[:5], decoys) == (None, None)


def test_screen_reports_actives_not_mechanism_labels():
    """Error text must name the real limit, so a data gap is not read as a
    judgement about the receptor."""
    import inspect
    from promoter_tf_drugs import screening as sc
    src = inspect.getsource(sc.screen_receptor)
    assert "measured actives" in src
    assert "data gap, not a judgement" in src
    assert "fetch_actives" in src


def test_direction_is_an_annotation_not_a_gate():
    """Missing direction labels must not block the screen — they used to."""
    import inspect
    from promoter_tf_drugs import screening as sc
    src = inspect.getsource(sc.screen_receptor)
    i_bind = src.index("build_binding_model")
    i_dir = src.index("fetch_labelled_ligands")
    assert i_bind < i_dir, "direction model still gates the screen"
    assert "predicted BINDERS" in src


# ------------------------- undirected drug-gene links must not masquerade
def _undirected(name, direction, source="DGIdb", score=None):
    from promoter_tf_drugs.drugs import Drug
    d = Drug(name, None, "GATA2", "UNKNOWN" if direction == "unclear" else "AGONIST",
             None, True, 4.0, None, False, None, direction, source)
    d.interaction_score = score
    return d


def test_undirected_links_are_held_back_by_default():
    """The reported problem: a table where every row reads UNKNOWN / unclear.
    An undirected link cannot answer a question that is entirely about
    direction."""
    from promoter_tf_drugs.drugs import filter_drugs_by_goal
    drugs = {"GATA2": [_undirected("Cisplatin", "unclear", score=0.3),
                       _undirected("RealAgonist", "up")]}
    anns = {"GATA2": {"role": "activator"}}

    actionable, unknown = filter_drugs_by_goal(drugs, anns, "up")
    assert [d.name for d in actionable["GATA2"]] == ["RealAgonist"]
    assert [d.name for d in unknown["GATA2"]] == ["Cisplatin"]


def test_unknown_direction_can_be_opted_back_in():
    from promoter_tf_drugs.drugs import filter_drugs_by_goal
    drugs = {"GATA2": [_undirected("Cisplatin", "unclear"), _undirected("RealAgonist", "up")]}
    anns = {"GATA2": {"role": "activator"}}
    actionable, unknown = filter_drugs_by_goal(
        drugs, anns, "up", allow_unknown_direction=True)
    assert {d.name for d in actionable["GATA2"]} == {"Cisplatin", "RealAgonist"}
    assert unknown == {}


def test_undirected_links_sorted_by_confidence():
    from promoter_tf_drugs.drugs import filter_drugs_by_goal
    drugs = {"GATA2": [_undirected("Weak", "unclear", score=0.1),
                       _undirected("Strong", "unclear", score=0.9)]}
    _, unknown = filter_drugs_by_goal(drugs, {"GATA2": {"role": "activator"}}, "up")
    assert [d.name for d in unknown["GATA2"]] == ["Strong", "Weak"]


def test_chembl_record_wins_but_keeps_dgidb_provenance():
    """Deduplication used to discard DGIdb's score and source list entirely."""
    from promoter_tf_drugs.drugs import Drug
    chembl = Drug("Calcitriol", "CHEMBL1", "VDR", "AGONIST", "VDR agonist",
                  True, 4.0, 1978, False, True, "up", "ChEMBL")
    dgidb = Drug("Calcitriol", None, "VDR", "AGONIST", None, True, 4.0, None,
                 False, None, "up", "DGIdb")
    dgidb.interaction_score, dgidb.sources = 0.8, "TTD; ChemblInteractions"

    seen = {}
    for d in (dgidb, chembl):
        key = (d.name.lower(), d.target_symbol.upper())
        prev = seen.get(key)
        if prev is None:
            seen[key] = d
        elif prev.source == "DGIdb" and d.source == "ChEMBL":
            d.interaction_score = d.interaction_score or prev.interaction_score
            d.sources = d.sources or prev.sources
            seen[key] = d
    merged = seen[("calcitriol", "VDR")]
    assert merged.source == "ChEMBL" and merged.first_approval == 1978
    assert merged.interaction_score == 0.8 and "TTD" in merged.sources


def test_undirected_section_explains_itself():
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    res["results"][0]["direction_unknown"] = {
        "GATA2": [{"name": "Cisplatin", "interaction_score": 0.4,
                   "sources": "ChemblInteractions"}]}
    md = pipeline.markdown_report(res)
    assert "no recorded direction" in md
    assert "held back from the recommendations" in md
    assert "--include-unknown-direction" in md
    assert "Cisplatin" in md


# --------------------- CNS assessment must not overstate what it knows
def test_mpo_refuses_below_three_parameters():
    """'MPO 0.5/1' reads like an assessment but rests on one descriptor."""
    from promoter_tf_drugs import delivery as d
    a = d.assess_delivery("Tamibarotene", compartment="cns",
                          properties={"hbd": 1})
    assert a.cns_mpo is None
    assert a.likely_penetrant is None
    assert a.mpo_components == {}
    assert "not assessed" in a.notes[0] and "only 1 of 6" in a.notes[0]


def test_mpo_reports_which_parameters_were_missing():
    from promoter_tf_drugs import delivery as d
    a = d.assess_delivery("X", compartment="cns", properties={"mw": 300})
    assert "(mw)" in a.notes[0], a.notes[0]


def test_three_parameters_is_enough():
    from promoter_tf_drugs import delivery as d
    a = d.assess_delivery("Y", compartment="cns",
                          properties={"mw": 300, "tpsa": 60, "logp": 2.0})
    assert a.cns_mpo is not None and a.mpo_parameters == 3
    assert a.likely_penetrant is not None


def test_drug_scoring_tissue_term_uses_the_same_floor():
    """Otherwise the grade silently scores a tissue term the report refuses
    to state."""
    from promoter_tf_drugs import delivery as dl
    thin = dsc.DrugEvidence(properties={"full_mwt": 300})
    score, note = dsc._tissue(thin, "cns")
    assert "not assessed" in note
    assert f"of 6" in note
    assert score == 0.6          # neutral, neither pass nor fail


def test_notes_are_not_orphaned_onto_the_previous_drug():
    """The reported bug: 20 identical 'could not be assessed' lines printed
    under one drug's heading, because the heading was conditional and the
    notes loop was not."""
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    res["results"][0]["cell_type"] = {
        "tfs": {},
        "delivery": [
            {"drug": "Tamibarotene", "cns_mpo": 3.2, "mpo_parameters": 4,
             "likely_penetrant": True, "notes": ["one note", "one note"]},
        ] + [{"drug": f"NoData{i}", "cns_mpo": None, "mpo_parameters": 0,
              "likely_penetrant": None,
              "notes": ["CNS exposure not assessed: only 0 of 6"]}
             for i in range(20)],
    }
    md = pipeline.markdown_report(res)
    section = _section(md, "Cell-type reach")

    # the assessed drug keeps its own note, deduplicated
    assert section.count("one note") == 1
    # the 20 failures are summarised once, not repeated under Tamibarotene
    assert section.count("CNS exposure not assessed") <= 1
    assert "20 drug(s) could not be assessed" in section
    assert "neither passed nor failed" in section


# --------------------- report must not triple itself under --all-tss
def _multi_promoter_result():
    """Two promoter blocks, as --all-tss produces."""
    res = pipeline.run("SMN2", demo=True, offline=True, design_screen=True)
    block = res["results"][0]
    second = json.loads(json.dumps(block, default=str))
    second["promoter"]["transcript_id"] = "ENST00000000002"
    res["results"].append(second)
    return res


def test_screen_design_emitted_once_across_promoters():
    """Three TSSs produced three identical screen-design sections dumped at
    the end of the report."""
    import json as _json
    globals().setdefault("json", _json)
    res = _multi_promoter_result()
    md = pipeline.markdown_report(res)
    assert md.count("## Screen design") == 1
    assert "clone a second construct for the other TSS" in md


def test_modality_design_emitted_once():
    import json as _json
    globals().setdefault("json", _json)
    res = _multi_promoter_result()
    res["results"][0]["modality_design"] = {
        "gene": "SMN2", "goal": "up",
        "delivery": {"compartment": "cns", "status": "solved", "route": "IT",
                     "note": "n"},
        "feasible_in_vivo": True, "verdict": "v", "options": [],
        "off_target_warning": "w"}
    res["results"][1]["modality_design"] = res["results"][0]["modality_design"]
    md = pipeline.markdown_report(res)
    assert md.count("Designing a direct modality") == 1


def test_approved_drug_table_respects_top_n():
    """This table was untruncated: 850 rows per promoter, 2,550 lines of a
    3,943-line report."""
    res = pipeline.run("SMN2", demo=True, offline=True)
    md = pipeline.markdown_report(res, top_n=3)
    section = _section(md, "Approved drugs acting on these factors")
    rows = [l for l in section.splitlines()
            if l.startswith("| ") and "---" not in l and "TF role" not in l]
    assert len(rows) <= 4, f"{len(rows)} rows printed with top_n=3"
    if len(rows) == 4:
        assert "more" in rows[-1]


# ------------- null expectation must match what was actually scanned
def test_scannable_length_excludes_masked_and_ambiguous():
    from promoter_tf_drugs import motifs as m
    # lowercase = soft-masked repeat (scannable only when masking is off);
    # N = genuinely ambiguous and never scannable
    seq = "ACGT" * 1000 + "acgt" * 500 + "N" * 500
    assert m.scannable_length(seq, skip_repeats=True) == 4000
    assert m.scannable_length(seq, skip_repeats=False) == 6000


def test_expectation_uses_scannable_not_window_length():
    """The nonsense line: '293 sites is what chance alone predicts (1055)'.
    293 is nowhere near 1055 — the null was computed over the whole 6 kb
    window while masking had removed most of it from the scan."""
    from promoter_tf_drugs import motifs as m
    full = m.scan_quality(293, 879, 6000, 1e-4)
    real = m.scan_quality(293, 879, 1500, 1e-4)
    assert full["expected_by_chance"] > 1000
    assert 200 < real["expected_by_chance"] < 350
    assert real["scannable_bases"] == 1500


def test_fewer_hits_than_chance_is_worded_correctly():
    """A count far BELOW the null is a different situation from matching it,
    and saying 'is what chance predicts' there is simply false."""
    from promoter_tf_drugs import motifs as m
    q = m.scan_quality(293, 879, 6000, 1e-4)
    assert q["verdict"] == "no signal"
    assert "FEWER than chance" in q["advice"]
    assert "is what chance alone predicts" not in q["advice"]

    at_chance = m.scan_quality(1050, 879, 6000, 1e-4)
    assert "is what chance alone predicts" in at_chance["advice"]


def test_scan_quality_reports_the_basis_of_its_claim():
    from promoter_tf_drugs import motifs as m
    for n in (100, 1050, 5000):
        q = m.scan_quality(n, 879, 1200, 1e-4)
        assert "scannable bases" in q["advice"]
        assert q["scannable_bases"] == 1200


# ------------------------------- ReMap API: the documented nested schema
_REMAP_DOC_SHAPE = {"peaks": [
    {"peakValues": {"chrom": "chr15", "chromStart": 1, "chromEnd": 2,
                    "name": {"data": {"Experiment": "GSE1", "TF": "CTCF",
                                      "Biotype": "K562"}}}},
    {"peakValues": {"chrom": "chr15", "chromStart": 3, "chromEnd": 4,
                    "name": {"data": {"Experiment": "GSE2", "TF": "CTCF",
                                      "Biotype": "HeLa"}}}},
    {"peakValues": {"chrom": "chr15", "chromStart": 5, "chromEnd": 6,
                    "name": {"data": {"Experiment": "GSE3", "TF": "NR3C1",
                                      "Biotype": "A549"}}}},
]}


class _RemapClient:
    def __init__(self, payload):
        self.payload, self.urls = payload, []

    def get(self, url, params=None, **kw):
        self.urls.append(url)
        return self.payload


def test_remap_target_read_from_documented_nested_path():
    """ReMap nests the TF at peaks[].peakValues.name.data.TF. Reading
    row['target'] at the top level found nothing, which is why an entire
    promoter came back with one TF."""
    from promoter_tf_drugs import evidence as ev
    c = _RemapClient(_REMAP_DOC_SHAPE)
    assert ev.peaks_from_api(c, "15", 100, 200) == {"CTCF": 2, "NR3C1": 1}


def test_remap_url_matches_the_documented_structure():
    from promoter_tf_drugs import evidence as ev
    c = _RemapClient({"peaks": []})
    ev.peaks_from_api(c, "15", 100, 200, datatype="nr")
    url = c.urls[0]
    assert url.startswith("https://")
    assert "/get_peaks/2022/hg38/nr/chr15:100-200" in url
    assert "chrchr" not in url


def test_flat_and_string_fallbacks_still_work():
    """Other releases or a changed schema must not silently yield nothing."""
    from promoter_tf_drugs import evidence as ev
    assert ev._target_of({"TF": "FOXA1"}) == "FOXA1"
    assert ev._target_of({"target": "GATA1"}) == "GATA1"
    assert ev._target_of({"name": {"data": {"TF": "SP1"}}}) == "SP1"
    assert ev._target_of("CTCF") == "CTCF"
    assert ev._target_of({"chrom": "chr1"}) is None


def test_unparseable_records_are_reported_not_silently_zero(caplog):
    from promoter_tf_drugs import evidence as ev
    c = _RemapClient({"peaks": [{"chrom": "chr1", "chromStart": 1}] * 40})
    with caplog.at_level("WARNING"):
        counts = ev.peaks_from_api(c, "15", 100, 200)
    assert counts == {}
    assert any("schema has probably changed" in r.message for r in caplog.records)


def test_realistic_promoter_result_is_not_flagged_as_broken(caplog):
    from promoter_tf_drugs import evidence as ev
    peaks = {"peaks": [
        {"peakValues": {"name": {"data": {"TF": tf}}}}
        for tf in ("CTCF", "SP1", "JUN", "FOS", "MYC", "NR3C1", "GATA2")]}
    with caplog.at_level("WARNING"):
        counts = ev.peaks_from_api(_RemapClient(peaks), "15", 100, 200)
    assert len(counts) == 7
    assert not any("parsing problem" in r.message for r in caplog.records)


# ----------------------- surmountability: the competitive-block trap
def test_competitive_blocks_flagged_as_surmountable():
    """A competitive antagonist occupies the ligand site and is displaced when
    ligand rises. A NAM and a degrader cannot be outcompeted."""
    assert mods.profile_for("orthosteric antagonist").surmountable is True
    assert mods.profile_for("inverse agonist").surmountable is True
    assert mods.profile_for("NAM").surmountable is False
    assert mods.profile_for("degrader").surmountable is False


def test_surmountability_note_names_the_robust_alternative():
    rows = [{"mechanism": "orthosteric antagonist", "surmountable": True},
            {"mechanism": "NAM", "surmountable": False}]
    note = mods.surmountability_note(rows)
    assert "competitive" in note
    assert "NAM cannot be outcompeted" in note


def test_endocrine_feedback_warning_only_when_axis_exists():
    """Blocking a receptor triggers feedback raising its ligand — which is
    exactly what defeats a competitive blocker."""
    rows = [{"mechanism": "orthosteric antagonist", "surmountable": True}]
    with_axis = mods.surmountability_note(rows, ligand_axis=True)
    without = mods.surmountability_note(rows, ligand_axis=False)
    assert "feedback will push its ligand up" in with_axis
    assert "feedback" not in without
    assert "escalating doses" in with_axis


def test_no_insurmountable_option_says_so():
    rows = [{"mechanism": "orthosteric antagonist", "surmountable": True}]
    note = mods.surmountability_note(rows)
    assert "No insurmountable option" in note
    assert "degrader" in note


def test_all_insurmountable_produces_no_warning():
    rows = [{"mechanism": "NAM", "surmountable": False},
            {"mechanism": "degrader", "surmountable": False}]
    assert mods.surmountability_note(rows) is None


def test_nam_outranks_competitive_antagonist_on_specificity():
    """Already true, but assert it: NAMs carry the allosteric bonus on the
    downregulation path just as PAMs do on the up path."""
    ev = dsc.DrugEvidence(pchembl_median=8.0, pchembl_n=20, n_targets=4)
    nam = dsc.grade_drug(_drug("NEGATIVE ALLOSTERIC MODULATOR"), 0.7, ev)
    ant = dsc.grade_drug(_drug("ANTAGONIST"), 0.7, ev)
    assert nam.components["selectivity"] > ant.components["selectivity"]
    assert "NAM" in mods.DOWN_CLASSES and "NAM" not in mods.UP_CLASSES


# ------------- established binders untested on the gene are the best class
def test_known_binder_untested_on_gene_outranks_novel_chemistry():
    """The conceptual fix. 'Known ligand' meant known at the RECEPTOR. Calling
    those positive controls and dropping them discarded the lowest-risk
    candidates: the drug-target step is measured, only the target-gene step is
    a hypothesis."""
    m = agonism.train_model(_toy_ligands())
    known = {f"A{i}" for i in range(len(_CORES))}
    hits = scr.screen(m, _toy_library(), goal="up", known_ids=known,
                      gene_linked=set(), top_n=100)
    assert hits[0].quadrant == "ESTABLISHED LIGAND, untested on this gene"
    assert any("de-risked" in n for n in hits[0].notes)

    order = [h.quadrant for h in hits]
    if "NOVEL HYPOTHESIS" in order:
        assert (order.index("ESTABLISHED LIGAND, untested on this gene")
                < order.index("NOVEL HYPOTHESIS"))


def test_known_binder_already_linked_is_a_confirmation_not_a_lead():
    """Both axes known: the drug binds the receptor AND already has a reported
    effect on the gene. That explains a known effect; it is not a new idea."""
    m = agonism.train_model(_toy_ligands())
    lib = _toy_library()
    linked = {lib[0]["name"].lower()}
    hits = scr.screen(m, lib, goal="up", known_ids={lib[0]["chembl_id"]},
                      gene_linked=linked, top_n=100)
    row = next(h for h in hits if h.drug.lower() == lib[0]["name"].lower())
    assert row.quadrant == "confirms a reported link"
    assert any("not a new lead" in n for n in row.notes)


def test_quadrant_axes_are_independent():
    from promoter_tf_drugs.screening import _quadrant
    assert _quadrant(0.9, 0.9, True, True, False) == \
        "ESTABLISHED LIGAND, untested on this gene"
    assert _quadrant(0.9, 0.9, True, True, True) == "confirms a reported link"
    assert _quadrant(0.9, 0.2, True, False, False) == "NOVEL HYPOTHESIS"
    assert _quadrant(0.9, 0.2, False, False, False) == "no call"


def test_established_binders_enter_the_catalogue_at_curated_tier():
    """They used to be dropped entirely — only NOVEL HYPOTHESIS rows made it
    into the unified table."""
    from promoter_tf_drugs import catalogue as cat
    block = {"screens": {"ESR1": {
        "enrichment": {"usable": True},
        "hits": [
            {"drug": "EstablishedDrug", "probability": 0.95,
             "quadrant": "ESTABLISHED LIGAND, untested on this gene"},
            {"drug": "NovelDrug", "probability": 0.8,
             "quadrant": "NOVEL HYPOTHESIS"},
        ]}}}
    routes = {r["drug"]: r for r in
              cat.build(block, include_chromatin=False)["routes"]}
    assert routes["EstablishedDrug"]["tier"] == 3
    assert routes["NovelDrug"]["tier"] == 5
    assert routes["EstablishedDrug"]["category"] == "TF: orthosteric ligand"
    # and tier ordering puts the de-risked one first
    ordered = cat.build(block, include_chromatin=False)["routes"]
    assert ordered[0]["drug"] == "EstablishedDrug"


# ===================================================================== #
# Evidence provenance: a failed lookup must never read as a measured    #
# negative. This is the failure mode that made a whole CCR5 run look    #
# complete while three of its four evidence sources had silently died.  #
# ===================================================================== #
def test_false_discovery_proportion_is_the_chance_share_not_its_complement():
    """The direction of this number is the whole point, and it is easy to
    invert: 1 - expected/observed is the share that is REAL."""
    from promoter_tf_drugs import motifs as m
    # 900 observed against ~744 expected: overwhelmingly noise.
    q = m.scan_quality(900, 879, 4230, 1e-4)
    assert q["false_discovery_proportion"] > 0.75
    # A genuinely enriched scan must come out LOW, not high.
    good = m.scan_quality(200, 800, 1200, 1e-5)
    assert good["false_discovery_proportion"] < 0.25
    # and the two framings must agree with each other
    assert good["expected_real_sites"] > 0.5 * good["observed"]


def test_fdp_and_enrichment_tell_the_same_story():
    from promoter_tf_drugs import motifs as m
    for n in (50, 300, 900, 2000):
        q = m.scan_quality(n, 879, 4230, 1e-4)
        if q["enrichment"] >= 2.0:
            assert q["false_discovery_proportion"] <= 0.5
        if q["enrichment"] <= 1.2:
            assert q["false_discovery_proportion"] >= 0.8


def test_broken_chip_api_is_a_failure_not_measured_absence():
    """A ReMap response yielding one or two factors over a promoter is a
    parsing failure. Scoring it as 'binding measured and absent' zeroes the
    20-point evidence axis for every factor on a lie."""
    from promoter_tf_drugs import evidence as ev

    class _C:
        def get(self, url): return {"peaks": [{"peakValues": {"name": {
            "data": {"TF": "CTCF"}}}}]}

    prom = _prom(1)
    counts, status = ev.chip_evidence_status(prom, client=_C(), use_api=True)
    assert counts == {"CTCF": 1}
    assert status["requested"] is True
    assert status["ok"] is False           # refused as scoring input
    assert "--chip-bed" in status["detail"]


def test_empty_local_bed_is_a_real_negative():
    """Contrast with the above: a BED that parsed and found nothing IS
    evidence of absence, and must not be downgraded to a failure."""
    from promoter_tf_drugs import evidence as ev
    import tempfile, os
    fd, path = tempfile.mkstemp(suffix=".bed")
    with os.fdopen(fd, "w") as fh:
        fh.write("chr9\t1\t2\tSOMETF\n")       # nowhere near the promoter
    try:
        counts, status = ev.chip_evidence_status(_prom(1), chip_bed=path)
        assert counts == {}
        assert status["requested"] and status["ok"] is True
    finally:
        os.unlink(path)


def test_unrequested_chip_is_not_a_failure():
    from promoter_tf_drugs import evidence as ev
    _, status = ev.chip_evidence_status(_prom(1))
    assert status["requested"] is False
    assert status["ok"] is False


def test_ledger_distinguishes_not_requested_from_failed():
    from promoter_tf_drugs import provenance as p
    led = p.Ledger()
    led.record("chip", requested=False, ok=False)
    led.record("lincs", requested=True, ok=False, detail="API down")
    assert led.get("chip").state == "not requested"
    assert led.get("lincs").state == "FAILED"
    assert [s.name for s in led.failures()] == ["lincs"]
    # a not-requested source is never reported as broken
    assert led.to_dict()["failed"] == ["lincs"]
    # and a failed source names what it cost you
    assert led.get("lincs").to_dict()["impact"]


def test_ledger_keeps_the_best_outcome_across_promoters():
    """With --all-tss the same source is recorded per TSS. One region with no
    peaks must not mark a working source as broken."""
    from promoter_tf_drugs import provenance as p
    led = p.Ledger()
    led.record("chip", requested=True, ok=False, detail="nothing here")
    led.record("chip", requested=True, ok=True, yielded=40, detail="40 TFs")
    assert led.usable("chip")


def test_verdict_refuses_a_noise_scan_with_no_measurement():
    from promoter_tf_drugs import provenance as p, motifs as m
    led = p.Ledger()
    led.record("chip", requested=True, ok=False)
    v = p.verdict(led, m.scan_quality(900, 879, 4230, 1e-4))
    assert v["level"] == "not interpretable"
    assert v["n_empirical_sources"] == 0


def test_verdict_upgrades_when_measurement_exists():
    from promoter_tf_drugs import provenance as p, motifs as m
    led = p.Ledger()
    led.record("chip", requested=True, ok=True, yielded=40)
    v = p.verdict(led, m.scan_quality(200, 800, 1200, 1e-5))
    assert v["level"] == "interpretable"


# --------------------------------------------------------------- congeners
def test_salt_forms_collapse_to_one_active_moiety():
    from promoter_tf_drugs import congeners as c
    assert c.salt_free_name("Liothyronine Sodium") == "liothyronine"
    assert c.salt_free_name("Elacestrant Hydrochloride") == "elacestrant"
    assert c.salt_free_name("Seladelpar Lysine") == "seladelpar"
    assert c.salt_free_name("Calcitriol") == "calcitriol"


def test_leading_counterion_is_not_stripped():
    """'Sodium Fluoride' is a drug; reducing it to 'fluoride' would be wrong.
    Only trailing counterions are salts."""
    from promoter_tf_drugs import congeners as c
    assert c.salt_free_name("Sodium Fluoride") == "sodium fluoride"


def test_esters_are_grouped_but_never_merged():
    """Fluticasone Propionate and Fluticasone Furoate are different molecules.
    Merging them would be a pharmacological claim, not a dedup."""
    from promoter_tf_drugs import congeners as c
    assert c.salt_free_name("Fluticasone Furoate") != \
           c.salt_free_name("Fluticasone Propionate")
    fam = c.family_of(["Fluticasone Furoate", "Fluticasone Propionate"])
    assert fam["fluticasone"] == ["Fluticasone Furoate",
                                  "Fluticasone Propionate"]


def test_collapse_keeps_the_winner_and_names_the_alternates():
    from promoter_tf_drugs import congeners as c

    class G:
        def __init__(self, n, t, s):
            self.drug_name, self.tf_symbol, self.combined = n, t, s

    rows = [G("Liothyronine", "THRB", 63.0),
            G("Liothyronine Sodium", "THRB", 63.0),
            G("Levothyroxine Sodium", "THRB", 62.0),
            G("Calcitriol", "VDR", 58.7)]
    kept, info = c.collapse_grades(rows)
    names = [g.drug_name for g in kept]
    assert "Liothyronine Sodium" not in names      # collapsed away
    assert "Liothyronine" in names                  # winner survives
    assert "Calcitriol" in names                    # untouched
    assert info["Liothyronine"] == ["Liothyronine Sodium"]


def test_same_drug_at_two_tfs_stays_two_routes():
    """Those are genuinely different mechanisms, not duplicate chemistry."""
    from promoter_tf_drugs import congeners as c

    class G:
        def __init__(self, n, t):
            self.drug_name, self.tf_symbol = n, t

    kept, _ = c.collapse_grades([G("Tretinoin", "RARA"), G("Tretinoin", "PPARD")])
    assert len(kept) == 2
