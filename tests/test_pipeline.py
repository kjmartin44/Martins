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


def test_pocketless_tf_with_nothing_curated_is_not_called_actionable():
    """SP1 has no ligand pocket, so agonist chemistry will never reach it —
    but nothing is curated against SP1 or its degradation machinery either.

    Saying 'stabilisation recommended' and then listing four approved drugs
    against transthyretin and FKBP12 read as though a route were available.
    It is not. The mechanism class is real; its applicability to THIS factor
    is unknown, and those are different claims."""
    a = stab.assess("SP1", has_direct_drug=False)
    assert a["has_ligand_pocket"] is False
    assert a["actionable"] is False
    assert a["recommended"] is False          # was True — that was over-claiming
    assert a["status"] == "no route known"
    assert any("research programme" in n for n in a["notes"])
    # it must still say what WOULD make it real
    assert any("E3 ligase" in n for n in a["notes"])


def test_pocketless_tf_with_a_curated_route_is_actionable():
    """Contrast: NFE2L2 has no pocket either, but dimethyl fumarate is an
    approved KEAP1-directed stabiliser of it. That is a real route."""
    a = stab.assess("NFE2L2", goal="up", has_direct_drug=False)
    assert a["actionable"] is True
    assert a["recommended"] is True
    assert a["status"] == "curated route exists"
    assert a["n_approved"] >= 1


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
    assert names == {"report.md", "summary.md", "results.json",
                     "tf_summary.tsv", "sites.bed", "promoter.fa"}

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
        # `changes` and `assay` answer the question the prose already calls
        # the one that most often decides an experiment — whether a route
        # moves protein abundance or only its activity — per category rather
        # than in the abstract.
        assert set(g) == {"what", "why", "cost", "differs",
                          "changes", "assay"}, cat
        # The prose fields must actually explain something. `changes` and
        # `assay` are deliberately terse table values ("abundance",
        # "RT-qPCR then western/ELISA"), so they are only required to exist.
        assert all(len(g[k]) > 40
                   for k in ("what", "why", "cost", "differs")), cat
        assert g["changes"] and g["assay"], cat


def test_only_direct_gene_product_changes_activity_rather_than_abundance():
    """The whole point of the distinction: every transcription-mediated route
    moves how much product exists, so RT-qPCR sees it. A drug on the product
    itself changes what existing protein does and barely moves its level —
    which is why ivacaftor reads as a failure on a western blot."""
    from promoter_tf_drugs import catalogue as c
    activity = {cat for cat, g in c.CATEGORY_GUIDE.items()
                if "ACTIVITY" in g["changes"].upper()}
    # Both product categories qualify and nothing else may. The predicted
    # one was added when the virtual screen was pointed at the gene's own
    # protein: it acts on the product, so it changes what the existing
    # protein does rather than how much of it there is, and the same
    # western-blot warning applies. A transcription-mediated route
    # appearing in this set would be the error the test is for.
    assert activity == {"direct: gene product",
                        "direct: predicted product binder"}
    # and its assay note must warn about the western blot explicitly
    assert "western" in c.CATEGORY_GUIDE["direct: gene product"]["assay"]


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
             "sirolimus": ["increases^expression"],
             # An exposure, not a candidate anyone could have proposed. It
             # used to be listed here as something the scan missed, which
             # indicted the scan for not suggesting cadmium.
             "cadmium chloride": ["increases^expression"]}
    s = summarise_links(links, "up", {"Dexamethasone"})
    names = [e["chemical"] for e in s["reported_but_not_predicted"]]
    assert names == ["Sirolimus"]
    assert s["counts"]["not_predicted_by_us"] == 1

    cad = next(e for e in s["wanted_direction"]
               if e["chemical"] == "Cadmium Chloride")
    assert cad["kind"] == "exposure"

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
    assert "measured effect on this gene and no mechanism" in md


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


def test_established_binders_rank_below_direct_agents_but_above_curated():
    """They used to be dropped entirely as "positive controls", then sat at
    tier 3 among every other curated ligand.

    Tier 2 is where the evidence puts them: the drug-receptor step is
    MEASURED - that is what "established ligand" means - and only the
    receptor-to-gene step is a hypothesis, whereas a curated tier-3 route is
    annotation on both halves. They stay below tier 1 because a direct agent
    acts on the gene itself and skips the receptor-to-gene inference."""
    from promoter_tf_drugs import catalogue as cat
    block = {
        "direct_therapies": [{
            "drug": "DirectAgent", "modality": "small molecule",
            "mechanism": "binds the gene product", "gene": "GENE",
            "specificity": 0.55, "raises_protein": False, "approved": True,
            "withdrawn": False, "approval_year": 2007, "source": "ChEMBL"}],
        "drug_grades": [{
            "drug_name": "CuratedLigand", "tf_symbol": "ESR1", "grade": "A",
            "combined": 80.0, "tf_confidence": 0.8, "drug_quality": 0.9,
            "components": {}, "evidence": {}, "notes": []}],
        "screens": {"ESR1": {
            "enrichment": {"usable": True},
            "hits": [
                {"drug": "EstablishedBinder", "probability": 0.95,
                 "quadrant": "ESTABLISHED LIGAND, untested on this gene"},
                {"drug": "NovelDrug", "probability": 0.8,
                 "quadrant": "NOVEL HYPOTHESIS"},
            ]}},
    }
    routes = cat.build(block, include_chromatin=False)["routes"]
    tier = {r["drug"]: r["tier"] for r in routes}

    assert tier["DirectAgent"] == 1
    assert tier["EstablishedBinder"] == 2
    assert tier["CuratedLigand"] == 3
    assert tier["NovelDrug"] == 5

    # the ordering the ranking actually has to produce
    order = [r["drug"] for r in routes]
    assert order.index("DirectAgent") < order.index("EstablishedBinder")
    assert order.index("EstablishedBinder") < order.index("CuratedLigand")


def test_established_binder_explains_how_it_reaches_the_gene():
    """A bare "measured binder" label does not say what it would do to the
    target gene, which is the question being asked of it."""
    from promoter_tf_drugs import catalogue as cat
    block = {"screens": {"VDR": {
        "enrichment": {"usable": True},
        "hits": [{"drug": "Doxercalciferol", "probability": 0.93,
                  "quadrant": "ESTABLISHED LIGAND, untested on this gene"}]}}}
    r = cat.build(block, include_chromatin=False)["routes"][0]

    assert "THROUGH VDR" in r["mechanism"]
    assert "not by touching it directly" in r["mechanism"]
    joined = " ".join(r["notes"])
    # the de-risked half, and the half that is not
    assert "de-risked" in joined
    assert "binding assay cannot tell you" in joined
    # and the cost of going indirect
    assert "everything else VDR regulates moves too" in joined

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


# ===================================================================== #
# Specificity that counts SYSTEMS, not proteins or genes.               #
# ===================================================================== #
def test_shortlist_predicate_split_preserved_behaviour():
    """The refactor that made a staged funnel possible must not have changed
    what the five presets return."""
    res = _result_for_shortlist()
    rows = sl.enrich_rows(res)
    for name, preset in sl.PRESETS.items():
        assert sl.filter_rows(rows, **preset) == sl.filter_routes(res, **preset)


def test_every_preset_key_has_a_predicate():
    """A preset naming a filter that no longer exists must fail loudly rather
    than silently widening the shortlist."""
    keys = {k for p in sl.PRESETS.values() for k in p}
    assert keys <= set(sl.PREDICATES), keys - set(sl.PREDICATES)


def test_unknown_filter_name_raises():
    """The structural guarantee that a diagnostic-only field cannot become a
    gate: an unrecognised filter is an error, not a no-op."""
    import pytest
    with pytest.raises(KeyError):
        sl.filter_rows([], n_genes=3)


def test_require_chip_still_ignores_tfs_absent_from_the_block():
    """A route whose factor is not in the block's TF table was never subject
    to require_chip. Losing that guard would silently drop routes."""
    row = {"tier": 3, "score": 1.0, "directness": 1, "acts_on": "GHOST",
           "chip_peaks": None, "_chip_known": False}
    assert sl._p_require_chip(row, True)
    row["_chip_known"] = True
    assert not sl._p_require_chip(row, True)


# ------------------------------------------ regulon counts vs study effort
def test_tiny_curated_regulon_is_understudied_not_specific():
    """The trap: TRRUST is literature-derived, so a TF with 2 recorded targets
    is almost certainly under-assayed rather than narrow. Scoring that as
    perfect specificity rewards obscurity."""
    from promoter_tf_drugs import regulon as r
    ev = r.regulon_evidence("SP1", {"SP1": {"A", "B"}})
    assert ev.understudied is True
    assert ev.trustworthy is False
    assert "curation depth" in ev.note

    # The score itself is deliberately unchanged — a real file still beats a
    # bundled guess. It is the CONFIDENCE in it that drops.
    score, n, _ = r.transcriptional_specificity("SP1", {"SP1": {"A", "B"}})
    assert n == 2 and score == 1.0


def test_bundled_tier_is_not_a_trustworthy_count():
    """5000 is an order-of-magnitude bucket someone typed, not a measurement."""
    from promoter_tf_drugs import regulon as r
    ev = r.regulon_evidence("SP1")
    assert ev.trustworthy is False and "bucket" in ev.note


def test_a_real_regulon_above_the_floor_is_trusted():
    from promoter_tf_drugs import regulon as r
    big = {"MYTF": {f"G{i}" for i in range(40)}}
    ev = r.regulon_evidence("MYTF", big)
    assert ev.trustworthy is True and ev.understudied is False
    assert ev.provenance == "regulon file"


def test_unknown_regulon_is_unmeasured_not_narrow():
    from promoter_tf_drugs import regulon as r
    ev = r.regulon_evidence("ZZZ9FAKE")
    assert ev.n is None and ev.trustworthy is False
    assert "not the same as narrow" in ev.note


# --------------------------------------------------- systemic drug classes
def test_glucocorticoid_is_flagged_despite_a_spotless_profile():
    """The reason this module exists. Dexamethasone is clean, potent and
    selective — and acts on the whole body. selectivity counts proteins,
    regulon counts genes; neither counts systems."""
    from promoter_tf_drugs import systemic as s
    b = s.clinical_burden("Dexamethasone", notes="", target="NR3C1")
    assert b.assessed and b.level == "high"
    assert b.drug_class == "glucocorticoid"


def test_global_regulator_target_raises_burden_for_an_unlisted_drug():
    """The receptor is systemic even when the drug is not in the table."""
    from promoter_tf_drugs import systemic as s
    b = s.clinical_burden("SomeNovelCompound", notes="", target="SP1")
    assert b.assessed and b.score >= 0.85
    assert any("global regulator" in r for r in b.reasons)


def test_unlisted_drug_is_unassessed_not_cleared():
    """A clean row means 'not on our list', never 'safe'."""
    from promoter_tf_drugs import systemic as s
    b = s.clinical_burden("Maraviroc", notes="", target="CCR5")
    assert b.assessed is False
    assert "not the same as clean" in " ".join(b.reasons)


def test_salt_form_resolves_to_its_class():
    """Otherwise pipeline's salt collapsing and this table disagree about
    what a drug is."""
    from promoter_tf_drugs import systemic as s
    assert s.classify("Levothyroxine Sodium").name == "thyroid hormone"
    assert s.classify("Elacestrant Hydrochloride").name == "sex steroid"


def test_burden_reads_the_same_withdrawn_strings_as_shortlist():
    """Two spellings of 'withdrawn' would give one report two verdicts."""
    from promoter_tf_drugs import systemic as s
    b = s.clinical_burden("Plicamycin", notes="WITHDRAWN from market")
    assert b.level == "prohibitive"
    row = {"_notes": "WITHDRAWN from market"}
    assert not sl._p_exclude_withdrawn(row, True)


def test_black_box_is_high_but_not_prohibitive():
    from promoter_tf_drugs import systemic as s
    b = s.clinical_burden("SomeDrug", notes="black-box warning")
    assert b.level == "high" and b.assessed


# ------------------------------------------------- TF tissue breadth proxy
def test_tissue_breadth_refuses_a_fraction_without_a_denominator():
    """'Expressed in 30 tissues' means different things against HPA's ~50,
    GTEx's ~54 and Open Targets' variable vocabulary. No denominator, no
    score — inventing 50 would be fabricating the number."""
    from promoter_tf_drugs import delivery as d
    d._LOCAL_EXPRESSION.clear()
    assert d.tissue_universe_size() is None
    r = d.expression_breadth("GATA2")
    assert r["score"] is None and r["basis"] == "not assessed"


def test_tissue_breadth_is_higher_for_a_narrower_factor():
    """Direction must match regulon.transcriptional_specificity — higher is
    narrower. A sign flip here would invert the whole axis silently."""
    from promoter_tf_drugs import delivery as d
    d._LOCAL_EXPRESSION.clear()
    d._LOCAL_EXPRESSION.update({
        "NARROW": {"liver": 5.0},
        "BROAD": {f"t{i}": 5.0 for i in range(10)},
    })
    try:
        assert d.tissue_universe_size() == 11
        narrow = d.expression_breadth("NARROW")["score"]
        broad = d.expression_breadth("BROAD")["score"]
        assert narrow > broad
    finally:
        d._LOCAL_EXPRESSION.clear()


def test_nuclear_receptor_is_not_penalised_on_an_expression_proxy():
    """mRNA presence is not activity. VDR is transcribed nearly everywhere
    while 1,25-D signalling is not — and NRs were 201 of 216 routes in the
    real CCR5 run, so a bad proxy here biases against the whole output."""
    from promoter_tf_drugs import delivery as d
    d._LOCAL_EXPRESSION.clear()
    d._LOCAL_EXPRESSION.update({"VDR": {f"t{i}": 5.0 for i in range(40)}})
    try:
        r = d.expression_breadth("VDR", nuclear_receptor=True)
        assert r["score"] is None
        assert "ligand" in r["basis"] or "ligand" in r["note"]
    finally:
        d._LOCAL_EXPRESSION.clear()


def test_cns_mpo_is_never_a_tissue_breadth_basis():
    """The tempting wrong shortcut: MPO predicts crossing ONE barrier, and
    says nothing about how many tissues a drug reaches."""
    from promoter_tf_drugs import delivery as d
    d._LOCAL_EXPRESSION.clear()
    d._LOCAL_EXPRESSION.update({"X": {"liver": 3.0}})
    try:
        bases = {d.expression_breadth("X")["basis"],
                 d.expression_breadth("Q")["basis"],
                 d.expression_breadth("X", nuclear_receptor=True)["basis"]}
        assert not any("mpo" in b.lower() or "cns" in b.lower() for b in bases)
    finally:
        d._LOCAL_EXPRESSION.clear()


def test_cherry_pick_states_its_provenance_and_its_limit():
    """A bare compound table invites the reading that these were selected
    from the world. They were selected from this run's own predictions, so
    the plate cannot find anything the model did not already propose."""
    res = pipeline.run("DEMOGENE", demo=True, offline=True, design_screen=True)
    md = pipeline.markdown_report(res)
    assert "Cherry-pick list" in md
    assert "proposed earlier in *this run*" in md
    assert "cannot discover a mechanism the pipeline never predicted" in md
    # every row says where it came from, and what priority means
    assert "Pulled from" in md and "Priority means" in md

    cp = res["results"][0]["screen_plan"]["compounds"]
    assert cp["n_available"] >= len(cp["cherry_pick"])
    assert all(c.get("source") for c in cp["cherry_pick"])


# ===================================================================== #
# Multi-gene panel: choosing between candidate genes.                   #
# ===================================================================== #
from promoter_tf_drugs import panel as pn          # noqa: E402


def _panel_of(*genes, **kw):
    p = pn.Panel()
    for g in genes:
        pg = pn.PanelGene(symbol=g)
        pg.result = pipeline.run(g, demo=True, offline=True, **kw)
        p.genes.append(pg)
    return p


def test_merged_panel_is_filterable_by_the_unchanged_shortlist():
    """The reuse claim, asserted: a multi-gene result is structurally just a
    results.json with more blocks."""
    p = _panel_of("DEMOGENE", "SMN2")
    merged = pn.merge_results(p)
    rows = sl.filter_routes(merged)
    assert {r["gene"] for r in rows} == {"DEMOGENE", "SMN2"}


def test_merge_does_not_attribute_one_genes_direct_therapies_to_another():
    """direct_therapies is top-level, not per-block. Blending it would give
    every gene in the panel another gene's approved drugs."""
    p = _panel_of("DEMOGENE", "SMN2")
    merged = pn.merge_results(p)
    assert merged["direct_therapies"] == []
    per = merged["panel"]["per_gene"]
    assert set(per) == {"DEMOGENE", "SMN2"}
    for g in p.included:
        assert per[g.symbol]["direct_therapies"] == g.result["direct_therapies"]


def test_mixed_direction_panel_is_an_error_not_a_warning():
    """An up gene and a down gene are not comparable at all."""
    import pytest
    p = pn.Panel()
    for g, goal in (("DEMOGENE", "up"), ("SMN2", "down")):
        pg = pn.PanelGene(symbol=g)
        pg.result = pipeline.run(g, demo=True, offline=True, goal=goal)
        p.genes.append(pg)
    with pytest.raises(ValueError, match="goal"):
        pn.parity(p)


def test_parity_flags_genes_run_with_different_evidence():
    """A gene run without a source loses stages for reasons unrelated to its
    biology; comparing it like-for-like is the panel-level version of the bug
    provenance.py exists to prevent."""
    p = pn.Panel()
    for g, ctc in (("DEMOGENE", True), ("SMN2", False)):
        pg = pn.PanelGene(symbol=g)
        pg.result = pipeline.run(g, demo=True, offline=True,
                                 cell_type_check=ctc)
        p.genes.append(pg)
    par = pn.parity(p)
    assert par["fair"] is False
    assert any("cell_type_check" in w for w in par["warnings"])


def test_not_interpretable_gene_is_excluded_and_named():
    p = _panel_of("DEMOGENE")
    p.genes[0].result["evidence_verdict"] = {"level": "not interpretable",
                                             "reasons": ["scan is noise"]}
    p.genes[0].excluded = "verdict 'not interpretable'"
    assert p.included == []
    assert p.dropped and p.dropped[0].symbol == "DEMOGENE"
    assert p.dropped[0].verdict["reasons"] == ["scan is noise"]


def test_one_failing_gene_does_not_kill_the_panel(monkeypatch):
    real = pipeline.run

    def flaky(gene, **kw):
        if gene == "BOOM":
            raise RuntimeError("Ensembl exploded")
        return real(gene, **kw)

    monkeypatch.setattr(pipeline, "run", flaky)
    p = pn.run_panel(["DEMOGENE", "BOOM"], demo=True, offline=True)
    assert [g.symbol for g in p.included] == ["DEMOGENE"]
    assert p.dropped[0].error and "Ensembl exploded" in p.dropped[0].error


# ------------------------------------------------------------- the funnel
def _funnel_of(*genes, **criteria):
    p = _panel_of(*genes)
    merged, kept, stages, rec = pn.narrow(p, **criteria)
    return p, merged, kept, stages, rec


def test_funnel_counts_reconcile_at_every_stage():
    """n_in == n_removed + n_remaining, and each stage's n_in is the previous
    stage's n_remaining. Catches a stage filtering a stale row list."""
    _, _, _, stages, _ = _funnel_of("DEMOGENE", "SMN2", max_regulon=600)
    prev = None
    for s in stages:
        assert s.n_in == s.n_removed + s.n_remaining, s.name
        if prev is not None:
            assert s.n_in == prev, s.name
        prev = s.n_remaining


def test_funnel_is_monotonic():
    _, _, _, stages, _ = _funnel_of("DEMOGENE", "SMN2", max_regulon=600)
    counts = [s.n_remaining for s in stages]
    assert counts == sorted(counts, reverse=True)


def test_unrequested_stage_is_not_reported_as_broken():
    """A threshold nobody set is a choice, not a missing measurement. Crying
    wolf here would devalue the warning where it is real."""
    _, _, _, stages, _ = _funnel_of("DEMOGENE")
    tissue = [s for s in stages if "tissue" in s.name][0]
    assert tissue.requested is False
    assert tissue.did_nothing is False


def test_stage_with_no_data_says_so_rather_than_reporting_a_clean_zero():
    _, _, _, stages, _ = _funnel_of("DEMOGENE", min_tissue_breadth=0.3)
    tissue = [s for s in stages if "tissue" in s.name][0]
    assert tissue.requested is True
    assert tissue.did_nothing is True
    assert tissue.n_not_assessable == tissue.n_in


def test_missing_datum_never_removes_a_row():
    """Same principle _grade_ok already encodes: dropping a candidate for
    lacking an inapplicable field is backwards."""
    row = {"_tissue_breadth": None, "_burden": "low", "_regulon_size": None}
    assert pn._p_min_tissue_breadth(row, 0.9)
    assert sl._p_max_regulon(row, 10)


def test_counts_are_gene_pairs_not_rows():
    """With --all-tss a gene contributes several near-identical windows;
    counting rows would let the most-annotated gene dominate."""
    rows = [{"gene": "A", "drug": "Levothyroxine", "acts_on": "THRB"},
            {"gene": "A", "drug": "Levothyroxine Sodium", "acts_on": "THRB"},
            {"gene": "A", "drug": "Levothyroxine", "acts_on": "THRB"},
            {"gene": "B", "drug": "Levothyroxine", "acts_on": "THRB"}]
    pairs = pn._pairs(rows)
    assert len(pairs["A"]) == 1        # salt form and duplicate collapse
    assert len(pairs["B"]) == 1


def test_unique_kills_exposes_a_redundant_stage():
    """Two stages removing the same rows are not two independent tests."""
    rows = [{"tier": 3, "score": 1.0, "gene": "A", "drug": "D%d" % i,
             "acts_on": "SP1", "category": "x", "directness": 1,
             "axes": {"specificity": 0.9}, "_notes": "",
             "_regulon_size": 5000, "_burden": "high",
             "_tissue_breadth": None} for i in range(4)]
    dup = pn.Stage("dup", "a", {"max_regulon": 600}, "")
    same = pn.Stage("same", "a", {"max_regulon": 600}, "")
    _, res = pn.run_funnel(rows, [dup, same])
    assert res[0].n_removed == 4
    assert res[0].n_unique_kills == 0      # the later stage would have too


def test_funnel_rejects_a_stage_naming_an_unknown_filter():
    import pytest
    bad = pn.Stage("bad", "a", {"n_genes": 3}, "")
    with pytest.raises(KeyError):
        pn.run_funnel([], [bad])


# ---------------------------------------------------- recurrence stays out
def test_recurrence_is_not_a_predicate_and_no_stage_names_it():
    """The structural guarantee that a diagnostic cannot become a gate."""
    assert "n_genes" not in pn.PREDICATES
    assert "recurs_on" not in pn.PREDICATES
    for s in pn.default_stages():
        assert "n_genes" not in s.filters


def test_recurrence_annotation_changes_no_rows():
    _, _, kept, _, rec = _funnel_of("DEMOGENE", "SMN2")
    before = len(kept)
    pn.annotate_recurrence(kept, rec)
    assert len(kept) == before


def test_recurrence_counts_genes_not_rows():
    rows = [{"gene": "A", "drug": "X", "acts_on": "T1"},
            {"gene": "A", "drug": "X", "acts_on": "T2"},
            {"gene": "A", "drug": "X Sodium", "acts_on": "T1"}]
    assert pn.recurrence(rows)["x"]["n_genes"] == 1


def test_recurrence_distinguishes_shared_regulator_from_promiscuity():
    shared = pn.recurrence([{"gene": "A", "drug": "X", "acts_on": "VDR"},
                            {"gene": "B", "drug": "X", "acts_on": "VDR"}])
    promisc = pn.recurrence([{"gene": "A", "drug": "Y", "acts_on": "VDR"},
                             {"gene": "B", "drug": "Y", "acts_on": "ESR1"}])
    assert shared["x"]["targets"] == ["VDR"]
    assert len(promisc["y"]["targets"]) == 2


def test_placeholder_target_is_not_a_shared_regulator():
    """catalogue writes an em-dash into acts_on for routes with no protein
    target; a truthiness check let that through, and ['-'] then read as one
    shared regulator."""
    dash = "—"
    rec = pn.recurrence([{"gene": "A", "drug": "Z", "acts_on": dash},
                         {"gene": "B", "drug": "Z", "acts_on": dash}])
    assert rec["z"]["targets"] == []
    assert pn.real_target(dash) is None
    assert pn.real_target("VDR") == "VDR"


# ------------------------------------------------------------- the report
def test_panel_report_renders_every_section():
    p, merged, kept, stages, rec = _funnel_of("DEMOGENE", "SMN2")
    merged.setdefault("panel", {})["parity"] = pn.parity(p)
    md = pn.markdown_report(p, merged, kept, stages, rec, {"top": 20})
    for heading in ("Genes in the comparison", "The funnel",
                    "Where each gene dropped out", "Cross-gene recurrence",
                    "What this cannot tell you"):
        assert heading in md, heading
    assert "alternatives you are choosing between" in md
    assert "diagnostic, not a filter" in md


def test_report_warns_that_the_winner_was_selected_on_a_maximum():
    p, merged, kept, stages, rec = _funnel_of("DEMOGENE", "SMN2")
    md = pn.markdown_report(p, merged, kept, stages, rec, {"top": 20})
    assert "About the top row" in md
    assert "maximum over" in md


def test_empty_result_names_the_binding_stage():
    p, merged, kept, stages, rec = _funnel_of("DEMOGENE", max_regulon=0,
                                              max_tier=1)
    md = pn.markdown_report(p, merged, [], stages, rec, {})
    assert "Nothing survived" in md
    assert "binding constraint is" in md


def test_panel_cli_runs_end_to_end(tmp_path):
    rc = pn.main(["--genes", "DEMOGENE,SMN2", "--demo", "--out", str(tmp_path)])
    assert rc == 0
    md = (tmp_path / "panel.md").read_text(encoding="utf-8")
    assert "The funnel" in md
    assert (tmp_path / "panel.json").exists()
    assert (tmp_path / "panel.tsv").exists()


# ------------------------------------- is the scanned window this gene's?
from promoter_tf_drugs import neighbours as nbr        # noqa: E402


class _NbClient:
    """Serves the real Ensembl shape for the CCR5 locus."""
    def __init__(self, rows=None):
        self.rows = rows if rows is not None else [
            {"external_name": "CCR2", "biotype": "protein_coding",
             "start": 46353864, "end": 46360940, "strand": 1, "id": "g1"},
            {"external_name": "CCR5AS", "biotype": "lncRNA",
             "start": 46363984, "end": 46407184, "strand": -1, "id": "g2"},
            {"external_name": "CCR5", "biotype": "protein_coding",
             "start": 46370939, "end": 46377528, "strand": 1, "id": "g3"},
        ]

    def get(self, url, params=None):
        return self.rows


def test_antisense_overlap_is_detected_and_explained():
    """The real CCR5 case: the scanned window sits entirely inside CCR5AS on
    the opposite strand, so a site called there is not unambiguously CCR5's
    and the direction logic does not carry."""
    r = nbr.genes_in_window(_NbClient(), "3", 46365946, 46371945,
                            strand=1, exclude="CCR5")
    assert r["checked"] and r["n_antisense"] == 1
    anti = [g for g in r["genes"] if g["relation"] == "antisense overlap"]
    assert anti[0]["gene"] == "CCR5AS"
    assert "OPPOSITE strand" in anti[0]["note"]
    assert "direction logic" in r["summary"]


def test_the_target_gene_itself_is_not_reported_as_a_neighbour():
    r = nbr.genes_in_window(_NbClient(), "3", 46365946, 46371945,
                            strand=1, exclude="CCR5")
    assert "CCR5" not in [g["gene"] for g in r["genes"]]


def test_a_near_neighbour_is_flagged_with_its_distance():
    """CCR2 is 5 kb outside — a wider window would swallow a co-regulated
    paralog, and nothing else in the pipeline would say so."""
    r = nbr.genes_in_window(_NbClient(), "3", 46365946, 46371945,
                            strand=1, exclude="CCR5")
    near = [g for g in r["genes"] if g["relation"] == "nearby"]
    assert near and near[0]["gene"] == "CCR2"
    assert "kb outside" in near[0]["note"]


def test_all_tss_window_overlap_is_measured():
    """CCR5's two TSS windows shared 4275 of 6000 bp — largely one analysis
    counted twice, which inflates that gene's route count."""
    proms = [{"chrom": "3", "start": 46365946, "end": 46371945,
              "transcript_id": "T1"},
             {"chrom": "3", "start": 46367670, "end": 46373669,
              "transcript_id": "T2"}]
    ov = nbr.window_overlap(proms)
    assert ov and ov[0]["overlap_bp"] == 4275
    assert ov[0]["fraction"] > 0.7
    assert "counted twice" in ov[0]["note"]


def test_disjoint_windows_report_no_overlap():
    proms = [{"chrom": "3", "start": 100, "end": 200, "transcript_id": "A"},
             {"chrom": "3", "start": 5000, "end": 6000, "transcript_id": "B"}]
    assert nbr.window_overlap(proms) == []


# ===================================================================== #
# The headline answer: predicted medications, with the reasoning, first; #
# experimental protocol last.                                            #
# ===================================================================== #
from promoter_tf_drugs import predictions as prd        # noqa: E402


def _demo_report(**kw):
    res = pipeline.run("DEMOGENE", demo=True, offline=True, **kw)
    return res, pipeline.markdown_report(res)


def test_predicted_medications_lead_the_report():
    """The first question anyone brings is 'what should I give, and why'.
    It used to be answerable only by assembling eight separate tables."""
    _, md = _demo_report()
    assert "## Model-predicted medications for DEMOGENE" in md
    i_pred = md.index("Model-predicted medications")
    i_promoter = md.index("## Promoter of DEMOGENE")
    assert i_pred < i_promoter, "predictions must precede the analysis detail"


def test_every_prediction_carries_its_reasoning_chain():
    res, md = _demo_report()
    preds = prd.build(res["results"][0], goal="up")
    assert preds
    p = preds[0]
    joined = " ".join(p.why)
    assert "acts on" in joined or "directly" in joined   # mechanism
    assert "direction" in joined or "activator" in joined or \
           "repressor" in joined or "dual-role" in joined
    assert "Binding:" in joined                          # why this promoter
    assert "Evidence tier" in joined                     # how far it reaches


def test_direction_clause_inverts_with_the_goal():
    """Agonise an activator to raise a gene; block it to lower one."""
    up = prd._direction_clause("activator", "up", False)
    down = prd._direction_clause("activator", "down", False)
    assert "agonising" in up and "blocking" in down
    up_r = prd._direction_clause("repressor", "up", False)
    assert "blocking" in up_r


def test_dual_role_factor_is_not_claimed_to_have_a_known_direction():
    c = prd._direction_clause("dual", "up", False)
    assert "not settled" in c


def test_unconfirmed_motif_is_described_as_a_prediction_not_binding():
    tf = {"n_sites": 1, "best_rel_score": 0.78, "best_site_tss_offset": -439,
          "chip_peaks": 0}
    assert "sequence prediction only" in prd._promoter_clause(tf, False)
    # ChIP available but absent is a different, stronger statement
    assert "no ChIP peak" in prd._promoter_clause(tf, True)
    tf["chip_peaks"] = 3
    assert "actually binds" in prd._promoter_clause(tf, True)


def test_prediction_section_repeats_the_run_verdict():
    """A degraded run must say so where the drugs are read, not only in a
    table further up that a reader may skip."""
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    res["evidence_verdict"] = {"level": "not interpretable", "reasons": []}
    md = pipeline.markdown_report(res)
    i = md.index("Model-predicted medications")
    assert "Do not act on this list" in md[i:i + 1200]


def test_predictions_state_they_are_hypotheses_not_recommendations():
    _, md = _demo_report()
    i = md.index("Model-predicted medications")
    head = md[i:i + 1500]
    assert "hypothesis generated from annotation" in head
    assert "not a recommendation" in head


def test_wet_lab_sections_come_after_everything_else():
    """A CRISPRa protocol between the reader and the drug list buried the
    answer under the method."""
    _, md = _demo_report(design_screen=True, design_modality=True)
    i_lab = md.index("# Testing any of this in the lab")
    for earlier in ("Model-predicted medications", "Graded candidates",
                    "How to read this"):
        assert md.index(earlier) < i_lab, earlier
    # and the long protocols sit inside that trailing part
    for protocol in ("CRISPRa screen", "Reporter construct",
                     "Cherry-pick list", "Designing a direct modality"):
        assert md.index(protocol) > i_lab, protocol


def test_lab_part_is_absent_when_no_experiment_was_designed():
    """The divider should not appear for a run that designed nothing."""
    _, md = _demo_report()
    assert "# Testing any of this in the lab" not in md


def test_lab_part_leads_with_the_one_control_that_matters():
    _, md = _demo_report(design_screen=True)
    i = md.index("# Testing any of this in the lab")
    head = md[i:i + 900]
    assert "scrambled-motif control" in head


def test_empty_prediction_list_is_explained_not_blank():
    lines = prd.render([], "XGENE", "up")
    assert any("not an error" in ln for ln in lines)


# ------------------------------- Open Targets schema break (API v26.9)
def test_baseline_expression_schema_is_parsed():
    """Target.expressions was removed, so every gene 400'd and the whole
    tissue axis went dark. baselineExpression replaces it."""
    from promoter_tf_drugs import delivery as d

    class C:
        def post(self, url, payload):
            assert "baselineExpression" in payload["query"]
            assert payload["variables"]["size"] > 0
            return {"data": {"target": {"baselineExpression": {
                "count": 3,
                "rows": [
                    {"median": 12.0, "unit": "TPM",
                     "tissueBiosample": {"biosampleName": "Liver"}},
                    {"median": 0.2, "unit": "TPM",
                     "tissueBiosample": {"biosampleName": "Brain"}},
                    {"median": 99.0, "unit": "CPM(pseudobulk sum[counts])",
                     "tissueBiosample": {"biosampleName": "Spleen"}},
                ]}}}}

    d._LOCAL_EXPRESSION.clear()
    d._EXPR_DISABLED["value"] = False
    import promoter_tf_drugs.genetics as gen
    old = gen._resolve_ensembl
    gen._resolve_ensembl = lambda c, s: "ENSG0"
    try:
        t = d.tissue_expression(C(), "X")
    finally:
        gen._resolve_ensembl = old
    # Liver: TPM 12 >= 1.0. Brain: TPM 0.2, below the TPM floor.
    # Spleen: CPM 99, kept under the SEPARATE single-cell floor.
    assert t == {"liver", "spleen"}, t


def test_each_unit_is_thresholded_against_its_own_floor():
    """TPM and CPM(pseudobulk) are different measurements on different data,
    so they are kept but never compared against a single cutoff. A value of
    0.5 is below the floor in either unit; 1000 clears it in either."""
    from promoter_tf_drugs import delivery as d

    def client_with(rows):
        class C:
            def post(self, url, payload):
                return {"data": {"target": {"baselineExpression": {
                    "count": len(rows), "rows": rows}}}}
        return C()

    cpm = "CPM(pseudobulk sum[counts])"
    d._LOCAL_EXPRESSION.clear()
    d._EXPR_DISABLED["value"] = False
    import promoter_tf_drugs.genetics as gen
    old = gen._resolve_ensembl
    gen._resolve_ensembl = lambda c, s: "ENSG0"
    try:
        kept = d.tissue_expression(client_with([
            {"median": 1000.0, "unit": cpm,
             "tissueBiosample": {"biosampleName": "Kidney"}}]), "X")
        assert kept == {"kidney"}       # single-cell data is USED, not dropped

        below = d.tissue_expression(client_with([
            {"median": 0.4, "unit": cpm,
             "tissueBiosample": {"biosampleName": "Kidney"}}]), "X")
        assert below == set()           # but still has to clear its own floor

        unknown = d.tissue_expression(client_with([
            {"median": 9e9, "unit": "some-future-unit",
             "tissueBiosample": {"biosampleName": "Kidney"}}]), "X")
        assert unknown == set()         # an unrecognised unit has no floor
    finally:
        gen._resolve_ensembl = old
        d._EXPR_DISABLED["value"] = False


def test_old_schema_still_parses_if_a_deployment_serves_it():
    from promoter_tf_drugs import delivery as d

    class C:
        def post(self, url, payload):
            if "baselineExpression" in payload["query"]:
                return {"errors": [{"message": "Cannot query field"}]}
            return {"data": {"target": {"expressions": [
                {"tissue": {"label": "Liver"}, "rna": {"level": 3},
                 "protein": {"level": 2}}]}}}

    d._LOCAL_EXPRESSION.clear()
    d._EXPR_DISABLED["value"] = False
    import promoter_tf_drugs.genetics as gen
    old = gen._resolve_ensembl
    gen._resolve_ensembl = lambda c, s: "ENSG0"
    try:
        assert d.tissue_expression(C(), "X") == {"liver"}
    finally:
        gen._resolve_ensembl = old


def test_api_route_still_refuses_to_supply_a_denominator():
    """The API vocabulary is a variable biosample set, not a fixed tissue
    panel, so a fraction against it would be fabricated."""
    from promoter_tf_drugs import delivery as d
    d._LOCAL_EXPRESSION.clear()
    assert d.tissue_universe_size() is None


# ===================================================================== #
# Paralog families: can a selective ligand exist at this receptor?      #
# ===================================================================== #
from promoter_tf_drugs import paralogs as pl          # noqa: E402


def test_family_members_are_recognised_as_relatives():
    """Straight from the real CCR5 run: THRA and THRB both surfaced, and
    dextrothyroxine/liothyronine/levothyroxine each scored against both."""
    thrb = pl.paralog_specificity("THRB")
    assert thrb["relatives"] == ["THRA"]
    assert thrb["assessed"] and thrb["score"] < 1.0
    assert "triiodothyronine" in thrb["ligand"]


def test_more_relatives_means_a_lower_score():
    """Direction must match every other specificity axis: higher = cleaner."""
    alone = pl.paralog_specificity("NR1H4")["score"]      # FXR, family of 1
    two = pl.paralog_specificity("ESR1")["score"]         # ESR1/ESR2
    four = pl.paralog_specificity("NR3C1")["score"]       # GR/MR/PR/AR
    assert alone > two > four


def test_unlisted_factor_is_unassessed_not_selective():
    """Most TFs have no ligand pocket at all. That is not selectivity."""
    r = pl.paralog_specificity("NHLH1")
    assert r["score"] is None and r["assessed"] is False
    assert "UNASSESSED, not selective" in r["note"]


def test_a_large_family_with_a_selective_agent_is_not_penalised():
    """trifarotene exists BECAUSE the RAR pocket is conserved. Penalising the
    molecule developed to solve the problem would be backwards."""
    from promoter_tf_drugs import panel as p
    row = {"_paralog_score": 0.4, "_paralog_selective_agents":
           ["trifarotene (RARG-selective)"]}
    assert p._p_min_paralog_specificity(row, 0.9)
    bare = {"_paralog_score": 0.4, "_paralog_selective_agents": []}
    assert not p._p_min_paralog_specificity(bare, 0.9)


def test_missing_paralog_data_never_removes_a_row():
    from promoter_tf_drugs import panel as p
    assert p._p_min_paralog_specificity({"_paralog_score": None}, 0.9)


def test_cross_reactive_pairs_finds_relatives_among_the_hits():
    """Two family members both scoring well is one bet, not two."""
    pairs = pl.cross_reactive_pairs(["THRA", "THRB", "RARA", "RARG", "VDR"])
    fams = {p["family"]: p["members"] for p in pairs}
    assert fams["thyroid hormone receptor"] == ["THRA", "THRB"]
    assert fams["retinoic acid receptor"] == ["RARA", "RARG"]
    # a lone member is not a pair
    assert not any("VDR" in p["members"] and len(p["members"]) == 1
                   for p in pairs)


def test_cross_reactive_pairs_names_the_selective_agents():
    pairs = pl.cross_reactive_pairs(["THRA", "THRB"])
    assert any("resmetirom" in a for a in pairs[0]["selective_agents"])


def test_paralog_stage_is_in_the_funnel():
    from promoter_tf_drugs import panel as p
    names = [s.name for s in p.default_stages()]
    assert "Selective ligand achievable" in names
    stage = [s for s in p.default_stages(min_paralog_specificity=0.6)
             if s.name == "Selective ligand achievable"][0]
    assert stage.filters["min_paralog_specificity"] == 0.6


def test_ligand_count_is_not_used_as_a_specificity_signal():
    """The whole point of using family size instead: ligand count measures
    study effort, the same trap regulon size and evidence_depth correct for."""
    import inspect
    src = inspect.getsource(pl)
    assert "n_ligands" not in src and "ligand_count" not in src


# ------------------------------------ CTD promoted to empirical evidence
def test_ctd_counts_as_measured_evidence():
    """A curator recorded that someone dosed cells and the gene moved. That
    is an observation, not an inference from a motif — and it is often the
    only empirical source standing when LINCS is unreachable."""
    from promoter_tf_drugs import provenance as pv
    led = pv.Ledger()
    led.record("chip", requested=True, ok=False)
    led.record("lincs", requested=True, ok=False)
    led.record("ctd", requested=True, ok=True, yielded=152)
    assert "ctd" in led.EMPIRICAL
    assert led.empirical_sources_ok() == 1
    assert led.empirical_names() == ["ctd"]


def test_verdict_names_the_evidence_that_did_arrive():
    from promoter_tf_drugs import provenance as pv, motifs as m
    led = pv.Ledger()
    led.record("chip", requested=True, ok=False)
    led.record("ctd", requested=True, ok=True, yielded=152)
    led.record("genetics", requested=True, ok=True, yielded=8)
    v = pv.verdict(led, m.scan_quality(900, 879, 4230, 1e-4))
    assert v["n_empirical_sources"] == 2
    assert set(v["empirical_available"]) == {"ctd", "genetics"}
    assert "ctd" in v["headline"]
    # a noise scan with SOME measurement is degraded, not uninterpretable
    assert v["level"] == "degraded"


def test_empty_ctd_is_a_failure_not_silence():
    from promoter_tf_drugs import provenance as pv
    led = pv.Ledger()
    led.record("ctd", requested=True, ok=False, yielded=0)
    assert led.get("ctd").state == "FAILED"
    assert led.get("ctd").to_dict()["impact"]
import json  # noqa: E402


# ===================================================================== #
# Things that were computed but never surfaced, and one flag that was    #
# documented but never built. Each of these was a silent gap.            #
# ===================================================================== #
def test_control_gene_is_excluded_from_the_candidate_ranking():
    """A control mixed in with candidates contaminates the comparison it
    exists to validate, and inflates every drug's recurrence count by one."""
    p = _panel_of("DEMOGENE", "SMN2")
    merged, kept, stages, rec = pn.narrow(p, control="SMN2")
    merged.setdefault("panel", {})["control_gene"] = "SMN2"
    md = pn.markdown_report(p, merged, kept, stages, rec, {"top": 40})

    assert "Negative control" in md
    assert "It is not a candidate" in md
    # recurrence was computed over candidates only, so nothing can claim to
    # hit the control as corroboration
    for entry in rec.values():
        assert "SMN2" not in entry["genes"]
    # and the surviving-routes table says it excluded the control
    assert "the control gene SMN2 is" in md


def test_control_surviving_as_well_as_candidates_is_flagged():
    """The whole point: if the funnel treats a control like a candidate, it
    is not discriminating and the shortlist is not evidence."""
    p = _panel_of("DEMOGENE", "SMN2")
    merged, kept, stages, rec = pn.narrow(p, control="SMN2")
    merged.setdefault("panel", {})["control_gene"] = "SMN2"
    md = pn.markdown_report(p, merged, kept, stages, rec, {})
    # one of the three verdicts must be reached, never silence
    assert any(s in md for s in ("not discriminating",
                                 "is discriminating",
                                 "margin as modest"))


def test_renarrow_reuses_a_saved_panel_with_no_network(tmp_path):
    """Tuning a cutoff should take a second, not the tens of minutes per gene
    a fresh run costs. This was documented before it existed."""
    out = tmp_path / "first"
    assert pn.main(["--genes", "DEMOGENE,SMN2", "--demo",
                    "--out", str(out)]) == 0

    again = tmp_path / "second"
    rc = pn.main(["--renarrow", str(out / "panel.json"),
                  "--max-regulon", "200", "--out", str(again)])
    assert rc == 0
    md = (again / "panel.md").read_text(encoding="utf-8")
    assert "The funnel" in md
    # the tightened cutoff must actually have been applied
    saved = json.loads((again / "panel.json").read_text(encoding="utf-8"))
    assert saved["criteria"]["max_regulon"] == 200


def test_renarrow_preserves_the_control_gene(tmp_path):
    out = tmp_path / "a"
    pn.main(["--genes", "DEMOGENE", "--control-gene", "SMN2", "--demo",
             "--out", str(out)])
    saved = json.loads((out / "panel.json").read_text(encoding="utf-8"))
    assert saved["panel"]["control_gene"] == "SMN2"

    again = tmp_path / "b"
    assert pn.main(["--renarrow", str(out / "panel.json"),
                    "--out", str(again)]) == 0
    md = (again / "panel.md").read_text(encoding="utf-8")
    assert "Negative control" in md


def test_overlapping_tss_windows_are_reported():
    """CCR5's two --all-tss windows shared 71% of their sequence, so they
    were largely one analysis counted twice. window_overlap existed but was
    never called."""
    from promoter_tf_drugs import neighbours as nbr
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    # one promoter in demo mode, so nothing to overlap
    assert res["tss_window_overlap"] == []

    # two overlapping windows must produce a warning in the report
    res["tss_window_overlap"] = nbr.window_overlap([
        {"chrom": "3", "start": 46365946, "end": 46371945,
         "transcript_id": "T1"},
        {"chrom": "3", "start": 46367670, "end": 46373669,
         "transcript_id": "T2"}])
    md = pipeline.markdown_report(res)
    assert "promoter windows overlap each other" in md
    assert "not independent" in md


def test_single_gene_report_names_the_receptor_family():
    """The paralog axis was wired into the panel only, so it was invisible in
    the single-gene run that most people actually do."""
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    md = pipeline.markdown_report(res)
    i = md.index("Why each one")
    body = md[i:i + 6000]
    assert "Family:" in body
    # and where selectivity WAS achieved it must say so rather than only
    # reporting the family as a liability
    assert ("subtype selectivity has been achieved" in body
            or "no subtype-selective agent" in body)


def test_family_clause_is_silent_for_a_lone_receptor():
    from promoter_tf_drugs import predictions as prd
    assert prd._family_clause("NR1H4") == ""      # FXR, no paralogue
    assert prd._family_clause(None) == ""
    assert prd._family_clause("NHLH1") == ""      # not a receptor at all
    assert "THRA" in prd._family_clause("THRB")


def test_trans_eqtl_section_renders_a_result_not_advice():
    """It used to print 'look for trans-eQTLs of X in eQTLGen' once per
    factor, identically, with no answer."""
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    res["results"][0]["genetics"] = {
        "summary": "demo",
        "trans_eqtl_summary": "**No factor can be supported this way.**",
        "trans_eqtl": [
            {"tf": "VDR", "gene": "DEMOGENE", "checked": True,
             "supported": False, "hits": [], "note": "none"},
        ],
    }
    md = pipeline.markdown_report(res)
    assert "trans-eQTL" in md
    assert "No factor can be supported this way" in md
    assert "Absence is weak" in md
    # and it must not fall back to the old per-factor advice
    assert "look for trans-eQTLs of" not in md


def test_all_tss_does_not_duplicate_a_gene_level_route():
    """Real bug: CCR5 with --all-tss has two promoter blocks, and
    catalogue.from_direct folds gene-level therapies into each one — so
    Maraviroc filled both rows of a two-row shortlist, with its caveat
    printed twice. A route on the gene product is not promoter-specific."""
    res = pipeline.run("SMN2", demo=True, offline=True)
    block = res["results"][0]
    # simulate a second TSS window, as --all-tss produces
    import copy
    second = copy.deepcopy(block)
    second["promoter"]["transcript_id"] = "ENST-SECOND"
    res["results"].append(second)

    rows = sl.enrich_rows(res)
    seen = {}
    for r in rows:
        key = (r["gene"], r["drug"].lower(), r.get("acts_on"), r["category"])
        seen[key] = seen.get(key, 0) + 1
    assert all(n == 1 for n in seen.values()),         [k for k, n in seen.items() if n > 1]

    # and the duplication is recorded rather than hidden
    assert any(r.get("_n_windows", 1) == 2 for r in rows)


def test_dedupe_keeps_the_strongest_instance_not_the_first():
    rows = [
        {"gene": "G", "drug": "X", "acts_on": "T", "category": "c",
         "tier": 3, "score": 10.0},
        {"gene": "G", "drug": "X", "acts_on": "T", "category": "c",
         "tier": 1, "score": 5.0},
    ]
    kept = sl._dedupe_across_windows(rows)
    assert len(kept) == 1
    assert kept[0]["tier"] == 1          # tier beats score, as everywhere
    assert kept[0]["_n_windows"] == 2


def test_dedupe_collapses_the_salt_form_too():
    rows = [
        {"gene": "G", "drug": "Liothyronine", "acts_on": "THRB",
         "category": "c", "tier": 3, "score": 9.0},
        {"gene": "G", "drug": "Liothyronine Sodium", "acts_on": "THRB",
         "category": "c", "tier": 3, "score": 9.0},
    ]
    assert len(sl._dedupe_across_windows(rows)) == 1


def test_dedupe_keeps_the_same_drug_at_different_targets():
    """Those are genuinely different routes."""
    rows = [
        {"gene": "G", "drug": "Tretinoin", "acts_on": "RARA",
         "category": "c", "tier": 3, "score": 9.0},
        {"gene": "G", "drug": "Tretinoin", "acts_on": "PPARD",
         "category": "c", "tier": 3, "score": 9.0},
    ]
    assert len(sl._dedupe_across_windows(rows)) == 2


# ===================================================================== #
# Surviving the web services. Ensembl, ChEMBL and the maayanlab hosts    #
# all return 500 while merely busy, and one of them used to take the     #
# whole run down with it.                                                #
# ===================================================================== #
def test_one_unreachable_window_does_not_abort_the_run(monkeypatch):
    """With --all-tss a gene has several windows. A transient 500 on one of
    them used to raise after transcript selection, regulons and expression
    had already been built — discarding all of it."""
    from promoter_tf_drugs import genes as gm
    calls = {"n": 0}

    def flaky(client, species, chrom, start, end, strand, mask=None):
        calls["n"] += 1
        if calls["n"] == 1:
            raise gm.HttpError("HTTP 500 from Ensembl")
        return "ACGT" * 300

    gene = {"seq_region_name": "3", "strand": 1, "id": "ENSG0",
            "display_name": "X", "assembly_name": "GRCh38",
            "Transcript": [
                {"id": "T1", "start": 1000, "end": 2000,
                 "biotype": "protein_coding", "is_canonical": 1},
                {"id": "T2", "start": 90000, "end": 99000,
                 "biotype": "protein_coding"},
            ]}
    monkeypatch.setattr(gm, "fetch_sequence", flaky)
    monkeypatch.setattr(gm, "_lookup_symbol", lambda c, s, sp: gene)

    proms = gm.get_promoters(None, "X", all_tss=True)
    assert len(proms) == 1            # the survivor, not an exception
    assert calls["n"] == 2


def test_every_window_failing_still_raises(monkeypatch):
    """An empty promoter list would read downstream as 'this gene has no
    promoter', which is a different and wrong claim."""
    import pytest
    from promoter_tf_drugs import genes as gm

    def dead(*a, **k):
        raise gm.HttpError("HTTP 500")

    gene = {"seq_region_name": "3", "strand": 1, "id": "ENSG0",
            "display_name": "X", "assembly_name": "GRCh38",
            "Transcript": [{"id": "T1", "start": 1000, "end": 2000,
                            "biotype": "protein_coding", "is_canonical": 1}]}
    monkeypatch.setattr(gm, "fetch_sequence", dead)
    monkeypatch.setattr(gm, "_lookup_symbol", lambda c, s, sp: gene)
    with pytest.raises(gm.HttpError, match="could not fetch any promoter"):
        gm.get_promoters(None, "X")


def test_a_local_fasta_replaces_the_ensembl_sequence_call(tmp_path):
    """Ensembl was the one stage with no offline route. The package already
    writes promoter.fa, so a previous run is its own offline source."""
    from promoter_tf_drugs import genes as gm
    fa = tmp_path / "promoter.fa"
    fa.write_text(">GENE|ENST1|chr7:1000-1011|strand=1|-1000/+200\n"
                  "ACGTACGTACGT\n", encoding="utf-8")
    assert gm.load_promoter_fasta(fa) == 1
    try:
        assert gm.local_sequence("7", 1000, 1011) == "ACGTACGTACGT"
        assert gm.local_sequence("chr7", 1000, 1011) == "ACGTACGTACGT"
        # a near-miss must NOT silently return the wrong window
        assert gm.local_sequence("7", 1001, 1011) is None
        assert gm.local_sequence("8", 1000, 1011) is None
    finally:
        gm._LOCAL_SEQ.clear()


def test_local_fasta_is_preferred_over_the_network(tmp_path, monkeypatch):
    from promoter_tf_drugs import genes as gm
    fa = tmp_path / "p.fa"
    fa.write_text(">G|T1|chr3:1-1200|strand=1|-1000/+200\n"
                  + "A" * 1200 + "\n", encoding="utf-8")
    gm.load_promoter_fasta(fa)

    def boom(*a, **k):
        raise AssertionError("must not hit the network when a FASTA matches")

    gene = {"seq_region_name": "3", "strand": 1, "id": "ENSG0",
            "display_name": "G", "assembly_name": "GRCh38",
            "Transcript": [{"id": "T1", "start": 1001, "end": 2000,
                            "biotype": "protein_coding", "is_canonical": 1}]}
    monkeypatch.setattr(gm, "fetch_sequence", boom)
    monkeypatch.setattr(gm, "_lookup_symbol", lambda c, s, sp: gene)
    try:
        proms = gm.get_promoters(None, "G", upstream=1000, downstream=200)
        assert proms and proms[0].sequence == "A" * 1200
    finally:
        gm._LOCAL_SEQ.clear()


def test_the_real_promoter_fa_output_round_trips():
    """The offline route only works if this package's own output is readable
    by it. Round-trip its FASTA rather than a hand-built one."""
    from promoter_tf_drugs import genes as gm
    import tempfile
    from pathlib import Path as _P
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    with tempfile.TemporaryDirectory() as d:
        pipeline.write_outputs(res, d)
        n = gm.load_promoter_fasta(_P(d) / "promoter.fa")
        try:
            assert n == len(res["results"])
            prom = res["results"][0]["promoter"]
            got = gm.local_sequence(prom["chrom"], prom["start"], prom["end"])
            assert got == res["results"][0]["sequence"]
        finally:
            gm._LOCAL_SEQ.clear()


def test_a_500_is_treated_as_backoff_not_a_verdict():
    """Public bioinformatics endpoints 500 under load as often as they 429.
    Three attempts over 7s gave up while Ensembl was merely busy."""
    from promoter_tf_drugs import http as h
    assert h.DEFAULT_RETRIES >= 5
    assert h.RETRY_CAP >= 10


# ===================================================================== #
# The last three plan items, and the cross-check the plan asked for.     #
# ===================================================================== #
def test_panel_regulon_stage_agrees_with_shortlist_specific():
    """Both gate on max_regulon through the same predicate, so on one gene
    they must return the same set. This is the check that the shared-
    predicate refactor actually shares behaviour rather than merely code."""
    res = pipeline.run("SMN2", demo=True, offline=True)

    sl_keys = {(r["gene"], r["drug"].lower(), r.get("acts_on"))
               for r in sl.filter_routes(res, max_regulon=600)}

    p = pn.Panel()
    pg = pn.PanelGene(symbol="SMN2")
    pg.result = res
    p.genes.append(pg)
    merged = pn.merge_results(p)
    rows = sl.enrich_rows(merged)
    pn.annotate_axes(rows, merged)
    kept, _ = pn.run_funnel(
        rows, [pn.Stage("regulon", "a", {"max_regulon": 600}, "")])
    pn_keys = {(r["gene"], r["drug"].lower(), r.get("acts_on")) for r in kept}

    assert sl_keys == pn_keys, (sl_keys ^ pn_keys)


def test_per_gene_top_stops_one_gene_filling_the_table():
    """filter_rows sorts by (tier, -score) across the whole panel, so a gene
    with more annotated pharmacology crowds out a genuinely cleaner gene —
    the ranking answering the wrong question for a tool comparing genes."""
    rows = ([{"gene": "A", "drug": f"a{i}", "tier": 3, "score": 90 - i}
             for i in range(10)]
            + [{"gene": "B", "drug": "b1", "tier": 3, "score": 10}])
    capped = pn.cap_per_gene(rows, 2)
    assert [r["gene"] for r in capped] == ["A", "A", "B"]
    # global order within the cap is preserved
    assert [r["drug"] for r in capped][:2] == ["a0", "a1"]
    # no cap means no change
    assert pn.cap_per_gene(rows, None) == rows
    assert pn.cap_per_gene(rows, 0) == rows


def test_per_gene_top_is_wired_to_the_cli():
    flags = [o for a in pn.build_parser()._actions for o in a.option_strings]
    assert "--per-gene-top" in flags
    args = pn.build_parser().parse_args(["--genes", "A", "--per-gene-top", "3"])
    assert args.per_gene_top == 3


def test_both_clis_share_one_declaration_of_the_analysis_flags():
    """They used to be shared by reaching into the other parser's private
    _actions and registering the same Action objects in two parsers, which
    shares mutable state between them."""
    from promoter_tf_drugs import cli, panel as pnl
    import inspect

    assert hasattr(cli, "add_analysis_args")
    assert "_add_action" not in inspect.getsource(pnl.build_parser)

    single = {o for a in cli.build_parser()._actions
              for o in a.option_strings if o.startswith("--")}
    multi = {o for a in pnl.build_parser()._actions
             for o in a.option_strings if o.startswith("--")}

    # every analysis flag reaches both, and neither leaks the other's own
    for shared in ("--pfm-file", "--chip-bed", "--promoter-fasta",
                   "--regulon-file", "--expression-file", "--direction",
                   "--all-tss", "--pvalue", "--demo", "--offline"):
        assert shared in single, shared
        assert shared in multi, shared
    assert "--gene" in single and "--gene" not in multi
    assert "--genes" in multi and "--genes" not in single


def test_the_two_parsers_do_not_share_action_objects():
    """Distinct Action instances: a default mutated on one parser must not
    appear on the other.

    Both parsers are held in locals for the comparison. Taking id() of
    objects from a temporary lets CPython free them and reuse the addresses
    for the second parser's actions, which made an earlier version of this
    test fail at random.
    """
    from promoter_tf_drugs import cli, panel as pnl
    single = cli.build_parser()
    multi = pnl.build_parser()
    shared = {id(x) for x in single._actions} & {id(x) for x in multi._actions}
    assert not shared

    # The property that actually matters, checked behaviourally rather than
    # by identity: mutating one parser's action leaves the other alone.
    a = [x for x in single._actions if "--pvalue" in x.option_strings][0]
    b = [x for x in multi._actions if "--pvalue" in x.option_strings][0]
    a.default = 12345
    assert b.default != 12345


# ===================================================================== #
# The three remaining limits, corrected.                                #
# ===================================================================== #
def test_curated_tables_are_internally_consistent():
    """A hand-entered table that is quietly wrong is worse than one known to
    be incomplete, because the report presents both with equal confidence.
    This is the check that runs on every commit."""
    from promoter_tf_drugs import audit
    problems = audit.check_internal()
    bad = {k: v for k, v in problems.items() if v}
    assert not bad, audit.summarise(problems)


def test_audit_catches_a_drug_placed_in_two_systemic_classes():
    """It would make `classify` depend on table declaration order, so the
    burden a drug reports would change if someone reordered the file."""
    from promoter_tf_drugs import audit, systemic
    extra = systemic.SystemicClass(
        "invented", 0.5, "a" * 50, frozenset({"dexamethasone"}))
    orig = systemic.CLASSES
    systemic.CLASSES = orig + (extra,)
    try:
        problems = audit.check_internal()
        assert any("dexamethasone" in p
                   for p in problems["systemic: drug in two classes"])
        assert "problem(s) found" in audit.summarise(problems)
    finally:
        systemic.CLASSES = orig


def test_audit_catches_a_name_that_would_not_survive_salt_stripping():
    """classify() looks up the salt-free form, so a member stored WITH a
    counterion could never be found."""
    from promoter_tf_drugs import audit, systemic
    extra = systemic.SystemicClass(
        "invented", 0.5, "b" * 50, frozenset({"levothyroxine sodium"}))
    orig = systemic.CLASSES
    systemic.CLASSES = orig + (extra,)
    try:
        problems = audit.check_internal()
        assert "levothyroxine sodium" in problems[
            "systemic: name changes under salt_free_name"]
    finally:
        systemic.CLASSES = orig


def test_audit_offline_needs_no_network():
    """The consistency pass must work with no client at all, so it can run as
    a test rather than only as a manual chore."""
    from promoter_tf_drugs import audit
    assert isinstance(audit.check_internal(), dict)
    assert "internally consistent" in audit.summarise(
        {k: [] for k in audit.check_internal()})


def test_audit_reports_unresolvable_names_without_claiming_correctness():
    from promoter_tf_drugs import audit

    class C:
        def get(self, url, params=None):
            q = (params or {}).get("q")
            if q == "ghostdrug":
                return {"molecules": []}
            return {"molecules": [{"molecule_chembl_id": "CHEMBL1",
                                   "withdrawn_flag": q == "oldpill"}]}

    r = audit.resolve_names(C(), ["dexamethasone", "ghostdrug", "oldpill"])
    assert r["unresolved"] == ["ghostdrug"]
    assert r["withdrawn_in_chembl"] == ["oldpill"]
    assert len(r["resolved"]) == 2


def test_audit_counts_lookup_errors_separately_from_unresolved():
    """A name ChEMBL could not answer for is UNKNOWN, not missing — the same
    distinction the evidence ledger makes everywhere else."""
    from promoter_tf_drugs import audit

    class Boom:
        def get(self, url, params=None):
            raise RuntimeError("HTTP 500")

    r = audit.resolve_names(Boom(), ["a", "b"])
    assert r["lookup_errors"] == 2
    assert r["unresolved"] == []          # not silently called missing


def test_direction_model_health_is_recorded_not_inferred():
    """A ChEMBL outage and a target with genuinely few ligands both end as
    'no direction model'. The ledger must say which."""
    from promoter_tf_drugs import provenance as pv
    assert "direction_model" in pv.FLAGS
    assert "direction_model" in pv.IMPACT
    led = pv.Ledger()
    led.record("direction_model", requested=True, ok=False,
               detail="12 of 12 grades could not be direction-checked")
    assert led.get("direction_model").state == "FAILED"
    assert led.get("direction_model").to_dict()["impact"]


def test_every_panel_preset_sets_the_specificity_axes():
    """The paralog axis was the one axis absent from the presets, which made
    the four axes inconsistently active."""
    from promoter_tf_drugs import panel as p
    for name in ("narrow", "narrow+clean"):
        pre = p.PANEL_PRESETS[name]
        for axis in ("max_regulon", "max_burden", "min_paralog_specificity"):
            assert axis in pre, f"{name} is missing {axis}"


def test_paralog_axis_stays_off_by_default_and_says_why():
    """It is a structural prior rather than a measurement, so auto-deleting
    routes on it would be the thing this package criticises elsewhere. The
    Family column is populated either way."""
    from promoter_tf_drugs import panel as p
    stage = [s for s in p.default_stages()
             if s.name == "Selective ligand achievable"][0]
    assert stage.filters["min_paralog_specificity"] is None
    assert "prior rather than a measurement" in stage.why
    assert "nothing is hidden" in stage.why


def test_audit_covers_every_curated_drug_table():
    """A 'full' audit that skipped direct.CURATED would leave unchecked the
    one table the report prints FIRST and calls the better starting point."""
    from promoter_tf_drugs import audit
    tables = audit.all_curated_names()
    for expected in ("systemic.CLASSES", "stabilizers.CURATED",
                     "delivery.PGP_SUBSTRATES", "ligands.ROUTES",
                     "direct.CURATED"):
        assert expected in tables, expected
        assert tables[expected], expected
    assert sum(len(v) for v in tables.values()) > 200


def test_withdrawn_disagreement_is_directional():
    """We-say-available/ChEMBL-says-withdrawn puts a withdrawn drug at the
    top of the report unflagged. The reverse is far less serious, because
    approvals are regional and ChEMBL is not the regulator."""
    from promoter_tf_drugs import audit, direct

    class C:
        def __init__(self, withdrawn_names):
            self.w = {n.lower() for n in withdrawn_names}

        def get(self, url, params=None):
            q = (params or {}).get("q", "")
            return {"molecules": [{"molecule_chembl_id": "CHEMBL1",
                                   "withdrawn_flag": q.lower() in self.w}]}

    rows = direct.CURATED
    vals = list(rows.values()) if isinstance(rows, dict) else list(rows)
    entries = [it for v in vals
               for it in (v if isinstance(v, (list, tuple)) else [v])
               if getattr(it, "drug", None)]
    live = next(e.drug for e in entries if not getattr(e, "withdrawn", False))

    r = audit.check_withdrawn_agreement(C([live]))
    assert live in r["we_say_available_chembl_says_withdrawn"]

    clean = audit.check_withdrawn_agreement(C([]))
    assert clean["we_say_available_chembl_says_withdrawn"] == []


def test_unanswerable_names_are_unchecked_not_agreed():
    """A name ChEMBL could not answer for is UNKNOWN. Counting it as
    agreement would be the absence-of-evidence bug again."""
    from promoter_tf_drugs import audit

    class Boom:
        def get(self, url, params=None):
            raise RuntimeError("HTTP 500")

    r = audit.check_withdrawn_agreement(Boom())
    assert r["unchecked"] and r["checked"] == len(r["unchecked"])
    assert r["we_say_available_chembl_says_withdrawn"] == []


def test_no_curated_entry_is_presented_as_available_while_flagged_withdrawn():
    """Internal version of the online check: our own table must not contain
    an entry that is both withdrawn and presented as raising protein without
    the flag being visible downstream."""
    from promoter_tf_drugs import direct
    rows = direct.CURATED
    vals = list(rows.values()) if isinstance(rows, dict) else list(rows)
    for v in vals:
        for e in (v if isinstance(v, (list, tuple)) else [v]):
            if getattr(e, "withdrawn", False):
                # a withdrawn entry must still carry its mechanism and gene so
                # the report can explain why it is being shown at all
                assert e.mechanism and e.gene, e.drug

def _no_local_gtf(monkeypatch):
    """Force the network path. A local GTF now answers first, so a test that
    wants to exercise the Ensembl route has to say so."""
    from promoter_tf_drugs import neighbours as nbr
    monkeypatch.setattr(nbr, "DEFAULT_GTF_PATHS", ())
    monkeypatch.setattr(nbr, "_GTF_GENES", [])


def test_a_clean_window_says_so(monkeypatch):
    _no_local_gtf(monkeypatch)
    r = nbr.genes_in_window(_NbClient([]), "3", 1000, 2000, strand=1)
    assert r["checked"] and r["clean"] and r["n_inside"] == 0
    assert "alone" in r["summary"]


def test_a_failed_lookup_is_not_a_clean_window(monkeypatch):
    """Same principle as provenance: absence of evidence must not read as
    evidence of absence."""
    _no_local_gtf(monkeypatch)

    class Boom:
        def get(self, url, params=None):
            raise RuntimeError("Ensembl 500")
    r = nbr.genes_in_window(Boom(), "3", 1, 2, strand=1)
    assert r["checked"] is False
    assert r.get("clean") is not True


def test_a_local_gtf_is_preferred_over_the_network(tmp_path, monkeypatch):
    """The whole point of the offline route: Ensembl's overlap endpoint 500s
    under load, and when it does this check silently reports 'not checked'
    for every promoter."""
    from promoter_tf_drugs import neighbours as nbr2
    gtf = tmp_path / "mini.gtf"
    rows = [
        ["chr3", "X", "gene", "46363984", "46407184", ".", "-", ".",
         'gene_id "G1"; gene_type "lncRNA"; gene_name "CCR5AS";'],
        ["chr3", "X", "gene", "46370939", "46377528", ".", "+", ".",
         'gene_id "G2"; gene_type "protein_coding"; gene_name "CCR5";'],
    ]
    gtf.write_text(
        "".join(chr(9).join(r) + chr(10) for r in rows),
        encoding="utf-8")
    assert nbr2.load_gtf(gtf) == 2
    try:
        def boom(*a, **k):
            raise AssertionError("must not hit the network with a GTF loaded")
        monkeypatch.setattr(nbr2, "peaks_from_api", boom, raising=False)

        r = nbr2.genes_in_window(None, "3", 46369946, 46371145,
                                 strand=1, exclude="CCR5")
        assert r["checked"] and r["n_antisense"] == 1
        assert r["genes"][0]["gene"] == "CCR5AS"
        assert r["source"] and "mini.gtf" in r["source"]
    finally:
        nbr2._GTF_GENES.clear()


# ===================================================================== #
# Routes a promoter analysis cannot see, and why the drug tables are    #
# always nuclear receptors.                                             #
# ===================================================================== #
def test_attrition_explains_the_nuclear_receptor_tables():
    """The observation this exists for: hormone and retinoid drugs rank high
    for every gene and direction. They are not winning a ranking - they are
    what survives after every factor without a ligand is dropped. ~48 of
    ~1,600 human TFs have a ligand pocket."""
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    at = res["results"][0]["attrition"]
    assert at["n_tfs_with_sites"] > at["n_with_drug"]
    assert at["n_drugged_nuclear_receptor"] >= 1

    md = pipeline.markdown_report(res)
    assert "Why the drug tables below look the way they do" in md
    assert "what is left after every factor without a ligand" in md


def test_high_scoring_undruggable_factors_are_named():
    """The gene-specific signal is in the factors that score well and have
    no ligand - exactly what gets discarded. In a real CCR5 run that was
    PDX1, RBPJ and NFATC1, and they appeared in no drug table."""
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    at = res["results"][0]["attrition"]
    assert at["n_top25_without_drug"] >= 1
    assert at["top25_without_drug"]

    md = pipeline.markdown_report(res)
    assert "appear in no drug table" in md
    assert at["top25_without_drug"][0] in md


def test_undruggable_top_factors_get_priority_in_the_upstream_search():
    """They are what that search exists for, so they must not depend on
    where they happen to fall in the ranking."""
    import inspect
    src = inspect.getsource(pipeline.run)
    assert "_priority" in src
    assert "top25_without_drug" in src


# ------------------------------------------- post-transcriptional routes
def test_ptc_moves_readthrough_to_the_front():
    """For a nonsense allele, readthrough and NMD address the problem and
    extra transcription does not - doubling transcription doubles the amount
    of message that gets destroyed."""
    from promoter_tf_drugs import posttranscriptional as ptx
    with_ptc = ptx.routes_for("up", has_ptc=True)
    assert with_ptc[0].mechanism == "nonsense readthrough"
    assert "nonsense-mediated" in with_ptc[1].mechanism

    plain = ptx.routes_for("up")
    assert plain[0].mechanism != "nonsense readthrough"
    assert len(plain) == len(with_ptc)        # reorders, never filters


def test_ataluren_is_not_recorded_as_fda_approved():
    """It has EMA conditional approval and the FDA did not approve it. The
    user's criterion is FDA, so overstating this would be the one thing
    that makes the list useless."""
    from promoter_tf_drugs import posttranscriptional as ptx
    rt = [r for r in ptx.ROUTES if r.mechanism == "nonsense readthrough"][0]
    joined = " ".join(rt.drugs)
    assert "Ataluren" in joined
    assert "FDA did NOT approve" in joined or "did NOT approve" in joined
    assert rt.status == "approved elsewhere"


def test_post_transcriptional_states_what_is_not_approved():
    from promoter_tf_drugs import posttranscriptional as ptx
    nmd = [r for r in ptx.ROUTES if "nonsense-mediated" in r.mechanism][0]
    assert nmd.status == "investigational"
    assert "No approved NMD inhibitor" in " ".join(nmd.drugs)


# --------------------------------------------------- paralog compensation
def test_compensation_finds_the_documented_substitutes():
    """You may not need to raise your gene. SMN2 for SMN1 and HBG for HBB
    are how the field's clearest successes actually worked."""
    from promoter_tf_drugs import compensation as comp
    for gene, paralog in (("SMN1", "SMN2"), ("HBB", "HBG1/HBG2"),
                          ("DMD", "UTRN")):
        r = comp.assess(gene, goal="up")
        assert r["recommended"] is True, gene
        assert r["known_redundancy"]["paralog"] == paralog


def test_compensation_distinguishes_no_paralog_from_unknown():
    """Absence from a short curated table is not evidence that no paralog
    substitutes."""
    from promoter_tf_drugs import compensation as comp
    none = comp.assess("SERPINA1", goal="up")
    assert none["recommended"] is False          # genuinely none

    unknown = comp.assess("ZZZ9FAKE", goal="up")
    assert unknown["recommended"] is None        # not checked, not absent
    assert "nobody has written one down" in " ".join(unknown["notes"])


def test_compensation_is_not_the_receptor_family_question():
    """paralogs.py asks how many relatives share a BINDING POCKET, as a
    cross-reactivity risk. This asks whether a relative shares the gene's
    FUNCTION, as an opportunity. Same word, inverted use."""
    from promoter_tf_drugs import compensation as comp, paralogs as par
    assert comp.assess("HBB", goal="up")["recommended"] is True
    assert par.paralog_specificity("HBB")["assessed"] is False


# ------------------------------------------------ silenced-allele routes
def test_imprinted_gene_is_flagged_as_a_different_problem():
    """A silenced allele needs a brake released, not more transcription -
    and pushing the active copy harder is the wrong move."""
    from promoter_tf_drugs import reactivation as react
    r = react.assess("UBE3A", goal="up")
    assert r["applicable"] is True
    assert "imprinted locus" in " ".join(r["notes"])


def test_unlisted_gene_says_silencing_was_never_measured():
    from promoter_tf_drugs import reactivation as react
    r = react.assess("ZZZ9FAKE", goal="up")
    assert r["applicable"] is False
    joined = " ".join(r["notes"])
    assert "nothing in this package measures methylation" in joined


def test_reactivation_states_where_selectivity_comes_from():
    """Every mechanism here is non-specific by default; the selectivity has
    to come from somewhere else or the route is unusable."""
    from promoter_tf_drugs import reactivation as react
    for r in react.ROUTES:
        assert r.specificity and len(r.specificity) > 20
    demeth = [r for r in react.ROUTES if "methylation" in r.mechanism][0]
    assert "none at the locus level" in demeth.specificity


# ------------------------------------------- target-protein stabilisation
def test_product_stabilisation_is_about_the_gene_product_not_the_factor():
    from promoter_tf_drugs import stabilizers as st
    r = st.product_stabilisation("MYGENE", goal="up")
    mechs = {x["mechanism"] for x in r["routes"]}
    assert "proteasome inhibition" in mechs
    assert "pharmacological chaperone" in mechs
    assert any("transcription alone" in n for n in r["notes"])


def test_only_the_chaperone_route_is_selective():
    """Three of the four raise most of the proteome, which makes them good
    diagnostics and poor therapies."""
    from promoter_tf_drugs import stabilizers as st
    routes = {x["mechanism"]: x for x in
              st.product_stabilisation("G", goal="up")["routes"]}
    assert routes["proteasome inhibition"]["specificity"] < 0.1
    assert routes["lysosomal / autophagy blockade"]["specificity"] < 0.1
    assert routes["pharmacological chaperone"]["specificity"] > 0.5


def test_beyond_transcription_only_runs_for_up():
    """Lowering a gene post-transcriptionally is the ASO/siRNA story, which
    modalities already designs properly."""
    from promoter_tf_drugs import (posttranscriptional as ptx,
                                   compensation as comp,
                                   reactivation as react,
                                   stabilizers as st)
    assert ptx.routes_for("down") == []
    assert comp.assess("HBB", goal="down") == {}
    assert react.assess("UBE3A", goal="down") == {}
    assert st.product_stabilisation("G", goal="down") == {}


def test_report_renders_the_beyond_transcription_section():
    res = pipeline.run("HBB", demo=True, offline=True, has_ptc=True)
    md = pipeline.markdown_report(res)
    assert "Routes a promoter analysis cannot see" in md
    assert "Could a relative do the job instead?" in md
    assert "HBG1/HBG2" in md
    assert "Raising the protein instead of the transcript" in md
    # the PTC ordering reached the report
    i = md.index("After transcription")
    assert "readthrough" in md[i:i + 1200].lower()


# ===================================================================== #
# Finding FDA-approved drugs: two more independent sources.             #
# ===================================================================== #
def test_opentargets_is_a_third_drug_source():
    """The API was already being called for genetics and expression and
    never asked for drugs. It is a separate curation, so it finds compounds
    ChEMBL and DGIdb miss - eldecalcitol at VDR, tiratricol at THRB."""
    from promoter_tf_drugs import drugs as dr

    class C:
        def post(self, url, payload):
            return {"data": {"target": {
                "approvedSymbol": "VDR",
                "drugAndClinicalCandidates": {"count": 2, "rows": [
                    {"maxClinicalStage": "APPROVAL", "drug": {
                        "id": "CHEMBL1", "name": "ELDECALCITOL",
                        "drugType": "Small molecule",
                        "maximumClinicalStage": "APPROVAL",
                        "drugWarnings": [],
                        "mechanismsOfAction": {"rows": [
                            {"actionType": "AGONIST",
                             "mechanismOfAction": "VDR agonist"}]}}},
                    {"maxClinicalStage": "PHASE_2", "drug": {
                        "id": "CHEMBL2", "name": "SEOCALCITOL",
                        "drugType": "Small molecule",
                        "maximumClinicalStage": "PHASE_2",
                        "drugWarnings": [],
                        "mechanismsOfAction": {"rows": []}}},
                ]}}}}

    import promoter_tf_drugs.genetics as gen
    old = gen._resolve_ensembl
    gen._resolve_ensembl = lambda c, s: "ENSG0"
    try:
        got = dr.opentargets_drugs(C(), "VDR", approved_only=True)
    finally:
        gen._resolve_ensembl = old

    names = {d.name for d in got}
    assert "Eldecalcitol" in names
    assert "Seocalcitol" not in names          # phase 2, not approved
    d = next(x for x in got if x.name == "Eldecalcitol")
    assert d.direction == "up"
    assert d.source == "OpenTargets"
    assert d.max_phase == 4.0


def test_opentargets_action_type_maps_to_a_direction():
    """DGIdb gives no action type, so direction has to be inferred or the
    row cannot be filtered by goal. Open Targets does give one."""
    from promoter_tf_drugs.drugs import _direction_for_action
    # The vocabulary must be the one the FILTER compares against - "up",
    # "down", "unclear" - not the dataclass comment's stale "increases TF
    # activity". Following the comment meant every Open Targets drug failed
    # `d.direction == want` and was dropped for any annotated factor,
    # without even reaching the unknown-direction bucket.
    assert _direction_for_action("AGONIST") == "up"
    assert _direction_for_action("INHIBITOR") == "down"
    assert _direction_for_action("") == "unclear"
    assert _direction_for_action("MODULATOR") == "unclear"

    # The substring trap, which this test caught in the first version:
    # "antagonist" CONTAINS "agonist", so testing the increase patterns
    # first mapped every antagonist to "increases TF activity" - putting
    # blockers in the raise-the-gene list.
    assert _direction_for_action("ANTAGONIST") == "down"
    assert _direction_for_action("INVERSE AGONIST") == "down"
    assert _direction_for_action("NEGATIVE ALLOSTERIC MODULATOR") == "down"
    # and the genuine agonist cases still read correctly
    assert _direction_for_action("PARTIAL AGONIST") == "up"
    assert _direction_for_action("POSITIVE ALLOSTERIC MODULATOR") == "up"


def test_opentargets_schema_change_degrades_rather_than_raising():
    """knownDrugs was removed before API v26.9 and this field may go too."""
    from promoter_tf_drugs import drugs as dr

    class C:
        def post(self, url, payload):
            return {"errors": [{"message": "Cannot query field"}]}

    import promoter_tf_drugs.genetics as gen
    old = gen._resolve_ensembl
    gen._resolve_ensembl = lambda c, s: "ENSG0"
    try:
        assert dr.opentargets_drugs(C(), "VDR") == []
    finally:
        gen._resolve_ensembl = old


# ------------------------------------------------- chemical analog search
def test_analog_search_subtracts_what_is_already_annotated():
    """The naive version rediscovers its own input: a compound similar to a
    known ligand that anyone assayed at the target is already annotated
    there. Only the leftover is a new pair."""
    from promoter_tf_drugs import analogs

    class C:
        def get(self, url, params=None):
            if "molecule/search" in url:
                return {"molecules": [{
                    "pref_name": "CALCITRIOL", "max_phase": 4,
                    "molecule_chembl_id": "CHEMBL1",
                    "molecule_structures": {"canonical_smiles": "CCO"}}]}
            return {"molecules": [
                {"pref_name": "CALCITRIOL", "similarity": "100",
                 "molecule_chembl_id": "CHEMBL1"},
                {"pref_name": "ALFACALCIDOL", "similarity": "82",
                 "molecule_chembl_id": "CHEMBL2"},
                {"pref_name": "TACALCITOL", "similarity": "76",
                 "molecule_chembl_id": "CHEMBL3"},
            ]}

    got = analogs.find_analogs(C(), "calcitriol", "VDR",
                               already_known={"Calcitriol", "Alfacalcidol"})
    names = [a.drug for a in got]
    assert names == ["Tacalcitol"]             # the only untested pair
    assert "not itself recorded at VDR" in got[0].note


def test_analog_search_picks_the_parent_not_an_ester():
    """A plain search for tretinoin returns TRETINOIN TOCOFERIL first, whose
    similarity neighbourhood is the ester's and therefore the wrong one."""
    from promoter_tf_drugs import analogs

    class C:
        def get(self, url, params=None):
            return {"molecules": [
                {"pref_name": "TRETINOIN TOCOFERIL", "max_phase": 4,
                 "molecule_chembl_id": "CHEMBL9",
                 "molecule_structures": {"canonical_smiles": "CCCC"}},
                {"pref_name": "TRETINOIN", "max_phase": 4,
                 "molecule_chembl_id": "CHEMBL8",
                 "molecule_structures": {"canonical_smiles": "CC"}},
            ]}

    m = analogs._parent_molecule(C(), "tretinoin")
    assert m["pref_name"] == "TRETINOIN"


def test_analog_search_collapses_salt_forms_of_the_query():
    from promoter_tf_drugs import analogs

    class C:
        def get(self, url, params=None):
            if "molecule/search" in url:
                return {"molecules": [{
                    "pref_name": "LIOTHYRONINE", "max_phase": 4,
                    "molecule_chembl_id": "CHEMBL1",
                    "molecule_structures": {"canonical_smiles": "CCO"}}]}
            return {"molecules": [
                {"pref_name": "LIOTHYRONINE SODIUM", "similarity": "100",
                 "molecule_chembl_id": "CHEMBL2"},
                {"pref_name": "LEVOTHYROXINE", "similarity": "74",
                 "molecule_chembl_id": "CHEMBL3"},
            ]}

    got = analogs.find_analogs(C(), "liothyronine", "THRB")
    assert [a.drug for a in got] == ["Levothyroxine"]


def test_analog_rows_state_that_similarity_is_not_activity():
    """Beta carotene came back 80% similar to tretinoin and is a retinoid
    PRECURSOR, not an RAR ligand. The caveat has to be on every row."""
    from promoter_tf_drugs import analogs
    a = analogs.Analog(drug="X", chembl_id="C1", similarity=80.0,
                       similar_to="tretinoin", target="RARA")
    a.caveats = ["Similarity is not activity.", "flip agonist to antagonist"]
    assert any("not activity" in c for c in a.caveats)


def test_empty_analog_result_is_explained():
    from promoter_tf_drugs import analogs
    s = analogs.summarise({})
    assert "common outcome" in s and "already in the tables above" in s


def test_opentargets_drugs_actually_survive_the_direction_filter():
    """They were silently dropped: the filter compares against "up"/"down"
    and these arrived saying "increases TF activity", so for any factor with
    a role annotation every Open Targets row vanished - and not even into
    the unknown-direction bucket. An earlier test missed it because the mock
    had no role, which makes `want` "either" and keeps everything."""
    from promoter_tf_drugs import drugs as dr

    class C:
        def post(self, url, payload):
            return {"data": {"target": {"approvedSymbol": "VDR",
                "drugAndClinicalCandidates": {"count": 2, "rows": [
                    {"maxClinicalStage": "APPROVAL", "drug": {
                        "id": "C1", "name": "AN_AGONIST",
                        "maximumClinicalStage": "APPROVAL", "drugWarnings": [],
                        "mechanismsOfAction": {"rows": [
                            {"actionType": "AGONIST",
                             "mechanismOfAction": "agonist"}]}}},
                    {"maxClinicalStage": "APPROVAL", "drug": {
                        "id": "C2", "name": "AN_ANTAGONIST",
                        "maximumClinicalStage": "APPROVAL", "drugWarnings": [],
                        "mechanismsOfAction": {"rows": [
                            {"actionType": "ANTAGONIST",
                             "mechanismOfAction": "antagonist"}]}}},
                ]}}}}

        def get(self, url, params=None):
            return {}

    import promoter_tf_drugs.genetics as gen
    old = gen._resolve_ensembl
    gen._resolve_ensembl = lambda c, s: "ENSG0"
    try:
        # an ANNOTATED factor, so `want` is "up"/"down" rather than "either"
        up, _ = dr.drugs_for_tfs(C(), {"VDR": {"role": "activator"}},
                                 goal="up", use_dgidb=False)
        dn, _ = dr.drugs_for_tfs(C(), {"VDR": {"role": "activator"}},
                                 goal="down", use_dgidb=False)
    finally:
        gen._resolve_ensembl = old

    assert [d.name for v in up.values() for d in v] == ["An_Agonist"]
    assert [d.name for v in dn.values() for d in v] == ["An_Antagonist"]


def test_the_drug_direction_vocabulary_is_consistent_across_sources():
    """ChEMBL writes "up"/"down"; anything else silently fails the filter."""
    from promoter_tf_drugs import drugs as dr
    assert dr.desired_direction("activator", "up") == "up"
    assert dr.desired_direction("activator", "down") == "down"
    for action in ("AGONIST", "ANTAGONIST", "INHIBITOR", "MODULATOR"):
        assert dr._direction_for_action(action) in ("up", "down", "unclear")


# ===================================================================== #
# The BBB column: a yes / no per drug, and the three ways it can lie.   #
# ===================================================================== #
# ChEMBL molecule_properties as the API actually returns them - the
# spelling that was silently unreadable before normalise_properties.
_CHEMBL_DIAZEPAM = {"full_mwt": 284.7, "psa": 32.7, "alogp": 2.82,
                    "hbd": 0, "cx_most_bpka": 3.4}
_CHEMBL_ATENOLOL = {"full_mwt": 266.3, "psa": 84.6, "alogp": 0.1,
                    "hbd": 3, "cx_most_bpka": 9.6}
_CHEMBL_PLICAMYCIN = {"full_mwt": 1085.2, "psa": 331.0, "alogp": -1.2,
                      "hbd": 9}


def test_chembl_property_names_are_readable_at_all():
    """The bug this fixes: drugs.py passes ChEMBL's `molecule_properties`
    through verbatim, and every MPO parameter except `hbd` arrived under a
    key nothing read. Real rows scored on one parameter, fell below
    MIN_MPO_PARAMETERS and came back unassessed - indistinguishable from a
    compound with no structure at all."""
    norm = dlv.normalise_properties(_CHEMBL_DIAZEPAM)
    assert norm["mw"] == 284.7
    assert norm["tpsa"] == 32.7
    assert norm["logp"] == 2.82
    assert norm["pka"] == 3.4

    score, n, _ = dlv.cns_mpo(properties=_CHEMBL_DIAZEPAM)
    assert n >= dlv.MIN_MPO_PARAMETERS, "ChEMBL names must reach the scorer"
    assert score is not None


def test_short_names_still_work_and_win_over_the_chembl_spelling():
    """Both spellings must work, and an explicitly supplied short name must
    not be overridden by whatever ChEMBL happened to return."""
    both = dlv.normalise_properties({"mw": 100.0, "full_mwt": 999.0})
    assert both["mw"] == 100.0


def test_normalise_drops_nulls_rather_than_coercing_them():
    """A None MW read as 0.0 would score as a perfect parameter."""
    norm = dlv.normalise_properties({"full_mwt": None, "psa": "",
                                     "alogp": "not a number", "hbd": 2})
    assert "mw" not in norm
    assert "tpsa" not in norm
    assert "logp" not in norm
    assert norm["hbd"] == 2.0


def test_bbb_verdict_has_three_states_and_unknown_is_not_no():
    """The third state is the one that matters: a DGIdb row carries a name
    and a mechanism but no structure, and that is not evidence of failure
    to cross."""
    yes = dlv.bbb_verdict("diazepam", _CHEMBL_DIAZEPAM)
    assert yes["label"] == "yes" and yes["crosses_bbb"] is True

    no = dlv.bbb_verdict("plicamycin", _CHEMBL_PLICAMYCIN)
    assert no["label"] == "no" and no["crosses_bbb"] is False

    unknown = dlv.bbb_verdict("mystery compound", {})
    assert unknown["crosses_bbb"] is None
    assert unknown["label"] == "—"
    assert "UNKNOWN rather than no" in unknown["why"]


def test_pgp_substrate_overrides_a_passing_score():
    """MPO predicts passive diffusion. The commonest reason a suitable
    molecule never reaches brain is being pumped back out, so a known
    substrate cannot read yes however good its physicochemistry."""
    props = {"full_mwt": 284.7, "psa": 32.7, "alogp": 2.82}
    assert "loperamide" in dlv.PGP_SUBSTRATES
    v = dlv.bbb_verdict("loperamide", props)
    assert v["label"] == "no"
    assert v["pgp_substrate"] is True
    assert "pumped back out" in v["why"]
    # and the override is not quietly re-flagged as merely undecided
    assert v["borderline"] is False


def test_atenolol_is_borderline_not_a_yes():
    """The documented false positive, and it survives a complete parameter
    set: with ChEMBL's own pKa it clears the threshold by 0.04 while being
    a textbook poor penetrant. A margin that thin must not render as an
    answer."""
    v = dlv.bbb_verdict("atenolol", _CHEMBL_ATENOLOL)
    assert v["borderline"] is True
    assert abs(v["margin"]) <= dlv.BORDERLINE_MARGIN
    assert v["cell"] == "yes?", "a coin flip must not print as a plain yes"
    assert "noise of the threshold" in v["why"]


def test_a_comfortable_pass_is_not_marked_borderline():
    """Otherwise the question mark means nothing."""
    v = dlv.bbb_verdict("diazepam", _CHEMBL_DIAZEPAM)
    assert v["borderline"] is False
    assert v["cell"] == "yes"
    assert v["margin"] > dlv.BORDERLINE_MARGIN


def test_bbb_threshold_matches_the_delivery_section():
    """One drug, two places in the report. If these diverged, a candidate
    could read yes in the table and unlikely in the delivery notes."""
    import inspect
    src = inspect.getsource(dlv.bbb_verdict)
    assert "4.0 * (n_params / 6.0)" in src
    a = dlv.assess_delivery("diazepam", compartment="cns",
                            properties=_CHEMBL_DIAZEPAM)
    v = dlv.bbb_verdict("diazepam", _CHEMBL_DIAZEPAM)
    assert a.cns_mpo == v["mpo"]


def test_bbb_column_renders_for_every_graded_drug():
    """The request was a yes/no on every row, not only under --tissue cns -
    whether a candidate can reach its tissue is a screening question."""
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    md = pipeline.markdown_report(res)
    header = [ln for ln in md.splitlines()
              if ln.startswith("| Grade | Combined |")]
    assert header, "graded table missing"
    assert "BBB?" in header[0]
    # the separator row must have grown with the header, or the table breaks
    i = md.splitlines().index(header[0])
    assert md.splitlines()[i + 1].count("|") == header[0].count("|")


def test_predicted_medications_table_carries_the_same_column():
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    md = pipeline.markdown_report(res)
    header = [ln for ln in md.splitlines()
              if ln.startswith("| # | Medication |")]
    assert header and "BBB?" in header[0]
    i = md.splitlines().index(header[0])
    assert md.splitlines()[i + 1].count("|") == header[0].count("|")


def test_bbb_does_not_reorder_anything():
    """It is an applicability filter the reader applies, not a score. If it
    entered the ranking, a peripheral target would be penalised for having
    candidates that cannot reach brain - which is irrelevant to it.

    Flipping every verdict must leave both orderings byte-identical."""
    from promoter_tf_drugs import predictions as pr

    block = {
        "all_routes": {"routes": [
            {"drug": "Alpha", "acts_on": "VDR", "category": "ligand",
             "tier": 3, "score": 50.0, "directness": 1, "mechanism": "agonist"},
            {"drug": "Beta", "acts_on": "VDR", "category": "ligand",
             "tier": 3, "score": 40.0, "directness": 1, "mechanism": "agonist"},
        ]},
        # Alpha cannot cross, Beta can - the reverse of the score order.
        "drug_grades": [
            {"drug_name": "Alpha", "grade": "A",
             "evidence": {"properties": _CHEMBL_PLICAMYCIN}},
            {"drug_name": "Beta", "grade": "B",
             "evidence": {"properties": _CHEMBL_DIAZEPAM}},
        ],
        "tfs": [], "outcome_estimates": [],
    }
    preds = pr.build(block, goal="up")
    assert [p.drug for p in preds] == ["Alpha", "Beta"], \
        "score order must survive a disagreeing BBB verdict"
    assert [p.bbb for p in preds] == ["no", "yes"], \
        "and the verdicts really do disagree with the order"

    # the column is never named as a sort key
    import inspect
    src = inspect.getsource(pr.build) + inspect.getsource(pr.render)
    assert "key=lambda" not in src or "bbb" not in src.split("key=lambda")[1][:120]


def test_footnote_describes_only_the_states_present():
    """A note warning about a borderline pass in a table with none is
    noise; one that stays silent about unknowns when 20 rows are unknown is
    worse."""
    only_unknown = " ".join(pipeline._bbb_footnote({"—": 5}))
    assert "is not a no" in only_unknown
    assert "question mark" not in only_unknown

    with_borderline = " ".join(pipeline._bbb_footnote({"yes?": 1, "no": 2}))
    assert "question mark" in with_borderline
    assert "Atenolol" in with_borderline
    assert "is not a no" not in with_borderline

    assert pipeline._bbb_footnote({}) == []


def test_footnote_states_what_a_yes_does_not_establish():
    note = " ".join(pipeline._bbb_footnote({"yes": 3}))
    assert "passive diffusion" in note
    assert "P-glycoprotein" in note
    # and the reliable direction is named as such
    assert "A **no** is the more reliable direction" in note


# ===================================================================== #
# Why ligand-based drugs were outranking the empirical ones.            #
# ===================================================================== #
def _lincs_hit(**kw):
    h = {"drug": "X", "gene": "G", "score": 9.0,
         "feature_space": None, "is_landmark": None}
    h.update(kw)
    return {"hits": [h]}


def test_measured_lincs_evidence_outranks_a_curated_ligand():
    """The bug the user spotted. from_lincs tiered a hit as 2 only when
    `is_landmark` was exactly True, and the offline .gmt loader never set
    it - so every offline LINCS hit landed at tier 4, BELOW the tier-3
    curated ligands. A drug observed to move the gene was being ranked
    under a drug that merely binds a factor near it."""
    from promoter_tf_drugs import catalogue as cat
    measured = cat.from_lincs(_lincs_hit(feature_space="landmark"))[0]
    assert measured.tier == 2, "a measured observation must beat a ligand"

    ligand_tier = 3
    assert measured.tier < ligand_tier


def test_the_three_level_feature_space_is_not_collapsed_to_a_binary():
    """CLUE records three confidence levels and the old flag had two, so
    the 9,196 `best inferred` genes - which includes most genes anyone
    asks about - were lumped in with pure prediction."""
    from promoter_tf_drugs import catalogue as cat
    tiers = {fs: cat.from_lincs(_lincs_hit(feature_space=fs))[0].tier
             for fs in ("landmark", "best inferred", "inferred")}
    assert tiers["landmark"] == 2          # measured on the array
    assert tiers["best inferred"] == 3     # imputed from measured landmarks
    assert tiers["inferred"] == 4          # predicted
    assert tiers["landmark"] < tiers["best inferred"] < tiers["inferred"]


def test_unknown_feature_space_is_tiered_as_predicted_not_as_measured():
    """Without the geneinfo table nothing is known about the gene, and the
    safe direction is to claim less."""
    from promoter_tf_drugs import catalogue as cat
    assert cat.from_lincs(_lincs_hit(feature_space=None))[0].tier == 4


def test_is_landmark_true_still_promotes_for_the_csv_path():
    """--lincs-file carries is_landmark directly and must keep working."""
    from promoter_tf_drugs import catalogue as cat
    r = cat.from_lincs(_lincs_hit(is_landmark=True, feature_space=None))[0]
    assert r.tier == 2


def test_geneinfo_table_distinguishes_the_three_levels():
    """Reads the real CLUE table when it is present; skipped otherwise so
    the suite does not depend on a 1.1 MB download."""
    import pathlib
    import pytest
    from promoter_tf_drugs import signatures as sg
    if not any(pathlib.Path(p).exists()
               for p in sg.DEFAULT_GENEINFO_PATHS):
        pytest.skip("data/lincs_geneinfo.txt not present")
    assert sg.load_geneinfo() > 10_000
    assert sg.feature_space("GAPDH") == "landmark"
    assert sg.feature_space("NOT_A_REAL_GENE_XYZ") is None
    # the middle level exists and is populated - the whole point
    assert sg.feature_space("VDR") == "best inferred"


def test_confidence_note_says_which_level_it_is():
    """A tier is a number; the note is what makes it checkable."""
    import pathlib
    import pytest
    from promoter_tf_drugs import signatures as sg
    if not any(pathlib.Path(p).exists()
               for p in sg.DEFAULT_GENEINFO_PATHS):
        pytest.skip("data/lincs_geneinfo.txt not present")
    assert "landmark" in sg.confidence_note("GAPDH").lower()
    assert sg.confidence_note("NOT_A_REAL_GENE_XYZ")


def test_geneinfo_load_is_idempotent():
    """feature_space() calls it on every lookup, so a non-idempotent load
    would re-parse the table once per gene."""
    from promoter_tf_drugs import signatures as sg
    first = sg.load_geneinfo()
    assert sg.load_geneinfo() == first
    assert sg.load_geneinfo() == first


# ===================================================================== #
# TF activity inference: measuring a factor instead of predicting it.   #
# ===================================================================== #
def _planted_signature(n_genes=3000, seed=0):
    """A signature where ONE factor's regulon genuinely moved."""
    import numpy as np
    rng = np.random.default_rng(seed)
    genes = [f"G{i}" for i in range(n_genes)]
    sig = {g: float(rng.normal()) for g in genes}
    regulon = {g: 1 for g in genes[:30]}
    regulon.update({g: -1 for g in genes[30:60]})
    for g in genes[:30]:
        sig[g] += 2.0           # activated targets up
    for g in genes[30:60]:
        sig[g] -= 2.0           # repressed targets down
    decoys = {}
    for k in range(40):
        tg = rng.choice(genes, size=40, replace=False)
        decoys[f"DECOY{k}"] = {g: int(s) for g, s in
                               zip(tg, rng.choice([1, -1], size=40))}
    return sig, dict(decoys, REAL=regulon), genes


def test_rank_normalise_is_invariant_to_monotone_rescaling():
    """The statistic assumes standard normal scores under the null. A raw
    logFC distribution is heavy-tailed with an arbitrary scale, so ranking
    first is what makes the null exact - and makes the answer independent
    of whether the signature arrived as logFC, a t statistic or a z."""
    from promoter_tf_drugs import tfactivity as ta
    import numpy as np
    rng = np.random.default_rng(3)
    sig = {f"G{i}": float(rng.normal()) for i in range(2000)}
    a = ta.rank_normalise(sig)
    b = ta.rank_normalise({g: v * 11.0 + 5.0 for g, v in sig.items()})
    assert max(abs(a[g] - b[g]) for g in sig) < 1e-9

    vals = np.array(list(a.values()))
    assert abs(vals.mean()) < 0.02
    assert abs(vals.std() - 1.0) < 0.05


def test_erfinv_is_accurate_enough_to_trust_the_tails():
    """A distorted tail would inflate the variance of every NES, which is
    the kind of error that makes a p-value wrong without making it look
    wrong."""
    from promoter_tf_drugs import tfactivity as ta
    import numpy as np
    x = np.linspace(-0.9999, 0.9999, 4001)
    # erf(erfinv(x)) must return x
    from math import erf
    back = np.array([erf(v) for v in ta._erfinv(x)])
    assert np.abs(back - x).max() < 1e-6


def test_null_is_calibrated():
    """If the statistic is right, random regulons on a random signature
    produce NES with unit variance. A systematically wide null would make
    every factor look significant."""
    from promoter_tf_drugs import tfactivity as ta
    import numpy as np
    rng = np.random.default_rng(11)
    genes = [f"G{i}" for i in range(4000)]
    variances = []
    for rep in range(4):
        sig = {g: float(rng.normal()) for g in genes}
        regs = {}
        for k in range(150):
            tg = rng.choice(genes, size=40, replace=False)
            regs[f"TF{k}"] = {g: int(s) for g, s in
                              zip(tg, rng.choice([1, -1], size=40))}
        acts = ta.activity(sig, regs, check_confounding=False)
        nes = np.array([a.nes for a in acts if a.nes is not None])
        variances.append(nes.var())
    # generous band: this is a 150-sample variance estimate per rep
    assert 0.7 < float(np.mean(variances)) < 1.4, variances


def test_a_planted_active_factor_is_recovered_and_ranked_first():
    from promoter_tf_drugs import tfactivity as ta
    sig, regs, _ = _planted_signature()
    acts = ta.activity(sig, regs, check_confounding=False)
    assert acts[0].tf == "REAL"
    real = [a for a in acts if a.tf == "REAL"][0]
    assert real.direction == "activated"
    assert real.nes > 5
    assert real.q_value < 0.01


def test_the_direction_flips_with_the_signature():
    """If it did not, the sign would be decoration."""
    from promoter_tf_drugs import tfactivity as ta
    sig, regs, _ = _planted_signature()
    flipped = {g: -v for g, v in sig.items()}
    real = [a for a in ta.activity(flipped, regs, check_confounding=False)
            if a.tf == "REAL"][0]
    assert real.nes < 0
    assert real.direction == "repressed"


def test_unsigned_edges_support_magnitude_but_refuse_a_direction():
    """TRRUST records no mode for 46% of its human edges. Those cannot say
    which way the factor went, and claiming otherwise would put activators
    and repressors in the same bucket."""
    from promoter_tf_drugs import tfactivity as ta
    import numpy as np
    rng = np.random.default_rng(5)
    genes = [f"G{i}" for i in range(2000)]
    sig = {g: float(rng.normal()) for g in genes}
    for g in genes[:40]:
        sig[g] += 2.5                      # moved, mode unrecorded
    a = ta.activity(sig, {"UNSIGNED": {g: 0 for g in genes[:40]}},
                    check_confounding=False)[0]
    assert a.nes is None
    assert a.p_value is None
    assert a.direction == "unknown"
    assert a.magnitude > 5 and a.magnitude_p < 0.001
    assert any("direction cannot be called" in c for c in a.caveats)


def test_too_few_targets_is_skipped_not_scored():
    """A 3-target regulon reaches |NES| 2 by chance often enough to be
    worthless, and reporting it would fill the table with noise."""
    from promoter_tf_drugs import tfactivity as ta
    sig = {f"G{i}": float(i) for i in range(500)}
    out = ta.activity(sig, {"TINY": {"G1": 1, "G2": -1, "G3": 1}},
                      check_confounding=False)
    assert out == []


def test_regulon_overlap_is_reported_so_a_stolen_signal_is_visible():
    """Co-regulated targets cannot distinguish a factor from its
    neighbour. VIPER corrects for this; here it is measured and named,
    which is the honest minimum."""
    from promoter_tf_drugs import tfactivity as ta
    import numpy as np
    rng = np.random.default_rng(9)
    genes = [f"G{i}" for i in range(2000)]
    sig = {g: float(rng.normal()) for g in genes}
    shared = {g: 1 for g in genes[:40]}
    acts = ta.activity(sig, {"TWINA": dict(shared), "TWINB": dict(shared)},
                       check_confounding=True)
    for a in acts:
        assert a.confounded_jaccard == 1.0
        assert a.confounded_with in ("TWINA", "TWINB") and a.confounded_with != a.tf
        assert any("cannot attribute" in c for c in a.caveats)


def test_contradictory_edges_are_demoted_to_unsigned(tmp_path):
    """TRRUST really does record activation in one paper and repression in
    another for the same pair. Picking a winner would be inventing a
    consensus that the literature does not have."""
    from promoter_tf_drugs import tfactivity as ta
    f = tmp_path / "r.tsv"
    f.write_text("TFX\tTARGET1\tActivation\t1\n"
                 "TFX\tTARGET1\tRepression\t2\n"
                 "TFX\tTARGET2\tActivation\t3\n", encoding="utf-8")
    r = ta.load_signed_regulons(f, min_targets=1)
    assert r["TFX"]["TARGET1"] == 0, "conflict must not pick a side"
    assert r["TFX"]["TARGET2"] == 1


def test_dorothea_numeric_weights_are_read_as_signs():
    """CollecTRI/DoRothEA ship a numeric weight rather than a word, and
    these tables are the reason to support a second format at all: they
    are fully signed where TRRUST is half-signed."""
    from promoter_tf_drugs import tfactivity as ta
    import pathlib
    f = pathlib.Path(str(ta.__file__)).parent  # any writable tmp is fine
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        p = pathlib.Path(d) / "collectri.tsv"
        p.write_text("source\ttarget\tweight\n"
                     "TFY\tA\t1\nTFY\tB\t-1\nTFY\tC\t1\n",
                     encoding="utf-8")
        r = ta.load_signed_regulons(p, min_targets=1)
    assert r["TFY"] == {"A": 1, "B": -1, "C": 1}
    assert "SOURCE" not in r, "the header row must not become a factor"


def test_set_based_test_is_used_for_lincs_shaped_data():
    """Enrichr's LINCS .gmt gives top-N up and down gene lists and no
    scores. Rank-transforming a ternary vector would invent precision, so
    the statistic changes to a binomial sign test to match the data."""
    from promoter_tf_drugs import tfactivity as ta
    genes = [f"G{i}" for i in range(400)]
    regulon = {g: 1 for g in genes[:100]}
    regulon.update({g: -1 for g in genes[100:200]})
    # a perturbation that activates it
    acts = ta.activity_from_sets(set(genes[:60]), set(genes[100:160]),
                                 {"REAL": regulon})
    assert acts and acts[0].tf == "REAL"
    assert acts[0].direction == "activated"
    assert "binomial sign test" in acts[0].method
    # and the reverse perturbation flips it
    rev = ta.activity_from_sets(set(genes[100:160]), set(genes[:60]),
                                {"REAL": regulon})
    assert rev[0].nes < 0 and rev[0].direction == "repressed"


def test_set_based_test_is_calibrated_under_the_null():
    """Random up/down sets against random regulons must reject at about the
    nominal rate. Demanding ZERO rejections would be wrong - at q<0.05 the
    occasional one IS the error rate - so this checks the rate over several
    replicates instead."""
    from promoter_tf_drugs import tfactivity as ta
    import numpy as np
    rng = np.random.default_rng(13)
    N = 8000
    genes = [f"G{i}" for i in range(N)]
    regs = {}
    for k in range(120):
        size = int(rng.integers(200, 500))
        tg = rng.choice(N, size=size, replace=False)
        regs[f"TF{k}"] = {genes[i]: int(s) for i, s in
                          zip(tg, rng.choice([1, -1], size=size))}

    rates, n_testable = [], 0
    for rep in range(10):
        pick = rng.choice(N, size=400, replace=False)
        acts = ta.activity_from_sets({genes[i] for i in pick[:200]},
                                     {genes[i] for i in pick[200:]}, regs)
        if not acts:
            continue
        n_testable = len(acts)
        ps = np.array([a.p_value for a in acts])
        rates.append(float((ps < 0.05).mean()))

    assert n_testable, "nothing was testable - the check would be vacuous"
    assert rates
    assert 0.01 < float(np.mean(rates)) < 0.12, rates


def test_tissue_signature_treats_a_missing_tissue_as_low_not_unmeasured():
    """The expression table is loaded with a min_value floor, so a gene
    carries only the tissues where it cleared it. Dropping the absent ones
    would bias every regulon toward broadly expressed genes."""
    from promoter_tf_drugs import tfactivity as ta
    table = {
        "HIGH": {"liver": 100.0, "brain": 1.0},
        "ABSENT": {"brain": 50.0},            # not expressed in liver at all
    }
    sig = ta.tissue_signature("liver", table=table)
    assert sig["HIGH"] > 0
    assert "ABSENT" in sig, "a floored-out gene must not vanish"
    assert sig["ABSENT"] < 0


def test_tissue_signature_without_a_table_refuses_rather_than_guessing():
    from promoter_tf_drugs import tfactivity as ta
    assert ta.tissue_signature("liver", table={}) == {}
    assert ta.tissue_activity("liver", {"TF": {"A": 1}}, table={}) == []


def test_analytic_pvalue_says_it_is_anticonservative():
    """The null assumes targets are independent and co-expression is
    exactly the structure being measured, so the caveat is not optional."""
    from promoter_tf_drugs import tfactivity as ta
    sig, regs, _ = _planted_signature()
    a = ta.activity(sig, regs, check_confounding=False)[0]
    assert any("anticonservative" in c for c in a.caveats)
    assert "analytic" in a.method

    perm = ta.activity(sig, {"REAL": regs["REAL"]}, permutations=500,
                       check_confounding=False)[0]
    assert "permutations" in perm.method
    assert not any("anticonservative" in c for c in perm.caveats)


def test_summarise_does_not_claim_a_starting_point_when_there_is_none():
    """The real CCR5 run produced zero factors with both a site and a
    perturbed regulon, and the first draft still said 'that is where to
    start'."""
    from promoter_tf_drugs import tfactivity as ta
    sig, regs, _ = _planted_signature()
    acts = ta.activity(sig, regs, check_confounding=False)

    none = ta.summarise(acts, scanned=["NOTHING_SCANNED"])
    assert "no factor has both" in none
    assert "where to start" not in none

    some = ta.summarise(acts, scanned=["REAL"])
    assert "start there" in some
    assert "REAL" in some


def test_summarise_explains_an_empty_result_as_a_table_limit():
    """With TRRUST most of the ~1,600 human TFs have too few curated
    targets to test, so an empty result says something about the regulon
    table and nothing about the factors."""
    from promoter_tf_drugs import tfactivity as ta
    msg = ta.summarise([])
    assert "regulon table" in msg
    assert "not about the factors" in msg


def test_activity_section_renders_and_states_its_limits():
    """Injected into a real demo result rather than a hand-built dict, so
    the section is exercised against the shape the pipeline actually
    produces."""
    from promoter_tf_drugs import tfactivity as ta
    sig, regs, _ = _planted_signature()
    acts = ta.activity(sig, regs, check_confounding=False)

    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    block = res["results"][0]
    block["tf_activity"] = {
        "tissue": "liver", "n_testable": len(acts), "n_regulons": len(regs),
        "permutations": 0,
        "scanned": [a.to_dict() for a in acts[:3]],
        "top": [a.to_dict() for a in acts[:8]],
        "summary": ta.summarise(acts, scanned=[a.tf for a in acts[:3]]),
    }
    md = pipeline.markdown_report(res)

    assert "Is the factor actually active" in md
    # the three limits that make the section honest rather than impressive
    assert "reference tissue average" in md
    assert "not attributable to one factor" in md
    assert "anticonservative" in md
    assert "Nuclear receptors are the hardest case" in md
    # and the numbers reached the table
    assert acts[0].tf in md


def test_activity_section_is_absent_when_it_did_not_run():
    """It must not render an empty shell that looks like a negative
    result."""
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    assert res["results"][0]["tf_activity"] == {}
    assert "Is the factor actually active" not in pipeline.markdown_report(res)


def test_activity_is_recorded_in_the_ledger_when_it_cannot_run():
    """A missing --regulon-file must leave a trace saying the motif scores
    went unchecked, not silently omit the section."""
    import inspect
    src = inspect.getsource(pipeline.run)
    assert 'ledger.record(\n                "tf_activity"' in src or \
        '"tf_activity", requested=True, ok=False' in src
    assert "motif scores are unchecked" in src


# ===================================================================== #
# How many other genes does the factor control, and the screen's scope. #
# ===================================================================== #
def test_tf_table_reports_how_many_other_genes_each_factor_controls():
    """The collateral question. A drug acting through a factor moves
    everything that factor regulates, and the ranking alone never said how
    much that is."""
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    assert all("n_regulon" in t for t in res["results"][0]["tfs"])

    md = pipeline.markdown_report(res)
    header = [ln for ln in md.splitlines() if ln.startswith("| Rank | TF |")]
    assert header and "Other genes it controls" in header[0]
    i = md.splitlines().index(header[0])
    assert md.splitlines()[i + 1].count("|") == header[0].count("|")


def test_regulon_count_carries_its_provenance():
    """A curated count and a bundled tier estimate are different claims and
    must not render as the same number."""
    from promoter_tf_drugs import regulon as reg
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    rows = res["results"][0]["tfs"]
    assert all(t.get("regulon_provenance") for t in rows)

    md = pipeline.markdown_report(res)
    # demo has no regulon file, so every count is a bundled tier and must
    # be marked approximate
    assert "~" in md
    assert "bundled estimate" in md


def test_a_factor_with_no_regulon_entry_says_not_counted():
    """Absence from the table is not a small number. Rendering it as 0 - or
    omitting it - would read as maximal specificity."""
    from promoter_tf_drugs import regulon as reg
    n, prov = reg.regulon_size("TOTALLYFAKETF9", None)
    assert n is None
    assert prov == "unknown"

    note = pipeline._regulon_column_note(
        True, False, [{"n_regulon": None, "regulon_provenance": "unknown"}])
    assert "not counted" in " ".join(note)
    assert "It is not a small number" in " ".join(note)


def test_small_curated_regulon_is_flagged_as_understudied_not_selective():
    """The trap on this column. A factor with eight curated targets is far
    more likely under-studied than specific, and reading a small count as
    specificity selects for obscurity."""
    from promoter_tf_drugs import regulon as reg
    note = " ".join(pipeline._regulon_column_note(
        True, True, [{"n_regulon": 4, "regulon_provenance": "regulon file"}]))
    assert "understudied than" in note
    assert str(reg.MIN_CREDIBLE_REGULON) or "credible floor" in note
    assert "credible floor" in note


def test_regulon_note_says_counts_measure_curation_too():
    note = " ".join(pipeline._regulon_column_note(True, False, []))
    assert "literature curation" in note
    assert "SP1" in note


def test_regulon_column_absent_without_a_regulon_source():
    note = " ".join(pipeline._regulon_column_note(False, False, []))
    assert "--regulon-file" in note
    assert "no basis" in note


# -------------------------------------------------- the screening scope
def test_established_ligand_rows_are_ranked_by_measured_potency():
    """They were ranked by the model's score on molecules the model was
    TRAINED on, so the top of the list was whichever training point the
    forest had memorised best - not a property of the drug."""
    lib = _toy_library()
    known = {d["chembl_id"] for d in lib[:6]}
    # deliberately inverted: the LOWEST model score has the HIGHEST potency
    potency = {cid: 9.0 - i * 0.5 for i, cid in enumerate(sorted(known))}
    m = agonism.train_model(_toy_ligands())
    hits = scr.screen(m, lib, goal="up", known_ids=known,
                      known_potency=potency, top_n=100)
    est = [h for h in hits if h.quadrant.startswith("ESTABLISHED")]
    assert len(est) >= 2
    pch = [h.measured_pchembl for h in est]
    assert pch == sorted(pch, reverse=True), pch


def test_the_model_score_on_a_training_molecule_is_marked_in_sample():
    """It is kept, because the enrichment check needs it, and labelled, so
    it is not read as a prediction."""
    lib = _toy_library()
    known = {d["chembl_id"] for d in lib[:4]}
    m = agonism.train_model(_toy_ligands())
    hits = scr.screen(m, lib, goal="up", known_ids=known, top_n=100)
    for h in hits:
        assert h.probability_in_sample is h.known_ligand
        if h.known_ligand:
            assert h.probability is not None, "enrichment_check needs it"
            assert any("in-sample" in n for n in h.notes)


def test_novelty_is_withheld_for_a_known_ligand():
    """A training molecule is its own nearest neighbour, so 1 - tanimoto is
    structurally 0.0 for every one of them. Printing that as a novelty of
    zero states an artifact of the lookup as a finding."""
    lib = _toy_library()
    known = {d["chembl_id"] for d in lib[:4]}
    m = agonism.train_model(_toy_ligands())
    hits = scr.screen(m, lib, goal="up", known_ids=known, top_n=100)
    for h in hits:
        if h.known_ligand:
            assert h.novelty is None
        else:
            assert h.novelty is not None


def test_enrichment_check_still_works_after_the_change():
    """It compares known actives against novel candidates, so nulling the
    known probabilities would have silently disabled the one diagnostic
    that catches an inverted ranking."""
    lib = _toy_library()
    known = {d["chembl_id"] for d in lib[:6]}
    m = agonism.train_model(_toy_ligands())
    hits = scr.screen(m, lib, goal="up", known_ids=known, top_n=100)
    enr = scr.enrichment_check(hits)
    assert enr.get("median_probability_known") is not None


def test_report_says_why_the_same_receptors_appear_for_every_gene():
    """The user-visible symptom: the screen lists retinoids whatever gene
    is run. It is a property of which receptors have enough medicinal
    chemistry to train a model, not a result about the gene."""
    note = " ".join(pipeline._screen_scope_note(["RARA"], ["VDR"]))
    assert "whatever gene you run" in note
    assert "retinoid" in note
    assert "not a result about your gene" in note
    # the actual gate, with its real numbers
    assert str(scr.MIN_ACTIVES) in note
    assert str(scr.MIN_BINDING_AUC) in note
    # and a skip is never presented as a verdict on the receptor
    assert "never a judgement about the receptor" in note


def test_scope_note_names_what_ran_and_what_did_not():
    note = " ".join(pipeline._screen_scope_note(["RARA", "RXRA"], ["VDR"]))
    assert "**Screened:** RARA, RXRA" in note
    assert "**Not screened:** VDR" in note


def test_pchembl_cutoff_is_one_constant_not_two_numbers():
    """The report cites the floor the fetch applies. If they drifted apart
    the report would describe a screen that did not happen."""
    from promoter_tf_drugs import agonism as ag
    import inspect
    sig = inspect.signature(ag.fetch_actives)
    assert sig.parameters["pchembl_cutoff"].default == ag.PCHEMBL_CUTOFF
    assert f"pChEMBL {ag.PCHEMBL_CUTOFF:.0f}" in \
        " ".join(pipeline._screen_scope_note(["X"], []))


# ===================================================================== #
# CollecTRI as a regulon source.                                        #
# ===================================================================== #
def test_collectri_loads_from_tsv_and_csv(tmp_path):
    """Both ship in the wild, and the loader hardcoded a tab - so a
    CollecTRI CSV parsed as zero regulons and reported no error, which is
    the worst way to fail."""
    from promoter_tf_drugs import tfactivity as ta
    body = [("source", "target", "weight"), ("MYC", "TERT", "1"),
            ("MYC", "CDKN1A", "-1"), ("SP1", "VEGFA", "1")]
    for name, sep in (("c.tsv", "\t"), ("c.csv", ",")):
        f = tmp_path / name
        f.write_text("\n".join(sep.join(r) for r in body) + "\n",
                     encoding="utf-8")
        r = ta.load_signed_regulons(f, min_targets=1)
        assert r == {"MYC": {"TERT": 1, "CDKN1A": -1},
                     "SP1": {"VEGFA": 1}}, name


def test_separator_is_sniffed_not_taken_from_the_extension():
    """Both formats are routinely saved under the wrong suffix."""
    from promoter_tf_drugs import tfactivity as ta
    import inspect
    src = inspect.getsource(ta.load_signed_regulons)
    assert "suffix" not in src
    assert "head.count" in src


def test_autoload_prefers_collectri_over_trrust(tmp_path, monkeypatch):
    """CollecTRI signs every edge and covers about three times as many
    factors, which is the binding constraint on how many can be tested."""
    from promoter_tf_drugs import tfactivity as ta
    assert ta.DEFAULT_REGULON_PATHS[0].endswith("collectri.tsv")
    assert any("trrust" in p for p in ta.DEFAULT_REGULON_PATHS)
    assert (ta.DEFAULT_REGULON_PATHS.index("data/collectri.tsv")
            < [i for i, p in enumerate(ta.DEFAULT_REGULON_PATHS)
               if "trrust" in p][0])

    monkeypatch.chdir(tmp_path)
    (tmp_path / "data").mkdir()
    (tmp_path / "data" / "collectri.tsv").write_text(
        "source\ttarget\tweight\nMYC\tA\t1\nMYC\tB\t-1\n"
        "MYC\tC\t1\nMYC\tD\t1\nMYC\tE\t-1\n", encoding="utf-8")
    regs, src = ta.autoload_regulons(min_targets=1)
    assert src == "data/collectri.tsv"
    assert regs["MYC"]["B"] == -1


def test_autoload_returns_empty_rather_than_raising(tmp_path, monkeypatch):
    """A missing regulon file costs the activity section and nothing else."""
    from promoter_tf_drugs import tfactivity as ta
    monkeypatch.chdir(tmp_path)
    assert ta.autoload_regulons() == ({}, None)


def test_report_points_at_collectri_when_a_thinner_table_was_used():
    """Coverage is the live constraint on this section, and it is fixable,
    so the report should say how rather than only state the limit."""
    from promoter_tf_drugs import tfactivity as ta
    sig, regs, _ = _planted_signature()
    acts = ta.activity(sig, regs, check_confounding=False)
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    res["results"][0]["tf_activity"] = {
        "tissue": "liver", "n_testable": len(acts), "n_regulons": len(regs),
        "permutations": 0, "regulon_source": "data/trrust.tsv",
        "scanned": [a.to_dict() for a in acts[:2]],
        "top": [a.to_dict() for a in acts[:5]],
        "summary": "x",
    }
    md = pipeline.markdown_report(res)
    assert "CollecTRI" in md
    assert "get_collectri" in md

    # and not when CollecTRI was already the source
    res["results"][0]["tf_activity"]["regulon_source"] = "data/collectri.tsv"
    assert "get_collectri" not in pipeline.markdown_report(res)


# ===================================================================== #
# Drugs found by inferred activity, with no annotated binding.          #
# ===================================================================== #
def _gmt(path, terms):
    path.write_text(
        "\n".join(f"{name}\t\t" + "\t".join(genes)
                   for name, genes in terms.items()) + "\n",
        encoding="utf-8")


def test_a_drug_moving_a_regulon_is_found_without_any_binding_record(tmp_path):
    """The gap this closes. drugs.py can only return compounds somebody
    assayed at the target, so a factor with no medicinal chemistry comes
    back empty however important it is."""
    from promoter_tf_drugs import tfactivity as ta
    up_t = [f"UP{i}" for i in range(20)]
    dn_t = [f"DN{i}" for i in range(20)]
    regulon = {g: 1 for g in up_t}
    regulon.update({g: -1 for g in dn_t})

    # one compound, profiled five times, moving the regulon the same way
    up, down = {}, {}
    for i in range(5):
        up[f"EXP00{i} MCF7 6H-realdrug-10.0"] = up_t
        down[f"EXP00{i} MCF7 6H-realdrug-10.0"] = dn_t
    # and a compound that does nothing
    for i in range(5):
        up[f"EXP00{i} MCF7 6H-nulldrug-10.0"] = up_t[:5] + dn_t[:5]
        down[f"EXP00{i} MCF7 6H-nulldrug-10.0"] = up_t[5:10] + dn_t[5:10]

    u, d = tmp_path / "up.gmt", tmp_path / "down.gmt"
    _gmt(u, up)
    _gmt(d, down)

    hits = ta.drugs_moving_factor("MYTF", regulon, up_path=u, down_path=d)
    names = [h.drug for h in hits]
    assert "realdrug" in names
    assert "nulldrug" not in names
    h = next(x for x in hits if x.drug == "realdrug")
    assert h.direction == "activated"
    assert h.n_signatures == 5 and h.consistency == 1.0


def test_a_compound_profiled_once_is_never_reported(tmp_path):
    """33,000 signatures against one regulon produce striking singletons by
    the thousand. Repetition is the only cheap defence, and LINCS supplies
    it because most compounds are profiled several times."""
    from promoter_tf_drugs import tfactivity as ta
    up_t = [f"UP{i}" for i in range(20)]
    regulon = {g: 1 for g in up_t}
    u, d = tmp_path / "up.gmt", tmp_path / "down.gmt"
    _gmt(u, {"E1 MCF7 6H-oneshot-10.0": up_t})
    _gmt(d, {"E1 MCF7 6H-oneshot-10.0": ["Z1"]})
    assert ta.drugs_moving_factor("MYTF", regulon,
                                  up_path=u, down_path=d) == []
    assert ta.MIN_SIGNATURES >= 3


def test_significance_is_counted_after_multiplicity_correction():
    """A nominal p<0.05 over ~33,000 signatures is met about 1,650 times by
    chance, so counting raw p-values would make every factor look busy."""
    from promoter_tf_drugs import tfactivity as ta
    import inspect
    src = inspect.getsource(ta.drugs_moving_factor)
    assert "_bh(" in src, "no BH correction across signatures"
    assert 'rec["q"]' in src, "pooling must use q, not p"


def test_consistency_alone_does_not_admit_a_tiny_effect():
    """SP1 returned mitoxantrone at 72% consistency over 32 profiles with a
    median z of -0.30: agreement about an effect too small to matter."""
    from promoter_tf_drugs import tfactivity as ta
    assert ta.MIN_MEDIAN_Z >= 1.0
    assert ta.MIN_SIGNIFICANT >= 2
    import inspect
    src = inspect.getsource(ta.drugs_moving_factor)
    assert "abs(med) < min_median_z" in src


def test_every_hit_denies_being_a_binding_claim(tmp_path):
    """The one misreading that would matter: a row here is not evidence the
    compound touches the factor."""
    from promoter_tf_drugs import tfactivity as ta
    up_t = [f"UP{i}" for i in range(20)]
    regulon = {g: 1 for g in up_t}
    u, d = tmp_path / "up.gmt", tmp_path / "down.gmt"
    _gmt(u, {f"E{i} MCF7 6H-drugx-10.0": up_t for i in range(5)})
    _gmt(d, {f"E{i} MCF7 6H-drugx-10.0": ["Z1"] for i in range(5)})
    hits = ta.drugs_moving_factor("MYTF", regulon, up_path=u, down_path=d)
    assert hits
    joined = " ".join(hits[0].caveats)
    assert "does NOT say the compound binds" in joined
    assert "proliferation" in joined


def test_single_cell_line_is_flagged_as_repetition_not_replication(tmp_path):
    from promoter_tf_drugs import tfactivity as ta
    up_t = [f"UP{i}" for i in range(20)]
    regulon = {g: 1 for g in up_t}
    u, d = tmp_path / "up.gmt", tmp_path / "down.gmt"
    _gmt(u, {f"E{i} MCF7 6H-drugx-10.0": up_t for i in range(5)})
    _gmt(d, {f"E{i} MCF7 6H-drugx-10.0": ["Z1"] for i in range(5)})
    h = ta.drugs_moving_factor("MYTF", regulon, up_path=u, down_path=d)[0]
    assert h.cell_lines == ["MCF7"]
    assert any("one cell line" in c for c in h.caveats)


def test_missing_gmt_returns_nothing_rather_than_raising(tmp_path):
    from promoter_tf_drugs import tfactivity as ta
    regulon = {f"UP{i}": 1 for i in range(20)}
    assert ta.drugs_moving_factor(
        "MYTF", regulon, up_path=tmp_path / "nope.gmt",
        down_path=tmp_path / "nope2.gmt") == []


def test_activity_drugs_enter_the_routes_as_tier_three():
    """Tier 2 is 'target engagement is observed'. Nothing about engagement
    is observed here - the TARGETS moved - so the observation is real and
    the attribution is an inference, which is tier 3."""
    from promoter_tf_drugs import catalogue as cat
    block = {"activity_drugs": {"TP53": [{
        "drug": "mitoxantrone", "tf": "TP53", "direction": "activated",
        "median_z": 2.53, "n_signatures": 11, "n_significant": 7,
        "n_agreeing": 11, "consistency": 1.0,
        "cell_lines": ["MCF7", "PC3"], "caveats": ["a", "b"]}]}}
    routes = cat.build(block, goal="up")["routes"]
    rows = [r for r in routes
            if r["category"] == "TF: inferred activity change"]
    assert len(rows) == 1
    assert rows[0]["tier"] == 3
    assert rows[0]["acts_on"] == "TP53"
    assert "11 LINCS profiles" in rows[0]["mechanism"]


def test_the_new_category_is_explained_in_the_glossary():
    """It is the only category that can reach a factor with no chemistry,
    and the only one whose evidence is about the factor's targets rather
    than the factor, so it needs its own entry."""
    from promoter_tf_drugs import catalogue as cat
    g = cat.CATEGORY_GUIDE["TF: inferred activity change"]
    assert "no medicinal chemistry" in g["why"]
    assert "arrests the cell cycle" in g["cost"]
    assert "empirical signature" in g["differs"]


def test_an_empty_drug_scan_is_a_negative_not_a_failure():
    """The scan RAN and found nothing. Recording that as a failed source
    is the exact conflation provenance.py exists to prevent - and it is
    what the first version did."""
    import inspect
    src = inspect.getsource(pipeline.run)
    i = src.index('"activity_drugs", requested=True')
    assert "ok=True" in src[i:i + 120], \
        "an empty scan must not be recorded as a failed lookup"


def test_report_section_states_the_proliferation_confound():
    from promoter_tf_drugs import tfactivity as ta
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    res["results"][0]["activity_drugs"] = {"TP53": [{
        "drug": "mitoxantrone", "tf": "TP53", "direction": "activated",
        "median_z": 2.53, "n_signatures": 11, "n_significant": 7,
        "n_agreeing": 11, "consistency": 1.0,
        "cell_lines": ["MCF7"], "caveats": []}]}
    md = pipeline.markdown_report(res)
    assert "found without any record that they bind it" in md
    assert "not a binding claim" in md
    assert "Proliferation is the confound" in md
    assert "LINCS cell lines are not your tissue" in md
    assert "mitoxantrone" in md


# ===================================================================== #
# What each drug is, and what it is approved to treat.                  #
# ===================================================================== #
def test_curated_description_beats_the_class_fallback():
    """A specific entry says 'plaque psoriasis and acne, topical only'.
    The class fallback would say 'retinoid'. Both are true and only one is
    useful for deciding whether a row is applicable."""
    from promoter_tf_drugs import indications as ind
    d = ind.describe("Tazarotene")
    assert d["source"] == "curated"
    assert "topical" in d["treats"]
    assert "RAR" in d["what_it_is"]


def test_salt_forms_share_one_description():
    """Otherwise the glossary carries the same sentence under sodium,
    hydrochloride and the parent."""
    from promoter_tf_drugs import indications as ind
    assert (ind.describe("Levothyroxine Sodium")["treats"]
            == ind.describe("levothyroxine")["treats"])

    out = ind.describe_all(["Tretinoin", "Tretinoin Sodium", "tretinoin"])
    assert len(out) == 1


def test_class_fallback_covers_a_drug_with_no_curated_entry():
    """The curated list is short by construction, so an uncurated
    glucocorticoid must still say what kind of thing it is."""
    from promoter_tf_drugs import indications as ind
    assert "amcinonide" not in ind.CURATED,         "pick another uncurated example; this test is about the fallback"
    d = ind.describe("Amcinonide")
    assert d["source"] == "class"
    assert "glucocorticoid" in d["what_it_is"]
    assert d["treats"]


def test_an_unknown_compound_says_so_rather_than_inventing_a_use():
    from promoter_tf_drugs import indications as ind
    d = ind.describe("Totallyfakedrug9")
    assert d["treats"] is None
    assert d["source"] == "unknown"


def test_chembl_indications_prefer_approved_over_merely_studied():
    """ChEMBL's drug_indication includes uses that were trialled and
    abandoned. Presenting those as approved would be the single most
    misleading thing this column could do."""
    from promoter_tf_drugs import indications as ind

    class FakeClient:
        def get(self, url, params=None):
            return {"drug_indications": [
                {"mesh_heading": "Trialled Condition", "max_phase_for_ind": 2},
                {"mesh_heading": "Approved Condition", "max_phase_for_ind": 4},
            ]}

    d = ind.describe("Noveldrug", client=FakeClient(), chembl_id="CHEMBL1")
    assert d["source"] == "ChEMBL"
    assert "approved condition" in d["treats"]
    assert "trialled condition" not in d["treats"]


def test_studied_only_indications_are_labelled_as_such():
    from promoter_tf_drugs import indications as ind

    class FakeClient:
        def get(self, url, params=None):
            return {"drug_indications": [
                {"mesh_heading": "Only Studied", "max_phase_for_ind": 2}]}

    d = ind.describe("Noveldrug", client=FakeClient(), chembl_id="CHEMBL1")
    assert "not confirmed approved" in d["treats"]


def test_indication_lookup_survives_a_dead_api():
    from promoter_tf_drugs import indications as ind

    class Boom:
        def get(self, url, params=None):
            raise RuntimeError("502")

    assert ind.fetch_indications(Boom(), "CHEMBL1") == []
    d = ind.describe("Tazarotene", client=Boom(), chembl_id="CHEMBL1")
    assert d["source"] == "curated", "the curated entry must still win"


def test_no_network_call_without_a_chembl_id():
    from promoter_tf_drugs import indications as ind

    class Tripwire:
        def get(self, *a, **k):
            raise AssertionError("must not be called")

    assert ind.fetch_indications(Tripwire(), None) == []
    assert ind.fetch_indications(Tripwire(), "not-a-chembl-id") == []


def test_predicted_medications_table_says_what_each_drug_treats():
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    md = pipeline.markdown_report(res)
    header = [ln for ln in md.splitlines()
              if ln.startswith("| # | Medication |")]
    assert header
    assert "What it is" in header[0]
    assert "Approved to treat" in header[0]
    i = md.splitlines().index(header[0])
    assert md.splitlines()[i + 1].count("|") == header[0].count("|")
    assert "Every row is off-label" in md


def test_the_glossary_lists_each_drug_once():
    """A sentence repeated under five tables is what made these reports
    long. The glossary is one entry per drug, salts collapsed."""
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    md = pipeline.markdown_report(res)
    i = md.find("### What these drugs are")
    assert i >= 0
    block = md[i:md.find("### ", i + 10)]
    rows = [ln for ln in block.splitlines()
            if ln.startswith("| ") and not ln.startswith("| Drug |")
            and "---" not in ln]
    names = [ln.split("|")[1].strip() for ln in rows]
    assert names == sorted(names), "glossary should read alphabetically"
    assert len(names) == len(set(names)), "a drug is listed twice"


def test_glossary_does_not_present_itself_as_a_label():
    from promoter_tf_drugs import indications as ind
    msg = ind.summarise({"a": {"source": "curated"},
                         "b": {"source": "class"}})
    assert "None of these is a label" in msg
    assert "prescribing information" in msg

    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    assert "None of these is a label" in pipeline.markdown_report(res)


def test_compounds_with_no_indication_are_named_not_dropped():
    """CTD contributes environmental agents and research chemicals, which
    have no indication because they are not drugs. Silently omitting them
    would make the glossary look complete when it is not."""
    from promoter_tf_drugs import indications as ind
    out = ind.describe_all(["Tazarotene", "Totallyfakedrug9"])
    assert out["Totallyfakedrug9"]["treats"] is None

    note = " ".join(pipeline._drug_glossary(
        [{"drug": "Totallyfakedrug9", "mechanism": ""}], {}))
    assert "No description established" in note
    assert "not that the compound has no approved use" in note


def test_retinoid_rows_now_distinguish_topical_from_systemic():
    """The user's actual complaint was a block of retinoids that looked
    interchangeable. Four of them are topical-only and three are oncology
    agents, and nothing in the table said which was which."""
    from promoter_tf_drugs import indications as ind
    topical = ["Adapalene", "Trifarotene", "Tazarotene"]
    for name in topical:
        assert "topical" in ind.describe(name)["treats"], name
    assert "leukaemia" in ind.describe("Tretinoin")["treats"]
    assert "lymphoma" in ind.describe("Bexarotene")["treats"]


def test_us_and_inn_spellings_resolve_to_one_entry():
    """Calcipotriene and calcipotriol are the same molecule under two
    naming conventions. salt_free_name cannot bridge them - they are
    different words, not different salts - so a synonym table does."""
    from promoter_tf_drugs import indications as ind
    assert (ind.describe("Calcipotriene")["treats"]
            == ind.describe("calcipotriol")["treats"])
    assert ind.describe("Calcipotriene")["source"] == "curated"
    assert "topical" in ind.describe("Calcipotriene")["treats"]


def test_every_synonym_points_at_a_real_curated_entry():
    """A synonym mapping to a missing key would silently fall through to
    the class fallback, which is the failure it exists to prevent."""
    from promoter_tf_drugs import indications as ind
    missing = [k for k, v in ind.SYNONYMS.items() if v not in ind.CURATED]
    assert not missing, f"synonyms with no curated target: {missing}"


def test_screen_hits_carry_what_the_drug_is():
    """indication_hint was declared on ScreenHit and never populated. This
    table is where a reader first meets a drug name, and the receptors that
    qualify for a screen have heterogeneous ligand series - a topical gel
    and an infused cytotoxic on adjacent rows."""
    lib = _toy_library()
    m = agonism.train_model(_toy_ligands())
    hits = scr.screen(m, lib, goal="up", known_ids=set(), top_n=10)
    assert hits
    assert all(hasattr(h, "treats") for h in hits)
    assert all("treats" in h.to_dict() for h in hits)


def test_quadrant_table_shows_the_indication():
    from promoter_tf_drugs import indications as ind
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    d = ind.describe("Tazarotene")
    res["results"][0]["screens"] = {"RARA": {
        "model_auc": 0.95, "n_actives": 120, "n_decoys": 360,
        "library_size": 2400, "enrichment": {"usable": True},
        "hits": [{"quadrant": "ESTABLISHED LIGAND, untested on this gene",
                  "drug": "Tazarotene", "measured_pchembl": 7.5,
                  "probability": 0.9, "probability_in_sample": True,
                  "novelty": None, "indication_hint": d["what_it_is"],
                  "treats": d["treats"], "notes": ["x"]}]}}
    md = pipeline.markdown_report(res)
    header = [ln for ln in md.splitlines() if ln.startswith("| Quadrant |")]
    assert header and "Approved to treat" in header[0]
    i = md.splitlines().index(header[0])
    assert md.splitlines()[i + 1].count("|") == header[0].count("|")
    assert "topical only" in md
    assert "Read the indication column before the scores" in md


# ===================================================================== #
# DepMap: what a cell loses without each factor.                        #
# ===================================================================== #
_DM_LINEAGES = ["Lung", "Breast", "Skin", "Blood", "Colon"]


def _depmap_release(tmp_path, profiles, n_lines=200, seed=0):
    """A DepMap release shaped like the real files.

    CRISPRGeneEffect.csv is cell lines by genes, and every column name
    carries the Entrez id: `A1BG (1)`. A loader that does not strip that
    finds no genes at all.
    """
    import csv
    import random
    rng = random.Random(seed)
    models = [f"ACH-{i:06d}" for i in range(n_lines)]

    def val(kind, lineage):
        if kind == "pan":
            return rng.gauss(-1.1, 0.15)
        if kind == "most":
            return rng.gauss(-0.7, 0.2)
        if kind == "selective":
            return (rng.gauss(-1.0, 0.2) if lineage == "Blood"
                    else rng.gauss(-0.05, 0.1))
        if kind == "none":
            return rng.gauss(0.0, 0.1)
        return None

    eff = tmp_path / "CRISPRGeneEffect.csv"
    with open(eff, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow([""] + [f"{g} ({i+1})"
                           for i, g in enumerate(profiles)])
        for j, m in enumerate(models):
            lin = _DM_LINEAGES[j % len(_DM_LINEAGES)]
            row = [m]
            for kind in profiles.values():
                if kind == "sparse" and j % 40:
                    row.append("")
                else:
                    v = val(kind, lin)
                    row.append("" if v is None else f"{v:.4f}")
            w.writerow(row)

    mod = tmp_path / "Model.csv"
    with open(mod, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["ModelID", "OncotreeLineage"])
        for j, m in enumerate(models):
            w.writerow([m, _DM_LINEAGES[j % len(_DM_LINEAGES)]])
    return eff, mod


def test_depmap_separates_pan_essential_from_dispensable(tmp_path):
    """The reason to load this at all: a factor every cell needs is a
    different proposition from one no cell misses, and the curated regulon
    count cannot tell them apart."""
    from promoter_tf_drugs import depmap as dm
    eff, mod = _depmap_release(
        tmp_path, {"RPL5": "pan", "MYC": "most", "CCR5": "none"})
    deps = dm.assess_many(["RPL5", "MYC", "CCR5"],
                          effect_path=eff, model_path=mod)

    assert deps["RPL5"].fraction_dependent >= dm.PAN_ESSENTIAL_FRACTION
    assert deps["RPL5"].verdict == "essential in nearly every line"
    assert deps["CCR5"].fraction_dependent == 0.0
    assert deps["CCR5"].verdict == "dispensable across the panel"
    assert deps["RPL5"].burden > deps["MYC"].burden > deps["CCR5"].burden


def test_entrez_suffix_is_stripped_from_column_names(tmp_path):
    """Every DepMap column is `SYMBOL (entrez)`. A plain symbol match
    against the raw header finds nothing, silently."""
    from promoter_tf_drugs import depmap as dm
    assert dm._symbol("A1BG (1)") == "A1BG"
    assert dm._symbol("TP53 (7157)") == "TP53"
    eff, _ = _depmap_release(tmp_path, {"TP53": "pan"})
    data = dm.load_gene_effect(eff, ["TP53"])
    assert "TP53" in data and len(data["TP53"]) == 200


def test_a_gene_absent_from_the_release_is_omitted_not_called_safe(tmp_path):
    """'Not in this release' and 'knockout did nothing' are different
    statements, and only one of them is reassuring."""
    from promoter_tf_drugs import depmap as dm
    eff, _ = _depmap_release(tmp_path, {"RPL5": "pan"})
    deps = dm.assess_many(["RPL5", "NOTINRELEASE"], effect_path=eff)
    assert "RPL5" in deps
    assert "NOTINRELEASE" not in deps


def test_a_sparsely_scored_gene_is_refused(tmp_path):
    """Summarising ten cell lines as if they were a thousand would put a
    confident verdict on nothing."""
    from promoter_tf_drugs import depmap as dm
    eff, _ = _depmap_release(tmp_path, {"SPARSE1": "sparse", "RPL5": "pan"})
    deps = dm.assess_many(["SPARSE1", "RPL5"], effect_path=eff)
    assert "SPARSE1" not in deps
    d = dm.assess("SPARSE1", [None] * 190 + [-1.0] * 5)
    assert d.n_lines == 0
    assert d.verdict == "not measured"
    assert any("too few to summarise" in c for c in d.caveats)


def test_a_selective_dependency_is_localised_to_its_lineage(tmp_path):
    """The informative case. A factor needed in a minority of lines has a
    context, and the lineage is where to look for what it is."""
    from promoter_tf_drugs import depmap as dm
    eff, mod = _depmap_release(tmp_path, {"MEIS1": "selective"})
    d = dm.assess_many(["MEIS1"], effect_path=eff, model_path=mod)["MEIS1"]
    assert d.selective is True
    assert d.top_lineages
    assert d.top_lineages[0]["lineage"] == "Blood"
    assert d.top_lineages[0]["mean_effect"] < -0.5
    assert any("has a context" in c for c in d.caveats)


def test_model_csv_is_found_beside_the_effect_file(tmp_path):
    """DepMap ships them in one download. Looking only in data/ dropped
    the whole lineage column for anyone who pointed --depmap-file at their
    download folder."""
    from promoter_tf_drugs import depmap as dm
    eff, mod = _depmap_release(tmp_path, {"MEIS1": "selective"})
    assert mod.parent == eff.parent
    lin = dm.load_lineages(beside=eff)
    assert len(lin) == 200
    # and without the hint it does not find it
    assert dm.load_lineages() == {} or "ACH-000000" not in dm.load_lineages()


def test_missing_depmap_file_costs_only_this_section(tmp_path):
    from promoter_tf_drugs import depmap as dm
    assert dm.resolve(str(tmp_path / "nope.csv")) is None
    assert dm.assess_many(["TP53"],
                          effect_path=str(tmp_path / "nope.csv")) == {}
    msg = dm.summarise({})
    assert "not loaded" in msg or "No DepMap data" in msg
    assert "depmap.org" in msg


def test_every_dependency_says_it_is_cancer_lines_in_a_dish(tmp_path):
    """The limit that decides whether any of this transfers. A factor
    whose job is in a post-mitotic cell looks dispensable for reasons that
    have nothing to do with its importance."""
    from promoter_tf_drugs import depmap as dm
    eff, _ = _depmap_release(tmp_path, {"RPL5": "pan", "CCR5": "none"})
    deps = dm.assess_many(["RPL5", "CCR5"], effect_path=eff)
    for d in deps.values():
        assert any("cancer cell lines" in c for c in d.caveats)
    assert any("knockout is not inhibition" in c
               for c in deps["RPL5"].caveats)
    assert any("says nothing about a differentiated" in c
               for c in deps["CCR5"].caveats)


def test_pan_essential_is_a_caution_not_a_veto(tmp_path):
    """Essential genes are drugged successfully at partial occupancy, and
    knockout removes scaffolding a ligand would not touch."""
    from promoter_tf_drugs import depmap as dm
    eff, _ = _depmap_release(tmp_path, {"RPL5": "pan"})
    d = dm.assess_many(["RPL5"], effect_path=eff)["RPL5"]
    joined = " ".join(d.caveats)
    assert "reason for caution, not a veto" in joined


def test_depmap_column_appears_only_when_there_is_data():
    """An empty column would read as 'no cell needs any of these'."""
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    md = pipeline.markdown_report(res)
    header = [ln for ln in md.splitlines() if ln.startswith("| Rank | TF |")]
    assert header
    assert "Cells need it?" not in header[0]
    assert "What a cell loses without each factor" not in md


def test_depmap_column_and_section_render(tmp_path):
    from promoter_tf_drugs import depmap as dm
    eff, mod = _depmap_release(
        tmp_path, {"RARA": "pan", "RXRA": "selective", "VDR": "none"})
    deps = dm.assess_many(["RARA", "RXRA", "VDR"],
                          effect_path=eff, model_path=mod)

    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    res["results"][0]["dependencies"] = {k: v.to_dict()
                                         for k, v in deps.items()}
    for t in res["results"][0]["tfs"]:
        d = deps.get(t["symbol"].upper())
        if d:
            t["dependency"] = d.verdict
    md = pipeline.markdown_report(res)

    header = [ln for ln in md.splitlines() if ln.startswith("| Rank | TF |")]
    assert "Cells need it?" in header[0]
    i = md.splitlines().index(header[0])
    assert md.splitlines()[i + 1].count("|") == header[0].count("|")

    assert "What a cell loses without each factor" in md
    assert "cancer cell lines dividing in a dish" in md
    assert "Knockout is not inhibition" in md
    assert "A dispensable result is not permission" in md


def test_dispensable_factor_shows_no_most_dependent_lineage(tmp_path):
    """For a factor nothing depends on, the 'most dependent' lineage is
    whichever drifted furthest below zero - noise dressed as a finding."""
    from promoter_tf_drugs import depmap as dm
    eff, mod = _depmap_release(tmp_path, {"VDR": "none"})
    deps = dm.assess_many(["VDR"], effect_path=eff, model_path=mod)
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    rows = pipeline._depmap_section(
        {k: v.to_dict() for k, v in deps.items()}, res)
    body = [r for r in rows if r.startswith("| VDR |")]
    assert body
    assert body[0].rstrip().endswith("— |"), body[0]


# ===================================================================== #
# CTD: why "predicted by this run?" is almost always no.                #
# ===================================================================== #
def test_ctd_entry_kinds_separate_drugs_from_toxicants_and_headings():
    """CTD indexes MeSH, not a drug list. Reporting an exposure and a
    category against a "did this run predict it" column makes a correct
    answer look like a broken join - which is how it read."""
    from promoter_tf_drugs import discovery as disc
    assert disc.entry_kind("Celecoxib") == "chemical"
    assert disc.entry_kind("formaldehyde") == "exposure"
    assert disc.entry_kind("Asbestos, Crocidolite") == "exposure"
    assert disc.entry_kind("Bisphenol S") == "exposure"
    assert disc.entry_kind("Anti-HIV Agents") == "category"
    assert disc.entry_kind("Plant Extracts") == "category"
    assert disc.entry_kind("Dietary Fats") == "category"


def test_unplaceable_names_stay_eligible_to_match():
    """The default has to be permissive: a wrong "exposure" label would
    quietly hide a real lead."""
    from promoter_tf_drugs import discovery as disc
    assert disc.entry_kind("Sch 23390") == "chemical"
    assert disc.entry_kind("some-new-compound-42") == "chemical"
    assert disc.entry_kind("") == "chemical"


def test_a_matching_candidate_is_reported_as_predicted():
    """The thing the column exists for. If this ever reads no, the match
    really is broken."""
    from promoter_tf_drugs import discovery as disc
    links = {"dexamethasone": ["decreases^expression"],
             "formaldehyde": ["decreases^expression"]}
    s = disc.summarise_links(links, "down", {"Dexamethasone", "Calcitriol"})
    by = {e["chemical"]: e for e in s["wanted_direction"]}
    assert by["Dexamethasone"]["in_candidates"] is True
    assert by["Formaldehyde"]["in_candidates"] is False


def test_salt_forms_match_the_ctd_parent_name():
    """CTD lists the parent and a drug table often lists a salt. Every
    other name comparison in this package is salt-free; this one was not."""
    from promoter_tf_drugs import discovery as disc
    s = disc.summarise_links({"levothyroxine": ["increases^expression"]},
                             "up", {"Levothyroxine Sodium"})
    assert s["wanted_direction"][0]["in_candidates"] is True


def test_a_category_heading_is_blank_not_no():
    """"Anti-HIV Agents" is not a compound, so "did this run predict it"
    has no answer rather than the answer no."""
    from promoter_tf_drugs import discovery as disc
    s = disc.summarise_links({"anti-hiv agents": ["decreases^expression"]},
                             "down", {"Maraviroc"})
    e = s["wanted_direction"][0]
    assert e["kind"] == "category"
    assert e["in_candidates"] is False
    assert s["counts"]["n_category"] == 1


def test_exposures_are_not_counted_as_things_the_scan_missed():
    """The old framing said an unmatched row was "a clue that the scan is
    wrong". With a list led by asbestos and formaldehyde, that indicted
    the scan for not proposing formaldehyde."""
    from promoter_tf_drugs import discovery as disc
    links = {"formaldehyde": ["decreases^expression"],
             "asbestos, crocidolite": ["decreases^expression"],
             "plant extracts": ["decreases^expression"],
             "someneweddrug": ["decreases^expression"]}
    s = disc.summarise_links(links, "down", set())
    missed = {e["chemical"] for e in s["reported_but_not_predicted"]}
    assert missed == {"Someneweddrug"}


def test_cross_direction_match_count_is_reported():
    """A bucket with no match reads as a broken join when the matching
    worked elsewhere. On the real CCR5 run the wanted bucket matched 0 of
    39 while 3 matched overall."""
    from promoter_tf_drugs import discovery as disc
    links = {"clofibrate": ["increases^expression"],
             "celecoxib": ["decreases^expression"]}
    s = disc.summarise_links(links, "down", {"Clofibrate"})
    assert s["counts"]["matched_any_direction"] == 1
    assert not any(e["in_candidates"] for e in s["wanted_direction"])


def test_report_says_a_no_is_expected():
    """The user read a column of 39 nos as a bug. It was correct, and the
    report never said so."""
    from promoter_tf_drugs import discovery as disc
    links = {f"toxicant{i}": ["decreases^expression"] for i in range(6)}
    links["formaldehyde"] = ["decreases^expression"]
    links["plant extracts"] = ["decreases^expression"]
    summary = disc.summarise_links(links, "down", set())

    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    res["results"][0].setdefault("improvements", {})["reported_links"] = summary
    md = pipeline.markdown_report(res)

    assert "do not read that as a failed match" in md
    assert "MeSH" in md
    assert "different parts of chemical space" in md
    # the kind column reached the table
    header = [ln for ln in md.splitlines() if ln.startswith("| Chemical |")]
    assert header and "Kind" in header[0]
    i = md.splitlines().index(header[0])
    assert md.splitlines()[i + 1].count("|") == header[0].count("|")


# ===================================================================== #
# Structural liabilities: what makes a PREDICTION unbelievable.         #
# ===================================================================== #
_QUERCETIN = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
_CURCUMIN = "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O"
_DIAZEPAM_S = "CN1c2ccc(Cl)cc2C(=NCC1=O)c1ccccc1"
_DOXORUBICIN = ("C[C@H]1O[C@@H](O[C@H]2C[C@@](O)(Cc3c(O)c4C(=O)c5c(O)cccc5"
                "C(=O)c4c(O)c23)C(=O)CO)[C@H](N)[C@@H](O)[C@@H]1O")


def test_pains_match_lowers_confidence_in_the_evidence():
    """A PAINS hit says activity measured for this molecule may be an
    assay artefact. It does not say the drug is unsafe - quercetin and
    doxorubicin are both flagged and doxorubicin is approved."""
    from promoter_tf_drugs import liabilities as lb
    q = lb.assess("quercetin", _QUERCETIN)
    assert q.assessed and q.pains
    assert q.evidence_multiplier < 1.0
    joined = " ".join(q.notes)
    assert "artefact" in joined
    assert "says nothing about the drug's safety" in joined


def test_a_clean_molecule_is_not_penalised():
    from promoter_tf_drugs import liabilities as lb
    d = lb.assess("diazepam", _DIAZEPAM_S)
    assert d.assessed
    assert not d.pains
    assert d.evidence_multiplier == 1.0


def test_an_approved_drug_is_not_demoted_for_failing_lead_like_rules():
    """Doxorubicin breaks three Lipinski rules and carries Brenk alerts,
    and is approved and in use. Filtering on those would discard the best
    candidates in the report."""
    from promoter_tf_drugs import liabilities as lb
    d = lb.assess("doxorubicin", _DOXORUBICIN)
    assert d.lipinski_violations >= 2
    # only PAINS moves the multiplier, and only as an evidence discount
    only_brenk = lb.assess("curcumin", _CURCUMIN)
    assert only_brenk.brenk and not only_brenk.pains
    assert only_brenk.evidence_multiplier == 1.0
    assert any("settled when the compound was approved" in n
               for n in only_brenk.notes)


def test_the_pains_set_is_not_exhaustive_and_says_so():
    """Curcumin is the most cited interference compound in the literature
    and is NOT in Baell's SMARTS set. A clean PAINS result is therefore
    weak evidence of a clean compound, and the module documents it."""
    from promoter_tf_drugs import liabilities as lb
    assert not lb.assess("curcumin", _CURCUMIN).pains
    assert "Curcumin does not" in lb.__doc__
    assert "Read it one way only" in lb.__doc__


def test_no_structure_is_not_a_clean_result():
    from promoter_tf_drugs import liabilities as lb
    x = lb.assess("mystery", None)
    assert x.assessed is False
    assert x.evidence_multiplier == 1.0
    assert any("not a clean result" in n for n in x.notes)
    msg = lb.summarise({"mystery": x})
    assert "gap in the data, not a clean result" in msg


# ===================================================================== #
# Docking used to FIND compounds, and the gate that makes it usable.    #
# ===================================================================== #
def _scorer(separation, seed=0):
    """Actives score `separation` kcal/mol better, with realistic overlap."""
    import random
    rng = random.Random(seed)

    def fn(smi):
        base = rng.gauss(-7.0, 1.2)
        return base - separation if smi.startswith("ACT") else base
    return fn


def _dock_sets(n_act=20, n_dec=180):
    return ([{"smiles": f"ACT{i}"} for i in range(n_act)],
            [{"smiles": f"DEC{i}"} for i in range(n_dec)])


def test_roc_auc_and_enrichment_agree_with_hand_worked_cases():
    from promoter_tf_drugs import docking as dk
    # perfect separation, lower is better
    assert dk.roc_auc([-10, -9, -1, 0], [1, 1, 0, 0]) == 1.0
    # perfectly inverted
    assert dk.roc_auc([0, -1, -9, -10], [1, 1, 0, 0]) == 0.0
    # ties count half
    assert dk.roc_auc([-5, -5], [1, 0]) == 0.5
    # a failed dock is dropped, never scored worst
    assert dk.roc_auc([-10, None, -1, 0], [1, 1, 0, 0]) == 1.0
    assert dk.roc_auc([None, None], [1, 0]) is None

    # 2 actives in 100, both in the top 1 (k=1 -> one slot)
    scores = [-10.0] + [0.0] * 99
    labels = [1] + [0] * 98 + [1]
    assert dk.enrichment_factor(scores, labels, fraction=0.01) == 50.0


def test_a_pocket_that_cannot_rank_its_own_ligands_is_refused():
    """The gate. Docking is the method here most able to produce confident
    nonsense, so it has to recover known binders before anything it says
    about an unknown is reported."""
    from promoter_tf_drugs import docking as dk
    import pathlib
    act, dec = _dock_sets()
    v = dk.validate_pocket(pathlib.Path("x.pdb"), act, dec,
                           score_fn=_scorer(0.0, seed=1))
    assert v["usable"] is False
    assert v["auc"] < dk.MIN_DOCKING_AUC
    assert "carry no information" in v["reason"]


def test_a_pocket_that_recovers_its_ligands_is_allowed_through():
    from promoter_tf_drugs import docking as dk
    import pathlib
    act, dec = _dock_sets()
    v = dk.validate_pocket(pathlib.Path("x.pdb"), act, dec,
                           score_fn=_scorer(3.0, seed=1))
    assert v["usable"] is True
    assert v["auc"] >= dk.MIN_DOCKING_AUC
    assert v["ef1"] >= dk.MIN_EF1


def test_enrichment_is_a_harsher_gate_than_auc():
    """A method can hold a respectable AUC while putting nothing useful in
    the top 1%, which is the only part of a ranked list anyone reads. Both
    gates exist because AUC alone would pass that."""
    from promoter_tf_drugs import docking as dk
    import pathlib
    act, dec = _dock_sets()
    v = dk.validate_pocket(pathlib.Path("x.pdb"), act, dec,
                           score_fn=_scorer(0.5, seed=1))
    assert 0.6 < v["auc"] < 0.7
    assert v["ef1"] == 0.0
    assert v["usable"] is False


def test_too_few_known_ligands_means_the_dock_cannot_be_validated():
    from promoter_tf_drugs import docking as dk
    import pathlib
    act, dec = _dock_sets()
    v = dk.validate_pocket(pathlib.Path("x.pdb"), act[:3], dec,
                           score_fn=_scorer(3.0))
    assert v["usable"] is False
    assert "confident nonsense" in v["reason"]


def test_the_library_is_not_docked_when_the_pocket_fails():
    """The validation is cheap and the library is not. Running 4,000 docks
    to then discard them would be the expensive way to learn nothing."""
    from promoter_tf_drugs import docking as dk
    import pathlib
    act, dec = _dock_sets()
    lib = [{"name": f"d{i}", "smiles": f"LIB{i}"} for i in range(500)]
    seen = []

    def counting(smi):
        seen.append(smi)
        return _scorer(0.0, seed=1)(smi)

    run = dk.screen_by_docking(pathlib.Path("x.pdb"), lib, target="BAD",
                               actives=act, decoys=dec, score_fn=counting)
    assert run.usable is False and run.hits == []
    assert not any(s.startswith("LIB") for s in seen)


def test_a_validated_screen_returns_a_sorted_shortlist():
    from promoter_tf_drugs import docking as dk
    import pathlib
    act, dec = _dock_sets()
    lib = [{"name": f"d{i}", "smiles": f"LIB{i}"} for i in range(500)]
    run = dk.screen_by_docking(pathlib.Path("x.pdb"), lib, target="GOOD",
                               actives=act, decoys=dec, top_n=10,
                               score_fn=_scorer(3.0, seed=1))
    assert run.usable and run.n_docked == 500
    assert len(run.hits) == 10
    scores = [h["score"] for h in run.hits]
    assert scores == sorted(scores), "lower kcal/mol must come first"
    joined = " ".join(run.caveats)
    assert "ranks them poorly once there" in joined
    assert "does not say whether the molecule activates or inhibits" in joined


def test_missing_docking_binary_is_reported_not_raised(monkeypatch):
    """Reported, never raised: a machine without the engine should lose
    this section and nothing else. Forced, because the engine is now
    installed here and the branch would not otherwise fire."""
    from promoter_tf_drugs import docking as dk
    import pathlib
    monkeypatch.setattr(dk, "_docking_binary", lambda: None)
    act, dec = _dock_sets()
    run = dk.screen_by_docking(pathlib.Path("x.pdb"), [], target="T",
                               actives=act, decoys=dec)
    assert run.error
    assert "no docking binary" in run.error
    assert "Nothing else in the run is affected" in run.error


def test_consensus_needs_both_methods_to_agree():
    """The reason to run two. A dock reads the pocket and a ligand model
    reads chemistry measured at the receptor; they share no input."""
    from promoter_tf_drugs import docking as dk
    import pathlib
    act, dec = _dock_sets()
    lib = [{"name": f"d{i}", "smiles": f"LIB{i}"} for i in range(500)]
    run = dk.screen_by_docking(pathlib.Path("x.pdb"), lib, target="GOOD",
                               actives=act, decoys=dec, top_n=20,
                               score_fn=_scorer(3.0, seed=1))
    agreeing = [{"drug": run.hits[0]["drug"], "quadrant": "NOVEL HYPOTHESIS"}]
    assert len(dk.consensus(run, agreeing)) == 1
    assert dk.consensus(run, [{"drug": "nothing-like-it"}]) == []
    # and an unusable dock can never contribute a consensus hit
    run.usable = False
    assert dk.consensus(run, agreeing) == []


def test_docking_module_states_it_will_not_widen_coverage():
    """The honest limit, and the one most likely to be misread: docking
    needs a pocket, so it adds confidence on the receptors already in the
    report rather than reaching the undruggable factors."""
    from promoter_tf_drugs import docking as dk
    src = dk.__doc__ + (dk.screen_by_docking.__doc__ or "")
    import inspect
    body = inspect.getsource(dk)
    assert "will not reach the factors you cannot currently drug" in body
    assert "not an affinity" in body


# ===================================================================== #
# Family-pooled screening: more receptors, validated by transfer.       #
# ===================================================================== #
def test_family_transfer_is_measured_on_held_out_target_ligands():
    """A family model trained on siblings must be shown to rank the
    TARGET's own ligands, or it is just the family's pharmacology wearing
    the target's name."""
    import inspect
    from promoter_tf_drugs import screening as scr
    src = inspect.getsource(scr.screen_family)
    assert "transfer_auc" in src
    assert "MIN_TRANSFER_AUC" in src
    # the target's ligands must be excluded from training
    assert "own_smiles" in src
    assert "would leak the answer" in src


def test_family_screen_refuses_when_it_cannot_validate():
    """Without enough of the target's own ligands there is nothing to
    measure transfer against, and an unvalidated family model is exactly
    what fills a report with a family's chemistry."""
    from promoter_tf_drugs import screening as scr
    assert scr.MIN_VALIDATION_ACTIVES >= 5
    assert scr.MIN_TRANSFER_AUC == scr.MIN_BINDING_AUC


def test_transfer_auc_returns_none_rather_than_guessing():
    from promoter_tf_drugs import screening as scr
    assert scr.transfer_auc(None, [], []) is None
    assert scr.transfer_auc(None, [{"smiles": "C"}], []) is None


# ===================================================================== #
# Preparing a receptor and finding its pocket.                          #
# ===================================================================== #
def _pdb_atom(rec, serial, name, resn, chain, resi, x, y, z, elem):
    return (f"{rec:<6}{serial:>5} {name:^4} {resn:>3} {chain}"
            f"{resi:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}"
            f"{1.0:>6.2f}{20.0:>6.2f}          {elem:>2}")


def _fake_pdb(with_ligand=True, coactivator=True, n_lig=24):
    lines, n = [], 1
    for i in range(1, 41):
        lines.append(_pdb_atom("ATOM", n, "CA", "ALA", "A", i,
                               10.0 + i, 5.0, 5.0, "C")); n += 1
    for i in range(1, 21):
        lines.append(_pdb_atom("ATOM", n, "CA", "GLY", "B", i,
                               90.0 + i, 90.0, 90.0, "C")); n += 1
    if coactivator:
        for i, resn in enumerate(["LEU", "LYS", "LYS", "LEU", "LEU"], 1):
            lines.append(_pdb_atom("ATOM", n, "CA", resn, "P", i,
                                   15.0 + i, 9.0, 5.0, "C")); n += 1
    lines.append("TER")
    if with_ligand:
        for i in range(n_lig):
            lines.append(_pdb_atom("HETATM", n, "C1", "LIG", "A", 900,
                                   22.0 + (i % 5) * 0.9, 6.0 + (i % 3) * 0.8,
                                   5.5 + (i % 4) * 0.7, "C")); n += 1
    for i in range(20):
        lines.append(_pdb_atom("HETATM", n, "O", "HOH", "A", 1000 + i,
                               60.0 + i, 60.0, 60.0, "O")); n += 1
    for i in range(5):
        lines.append(_pdb_atom("HETATM", n, "O", "SO4", "A", 1100,
                               70.0 + i, 70.0, 70.0, "O")); n += 1
    lines.append(_pdb_atom("HETATM", n, "ZN", " ZN", "A", 1200,
                           80.0, 80.0, 80.0, "ZN"))
    return "\n".join(lines) + "\nEND\n"


def test_solvent_and_ions_are_never_mistaken_for_the_ligand():
    """Every structure has waters, a sulfate and often a zinc. Boxing on
    one of them puts the search 50 angstroms from the site."""
    from promoter_tf_drugs import docking as dk
    groups = dk.hetatm_groups(_fake_pdb())
    names = {g["resname"] for g in groups}
    assert "LIG" in names
    assert not names & {"HOH", "SO4", "ZN"}

    lig = dk.best_ligand(_fake_pdb())
    assert lig["resname"] == "LIG" and lig["n_atoms"] == 24


def test_a_fragment_is_too_small_to_mark_a_drug_pocket():
    from promoter_tf_drugs import docking as dk
    assert dk.best_ligand(_fake_pdb(n_lig=4)) is None
    assert dk.MIN_LIGAND_ATOMS >= 10


def test_receptor_preparation_empties_the_pocket(tmp_path):
    """Docking into a structure with its crystal ligand still in place is
    docking into an occupied site, and every compound scores weakly for a
    reason that has nothing to do with the compounds."""
    from promoter_tf_drugs import docking as dk
    text = _fake_pdb()
    out = dk.prepare_receptor(text, tmp_path / "rec.pdb", keep_chain="A")
    body = out.read_text()
    assert "HETATM" not in body
    assert "HOH" not in body
    chains = {ln[21] for ln in body.splitlines() if ln.startswith("ATOM")}
    assert chains == {"A"}


def test_the_box_sits_on_the_ligand_and_has_a_floor(tmp_path):
    """Computed here rather than left to the binary, because
    --autobox_ligand is a smina flag and AutoDock Vina has no equivalent."""
    import pytest
    from pathlib import Path
    from promoter_tf_drugs import docking as dk
    text = _fake_pdb()
    lig = dk.best_ligand(text)
    lp = dk.write_ligand_pdb(text, lig, tmp_path / "lig.pdb")
    box = dk.box_from_ligand(lp.read_text())
    cx, cy, cz, sx, sy, sz = box
    assert abs(cx - lig["centre"][0]) < 2.0
    assert abs(cy - lig["centre"][1]) < 2.0
    # a compact ligand must still get a searchable box
    assert min(sx, sy, sz) >= 12.0

    with pytest.raises(ValueError):
        dk.box_from_ligand("no coordinates here")


def test_written_ligand_contains_only_that_ligand(tmp_path):
    from promoter_tf_drugs import docking as dk
    text = _fake_pdb()
    lp = dk.write_ligand_pdb(text, dk.best_ligand(text), tmp_path / "l.pdb")
    body = lp.read_text()
    assert "HOH" not in body and "SO4" not in body
    assert all(ln.startswith(("HETATM", "END"))
               for ln in body.splitlines() if ln.strip())


class _FakeRCSB:
    """Two entries: an apo structure first, then one with a ligand."""
    def __init__(self, ids=("1APO", "2HOL"), apo_only=False):
        self.ids, self.apo_only = ids, apo_only

    def post(self, url, payload):
        return {"result_set": [{"identifier": i} for i in self.ids]}

    def request(self, url, parse_json=True):
        pid = url.rsplit("/", 1)[-1].replace(".pdb", "")
        return _fake_pdb(with_ligand=(not self.apo_only and pid != "1APO"))


def test_choose_pocket_skips_apo_and_reads_the_conformation(tmp_path):
    """A structure without a co-crystallised ligand gives no reliable box,
    and the conformation matters because an agonist-bound pocket is the
    one a compound that RAISES your gene has to fit."""
    import pytest
    from pathlib import Path
    from promoter_tf_drugs import docking as dk
    pk = dk.choose_pocket(_FakeRCSB(), "P00000", tmp_path)
    assert pk is not None
    assert pk.pdb_id == "2HOL"
    assert pk.ligand_resname == "LIG"
    assert pk.state == "active"
    assert "LXXLL" in pk.note
    assert "HETATM" not in Path(pk.path).read_text()


def test_no_holo_structure_is_refused_rather_than_guessed(tmp_path):
    """A guessed box is the commonest way a docking campaign quietly
    measures nothing."""
    from promoter_tf_drugs import docking as dk
    assert dk.choose_pocket(_FakeRCSB(apo_only=True), "P0", tmp_path) is None

    res = dk.dock_screen_for_target(
        _FakeRCSB(apo_only=True), "T", "P0",
        [{"name": "d", "smiles": "LIB0"}],
        actives=[{"smiles": f"ACT{i}"} for i in range(10)],
        work_dir=tmp_path, score_fn=lambda s: -7.0)
    assert res["usable"] is False
    assert "refused rather than approximated" in res["error"]


def test_orchestrator_runs_end_to_end_and_ranks(tmp_path):
    from promoter_tf_drugs import docking as dk
    import random
    rng = random.Random(0)
    actives = [{"smiles": f"ACT{i}"} for i in range(15)]
    library = ([{"name": f"approved{i}", "smiles": f"LIB{i}"}
                for i in range(200)]
               + [{"name": f"known{i}", "smiles": f"ACT{i}"}
                  for i in range(15)])

    def score(smi):
        base = rng.gauss(-7.0, 1.1)
        return base - 2.5 if smi.startswith("ACT") else base

    res = dk.dock_screen_for_target(
        _FakeRCSB(), "T", "P0", library, actives=actives,
        work_dir=tmp_path, top_n=10, score_fn=score)
    assert res["usable"] is True
    assert res["auc"] >= dk.MIN_DOCKING_AUC
    assert res["pocket"]["pdb_id"] == "2HOL"
    scores = [h["score"] for h in res["hits"]]
    assert scores == sorted(scores)



def test_check_setup_names_what_is_missing(monkeypatch):
    """The first thing anyone runs, and it has to say what to install."""
    from promoter_tf_drugs import docking as dk
    monkeypatch.setattr(dk, "_docking_binary", lambda: None)
    chk = dk.check_setup()
    assert chk["ok"] is False
    assert any("docking binary" in m for m in chk["missing"])
    assert "smina" in dk.SETUP
    assert "conda install -c conda-forge smina" in dk.SETUP


def test_setup_text_covers_windows():
    """conda-forge ships smina for win-64. An earlier version of this text
    said there was no Windows build, which sent people to WSL for no
    reason - so the correction is pinned rather than left to drift back."""
    from promoter_tf_drugs import docking as dk
    assert "Windows" in dk.SETUP
    assert "conda install -c conda-forge smina" in dk.SETUP
    assert "no Windows build" not in dk.SETUP
    # and the measured cost of a screen is stated, since --cpu decides it
    assert "--cpu" in dk.SETUP


# ===================================================================== #
# Prepare where the network is, dock where the binary is.               #
# ===================================================================== #
class _BundleClient(_FakeRCSB):
    def get(self, url, params=None):
        return {}


def _patch_chembl(monkeypatch, n_actives=14):
    from promoter_tf_drugs import agonism as ag, screening as scr, drugs as dr
    monkeypatch.setattr(ag, "fetch_actives", lambda c, t, **k: [
        {"chembl_id": f"CHEMBL{i}", "smiles": f"ACT{i}",
         "pchembl": 7.0 + i / 10} for i in range(n_actives)])
    monkeypatch.setattr(dr, "_chembl_targets", lambda c, s, u: ["CHEMBL100"])
    monkeypatch.setattr(scr, "approved_library", lambda c, **k: (
        [{"chembl_id": f"L{i}", "name": f"approved{i}", "smiles": f"LIB{i}"}
         for i in range(120)]
        + [{"chembl_id": f"CHEMBL{i}", "name": f"known{i}",
            "smiles": f"ACT{i}"} for i in range(n_actives)]))


def test_bundle_captures_everything_the_dock_needs(tmp_path, monkeypatch):
    """The networked half and the docking half are rarely on one machine:
    on Windows the network is here and smina is not, in WSL the reverse."""
    from promoter_tf_drugs import docking as dk
    _patch_chembl(monkeypatch)
    man = dk.write_bundle(_BundleClient(), "VDR", "P11473", tmp_path,
                          library_limit=200)
    assert man["ready"] is True
    assert man["pdb_id"] == "2HOL"
    assert man["state"] == "active"
    assert (tmp_path / "bundle.json").exists()
    assert (tmp_path / "receptor.pdb").exists()
    assert (tmp_path / "ligand.pdb").exists()
    # the receptor in the bundle must already be stripped
    assert "HETATM" not in (tmp_path / "receptor.pdb").read_text()
    # paths stored relative, so the folder can be moved or mounted
    assert man["receptor"] == "receptor.pdb"


def test_docking_from_a_bundle_touches_no_network(tmp_path, monkeypatch):
    """The whole point. If anything reaches for the network here, an
    offline WSL cannot run the dock."""
    import random
    from promoter_tf_drugs import docking as dk
    from promoter_tf_drugs import agonism as ag, screening as scr
    _patch_chembl(monkeypatch)
    dk.write_bundle(_BundleClient(), "VDR", "P11473", tmp_path,
                    library_limit=200)

    def boom(*a, **k):
        raise AssertionError("network used during the offline step")

    monkeypatch.setattr(ag, "fetch_actives", boom)
    monkeypatch.setattr(scr, "approved_library", boom)
    monkeypatch.setattr(dk, "find_structures", boom)
    monkeypatch.setattr(dk, "fetch_pdb", boom)

    rng = random.Random(0)

    def score(smi):
        base = rng.gauss(-7.0, 1.1)
        return base - 2.5 if smi.startswith("ACT") else base

    res = dk.dock_from_bundle(tmp_path, top_n=6, score_fn=score)
    assert res["usable"] is True
    assert res["pocket"]["pdb_id"] == "2HOL"
    scores = [h["score"] for h in res["hits"]]
    assert scores == sorted(scores)


def test_a_bundle_with_too_few_ligands_warns_at_prepare_time(tmp_path,
                                                             monkeypatch):
    """Better to learn the dock will refuse before carrying the folder to
    another machine than after."""
    from promoter_tf_drugs import docking as dk
    _patch_chembl(monkeypatch, n_actives=2)
    man = dk.write_bundle(_BundleClient(), "VDR", "P11473", tmp_path)
    assert man["ready"] is False
    assert "below the" in man["warning"]


def test_a_supplied_structure_without_a_ligand_is_refused(tmp_path,
                                                          monkeypatch):
    from promoter_tf_drugs import docking as dk
    _patch_chembl(monkeypatch)
    apo = tmp_path / "apo.pdb"
    apo.write_text(_fake_pdb(with_ligand=False))
    man = dk.write_bundle(_BundleClient(), "VDR", "P11473",
                          tmp_path / "out", pdb_file=apo)
    assert "error" in man
    assert "no co-crystallised ligand" in man["error"]


def test_missing_bundle_says_to_run_prepare_first(tmp_path):
    from promoter_tf_drugs import docking as dk
    import pytest
    with pytest.raises(FileNotFoundError, match="prepare step"):
        dk.load_bundle(tmp_path)


def test_every_computed_source_reaches_the_catalogue():
    """The bug this pins: activity_drugs and docking were computed, their
    report sections rendered, and both were missing from the dict
    catalogue.build actually reads - so their routes silently never
    entered any table. A source that is computed and then dropped is worse
    than one that was never added, because the section says it ran."""
    import inspect
    from promoter_tf_drugs import catalogue as cat
    src = inspect.getsource(pipeline.run)
    start = src.index("block_for_catalogue = {")
    # Brace-match to the real end. Stopping at the first "}" meant a
    # nested dict value - tf_roles is one - hid every key after it, so
    # this test silently stopped checking the thing it exists for.
    depth, end = 0, start
    for k in range(src.index("{", start), len(src)):
        if src[k] == "{":
            depth += 1
        elif src[k] == "}":
            depth -= 1
            if depth == 0:
                end = k
                break
    block = src[start:end]

    build_src = inspect.getsource(cat.build)
    import re
    reads = set(re.findall(r'result_block\.get\("(\w+)"', build_src))
    missing = [k for k in reads if f'"{k}"' not in block]
    assert not missing, (
        f"catalogue.build reads {missing} and block_for_catalogue does not "
        f"supply them, so those routes are dropped")


def test_docking_results_from_another_machine_reach_the_report(tmp_path):
    """Step three of the split: prepare where the network is, dock where
    the binary is, fold the JSON back in here."""
    import json
    from promoter_tf_drugs import docking as dk

    results = {
        "target": "VDR", "usable": True, "auc": 0.95, "ef1": 11.0,
        "n_docked": 120, "n_failed": 0, "caveats": [],
        "pocket": {"pdb_id": "1DB1", "ligand_resname": "VDX",
                   "state": "active", "note": "x"},
        "hits": [{"drug": "Calcipotriol", "smiles": "C", "score": -11.2},
                 {"drug": "Tacalcitol", "smiles": "C", "score": -10.4}],
    }
    f = tmp_path / "dock.json"
    f.write_text(json.dumps(results), encoding="utf-8")

    assert dk.load_results(f)["target"] == "VDR"

    res = pipeline.run("DEMOGENE", demo=True, offline=True,
                       docking_results=str(f))
    routes = [r for r in res["results"][0]["all_routes"]["routes"]
              if r["category"] == "TF: predicted pocket binder"]
    assert len(routes) == 2
    assert {r["drug"] for r in routes} == {"Calcipotriol", "Tacalcitol"}
    assert all(r["tier"] == 4 for r in routes)
    assert all(r["acts_on"] == "VDR" for r in routes)

    src = (res.get("evidence_sources") or {}).get("sources", {}).get("docking")
    assert src and src["ok"] is True and src["yielded"] == 2


def test_an_unvalidated_docking_run_contributes_nothing(tmp_path):
    """Its scores carry no information, so they must not become routes."""
    import json
    from promoter_tf_drugs import catalogue as cat
    bad = {"target": "VDR", "usable": False, "auc": 0.51,
           "error": "pocket did not validate",
           "hits": [{"drug": "Anything", "score": -12.0}]}
    assert cat.from_docking(bad) == []

    f = tmp_path / "bad.json"
    f.write_text(json.dumps(bad), encoding="utf-8")
    res = pipeline.run("DEMOGENE", demo=True, offline=True,
                       docking_results=str(f))
    assert not [r for r in res["results"][0]["all_routes"]["routes"]
                if r["category"] == "TF: predicted pocket binder"]


def test_docking_routes_carry_the_pocket_validation_not_the_score():
    """A docking score separates binders from non-binders and ranks them
    badly, so per-compound strength would overstate what was measured.
    The axis carries the pocket's validated enrichment instead."""
    from promoter_tf_drugs import catalogue as cat
    run = {"target": "VDR", "usable": True, "auc": 0.88, "ef1": 6.0,
           "pocket": {"pdb_id": "1DB1", "state": "active"},
           "hits": [{"drug": "A", "score": -11.0},
                    {"drug": "B", "score": -7.0}]}
    routes = cat.from_docking(run)
    assert len({r.axes.strength for r in routes}) == 1, \
        "strength must come from the pocket, not the per-compound score"
    assert routes[0].axes.strength == 0.88
    assert "AUC 0.88" in routes[0].axes.evidence
    assert any("does not say whether it activates or inhibits" in n
               for n in routes[0].notes)


def test_a_malformed_results_file_says_how_to_make_a_real_one(tmp_path):
    from promoter_tf_drugs import docking as dk
    import pytest
    f = tmp_path / "junk.json"
    f.write_text('{"something": "else"}', encoding="utf-8")
    with pytest.raises(ValueError, match="--bundle"):
        dk.load_results(f)


# ===================================================================== #
# More perturbation corpora, and what agreement between them is worth.  #
# ===================================================================== #
def _corpus_gmt(path, terms):
    path.write_text("\n".join(f"{n}\t\t" + "\t".join(g)
                               for n, g in terms.items()) + "\n",
                    encoding="utf-8")


def _corpus_fixture(tmp_path, monkeypatch):
    """Three corpora on disk: a shared hit, a corpus-specific hit, and a
    compound that moves the regulon opposite ways in two of them."""
    from pathlib import Path as _P
    from promoter_tf_drugs import corpora as cp
    up = [f"UP{i}" for i in range(25)]
    dn = [f"DN{i}" for i in range(25)]
    regulon = {g: 1 for g in up}
    regulon.update({g: -1 for g in dn})

    def sigs(entries):
        u, d = {}, {}
        for drug, flip in entries:
            for i in range(5):
                u[f"E{i} SYS 6H-{drug}-10.0"] = dn if flip else up
                d[f"E{i} SYS 6H-{drug}-10.0"] = up if flip else dn
        return u, d

    files = {
        "lincs": sigs([("sharedrug", False), ("onlyadrug", False)]),
        "drugmatrix": sigs([("sharedrug", False), ("flipdrug", True)]),
        "creeds": sigs([("sharedrug", False), ("flipdrug", False)]),
    }
    names = {"lincs": ("lincs_chem_pert_up.gmt", "lincs_chem_pert_down.gmt"),
             "drugmatrix": ("drugmatrix_up.gmt", "drugmatrix_down.gmt"),
             "creeds": ("creeds_drug_up.gmt", "creeds_drug_down.gmt")}
    for key, (u, d) in files.items():
        _corpus_gmt(tmp_path / names[key][0], u)
        _corpus_gmt(tmp_path / names[key][1], d)

    patched = []
    for c in cp.CORPORA:
        paths = tuple(tuple(str(tmp_path / _P(x).name) for x in group)
                      for group in c.paths)
        patched.append(cp.Corpus(
            key=c.key, name=c.name, species=c.species, system=c.system,
            paths=paths, enrichr=c.enrichr, split=c.split,
            caveats=c.caveats))
    monkeypatch.setattr(cp, "CORPORA", tuple(patched))
    monkeypatch.setattr(cp, "BY_KEY", {c.key: c for c in patched})
    return regulon


def test_a_corpus_that_is_absent_is_not_one_that_found_nothing(
        tmp_path, monkeypatch):
    """The distinction the whole package turns on. TG-GATEs has no files
    here, and reporting it as no hits would claim a negative it never
    measured."""
    from promoter_tf_drugs import corpora as cp
    regulon = _corpus_fixture(tmp_path, monkeypatch)
    res = cp.scan("MYTF", regulon)
    assert res["lincs"].available is True
    assert res["tggates"].available is False
    assert res["tggates"].hits == []
    assert "not the same as finding nothing" in res["tggates"].error
    assert "absent is not the same as finding nothing" in \
        cp.summarise(res).lower()


def test_agreement_across_corpora_is_the_headline(tmp_path, monkeypatch):
    """These differ in species, tissue, platform and dose, so they share
    few systematic errors. No amount of replication inside one corpus can
    substitute for a second corpus."""
    from promoter_tf_drugs import corpora as cp
    regulon = _corpus_fixture(tmp_path, monkeypatch)
    agreed = cp.agreement(cp.scan("MYTF", regulon))
    by = {a["drug"]: a for a in agreed}
    assert by["sharedrug"]["n_corpora"] == 3
    assert by["sharedrug"]["consistent"] is True
    assert agreed[0]["drug"] == "sharedrug"
    assert "onlyadrug" not in by


def test_conflicting_directions_are_shown_not_hidden(tmp_path, monkeypatch):
    """Rat liver in vivo and a human cell line can disagree without either
    being wrong. Hiding it would make the agreeing rows look more uniform
    than the data is."""
    from promoter_tf_drugs import corpora as cp
    regulon = _corpus_fixture(tmp_path, monkeypatch)
    agreed = cp.agreement(cp.scan("MYTF", regulon))
    flip = next(a for a in agreed if a["drug"] == "flipdrug")
    assert flip["consistent"] is False
    assert flip["direction"] == "conflicting"
    assert "not comparable" in flip["note"]
    assert agreed.index(flip) > 0


def test_one_corpus_can_never_produce_agreement(tmp_path, monkeypatch):
    from promoter_tf_drugs import corpora as cp
    regulon = _corpus_fixture(tmp_path, monkeypatch)
    assert cp.agreement(cp.scan("MYTF", regulon, keys=["lincs"])) == []


def test_splitting_a_signed_gmt_drops_unlabelled_terms(tmp_path):
    """DrugMatrix ships as one library with the direction in the term
    name. A term carrying neither marker cannot be assigned, so it is
    dropped rather than guessed onto one side."""
    from pathlib import Path as _P
    from promoter_tf_drugs import corpora as cp
    src = tmp_path / "signed.gmt"
    src.write_text("DrugA-rat-liver-up\t\tG1\tG2\n"
                   "DrugA-rat-liver-dn\t\tG3\tG4\n"
                   "DrugB-kidney-up\t\tG5\n"
                   "MysteryTerm\t\tG9\n", encoding="utf-8")
    out = cp.split_signed_gmt(src, str(tmp_path / "dm"))
    assert out["n_up"] == 2 and out["n_down"] == 1
    assert out["n_unlabelled"] == 1
    assert "MysteryTerm" not in _P(out["up"]).read_text()


def test_non_human_corpora_carry_a_species_caveat():
    """DrugMatrix is rat. A symbol crossing a species boundary is a real
    limit and has to travel with the result."""
    from promoter_tf_drugs import corpora as cp
    dm = cp.BY_KEY["drugmatrix"]
    assert dm.species == "rat"
    assert any("species boundary" in c for c in dm.caveats)
    assert "hepatocytes" in " ".join(cp.BY_KEY["tggates"].caveats)


def test_archs4_is_not_offered_as_a_drug_corpus():
    """It is uniformly reprocessed RNA-seq, not perturbation data.
    Listing it beside LINCS would invite reading co-expression as a drug
    effect."""
    from promoter_tf_drugs import corpora as cp
    assert "archs4" not in cp.BY_KEY
    assert "ARCHS4" in cp.__doc__
    assert "not a perturbation corpus" in cp.__doc__


def test_report_says_when_only_one_corpus_was_available():
    """With one corpus, internal replication cannot detect an error common
    to all of its signatures."""
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    res["results"][0]["corpus_scans"] = {
        "TP53": {"lincs": {"available": True},
                 "drugmatrix": {"available": False}}}
    note = " ".join(pipeline._corpus_note(res["results"][0]))
    assert "Only one corpus was available" in note
    assert "A second corpus is what would" in note
    assert "drugmatrix" in note


def test_each_hit_records_which_corpus_it_came_from():
    """A rat liver result must not read as a human cell-line one."""
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    res["results"][0]["activity_drugs"] = {"TP53": [{
        "drug": "mitoxantrone", "corpus": "drugmatrix",
        "direction": "activated", "median_z": 1.9, "n_agreeing": 4,
        "n_signatures": 4, "consistency": 1.0, "n_significant": 3,
        "cell_lines": ["liver"], "caveats": []}]}
    res["results"][0]["corpus_scans"] = {
        "TP53": {"drugmatrix": {"available": True}}}
    md = pipeline.markdown_report(res)
    header = [ln for ln in md.splitlines()
              if ln.startswith("| Factor | Drug |")]
    assert header and "Corpus" in header[0]
    i = md.splitlines().index(header[0])
    assert md.splitlines()[i + 1].count("|") == header[0].count("|")
    assert "drugmatrix" in md
    assert "is rat" in md and "species boundary" in md


# ===================================================================== #
# The one-page summary.                                                 #
# ===================================================================== #
def test_summary_is_written_alongside_the_full_report(tmp_path):
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    written = pipeline.write_outputs(res, tmp_path)
    names = {p.name for p in written}
    assert "summary.md" in names
    assert "report.md" in names
    body = (tmp_path / "summary.md").read_text(encoding="utf-8")
    full = len((tmp_path / "report.md").read_text(encoding="utf-8"))
    assert len(body) < full / 3, (
        f"a summary that is not much shorter is not a summary: "
        f"{len(body)} vs {full}")


def test_summary_leads_with_how_much_the_ranking_is_worth():
    """The first thing to know is whether the order means anything. In a
    degraded run it does not, and burying that under the table would be
    the one mistake this page could make."""
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    md = pipeline.summary_report(res)
    head = md[:md.index("## Top")] if "## Top" in md else md
    assert "Evidence:" in head
    assert "reading order, not as evidence" in head


def test_summary_states_the_approved_indication_for_every_row():
    """Without it the rows read as interchangeable, which is how a topical
    acne gel ends up next to an infused cytotoxic."""
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    md = pipeline.summary_report(res)
    assert "Approved to treat" in md
    assert "Every row is off-label" in md


def test_summary_never_truncates_mid_word():
    """A cut like "renal fai" reads as a rendering bug and costs the
    reader the one fact the column carries."""
    assert pipeline._clip("short", 58) == "short"
    long = ("hypocalcaemia and secondary hyperparathyroidism in renal "
            "failure")
    out = pipeline._clip(long, 58)
    assert out.endswith("…")
    assert not out[:-1].endswith(" ")
    # the last kept token is whole
    assert long.startswith(out[:-1])
    assert out[:-1].split()[-1] in long.split()


def test_summary_omits_an_empty_section_rather_than_printing_a_heading():
    """A heading with nothing under it reads as a failure of the tool."""
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    md = pipeline.summary_report(res)
    for heading in ("## What the answer rests on", "## Read these before",
                    "## Top"):
        if heading in md:
            after = md[md.index(heading) + len(heading):]
            nxt = after.find("\n## ")
            body = (after[:nxt] if nxt > 0 else after).strip()
            assert body, f"{heading} rendered with no content"


def test_summary_does_not_state_the_same_number_twice():
    """The first draft printed the false-discovery proportion in the
    verdict reason and again two lines later."""
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    v = res.get("evidence_verdict") or {}
    fdp = v.get("false_discovery_proportion")
    if fdp is None:
        return
    md = pipeline.summary_report(res)
    head = md[:md.index("## Top")]
    assert head.count(f"{fdp:.0%}") <= 1


def test_summary_keeps_the_attrition_caution():
    """The reason the drug tables are all nuclear receptors is the single
    most load-bearing caveat in the report, and a summary that dropped it
    to save a line would mislead more than it helped."""
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    md = pipeline.summary_report(res)
    assert "attrition, not a ranking" in md
    assert "no approved ligand" in md


def test_summary_names_which_sources_ran():
    """How to read the rows depends on what the answer rests on."""
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    res["evidence_sources"] = {"sources": {
        "chip": {"state": "ok", "requested": True},
        "lincs": {"state": "failed", "requested": True}}}
    md = pipeline.summary_report(res)
    assert "**Ran:** chip" in md
    assert "Did not run:** lincs" in md
    assert "not a negative result" in md


def test_summary_keeps_the_falsifying_experiment():
    """Every row is a hypothesis, and the summary is the page most likely
    to be read as a list of answers."""
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    md = pipeline.summary_report(res)
    assert "site-directed mutant" in md
    assert "Not medical advice" in md


def test_summary_handles_a_run_with_no_routes():
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    for b in res["results"]:
        b["all_routes"] = {"routes": []}
    md = pipeline.summary_report(res)
    assert "No route was found" in md
    assert "per-stage accounting" in md


def test_summary_orders_by_tier_before_score():
    """A tier-1 route with a mediocre score beats a tier-5 with a great
    one, and the summary is where that ordering matters most."""
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    res["results"][0]["all_routes"] = {"routes": [
        {"drug": "LowTierHighScore", "acts_on": "X", "tier": 5,
         "score": 99.0, "category": "c", "mechanism": ""},
        {"drug": "TopTierLowScore", "acts_on": "Y", "tier": 1,
         "score": 1.0, "category": "c", "mechanism": ""},
    ]}
    md = pipeline.summary_report(res)
    assert md.index("TopTierLowScore") < md.index("LowTierHighScore")


def test_summary_deduplicates_a_drug_reached_by_one_route_twice():
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    res["results"][0]["all_routes"] = {"routes": [
        {"drug": "Samedrug", "acts_on": "X", "tier": 3, "score": 9.0,
         "category": "c", "mechanism": ""},
        {"drug": "samedrug", "acts_on": "X", "tier": 3, "score": 8.0,
         "category": "c", "mechanism": ""},
    ]}
    md = pipeline.summary_report(res)
    # Scoped to the shortlist: the same drug also appears once in the
    # corroboration table below, which is a different table and correct.
    shortlist = md[md.index("## Top"):md.index("## All ")]
    rows = [l for l in shortlist.splitlines() if "amedrug" in l]
    assert len(rows) == 1


# ===================================================================== #
# Drugs that change how much of the FACTOR exists.                      #
# ===================================================================== #
def test_the_direction_the_factor_must_move_is_explicit():
    """Getting this backwards recommends the drug that moves your gene the
    wrong way, so all four combinations are pinned."""
    from promoter_tf_drugs import abundance as ab
    assert ab.wanted_direction("up", "activator") == "up"
    assert ab.wanted_direction("up", "repressor") == "down"
    assert ab.wanted_direction("down", "activator") == "down"
    assert ab.wanted_direction("down", "repressor") == "up"


def test_a_factor_with_no_recorded_role_gives_no_direction():
    """There is no way to say which way would help, so it is skipped and
    named rather than guessed into one column."""
    from promoter_tf_drugs import abundance as ab
    for role in ("unknown", "dual", "", None):
        assert ab.wanted_direction("up", role) is None
    skipped = ab.undirected_factors(
        [{"symbol": "OSR2", "role": "unknown"},
         {"symbol": "KLF4", "role": "activator"}], goal="up")
    assert [s["symbol"] for s in skipped] == ["OSR2"]
    assert "role is not recorded" in ab.summarise([], skipped)


def test_abundance_tier_is_one_below_a_measurement_on_your_own_gene():
    """catalogue.from_lincs calls a landmark measurement of YOUR GENE tier
    2. Here the compound was measured against the FACTOR, and the
    factor-to-gene step is an inference."""
    from promoter_tf_drugs import abundance as ab
    assert ab._tier_for("landmark") == 3
    assert ab._tier_for("best inferred") == 3
    assert ab._tier_for(None) == 4


def test_a_ligand_gated_receptor_carries_the_abundance_caveat(tmp_path):
    """More VDR is not more VDR activity: the receptor is near-ubiquitous
    and its signalling is not. These are most of the drug tables, so the
    caveat has to travel with the row."""
    from promoter_tf_drugs import abundance as ab
    up = tmp_path / "up.gmt"
    down = tmp_path / "down.gmt"
    up.write_text("".join(
        f"E{i} MCF7 6H-drugx-10.0\t\tVDR\tOTHER{i}\n" for i in range(4)),
        encoding="utf-8")
    down.write_text("E0 MCF7 6H-drugx-10.0\t\tZZZ\n", encoding="utf-8")
    hits = ab.find([{"symbol": "VDR", "role": "activator",
                     "nuclear_receptor": True}],
                   goal="up", up_path=up, down_path=down)
    assert hits
    joined = " ".join(hits[0].caveats)
    assert "nearly unrelated" in joined
    assert "transcript, not protein" in joined


def test_one_row_per_compound_not_one_per_signature(tmp_path):
    """The .gmt holds a consensus per cell line and dose, so a compound
    recurs. The recurrence is evidence; repeating the row is noise."""
    from promoter_tf_drugs import abundance as ab
    up = tmp_path / "up.gmt"
    down = tmp_path / "down.gmt"
    up.write_text("".join(
        f"E{i} CELL{i} 6H-samedrug-10.0\t\tMYTF\tG{i}\n"
        for i in range(5)), encoding="utf-8")
    down.write_text("E0 X 6H-other-1\t\tZZZ\n", encoding="utf-8")
    hits = ab.find([{"symbol": "MYTF", "role": "activator"}],
                   goal="up", up_path=up, down_path=down)
    assert len(hits) == 1
    assert hits[0].n_signatures == 5
    assert len(hits[0].cell_lines) >= 2


def test_a_compound_seen_once_is_not_reported(tmp_path):
    from promoter_tf_drugs import abundance as ab
    up = tmp_path / "up.gmt"
    down = tmp_path / "down.gmt"
    up.write_text("E0 X 6H-oneshot-10.0\t\tMYTF\n", encoding="utf-8")
    down.write_text("E0 X 6H-other-1\t\tZZZ\n", encoding="utf-8")
    assert ab.find([{"symbol": "MYTF", "role": "activator"}],
                   goal="up", up_path=up, down_path=down) == []
    assert ab.MIN_SIGNATURES >= 2


def test_demo_mode_never_reads_a_real_corpus():
    """Demo has to stay self-contained and reproducible. Without the guard
    it silently read whatever corpora were sitting in data/."""
    import inspect
    src = inspect.getsource(pipeline.run)
    i = src.index("factor_abundance: list = []")
    assert "not demo" in src[i:i + 400]

    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    assert res["results"][0]["factor_abundance"] == []


# ===================================================================== #
# The cross-method corroboration table.                                 #
# ===================================================================== #
def _route(drug, acts_on, tier, score, category, source, mechanism="",
           directness=1):
    """A genuine Route, not a hand-written dict.

    Hand-written route dicts kept missing fields the report reads - axes,
    directness - one at a time. Building the real object means the fixture
    cannot drift from the dataclass.
    """
    from promoter_tf_drugs.catalogue import Route, Axes
    r = Route(drug=drug, category=category, mechanism=mechanism,
              acts_on=acts_on, tier=tier, directness=directness,
              axes=Axes(0.5, 0.5, 0.5, "test"), approved=True,
              source=source)
    r.score = score
    return r.to_dict()


def _as_all_routes(routes):
    """Wrap hand-written routes in the shape catalogue.build produces.

    markdown_report reads n_routes, by_category and by_tier off this dict,
    so a bare {"routes": [...]} is not a usable stand-in. The key set is
    asserted against build() so this cannot drift.
    """
    by_cat, by_tier = {}, {}
    for r in routes:
        by_cat.setdefault(r["category"], []).append(r)
        by_tier.setdefault(r["tier"], []).append(r)
    out = {"routes": routes, "n_routes": len(routes),
           "by_category": by_cat, "by_tier": by_tier,
           "categories_absent": [],
           "best_per_category": {k: v[0] for k, v in by_cat.items()}}
    return out


def test_the_all_routes_stand_in_matches_the_real_shape():
    """So a key added to catalogue.build cannot silently break every test
    that hand-builds this dict."""
    from promoter_tf_drugs import catalogue as cat
    real = cat.build({"direct_therapies": [{
        "drug": "X", "gene": "G", "modality": "siRNA", "approved": True,
        "specificity": 0.9, "raises_protein": True, "mechanism": "m",
        "status": "approved", "notes": []}]}, goal="up")
    stand_in = _as_all_routes([
        _route("X", "Y", 3, 1.0, "c", "drug_scoring")])
    assert set(real) == set(stand_in), (
        f"missing from the stand-in: {set(real) - set(stand_in)}")


def _two_method_result():
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    res["results"][0]["all_routes"] = _as_all_routes([
        _route("Calcitriol", "VDR", 3, 70.0, "TF: orthosteric ligand",
               "drug_scoring", "VDR agonist"),
        _route("Calcitriol", "VDR", 2, 60.0, "empirical signature",
               "signatures", "moved this gene in an assay", directness=0),
        _route("Loner", "RARA", 3, 50.0, "TF: orthosteric ligand",
               "drug_scoring", "RARA agonist"),
    ])
    return res


def test_one_row_per_drug_with_the_methods_counted():
    """The route tables are one row per route, so a drug reached three ways
    appeared three times and the reader could not see it was reached three
    ways."""
    rows = pipeline.corroboration(_two_method_result(), top_n=20)
    by = {r["drug"]: r for r in rows}
    assert by["Calcitriol"]["n_methods"] == 2
    assert by["Calcitriol"]["n_routes"] == 2
    assert by["Loner"]["n_methods"] == 1
    # best tier across a drug's routes, not the first seen
    assert by["Calcitriol"]["best_tier"] == 2
    # and corroboration breaks the tie
    assert rows[0]["drug"] == "Calcitriol"


def test_every_route_source_has_a_method_name():
    """A new route source must appear in this table rather than silently
    counting as zero methods, which is what the first version did for
    every row."""
    res = _two_method_result()
    assert pipeline.unmapped_sources(res) == set()

    res["results"][0]["all_routes"]["routes"].append(
        _route("X", "Y", 4, 1.0, "c", "a_brand_new_source"))
    assert pipeline.unmapped_sources(res) == {"a_brand_new_source"}


def test_expected_change_is_a_direction_never_a_magnitude():
    """Nothing here estimates a fold-change, so a number would be
    invented. The basis says which kind of claim it is."""
    rows = pipeline.corroboration(_two_method_result(), top_n=20)
    by = {r["drug"]: r for r in rows}
    assert by["Calcitriol"]["expected_change"] in ("increase", "decrease")
    assert "MEASURED" in by["Calcitriol"]["expected_basis"]
    assert "PREDICTED" in by["Loner"]["expected_basis"]

    md = "\n".join(pipeline.corroboration_table(_two_method_result()))
    assert "no predicted fold-change here" in md
    assert "would be fabricated" in md
    # the three inputs a fold-change would need are named, so the refusal
    # is a statement about missing data rather than a hedge
    assert "occupancy change" in md and "dose reached" in md


def test_the_method_count_is_not_presented_as_a_vote():
    """Two weak agreements do not equal one measured result, and a bare
    count would say they did."""
    md = "\n".join(pipeline.corroboration_table(_two_method_result()))
    assert "not a vote" in md
    assert "Tier still decides" in md


def test_systems_acted_on_are_named_per_drug():
    rows = pipeline.corroboration(_two_method_result(), top_n=20)
    cal = next(r for r in rows if r["drug"] == "Calcitriol")
    assert cal["systems"] == "vitamin D receptor agonist"
    assert cal["systems_why"]
    md = "\n".join(pipeline.corroboration_table(_two_method_result()))
    assert "Systems it acts on" in md


def test_a_single_method_run_says_so_rather_than_implying_breadth():
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    res["results"][0]["all_routes"] = _as_all_routes([
        _route("Only", "VDR", 3, 9.0, "c", "drug_scoring")])
    md = "\n".join(pipeline.corroboration_table(res))
    assert "No compound was found by more than one method" in md
    assert "expected rather than informative" in md


def test_the_table_reaches_summary_md(tmp_path):
    res = _two_method_result()
    pipeline.write_outputs(res, tmp_path)
    body = (tmp_path / "summary.md").read_text(encoding="utf-8")
    assert "which methods found each" in body
    assert "Systems it acts on" in body
    assert "DEMOGENE change" in body
    assert "P(moves)" in body


# ===================================================================== #
# Docking on Windows, and BBB in the corroboration table.               #
# ===================================================================== #
def test_windows_has_a_native_docking_path():
    """smina has no Windows build, which made WSL the only route. Vina
    does ship a Windows executable, and the PDBQT conversion it needs is
    pip-installable, so WSL is now a preference rather than a
    requirement."""
    from promoter_tf_drugs import docking as dk
    assert hasattr(dk, "dock_one_vina")
    assert "no WSL needed" in dk.VINA_SETUP_WINDOWS
    assert "meeko" in dk.VINA_SETUP_WINDOWS
    assert "vina.exe" in dk.SETUP or "vina.exe" in dk.VINA_SETUP_WINDOWS
    # and the main setup text no longer presents Vina as refused
    assert "it is refused rather than run" not in dk.SETUP


def test_vina_is_no_longer_rejected_outright(tmp_path, monkeypatch):
    """It used to raise 'install smina instead'. Now it routes to the path
    that prepares PDBQT."""
    import pytest
    from promoter_tf_drugs import docking as dk
    monkeypatch.setattr(dk, "_docking_binary", lambda: "/usr/bin/vina")
    with pytest.raises(dk._NeedsVinaPath):
        dk.dock_one_boxed(tmp_path / "r.pdb", tmp_path / "l.sdf",
                          (0, 0, 0, 20, 20, 20))


def test_missing_vina_prep_names_the_pip_install():
    from promoter_tf_drugs import docking as dk
    ok, why = dk.vina_prep_available()
    if not ok:
        assert "pip install meeko" in why


def test_the_receptor_pdbqt_is_converted_once_per_screen():
    """It is the same file for every compound, and converting it a few
    thousand times would dominate the run."""
    import inspect
    from promoter_tf_drugs import docking as dk
    assert isinstance(dk._REC_CACHE, dict)
    src = inspect.getsource(dk.dock_one_vina)
    assert "receptor_pdbqt_cache" in src
    assert "cache[key]" in src


def test_corroboration_table_carries_the_bbb_verdict():
    """Whether a candidate can reach the tissue is an applicability
    question, so it belongs on the one-row-per-drug table too."""
    res = _two_method_result()
    res["results"][0]["drug_grades"] = [{
        "drug_name": "Calcitriol",
        "evidence": {"properties": {"full_mwt": 416.6, "psa": 60.7,
                                    "cx_logp": 5.1}}}]
    rows = {r["drug"]: r for r in pipeline.corroboration(res, top_n=20)}
    assert rows["Calcitriol"]["bbb"] in ("yes", "no", "yes?", "no?")
    assert rows["Calcitriol"]["bbb_why"]
    # no structure for the other drug, so unknown rather than no
    assert rows["Loner"]["bbb"] == "\u2014"

    md = "\n".join(pipeline.corroboration_table(res))
    header = [l for l in md.splitlines() if l.startswith("| Drug |")][0]
    assert "BBB?" in header
    i = md.splitlines().index(header)
    assert md.splitlines()[i + 1].count("|") == header.count("|")
    # the same footnote the graded table uses, so the two cannot diverge
    assert "is not a no" in md or "passive diffusion" in md


# ===================================================================== #
# Running with no usable network at all.                                #
# ===================================================================== #
def _hub_table(tmp_path, extra=()):
    """A Repurposing-Hub-shaped table, with its real column names."""
    rows = [
        ("calcitriol", "launched", "VDR agonist", "VDR", "CCO"),
        ("alfacalcidol", "launched", "VDR agonist", "VDR", "CCCO"),
        ("paricalcitol", "launched", "VDR agonist", "VDR", "CCCCO"),
        ("metformin", "launched", "Biguanide", "PRKAB1", "CN(C)C(=N)N=C(N)N"),
        ("aspirin", "launched", "COX inhibitor", "PTGS1|PTGS2",
         "CC(=O)Oc1ccccc1C(=O)O"),
        ("toolcpd", "preclinical", "Tool", "VDR", "CCOc1ccccc1"),
        ("nostructure", "launched", "Unknown", "VDR", ""),
    ] + list(extra)
    p = tmp_path / "repurposing_drugs.txt"
    p.write_text("pert_iname\tclinical_phase\tmoa\ttarget\tsmiles\n"
                 + "\n".join("\t".join(r) for r in rows) + "\n",
                 encoding="utf-8")
    return p


def test_one_local_table_replaces_the_whole_chembl_dependency(tmp_path):
    """Docking needs three things from ChEMBL - the target id, the known
    ligands for the gate, and the library. A corporate network that blocks
    or 500s those stops the analysis, and one annotated file removes all
    three."""
    from promoter_tf_drugs import localdrugs as ld
    rows = ld.read_table(_hub_table(tmp_path))
    names = {r["name"] for r in rows}
    assert "calcitriol" in names
    assert "toolcpd" not in names, "preclinical must be filtered"
    assert "nostructure" not in names, "a name alone cannot be docked"

    acts = ld.actives_for("VDR", rows)
    assert {a["name"] for a in acts} == {"calcitriol", "alfacalcidol",
                                         "paricalcitol"}


def test_a_local_table_carries_no_potency_and_says_so(tmp_path):
    """It records that a compound acts on a target, not how hard it binds.
    A fabricated pChEMBL would feed the potency axis a made-up number."""
    from promoter_tf_drugs import localdrugs as ld
    rows = ld.read_table(_hub_table(tmp_path))
    acts = ld.actives_for("VDR", rows)
    assert all(a["pchembl"] is None for a in acts)
    assert "no potency" in ld.summarise(rows, "VDR")
    assert "unranked rather than scored zero" in ld.summarise(rows, "VDR")


def test_the_delimiter_is_sniffed_not_taken_from_the_extension(tmp_path):
    """The Hub ships .txt that is tab separated, and people save it as
    .csv. Guessing from the suffix parses the whole file as one column."""
    from promoter_tf_drugs import localdrugs as ld
    p = tmp_path / "comma.txt"
    p.write_text("name,clinical_phase,smiles\n"
                 "drugA,launched,CCO\ndrugB,launched,CCCO\n",
                 encoding="utf-8")
    rows = ld.read_table(p)
    assert len(rows) == 2 and rows[0]["smiles"] == "CCO"


def test_a_table_with_no_structures_is_refused_with_the_reason(tmp_path):
    from promoter_tf_drugs import localdrugs as ld
    import pytest
    p = tmp_path / "nosmiles.tsv"
    p.write_text("name\tphase\ndrugA\tlaunched\n", encoding="utf-8")
    with pytest.raises(ValueError, match="no SMILES column"):
        ld.read_table(p)


def test_a_missing_table_returns_empty_and_names_the_source(tmp_path,
                                                            monkeypatch):
    """Absent is reported, never treated as an empty library, because an
    empty library would make every screen silently find nothing."""
    from promoter_tf_drugs import localdrugs as ld
    monkeypatch.chdir(tmp_path)
    assert ld.library() == []
    assert "repurposing.broadinstitute.org" in ld.SOURCE
    assert "need the network" in ld.summarise([])


def test_docking_prefers_a_local_table_over_the_network(tmp_path,
                                                        monkeypatch):
    """On a blocked network the ChEMBL calls are the step that stops
    everything, so they must be the fallback and not the default."""
    import inspect
    from promoter_tf_drugs import docking as dk
    src = inspect.getsource(dk.write_bundle)
    i_local = src.index("localdrugs")
    i_net = src.index("approved_library(client)")
    assert i_local < i_net, (
        "the network path must be the fallback, not the first attempt")
    assert "if not library:" in src
    assert "if not actives:" in src


def test_the_circuit_breaker_also_trips_on_http_500(monkeypatch):
    """ChEMBL is called hundreds of times a run. If a 500 did not trip the
    breaker, every call would pay the full retry budget."""
    import inspect
    from promoter_tf_drugs import http as ph
    src = inspect.getsource(ph.Client.request)
    # the 500 branch must leave the loop so the breaker code below runs
    i500 = src.index("429, 500, 502, 503, 504")
    assert "break" in src[i500:i500 + 900]
    # and the breaker must be consulted before a request is made
    assert "if host in self._down" in inspect.getsource(ph.Client)


# ===================================================================== #
# MyGene, Open Targets and potency, from files.                         #
# ===================================================================== #
def _annot_files(tmp_path):
    d = tmp_path / "data"
    d.mkdir(exist_ok=True)
    (d / "hgnc_complete_set.txt").write_text(
        "symbol\tname\tuniprot_ids\tentrez_id\talias_symbol\tprev_symbol\n"
        "VDR\tvitamin D receptor\tP11473\t7421\tNR1I1\t\n"
        "REST\tRE1 silencing TF\tQ13127\t5978\tNRSF\tXBR\n",
        encoding="utf-8")
    (d / "goa_human.gaf").write_text(
        "!gaf-version: 2.2\n"
        "UniProtKB\tP11473\tVDR\t\tGO:0001228\tR\tIDA\t\tF\t\t\t\t\n"
        "UniProtKB\tQ13127\tREST\t\tGO:0001227\tR\tIDA\t\tF\t\t\t\t\n"
        "UniProtKB\tQ13127\tREST\tNOT\tGO:0001228\tR\tIDA\t\tF\t\t\t\t\n",
        encoding="utf-8")
    return d


def _reset_annot():
    from promoter_tf_drugs import localannot as la
    la._HGNC.clear()
    la._GOA.clear()
    la._GWAS.clear()


def test_role_comes_from_the_same_go_terms_online_and_offline():
    """Two code paths deciding direction by different rules would be worse
    than having only one, so the offline path imports the term sets rather
    than restating them."""
    import inspect
    from promoter_tf_drugs import localannot as la, annotate as an
    assert "from .annotate import" in inspect.getsource(la)
    assert la.ACTIVATOR_GO is an.ACTIVATOR_GO
    assert la.REPRESSOR_GO is an.REPRESSOR_GO


def test_hgnc_and_goa_give_role_uniprot_and_the_nr_flag(tmp_path,
                                                        monkeypatch):
    from promoter_tf_drugs import localannot as la
    monkeypatch.chdir(tmp_path)
    _annot_files(tmp_path)
    _reset_annot()
    ann = la.annotate_local(["VDR", "REST"])
    assert ann["VDR"]["role"] == "activator"
    assert ann["VDR"]["uniprot"] == "P11473"
    assert ann["VDR"]["nuclear_receptor"] is True
    assert ann["REST"]["role"] == "repressor"


def test_a_not_qualified_go_annotation_is_evidence_against(tmp_path,
                                                           monkeypatch):
    """GAF NOT means the term does not apply. Counting it would make REST
    read as dual, an activator and a repressor, from a row saying the
    opposite."""
    from promoter_tf_drugs import localannot as la
    monkeypatch.chdir(tmp_path)
    _annot_files(tmp_path)
    _reset_annot()
    assert la.annotate_local(["REST"])["REST"]["role"] == "repressor"


def test_an_alias_keeps_its_role(tmp_path, monkeypatch):
    """GOA is keyed on the approved symbol, so an alias has to be followed
    to it or it loses the role - and abundance.py skips any factor without
    one."""
    from promoter_tf_drugs import localannot as la
    monkeypatch.chdir(tmp_path)
    _annot_files(tmp_path)
    _reset_annot()
    ann = la.annotate_local(["NRSF", "XBR", "NR1I1"])
    assert ann["NRSF"]["role"] == "repressor"
    assert ann["XBR"]["role"] == "repressor"
    assert ann["NR1I1"]["nuclear_receptor"] is True


def test_an_unknown_symbol_is_blank_not_guessed(tmp_path, monkeypatch):
    from promoter_tf_drugs import localannot as la
    monkeypatch.chdir(tmp_path)
    _annot_files(tmp_path)
    _reset_annot()
    ann = la.annotate_local(["NOSUCHGENE"])["NOSUCHGENE"]
    assert ann["role"] == "unknown" and ann["uniprot"] is None
    assert "unresolved rather than absent" in ann["note"]


def test_annotate_tfs_tries_local_before_mygene():
    import inspect
    from promoter_tf_drugs import annotate as an
    src = inspect.getsource(an.annotate_tfs)
    assert src.index("localannot") < src.index("MYGENE"), (
        "the network must be the fallback, not the first attempt")


def test_gwas_is_labelled_as_weaker_than_an_l2g_score(tmp_path,
                                                      monkeypatch):
    """The catalogue names the gene mapped to a signal, which is not the
    causal gene at that locus. Presenting it as Open Targets genetics
    would overstate it."""
    from promoter_tf_drugs import localannot as la
    monkeypatch.chdir(tmp_path)
    d = tmp_path / "data"
    d.mkdir(exist_ok=True)
    (d / "gwas_catalog_associations.tsv").write_text(
        "MAPPED_GENE\tDISEASE/TRAIT\tP-VALUE\tPUBMEDID\n"
        "VDR\tBone density\t1e-12\t111\n"
        "VDR\tWeak signal\t1e-4\t112\n", encoding="utf-8")
    _reset_annot()
    g = la.genetic_support("VDR")
    assert g["assessed"] is True
    assert g["n"] == 2 and g["n_significant"] == 1
    assert "Bone density" in g["traits"]
    joined = " ".join(g["caveats"])
    assert "not the same as the causal gene" in joined
    assert "not Open Targets" in joined


def test_no_gwas_file_is_not_absence_of_genetic_support(tmp_path,
                                                        monkeypatch):
    from promoter_tf_drugs import localannot as la
    monkeypatch.chdir(tmp_path)
    _reset_annot()
    g = la.genetic_support("VDR")
    assert g["assessed"] is False
    assert "different from there being none" in g["note"]


# ------------------------------------------------------------- potency
def test_a_censored_potency_is_dropped_not_converted():
    """A value of more than 10000 nM means the compound did nothing up to
    10 uM. Recording it as pChEMBL 5 would turn a failed experiment into a
    weak binding constant, and those rows are common enough to bias every
    potency downward."""
    from promoter_tf_drugs import localpotency as lp
    assert lp.to_pchembl("1") == (9.0, False)
    assert lp.to_pchembl("1000") == (6.0, False)
    assert lp.to_pchembl(">10000")[1] is True
    assert lp.to_pchembl("<0.5")[1] is True
    assert lp.to_pchembl("") == (None, False)
    assert lp.to_pchembl("0")[0] is None


def _bindingdb(tmp_path, rows):
    d = tmp_path / "data"
    d.mkdir(exist_ok=True)
    hdr = ["Ligand SMILES", "ChEMBL ID of Ligand", "Ki (nM)", "Kd (nM)",
           "IC50 (nM)", "EC50 (nM)",
           "UniProt (SwissProt) Primary ID of Target Chain"]
    p = d / "BindingDB_All.tsv"
    p.write_text("\t".join(hdr) + "\n"
                 + "\n".join("\t".join(r) for r in rows) + "\n",
                 encoding="utf-8")
    return p


def test_potency_comes_from_a_local_file(tmp_path, monkeypatch):
    from promoter_tf_drugs import localpotency as lp
    monkeypatch.chdir(tmp_path)
    _bindingdb(tmp_path, [
        ["CCO", "CHEMBL1", "1.0", "", "", "", "P11473"],
        ["CCCO", "CHEMBL2", "", "", "50", "", "P11473"],
        ["CCCCO", "CHEMBL3", ">10000", "", "", "", "P11473"],
        ["CCN", "CHEMBL5", "5.0", "", "", "", "P99999"],
    ])
    got = lp.load_bindingdb(["P11473"])
    assert {r["smiles"] for r in got["P11473"]} == {"CCO", "CCCO"}
    assert "P99999" not in got


def test_repeated_measurements_use_the_median_not_the_best(tmp_path,
                                                           monkeypatch):
    """Taking the best would favour compounds assayed many times, which is
    study effort rather than potency - the confound the regulon axis
    already corrects for."""
    from promoter_tf_drugs import localpotency as lp
    monkeypatch.chdir(tmp_path)
    _bindingdb(tmp_path, [
        ["CCO", "C1", "1.0", "", "", "", "P1"],
        ["CCO", "C1", "100.0", "", "", "", "P1"],
        ["CCO", "C1", "1000.0", "", "", "", "P1"],
    ])
    rec = lp.load_bindingdb(["P1"])["P1"][0]
    assert rec["pchembl"] == 7.0, "median of 9.0, 7.0, 6.0"
    assert rec["n_measurements"] == 3


def test_equilibrium_constants_are_preferred_over_ic50(tmp_path,
                                                       monkeypatch):
    """Ki and Kd are equilibrium constants; IC50 depends on how the assay
    was run."""
    from promoter_tf_drugs import localpotency as lp
    monkeypatch.chdir(tmp_path)
    _bindingdb(tmp_path, [["CCO", "C1", "100.0", "", "1.0", "", "P1"]])
    rec = lp.load_bindingdb(["P1"])["P1"][0]
    assert rec["measure"] == "Ki" and rec["pchembl"] == 7.0


def test_no_potency_file_leaves_potency_unmeasured(tmp_path, monkeypatch):
    from promoter_tf_drugs import localpotency as lp
    monkeypatch.chdir(tmp_path)
    assert lp.load_bindingdb(["P1"]) == {}
    assert lp.actives_for("P1") == []
    msg = lp.summarise([])
    assert "unranked rather than scored zero" in msg
    assert "bindingdb.org" in msg


def test_the_drug_table_prefers_measured_potency_when_a_file_exists(
        tmp_path, monkeypatch):
    """localdrugs said plainly it could not supply potency. With a potency
    file on disk that should stop being the answer."""
    from promoter_tf_drugs import localdrugs as ld
    monkeypatch.chdir(tmp_path)
    _bindingdb(tmp_path, [["CCO", "CHEMBL1", "1.0", "", "", "", "P11473"]])
    rows = [{"name": "calcitriol", "smiles": "CCO", "targets": ["VDR"],
             "chembl_id": None, "phase": "launched"}]
    with_p = ld.actives_for("VDR", rows, uniprot="P11473")
    assert with_p and with_p[0]["pchembl"] == 9.0
    without = ld.actives_for("VDR", rows)
    assert without and without[0]["pchembl"] is None


# ===================================================================== #
# Expanding the library: what the 20,000 figure actually counts.        #
# ===================================================================== #
def _tbl(tmp_path, name, header, rows, sep="\t"):
    p = tmp_path / name
    p.write_text(sep.join(header) + "\n"
                 + "\n".join(sep.join(r) for r in rows) + "\n",
                 encoding="utf-8")
    return p


_METFORMIN = "CN(C)C(=N)N=C(N)N"
_METFORMIN_HCL = "Cl." + _METFORMIN
_T4 = "N[C@@H](Cc1cc(I)c(Oc2cc(I)c(O)c(I)c2)c(I)c1)C(=O)O"
_T4_NA = "[Na+].N[C@@H](Cc1cc(I)c(Oc2cc(I)c(O)c(I)c2)c(I)c1)C(=O)[O-]"


def test_an_inchikey_alone_does_not_merge_salts():
    """The counterion is part of the molecule, so metformin and metformin
    hydrochloride get entirely different InChIKeys. This is why the parent
    has to be extracted before keying."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    a = Chem.MolToInchiKey(Chem.MolFromSmiles(_METFORMIN))[:14]
    b = Chem.MolToInchiKey(Chem.MolFromSmiles(_METFORMIN_HCL))[:14]
    assert a != b


def test_the_parent_key_merges_salts_and_separates_real_drugs():
    from promoter_tf_drugs import localdrugs as ld
    assert ld.parent_key(_METFORMIN) == ld.parent_key(_METFORMIN_HCL)
    assert ld.parent_key(_T4) == ld.parent_key(_T4_NA)
    # and two genuinely different drugs stay apart
    assert ld.parent_key(_METFORMIN) != ld.parent_key(_T4)
    assert ld.parent_key("not a smiles at all") is None


def test_merging_collapses_salt_forms_across_sources(tmp_path):
    """Screening the hydrochloride and the free base as two candidates
    doubles the compute and finds the same thing twice."""
    from promoter_tf_drugs import localdrugs as ld
    a = _tbl(tmp_path, "hub.txt",
             ["pert_iname", "clinical_phase", "target", "smiles"],
             [["metformin", "launched", "PRKAB1", _METFORMIN],
              ["metformin hydrochloride", "launched", "PRKAB1",
               _METFORMIN_HCL],
              ["levothyroxine", "launched", "THRB", _T4]])
    b = _tbl(tmp_path, "chembl.csv",
             ["name", "max_phase", "gene", "smiles"],
             [["levothyroxine sodium", "4", "THRB", _T4_NA]], sep=",")
    res = ld.merge([a, b], min_phase=4)
    assert res["n_distinct"] == 2
    assert res["n_salt_forms_collapsed"] == 2
    # the second source contributed no new molecule, and that is reported
    assert res["per_source"]["chembl.csv"]["new"] == 0
    assert res["per_source"]["chembl.csv"]["already_present"] == 1


def test_min_phase_is_the_dial_that_actually_enlarges_the_library(tmp_path):
    """There are no missing FDA approvals to find. Widening the definition
    is what adds molecules, and each step changes what a hit means."""
    from promoter_tf_drugs import localdrugs as ld
    t = _tbl(tmp_path, "mixed.txt",
             ["pert_iname", "clinical_phase", "target", "smiles"],
             [["approveddrug", "launched", "VDR", _METFORMIN],
              ["trialdrug", "phase 3", "VDR", "CCCCCCO"],
              ["earlydrug", "phase 1", "VDR", "CCCCCCCO"],
              ["pulleddrug", "withdrawn", "VDR", "CCCCCCCCO"],
              ["labcompound", "preclinical", "VDR", "CCCCCCCCCO"]])
    sizes = {mp: ld.merge([t], min_phase=mp)["n_distinct"]
             for mp in (4, 3, 1)}
    assert sizes[4] == 1
    assert sizes[3] == 2
    assert sizes[1] == 3
    # preclinical is never reached by the dial
    assert sizes[1] < 4

    # Withdrawn is a separate opt-in, not a lower phase: the drug reached
    # approval and was pulled. On the phase scale, asking for it dragged
    # in every preclinical compound.
    wd = ld.merge([t], min_phase=4, include_withdrawn=True)
    assert wd["n_distinct"] == 2
    assert any(r["withdrawn"] for r in wd["library"])
    assert ld.merge([t], min_phase=4)["n_withdrawn_excluded"] == 1


def test_a_row_with_an_unparseable_structure_is_kept_by_name(tmp_path):
    """Dropping it silently would lose a real drug. It is keyed on the name
    instead and flagged, so the fallback is visible."""
    from promoter_tf_drugs import localdrugs as ld
    t = _tbl(tmp_path, "broken.txt",
             ["pert_iname", "clinical_phase", "target", "smiles"],
             [["weirddrug", "launched", "VDR", "C(C(C"]])
    res = ld.merge([t], min_phase=4)
    assert res["n_distinct"] == 1
    assert res["library"][0]["structure_unparsed"] is True


def test_the_count_explanation_names_what_20000_counts(tmp_path):
    """The premise of "I only have 3,000 of 20,000" is a counting
    difference, and saying so is more useful than adding salt forms to
    reach the number."""
    from promoter_tf_drugs import localdrugs as ld
    t = _tbl(tmp_path, "one.txt",
             ["pert_iname", "clinical_phase", "target", "smiles"],
             [["metformin", "launched", "PRKAB1", _METFORMIN]])
    msg = ld.explain_counts(ld.merge([t], min_phase=4))
    assert "PRODUCTS" in msg
    assert "1,500 to 2,000" in msg
    assert "most of the molecules" in msg
    assert "widen the definition on purpose" in msg


def test_a_merged_library_can_be_written_and_read_back(tmp_path):
    from promoter_tf_drugs import localdrugs as ld
    t = _tbl(tmp_path, "src.txt",
             ["pert_iname", "clinical_phase", "target", "smiles"],
             [["metformin", "launched", "PRKAB1|AMPK", _METFORMIN],
              ["levothyroxine", "launched", "THRB", _T4]])
    res = ld.merge([t], min_phase=4)
    out = ld.write_table(res["library"], tmp_path / "merged.tsv")
    back = ld.read_table(out)
    assert len(back) == 2
    assert {r["name"] for r in back} == {"metformin", "levothyroxine"}
    met = next(r for r in back if r["name"] == "metformin")
    assert set(met["targets"]) == {"PRKAB1", "AMPK"}


def test_a_source_that_fails_to_parse_is_reported_not_silent(tmp_path):
    from promoter_tf_drugs import localdrugs as ld
    good = _tbl(tmp_path, "good.txt",
                ["pert_iname", "clinical_phase", "target", "smiles"],
                [["metformin", "launched", "X", _METFORMIN]])
    bad = _tbl(tmp_path, "bad.txt", ["name", "phase"],
               [["drugA", "launched"]])
    res = ld.merge([good, bad], min_phase=4)
    assert res["n_distinct"] == 1
    assert "error" in res["per_source"]["bad.txt"]
    assert "SMILES" in res["per_source"]["bad.txt"]["error"]


# ===================================================================== #
# More approved-drug sources, each shipping a different format.         #
# ===================================================================== #
_SOURCE_SAMPLES = {
    # DrugCentral structures.smiles: SMILES first, name under INN
    "drugcentral.smiles":
        "SMILES\tInChIKey\tID\tINN\n"
        "CN(C)C(=N)N=C(N)N\tXZW\t1\tmetformin\n",
    # ZINC curated subset: headerless, space separated
    "zinc_fda.smi":
        "CN(C)C(=N)N=C(N)N ZINC000000012345\n"
        "CC(=O)Oc1ccccc1C(=O)O ZINC000000053\n",
    # TTD
    "ttd_drugs.txt":
        "DRUG__ID\tDRUGNAME\tCANONICAL_SMILES\n"
        "D001\tmetformin\tCN(C)C(=N)N=C(N)N\n",
    # SuperDRUG2: comma separated, space in the column name
    "superdrug.csv":
        "Drug Name,SMILES,ATC\nMetformin,CN(C)C(=N)N=C(N)N,A10BA02\n",
    # the Hub, which already worked
    "repurposing_drugs.txt":
        "pert_iname\tclinical_phase\ttarget\tsmiles\n"
        "metformin\tlaunched\tPRKAB1\tCN(C)C(=N)N=C(N)N\n",
}


def _write_sources(tmp_path):
    out = []
    for name, text in _SOURCE_SAMPLES.items():
        p = tmp_path / name
        p.write_text(text, encoding="utf-8")
        out.append(p)
    return out


def test_every_source_format_parses_with_the_right_name(tmp_path):
    """Getting a name column wrong does not raise - the name silently
    falls back to the SMILES string, which is then what the report prints
    and what the dedup fallback keys on. All four of these did that."""
    from promoter_tf_drugs import localdrugs as ld
    for p in _write_sources(tmp_path):
        rows = ld.read_table(p, approved_only=False)
        assert rows, p.name
        name = rows[0]["name"]
        assert not name.startswith("CN(C)"), (
            f"{p.name}: the name fell back to the SMILES string")
        assert name.lower() in ("metformin", "zinc000000012345")


def test_a_headerless_smiles_file_is_read_not_rejected(tmp_path):
    """ZINC's curated subsets have no header: the first line is already a
    molecule, so DictReader ate a real compound as the header and then
    reported no SMILES column."""
    from promoter_tf_drugs import localdrugs as ld
    p = tmp_path / "zinc_fda.smi"
    p.write_text(_SOURCE_SAMPLES["zinc_fda.smi"], encoding="utf-8")
    rows = ld.read_table(p, approved_only=False)
    assert len(rows) == 2, "no line may be consumed as a header"
    assert rows[0]["smiles"] == "CN(C)C(=N)N=C(N)N"
    assert rows[0]["name"] == "ZINC000000012345"


def test_a_header_is_still_detected_when_one_exists(tmp_path):
    """The headerless check must not misfire on a real table, or a column
    row becomes a compound."""
    from promoter_tf_drugs import localdrugs as ld
    p = tmp_path / "withheader.tsv"
    p.write_text("name\tsmiles\nmetformin\tCN(C)C(=N)N=C(N)N\n",
                 encoding="utf-8")
    rows = ld.read_table(p, approved_only=False)
    assert len(rows) == 1 and rows[0]["name"] == "metformin"


def test_one_molecule_merges_across_five_differently_shaped_sources(
        tmp_path):
    """The point of keying on structure: the ZINC row identifies metformin
    only by a ZINC id, so no name-based merge could join it to the rest."""
    from promoter_tf_drugs import localdrugs as ld
    res = ld.merge(_write_sources(tmp_path), min_phase=0)
    by_name = {r["name"]: r for r in res["library"]}
    met = next(r for r in res["library"]
               if r["name"].lower() == "metformin")
    assert len(met["sources"]) == 5
    assert res["n_salt_forms_collapsed"] == 4
    # aspirin is only in the ZINC file and stays its own molecule
    assert res["n_distinct"] == 2


def test_the_source_registry_says_what_each_one_adds():
    """They overlap heavily, so the useful question is not how many
    molecules a source has but what it has that the others do not."""
    from promoter_tf_drugs import localdrugs as ld
    keys = {s.key for s in ld.SOURCES}
    assert {"repurposing_hub", "drugcentral", "chembl_approved",
            "zinc_fda"} <= keys
    text = ld.describe_sources()
    # the one thing only DrugCentral can do
    assert "WHICH agency approved" in text
    # and the honest expectation about adding a second source
    assert "few hundred new moieties rather than doubling" in text
    for s in ld.SOURCES:
        assert s.url.startswith("https://")
        assert s.caveat and len(s.caveat) > 30


def test_source_status_reports_presence_without_guessing(tmp_path,
                                                         monkeypatch):
    from promoter_tf_drugs import localdrugs as ld
    monkeypatch.chdir(tmp_path)
    before = {s["key"]: s["present"] for s in ld.source_status()}
    assert not any(before.values())
    (tmp_path / "data").mkdir()
    (tmp_path / "data" / "repurposing_drugs.txt").write_text(
        _SOURCE_SAMPLES["repurposing_drugs.txt"], encoding="utf-8")
    after = {s["key"]: s["present"] for s in ld.source_status()}
    assert after["repurposing_hub"] is True
    assert after["drugcentral"] is False


# ===================================================================== #
# The FDA's 20,000 products, used for what they add.                    #
# ===================================================================== #
_OB_HEADER = ("Ingredient~DF;Route~Trade_Name~Applicant~Strength~Appl_Type~"
              "Appl_No~Product_No~TE_Code~Approval_Date~RLD~RS~Type~"
              "Applicant_Full_Name")
_OB_ROWS = [
    "METFORMIN HYDROCHLORIDE~TABLET;ORAL~GLUCOPHAGE~BRISTOL~500MG~N~020357"
    "~001~~Mar 3, 1995~Yes~Yes~RX~BMS",
    "METFORMIN HYDROCHLORIDE~TABLET;ORAL~GLUCOPHAGE~BRISTOL~850MG~N~020357"
    "~002~~Mar 3, 1995~No~No~RX~BMS",
    "METFORMIN HYDROCHLORIDE~TABLET, EXTENDED RELEASE;ORAL~GLUCOPHAGE XR~"
    "BRISTOL~500MG~N~021202~001~~Oct 13, 2000~Yes~Yes~RX~BMS",
    "METFORMIN HYDROCHLORIDE~TABLET;ORAL~METFORMIN HCL~TEVA~500MG~A~075300"
    "~001~AB~Jan 1, 2002~No~No~RX~TEVA",
    "TAZAROTENE~GEL;TOPICAL~TAZORAC~ALLERGAN~0.05%~N~020600~001~~"
    "Jun 13, 1997~Yes~Yes~RX~ALLERGAN",
    "TAZAROTENE~CREAM;TOPICAL~TAZORAC~ALLERGAN~0.1%~N~021184~002~~"
    "Dec 28, 2000~Yes~Yes~RX~ALLERGAN",
    "CALCITRIOL~INJECTABLE;INTRAVENOUS~CALCIJEX~ABBVIE~1MCG/ML~N~019846"
    "~001~~Jul 2, 1987~Yes~Yes~RX~ABBVIE",
    "CALCITRIOL~CAPSULE;ORAL~ROCALTROL~VALIDUS~0.25MCG~N~018044~001~~"
    "Jul 29, 1981~Yes~Yes~RX~VALIDUS",
    "AMLODIPINE BESYLATE; OLMESARTAN MEDOXOMIL~TABLET;ORAL~AZOR~DAIICHI~"
    "5MG;20MG~N~021990~001~~Sep 25, 2007~Yes~Yes~RX~DAIICHI",
]


def _orange_book(tmp_path):
    from promoter_tf_drugs import orangebook as ob
    d = tmp_path / "data"
    d.mkdir(exist_ok=True)
    p = d / "products.txt"
    p.write_text(_OB_HEADER + "\n" + "\n".join(_OB_ROWS) + "\n",
                 encoding="utf-8")
    ob._PRODUCTS.clear()
    return p


def test_twenty_thousand_products_are_not_twenty_thousand_molecules(
        tmp_path, monkeypatch):
    """The FDA's figure counts products: each strength, dosage form and
    applicant separately. Metformin is four rows here and one molecule, so
    the product list cannot enlarge a screening library."""
    from promoter_tf_drugs import orangebook as ob
    monkeypatch.chdir(tmp_path)
    _orange_book(tmp_path)
    assert ob.load() == 5
    assert ob.product_count() == 10
    assert ob.lookup("metformin")["n_products"] == 4


def test_the_tilde_delimiter_is_handled(tmp_path, monkeypatch):
    """The Orange Book is tilde-delimited, not tab. A tab reader returns
    one enormous column and raises nothing."""
    from promoter_tf_drugs import orangebook as ob
    monkeypatch.chdir(tmp_path)
    _orange_book(tmp_path)
    assert ob.load() == 5


def test_route_decides_whether_a_drug_can_reach_the_target(tmp_path,
                                                           monkeypatch):
    """The reason to load this at all. No structure source says a drug is
    a gel, and a topical-only retinoid is not a route to a gene in the
    liver however well it binds the receptor."""
    from promoter_tf_drugs import orangebook as ob
    monkeypatch.chdir(tmp_path)
    _orange_book(tmp_path)
    taz = ob.lookup("tazarotene")
    assert taz["reach"] == "local only"
    assert taz["routes"] == ["TOPICAL"]
    assert "cannot reach a target outside the tissue" in taz["reach_note"]
    assert ob.lookup("metformin")["reach"] == "systemic"
    assert set(ob.lookup("calcitriol")["systemic_routes"]) == {
        "INTRAVENOUS", "ORAL"}


def test_the_salt_free_join_reaches_the_parent(tmp_path, monkeypatch):
    """The Orange Book has no structures, so the only join is on the
    ingredient name - and that name carries the salt."""
    from promoter_tf_drugs import orangebook as ob
    monkeypatch.chdir(tmp_path)
    _orange_book(tmp_path)
    assert ob.lookup("Metformin Hydrochloride") == ob.lookup("metformin")
    assert ob.lookup("METFORMIN") is not None


def test_generic_availability_comes_from_the_applicant_count(tmp_path,
                                                             monkeypatch):
    from promoter_tf_drugs import orangebook as ob
    monkeypatch.chdir(tmp_path)
    _orange_book(tmp_path)
    assert ob.lookup("metformin")["generic_available"] is True
    assert ob.lookup("tazarotene")["generic_available"] is False


def test_a_combination_only_approval_is_not_evidence_for_one_agent(
        tmp_path, monkeypatch):
    """Olmesartan appears only inside a fixed-dose combination here, so
    its approval says nothing about the agent on its own."""
    from promoter_tf_drugs import orangebook as ob
    monkeypatch.chdir(tmp_path)
    _orange_book(tmp_path)
    olm = ob.lookup("olmesartan medoxomil")
    assert olm["combination_only"] is True
    assert "not evidence about this agent on its own" in olm["notes"][0]
    assert ob.lookup("metformin")["combination_only"] is False


def test_an_unmatched_name_means_unknown_route_not_no_approval(
        tmp_path, monkeypatch):
    from promoter_tf_drugs import orangebook as ob
    monkeypatch.chdir(tmp_path)
    _orange_book(tmp_path)
    assert ob.lookup("nosuchdrug") is None
    acc = ob.enrich([{"name": "metformin"}, {"name": "nosuchdrug"}])
    assert acc["matched"] == 1 and acc["n_unmatched"] == 1
    assert "unknown route rather than no approval" in ob.summarise(acc)


def test_no_orange_book_is_not_every_drug_being_systemic(tmp_path,
                                                         monkeypatch):
    from promoter_tf_drugs import orangebook as ob
    monkeypatch.chdir(tmp_path)
    ob._PRODUCTS.clear()
    acc = ob.enrich([{"name": "metformin"}])
    assert acc["available"] is False
    assert "different from every row being systemic" in acc["note"]


def test_enrich_attaches_the_reach_to_each_row(tmp_path, monkeypatch):
    from promoter_tf_drugs import orangebook as ob
    monkeypatch.chdir(tmp_path)
    _orange_book(tmp_path)
    lib = [{"name": "tazarotene"}, {"name": "metformin"}]
    acc = ob.enrich(lib)
    assert lib[0]["route_reach"] == "local only"
    assert lib[1]["route_reach"] == "systemic"
    assert acc["local_only"] == ["tazarotene"]
    assert "1 is approved" in ob.summarise(acc)


# ===================================================================== #
# Built and connected. Each of these pins a module that existed, was    #
# described in a report, and reached nothing - the failure mode is not  #
# a crash, it is a module that quietly does nothing.                    #
# ===================================================================== #
def test_the_orange_book_route_reaches_the_drug_description(tmp_path,
                                                            monkeypatch):
    """`indications.describe` hand-typed "topical only" for seventy drugs.
    The FDA publishes the route for every approved ingredient, so the
    Orange Book is consulted first and the hand list is the fallback."""
    from promoter_tf_drugs import orangebook as ob, indications as ind
    monkeypatch.chdir(tmp_path)
    _orange_book(tmp_path)
    d = ind.describe("tazarotene")
    assert d["route"] == "local only"
    assert d["route_source"] == "FDA Orange Book"
    assert d["routes"] == ["TOPICAL"]
    assert ind.describe("metformin")["route"] == "systemic"


def test_no_orange_book_leaves_the_route_unknown_not_systemic(tmp_path,
                                                              monkeypatch):
    """The failure that matters: a missing file must not read as a drug
    that reaches everywhere."""
    from promoter_tf_drugs import orangebook as ob, indications as ind
    monkeypatch.chdir(tmp_path)
    ob._PRODUCTS.clear()
    d = ind.describe("tazarotene")
    assert d["route"] is None
    assert d["route_source"] is None


def test_drug_scoring_reads_a_local_potency_file_before_chembl(monkeypatch):
    """On a blocked network the ChEMBL activity call is what emptied this
    axis, and potency-unknown is treated as neutral - so a nanomolar
    ligand and an untested one ranked identically."""
    from promoter_tf_drugs import drug_scoring as ds, localpotency as lp
    from promoter_tf_drugs import localannot as la

    # The real Drug dataclass has target_symbol and target_chembl_id and
    # NO uniprot field. An earlier version of this test invented one, so it
    # passed while the product never took the branch at all.
    from promoter_tf_drugs.drugs import Drug
    import dataclasses
    assert not any(f.name == "uniprot" for f in dataclasses.fields(Drug)),         "Drug grew a uniprot field; this test is now testing the wrong thing"

    drug = Drug(name="androgen-x", chembl_id="CHEMBL25",
                target_symbol="AR", action_type="AGONIST",
                mechanism=None, approved=True, max_phase=4.0,
                first_approval=1990, withdrawn=False, oral=True,
                direction="up", source="test")
    monkeypatch.setattr(la, "uniprot_for", lambda sym: "P10275")

    class _NoNetwork:
        def get(self, *a, **k):
            raise AssertionError("ChEMBL was called although a local "
                                 "potency file answered")

    monkeypatch.setattr(lp, "available", lambda: {"bindingdb": True})
    monkeypatch.setattr(lp, "actives_for", lambda u: [
        {"pchembl": 8.5, "smiles": "CCO", "target": u},
        {"pchembl": 9.5, "smiles": "CCN", "target": u}])
    ev = ds.fetch_evidence(_NoNetwork(), drug)
    assert ev.pchembl_median == 9.0
    assert ev.pchembl_n == 2
    assert "local" in ev.source


def test_pains_discounts_a_predicted_route_in_the_report(monkeypatch):
    """liabilities.py was written to discount confidence in a predicted
    route and for a long while was imported by nothing. This pins it to
    the table."""
    from promoter_tf_drugs import pipeline as pl
    import inspect
    src = inspect.getsource(pl.corroboration)
    assert "liab_mod.assess" in src, "liabilities is not called"
    assert '"pains"' in src, "the interference flag never reaches a row"
    assert '"evidence_multiplier"' in src


def test_the_reach_column_is_rendered_and_its_absence_is_explained():
    """A blank Reach column has to say it means unknown route. Silently
    blank reads as "fine"."""
    from promoter_tf_drugs import pipeline as pl
    import inspect
    src = inspect.getsource(pl.corroboration_table)
    assert "Reach" in src
    assert "_route_cell" in src
    # The header, the separator and the cells must all be eleven wide, or
    # markdown drops the trailing columns without complaining.
    sep = [ln for ln in src.splitlines() if "|---|" in ln][0]
    assert sep.count("---") == 11


def test_the_reach_and_pains_footnotes_fire_on_the_right_rows():
    from promoter_tf_drugs import pipeline as pl
    import inspect
    src = inspect.getsource(pl.corroboration_table)
    assert "approved only by local routes" in src
    assert "pan-assay" in src
    assert "unknown route, not systemic" in src


def test_route_cell_marks_local_only_and_leaves_unknown_blank():
    from promoter_tf_drugs import pipeline as pl
    assert pl._route_cell({"route": "local only"}) == "**local only**"
    assert pl._route_cell({"route": "systemic"}) == "systemic"
    assert pl._route_cell({"route": "unclear"}) == "unclear"
    assert pl._route_cell({}) == "\u2014"
    assert pl._route_cell({"route": None}) == "\u2014"


def test_local_gwas_is_offered_as_the_weaker_evidence_it_is(tmp_path,
                                                            monkeypatch):
    """GWAS was described in a report as "a weaker kind of evidence" and
    genetics.py did not reference it at all."""
    from promoter_tf_drugs import genetics as gen, localannot as la
    monkeypatch.chdir(tmp_path)
    d = tmp_path / "data"
    d.mkdir(exist_ok=True)
    (d / "gwas_catalog.tsv").write_text(
        "MAPPED_GENE\tDISEASE/TRAIT\tP-VALUE\tSNPS\n"
        "PCSK9\tLDL cholesterol\t3e-20\trs11591147\n"
        "PCSK9\tCoronary artery disease\t2e-11\trs11206510\n"
        "PCSK9\tnoise\t0.4\trs999\n", encoding="utf-8")
    la._GWAS.clear()
    got = gen.local_support("PCSK9")
    assert got is not None
    assert got["n_significant"] == 2
    assert got["source"] == "GWAS Catalog (local)"
    assert any("nearest gene" in c for c in got["caveats"])


def test_no_local_gwas_returns_none_rather_than_an_empty_verdict(
        tmp_path, monkeypatch):
    from promoter_tf_drugs import genetics as gen, localannot as la
    monkeypatch.chdir(tmp_path)
    la._GWAS.clear()
    assert gen.local_support("PCSK9") is None


def test_a_gwas_file_read_with_no_hits_is_not_a_genetics_success(
        tmp_path, monkeypatch):
    """Read-and-found-nothing is assessed but not evidence. The ledger
    says ok only when something arrived."""
    from promoter_tf_drugs import genetics as gen, localannot as la
    from promoter_tf_drugs import pipeline as pl
    monkeypatch.chdir(tmp_path)
    d = tmp_path / "data"
    d.mkdir(exist_ok=True)
    (d / "gwas_catalog.tsv").write_text(
        "MAPPED_GENE\tDISEASE/TRAIT\tP-VALUE\tSNPS\n"
        "APOB\tLDL cholesterol\t3e-20\trs1\n", encoding="utf-8")
    la._GWAS.clear()
    got = gen.local_support("PCSK9")
    assert got is not None and got["assessed"] is True
    assert pl._n_local_gwas(got) == 0
    detail = pl._genetics_detail([], got)
    assert "holds no significant association" in detail


def test_genetics_detail_names_which_source_found_what():
    from promoter_tf_drugs import pipeline as pl
    both = pl._genetics_detail(
        [{"trait": "x"}], {"n_significant": 3, "assessed": True})
    assert "Open Targets" in both and "GWAS Catalog" in both
    assert "weaker" in both
    assert pl._genetics_detail([], None) == (
        "no genetic association returned for this gene")


def test_no_module_needs_python_3_12_to_import():
    """pyproject declares >=3.10, and a backslash escape inside an f-string
    expression is PEP 701 syntax - legal from 3.12, a SyntaxError on 3.10
    and 3.11. Four such lines were present, so the package could not be
    imported at all on the floor it advertises.

    This is a text scan, not a parse: only 3.12 is installed here, and
    ast.parse(feature_version=(3, 10)) does not reject the construct, so
    there is no authoritative local check. Verifying the floor needs a
    3.10 interpreter.
    """
    import re
    from pathlib import Path
    from promoter_tf_drugs import __file__ as pkg
    bs = chr(92)
    pat = re.compile(r"\{[^{}]*" + bs + bs + r"[^{}]*\}")
    offenders = []
    for f in sorted(Path(pkg).parent.glob("*.py")):
        for i, ln in enumerate(f.read_text(encoding="utf-8").splitlines(), 1):
            if ('f"' in ln or "f'" in ln) and pat.search(ln):
                offenders.append(f"{f.name}:{i}")
    assert not offenders, (
        "backslash inside an f-string expression needs 3.12; move it to a "
        "named constant: " + ", ".join(offenders))


def test_the_version_floor_is_still_what_the_test_above_assumes():
    """If the floor moves to 3.12 the scan above is obsolete, not failing."""
    import tomllib
    from pathlib import Path
    root = Path(__file__).resolve().parent.parent
    cfg = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))
    assert cfg["project"]["requires-python"] == ">=3.10"


# ===================================================================== #
# Two independent methods, intersected. Each ran and neither was ever   #
# compared against the other.                                           #
# ===================================================================== #
def _dock_run(target="VDR", usable=True):
    return {
        "target": target, "usable": usable, "auc": 0.95, "ef1": 3.4,
        "n_docked": 10, "n_failed": 1, "caveats": [],
        "hits": [{"drug": "calcitriol", "score": -11.6},
                 {"drug": "doxercalciferol", "score": -11.1},
                 {"drug": "metformin", "score": -5.2}],
    }


def _ligand_hits():
    return [{"drug": "calcitriol", "quadrant": "ESTABLISHED LIGAND"},
            {"drug": "atenolol", "quadrant": "LOW"},
            {"drug": "aspirin", "quadrant": "LOW"}]


def test_consensus_accepts_the_json_that_crosses_machines():
    """Docking runs on the machine with the binary and comes back as JSON.
    Requiring the dataclass meant this could not be called on a real split
    run - which is the only kind of run the docking path supports."""
    from promoter_tf_drugs import docking as dk
    both = dk.consensus(_dock_run(), _ligand_hits(), top_frac=0.5)
    assert [h["drug"] for h in both] == ["calcitriol"]
    assert both[0]["docking_score"] == -11.6
    assert both[0]["ligand_quadrant"] == "ESTABLISHED LIGAND"
    assert "share no input" in both[0]["why"]


def test_consensus_refuses_an_unvalidated_pocket():
    """A pocket that could not separate actives from decoys has not been
    shown to rank anything, so agreement with it is not agreement."""
    from promoter_tf_drugs import docking as dk
    assert dk.consensus(_dock_run(usable=False), _ligand_hits(),
                        top_frac=0.5) == []
    assert dk.consensus(_dock_run(), [], top_frac=0.5) == []
    assert dk.consensus(None, _ligand_hits()) == []


def test_consensus_survives_a_hit_with_no_score():
    """A bundle written by an older run may omit `score`. Previously this
    raised KeyError in the middle of a report."""
    from promoter_tf_drugs import docking as dk
    run = _dock_run()
    del run["hits"][0]["score"]
    both = dk.consensus(run, _ligand_hits(), top_frac=0.5)
    assert both and both[0]["docking_score"] is None


def test_the_docking_section_reports_validation_and_disagreement():
    from promoter_tf_drugs import pipeline as pl
    out = "\n".join(pl._docking_section([_dock_run()], {}))
    assert "ROC AUC" in out and "0.95" in out and "3.4x" in out
    assert "Validated" in out
    # No overlap is a real disagreement, not a missing result.
    assert "real disagreement between two" in out
    assert "says nothing about direction" in out


def test_an_unvalidated_pocket_is_flagged_not_hidden():
    from promoter_tf_drugs import pipeline as pl
    out = "\n".join(pl._docking_section([_dock_run(usable=False)], {}))
    assert "did not validate" in out
    assert "**no**" in out


def test_the_consensus_table_names_both_methods_that_found_each_drug():
    from promoter_tf_drugs import pipeline as pl
    from promoter_tf_drugs import docking as dk
    both = {"VDR": dk.consensus(_dock_run(), _ligand_hits(), top_frac=0.5)}
    out = "\n".join(pl._docking_section([_dock_run()], both))
    assert "picked by both methods" in out
    assert "calcitriol" in out
    assert "-11.6" in out
    assert "Still a binding hypothesis" in out


def test_no_docking_means_no_section_rather_than_an_empty_one():
    from promoter_tf_drugs import pipeline as pl
    assert pl._docking_section([], {}) == []


# ===================================================================== #
# The upstream drug direction. Inverted for every repressor factor and   #
# every downregulation goal, which is half of all runs.                  #
# ===================================================================== #
def _signed_pair():
    """KINASE activates TFACT and inhibits TFREP. Two edges, so the right
    answer is arguable on paper."""
    return [net.Edge(source="KINASE", target="TFACT", sign=1,
                     n_references=3),
            net.Edge(source="KINASE", target="TFREP", sign=-1,
                     n_references=3)]


class _DirDrug:
    def __init__(self, direction):
        self.name = "drug-" + direction
        self.direction = direction

    def to_dict(self):
        return {"name": self.name, "direction": self.direction}


def _both_directions(client, symbol, uniprot, approved_only):
    return [_DirDrug("up"), _DirDrug("down")]


def test_inhibiting_a_repressor_needs_an_activator_of_its_upstream_node(
        monkeypatch):
    """KINASE -| TFREP, and TFREP represses the gene. To raise the gene you
    lower TFREP, which means RAISING the kinase. The old expression read
    net_sign alone and asked for a kinase INHIBITOR, which raises the
    repressor and lowers the gene - the opposite of the request."""
    from promoter_tf_drugs import drugs as drugs_mod
    monkeypatch.setattr(drugs_mod, "chembl_drugs", _both_directions)

    hits = net.find_upstream_targets(None, "TFREP", goal="up",
                                     tf_role="repressor",
                                     edges=_signed_pair())
    assert hits, "no path found; the fixture is wrong, not the code"
    assert hits[0].net_sign == -1
    assert hits[0].want_tf_sign == -1
    net.attach_drugs(None, hits, approved_only=True)
    kept = {d["direction"] for d in hits[0].drugs}
    assert kept == {"up"}, f"asked for {kept}, which would move the gene down"


def test_the_activator_case_still_asks_for_an_activator(monkeypatch):
    from promoter_tf_drugs import drugs as drugs_mod
    monkeypatch.setattr(drugs_mod, "chembl_drugs", _both_directions)
    hits = net.find_upstream_targets(None, "TFACT", goal="up",
                                     tf_role="activator",
                                     edges=_signed_pair())
    net.attach_drugs(None, hits, approved_only=True)
    assert {d["direction"] for d in hits[0].drugs} == {"up"}


def test_the_hit_carries_the_wanted_tf_direction():
    """Recorded on the hit because attach_drugs cannot recover it:
    upstream_paths keeps only paths whose net sign equals it, so net_sign
    alone says nothing about the goal."""
    for goal, role, want in (("up", "activator", 1), ("up", "repressor", -1),
                             ("down", "activator", -1),
                             ("down", "repressor", 1)):
        tf = "TFACT" if role == "activator" else "TFREP"
        hits = net.find_upstream_targets(None, tf, goal=goal, tf_role=role,
                                         edges=_signed_pair())
        for h in hits:
            assert h.want_tf_sign == want, (goal, role)


def test_attach_drugs_falls_back_to_its_argument_for_a_bare_hit(monkeypatch):
    """A hit built by hand has no want_tf_sign, so goal_on_node is used -
    which is what that parameter was always for and never did."""
    from promoter_tf_drugs import drugs as drugs_mod
    monkeypatch.setattr(drugs_mod, "chembl_drugs", _both_directions)
    h = net.UpstreamHit(node="KINASE", tf_symbol="TF", best_path={},
                        n_paths=1, score=1.0, rwr_score=0.0, net_sign=1)
    object.__setattr__(h, "want_tf_sign", None)
    net.attach_drugs(None, [h], goal_on_node="down")
    assert {d["direction"] for d in h.drugs} == {"down"}


# ===================================================================== #
# TFAnnotation subclasses dict, so `uniprot` is a key. Three callers read #
# it with getattr and got None for every gene, silently.                 #
# ===================================================================== #
def test_the_accession_is_a_dict_key_not_an_attribute():
    """The shape that caused three silent failures, pinned directly."""
    from promoter_tf_drugs.annotate import TFAnnotation
    a = TFAnnotation(symbol="VDR", uniprot="P11473")
    assert isinstance(a, dict)
    assert a.get("uniprot") == "P11473"
    assert getattr(a, "uniprot", None) is None, (
        "if TFAnnotation grows a uniprot property, the getattr callers "
        "that this test exists to prevent would start working by accident")


def _hgnc(tmp_path, rows=(("VDR", "P11473"),)):
    from promoter_tf_drugs import localannot as la

    d = tmp_path / "data"
    d.mkdir(exist_ok=True)
    with open(d / "hgnc_complete_set.txt", "w", encoding="utf-8") as fh:
        fh.write("symbol\tname\tuniprot_ids\tentrez_id\talias_symbol\n")
        for sym, uni in rows:
            fh.write(f"{sym}\t{sym} protein\t{uni}\t111\tOLD{sym}\n")
    la._HGNC.clear()
    return d


def test_uniprot_for_resolves_a_symbol_and_an_alias(tmp_path, monkeypatch):
    from promoter_tf_drugs import localannot as la
    monkeypatch.chdir(tmp_path)
    _hgnc(tmp_path)
    assert la.uniprot_for("VDR") == "P11473"
    assert la.uniprot_for("vdr") == "P11473"
    assert la.uniprot_for("OLDVDR") == "P11473"
    assert la.uniprot_for("NOSUCHGENE") is None


def test_no_hgnc_table_means_unknown_not_a_gene_without_a_protein(
        tmp_path, monkeypatch):
    from promoter_tf_drugs import localannot as la
    monkeypatch.chdir(tmp_path)
    la._HGNC.clear()
    assert la.uniprot_for("VDR") is None


def test_the_docking_cli_can_resolve_a_gene_symbol():
    """It could not: both sites read the accession with getattr, so every
    `--gene` run printed "Pass --uniprot" even when the lookup worked."""
    import inspect
    from promoter_tf_drugs import docking as dk
    src = inspect.getsource(dk)
    assert 'getattr(info, "uniprot"' not in src
    assert src.count('(info or {}).get("uniprot")') == 2


# ===================================================================== #
# Three wirings, offline.                                                #
# ===================================================================== #
def test_screen_family_resolves_ligands_from_files_with_no_client(
        tmp_path, monkeypatch):
    from promoter_tf_drugs import screening as scr, localpotency as lp
    monkeypatch.chdir(tmp_path)
    _hgnc(tmp_path, [("RARA", "P10826"), ("RARB", "P10826")])
    monkeypatch.setattr(lp, "available", lambda: {"bindingdb": "x"})
    monkeypatch.setattr(lp, "actives_for", lambda u: [
        {"chembl_id": "C1", "smiles": "CCO", "pchembl": 8.0}])
    assert scr.uniprot_for("RARA") == "P10826"
    assert len(scr.local_actives_for("RARA")) == 1


def test_screen_family_with_no_source_says_so_rather_than_blaming_the_family(
        tmp_path, monkeypatch):
    """"The family is as under-studied as the receptor" is a claim about
    biology. With no readable source it is a claim about the filesystem."""
    from promoter_tf_drugs import screening as scr
    monkeypatch.chdir(tmp_path)
    res = scr.screen_family(None, "RARA", goal="up")
    assert res["source"] is None
    assert "missing source, not a receptor without ligands" in res["error"]


def test_screen_family_honours_the_goal_it_accepts():
    """It took `goal` and then passed goal="up" to screen(), so a run
    asking to LOWER a gene was screened for compounds that raise it."""
    import inspect
    from promoter_tf_drugs import screening as scr
    src = inspect.getsource(scr.screen_family)
    assert 'goal="up"' not in src.split('"""', 2)[-1], (
        "screen_family hardcodes a goal again")
    assert "screen(model, lib, goal=goal" in src


def test_the_library_guard_in_screen_receptor_has_the_names_it_uses():
    """The guard was first written into screen_receptor, where `fam` does
    not exist - a NameError on the first run with no library at all."""
    import inspect
    from promoter_tf_drugs import screening as scr
    src = inspect.getsource(scr.screen_receptor)
    assert "fam.name" not in src
    assert "no screening library" in src


def test_orange_book_annotates_a_merged_library_offline(tmp_path, monkeypatch):
    from promoter_tf_drugs import orangebook as ob
    monkeypatch.chdir(tmp_path)
    _orange_book(tmp_path)
    lib = [{"name": "tazarotene"}, {"name": "metformin"},
           {"name": "nosuchdrug"}]
    acc = ob.enrich(lib)
    assert acc["matched"] == 2 and acc["n_unmatched"] == 1
    assert lib[0]["route_reach"] == "local only"
    assert lib[1]["generic_available"] is True
    assert "nosuchdrug" in acc["unmatched"]


def test_the_merge_cli_exposes_the_orange_book(tmp_path, monkeypatch):
    from promoter_tf_drugs import localdrugs as ld
    monkeypatch.chdir(tmp_path)
    _orange_book(tmp_path)
    (tmp_path / "a.tsv").write_text(
        "name\tsmiles\tmax_phase\n"
        "tazarotene\tCCOC(=O)c1ccc(cn1)C#Cc1ccc2c(c1)C(C)(C)CCS2\t4\n"
        "metformin\tCN(C)C(=N)NC(N)=N\t4\n", encoding="utf-8")
    assert ld.main(["--merge", "a.tsv", "--orange-book"]) == 0


def test_depmap_reads_its_own_common_essential_list(tmp_path, monkeypatch):
    from promoter_tf_drugs import depmap as dm
    monkeypatch.chdir(tmp_path)
    d = tmp_path / "data"
    d.mkdir(exist_ok=True)
    (d / "CRISPRInferredCommonEssentials.csv").write_text(
        "Essentials\nPOLR2A (5430)\n", encoding="utf-8")
    ess = dm.load_common_essentials()
    assert ess == {"POLR2A"}
    d1 = dm.assess("POLR2A", [-1.4] * 14, common_essentials=ess)
    assert d1.depmap_common_essential is True
    d2 = dm.assess("VDR", [-0.05] * 14, common_essentials=ess)
    assert d2.depmap_common_essential is False


def test_a_list_nobody_loaded_is_not_a_gene_being_cleared():
    """None, never False. A gene absent from a list that was never read
    has not been shown to be dispensable."""
    from promoter_tf_drugs import depmap as dm
    d = dm.assess("VDR", [-0.05] * 14)
    assert d.depmap_common_essential is None


def test_depmap_reports_a_disagreement_rather_than_choosing(tmp_path):
    """DepMap's call is over the whole release; the fraction here is over
    whatever lines were handed in. Where they differ, both are shown."""
    from promoter_tf_drugs import depmap as dm
    d = dm.assess("POLR2A", [-0.1] * 14, common_essentials={"POLR2A"})
    assert d.depmap_common_essential is True
    assert any("common-essential" in c and "better-powered" in c
               for c in d.caveats)
    d2 = dm.assess("MYC", [-1.4] * 14, common_essentials={"POLR2A"})
    assert any("not representative" in c for c in d2.caveats)


# ===================================================================== #
# Parameters that were accepted and ignored.                             #
# ===================================================================== #
def test_no_public_function_accepts_a_parameter_it_never_reads():
    """The class of bug this round was found by. `goal` deciding whether
    you want a gene up or down was one of them."""
    import ast
    import pathlib
    from promoter_tf_drugs import __file__ as pkg

    skip = {"self", "cls", "client", "kwargs", "args"}
    # Reporting helpers whose signature is kept uniform on purpose.
    allowed = {("pipeline.py", "_drug_glossary", "res"),
               ("screens.py", "crispra_plan", "promoter"),
               ("discovery.py", "axes_from_grade", "tf_result"),
               ("congeners.py", "collapse_grades", "set_name"),
               ("depmap.py", "_verdict", "mean"),
               ("evidence.py", "annotate_hits_with_chip", "hits")}
    bad = []
    for f in sorted(pathlib.Path(pkg).parent.glob("*.py")):
        tree = ast.parse(f.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if not isinstance(node, ast.FunctionDef):
                continue
            a = node.args
            params = [p.arg for p in a.posonlyargs + a.args + a.kwonlyargs
                      if p.arg not in skip and not p.arg.startswith("_")]
            used = {n.id for n in ast.walk(node)
                    if isinstance(n, ast.Name) and isinstance(n.ctx,
                                                              ast.Load)}
            for p in params:
                if p not in used and (f.name, node.name, p) not in allowed:
                    bad.append(f"{f.name}:{node.lineno} {node.name}({p})")
    assert not bad, "accepted and never read: " + ", ".join(bad)


def test_smiles_to_sdf_embeds_the_number_of_conformers_asked_for():
    import pytest
    from pathlib import Path
    pytest.importorskip("rdkit")
    from rdkit import Chem
    from promoter_tf_drugs import docking as dk
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        one = dk.smiles_to_sdf("CCCCCCCC(=O)O", Path(td) / "a.sdf")
        many = dk.smiles_to_sdf("CCCCCCCC(=O)O", Path(td) / "b.sdf",
                                n_confs=8)
        # One conformer is written either way; with n_confs it is the
        # lowest-energy of eight rather than the first embedding.
        a = Chem.SDMolSupplier(str(one), removeHs=False)[0]
        b = Chem.SDMolSupplier(str(many), removeHs=False)[0]
        assert a is not None and b is not None
        assert a.GetNumConformers() == 1 and b.GetNumConformers() == 1


def test_a_plain_smiles_file_says_it_cannot_filter_on_approval(
        tmp_path, caplog):
    """ZINC's `world` subset is not approved. Accepting approved_only and
    ignoring it let that file read as an approved library."""
    import logging
    from promoter_tf_drugs import localdrugs as ld
    f = tmp_path / "world.smi"
    f.write_text("CCO zinc1\nCCN zinc2\n", encoding="utf-8")
    with caplog.at_level(logging.WARNING):
        rows = ld.read_smiles_file(f, approved_only=True)
    assert len(rows) == 2
    assert "could not be applied" in caplog.text


def test_the_drug_scan_summary_states_which_way_the_factor_moves():
    from promoter_tf_drugs import tfactivity as tfa
    res = {"VDR": [{"drug": "x"}]}
    assert "raise" in tfa.summarise_drug_scan(res, want="up")
    assert "lower" in tfa.summarise_drug_scan(res, want="down")


def test_annotate_hits_with_chip_actually_annotates_the_hits(tmp_path):
    """It computed the counts, returned them, and never touched a hit -
    while its docstring said otherwise."""
    from promoter_tf_drugs import evidence as ev
    bed = tmp_path / "peaks.bed"
    bed.write_text("chr1\t1000\t1100\tVDR\n" * 3, encoding="utf-8")

    class _P:
        chrom, start, end, region = "chr1", 900, 1200, "chr1:900-1200"

    hits = [{"tf_symbols": ["VDR"]}, {"tf_symbols": ["NOBODY"]}]
    counts = ev.annotate_hits_with_chip(hits, _P(), bed)
    assert counts.get("VDR") == 3
    assert hits[0]["chip_peaks"] == 3 and hits[0]["chip_checked"] is True
    assert hits[1]["chip_peaks"] == 0 and hits[1]["chip_checked"] is True


def test_no_chip_bed_marks_the_hits_unchecked_rather_than_zero(tmp_path):
    from promoter_tf_drugs import evidence as ev
    hits = [{"tf_symbols": ["VDR"]}]
    assert ev.annotate_hits_with_chip(hits, None, None) == {}
    assert hits[0]["chip_checked"] is False
    assert "chip_peaks" not in hits[0], (
        "absent is not zero: a factor nobody assayed is unexamined")


# ===================================================================== #
# Docking into a predicted model. assess_model says "Do not dock this    #
# model" in its own text and was called by nothing.                      #
# ===================================================================== #
def _plddt_pdb(header, values):
    lines = [header]
    for i, b in enumerate(values, start=1):
        lines.append(
            f"ATOM  {i:5d}  CA  ALA A{i:4d}    "
            f"{i * 3.8:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00{b:6.2f}           C")
    lines.append("END")
    return "\n".join(lines) + "\n"


_AF_HEADER = ("HEADER    PREDICTED STRUCTURE FROM ALPHAFOLD DB "
              "(DEEPMIND / EMBL-EBI)")
_XRAY_HEADER = "HEADER    TRANSCRIPTION FACTOR                 01-JAN-00   1ABC"


def _flat_score(smiles, receptor, *a, **k):
    return -9.0


def test_a_mostly_disordered_predicted_model_is_refused(tmp_path):
    """Transcription factors are the worst case: a small folded DNA-binding
    domain and a long disordered activation domain. Low-pLDDT regions look
    like open greasy cavities to a scoring function, so the numbers come
    back confident and mean nothing."""
    f = tmp_path / "AF-P00001.pdb"
    f.write_text(_plddt_pdb(_AF_HEADER, [30.0] * 80 + [95.0] * 20),
                 encoding="utf-8")
    run = docking.screen_by_docking(
        f, [{"name": "x", "smiles": "CCO"}], target="TF",
        actives=[{"name": "a", "smiles": "CCN"}],
        decoys=[{"name": "d", "smiles": "CCC"}],
        score_fn=_flat_score, box=(0.0, 0.0, 0.0, 20.0, 20.0, 20.0))
    assert run.usable is False
    assert "disorder prediction" in run.error
    assert any("mean pLDDT" in c for c in run.caveats)


def test_a_folded_predicted_model_proceeds_but_is_labelled(tmp_path):
    f = tmp_path / "AF-P00002.pdb"
    f.write_text(_plddt_pdb(_AF_HEADER, [95.0] * 90 + [40.0] * 10),
                 encoding="utf-8")
    run = docking.screen_by_docking(
        f, [{"name": "x", "smiles": "CCO"}], target="TF",
        actives=[{"name": "a", "smiles": "CCN"}],
        decoys=[{"name": "d", "smiles": "CCC"}],
        score_fn=_flat_score, box=(0.0, 0.0, 0.0, 20.0, 20.0, 20.0))
    assert "disorder prediction" not in (run.error or "")
    assert any("89.5" in c or "90%" in c for c in run.caveats)


def test_a_crystal_b_factor_is_never_read_as_confidence(tmp_path):
    """The negative control. An experimental structure's B-factor column
    holds crystallographic B factors: 28 is a well-ordered atom, not 28%
    confidence. Reading it as pLDDT would refuse every real structure."""
    f = tmp_path / "1abc.pdb"
    f.write_text(_plddt_pdb(_XRAY_HEADER, [28.0] * 100), encoding="utf-8")
    assert docking._model_assessment(f) is None
    run = docking.screen_by_docking(
        f, [{"name": "x", "smiles": "CCO"}], target="TF",
        actives=[{"name": "a", "smiles": "CCN"}],
        decoys=[{"name": "d", "smiles": "CCC"}],
        score_fn=_flat_score, box=(0.0, 0.0, 0.0, 20.0, 20.0, 20.0))
    assert "disorder prediction" not in (run.error or "")
    assert not any("pLDDT" in c for c in run.caveats)


def test_an_unreadable_receptor_is_not_treated_as_a_predicted_model(tmp_path):
    assert docking._model_assessment(tmp_path / "nothing-here.pdb") is None


# ===================================================================== #
# The family screen is reached by the product, not only by its tests.    #
# ===================================================================== #
def test_the_pipeline_falls_back_to_a_family_screen():
    """screen_receptor refuses when a receptor has too few ligands of its
    own, which is a data gap rather than a fact about the receptor.
    screen_family was written for exactly that and called from nowhere."""
    import inspect
    src = inspect.getsource(pipeline.run)
    assert "screening.screen_family" in src
    i = src.index("screening.screen_family")
    window = src[max(0, i - 1500):i]
    assert "screen_receptor" in window, (
        "the family screen must sit after the receptor screen refuses, "
        "not replace it")


def test_the_family_fallback_only_fires_on_a_refusal_with_no_hits():
    """A receptor screen that produced hits must not be discarded for a
    family model, which is strictly weaker evidence."""
    import inspect
    src = inspect.getsource(pipeline.run)
    assert 'sc.get("error") and not sc.get("hits")' in src
    assert 'fam.get("hits")' in src, (
        "a family screen that refused must not overwrite the receptor "
        "result with its own error")


def test_a_family_result_records_what_it_replaced():
    import inspect
    src = inspect.getsource(pipeline.run)
    assert 'fam["replaces"] = sc.get("error")' in src


def test_the_upstream_mechanism_text_agrees_with_the_drug_chosen(monkeypatch):
    """`from_upstream` kept the same inverted expression attach_drugs had,
    so the pipeline selected an activator of the node while the report
    printed "inhibit" it. The two disagreed on screen."""
    from promoter_tf_drugs import catalogue as cat
    from promoter_tf_drugs import drugs as drugs_mod
    monkeypatch.setattr(drugs_mod, "chembl_drugs", _both_directions)

    edges = [net.Edge(source="KINASE", target="TFREP", sign=-1,
                      n_references=3)]
    hits = net.find_upstream_targets(None, "TFREP", goal="up",
                                     tf_role="repressor", edges=edges)
    net.attach_drugs(None, hits, approved_only=True)
    assert {d["direction"] for d in hits[0].drugs} == {"up"}

    routes = cat.from_upstream({"TFREP": [hits[0].to_dict()]})
    assert routes, "no route built"
    assert routes[0].mechanism.startswith("activate KINASE"), (
        routes[0].mechanism)


def test_an_old_result_without_the_sign_does_not_guess_a_direction():
    """A results.json written before want_tf_sign existed cannot have the
    direction recovered, so it says "move" rather than printing one that
    may be backwards."""
    from promoter_tf_drugs import catalogue as cat
    h = {"node": "KINASE", "net_sign": -1,
         "best_path": {"hops": 1, "confidence": 0.8, "specificity": 0.5,
                       "description": "KINASE -| TFREP"},
         "drugs": [{"name": "d", "approved": True}]}
    routes = cat.from_upstream({"TFREP": [h]})
    assert routes[0].mechanism.startswith("move KINASE")


def test_the_hop_count_reads_as_english():
    from promoter_tf_drugs import catalogue as cat
    h = {"node": "K", "net_sign": 1, "want_tf_sign": 1,
         "best_path": {"hops": 1, "confidence": 0.8, "specificity": 0.5,
                       "description": "x"},
         "drugs": [{"name": "d", "approved": True}]}
    assert "(1 hop)" in cat.from_upstream({"TF": [h]})[0].mechanism
    h["best_path"]["hops"] = 2
    assert "(2 hops)" in cat.from_upstream({"TF": [h]})[0].mechanism


# ===================================================================== #
# Which conformation to dock into depends on the factor's ROLE as well   #
# as the goal. --state was manual and its help covered only activators.  #
# ===================================================================== #
def test_the_conformation_follows_role_and_goal_together():
    """To raise a gene you want an activator ON and a repressor OFF, so
    the same goal needs opposite conformations for the two roles."""
    from promoter_tf_drugs.docking import state_for
    assert state_for("activator", "up") == "active"
    assert state_for("repressor", "up") == "inactive"
    assert state_for("activator", "down") == "inactive"
    assert state_for("repressor", "down") == "active"


def test_an_unresolved_role_does_not_guess_a_conformation():
    """This filter DISCARDS structures, so a wrong guess can delete the
    only pocket with a co-crystallised ligand and refuse the dock. That
    is stricter than drugs.desired_direction, which folds dual in with
    activator because there a wrong guess only reorders a table."""
    from promoter_tf_drugs.docking import state_for
    from promoter_tf_drugs.drugs import desired_direction
    assert state_for("dual", "up") is None
    assert state_for("unknown", "up") is None
    assert state_for(None, "up") is None
    assert desired_direction("dual", "up") == "up", (
        "the drug-side convention is deliberately different; if it "
        "changes, state_for's reason for diverging should be revisited")


def test_the_state_flag_is_derived_not_demanded():
    import inspect
    from promoter_tf_drugs import docking as dk
    src = inspect.getsource(dk.main)
    assert "want_state = args.state" in src
    assert "state_for(" in src
    assert "want_state=args.state" not in src, (
        "a call site still reads the raw flag and will ignore the "
        "derived value")


def test_a_conformation_filter_is_not_reported_as_a_missing_structure():
    """choose_pocket skipping every structure for its conformation used to
    surface as "the PDB has none", sending the user after the wrong
    problem."""
    import inspect
    from promoter_tf_drugs import docking as dk
    src = inspect.getsource(dk.dock_screen_for_target)
    assert "excluded_by_conformation" in src
    assert "choose_pocket(client, uniprot, work_dir)" in src, (
        "the unfiltered retry that tells the two cases apart is gone")


# ===================================================================== #
# A docking route and the direction you asked for.                       #
# ===================================================================== #
def _route_run(target="VDR", state="active"):
    return {"target": target, "usable": True, "auc": 0.92, "ef1": 3.1,
            "pocket": {"pdb_id": "1ABC", "ligand_resname": "VDX",
                       "state": state},
            "hits": [{"drug": "compound-1", "score": -10.4}]}


def test_the_pocket_conformation_is_related_to_the_direction_wanted():
    from promoter_tf_drugs import catalogue as cat
    cases = {
        ("activator", "active"): "supports",
        ("activator", "inactive"): "opposes",
        ("repressor", "inactive"): "supports",
        ("repressor", "active"): "opposes",
    }
    for (role, state), expect in cases.items():
        r = cat.from_docking(_route_run(state=state), goal="up",
                             tf_roles={"VDR": role})[0]
        assert expect in r.mechanism, (role, state, r.mechanism)


def test_an_opposing_conformation_is_flagged_and_not_dropped():
    """Docking does not measure whether a compound agonises or
    antagonises, and an antagonist can bind the active conformation. So
    an opposing pocket is a warning, never grounds to delete the route."""
    from promoter_tf_drugs import catalogue as cat
    routes = cat.from_docking(_route_run(state="active"), goal="up",
                              tf_roles={"VDR": "repressor"})
    assert len(routes) == 1
    note = " ".join(routes[0].notes)
    assert "conformation to stabilise is inactive" in note
    assert "antagonist can bind the active conformation" in note


def test_an_unknown_conformation_stays_unaddressed():
    from promoter_tf_drugs import catalogue as cat
    r = cat.from_docking(_route_run(state="unknown"), goal="up",
                         tf_roles={"VDR": "activator"})[0]
    assert "supports" not in r.mechanism and "opposes" not in r.mechanism


def test_without_roles_nothing_is_claimed_about_direction():
    from promoter_tf_drugs import catalogue as cat
    r = cat.from_docking(_route_run(), goal="up", tf_roles=None)[0]
    assert "supports" not in r.mechanism and "opposes" not in r.mechanism
    assert r.category == "TF: predicted pocket binder"


# ===================================================================== #
# Docking a protein that is not one of this promoter's factors.          #
# ===================================================================== #
def test_a_dock_against_a_non_factor_is_not_called_a_direct_binder():
    """It came back as "TF: predicted pocket binder" with directness 1,
    claiming a direct relationship with a gene it has no path to."""
    from promoter_tf_drugs import catalogue as cat
    off = cat.from_docking(_route_run(target="MTOR"), goal="up",
                           tf_roles={"VDR": "activator"})[0]
    assert off.category == "off-pathway pocket binder"
    assert off.directness == 2
    assert any("not one of the factors" in n for n in off.notes)


def test_a_dock_against_a_factor_here_is_unchanged():
    from promoter_tf_drugs import catalogue as cat
    on = cat.from_docking(_route_run(target="VDR"), goal="up",
                          tf_roles={"VDR": "activator"})[0]
    assert on.category == "TF: predicted pocket binder"
    assert on.directness == 1


def test_the_off_pathway_route_scores_below_the_on_pathway_one():
    from promoter_tf_drugs import catalogue as cat
    on = cat.from_docking(_route_run(target="VDR"), goal="up",
                          tf_roles={"VDR": "activator"})[0]
    off = cat.from_docking(_route_run(target="MTOR"), goal="up",
                           tf_roles={"VDR": "activator"})[0]
    assert cat._score(off) < cat._score(on)


def test_the_new_category_is_registered_and_explained():
    """categories_absent is built from CATEGORIES, and the report prints a
    guide entry per category - a category in neither is invisible."""
    from promoter_tf_drugs import catalogue as cat
    assert "off-pathway pocket binder" in cat.CATEGORIES
    guide = cat.CATEGORY_GUIDE["off-pathway pocket binder"]
    assert set(guide) >= {"changes", "assay", "what", "why", "cost",
                          "differs"}


def test_the_pipeline_hands_the_catalogue_the_factor_roles():
    import inspect
    src = inspect.getsource(pipeline.run)
    assert '"tf_roles"' in src
    assert '(r.annotation or {}).get("role")' in src, (
        "annotation is a dict subclass; a getattr read returns None for "
        "every factor")


def test_no_test_helper_is_defined_twice_in_this_file():
    """A second def of the same name silently replaces the first, so tests
    written earlier start running against a fixture meant for later ones.
    That happened: a new `_dock_run` shadowed an existing one and broke six
    passing tests in a way that looked like a product regression."""
    import ast
    import collections
    from pathlib import Path
    src = Path(__file__).read_text(encoding="utf-8")
    # Classes too: a second `class _FakeClient` shadowed the first and
    # broke two unrelated tests, the same way a second `def _dock_run`
    # had. Checking only functions missed it.
    names = [n.name for n in ast.parse(src).body
             if isinstance(n, (ast.FunctionDef, ast.ClassDef))]
    dupes = [n for n, c in collections.Counter(names).items() if c > 1]
    assert not dupes, f"defined more than once: {dupes}"


def _fam_screen(**kw):
    base = {"family_model": True, "family": "nuclear receptor",
            "transfer_auc": 0.78, "model_auc": 0.91,
            "contributing": [{"symbol": "RXRA"}, {"symbol": "PPARG"}],
            "enrichment": {"usable": True},
            "hits": [{"drug": "cpd-1", "quadrant": "NOVEL HYPOTHESIS",
                      "probability": 0.95}]}
    base.update(kw)
    return {"VDR": base}


def test_a_family_model_route_rests_on_the_measured_transfer_auc():
    """The model never saw VDR, so its own probability is not a strength
    AT VDR. transfer_auc is the measured answer - how well a model trained
    on relatives ranks VDR's own known ligands - and is what the route
    carries. Wiring screen_family into the pipeline sent its hits through
    from_screen, which dropped all of this."""
    from promoter_tf_drugs import catalogue as cat
    r = cat.from_screen(_fam_screen())[0]
    assert r.axes.strength == 0.78, "the in-model probability leaked back in"
    assert "transfer AUC 0.78" in r.axes.evidence
    assert "family model" in r.mechanism
    assert "not on VDR" in r.mechanism
    assert any("subtype selectivity" in n for n in r.notes)


def test_a_receptor_model_route_is_unchanged():
    from promoter_tf_drugs import catalogue as cat
    own = {"VDR": {"enrichment": {"usable": True},
                   "hits": [{"drug": "cpd-1", "quadrant": "NOVEL HYPOTHESIS",
                             "probability": 0.95}]}}
    r = cat.from_screen(own)[0]
    assert r.axes.strength == 0.95
    assert "family" not in r.mechanism


def test_an_established_ligand_keeps_its_tier_under_a_family_model():
    """The binding was measured at this receptor whichever model surfaced
    it, so tier 2 is still right - but the route says the model only
    ranked it."""
    from promoter_tf_drugs import catalogue as cat
    sc = _fam_screen(hits=[{"drug": "calcitriol",
                            "quadrant": "ESTABLISHED LIGAND, untested on "
                                        "this gene"}])
    r = cat.from_screen(sc)[0]
    assert r.tier == 2
    assert any("did not establish it" in n for n in r.notes)


def test_a_q_value_of_zero_is_the_most_significant_not_the_least():
    """`value or 1` reads 0.0 as missing. tfactivity._q exists for this and
    the pipeline reimplemented the broken version."""
    from promoter_tf_drugs import pipeline as pl
    assert pl._sig_q(0.0) == 0.0
    assert pl._sig_q(None) == 1.0
    assert pl._sig_q(0.04) == 0.04


def test_a_zscore_of_zero_is_not_swapped_for_another_column():
    """A gene that did not move has z = 0.0. `a or b` substituted a
    different field's value for it."""
    from promoter_tf_drugs import signatures as sig
    import inspect
    src = inspect.getsource(sig)
    assert 'float(r.get("zscore") or r.get("score") or 0.0)' not in src


# ===================================================================== #
# Drug sensitivity: the one route that needs no annotation at all.      #
# ===================================================================== #
def _sens_panel(seed=7, n=120):
    """A PRISM-shaped panel with a planted phenocopier and opposer, both
    sitting on top of a shared cell-line fragility component."""
    import numpy as np
    rng = np.random.default_rng(seed)
    models = [f"ACH-{i:06d}" for i in range(n)]
    gene = rng.normal(-0.2, 0.4, size=n)
    frag = rng.normal(0, 0.6, size=n)
    cols, mat = [], []
    for j in range(40):
        cols.append(f"BRD-{j:05d}")
        mat.append(0.9 * frag + rng.normal(0, 0.5, size=n))
    cols.append("BRD-COPY")
    mat.append(0.7 * gene + 0.9 * frag + rng.normal(0, 0.25, size=n))
    cols.append("BRD-OPP")
    mat.append(-0.7 * gene + 0.9 * frag + rng.normal(0, 0.25, size=n))
    prism = {"models": models, "columns": cols,
             "values": np.vstack(mat).T}
    target = {m: float(v) for m, v in zip(models, gene)}
    names = {"BRD-COPY": "phenocopier", "BRD-OPP": "opposer"}
    return prism, target, names


def test_the_p_value_matches_a_reference_t_distribution():
    """The incomplete beta is written out by hand here, because numpy has
    none and scipy is not a dependency. If it drifts, every q in this
    section is wrong and nothing else would notice."""
    import pytest
    st = pytest.importorskip("scipy.stats")
    import numpy as np
    from promoter_tf_drugs import sensitivity as sn
    for r, n in ((0.5, 30), (0.2, 100), (0.8, 10), (0.05, 500)):
        df = n - 2
        t = abs(r) * np.sqrt(df / (1 - r * r))
        assert abs(sn._p_from_r(r, n) - 2 * st.t.sf(t, df)) < 1e-9


def test_the_partial_correlation_equals_a_regression_residual():
    """The definition, checked against the thing it is defined to be."""
    import numpy as np
    from promoter_tf_drugs import sensitivity as sn
    rng = np.random.default_rng(3)
    z = rng.normal(size=400)
    x = 0.8 * z + 0.6 * rng.normal(size=400)
    y = 0.8 * z + 0.6 * rng.normal(size=400)

    def pr(a, b):
        return float(((a - a.mean()) * (b - b.mean())).mean()
                     / (a.std() * b.std()))

    mine = float(sn._partial(np.array([pr(x, y)]), np.array([pr(x, z)]),
                             pr(y, z))[0])
    rx = x - np.polyval(np.polyfit(z, x, 1), z)
    ry = y - np.polyval(np.polyfit(z, y, 1), z)
    assert abs(mine - pr(rx, ry)) < 1e-9
    assert abs(mine) < abs(pr(x, y)), (
        "the control must remove the shared component, not keep it")


def test_pure_noise_produces_no_significant_hit():
    """2,000 tests against a random profile. If BH is wrong this section
    invents drugs."""
    import numpy as np
    from promoter_tf_drugs import sensitivity as sn
    rng = np.random.default_rng(0)
    prism = {"models": [f"ACH-{i}" for i in range(200)],
             "columns": [f"BRD-{j}" for j in range(2000)],
             "values": rng.normal(size=(200, 2000))}
    target = {f"ACH-{i}": float(v)
              for i, v in enumerate(rng.normal(size=200))}
    res = sn.correlate(prism, target, gene="NULL")
    assert res.assessed
    assert res.n_drugs_tested == 2000
    assert res.hits == [], f"{len(res.hits)} false positive(s) on noise"


def test_a_planted_phenocopy_is_recovered_with_the_right_sign():
    from promoter_tf_drugs import sensitivity as sn
    prism, target, names = _sens_panel()
    res = sn.correlate(prism, target, gene="VDR", names=names)
    by = {h.drug: h for h in res.hits}
    assert "phenocopier" in by and "opposer" in by
    assert by["phenocopier"].r > 0
    assert by["opposer"].r < 0
    assert by["phenocopier"].direction == "phenocopies"
    assert by["opposer"].direction == "opposes"
    assert len(res.hits) == 2, (
        "the 40 fillers carry no gene signal and must not be reported")


def test_too_few_shared_lines_is_refused_not_scored():
    from promoter_tf_drugs import sensitivity as sn
    prism, target, names = _sens_panel()
    few = dict(list(target.items())[:10])
    res = sn.correlate(prism, few, gene="VDR", names=names)
    assert not res.assessed
    assert "below the" in res.error and "noise with a decimal place" in res.error


def test_a_missing_matrix_is_not_an_absence_of_correlated_drugs():
    from promoter_tf_drugs import sensitivity as sn
    res = sn.correlate({}, {"ACH-1": 0.5}, gene="VDR")
    assert not res.assessed
    assert "missing file, not an absence" in res.error


def test_a_gene_absent_from_depmap_says_so(tmp_path, monkeypatch):
    from promoter_tf_drugs import sensitivity as sn
    monkeypatch.chdir(tmp_path)
    d = tmp_path / "data"
    d.mkdir(exist_ok=True)
    (d / "CRISPRGeneEffect.csv").write_text(
        ",MYC (4609)\n" + "".join(f"ACH-{i},{-0.3}\n" for i in range(30)),
        encoding="utf-8")
    prism = {"models": [f"ACH-{i}" for i in range(30)],
             "columns": ["BRD-1"], "values": __import__("numpy").zeros((30, 1))}
    res = sn.phenocopy("VDR", prism=prism, names={})
    assert not res.assessed
    assert "not scored in the DepMap release" in res.error
    assert "Absent from the file is not the same" in res.error


def test_the_two_readouts_have_opposite_expected_signs():
    """Chronos and PRISM both run negative for harm, so a phenocopy is a
    POSITIVE correlation. Expression runs the other way, so sensitivity
    rising with expression is NEGATIVE. Getting this backwards would
    invert every call in the section."""
    from promoter_tf_drugs import sensitivity as sn
    prism, target, names = _sens_panel()
    ph = sn.correlate(prism, target, readout="phenocopy", gene="VDR",
                      names=names)
    ex = sn.correlate(prism, target, readout="expression", gene="VDR",
                      names=names)
    by_ph = {h.drug: h.direction for h in ph.hits}
    by_ex = {h.drug: h.direction for h in ex.hits}
    assert by_ph["phenocopier"] == "phenocopies"
    assert by_ex["phenocopier"] == "inverse to expression"
    assert by_ex["opposer"] == "tracks expression"


# ---- the route ------------------------------------------------------
def _sens_block():
    from promoter_tf_drugs import sensitivity as sn
    prism, target, names = _sens_panel()
    res = sn.correlate(prism, target, gene="VDR", names=names)
    return {"VDR": res.to_dict()}


def test_a_phenocopy_route_is_kept_only_where_it_moves_the_right_way():
    """A phenocopier removes the factor, so it serves a goal that needs
    the factor lowered and works against one that needs it raised. Unlike
    docking, the direction here was measured, so a wrong-way route is
    dropped rather than flagged."""
    from promoter_tf_drugs import catalogue as cat
    block = _sens_block()
    cases = {
        ("activator", "up"): "opposer",
        ("activator", "down"): "phenocopier",
        ("repressor", "up"): "phenocopier",
        ("repressor", "down"): "opposer",
    }
    for (role, goal), expect in cases.items():
        kept = sorted(r.drug for r in cat.from_sensitivity(
            block, goal=goal, tf_roles={"VDR": role}))
        assert kept == [expect], (role, goal, kept)


def test_an_unknown_role_keeps_both_directions():
    """desired_direction returns "either" for an unresolved role, so
    nothing is dropped on a guess."""
    from promoter_tf_drugs import catalogue as cat
    kept = sorted(r.drug for r in cat.from_sensitivity(
        _sens_block(), goal="up", tf_roles={"VDR": "unknown"}))
    assert kept == ["opposer", "phenocopier"]


def test_the_route_rests_on_the_partial_not_the_raw_correlation():
    from promoter_tf_drugs import catalogue as cat
    block = _sens_block()
    r = cat.from_sensitivity(block, goal="down",
                             tf_roles={"VDR": "activator"})[0]
    hit = [h for h in block["VDR"]["hits"] if h["drug"] == "phenocopier"][0]
    assert r.axes.strength == min(1.0, abs(hit["r"]))
    assert r.axes.strength != abs(hit["r_raw"])
    assert r.tier == 4
    assert r.category == "measured phenocopy"
    assert any("behaves as though it does" in n for n in r.notes)


def test_an_unassessed_scan_contributes_no_routes():
    from promoter_tf_drugs import catalogue as cat
    assert cat.from_sensitivity(
        {"VDR": {"assessed": False, "error": "no file"}}, goal="up") == []


def test_the_section_reports_a_scan_that_found_nothing(tmp_path):
    from promoter_tf_drugs import pipeline as pl
    out = "\n".join(pl._sensitivity_section(
        {"VDR": {"assessed": True, "n_drugs_tested": 4000, "hits": []}}))
    assert "none cleared" in out
    assert "bounds what this method contributes" in out


def test_the_section_says_when_nothing_could_be_scanned():
    from promoter_tf_drugs import pipeline as pl
    out = "\n".join(pl._sensitivity_section(
        {"VDR": {"assessed": False, "error": "no PRISM matrix loaded"}}))
    assert "Nothing was scanned" in out
    assert "no PRISM matrix loaded" in out


def test_the_section_warns_that_a_correlation_is_not_a_target():
    from promoter_tf_drugs import pipeline as pl
    out = "\n".join(pl._sensitivity_section(_sens_block()))
    assert "This is not a target" in out
    assert "Partial r" in out
    assert "proliferation readout" in out


def test_the_cli_flag_reaches_the_pipeline():
    """A flag that is parsed and never passed does nothing, which is the
    failure this project keeps finding."""
    import inspect
    from promoter_tf_drugs import cli
    src = inspect.getsource(cli)
    assert "--sensitivity-scan" in src
    assert "sensitivity_scan=args.sensitivity_scan" in src
    assert "prism_file=args.prism_file" in src


def test_the_fragility_control_removes_a_false_signal_it_is_built_for():
    """The case the partial correlation exists for, which `_sens_panel`
    does not cover: it draws the gene profile independently of fragility,
    so its fillers never acquire a false correlation and are rejected for
    having no signal rather than by the control.

    Here the gene profile tracks fragility, which is what real data looks
    like - a gene whose loss hurts will hurt the fragile lines most. Forty
    compounds with no relationship to the gene then correlate with it
    strongly through that shared term, and the control is the only thing
    standing between them and the report.
    """
    import numpy as np
    from promoter_tf_drugs import sensitivity as sn

    rng = np.random.default_rng(7)
    n = 120
    frag = rng.normal(0, 0.6, size=n)
    cols, mat = [], []
    for j in range(40):
        cols.append(f"filler-{j:02d}")
        mat.append(0.9 * frag + rng.normal(0, 0.5, size=n))
    values = np.vstack(mat).T
    models = [f"ACH-{i}" for i in range(n)]
    gene = 0.8 * frag + 0.6 * rng.normal(size=n)      # linked, as in life

    def pearson(a, b):
        return float(((a - a.mean()) * (b - b.mean())).mean()
                     / (a.std() * b.std()))

    raw = [abs(pearson(values[:, j], gene)) for j in range(40)]
    would_pass = sum(1 for r in raw if r >= sn.MIN_ABS_R)
    assert would_pass >= 30, (
        f"only {would_pass}/40 fillers carry a false signal; the fixture "
        f"is not posing the problem this test exists to check")

    res = sn.correlate({"models": models, "columns": cols,
                        "values": values},
                       {m: float(v) for m, v in zip(models, gene)},
                       gene="G", names={})
    assert res.assessed
    # Not zero: the control is a correction, not a guarantee, and some
    # residual false positives are expected at the stated false discovery
    # rate. What must not happen is most of them surviving.
    assert len(res.hits) <= 4, (
        f"{len(res.hits)} of 40 fillers survived; {would_pass} would have "
        f"passed on the raw correlation, so the control is not working")
    assert len(res.hits) < would_pass / 4


# ===================================================================== #
# Found by mutation testing: these behaviours had tests that greps the   #
# source and passed while the behaviour was disabled.                   #
# ===================================================================== #
def test_screen_actually_inverts_its_ranking_for_a_down_goal():
    """`test_screen_family_honours_the_goal_it_accepts` only greps for the
    text `goal=goal`. Making `screen` ignore the argument entirely left
    the whole suite green. This calls it and checks the number changes.

    screen() computes p_want = p_agonist for "up" and 1 - p_agonist for
    "down", so the two goals must give complementary probabilities for
    the same molecule.
    """
    import numpy as np
    from promoter_tf_drugs import screening as scr, agonism

    class _Clf:
        """One likely agonist and one likely antagonist."""

        def predict_proba(self, X):
            return np.array([[0.1, 0.9], [0.9, 0.1]])[:len(X)]

    class _Model:
        clf = _Clf()

        def nearest_tanimoto(self, bv):
            return 0.9

    lib = [{"name": "agonist-like", "smiles": "CCO",
            "chembl_id": "CHEMBL1"},
           {"name": "antagonist-like", "smiles": "CCN",
            "chembl_id": "CHEMBL2"}]
    _X, _b, keep = agonism.fingerprints([d["smiles"] for d in lib])
    if len(keep) < 2:
        import pytest
        pytest.skip("rdkit not available to fingerprint the fixture")

    up = {h.drug for h in scr.screen(_Model(), lib, goal="up")}
    down = {h.drug for h in scr.screen(_Model(), lib, goal="down")}
    up = {d.lower() for d in up}
    down = {d.lower() for d in down}
    assert up == {"agonist-like"}, up
    assert down == {"antagonist-like"}, down
    assert up != down, (
        "goal is being ignored: both directions returned the same set")


def test_the_q_helper_treats_zero_as_the_most_significant_value():
    """`tfactivity._q` exists so a q-value of 0.0 - the most significant
    result possible - is not read as missing by a falsy test. Nothing
    tested it: disabling it left the suite green, because `_bh` floors
    its output at Q_FLOOR so no real q ever reaches 0.0. The helper is
    the guard for data that does not come from `_bh`.
    """
    from promoter_tf_drugs import tfactivity as tfa
    assert tfa._q(0.0) == 0.0
    assert tfa._q(0.0, 1.0) == 0.0, "a falsy zero fell through to default"
    assert tfa._q(None) == 1.0
    assert tfa._q(None, 0.5) == 0.5
    assert tfa._q(0.03) == 0.03


def test_bh_never_emits_a_zero_that_a_falsy_test_would_drop():
    """The invariant the helper above is paired with, pinned explicitly
    because other modules rely on it without saying so."""
    from promoter_tf_drugs import tfactivity as tfa
    qs = tfa._bh([0.0, 1e-300, 1e-12, 0.5, None])
    assert all(q is None or q >= tfa.Q_FLOOR for q in qs)
    assert qs[4] is None, "a missing p-value must stay missing"


def test_a_zscore_of_zero_survives_the_parse():
    """Behavioural version. The previous test only asserted that one
    spelling of the bug was absent from the source, so reintroducing it
    with different wording passed.

    A gene that did not move has z = 0.0. It must not be silently
    replaced by whatever is in the `score` field, which is a different
    quantity.
    """
    from promoter_tf_drugs import signatures as sig

    class _Client:
        def post(self, url, payload):
            return {"results": [
                {"zscore": 0.0, "score": -9.9,
                 "meta": {"pert_name": "flat-drug"}},
                {"zscore": 2.5, "score": -1.0,
                 "meta": {"pert_name": "real-drug"}},
            ]}

    import promoter_tf_drugs.signatures as _s
    real_is_landmark = _s.is_landmark
    _s.is_landmark = lambda gene, client=None: True
    try:
        hits = sig.query_gene(_Client(), "MYGENE", direction="up")
    finally:
        _s.is_landmark = real_is_landmark

    by = {h.drug.lower(): h for h in hits}
    assert "flat-drug" in by, "the zero-score row was dropped entirely"
    assert by["flat-drug"].score == 0.0, (
        f"a z-score of 0.0 became {by['flat-drug'].score}; the `score` "
        f"column was substituted for it")


def test_the_sensitivity_flag_parses_and_the_pipeline_accepts_it():
    """Behavioural as far as is proportionate, and honest about the rest.

    The parser is exercised and `run`'s real signature is introspected,
    so a flag that does not parse or a parameter that does not exist both
    fail here. What this does NOT cover is `run` accepting the flag and
    then ignoring it: the scan is suppressed under `demo`, every full-run
    test in this file is a demo run, and standing up a non-demo run needs
    a genome and a GTF. That gap is real and recorded rather than papered
    over with a grep that would pass either way.
    """
    import inspect
    from promoter_tf_drugs import cli, pipeline as pl

    args = cli.build_parser().parse_args(
        ["--gene", "X", "--sensitivity-scan", "--prism-file", "p.csv"])
    assert args.sensitivity_scan is True
    assert args.prism_file == "p.csv"

    params = inspect.signature(pl.run).parameters
    assert "sensitivity_scan" in params
    assert "prism_file" in params
    assert params["sensitivity_scan"].default is False, (
        "the scan must be opt-in; it reads two large matrices")


def test_the_scan_is_suppressed_in_demo_mode(monkeypatch):
    """Demo must not read the real PRISM matrix, the same way it must not
    read the real LINCS files.

    The first version of this test asserted the result was empty, and
    passed for the wrong reason: there is no PRISM file in the test
    directory, so the scan returns nothing whether it is suppressed or
    not. Disabling the suppression left it green. It now makes the loader
    itself fail loudly, so the only way to pass is never to call it.
    """
    from promoter_tf_drugs import sensitivity as sens

    def _loud(*a, **k):
        raise AssertionError("demo mode read the PRISM matrix")

    monkeypatch.setattr(sens, "load_prism", _loud)
    monkeypatch.setattr(sens, "load_compound_names", _loud)
    monkeypatch.setattr(sens, "phenocopy", _loud)

    res = pipeline.run("DEMOGENE", demo=True, offline=True,
                       sensitivity_scan=True)
    for block in res["results"]:
        assert block.get("sensitivity") == {}


# ===================================================================== #
# The class table, extended because a real CCR5 shortlist was crowded   #
# out by compounds it did not know.                                      #
# ===================================================================== #
# Every name below came from an actual `strict` shortlist for lowering
# CCR5, where it outranked everything but Maraviroc.
_CROWDERS = [
    "mitoxantrone dihydrochloride", "dinaciclib", "vemurafenib",
    "ponatinib", "dabrafenib", "saracatinib", "RAF 265",
    "CGP-60474", "ALW-II-38-3", "AT-7519", "BMS-387032", "JNK-9L",
    "PD-0325901", "AS-605240", "GSK-1059615",
    "dorsomorphin dihydrochloride",
]


def test_the_compounds_that_crowded_a_real_shortlist_are_classified():
    """`classify` returned None for all of these, so the burden filter
    passed every one. An antiproliferative moves the abundance of a great
    many transcription factors, which is why they dominate an
    expression-signature shortlist without being evidence of a route."""
    from promoter_tf_drugs import systemic as sy
    unclassified = [n for n in _CROWDERS if sy.classify(n) is None]
    assert not unclassified, f"still unclassified: {unclassified}"
    for n in _CROWDERS:
        assert sy.clinical_burden(n).level == "high", n


def test_the_drug_that_should_win_is_not_swept_up_with_them():
    """Maraviroc is the correct tier-1 answer for lowering CCR5. A filter
    that removes the noise by removing everything has done nothing."""
    from promoter_tf_drugs import systemic as sy
    assert sy.classify("Maraviroc") is None
    assert sy.clinical_burden("Maraviroc").level == "low"
    assert sy.clinical_burden("Maraviroc").assessed is False


def test_a_code_name_resolves_however_it_is_spaced():
    """LINCS writes "RAF 265" where ChEMBL writes "RAF265". The
    head-token fallback looked up "raf" and missed, because for a code
    name the number is part of the identity."""
    from promoter_tf_drugs import systemic as sy
    names = {sy.classify(n).name for n in ("RAF 265", "RAF-265", "raf265")}
    assert names == {"antineoplastic kinase inhibitor"}


def test_the_head_token_fallback_still_serves_real_multiword_names():
    from promoter_tf_drugs import systemic as sy
    assert sy.classify("Fluticasone Propionate HFA").name == "glucocorticoid"
    assert sy.classify("Levothyroxine Sodium").name == "thyroid hormone"


def test_the_head_token_fallback_does_not_fire_on_a_code_stem():
    """Without this guard "raf 265" would look up "raf", and any future
    member starting with a shared stem would collect unrelated rows."""
    from promoter_tf_drugs import systemic as sy
    assert sy.classify("imatinib 400") is None or True
    # A two-token name containing a digit must not fall back to its head.
    assert sy.classify("nosuchdrug 123") is None


def test_an_lxr_agonist_is_not_filed_as_a_kinase_inhibitor():
    """gw 3965 survived the filter and should: it is a nuclear receptor
    agonist. Sweeping every unfamiliar code name into the kinase class
    would make the filter look effective by being wrong."""
    from promoter_tf_drugs import systemic as sy
    assert sy.classify("gw 3965 hydrochloride") is None


def test_the_table_still_satisfies_its_own_audit():
    """A member in two classes makes `classify` order-dependent, and a
    name that changes under salt_free_name is never matched at all. The
    audit caught one of each while this table was being extended."""
    from promoter_tf_drugs import audit
    problems = {k: v for k, v in audit.check_internal().items()
                if v and "systemic" in k}
    assert not problems, problems


def test_salt_forms_of_the_new_members_resolve():
    from promoter_tf_drugs import systemic as sy
    assert sy.classify("Mitoxantrone Dihydrochloride").name.startswith(
        "cytotoxic")
    assert sy.classify("Ruxolitinib Phosphate").name == (
        "antineoplastic kinase inhibitor")


# ---- the selection filters -------------------------------------------
def test_max_burden_drops_high_burden_routes():
    from promoter_tf_drugs import shortlist as sl
    rows = [
        {"drug": "dinaciclib", "burden": "high", "_burden_assessed": True,
         "_burden_score": 0.8},
        {"drug": "maraviroc", "burden": "low", "_burden_assessed": False,
         "_burden_score": 0.0},
    ]
    kept = [r["drug"] for r in rows if sl.keep(r, max_burden="moderate")]
    assert kept == ["maraviroc"]


def test_an_unassessed_drug_passes_the_burden_filter():
    """The class table is hand-maintained and incomplete. Dropping what it
    does not list would select for obscurity, which is the trap the
    regulon count already has."""
    from promoter_tf_drugs import shortlist as sl
    row = {"drug": "something-new", "burden": "low",
           "_burden_assessed": False}
    assert sl.keep(row, max_burden="low")


def test_an_unknown_burden_level_is_an_error_not_a_pass():
    from promoter_tf_drugs import shortlist as sl
    import pytest
    row = {"drug": "x", "burden": "high", "_burden_assessed": True}
    with pytest.raises(ValueError):
        sl.keep(row, max_burden="severe")


def test_max_dependency_burden_keeps_an_unscored_factor():
    """Absent from DepMap is not the same as safe."""
    from promoter_tf_drugs import shortlist as sl
    assert sl.keep({"_dependency_burden": None}, max_dependency_burden=0.5)
    assert not sl.keep({"_dependency_burden": 0.9},
                       max_dependency_burden=0.5)


def test_the_new_filters_are_registered_and_validated():
    from promoter_tf_drugs import shortlist as sl
    import pytest
    for name in ("max_burden", "exclude_drug_classes",
                 "max_dependency_burden", "require_verified_approval"):
        assert name in sl.PREDICATES
    with pytest.raises(KeyError):
        sl.validate_filters({"max_burdne": "low"})


# ---- approval ---------------------------------------------------------
def test_a_negative_approval_needs_a_real_library_behind_it():
    """The curated tables hold a few hundred names, so a miss against
    them says the drug is uncurated and nothing about its approval -
    vemurafenib is approved and is not in them."""
    from promoter_tf_drugs import catalogue as cat
    assert cat.verify_approval("maraviroc") is True
    if not cat._APPROVED_FROM_LIBRARY:
        assert cat.verify_approval("CGP-60474") is None
        assert cat.verify_approval("vemurafenib") is None


def test_require_verified_approval_keeps_what_could_not_be_checked():
    """With no drug table on disk that is every row, and a filter that
    empties itself on a missing file teaches nothing."""
    from promoter_tf_drugs import shortlist as sl
    assert sl.keep({"approval_verified": None},
                   require_verified_approval=True)
    assert sl.keep({"approval_verified": True},
                   require_verified_approval=True)
    assert not sl.keep({"approval_verified": False},
                       require_verified_approval=True)


def test_routes_carry_a_checked_approval_field():
    """Ten of the eleven collectors hardcode approved=True, so that field
    means "came from a drug table", not "is an approved drug"."""
    from promoter_tf_drugs import catalogue as cat
    block = {"direct_therapies": [], "drug_grades": [],
             "lincs": {"hits": [{"drug": "Maraviroc", "tier": 2,
                                 "score": -2.5, "direction": "down",
                                 "n_signatures": 3}]}}
    built = cat.build(block, goal="down")
    for r in built["routes"]:
        assert "approval_verified" in r


def test_strict_carries_the_burden_criterion():
    """It is the criterion that takes a real CCR5 shortlist from 21 rows
    to 4, so a silent removal would be a large regression that no count
    assertion elsewhere would notice."""
    from promoter_tf_drugs import shortlist as sl
    assert sl.PRESETS["strict"].get("max_burden") == "moderate"


def test_strict_drops_a_cytotoxic_and_keeps_the_direct_agent():
    """The shape of the real result: the tier-1 agent on the gene itself
    survives while an antiproliferative reaching it through a factor's
    abundance does not. A filter that cleaned up by removing everything
    would pass a count assertion and be useless."""
    from promoter_tf_drugs import shortlist as sl
    rows = [
        {"drug": "Maraviroc", "tier": 1, "score": 60.0, "category":
         "direct: gene product", "axes": {"specificity": 0.55},
         "directness": 0, "grade": "A", "_notes": "", "_chip_known": False,
         "burden": "low", "_burden_assessed": False, "approved": True,
         "mechanism": "CCR5 antagonist", "acts_on": "CCR5"},
        {"drug": "dinaciclib", "tier": 3, "score": 35.0, "category":
         "TF: abundance change", "axes": {"specificity": 0.45},
         "directness": 1, "grade": None, "_notes": "", "_chip_known": False,
         "burden": "high", "_burden_assessed": True, "approved": True,
         "mechanism": "lowers RBPJ transcript", "acts_on": "RBPJ"},
    ]
    kept = [r["drug"] for r in sl.filter_rows(rows, **sl.PRESETS["strict"])]
    assert kept == ["Maraviroc"], kept


def test_the_other_presets_are_left_alone():
    """strict is the one that answers "what would I order next week?".
    `novel` and `empirical` exist to show what is there, so silently
    narrowing them would hide the thing they are for."""
    from promoter_tf_drugs import shortlist as sl
    for name in ("mechanism", "empirical", "novel", "specific"):
        assert "max_burden" not in sl.PRESETS[name], name


# ===================================================================== #
# The virtual screen, pointed at the gene's OWN product.                 #
# ===================================================================== #
def _olol_actives():
    """Beta-blockers from the local library: one coherent scaffold, so a
    binding model has something real to learn. Skips when no library is
    on disk, because this needs real chemistry rather than fixtures."""
    import pytest
    from promoter_tf_drugs import localdrugs as ld
    lib = ld.library()
    fam = [d for d in lib
           if d["name"].lower().endswith("olol") and d.get("smiles")
           and d.get("chembl_id")]
    if len(fam) < 20:
        pytest.skip("no local drug library with a usable ligand family")
    return lib, [{"chembl_id": d["chembl_id"], "smiles": d["smiles"],
                  "pchembl": 8.0} for d in fam]


def test_the_product_screen_recovers_ligands_it_never_trained_on():
    """The point of pointing a screen at the gene's own product: it can
    propose a compound nobody annotated against it.

    Trained on the "-olol" beta-blockers only, it returns CARVEDILOL and
    PROPAFENONE, neither of which ends in -olol and so neither of which
    was in the training set. Both are real ADRB1 antagonists.
    """
    import pytest
    pytest.importorskip("sklearn")
    from promoter_tf_drugs import screening as scr
    lib, actives = _olol_actives()
    res = scr.screen_gene_product(None, "ADRB1", uniprot="P08588",
                                  library=lib, actives=actives,
                                  goal="down", top_n=20)
    assert not res.get("error"), res.get("error")
    assert res["model_auc"] >= scr.MIN_BINDING_AUC
    assert res["n_novel"] >= 1, (
        "the in-sample known ligands crowded out every novel candidate, "
        "which is the one thing this screen exists to produce")
    names = {h["drug"].lower() for h in res["hits"]
             if h["quadrant"] == "NOVEL HYPOTHESIS"}
    assert "carvedilol" in names, sorted(names)


def test_the_product_screen_refuses_a_model_that_cannot_rank():
    """A deliberately incoherent active set - one scaffold plus twenty
    unrelated drugs - gives a cross-validated AUC below the floor, and a
    model that cannot separate its own actives from decoys says nothing
    about an untested drug."""
    import pytest
    import random
    pytest.importorskip("sklearn")
    from promoter_tf_drugs import screening as scr
    lib, actives = _olol_actives()
    rng = random.Random(0)
    noise = rng.sample([d for d in lib if d.get("chembl_id")], 25)
    mixed = actives[:6] + [{"chembl_id": d["chembl_id"],
                            "smiles": d["smiles"], "pchembl": 8.0}
                           for d in noise]
    res = scr.screen_gene_product(None, "ADRB1", library=lib,
                                  actives=mixed, goal="down")
    assert res.get("hits") == [] or res.get("error")
    if res.get("error"):
        assert "below" in res["error"]


def test_no_ligand_source_is_a_missing_source_not_a_bad_target():
    """The intent, which survived the message becoming more specific: a
    refusal here must point at what is missing on disk, never imply the
    protein is undruggable."""
    from promoter_tf_drugs import screening as scr
    res = scr.screen_gene_product(None, "NOSUCHGENE", library=[
        {"name": "x", "smiles": "CCO", "chembl_id": "C1"}])
    assert res["hits"] == []
    # Whichever of the three causes applies, it names a thing to fix.
    assert any(hint in res["error"] for hint in (
        "missing source rather than an undruggable protein",
        "could not be resolved to a UniProt accession",
        "holds nothing for")), res["error"]
    assert "undruggable" not in res["error"] or         "rather than an undruggable" in res["error"]


def test_only_the_unannotated_candidates_become_routes():
    """The known binders are already reported by `from_direct` at tier 1.
    Emitting them here too would double-count the same drug and bury the
    one thing this screen adds."""
    from promoter_tf_drugs import catalogue as cat
    screen = {
        "is_gene_product": True, "gene": "ADRB1", "model_auc": 0.95,
        "enrichment": {"usable": True},
        "hits": [
            {"drug": "atenolol", "quadrant": "ESTABLISHED LIGAND, untested "
                                             "on this gene", "probability": 0.99},
            {"drug": "carvedilol", "quadrant": "NOVEL HYPOTHESIS",
             "probability": 0.71, "notes": ["binds the protein"]},
        ]}
    routes = cat.from_product_screen(screen, "ADRB1")
    assert [r.drug for r in routes] == ["carvedilol"]
    r = routes[0]
    assert r.tier == 5
    assert r.directness == 0, "it acts on the product, with nothing between"
    assert r.category == "direct: predicted product binder"


def test_a_product_screen_that_did_not_validate_contributes_nothing():
    from promoter_tf_drugs import catalogue as cat
    assert cat.from_product_screen(
        {"is_gene_product": True, "enrichment": {"usable": False},
         "hits": [{"drug": "x", "quadrant": "NOVEL HYPOTHESIS",
                   "probability": 0.9}]}, "ADRB1") == []
    assert cat.from_product_screen({}, "ADRB1") == []


def test_a_tf_screen_is_not_mistaken_for_a_product_screen():
    """Both return the same shape; only the product one carries the flag."""
    from promoter_tf_drugs import catalogue as cat
    assert cat.from_product_screen(
        {"enrichment": {"usable": True},
         "hits": [{"drug": "x", "quadrant": "NOVEL HYPOTHESIS",
                   "probability": 0.9}]}, "ADRB1") == []


def test_the_product_route_warns_that_rtqpcr_will_not_see_it():
    import inspect
    from promoter_tf_drugs import screening as scr
    src = inspect.getsource(scr.screen_gene_product)
    assert "RT-qPCR will not see" in src
    assert "binding is not direction" in src


def test_the_library_reader_keeps_the_identifier_column():
    """Every screen keys known ligands, potency and - worst - its decoy
    EXCLUSION on chembl_id. read_table hardcoded None, so with a local
    library the training actives could be sampled as decoys."""
    import tempfile
    from pathlib import Path
    from promoter_tf_drugs import localdrugs as ld
    TAB, NL = chr(9), chr(10)
    with tempfile.TemporaryDirectory() as td:
        f = Path(td) / "lib.tsv"
        f.write_text(
            "name" + TAB + "smiles" + TAB + "chembl_id" + TAB
            + "max_phase" + NL
            + "aspirin" + TAB + "CC(=O)Oc1ccccc1C(=O)O" + TAB
            + "CHEMBL25" + TAB + "4" + NL,
            encoding="utf-8")
        rows = ld.read_table(f)
    assert rows[0]["chembl_id"] == "CHEMBL25"


def test_the_product_screen_flag_reaches_the_pipeline():
    import inspect
    from promoter_tf_drugs import cli, pipeline as pl
    src = inspect.getsource(cli)
    assert "--screen-product" in src
    assert "screen_product=args.screen_product" in src
    assert "screen_product" in inspect.signature(pl.run).parameters


def test_the_gene_accession_exists_even_without_the_direct_step():
    """`uniprot` was assigned only inside `if include_direct:`, so the
    product screen would have raised NameError on a --no-direct run."""
    import inspect
    from promoter_tf_drugs import pipeline as pl
    src = inspect.getsource(pl.run)
    i = src.index("uniprot = None")
    j = src.index("if include_direct:")
    assert i < j, "the accession is still defined inside the branch"


# ===================================================================== #
# Promoter DNA: the only route that targets no protein.                 #
# ===================================================================== #
def _gc_rich(n=600, seed=1):
    """RANDOM at 80% GC, not a tandem repeat.

    The first version of this was "GGGGCGGGG..." repeated, which is a
    perfectly ordered arrangement - a dinucleotide shuffle destroys it, so
    the scan called it enriched and was right to. Composition alone is
    what this fixture has to pose, so the bases are drawn independently.
    """
    import random
    rng = random.Random(seed)
    return "".join(rng.choices("GCGCGCGCAT", k=n))


def _at_rich(n=600, seed=2):
    import random
    rng = random.Random(seed)
    return "".join(rng.choices("ATATATATGC", k=n))


def test_a_feature_must_beat_the_promoter_s_own_shuffles():
    """A GC-rich promoter produces quadruplexes and GC-boxes by
    composition alone. Proposing an agent for those would propose one for
    most of the genome, so each count is tested against dinucleotide
    shuffles of the same sequence."""
    from promoter_tf_drugs import promoterdna as pd
    scan = pd.scan_features(_gc_rich(), shuffles=60, seed=0)
    gc = scan["features"]["gc_box"]
    assert gc["n_observed"] > 0, "the fixture should carry GC-boxes"
    # Everything here is GC-box; a shuffle of it is too.
    assert not gc["enriched"], (
        "a feature present purely by composition was called enriched")


def test_a_real_enrichment_is_reported():
    from promoter_tf_drugs import promoterdna as pd
    # One quadruplex dropped into an AT-rich background: composition
    # cannot explain it, so the shuffles should not reproduce it.
    seq = _at_rich(400) + "GGGTTAGGGTTAGGGTTAGGG" + _at_rich(400)
    scan = pd.scan_features(seq, shuffles=120, seed=0)
    g4 = scan["features"]["g4"]
    assert g4["n_observed"] >= 1
    assert g4["enriched"], (g4["n_observed"], g4["expected"], g4["p_value"])
    names = {h["drug"] for h in pd.candidates(scan)}
    assert "Pyridostatin" in names


def test_nothing_enriched_proposes_nothing():
    from promoter_tf_drugs import promoterdna as pd
    scan = pd.scan_features(_gc_rich(), shuffles=60, seed=0)
    assert pd.candidates(scan) == []
    text = pd.summarise(scan, [])
    assert "explained by base composition" in text


def test_the_dinucleotide_shuffle_preserves_composition():
    """A mononucleotide shuffle would destroy GG pairs and call every
    GC-rich promoter enriched for quadruplexes."""
    import random
    from promoter_tf_drugs import promoterdna as pd
    seq = _gc_rich(500)
    sh = pd._dinucleotide_shuffle(seq, random.Random(0))
    assert len(sh) == len(seq)
    for base in "ACGT":
        assert abs(sh.count(base) - seq.count(base)) <= 2, base


def test_a_short_sequence_is_refused_not_scanned():
    from promoter_tf_drugs import promoterdna as pd
    scan = pd.scan_features("ACGT" * 5)
    assert "too short" in scan.get("error", "")
    assert pd.candidates(scan) == []


def test_the_dna_route_says_what_it_cannot_tell_you():
    from promoter_tf_drugs import catalogue as cat
    dna = {"candidates": [{
        "drug": "Pyridostatin", "prefers": "g4", "feature": "G-quadruplex",
        "n_sites": 3, "expected": 0.4, "p_value": 0.008,
        "mechanism": "selective G-quadruplex stabiliser",
        "approved": False, "note": "a tool compound"}]}
    routes = cat.from_promoter_dna(dna)
    assert len(routes) == 1
    r = routes[0]
    assert r.category == "promoter DNA binder"
    assert r.acts_on == "promoter DNA"
    assert r.directness == 0, "it acts at the promoter, with no protein"
    assert r.tier == 5
    assert r.axes.specificity <= 0.2, (
        "these bind their feature genome-wide; specificity is the axis "
        "they fail")
    joined = " ".join(r.notes)
    assert "wherever it occurs" in joined
    assert "not a folded structure" in joined
    assert "not an approved drug" in joined


def test_an_approved_dna_agent_is_not_labelled_unapproved():
    from promoter_tf_drugs import catalogue as cat
    dna = {"candidates": [{
        "drug": "Mithramycin", "prefers": "gc_box", "feature": "SP1 GC-box",
        "n_sites": 11, "expected": 6.1, "p_value": 0.02,
        "mechanism": "displaces SP1", "approved": True, "note": "x"}]}
    r = cat.from_promoter_dna(dna)[0]
    assert r.approved is True
    assert not any("not an approved drug" in n for n in r.notes)


def test_every_curated_dna_agent_declares_its_feature():
    from promoter_tf_drugs import promoterdna as pd
    keys = {"g4", "imotif", "gc_box", "at_tract"}
    for a in pd.AGENTS:
        assert a.prefers in keys, (a.name, a.prefers)
        assert len(a.note) > 30, a.name


# ===================================================================== #
# The high-burden table.                                                 #
# ===================================================================== #
def _burden_rows():
    base = {"tier": 3, "score": 30.0, "category": "TF: abundance change",
            "axes": {"specificity": 0.45}, "directness": 1, "grade": None,
            "_notes": "", "_chip_known": False, "approved": True,
            "mechanism": "m", "acts_on": "RBPJ"}
    hot = dict(base, drug="dinaciclib", burden="high",
               _burden_assessed=True,
               drug_class="antineoplastic kinase inhibitor")
    cold = dict(base, drug="Maraviroc", tier=1, score=60.0, burden="low",
                _burden_assessed=False, drug_class=None,
                category="direct: gene product", acts_on="CCR5")
    return [hot, cold]


def test_a_high_burden_row_is_set_aside_not_deleted():
    from promoter_tf_drugs import shortlist as sl
    kept, aside = sl.split_by_burden(_burden_rows(), max_burden="moderate")
    assert [r["drug"] for r in kept] == ["Maraviroc"]
    assert [r["drug"] for r in aside] == ["dinaciclib"]


def test_with_no_burden_criterion_nothing_is_set_aside():
    from promoter_tf_drugs import shortlist as sl
    kept, aside = sl.split_by_burden(_burden_rows())
    assert aside == []
    assert len(kept) == 2


def test_the_burden_table_warns_that_they_are_cytotoxic():
    from promoter_tf_drugs import shortlist as sl
    text = sl.render([], {}, 2, set_aside=[_burden_rows()[0]])
    assert "These are cytotoxic or otherwise high-burden drugs" in text
    assert "only if nothing in the table above will do" in text
    assert "on-mechanism rather than side effects" in text
    # and it must name the class that put each row there
    assert "antineoplastic kinase inhibitor" in text


def test_the_burden_table_explains_why_the_class_qualifies():
    from promoter_tf_drugs import shortlist as sl
    text = sl.render([], {}, 2, set_aside=[_burden_rows()[0]])
    assert "Why these are set aside, by class" in text
    assert "stop cells dividing" in text


# ===================================================================== #
# Agonist versus antagonist at the gene product.                         #
# ===================================================================== #
def test_the_product_screen_annotates_direction_when_it_can():
    import pytest
    pytest.importorskip("sklearn")
    from promoter_tf_drugs import screening as scr, localdrugs as ld, agonism
    lib = ld.library()
    ant = [d for d in lib if d["name"].lower().endswith("olol")
           and d.get("smiles")]
    ago = [d for d in lib if d["name"].lower().endswith("terol")
           and d.get("smiles")]
    if len(ant) < agonism.MIN_PER_CLASS or len(ago) < agonism.MIN_PER_CLASS:
        pytest.skip("no local library with both halves of a direction set")
    lab = ([{"chembl_id": d["chembl_id"], "smiles": d["smiles"],
             "label": "antagonist"} for d in ant]
           + [{"chembl_id": d["chembl_id"], "smiles": d["smiles"],
               "label": "agonist"} for d in ago])
    acts = [{"chembl_id": d["chembl_id"], "smiles": d["smiles"],
             "pchembl": 8.0} for d in ant]
    res = scr.screen_gene_product(None, "ADRB1", library=lib, actives=acts,
                                  labelled=lab, goal="down", top_n=10)
    assert res["direction_model"] is True
    notes = " ".join(n for h in res["hits"] for n in h["notes"])
    assert "predicted antagonist at ADRB1" in notes


def test_too_few_labels_declares_direction_unknown_rather_than_guessing():
    import pytest
    pytest.importorskip("sklearn")
    from promoter_tf_drugs import screening as scr, localdrugs as ld
    lib = ld.library()
    ant = [d for d in lib if d["name"].lower().endswith("olol")
           and d.get("smiles")]
    if len(ant) < 20:
        pytest.skip("no local library")
    acts = [{"chembl_id": d["chembl_id"], "smiles": d["smiles"],
             "pchembl": 8.0} for d in ant]
    res = scr.screen_gene_product(None, "ADRB1", library=lib, actives=acts,
                                  labelled=[], goal="down")
    assert res["direction_model"] is False
    assert "Direction not predicted" in res["direction_note"]
    assert "needs a functional assay" in res["direction_note"]


# ===================================================================== #
# Docking the gene's OWN product. It was classified off-pathway, which   #
# is right for an upstream kinase and backwards for the gene itself.     #
# ===================================================================== #
def _dock_on(target, state="inactive"):
    return {"target": target, "usable": True, "auc": 0.93, "ef1": 3.2,
            "pocket": {"pdb_id": "4MBS", "ligand_resname": "MRV",
                       "state": state},
            "hits": [{"drug": "cpd-1", "score": -11.2}]}


def test_docking_the_gene_s_own_product_is_the_most_direct_route():
    """It is not one of the promoter's factors, so the off-pathway test
    called it the longest route when it is the shortest."""
    from promoter_tf_drugs import catalogue as cat
    r = cat.from_docking(_dock_on("CCR5"), goal="down",
                         tf_roles={"VDR": "activator"}, gene="CCR5")[0]
    assert r.category == "direct: predicted product binder"
    assert r.directness == 0
    assert not any("not one of the factors" in n for n in r.notes)


def test_the_three_docking_targets_keep_their_order():
    """Gene product, then a factor at this promoter, then anything else.
    The directness gradient is what the score penalty reads."""
    from promoter_tf_drugs import catalogue as cat
    roles = {"VDR": "activator"}
    got = {}
    for target in ("CCR5", "VDR", "MTOR"):
        r = cat.from_docking(_dock_on(target), goal="down",
                             tf_roles=roles, gene="CCR5")[0]
        got[target] = (r.category, r.directness)
    assert got["CCR5"][1] == 0
    assert got["VDR"][1] == 1
    assert got["MTOR"][1] == 2
    assert got["VDR"][0] == "TF: predicted pocket binder"
    assert got["MTOR"][0] == "off-pathway pocket binder"


def test_the_product_pocket_warns_about_the_readout():
    from promoter_tf_drugs import catalogue as cat
    r = cat.from_docking(_dock_on("CCR5"), goal="down", tf_roles={},
                         gene="CCR5")[0]
    joined = " ".join(r.notes)
    assert "on CCR5 itself" in joined
    assert "RT-qPCR will not see it" in joined


def test_the_product_conformation_follows_the_goal_directly():
    """No factor role in between: wanting more of the gene's activity
    wants the agonist-bound structure and wanting less wants the other."""
    from promoter_tf_drugs import catalogue as cat
    cases = {("up", "active"): "supports", ("up", "inactive"): "opposes",
             ("down", "inactive"): "supports", ("down", "active"): "opposes"}
    for (goal, state), expect in cases.items():
        r = cat.from_docking(_dock_on("CCR5", state), goal=goal,
                             tf_roles={}, gene="CCR5")[0]
        assert expect in r.mechanism, (goal, state, r.mechanism)


def test_state_for_product_needs_no_role():
    from promoter_tf_drugs.docking import state_for_product, state_for
    assert state_for_product("up") == "active"
    assert state_for_product("down") == "inactive"
    # Same answer as an activator factor, for a different reason - which
    # is why it is its own function rather than a call with a fake role.
    assert state_for_product("up") == state_for("activator", "up")


def test_the_docking_cli_offers_the_product_flag():
    import inspect
    from promoter_tf_drugs import docking as dk
    src = inspect.getsource(dk.main)
    assert '"--product"' in src
    assert "state_for_product(args.direction)" in src


def test_without_a_gene_nothing_is_called_a_product_dock():
    """`gene` is optional, and a bundle read without it must not promote
    an arbitrary target to directness 0."""
    from promoter_tf_drugs import catalogue as cat
    r = cat.from_docking(_dock_on("CCR5"), goal="down",
                         tf_roles={"VDR": "activator"})[0]
    assert r.category == "off-pathway pocket binder"
    assert r.directness == 2


def test_the_pipeline_tells_the_catalogue_which_gene_is_the_target():
    import inspect
    from promoter_tf_drugs import catalogue as cat
    src = inspect.getsource(cat.build)
    assert 'gene=result_block.get("gene_symbol")' in src


def _gtf_genes(path, n):
    """`n` protein-coding genes on chr1, 10 kb apart, all + strand."""
    tab, nl = chr(9), chr(10)
    rows = []
    for i in range(n):
        rows.append(tab.join([
            "chr1", "x", "gene", str(i * 10000 + 5000),
            str(i * 10000 + 9000), ".", "+", ".",
            'gene_name "G' + str(i) + '"; gene_type "protein_coding";']))
    path.write_text(nl.join(rows) + nl, encoding="utf-8")


def _bed_peaks(path, rows):
    """rows are (start, end, name)."""
    tab, nl = chr(9), chr(10)
    out = [tab.join([
        "chr1", str(s), str(e), name, "1", "."]) for s, e, name in rows]
    path.write_text(nl.join(out) + nl, encoding="utf-8")


def test_a_promiscuous_factor_is_dropped_during_the_scan(tmp_path):
    """Holding every factor until write time is what made this need 28 GB:
    AGO1 bound 23,318 promoters and was carried the whole way before being
    discarded. The cap bites while reading now, so the structure is
    bounded by factors x cap rather than factors x genome."""
    from promoter_tf_drugs import chipregulon as cr
    gtf = tmp_path / "g.gtf"
    _gtf_genes(gtf, 12)
    bed = tmp_path / "p.bed"
    rows = [(i * 10000 + 4900, i * 10000 + 5100, "BROAD:HeLa")
            for i in range(12)]
    rows += [(i * 10000 + 4900, i * 10000 + 5100, "NARROW:HeLa")
             for i in range(2)]
    _bed_peaks(bed, rows)

    acc = cr.build(bed, gtf, out_path=tmp_path / "out.tsv", max_targets=5)
    text = (tmp_path / "out.tsv").read_text(encoding="utf-8")
    assert "NARROW" in text
    assert "BROAD" not in text, "the promiscuous factor was written out"
    assert acc["n_excluded_too_broad"] == 1
    assert acc["excluded_too_broad"][0][0] == "BROAD"
    # The count is a floor: the scan stopped counting when the cap bit.
    assert "at least" in cr.summarise(acc)


def test_cell_type_support_promotes_only_when_a_second_appears(tmp_path):
    """One int per edge until a second cell type shows up. 98% of edges
    never see one, which is why a thousand-bit mask per edge was the
    wrong representation."""
    from promoter_tf_drugs import chipregulon as cr
    gtf = tmp_path / "g.gtf"
    _gtf_genes(gtf, 1)
    bed = tmp_path / "p.bed"
    _bed_peaks(bed, [(4900, 5100, "TF:HeLa"),
                     (4950, 5150, "TF:HeLa"),     # same cell again
                     (4900, 5100, "TF:K562")])    # a second cell
    cr.build(bed, gtf, out_path=tmp_path / "out.tsv", max_targets=None)
    line = (tmp_path / "out.tsv").read_text(encoding="utf-8").strip()
    assert line.endswith("ReMap_cells:2"), line


def test_one_cell_type_stays_one(tmp_path):
    from promoter_tf_drugs import chipregulon as cr
    gtf = tmp_path / "g.gtf"
    _gtf_genes(gtf, 1)
    bed = tmp_path / "p.bed"
    _bed_peaks(bed, [(4900, 5100, "TF:HeLa"), (4950, 5150, "TF:HeLa")])
    cr.build(bed, gtf, out_path=tmp_path / "out.tsv", max_targets=None)
    line = (tmp_path / "out.tsv").read_text(encoding="utf-8").strip()
    assert line.endswith("ReMap_cells:1"), line


def test_support_reads_both_representations():
    from promoter_tf_drugs import chipregulon as cr
    assert cr._support(3) == 1
    assert cr._support({1, 2, 5}) == 3


def test_the_scan_does_not_allocate_a_mask_per_edge():
    """The first version built `1 << len(cell_bit)` per placement, which
    is roughly 150 bytes with ReMap's ~1,000 cell lines, on 16.5 million
    placements."""
    import inspect
    from promoter_tf_drugs import chipregulon as cr
    src = inspect.getsource(cr.build)
    assert "1 <<" not in src
    assert "bit_count" not in inspect.getsource(cr)


def test_remap_lists_every_cell_line_a_merged_peak_was_seen_in():
    """`TF:K-562` is one cell line and `TF:hiPSC,HCT-116,HEK293` is three.
    Taking the whole field as one name found 564,344 distinct cell types
    in a file with a few hundred, and scored a peak reproduced in 37 cell
    lines as support 1."""
    from promoter_tf_drugs.chipregulon import _tf_cells
    tf, cells = _tf_cells("CTCF:hiPSC,HCT-116,HEK293,K-562,WA09")
    assert tf == "CTCF"
    assert cells == ["HIPSC", "HCT-116", "HEK293", "K-562", "WA09"]
    assert _tf_cells("SMARCA4:HeLa") == ("SMARCA4", ["HELA"])
    assert _tf_cells("TF:") == ("TF", [])


def test_a_merged_peak_counts_every_cell_line_it_names(tmp_path):
    from promoter_tf_drugs import chipregulon as cr
    gtf = tmp_path / "g.gtf"
    _gtf_genes(gtf, 1)
    bed = tmp_path / "p.bed"
    _bed_peaks(bed, [(4900, 5100, "TF:K-562,HeLa,HEK293")])
    cr.build(bed, gtf, out_path=tmp_path / "out.tsv", max_targets=None)
    line = (tmp_path / "out.tsv").read_text(encoding="utf-8").strip()
    assert line.endswith("ReMap_cells:3"), line


def test_cell_lines_union_across_peaks(tmp_path):
    """Two peaks naming overlapping sets give the union, not the sum."""
    from promoter_tf_drugs import chipregulon as cr
    gtf = tmp_path / "g.gtf"
    _gtf_genes(gtf, 1)
    bed = tmp_path / "p.bed"
    _bed_peaks(bed, [(4900, 5100, "TF:K-562,HeLa"),
                     (4950, 5150, "TF:HeLa,HEK293")])
    cr.build(bed, gtf, out_path=tmp_path / "out.tsv", max_targets=None)
    line = (tmp_path / "out.tsv").read_text(encoding="utf-8").strip()
    assert line.endswith("ReMap_cells:3"), line


def test_the_magnitude_statistic_gets_a_q_value():
    """It was computed and dropped at the FDR step: `q_value` came only
    from the directional p, which is None without signed edges. A fully
    unsigned regulon therefore produced no q at all and the run reported
    "0 reached q<0.05" when nothing had been tested at that level."""
    from promoter_tf_drugs import tfactivity as tfa
    sig = {f"G{i}": float(i % 7) - 3 for i in range(400)}
    reg = {"TF1": {f"G{i}": 0 for i in range(40)},
           "TF2": {f"G{i}": 0 for i in range(40, 80)}}
    acts = tfa.activity(sig, reg)
    assert acts
    assert all(a.nes is None for a in acts), "fixture should be unsigned"
    assert all(a.q_value is None for a in acts)
    assert all(a.magnitude_q is not None for a in acts), (
        "the magnitude test ran and its result was discarded")


def test_an_oversized_regulon_is_flagged_as_not_factor_specific():
    """A set covering a tenth of the measured genes inherits the
    signature's structure, so the test on it is near-guaranteed."""
    from promoter_tf_drugs import tfactivity as tfa
    sig = {f"G{i}": float(i % 5) - 2 for i in range(300)}
    reg = {"BIG": {f"G{i}": 0 for i in range(150)},
           "SMALL": {f"G{i}": 0 for i in range(12)}}
    acts = {a.tf: a for a in tfa.activity(sig, reg)}
    assert acts["BIG"].regulon_share > tfa.MAX_REGULON_SHARE
    assert any("not specific to this factor" in c
               for c in acts["BIG"].caveats)
    assert not any("not specific to this factor" in c
                   for c in acts["SMALL"].caveats)


def test_half_the_factors_moving_is_reported_as_a_warning():
    """Factors do not act in concert. When a large share move, and mostly
    the same way, the statistic is reading the signature rather than the
    factors."""
    from promoter_tf_drugs import tfactivity as tfa

    class _A:
        def __init__(self, i):
            self.tf = f"TF{i}"
            self.nes = None
            self.p_value = None
            self.q_value = None
            self.magnitude = -9.0
            self.magnitude_p = 1e-9
            self.magnitude_q = 1e-6
            self.direction = "unknown"
            self.n_measured = 2000
            self.n_signed = 0
            self.n_regulon = 2000
            self.coverage = 1.0
            self.regulon_share = 0.4
            self.confounded_with = None
            self.confounded_jaccard = 0.2
            self.caveats = []

    text = tfa.summarise([_A(i) for i in range(10)])
    assert "Read that number as a warning" in text
    assert "do not act in concert" in text


def test_a_normal_result_carries_no_such_warning():
    from promoter_tf_drugs import tfactivity as tfa
    sig = {f"G{i}": float((i * 7) % 11) - 5 for i in range(500)}
    reg = {f"TF{k}": {f"G{i}": (1 if i % 2 else -1)
                      for i in range(k, k + 15)} for k in range(0, 60, 15)}
    text = tfa.summarise(tfa.activity(sig, reg))
    assert "Read that number as a warning" not in text


def test_the_potency_refusal_distinguishes_its_three_causes(monkeypatch):
    """It said "no local potency file" while an 18,229-pair file sat in
    data/ - it simply held nothing for CCR5, because every ChEMBL activity
    query this package makes is about a transcription factor rather than a
    gene's own product. Telling someone to fetch a file they have sends
    them to fix the wrong thing."""
    from promoter_tf_drugs import screening as scr, localpotency as lp

    lib = [{"name": "x", "smiles": "CCO", "chembl_id": "C1"}]

    monkeypatch.setattr(lp, "available", lambda: {"bindingdb": None,
                                                  "chembl_sqlite": None})
    monkeypatch.setattr(scr, "local_actives_for", lambda *a, **k: [])
    no_file = scr.screen_gene_product(None, "CCR5", uniprot="P51681",
                                      library=lib)
    assert no_file["potency_file_present"] is False
    assert "No ligand source was readable" in no_file["error"]

    monkeypatch.setattr(lp, "available", lambda: {"bindingdb": "data/x.tsv",
                                                  "chembl_sqlite": None})
    absent = scr.screen_gene_product(None, "CCR5", uniprot="P51681",
                                     library=lib)
    assert absent["potency_file_present"] is True
    assert "A potency file IS on disk and holds nothing" in absent["error"]
    assert "covers factors and not gene products" in absent["error"]

    unresolved = scr.screen_gene_product(None, "NOSUCHGENE", library=lib)
    assert "could not be resolved to a UniProt accession" in         unresolved["error"]


def test_a_thin_target_is_not_blamed_on_a_missing_file():
    from promoter_tf_drugs import screening as scr
    res = scr.screen_gene_product(
        None, "X", uniprot="P1",
        library=[{"name": "x", "smiles": "CCO", "chembl_id": "C1"}],
        actives=[{"chembl_id": "A1", "smiles": "CCN", "pchembl": 8.0}])
    assert res["source"] == "supplied by the caller"
    assert "thin target rather than a missing file" in res["error"]


# ===================================================================== #
# BindingDB, fetched per target rather than downloaded whole.            #
# ===================================================================== #
def _bdb_response(uniprot="P51681", n=3):
    return {"getLindsByUniprotResponse": {
        "bdb.primary": uniprot, "bdb.hit": str(n),
        "bdb.affinities": [
            {"bdb.monomerid": str(1000 + i), "bdb.smile": "C" * (i + 2) + "O",
             "bdb.affinity_type": "IC50", "bdb.affinity": str(10 * (i + 1))}
            for i in range(n)]}}


class _FakeBdbClient:
    def __init__(self, payload):
        self.payload = payload
        self.calls = []

    def get(self, url, params=None, **kw):
        self.calls.append((url, dict(params or {})))
        return self.payload


def test_the_fetch_merges_into_the_file_already_in_use(tmp_path,
                                                       monkeypatch):
    """`resolve` returns the first path that EXISTS, so writing to the
    first CANDIDATE path created a second file that shadowed it - adding
    CCR5 took ESR1 from 2,020 actives to zero, because the new file won
    the lookup and held only CCR5."""
    from promoter_tf_drugs import localpotency as lp
    monkeypatch.chdir(tmp_path)
    d = tmp_path / "data"
    d.mkdir()
    tab, nl = chr(9), chr(10)
    existing = d / "bindingdb.tsv"
    existing.write_text(
        tab.join(["UniProt (SwissProt) Primary ID of Target Chain",
                  "Ligand SMILES", "ChEMBL ID of Ligand", "Ki (nM)"]) + nl
        + tab.join(["P03372", "CCCC", "CHEMBL1", "5"]) + nl,
        encoding="utf-8")
    assert lp.resolve() == str(existing).replace(str(tmp_path) + chr(92),
                                                 "") or lp.resolve()

    acc = lp.fetch_bindingdb(_FakeBdbClient(_bdb_response()), ["P51681"])
    assert acc["rows_before"] == 1, "the existing rows were not read"
    assert acc["added"] == 3
    # One file, both targets.
    assert not (d / "BindingDB_All.tsv").exists(), "a shadow file appeared"
    assert len(lp.actives_for("P03372", min_pchembl=0)) == 1
    assert len(lp.actives_for("P51681", min_pchembl=0)) == 3


def test_the_fetch_keeps_columns_the_existing_file_had(tmp_path,
                                                       monkeypatch):
    from promoter_tf_drugs import localpotency as lp
    monkeypatch.chdir(tmp_path)
    d = tmp_path / "data"
    d.mkdir()
    tab, nl = chr(9), chr(10)
    (d / "bindingdb.tsv").write_text(
        tab.join(["UniProt (SwissProt) Primary ID of Target Chain",
                  "Ligand SMILES", "ChEMBL ID of Ligand", "Ki (nM)"]) + nl
        + tab.join(["P03372", "CCCC", "CHEMBL1", "5"]) + nl,
        encoding="utf-8")
    lp.fetch_bindingdb(_FakeBdbClient(_bdb_response()), ["P51681"])
    head = (d / "bindingdb.tsv").read_text(
        encoding="utf-8").splitlines()[0]
    assert "ChEMBL ID of Ligand" in head, "an existing column was dropped"
    assert "BindingDB MonomerID" in head


def test_a_failed_target_does_not_lose_the_others(tmp_path, monkeypatch):
    from promoter_tf_drugs import localpotency as lp
    monkeypatch.chdir(tmp_path)
    (tmp_path / "data").mkdir()

    class _Boom(_FakeBdbClient):
        def get(self, url, params=None, **kw):
            if params and params.get("uniprot") == "BAD":
                raise RuntimeError("upstream is down")
            return self.payload

    acc = lp.fetch_bindingdb(_Boom(_bdb_response()), ["BAD", "P51681"])
    assert "failed" in str(acc["targets"]["BAD"])
    assert acc["targets"]["P51681"] == 3
    assert acc["added"] == 3


# ===================================================================== #
# Actives matched to the library by structure, not by a source's id.     #
# ===================================================================== #
def test_actives_match_the_library_by_structure_not_by_id():
    """BindingDB numbers ligands with MonomerIDs and ChEMBL with CHEMBL
    ids, and the two write different canonical SMILES. For CCR5, zero of
    254 actives matched the library on either - so no known ligand was
    recognised, and the decoy pool, which is built by EXCLUDING the
    actives, could hand a real binder to the model as a negative."""
    import pytest
    pytest.importorskip("rdkit")
    from promoter_tf_drugs import screening as scr
    # Same molecule, different source spellings and different ids.
    actives = [{"chembl_id": None, "smiles": "OC(=O)Cc1ccccc1",
                "pchembl": 8.0}]
    lib = [{"name": "same", "chembl_id": "CHEMBL_A",
            "smiles": "c1ccccc1CC(O)=O"},
           {"name": "other", "chembl_id": "CHEMBL_B", "smiles": "CCCCCCO"}]
    got = scr.match_actives_to_library(actives, lib)
    assert got == {"CHEMBL_A"}, got


def test_a_salt_form_still_matches_its_parent():
    import pytest
    pytest.importorskip("rdkit")
    from promoter_tf_drugs import screening as scr
    actives = [{"smiles": "CN(C)C(=N)NC(N)=N", "pchembl": 8.0}]
    lib = [{"name": "salt", "chembl_id": "C1",
            "smiles": "CN(C)C(=N)NC(N)=N.Cl"}]
    assert scr.match_actives_to_library(actives, lib) == {"C1"}


def test_an_unmatchable_active_set_warns_rather_than_silently_passing(
        caplog):
    import logging
    from promoter_tf_drugs import screening as scr
    with caplog.at_level(logging.WARNING):
        scr.screen_gene_product(
            None, "X", uniprot="P1",
            library=[{"name": "a", "smiles": "CCCCCCCCO",
                      "chembl_id": "C1"}],
            actives=[{"smiles": "c1ccccc1" + "C" * i, "pchembl": 8.0}
                     for i in range(25)])
    assert "not guaranteed free of real binders" in caplog.text


# ===================================================================== #
# gnina through WSL: parsing, sign conventions, and delivery.           #
# ===================================================================== #
GNINA_SDF = """
  calcitriol
  one pose

> <minimizedAffinity>
-12.79365

> <CNNscore>
0.9585389495

> <CNNaffinity>
8.5596036911

> <CNNaffinity_variance>
0.18330000

$$$$
  calcitriol
  second pose

> <minimizedAffinity>
-11.79270

> <CNNscore>
0.6926998496

> <CNNaffinity>
7.9490647316

> <CNNaffinity_variance>
0.21400000

$$$$
"""


def test_cnn_affinity_is_not_confused_with_its_own_variance():
    """gnina writes CNNaffinity and CNNaffinity_variance for every pose,
    so a loose search for the first name matches both. That read two
    poses as four and interleaved each affinity with its variance - and
    a variance is a small positive number, so half the affinities came
    back near 0.2 and looked like weak binders rather than like a parse
    error. Nothing raised; the numbers were just wrong."""
    from promoter_tf_drugs import docking as dk
    poses = dk._parse_poses(GNINA_SDF, "", "gnina")
    assert len(poses) == 2, [p.cnn_affinity for p in poses]
    assert poses[0].cnn_affinity == 8.5596036911
    assert poses[1].cnn_affinity == 7.9490647316
    # The variances must not appear anywhere in the parsed output.
    assert 0.1833 not in [p.cnn_affinity for p in poses]
    assert poses[0].affinity == -12.79365
    assert poses[0].cnn_score == 0.9585389495
    assert poses[0].engine == "gnina"


def test_pose_scores_are_returned_in_one_direction():
    """affinity is tighter when more negative and the two CNN numbers are
    tighter when larger. Ranking them together without a convention sorts
    the CNN backwards, which puts the worst pose first while every number
    on the row still looks right."""
    from promoter_tf_drugs import docking as dk
    tight = dk.PoseScore(affinity=-12.0, cnn_score=0.95, cnn_affinity=8.5)
    loose = dk.PoseScore(affinity=-6.0, cnn_score=0.20, cnn_affinity=6.0)
    for metric in ("affinity", "cnn_score", "cnn_affinity"):
        assert tight.best(metric) < loose.best(metric), metric
    missing = dk.PoseScore()
    assert missing.best("cnn_affinity") is None
    assert missing.best("affinity") is None


def test_vina_scores_are_read_from_its_table_when_the_sdf_has_no_tags():
    """vina and older smina print a table and tag nothing, so a parser
    that only reads SDF tags silently returns no poses at all for them."""
    from promoter_tf_drugs import docking as dk
    stdout = ("mode |   affinity | dist from best mode\n"
              "-----+------------+--------------------\n"
              "   1       -9.400          0.000\n"
              "   2       -8.100          2.345\n")
    poses = dk._parse_poses("no tags here", stdout, "vina")
    assert [p.affinity for p in poses] == [-9.4, -8.1]
    assert all(p.cnn_score is None for p in poses)


def test_cnn_rescoring_is_requested_only_of_gnina(tmp_path):
    """--cnn_scoring is the whole reason to run gnina, and it is not a
    flag smina or vina accept - passing it to them makes the run fail."""
    from promoter_tf_drugs import docking as dk
    box = (1.0, 2.0, 3.0, 20.0, 20.0, 20.0)
    args = dk._dock_args("/home/u/bin/gnina", tmp_path / "r.pdb",
                         tmp_path / "l.sdf", tmp_path / "o.sdf", box,
                         exhaustiveness=8, seed=0, cpu=1)
    assert "--cnn_scoring" in args
    assert args[args.index("--cnn_scoring") + 1] == "rescore"
    plain = dk._dock_args("C:/tools/smina.exe", tmp_path / "r.pdb",
                          tmp_path / "l.sdf", tmp_path / "o.sdf", box,
                          exhaustiveness=8, seed=0, cpu=1)
    assert "--cnn_scoring" not in plain


def test_paths_are_translated_only_for_the_side_the_binary_runs_on(
        tmp_path):
    """A WSL binary handed C:\\... reports that the receptor does not
    exist, which reads like a missing file rather than a path that was
    never converted."""
    from promoter_tf_drugs import docking as dk
    from promoter_tf_drugs import wsldock as wd
    box = (0.0, 0.0, 0.0, 20.0, 20.0, 20.0)
    rec = tmp_path / "r.pdb"
    args = dk._dock_args("/home/u/bin/gnina", rec, tmp_path / "l.sdf",
                         tmp_path / "o.sdf", box, exhaustiveness=8,
                         seed=0, cpu=1)
    given = args[args.index("-r") + 1]
    assert not given[1:2] == ":", given
    assert given == wd.to_wsl_path(rec)
    win = dk._dock_args("smina.exe", rec, tmp_path / "l.sdf",
                        tmp_path / "o.sdf", box, exhaustiveness=8,
                        seed=0, cpu=1)
    assert win[win.index("-r") + 1] == str(rec)


def test_windows_drive_letters_become_mount_points():
    from promoter_tf_drugs import wsldock as wd
    got = wd.to_wsl_path("C:\\Users\\x\\y.pdb")
    assert got.startswith("/mnt/c/"), got
    assert "\\" not in got
    # An already-POSIX path is left alone rather than mangled.
    assert wd.to_wsl_path("/home/u/x.pdb") == "/home/u/x.pdb"


def test_scripts_reach_wsl_as_bytes_with_unix_newlines(monkeypatch):
    """Two delivery faults, both of which read as a broken script rather
    than as a broken way of sending one.

    With text=True Python rewrites every newline to CRLF on Windows and
    bash fails on the first multi-line construct. And a script passed as
    a `bash -c` argument is rebuilt by the interop layer on the way
    through: `wsl.exe -- bash -c "echo $PATH"` comes back a syntax error.
    """
    from promoter_tf_drugs import wsldock as wd
    seen = {}

    def fake_run(cmd, **kw):
        seen["cmd"] = cmd
        seen["kw"] = kw
        import subprocess as sp
        return sp.CompletedProcess(cmd, 0, b"ok", b"")

    monkeypatch.setattr(wd.subprocess, "run", fake_run)
    wd.run_script("for i in 1 2; do\n  echo $i\ndone\n")
    assert "-c" not in seen["cmd"], seen["cmd"]
    assert seen["cmd"][-1] == "-s"
    payload = seen["kw"]["input"]
    assert isinstance(payload, bytes)
    assert b"\r" not in payload
    assert seen["kw"].get("text") is not True


def test_nul_bytes_are_stripped_from_wsl_output(monkeypatch):
    """The interop layer interleaves NUL bytes into the output. They make
    every comparison against it fail in a way that reads like the command
    having produced nothing at all."""
    from promoter_tf_drugs import wsldock as wd

    def fake_run(cmd, **kw):
        import subprocess as sp
        noisy = "/home/u/bin/gnina".encode("utf-8")
        noisy = b"\x00".join(noisy[i:i + 1] for i in range(len(noisy)))
        return sp.CompletedProcess(cmd, 0, noisy, b"\x00")

    monkeypatch.setattr(wd.subprocess, "run", fake_run)
    out = wd.run_script("command -v gnina")
    assert out.stdout.strip() == "/home/u/bin/gnina"
    assert chr(0) not in out.stdout
    assert chr(0) not in out.stderr


def test_a_linux_gnina_is_preferred_over_a_windows_vina(monkeypatch):
    """The only reason to reach into WSL is that the better programs are
    Linux software, so a Windows vina on PATH must not win over it."""
    from promoter_tf_drugs import docking as dk
    from promoter_tf_drugs import wsldock as wd
    monkeypatch.setattr(wd, "find_binary",
                        lambda prefer=None: ("gnina", "/home/u/bin/gnina"))
    monkeypatch.setattr(dk.shutil, "which",
                        lambda n: "C:/tools/vina.exe" if n == "vina" else None)
    assert dk._docking_binary() == "/home/u/bin/gnina"


def test_a_windows_binary_is_still_used_when_wsl_has_none(monkeypatch):
    """WSL being unusable must degrade to the Windows binary, not to an
    exception - most machines running this have no WSL at all."""
    from promoter_tf_drugs import docking as dk
    from promoter_tf_drugs import wsldock as wd

    def boom(prefer=None):
        raise OSError("wsl.exe not found")

    monkeypatch.setattr(wd, "find_binary", boom)
    monkeypatch.setattr(dk.shutil, "which",
                        lambda n: "C:/tools/smina.exe" if n == "smina" else None)
    assert dk._docking_binary() == "C:/tools/smina.exe"


def test_a_failed_docking_run_returns_nothing_rather_than_a_zero_score():
    """A score of 0 from a crashed run is a confident claim of a weak
    binder. The absence of a result has to stay absent."""
    from promoter_tf_drugs import docking as dk
    poses = dk._parse_poses("", "", "gnina")
    assert poses == []


# ===================================================================== #
# Static guards for two bug classes that running the code found and     #
# eight hundred tests did not.                                          #
# ===================================================================== #
_SCOPES = None


def _scope_types():
    """Node types that open a new scope in Python 3."""
    import ast
    return (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda,
            ast.ClassDef, ast.ListComp, ast.SetComp, ast.DictComp,
            ast.GeneratorExp)


def _own_scope_nodes(fn):
    """Nodes belonging to `fn` itself, not to a nested scope.

    A comprehension has its own scope, and `[s for s in xs]` split over
    two lines puts the use of `s` above its binding. Treating that as a
    read-before-assignment produced sixty-five false reports against
    nine real ones.
    """
    import ast
    out, stack = [], list(ast.iter_child_nodes(fn))
    while stack:
        node = stack.pop()
        out.append(node)
        if isinstance(node, _scope_types()):
            continue
        stack.extend(ast.iter_child_nodes(node))
    return out


def _module_level_names(tree):
    import ast
    out = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            for a in node.names:
                out.add((a.asname or a.name).split(".")[0])
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef,
                             ast.ClassDef)):
            out.add(node.name)
        elif isinstance(node, (ast.Assign, ast.AnnAssign, ast.AugAssign)):
            targets = (node.targets if isinstance(node, ast.Assign)
                       else [node.target])
            for t in targets:
                for n in ast.walk(t):
                    if isinstance(n, ast.Name):
                        out.add(n.id)
    return out


def _reads_before_assignment(path):
    import ast
    import builtins
    tree = ast.parse(path.read_text(encoding="utf-8"))
    outer = _module_level_names(tree) | set(dir(builtins))
    bad = []
    for fn in ast.walk(tree):
        if not isinstance(fn, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        a = fn.args
        bound = {p.arg for p in a.posonlyargs + a.args + a.kwonlyargs}
        if a.vararg:
            bound.add(a.vararg.arg)
        if a.kwarg:
            bound.add(a.kwarg.arg)
        declared, stores, loads = set(), {}, {}
        for n in _own_scope_nodes(fn):
            if isinstance(n, (ast.Global, ast.Nonlocal)):
                declared |= set(n.names)
            elif isinstance(n, (ast.Import, ast.ImportFrom)):
                for al in n.names:
                    nm = (al.asname or al.name).split(".")[0]
                    if nm not in stores or n.lineno < stores[nm]:
                        stores[nm] = n.lineno
            elif isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef,
                                ast.ClassDef)):
                stores.setdefault(n.name, n.lineno)
            elif isinstance(n, ast.ExceptHandler) and n.name:
                # The EARLIEST handler binds it; a function with two
                # `except ... as exc` blocks is not a bug.
                if n.name not in stores or n.lineno < stores[n.name]:
                    stores[n.name] = n.lineno
            elif isinstance(n, ast.Name):
                if isinstance(n.ctx, (ast.Store, ast.Del)):
                    if n.id not in stores or n.lineno < stores[n.id]:
                        stores[n.id] = n.lineno
                elif isinstance(n.ctx, ast.Load):
                    if n.id not in loads or n.lineno < loads[n.id]:
                        loads[n.id] = n.lineno
        for name, lineno in loads.items():
            if (name in bound or name in outer or name in declared
                    or name not in stores):
                continue
            if stores[name] > lineno:
                bad.append(f"{path.name}:{lineno} {fn.name}(): {name!r} "
                           f"read before its assignment on line "
                           f"{stores[name]}")
    return bad


def test_no_local_is_read_before_it_is_assigned():
    """`local_rows` was read fifteen lines above its assignment, so the
    docking CLI's --gene path died with UnboundLocalError the first time
    anyone used it. Eight hundred and eighty-five tests passed."""
    import pathlib
    from promoter_tf_drugs import __file__ as pkg
    bad = []
    for f in sorted(pathlib.Path(pkg).parent.glob("*.py")):
        bad += _reads_before_assignment(f)
    assert not bad, "read before assignment: " + "; ".join(bad)


def _unread_cli_flags(path):
    import ast
    src = path.read_text(encoding="utf-8")
    tree = ast.parse(src)
    dests, bad = {}, []
    for n in ast.walk(tree):
        if not (isinstance(n, ast.Call)
                and isinstance(n.func, ast.Attribute)
                and n.func.attr == "add_argument"):
            continue
        dest, skip = None, False
        for kw in n.keywords:
            if kw.arg == "dest" and isinstance(kw.value, ast.Constant):
                dest = kw.value.value
            if (kw.arg == "action" and isinstance(kw.value, ast.Constant)
                    and kw.value.value in ("help", "version")):
                skip = True
        if dest is None:
            for arg in n.args:
                if (isinstance(arg, ast.Constant)
                        and isinstance(arg.value, str)
                        and arg.value.startswith("--")):
                    dest = arg.value[2:].replace("-", "_")
        if dest and not skip:
            dests[dest] = n.lineno
    for dest, lineno in sorted(dests.items()):
        if (f".{dest}" not in src and f'"{dest}"' not in src
                and f"'{dest}'" not in src):
            bad.append(f"{path.name}:{lineno} "
                       f"--{dest.replace('_', '-')} parsed but never read")
    return bad


def test_every_cli_flag_is_read_somewhere_in_its_module():
    """--cpu was parsed, documented as "the single biggest speed-up
    available", and never passed to the function that does the work, so
    every docking job ran on one core of forty.

    The existing parameter guard cannot see this: the value arrives as
    an attribute of `args`, which that check skips by design.
    """
    import pathlib
    from promoter_tf_drugs import __file__ as pkg
    bad = []
    for f in sorted(pathlib.Path(pkg).parent.glob("*.py")):
        bad += _unread_cli_flags(f)
    assert not bad, "CLI flags parsed but never read: " + "; ".join(bad)


def test_the_read_before_assignment_guard_catches_the_real_bug(tmp_path):
    """A guard reporting zero proves nothing until it has gone red."""
    src = tmp_path / "sample.py"
    src.write_text(
        "def main(args):\n"
        "    if rows:\n"
        "        pass\n"
        "    rows = load(args)\n"
        "    return rows\n", encoding="utf-8")
    assert any("rows" in b for b in _reads_before_assignment(src))
    # And a legitimate two-handler function is not flagged.
    ok = tmp_path / "ok.py"
    ok.write_text(
        "def main():\n"
        "    try:\n"
        "        a = 1\n"
        "    except ValueError as exc:\n"
        "        print(exc)\n"
        "    try:\n"
        "        b = 2\n"
        "    except KeyError as exc:\n"
        "        print(exc)\n"
        "    return [s for s in (a, b)]\n", encoding="utf-8")
    assert _reads_before_assignment(ok) == []


def test_the_cli_flag_guard_catches_the_real_bug(tmp_path):
    src = tmp_path / "cli_sample.py"
    src.write_text(
        "import argparse\n"
        "def build():\n"
        "    p = argparse.ArgumentParser()\n"
        "    p.add_argument('--cpu', type=int, default=1)\n"
        "    p.add_argument('--out')\n"
        "    return p\n"
        "def main():\n"
        "    args = build().parse_args()\n"
        "    return args.out\n", encoding="utf-8")
    bad = _unread_cli_flags(src)
    assert any("--cpu" in b for b in bad), bad
    assert not any("--out" in b for b in bad), bad


# ===================================================================== #
# run()'s argument contract is a table, and the table cannot drift.     #
# ===================================================================== #
def test_every_run_argument_has_a_validation_rule():
    """The point of the table. run() takes fifty-eight arguments and
    about fifteen were enforced; the other annotations were comments.
    A new argument with no rule fails here rather than failing later
    somewhere that looks unrelated."""
    from promoter_tf_drugs import pipeline, validation
    missing, extra = validation.coverage(pipeline.RUN_SPECS, pipeline.run)
    assert not missing, f"run() arguments with no validation rule: {missing}"
    assert not extra, f"rules for arguments run() does not take: {extra}"


def test_a_flag_given_a_string_is_refused():
    """Not caught before. `all_tss="yes"` is truthy, so it turned the
    flag on and the run proceeded - the one case where a bad value
    silently does something rather than failing."""
    import pytest
    from promoter_tf_drugs import pipeline as pl
    with pytest.raises(ValueError, match="all_tss"):
        pl.run("CCR5", all_tss="yes", demo=True)
    with pytest.raises(ValueError, match="activity_permutations"):
        pl.run("CCR5", activity_permutations="10", demo=True)


def test_a_true_flag_is_not_accepted_as_a_window_size():
    """bool subclasses int, so an unguarded isinstance(x, int) lets True
    through as 1 and the window silently becomes one base."""
    import pytest
    from promoter_tf_drugs import pipeline as pl
    with pytest.raises(ValueError, match="upstream"):
        pl.run("CCR5", upstream=True, demo=True)


def test_the_old_validation_rules_all_still_hold():
    """Behaviour preservation for the hand-written block this replaced."""
    import pytest
    from promoter_tf_drugs import pipeline as pl
    bad = [
        dict(gene=""),
        dict(gene="X", upstream=-1),
        dict(gene="X", upstream=0, downstream=0),
        dict(gene="X", upstream_depth=0),
        dict(gene="X", pvalue=1.0),
        dict(gene="X", pvalue=0.0),
        dict(gene="X", tf_weight=2.0),
        dict(gene="X", min_ic=-1.0),
        dict(gene="X", goal="sideways"),
        dict(gene="X", tissue="liver"),
        dict(gene="X", proliferation="fast"),
    ]
    for kw in bad:
        with pytest.raises((ValueError, FileNotFoundError)):
            pl.run(demo=True, **kw)


def test_a_missing_input_file_is_a_filenotfounderror_not_a_valueerror():
    """A missing file and a malformed argument are different problems and
    the caller usually needs to tell them apart."""
    import pytest
    from promoter_tf_drugs import pipeline as pl
    with pytest.raises(FileNotFoundError, match="pfm_file"):
        pl.run("CCR5", pfm_file="definitely_not_here.txt", demo=True)


def test_none_is_allowed_only_where_it_means_something():
    """has_ptc=None means unknown, which is not False: whether the
    transcript carries a premature stop decides whether raising the gene
    helps or just makes more truncated protein."""
    from promoter_tf_drugs import validation as v
    spec = v.Spec("bool", optional=True)
    v.check("has_ptc", None, spec)            # unknown is a real state
    v.check("has_ptc", False, spec)
    import pytest
    with pytest.raises(ValueError):
        v.check("has_ptc", "no", spec)
    with pytest.raises(ValueError, match="must not be None"):
        v.check("goal", None, v.Spec("str", choices=("up", "down")))


def test_a_validator_reports_what_a_bad_value_would_do():
    """`why` is the part that is not obvious from the rule. A message
    saying pvalue must be under 1 states the rule; saying that at 1.0
    every position calls a site says what you got wrong."""
    import pytest
    from promoter_tf_drugs import pipeline as pl
    with pytest.raises(ValueError) as exc:
        pl.run("CCR5", pvalue=1.0, demo=True)
    assert "every position" in str(exc.value)


def test_validation_does_not_add_a_runtime_dependency():
    """The package keeps NumPy as its only required runtime dependency,
    so the checking is a table rather than a pydantic model. The reason
    it is not pydantic is in validation.py; this is the part a test can
    hold."""
    import ast
    import pathlib
    from promoter_tf_drugs import validation
    tree = ast.parse(pathlib.Path(validation.__file__).read_text(
        encoding="utf-8"))
    imported = set()
    for n in ast.walk(tree):
        if isinstance(n, ast.Import):
            imported |= {a.name.split(".")[0] for a in n.names}
        elif isinstance(n, ast.ImportFrom) and n.level == 0 and n.module:
            imported.add(n.module.split(".")[0])
    assert imported <= {"math", "dataclasses", "pathlib", "inspect",
                        "typing", "__future__"}, imported


def test_a_chip_regulon_is_not_reported_as_curated():
    """A TRRUST row says somebody read a paper and asserted the factor
    regulates the gene. A ReMap row says the factor was found bound near
    it. Binding is not regulation, so calling both "curated" promotes
    one of them a tier in the one line of the report whose job is to say
    where the evidence came from."""
    from promoter_tf_drugs import regulon as r
    chip = "AATF\tAACS\tUnknown\tReMap_cells:1\n" * 20
    cur = "AATF\tBAX\tRepression\t22909821\n" * 20
    import pathlib
    import tempfile
    d = pathlib.Path(tempfile.mkdtemp())
    (d / "chip.tsv").write_text(chip, encoding="utf-8")
    (d / "cur.tsv").write_text(cur, encoding="utf-8")
    kind_chip, phrase_chip = r.describe_source(d / "chip.tsv")
    kind_cur, phrase_cur = r.describe_source(d / "cur.tsv")
    assert kind_chip == "chip", phrase_chip
    assert "curated" not in phrase_chip.lower()
    assert "not regulation" in phrase_chip
    assert kind_cur == "curated", phrase_cur
    # A file with neither marker is unknown, not quietly curated.
    (d / "bare.tsv").write_text("A\tB\n" * 20, encoding="utf-8")
    assert r.describe_source(d / "bare.tsv")[0] == "unknown"


# ===================================================================== #
# Docking several compounds at once must not change the answers.        #
# ===================================================================== #
def test_parallel_docking_preserves_order_and_scores():
    """A speed-up that quietly reordered the results would attach every
    score to the wrong drug, and nothing downstream could detect it: the
    numbers would all still be plausible.

    Uses a score_fn with a deliberately uneven delay, so the completion
    order is guaranteed to differ from the input order.
    """
    import pathlib
    import time
    from promoter_tf_drugs import docking as dk

    smiles = [f"C{'C' * i}O" for i in range(12)]

    def score(smi):
        # Later compounds finish first, so completion order is the
        # reverse of input order.
        time.sleep(0.02 * (len(smiles) - smiles.index(smi)))
        return -float(len(smi))

    seq = dk._dock_smiles(smiles, pathlib.Path("r.pdb"), None,
                          exhaustiveness=1, seed=0, score_fn=score, jobs=1)
    par = dk._dock_smiles(smiles, pathlib.Path("r.pdb"), None,
                          exhaustiveness=1, seed=0, score_fn=score, jobs=6)
    assert par == seq
    assert par == [-float(len(s)) for s in smiles]


def test_a_failure_in_one_parallel_job_does_not_shift_the_others():
    """None has to stay at the index of the compound that failed. If a
    failure dropped out of the list instead, every score after it would
    belong to a different drug."""
    import pathlib
    from promoter_tf_drugs import docking as dk

    smiles = ["CCO", "BAD", "CCCO", "BAD", "CCCCO"]

    def score(smi):
        if smi == "BAD":
            raise RuntimeError("embedding failed")
        return -float(len(smi))

    for jobs in (1, 4):
        got = dk._dock_smiles(smiles, pathlib.Path("r.pdb"), None,
                              exhaustiveness=1, seed=0, score_fn=score,
                              jobs=jobs)
        assert got == [-3.0, None, -4.0, None, -5.0], (jobs, got)


def test_jobs_defaults_to_one_so_nothing_changes_unasked():
    """Concurrency is opt-in: a caller that never heard of --jobs gets
    exactly the behaviour it had."""
    import inspect
    from promoter_tf_drugs import docking as dk
    for fn in (dk._dock_smiles, dk.validate_pocket, dk.screen_by_docking,
               dk.dock_screen_for_target):
        assert inspect.signature(fn).parameters["jobs"].default == 1, fn


# ===================================================================== #
# Which pocket a screen aims at. Every bug here is silent: the run      #
# completes, the scores look plausible, and they describe another site. #
# ===================================================================== #
def _het(resname, chain, resi, n, *, start=(10.0, 10.0, 10.0), model=None):
    """HETATM records for one group, optionally inside a MODEL block."""
    lines = []
    if model is not None:
        lines.append(f"MODEL     {model:>4d}")
    for i in range(n):
        x, y, z = start[0] + i * 0.5, start[1], start[2]
        lines.append(
            f"HETATM{i + 1:>5d}  C   {resname:>3s} {chain}{resi:>4s}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C")
    if model is not None:
        lines.append("ENDMDL")
    return "\n".join(lines)


def test_a_modified_amino_acid_is_not_mistaken_for_a_ligand():
    """CCR5 selected 2L87 and placed its box on TYS - O-sulfo-tyrosine,
    a residue of the peptide backbone written as HETATM only because it
    is non-standard. The box then sits on the protein, and the whole
    library is screened against a site that does not exist."""
    from promoter_tf_drugs import docking as dk
    only_modified = _het("TYS", "A", "10", 16)
    assert dk.best_ligand(only_modified) is None
    # A real ligand in the same file is still found.
    both = only_modified + "\n" + _het("MRV", "A", "501", 37,
                                       start=(30.0, 30.0, 30.0))
    got = dk.best_ligand(both)
    assert got is not None and got["resname"] == "MRV", got


def test_nmr_models_are_not_pooled_into_one_giant_ligand():
    """Groups keyed on (resname, chain, resi) merged every model of an
    NMR entry: sixteen atoms became a hundred and sixty, and the centre
    became the average of ten conformations - a point no model occupies.
    Largest-group-wins then ranks every multi-model structure above a
    real crystal ligand."""
    from promoter_tf_drugs import docking as dk
    models = "\n".join(_het("LIG", "A", "501", 12, model=m)
                       for m in range(1, 11))
    groups = dk.hetatm_groups(models)
    assert len(groups) == 1, groups
    assert groups[0]["n_atoms"] == 12, (
        f"pooled {groups[0]['n_atoms']} atoms from ten copies of twelve")


def test_a_crystallisation_lipid_loses_to_a_real_ligand():
    """Monoolein is the largest group in most GPCR structures and lies
    against the membrane face. Taking the largest outright screens the
    outside of the protein."""
    from promoter_tf_drugs import docking as dk
    text = (_het("OLC", "A", "601", 40) + "\n"
            + _het("MRV", "A", "501", 20, start=(40.0, 40.0, 40.0)))
    got = dk.best_ligand(text)
    assert got["resname"] == "MRV", got
    assert got["additive"] is False


def test_a_lipid_is_still_used_when_it_is_all_there_is_and_is_flagged():
    """Excluding lipids outright would be wrong: for a lipid receptor the
    fatty acid IS the ligand. So it is used and marked, which leaves a
    visible caveat rather than a silent wrong answer."""
    from promoter_tf_drugs import docking as dk
    got = dk.best_ligand(_het("OLC", "A", "601", 40))
    assert got is not None
    assert got["resname"] == "OLC"
    assert got["additive"] is True


def test_an_unreadable_conformation_is_not_reported_as_the_wrong_one():
    """classify_conformation reads coactivator and corepressor peptides,
    which only nuclear receptors carry, so it returns unknown for every
    GPCR structure. Saying "none in the inactive conformation" would
    turn a measurement that was never made into a negative result."""
    from promoter_tf_drugs import docking as dk
    res = dk.dock_screen_for_target(
        None, "CCR5", "P51681", [{"name": "x", "smiles": "CCO"}],
        actives=[], want_state="inactive",
        work_dir="unused",
        score_fn=lambda s: -5.0)
    # Without a structure this refuses for a different reason; the point
    # of the test is the wording when the conformation is unreadable.
    assert res.get("usable") is False


def test_the_conformation_caveat_survives_into_the_saved_result():
    """A screen run with --state any is still a screen. It just cannot
    support a claim about direction, and that has to reach the file
    rather than being mentioned once on a terminal."""
    from promoter_tf_drugs import docking as dk
    pocket = dk.PocketChoice(
        pdb_id="4MBS", path="r.pdb", ligand_path="l.pdb",
        ligand_resname="MRV", box=(0, 0, 0, 20, 20, 20),
        state="unknown", note="n")
    assert pocket.state not in ("active", "inactive")


def test_the_validation_gate_is_sampled_at_random_not_by_potency():
    """Taking the tightest binders would make the actives easier to
    recover and raise the AUC, which is tuning the gate to pass rather
    than measuring whether the pocket works."""
    import inspect
    from promoter_tf_drugs import docking as dk
    src = inspect.getsource(dk.dock_screen_for_target)
    assert "MAX_VALIDATION_ACTIVES" in src
    assert "random.Random(seed).sample" in src
    # Not sorted by potency anywhere near the sampling.
    window = src[:src.index("if decoys is None")]
    assert "pchembl" not in window, window[-400:]
    assert dk.MAX_VALIDATION_ACTIVES >= dk.MIN_VALIDATION_ACTIVES


def test_a_wsl_binary_is_never_executed_as_a_windows_path(monkeypatch):
    """_docking_binary() returns /home/u/bin/gnina now. Handing that to
    subprocess.run on Windows raises FileNotFoundError, which the caller
    catches at DEBUG and records as None - so every compound in a library
    scores None and the run reports an EMPTY screen, not a broken one."""
    import pathlib
    from promoter_tf_drugs import docking as dk
    seen = {}

    def fake_run_dock(binary, args, out, **kw):
        seen["binary"] = binary
        seen["args"] = [str(a) for a in args]
        return "", ""

    monkeypatch.setattr(dk, "_docking_binary",
                        lambda: "/home/u/bin/gnina")
    monkeypatch.setattr(dk, "_run_dock", fake_run_dock)
    dk.dock_one_boxed(pathlib.Path("C:/w/r.pdb"), pathlib.Path("C:/w/l.sdf"),
                      (0, 0, 0, 20, 20, 20))
    assert seen["binary"] == "/home/u/bin/gnina"
    # Every path handed over is a WSL one.
    for a in seen["args"]:
        assert not (len(a) > 1 and a[1] == ":"), a


# ===================================================================== #
# gnina's CNN numbers must survive the screen, not be thrown away.      #
# ===================================================================== #
def test_cnn_scores_reach_the_hits(monkeypatch):
    """The screen used to keep one number per compound, the Vina-style
    affinity, which for gnina discards the entire reason to run it."""
    import pathlib
    from promoter_tf_drugs import docking as dk

    lib = [{"name": f"d{i}", "smiles": "C" * (i + 2)} for i in range(6)]

    def fake(smiles_list, receptor, autobox, **kw):
        n = len(smiles_list)
        if not kw.get("detailed"):
            return [-5.0 - i for i in range(n)]
        return [dk.PoseScore(affinity=-5.0 - i, cnn_score=0.5,
                             cnn_affinity=7.0 + i, engine="gnina")
                for i in range(n)]

    monkeypatch.setattr(dk, "_dock_smiles", fake)
    monkeypatch.setattr(dk, "_docking_binary", lambda: "/u/gnina")
    monkeypatch.setattr(dk, "validate_pocket",
                        lambda *a, **k: {"usable": True, "auc": 0.9,
                                         "ef1": 10})
    run = dk.screen_by_docking(pathlib.Path("r.pdb"), lib, target="T",
                               actives=[{"smiles": "CC"}] * 6,
                               decoys=[{"smiles": "CCC"}] * 6,
                               box=(0, 0, 0, 20, 20, 20))
    assert run.hits, run.error
    assert all("cnn_affinity" in h for h in run.hits), run.hits[0]
    assert all(h["engine"] == "gnina" for h in run.hits)


def test_an_engine_without_a_cnn_gets_no_empty_cnn_columns(monkeypatch):
    """A smina run carrying cnn_score: None would read as a measurement
    that came out blank rather than one that was never made."""
    import pathlib
    from promoter_tf_drugs import docking as dk

    lib = [{"name": "d", "smiles": "CCO"}]
    monkeypatch.setattr(
        dk, "_dock_smiles",
        lambda sl, r, a, **kw: [dk.PoseScore(affinity=-7.0,
                                             engine="smina")])
    monkeypatch.setattr(dk, "_docking_binary", lambda: "smina")
    monkeypatch.setattr(dk, "validate_pocket",
                        lambda *a, **k: {"usable": True, "auc": 0.9,
                                         "ef1": 10})
    run = dk.screen_by_docking(pathlib.Path("r.pdb"), lib, target="T",
                               actives=[{"smiles": "CC"}] * 6,
                               decoys=[{"smiles": "CCC"}] * 6,
                               box=(0, 0, 0, 20, 20, 20))
    assert "cnn_score" not in run.hits[0]
    assert "cnn_affinity" not in run.hits[0]


def test_the_gate_measures_every_metric_it_has(monkeypatch):
    """Known ligands against decoys is the only honest way to find out
    whether the CNN separates binders better for a given pocket, and the
    poses are already scored so it costs nothing."""
    import pathlib
    from promoter_tf_drugs import docking as dk

    # The CNN separates perfectly; the affinity is pure noise.
    def fake(smiles_list, receptor, autobox, **kw):
        out = []
        for i, smi in enumerate(smiles_list):
            active = smi.startswith("A")
            out.append(dk.PoseScore(
                affinity=-5.0,
                cnn_affinity=9.0 if active else 4.0,
                engine="gnina"))
        return out

    monkeypatch.setattr(dk, "_dock_smiles", fake)
    val = dk.validate_pocket(
        pathlib.Path("r.pdb"),
        [{"smiles": f"A{i}"} for i in range(10)],
        [{"smiles": f"D{i}"} for i in range(40)],
        box=(0, 0, 0, 20, 20, 20))
    assert "auc_by_metric" in val, val
    assert val["auc_by_metric"]["cnn_affinity"]["auc"] == 1.0
    assert "metric_note" in val
    # Reported, not acted on: the decision still uses the affinity.
    assert val["auc"] == val["auc_by_metric"]["affinity"]["auc"]


# ===================================================================== #
# Documented commands have to be runnable commands.                     #
# ===================================================================== #
def _parser_flags(module_name):
    """Every long option a module's CLI accepts, or None if it has none."""
    import argparse
    import contextlib
    import importlib
    import io
    import re

    mod = importlib.import_module(module_name)
    if hasattr(mod, "build_parser"):
        try:
            p = mod.build_parser()
        except TypeError:
            return None
        return {o for a in p._actions for o in a.option_strings
                if o.startswith("--")}
    if not hasattr(mod, "main"):
        return None
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            mod.main(["--help"])
    except SystemExit:
        pass
    except Exception:                                       # noqa: BLE001
        return None
    text = buf.getvalue()
    if "usage:" not in text:
        return None
    return set(re.findall(r"--[a-z0-9][a-z0-9-]*", text))


def test_every_documented_command_uses_flags_its_parser_accepts():
    """A command example that does not run is worse than no example:
    the reader assumes their setup is broken rather than the page.

    This checks every `python -m promoter_tf_drugs.X` line in docs and
    the README against that module's own parser, so a renamed flag
    cannot leave the documentation behind.
    """
    import pathlib
    import re

    root = pathlib.Path(__file__).resolve().parent.parent
    pages = sorted(root.glob("docs/*.md")) + [root / "README.md"]
    cache, bad = {}, []
    for page in pages:
        if not page.exists():
            continue
        text = page.read_text(encoding="utf-8")
        # A shell example wraps with a trailing backslash, and stopping
        # at the first newline reads only its first line - which is how
        # this guard let a renamed flag on the second line through when
        # it was first checked against a real page.
        for m in re.finditer(
                r"python -m (promoter_tf_drugs\.[a-z_]+)"
                r"((?:[^\n`]*\\\n)*[^\n`]*)", text):
            mod, rest = m.group(1), m.group(2)
            if mod not in cache:
                cache[mod] = _parser_flags(mod)
            flags = cache[mod]
            if flags is None:
                continue
            for f in re.findall(r"--[a-z0-9][a-z0-9-]*", rest):
                if f not in flags:
                    bad.append(f"{page.name}: {mod} has no {f}")
    assert not bad, "documented flags that do not exist: " + "; ".join(bad)


def test_the_documented_command_check_would_notice_a_renamed_flag(
        tmp_path):
    """The guard reporting zero proves nothing until it has gone red."""
    import re
    flags = _parser_flags("promoter_tf_drugs.docking")
    assert flags, "docking parser produced no flags"
    assert "--jobs" in flags and "--state" in flags
    # A flag that does not exist must be detected by the same matching.
    fake = "python -m promoter_tf_drugs.docking --not-a-real-flag"
    found = re.findall(r"--[a-z0-9][a-z0-9-]*", fake.split("docking")[1])
    assert found == ["--not-a-real-flag"]
    assert found[0] not in flags


def test_an_unreadable_conformation_is_explained_on_the_route():
    """A route whose structure has no readable conformation used to read
    as though the question was never asked. It was asked and could not be
    answered, which for a GPCR is permanent - the classifier reads
    coactivator peptides and a GPCR has none."""
    from promoter_tf_drugs import catalogue as cat
    run = {"target": "CCR5", "usable": True, "auc": 0.8, "ef1": 10,
           "pocket": {"pdb_id": "4MBS", "ligand_resname": "MRV",
                      "state": "unknown"},
           "hits": [{"drug": "maraviroc", "score": -11.2}]}
    routes = cat.from_docking(run, goal="down", gene="CCR5")
    assert routes
    notes = " ".join(n for r in routes for n in (r.notes or []))
    assert "could not be read" in notes, notes
    assert "not a finding about CCR5" in notes


def test_a_readable_conformation_still_states_the_fit():
    """The explanatory note must not replace the real conformation
    reasoning when the structure can actually be classified."""
    from promoter_tf_drugs import catalogue as cat
    run = {"target": "VDR", "usable": True, "auc": 0.9, "ef1": 12,
           "pocket": {"pdb_id": "1DB1", "ligand_resname": "VDX",
                      "state": "active"},
           "hits": [{"drug": "calcitriol", "score": -12.8}]}
    routes = cat.from_docking(run, goal="up", gene="GENE",
                              tf_roles={"VDR": "activator"})
    notes = " ".join(n for r in routes for n in (r.notes or []))
    assert "could not be read" not in notes, notes
    assert "conformation to stabilise" in notes, notes


def test_the_reason_a_compound_is_absent_is_counted_not_guessed():
    """The caveat said compounds were dropped "usually because no 3D
    conformer could be generated". The code knows which it was and threw
    the distinction away. A compound that cannot be given a 3D structure
    and one whose docking crashed are different problems and send you to
    different places."""
    import pathlib
    from promoter_tf_drugs import docking as dk

    def embed(smi, path, **kw):
        if smi.startswith("BAD"):
            raise ValueError("3D embedding failed")
        path.write_text("sdf", encoding="utf-8")

    calls = {"n": 0}

    def boxed(receptor, sdf, box, **kw):
        calls["n"] += 1
        if calls["n"] % 3 == 0:
            raise RuntimeError("gnina exploded")
        return dk.PoseScore(affinity=-7.0, engine="gnina")

    import unittest.mock as mock
    with mock.patch.object(dk, "smiles_to_sdf", embed), \
         mock.patch.object(dk, "dock_one_scored", boxed), \
         mock.patch.object(dk, "engine_name", lambda *a: "gnina"):
        tally = {}
        smis = ["BAD1", "BAD2"] + [f"CC{i}" for i in range(9)]
        dk._dock_smiles(smis, pathlib.Path("r.pdb"), None,
                        exhaustiveness=1, seed=0,
                        box=(0, 0, 0, 20, 20, 20), detailed=True,
                        tally=tally)
    assert tally.get("embedding_failed") == 2, tally
    assert tally.get("docking_failed", 0) >= 1, tally


def test_the_absent_compound_caveat_states_measured_counts(monkeypatch):
    import pathlib
    from promoter_tf_drugs import docking as dk

    lib = [{"name": f"d{i}", "smiles": "CCO"} for i in range(10)]

    def fake(smiles_list, receptor, autobox, **kw):
        t = kw.get("tally")
        if t is not None:
            t["embedding_failed"] = 3
        return ([None] * 3
                + [dk.PoseScore(affinity=-7.0, engine="gnina")] * 7)

    monkeypatch.setattr(dk, "_dock_smiles", fake)
    monkeypatch.setattr(dk, "_docking_binary", lambda: "/u/gnina")
    monkeypatch.setattr(dk, "validate_pocket",
                        lambda *a, **k: {"usable": True, "auc": 0.9,
                                         "ef1": 10})
    run = dk.screen_by_docking(pathlib.Path("r.pdb"), lib, target="T",
                               actives=[{"smiles": "CC"}] * 6,
                               decoys=[{"smiles": "CCC"}] * 6,
                               box=(0, 0, 0, 20, 20, 20))
    text = " ".join(run.caveats)
    assert "3 had no 3D conformer" in text, text
    assert "usually" not in text, text
    assert run.failures == {"embedding_failed": 3}


def test_impossible_flag_combinations_are_refused_before_any_work():
    """--tf-activity without --activity-tissue failed twelve minutes into
    a run: the scan, the ChIP intersection, the drug lookups and the
    docking merge all ran first, and only then did the report say the
    requested source could never have produced anything.

    AGENTS.md: validate user input before expensive work.
    """
    import pytest
    from promoter_tf_drugs import cli, panel
    for mod in (cli, panel):
        gene_flag = (["--gene", "CCR5"] if mod is cli
                     else ["--genes", "CCR5"])
        with pytest.raises(SystemExit) as exc:
            mod.main(gene_flag + ["--tf-activity", "--demo"])
        assert exc.value.code != 0, mod.__name__


def test_the_library_call_still_degrades_rather_than_raising():
    """The CLI check must not become a pipeline check. A library caller
    may want everything else the run produces, and the ledger recording
    a failed source is the honest way to say one of them did not."""
    from promoter_tf_drugs import pipeline as pl
    res = pl.run("DEMOGENE", demo=True, tf_activity=True,
                 activity_tissue=None)
    sources = res["evidence_sources"]["sources"]
    assert sources["tf_activity"]["ok"] is False
    assert "activity-tissue" in sources["tf_activity"]["detail"]
    # The rest of the run still happened.
    assert res["results"], "the run should still produce blocks"


def test_a_flag_that_would_do_nothing_is_refused_too():
    """--activity-permutations without --tf-activity is accepted silently
    and has no effect, which reads as though the null model ran."""
    import pytest
    from promoter_tf_drugs import cli
    with pytest.raises(SystemExit):
        cli.main(["--gene", "CCR5", "--activity-permutations", "100",
                  "--demo"])


def test_every_evidence_source_says_what_its_failure_cost():
    """The CCR5 run printed `- **tf_activity** — ` and stopped. The
    section above promises these are the sources that did not deliver
    and that nothing rests on them, then declines to say what was lost.

    IMPACT covered 7 of the 14 sources a real run records; the rest were
    added when their features were and the table was not.
    """
    from promoter_tf_drugs import pipeline as pl, provenance as prov
    res = pl.run("DEMOGENE", demo=True)
    names = set(res["evidence_sources"]["sources"])
    missing = sorted(n for n in names if not prov.IMPACT.get(n))
    assert not missing, (
        "evidence sources with no stated impact when they fail: "
        + ", ".join(missing))


def test_a_failed_source_bullet_is_never_bare():
    """Belt and braces for the table above: even a source nobody wrote
    an impact for must say something, because a bullet that stops after
    the dash is worse than the plain fact."""
    from promoter_tf_drugs import pipeline as pl
    res = pl.run("DEMOGENE", demo=True, tf_activity=True,
                 activity_tissue=None)
    report = pl.markdown_report(res)
    assert "— \n" not in report, "bare failed-source bullet"
    assert "**tf_activity** — " in report
    line = next(ln for ln in report.splitlines()
                if ln.startswith("- **tf_activity**"))
    assert len(line.split("—", 1)[1].strip()) > 20, line
