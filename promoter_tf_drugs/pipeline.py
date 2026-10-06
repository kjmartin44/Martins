"""Orchestration: gene -> promoter -> motifs -> activators -> approved drugs."""

from __future__ import annotations

import json
import logging
import math
from dataclasses import dataclass, field
from pathlib import Path

from . import (annotate, catalogue as cat_mod, cellcycle, congeners as cong_mod,
               delivery as dlv_mod,
               liabilities as liab_mod,
               sensitivity as sens_mod,
               promoterdna as dna_mod,
               tfactivity as tfa_mod,
               depmap as depmap_mod,
               corpora as corpora_mod,
               abundance as abund_mod,
               demo_data, direct, discovery, drug_scoring, eqtl as eqtl_mod,
               modalities as modal_mod,
               drugs as drugs_mod, evidence, genes, motifs as motifs_mod,
               genetics as gen_mod, ligands as lig_mod, modulators as mod_mod,
               neighbours as nbr_mod, network as net_mod,
               predictions as pred_mod, systemic as sys_mod,
               analogs as analog_mod,
               compensation as comp_mod,
               posttranscriptional as ptx_mod,
               reactivation as react_mod,
               outcome as outcome_mod, provenance as prov_mod,
               regulon as regulon_mod,
               screens as screens_mod, signatures as sig_mod,
               stabilizers as stab_mod)
from .http import Client

log = logging.getLogger(__name__)

# Transparent, tunable. These are heuristics, not a validated model.
WEIGHTS = {
    "motif_quality": 30,    # best relative PWM score for this TF
    "site_count": 10,       # saturating bonus for multiple independent sites
    "proximity": 15,        # closeness of best site to the TSS
    "activator_role": 15,   # GO-annotated activating function
    "chip_evidence": 20,    # measured binding over this promoter
    "druggable": 10,        # an approved drug pushes it the right way
}


@dataclass
class TFResult:
    symbol: str
    annotation: dict
    hits: list = field(default_factory=list)
    chip_peaks: int = 0
    drugs: list = field(default_factory=list)
    score: float = 0.0
    # How many other genes this factor controls. Set after ranking, because
    # it needs the regulon table. None means nobody has counted, which is a
    # different statement from a small number.
    n_regulon: int | None = None
    regulon_provenance: str | None = None
    regulon_understudied: bool = False
    # What a cell loses without this factor, from DepMap CRISPR knockout.
    # The same collateral question as n_regulon, measured instead of
    # curated, and available for factors no curator has written up.
    dependency: str | None = None
    dependency_burden: float | None = None

    @property
    def best_hit(self):
        return max(self.hits, key=lambda h: h.score) if self.hits else None

    def to_dict(self) -> dict:
        b = self.best_hit
        return {
            "symbol": self.symbol,
            "name": self.annotation.get("name"),
            "role": self.annotation.get("role"),
            "nuclear_receptor": self.annotation.get("nuclear_receptor"),
            "n_sites": len(self.hits),
            "best_rel_score": round(b.rel_score, 3) if b else None,
            "best_score_bits": round(b.score, 2) if b else None,
            "best_site_tss_offset": b.tss_offset if b else None,
            "best_site_sequence": b.matched_sequence if b else None,
            "matrix_id": b.matrix_id if b else None,
            "chip_peaks": self.chip_peaks,
            "score": round(self.score, 1),
            "n_regulon": self.n_regulon,
            "regulon_provenance": self.regulon_provenance,
            "regulon_understudied": self.regulon_understudied,
            "dependency": self.dependency,
            "dependency_burden": self.dependency_burden,
            "drugs": [d.to_dict() for d in self.drugs],
        }


# How many ranked factors the LINCS drug scan covers. Each is a pass over
# roughly 33,000 signatures, about a second, so this is a wall-clock cap
# rather than a claim that the rest are uninteresting.
# A blank table cell. Named rather than inlined because a
# backslash escape inside an f-string expression is PEP 701
# syntax: legal from 3.12 and a SyntaxError on 3.10 and 3.11,
# which pyproject declares support for. The module would not
# import at all on those versions.
EMDASH = "\u2014"

# How many ranked factors get a sensitivity scan. Each one is a
# correlation over the whole PRISM matrix, so this is the expensive
# step in that section rather than the file read, which happens once.
SENSITIVITY_FACTORS = 6

ACTIVITY_DRUG_FACTORS = 8


def _optional(name: str, fn, default):
    """Run an optional analysis; never let its failure abort the run.

    Everything gated behind a flag is an enhancement, not a requirement. A
    dead third-party API should cost you that section, not thirty minutes of
    completed work.
    """
    try:
        return fn()
    except Exception as exc:                              # noqa: BLE001
        log.warning("%s unavailable (%s: %s) — skipping that section; the rest "
                    "of the run is unaffected.", name, type(exc).__name__,
                    str(exc)[:160])
        return default


def _proximity(offset: int, scale: float = 400.0) -> float:
    """Sites near the TSS are likelier to be functional core/proximal elements."""
    return math.exp(-abs(offset) / scale)


# How much credit an annotated role earns, per goal. For downregulation a
# bound activator is the *most* actionable thing in the list, because blocking
# it is the easy half of pharmacology.
ROLE_CREDIT = {
    "up":   {"activator": 1.0, "dual": 0.6, "repressor": 0.3},
    "down": {"activator": 1.0, "dual": 0.6, "repressor": 0.9},
}


def score_tf(res: TFResult, has_chip: bool, goal: str = "up") -> float:
    b = res.best_hit
    if b is None:
        return 0.0
    w = WEIGHTS
    s = w["motif_quality"] * b.rel_score
    s += w["site_count"] * min(len(res.hits), 3) / 3
    s += w["proximity"] * _proximity(b.tss_offset)

    role = res.annotation.get("role", "unknown")
    s += w["activator_role"] * ROLE_CREDIT[goal].get(role, 0.0)

    if has_chip:
        s += w["chip_evidence"] * (1.0 if res.chip_peaks else 0.0)
    else:
        s += w["chip_evidence"] * 0.35 * min(b.rel_score, 1.0)  # unverified discount

    if res.drugs:
        s += w["druggable"]
    return s


def run(gene: str, *, upstream: int = 1000, downstream: int = 200,
        pvalue: float = 1e-4, species: str = "homo_sapiens",
        pfm_file: str | None = None, chip_bed: str | None = None,
        matrix_ids: list[str] | None = None, all_tss: bool = False,
        mask_repeats: bool = True, goal: str = "up",
        grade_drugs: bool = True, tissue: str | None = None,
        tf_weight: float = 0.5, predict_direction: bool = False,
        regulon_file: str | None = None, include_direct: bool = True,
        estimate_outcome: bool = True, screen_receptors: bool = False,
        find_improvements: bool = True, chip_api: bool = False,
        compare_mechanisms: bool = True, proliferation: str = "unknown",
        allow_unknown_direction: bool = False,
        cell_type_check: bool = False, include_stabilizers: bool = True,
        design_screen: bool = False, design_modality: bool = False,
        use_lincs: bool = False, signature_file: str | None = None,
        lincs_file: str | None = None,
        lincs_geneinfo=None, tf_activity: bool = False,
        depmap_file: str | None = None,
        prism_file: str | None = None,
        sensitivity_scan: bool = False,
        screen_product: bool = False,
        promoter_dna: bool = False,
        docking_results: str | None = None,
        activity_tissue: str | None = None,
        activity_permutations: int = 0, lincs_gmt: tuple | None = None, expression_file: str | None = None,
        promoter_fasta: str | None = None,
        has_ptc: bool | None = None,
        gtf: str | None = None,
        genome: str | None = None,
        use_genetics: bool = False,
        upstream_network: bool = False,
        upstream_depth: int = 3,
        approved_only: bool = True, use_dgidb: bool = True,
        use_opentargets: bool = True,
        find_analogs: bool = False,
        demo: bool = False, offline: bool = False,
        cache_dir: str | None = None, min_ic: float = 6.0) -> dict:
    """Run the full pipeline and return a JSON-serialisable result dict."""
    if not isinstance(gene, str) or not gene.strip():
        raise ValueError("gene must be a non-empty symbol or Ensembl ID")
    gene = gene.strip()
    for name, value in (("upstream", upstream), ("downstream", downstream),
                        ("upstream_depth", upstream_depth)):
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError(f"{name} must be a non-negative integer")
    if upstream + downstream == 0:
        raise ValueError("upstream and downstream must define a non-empty window")
    if upstream_depth < 1:
        raise ValueError("upstream_depth must be at least 1")
    if not math.isfinite(pvalue) or not 0 < pvalue < 1:
        raise ValueError("pvalue must be finite and between 0 and 1 (exclusive)")
    if not math.isfinite(tf_weight) or not 0 <= tf_weight <= 1:
        raise ValueError("tf_weight must be finite and between 0 and 1")
    if not math.isfinite(min_ic) or min_ic < 0:
        raise ValueError("min_ic must be finite and non-negative")
    if goal not in ("up", "down"):
        raise ValueError("goal must be 'up' or 'down'")
    if tissue not in (None, "cns"):
        raise ValueError("tissue must be None or 'cns'")
    if proliferation not in ("cycling", "post-mitotic", "mixed", "unknown"):
        raise ValueError("invalid proliferation state")
    for name, value in (("pfm_file", pfm_file), ("chip_bed", chip_bed),
                        ("regulon_file", regulon_file), ("lincs_file", lincs_file),
                        ("signature_file", signature_file),
                        ("expression_file", expression_file)):
        if value is not None and not Path(value).is_file():
            raise FileNotFoundError(f"{name}: file not found: {value}")
    client = (Client(cache_dir=cache_dir, offline=offline or demo)
              if cache_dir else Client(offline=offline or demo))

    _screen_cache: dict = {}   # receptor -> screen, shared across promoters
    _gene_linked: set = set()  # drugs CTD already reports as moving this gene
    _library: list = []        # approved-drug library, fetched once per run

    # Which evidence sources actually produced evidence. Populated as the run
    # goes; consulted by the scoring below and printed at the top of the
    # report, so a failed lookup can never be read as a measured negative.
    ledger = prov_mod.Ledger()

    regulons = regulon_mod.load_regulons(regulon_file) if regulon_file else None
    ledger.record("regulon", requested=bool(regulon_file),
                  ok=bool(regulons), yielded=len(regulons or {}),
                  detail=(f"{len(regulons or {})} TFs with curated targets"
                          if regulons else
                          "no regulon file; coarse bundled tiers in use"))
    if expression_file:
        dlv_mod.load_expression_file(expression_file)

    # Signed regulons are a different parse of the same file regulon_mod
    # reads: it returns bare sets and discards the Activation/Repression
    # column, which is exactly what a direction needs.
    # Docking happens on another machine when the binary lives
    # there; this is where its output comes back in.
    docking_runs: dict = {}
    if docking_results:
        from . import docking as dock_mod
        loaded = _optional(
            "docking results",
            lambda: dock_mod.load_results(docking_results), None)
        if loaded is not None:
            runs = loaded if isinstance(loaded, list) else [loaded]
            docking_runs = {r.get("target") or "?": r for r in runs}
            ok = [r for r in runs if r.get("usable")]
            ledger.record(
                "docking", requested=True, ok=True,
                yielded=sum(len(r.get("hits") or []) for r in ok),
                detail=(f"{len(ok)} of {len(runs)} pocket(s) validated"
                        if ok else
                        "no pocket validated, so no docking hit is "
                        "reported"))

    signed_regulons = None
    regulon_source = None
    if tf_activity:
        if regulon_file:
            signed_regulons = _optional(
                "signed regulon parse",
                lambda: tfa_mod.load_signed_regulons(regulon_file), None)
            regulon_source = regulon_file
        else:
            # Fall back to whatever is in data/, preferring CollecTRI. Asking
            # for --regulon-file when a usable table is sitting in the
            # default location is friction with no purpose.
            signed_regulons, regulon_source = _optional(
                "signed regulon autoload",
                tfa_mod.autoload_regulons, ({}, None))
            if not signed_regulons:
                log.warning(
                    "--tf-activity found no regulon table. Pass "
                    "--regulon-file, or put one at %s. %s",
                    " or ".join(tfa_mod.DEFAULT_REGULON_PATHS),
                    tfa_mod.COLLECTRI_SOURCE)
    if lincs_geneinfo:
        # Explicit path wins over the data/ autoload. Loaded before
        # any signature lookup, because the feature space decides the
        # evidence tier and a late load would tier the early hits as
        # predicted.
        sig_mod.load_geneinfo(lincs_geneinfo)
    if lincs_gmt:
        sig_mod.DEFAULT_GMT_PATHS['up'] = (lincs_gmt[0],)
        sig_mod.DEFAULT_GMT_PATHS['down'] = (lincs_gmt[1],)
    if gtf:
        # One GTF serves two purposes: the neighbour check, and the gene
        # models that used to require an Ensembl lookup.
        nbr_mod.load_gtf(gtf)
        genes.load_transcripts_gtf(gtf)
    if genome:
        genes.load_genome(genome)
    if not (gtf or genome) and not demo:
        # Pick up conventionally-placed files without being asked, so a
        # local corpus is used whenever it is present.
        found = genes.autoload_offline()
        if found.get("gtf") or found.get("genome"):
            log.info("offline sources found: %s", found)
    if promoter_fasta:
        # Removes the one stage that had no offline route.
        genes.load_promoter_fasta(promoter_fasta)

    # ------------------------------------- 0. drugs that skip the TF route
    # Done first because if one of these exists it outranks everything below.
    direct_hits = []
    # Defined before the branch: it was only ever assigned inside it, so
    # anything later that wanted the gene's own accession raised
    # NameError on a --no-direct run.
    uniprot = None
    if include_direct:
        if not demo:
            try:
                ann = annotate.annotate_tfs(client, [gene])
                uniprot = (ann.get(gene.upper()) or {}).get("uniprot")
            except Exception:                              # noqa: BLE001
                pass
        direct_hits = direct.find_direct(
            client, gene, uniprot, goal=goal,
            approved_only=approved_only, use_dgidb=use_dgidb, demo=demo)

    # ---------------------------------------------------------- 1. promoter
    if demo:
        seq = demo_data.demo_promoter_sequence(upstream, downstream)
        promoters = [genes.Promoter(
            gene_symbol=gene or "DEMOGENE", gene_id="ENSG00000000000",
            transcript_id="ENST00000000000", transcript_selection="demo fixture",
            chrom="1", strand=1, tss=1_000_000,
            start=1_000_000 - upstream, end=1_000_000 + downstream - 1,
            upstream=upstream, downstream=downstream,
            assembly="DEMO", sequence=seq)]
    else:
        promoters = genes.get_promoters(client, gene, species=species,
                                        upstream=upstream, downstream=downstream,
                                        all_tss=all_tss,
                                        mask="soft" if mask_repeats else None)

    # ------------------------------------------------------------ 2. motifs
    motif_set = (demo_data.demo_motifs() if demo
                 else motifs_mod.load_motifs(client, pfm_file, matrix_ids))
    log.info("scanning with %d motifs at p <= %g", len(motif_set), pvalue)

    outputs = []
    for prom in promoters:
        hits = motifs_mod.scan_promoter(
            prom.sequence, motif_set, pvalue=pvalue, min_ic=min_ic,
            offset_fn=prom.offset_to_tss, skip_repeats=mask_repeats and not demo)
        hits = motifs_mod.collapse_overlapping(hits)
        scannable = motifs_mod.scannable_length(
            prom.sequence, skip_repeats=mask_repeats and not demo)
        quality = motifs_mod.scan_quality(len(hits), len(motif_set),
                                          scannable, pvalue)
        log.info("%s: %d non-redundant sites", prom.transcript_id, len(hits))
        if quality["verdict"] == "no signal":
            log.warning("SCAN FOUND NO SIGNAL: %s", quality["advice"])
        elif quality["verdict"] == "weak":
            log.warning("%s", quality["advice"])

        # -------------------------------------------------- 3. annotate TFs
        symbols = sorted({s for h in hits for s in h.tf_symbols})
        anns = (demo_data.demo_annotations(symbols) if demo
                else annotate.annotate_tfs(client, symbols))

        by_tf: dict[str, TFResult] = {
            s: TFResult(symbol=s, annotation=dict(anns.get(s, {"symbol": s})))
            for s in symbols
        }
        for h in hits:
            for s in h.tf_symbols:
                by_tf[s].hits.append(h)

        # --------------------------------------------------- 4. ChIP evidence
        # ------------------- 3b. is this window uniquely this gene's?
        # Asked before any ranking, because a window containing another
        # transcription unit makes every site call below ambiguous.
        neighbours = ({} if demo else _optional(
            "neighbouring-gene check",
            lambda: nbr_mod.genes_in_window(
                client, prom.chrom, prom.start, prom.end,
                strand=prom.strand, species=species,
                exclude=prom.gene_symbol),
            {}))

        chip_counts, chip_status = evidence.chip_evidence_status(
            prom, chip_bed=chip_bed, client=client, use_api=chip_api)
        chip_src = ledger.record("chip", **chip_status)
        for s, r in by_tf.items():
            r.chip_peaks = chip_counts.get(s.upper(), 0)

        # ---------------------------------------------- 5. actionable subset
        actionable = {s: r.annotation for s, r in by_tf.items()
                      if r.annotation.get("role") in ("activator", "dual", "repressor")
                      or r.annotation.get("nuclear_receptor")}
        if demo:
            drug_map, unknown_direction = demo_data.demo_drugs(
                actionable, goal=goal), {}
        else:
            drug_map, unknown_direction = drugs_mod.drugs_for_tfs(
                client, actionable, approved_only=approved_only,
                use_dgidb=use_dgidb, use_opentargets=use_opentargets,
                goal=goal,
                allow_unknown_direction=allow_unknown_direction)
        for s, ds in drug_map.items():
            if s in by_tf:
                by_tf[s].drugs = ds

        # Only *usable* ChIP counts as measured binding. Asking for ChIP and
        # getting a broken response previously scored every factor a hard zero
        # on the 20-point evidence axis — presenting a failed lookup as though
        # binding had been measured and ruled out. A failed source falls back
        # to the unverified discount instead.
        # Why the drug tables are always nuclear receptors. 62 factors had
        # sites in a real CCR5 run; 5 had an approved drug; all 5 were NRs -
        # and the highest-SCORING factors (PDX1, RBPJ, NFATC1) had none and
        # so appeared in no drug table at all. The ranking below is therefore
        # a filter with one surviving class, not a ranking of your gene's
        # biology, and the report has to say so rather than letting the
        # output look like a finding.
        _drugged = [r for r in by_tf.values() if r.drugs]
        _nr_drugged = [r for r in _drugged
                       if r.annotation.get("nuclear_receptor")]
        attrition = {
            "n_tfs_with_sites": len(by_tf),
            "n_with_drug": len(_drugged),
            "n_drugged_nuclear_receptor": len(_nr_drugged),
            "all_drugged_are_nr": bool(_drugged) and
                                  len(_nr_drugged) == len(_drugged),
        }

        for r in by_tf.values():
            r.score = score_tf(r, has_chip=chip_src.usable, goal=goal)

        ranked = sorted(by_tf.values(), key=lambda r: -r.score)

        # The gene-specific signal lives in the high-scoring factors, and it
        # is exactly what gets discarded for having no ligand. Name them, and
        # route them into the searches that exist for this case.
        _top_nodrug = [r for r in ranked[:25] if not r.drugs]
        attrition["top25_without_drug"] = [r.symbol for r in _top_nodrug]
        attrition["n_top25_without_drug"] = len(_top_nodrug)

        # ------------------------- 5b. is this factor active in these cells?
        cycle = cellcycle.filter_and_annotate(ranked, proliferation)
        cycle_summary = cellcycle.summarise(cycle, proliferation)
        # Routes that cannot work in the target tissue are removed, not demoted.
        blocked = {s for s, a in cycle.items() if not a["compatible"]}
        if blocked:
            log.info("removed %d cell-cycle-incompatible TFs", len(blocked))
            ranked = [r for r in ranked if r.symbol not in blocked]

        # ------------------------------------------- 6. grade each TF/drug pair
        grades = []
        collapsed: dict = {}
        if grade_drugs:
            max_score = sum(WEIGHTS.values())
            grades = drug_scoring.grade_all(
                client, ranked, tissue=tissue, tf_weight=tf_weight,
                max_tf_score=max_score,
                evidence_lookup=demo_data.demo_evidence if demo else None,
                predict_direction=predict_direction and not demo, goal=goal,
                regulons=regulons)

            # A ChEMBL outage and a target with genuinely few ligands both
            # end as "no direction model", and the grade rows look identical.
            # Record which it was: the whole point of the ledger is that a
            # source which was asked and did not answer cannot pass as a
            # measured negative.
            if predict_direction and not demo:
                broke = sum(1 for g in grades if any(
                    "direction model unavailable" in n
                    or "direction prediction failed" in n
                    for n in (g.notes or [])))
                ledger.record(
                    "direction_model", requested=True,
                    ok=broke == 0 and bool(grades),
                    yielded=len(grades) - broke,
                    detail=(f"{broke} of {len(grades)} grade(s) could not be "
                            f"direction-checked" if broke else
                            f"{len(grades)} grade(s) direction-checked"))

            # Collapse salt forms. One pharmacology was occupying five rows
            # (Liothyronine / Liothyronine Sodium / Levothyroxine Sodium ...),
            # pushing distinct chemistry off the bottom of every table and
            # inflating the candidate count the Pareto accounting reports.
            # Done here so every downstream consumer — catalogue, outcome
            # chain, discovery — inherits one row per active moiety.
            grades, collapsed = cong_mod.collapse_grades(grades)
            for g in grades:
                alts = collapsed.get(g.drug_name)
                if alts:
                    g.notes.append(
                        "same active moiety as " + ", ".join(alts[:4])
                        + (f" (+{len(alts) - 4} more)" if len(alts) > 4 else ""))

        # ------------------------------ 7. propagate to a transcriptional call
        estimates = []
        if estimate_outcome:
            by_symbol = {r.symbol: r for r in ranked}
            for g in grades:
                est = outcome_mod.estimate(
                    g.drug_name, g.tf_symbol, prom.gene_symbol, goal=goal,
                    grade=g, tf_result=by_symbol.get(g.tf_symbol))
                d = est.to_dict()
                d["verdict"] = est.verdict
                d["next_experiment"] = outcome_mod.narrowing_experiment(est)
                estimates.append(d)

            # Order by the chain's own result. Inheriting grade order meant the
            # discovery winner could fall outside the printed rows while worse
            # candidates were shown — two rankings disagreeing silently.
            estimates.sort(key=lambda e: -e["p_point"])

        # ------------------------------------- 8. optional library screen
        screens = {}
        if screen_receptors and not demo:
            from . import screening
            seen_targets = set(_screen_cache)
            for r in ranked:
                for drug in r.drugs:
                    tid = getattr(drug, "target_chembl_id", None)
                    if tid and tid not in seen_targets:
                        seen_targets.add(tid)
                        if not _library:
                            _library.extend(_optional(
                                "approved drug library",
                                lambda: screening.approved_library(client), []))
                        sc = _optional(
                            f"virtual screen for {r.symbol}",
                            lambda t=tid: screening.screen_receptor(
                                client, t, goal=goal,
                                library=_library or None,
                                gene_linked=_gene_linked), None)

                        # A receptor with too few ligands of its own cannot
                        # train a model, which is a data gap rather than a
                        # fact about the receptor. `screen_family` borrows
                        # its relatives' chemistry and MEASURES whether that
                        # chemistry transfers before reporting anything. It
                        # existed for this and was called from nowhere.
                        if sc and sc.get("error") and not sc.get("hits"):
                            fam = _optional(
                                f"family screen for {r.symbol}",
                                lambda sym=r.symbol: screening.screen_family(
                                    client, sym, goal=goal,
                                    library=_library or None,
                                    gene_linked=_gene_linked), None)
                            if fam and fam.get("hits"):
                                fam["replaces"] = sc.get("error")
                                sc = fam
                        if sc:
                            _screen_cache[tid] = sc
                            screens[r.symbol] = sc
                        break
                    elif tid in _screen_cache:
                        # Already screened for an earlier TSS of this gene —
                        # the receptor has not changed, so reuse the result
                        # instead of repeating the whole screen and printing
                        # it again.
                        screens.setdefault(r.symbol, _screen_cache[tid])
                        break

        # ------------------- 8a-ter. the promoter sequence itself
        # The only route that targets no protein. Gated on the sequence
        # beating its own dinucleotide shuffles, because a GC-rich
        # promoter produces quadruplexes and GC-boxes by composition and
        # proposing an agent for those would be proposing one for most of
        # the genome.
        promoter_dna_scan: dict = {}
        if promoter_dna and prom.sequence:
            scan = _optional(
                "promoter DNA feature scan",
                lambda: dna_mod.scan_features(prom.sequence), {}) or {}
            if scan and not scan.get("error"):
                cands = dna_mod.candidates(scan)
                promoter_dna_scan = dict(scan)
                promoter_dna_scan["candidates"] = cands
                promoter_dna_scan["summary"] = dna_mod.summarise(
                    scan, cands)
                n_enr = sum(1 for f in (scan.get("features") or {}).values()
                            if f.get("enriched"))
                ledger.record(
                    "promoter_dna", requested=True, ok=bool(cands),
                    yielded=len(cands),
                    detail=(f"{n_enr} sequence feature(s) enriched over "
                            f"{scan.get('shuffles')} shuffles; "
                            f"{len(cands)} agent(s) prefer them"
                            if cands else
                            "no sequence feature beat the shuffled null, "
                            "so no DNA-directed agent is proposed"))
            else:
                ledger.record("promoter_dna", requested=True, ok=False,
                              yielded=0,
                              detail=str(scan.get("error") or
                                         "no promoter sequence"))

        # ------------------- 8a-bis. the gene's OWN product
        # Every screen above points at a transcription factor. This one
        # points at the gene's own protein, which `direct` only ever
        # reached by looking up what was already annotated - so the
        # package could report Maraviroc for CCR5 and could not propose a
        # compound that might bind CCR5 and is recorded against nothing.
        product_screen: dict = {}
        if screen_product and not demo:
            from . import screening as _scr
            product_screen = _optional(
                f"gene-product screen for {prom.gene_symbol}",
                lambda: _scr.screen_gene_product(
                    client, prom.gene_symbol, uniprot=uniprot,
                    goal=goal, library=_library or None),
                {}) or {}
            ledger.record(
                "product_screen", requested=True,
                ok=bool(product_screen.get("hits")),
                yielded=len(product_screen.get("hits") or []),
                detail=(f"{product_screen.get('n_novel', 0)} unannotated "
                        f"candidate(s) at the {prom.gene_symbol} protein, "
                        f"cross-validated AUC "
                        f"{product_screen.get('model_auc')}"
                        if product_screen.get("hits") else
                        str(product_screen.get("error") or
                            "no result")[:200]))

        # ------------------- 8b. orthosteric vs allosteric at the top targets
        mechanisms = {}
        if compare_mechanisms and not demo:
            for r in ranked[:5]:
                if not r.drugs:
                    continue
                try:
                    mod_hits = mod_mod.find_modulators(
                        client, r.symbol, r.annotation.get("uniprot"), goal=goal,
                        approved_only=approved_only, use_dgidb=use_dgidb)
                except Exception as exc:                      # noqa: BLE001
                    log.warning("modulator search failed for %s: %s", r.symbol, exc)
                    continue
                if mod_hits:
                    cmp = mod_mod.compare_mechanisms(
                        mod_hits, ligand_axis=lig_mod.has_endogenous_ligand(
                            r.symbol))
                    cmp["hits"] = [h.to_dict() for h in mod_hits]
                    mechanisms[r.symbol] = cmp
        elif compare_mechanisms and demo:
            for r in ranked[:5]:
                if not r.drugs:
                    continue
                mod_hits = [mod_mod.ModulatorHit(
                    drug=d.name, target=r.symbol, mechanism=mod_mod.classify(d),
                    action_type=d.action_type, direction=d.direction,
                    approved=d.approved, withdrawn=d.withdrawn,
                    first_approval=d.first_approval,
                    is_allosteric=mod_mod.classify(d) in mod_mod.ALLOSTERIC,
                    profile=mod_mod.profile_for(mod_mod.classify(d)).to_dict(),
                ) for d in r.drugs]
                cmp = mod_mod.compare_mechanisms(
                    mod_hits, ligand_axis=lig_mod.has_endogenous_ligand(r.symbol))
                cmp["hits"] = [h.to_dict() for h in mod_hits]
                mechanisms[r.symbol] = cmp

        # ------------------- 8c. multi-hop upstream routes to undruggable TFs
        upstream_hits = {}
        if upstream_network and not demo:
            edges = net_mod.fetch_network(client)
            # The high-scoring factors with no ligand are precisely what this
            # search exists for, so they go in first rather than depending on
            # where they happen to fall in the ranking.
            _priority = [r for r in ranked if r.symbol in
                         set(attrition.get("top25_without_drug") or [])]
            _rest = [r for r in ranked[:6] if r not in _priority]
            for r in (_priority + _rest)[:8]:
                # The point of this step is TFs with no direct pharmacology.
                if r.drugs:
                    continue
                try:
                    hs = net_mod.find_upstream_targets(
                        client, r.symbol, goal=goal,
                        tf_role=r.annotation.get("role", "activator"),
                        max_depth=upstream_depth, edges=edges)
                    hs = net_mod.attach_drugs(client, hs,
                                              approved_only=approved_only)
                except Exception as exc:                      # noqa: BLE001
                    log.warning("upstream search failed for %s: %s", r.symbol, exc)
                    continue
                hs = [h for h in hs if h.drugs]
                if hs:
                    upstream_hits[r.symbol] = [h.to_dict() for h in hs[:10]]
            ledger.record("upstream", requested=True, ok=bool(upstream_hits),
                          yielded=len(upstream_hits),
                          detail=(f"druggable routes for "
                                  f"{len(upstream_hits)} factor(s)"
                                  if upstream_hits else
                                  "no signed causal path reached a druggable "
                                  "node"))

        # ------------- 8b2. stabilisation routes for pocket-less factors
        stabilisation = {}
        if include_stabilizers:
            for r in ranked[:8]:
                a = stab_mod.assess(r.symbol, goal=goal,
                                    annotation=r.annotation,
                                    has_direct_drug=bool(r.drugs))
                if not demo:
                    extra = stab_mod.chembl_stabilisers(
                        client, r.symbol, r.annotation.get("uniprot"),
                        goal=goal, approved_only=approved_only)
                    a["stabilisers"] += [e.to_dict() for e in extra]
                # Keep the no-route factors too. They are the majority, and
                # dropping them silently means a reader never learns that
                # stabilisation was considered for them and found nothing —
                # which is exactly the question "can we stabilise this?"
                if (a["stabilisers"] or a["recommended"]
                        or not a["has_ligand_pocket"]):
                    stabilisation[r.symbol] = a

        # ------------- 8b3. change the natural ligand instead of the receptor
        ligand_axes = {}
        if include_stabilizers:      # same "beyond direct ligands" switch
            for r in ranked[:8]:
                a = lig_mod.assess(r.symbol, goal=goal)
                if a["routes"]:
                    ligand_axes[r.symbol] = a

        # ------------- 8c2. does the drug reach the cells expressing the gene?
        celltype = {}
        if cell_type_check and not demo:
            gene_tissues = _optional(
                "tissue expression lookup",
                lambda: dlv_mod.tissue_expression(client, prom.gene_symbol),
                set())
            expr_status = dlv_mod.expression_source_status()
            # A gene with no tissues AND a working source is a real result;
            # with a broken source it is nothing at all.
            ledger.record("expression", requested=expr_status["requested"],
                          ok=expr_status["ok"] and bool(gene_tissues),
                          yielded=len(gene_tissues),
                          detail=(expr_status["detail"] if not expr_status["ok"] else
                                  f"{len(gene_tissues)} tissue(s) for "
                                  f"{prom.gene_symbol}" if gene_tissues else
                                  f"source available but returned no tissue "
                                  f"for {prom.gene_symbol} — co-expression "
                                  f"cannot be checked for any factor"))
            celltype["gene_tissues"] = sorted(gene_tissues)[:20]
            # The denominator a tissue count needs to mean anything. None
            # when no local table was loaded — see tissue_universe_size.
            celltype["n_tissues_total"] = dlv_mod.tissue_universe_size()
            celltype["tfs"] = {}
            for r in ranked[:10]:
                tf_t = dlv_mod.tissue_expression(client, r.symbol)
                ov = dlv_mod.expression_overlap(gene_tissues, tf_t)
                celltype["tfs"][r.symbol] = ov
                if ov.get("compatible") is False:
                    r.score = round(r.score * 0.2, 1)
            ranked.sort(key=lambda r: -r.score)

            celltype["delivery"] = []
            for g in grades[:20]:
                a = dlv_mod.assess_delivery(
                    g.drug_name, compartment=(tissue or "unknown"),
                    properties=(g.evidence or {}).get("properties"))
                if a.multiplier != 1.0 or a.notes:
                    celltype["delivery"].append(a.to_dict())

        # ---- how many other genes does each factor control?
        # The collateral question, and the reason a high-scoring route can
        # still be a bad idea: a drug acting through a factor moves
        # everything that factor controls, not only your gene. Counted here
        # for every ranked factor rather than only the ones that acquired a
        # drug, so the column is populated even when the drug lookup is
        # empty.
        for r in ranked:
            n, prov = regulon_mod.regulon_size(r.symbol, regulons)
            r.n_regulon, r.regulon_provenance = n, prov
            # A low curated count usually means nobody looked, not that the
            # factor is selective. regulon.py already draws this line; the
            # flag carries it onto the row so the table can show it.
            r.regulon_understudied = bool(
                n is not None and prov == "regulon file"
                and n < regulon_mod.MIN_CREDIBLE_REGULON)

        # ---- what does the cell lose without this factor?
        # The regulon count says how many genes a factor controls, from
        # literature curation. This says what a cell cannot do without it,
        # measured by CRISPR knockout across ~1,100 lines. Same collateral
        # question, answered rather than approximated - and it covers the
        # factors the curated tables leave at "not counted".
        dependencies: dict = {}
        _dm_path = depmap_mod.resolve(depmap_file)
        if _dm_path:
            # Every ranked factor, not a slice: the file is read in one
            # streaming pass whatever is asked for, so narrowing the
            # request saves nothing and loses the factors further down -
            # which are exactly the ones with no curated regulon.
            wanted = [r.symbol for r in ranked]
            deps = _optional(
                "DepMap dependency",
                lambda: depmap_mod.assess_many(wanted, effect_path=_dm_path),
                {})
            dependencies = {k: v.to_dict() for k, v in (deps or {}).items()}
            for r in ranked:
                d = (deps or {}).get(r.symbol.upper())
                if d is not None:
                    r.dependency = d.verdict
                    r.dependency_burden = d.burden
            ledger.record(
                "depmap", requested=True, ok=True, yielded=len(dependencies),
                detail=(f"{len(dependencies)} of {len(wanted)} factor(s) "
                        f"scored in {Path(_dm_path).name}"
                        if dependencies else
                        f"no requested factor appears in "
                        f"{Path(_dm_path).name}"))
        elif depmap_file:
            ledger.record("depmap", requested=True, ok=False, yielded=0,
                          detail=f"{depmap_file} not found")

        # ---- 8c-ter. which approved drugs behave like losing the factor
        # Every other route here begins with something somebody wrote
        # down, so none of them can surface a drug nobody has connected to
        # the factor. This one reads a measured phenotype instead: if a
        # compound hurts the same cell lines that knocking the factor out
        # hurts, it behaves as though it removes the factor.
        sensitivity: dict = {}
        if sensitivity_scan and not demo:
            _prism = sens_mod.load_prism(prism_file)
            if not _prism:
                ledger.record(
                    "sensitivity", requested=True, ok=False, yielded=0,
                    detail=("no PRISM matrix on disk, so no drug was "
                            "correlated. " + sens_mod.SOURCE.splitlines()[0]))
            else:
                _names = sens_mod.load_compound_names()
                for r in ranked[:SENSITIVITY_FACTORS]:
                    got = _optional(
                        f"sensitivity scan for {r.symbol}",
                        lambda sym=r.symbol: sens_mod.phenocopy(
                            sym, prism=_prism, effect_path=depmap_file,
                            names=_names),
                        None)
                    if got is not None:
                        sensitivity[r.symbol] = got.to_dict()
                n_hits = sum(len(v.get("hits") or [])
                             for v in sensitivity.values())
                n_ok = sum(1 for v in sensitivity.values()
                           if v.get("assessed"))
                ledger.record(
                    "sensitivity", requested=True, ok=bool(n_hits),
                    yielded=n_hits,
                    detail=(f"{n_hits} correlated compound(s) across "
                            f"{n_ok} of {len(sensitivity)} factor(s) "
                            f"scanned"
                            if n_ok else
                            "no factor could be scanned: none is scored in "
                            "both the DepMap and PRISM panels"))

        # ---- 8c-bis. is the factor actually ACTIVE, not merely motif-matched?
        # The independent check on everything above. A PWM scan of a 6 kb
        # window returns mostly noise at p=1e-4, so a site is evidence that
        # a factor COULD bind. This asks whether its known targets actually
        # moved - a different question, from different data, and the two
        # agreeing is worth far more than either alone.
        activity: dict = {}
        if signed_regulons and activity_tissue:
            acts = _optional(
                "TF activity inference",
                lambda: tfa_mod.tissue_activity(
                    activity_tissue, signed_regulons,
                    permutations=activity_permutations), [])
            scanned = [r.symbol for r in ranked]
            scanned_up = {s.upper() for s in scanned}
            by_tf = {a.tf.upper(): a for a in acts}

            activity = {
                "tissue": activity_tissue,
                "n_testable": len(acts),
                "n_regulons": len(signed_regulons),
                "regulon_source": regulon_source,
                "permutations": activity_permutations,
                # Only the factors this run actually found sites for, plus
                # the overall top calls for context.
                "scanned": [a.to_dict() for a in acts
                            if a.tf.upper() in scanned_up],
                "top": [a.to_dict() for a in acts[:15]],
                "summary": tfa_mod.summarise(acts, scanned=scanned),
            }

            # Annotate every scanned factor, and say so rather than silently
            # moving the score: the activity call is a reference-tissue
            # average and should inform a reader, not overrule a scan.
            for r in ranked:
                a = by_tf.get(r.symbol.upper())
                if a is None:
                    r.annotation.setdefault(
                        "activity", "not testable - too few curated targets "
                                    "in the regulon table to infer activity")
                    continue
                r.annotation["activity"] = a.direction
                r.annotation["activity_nes"] = a.nes
                r.annotation["activity_q"] = a.q_value

            ledger.record("tf_activity", requested=True, ok=bool(acts),
                          yielded=len(acts),
                          detail=(f"{len(acts)} factor(s) testable in "
                                  f"{activity_tissue}; "
                                  f"{len(activity['scanned'])} of them have a "
                                  f"site in this promoter"
                                  if acts else
                                  f"no factor had enough curated targets "
                                  f"measured in {activity_tissue}"))
        elif tf_activity:
            # Name only what is actually missing. The fixed string here
            # listed both flags, so a run that passed --regulon-file was
            # told to pass --regulon-file, which sends you to check the
            # one thing that was already right.
            missing = []
            if not signed_regulons:
                missing.append("a regulon table (--regulon-file, or put "
                               "CollecTRI/TRRUST in data/)")
            if not activity_tissue:
                missing.append("--activity-tissue")
            ledger.record(
                "tf_activity", requested=True, ok=False, yielded=0,
                detail=("--tf-activity needs " + " and ".join(missing)
                        + "; motif scores are unchecked"
                        if missing else
                        "--tf-activity found no testable factor; motif "
                        "scores are unchecked"))

        # ---- which drugs change how much of the factor EXISTS?
        # Everything else here acts on a factor's activity, which needs it
        # to be druggable. The corpora already record what moves ~20,000
        # transcripts, and the pipeline only ever asked that about the
        # target gene. Asking it about the factors reaches the ones with
        # no ligand at all.
        factor_abundance: list = []
        undirected_factors: list = []
        # `not demo` matters: demo mode must stay self-contained and
        # reproducible, and without it this silently read whatever real
        # corpora happened to be sitting in data/.
        if not demo and sig_mod.autoload_gmt_available("up") \
                and sig_mod.autoload_gmt_available("down"):
            rows = [r.to_dict() for r in ranked[:ACTIVITY_DRUG_FACTORS]]
            found = _optional(
                "factor abundance scan",
                lambda: abund_mod.find(rows, goal=goal), [])
            factor_abundance = [h.to_dict() for h in (found or [])]
            undirected_factors = abund_mod.undirected_factors(rows, goal=goal)
            ledger.record(
                "factor_abundance", requested=True, ok=True,
                yielded=len(factor_abundance),
                detail=(f"{len(factor_abundance)} compound(s) move a "
                        f"factor's own transcript the needed way; "
                        f"{len(undirected_factors)} factor(s) skipped for "
                        f"an unrecorded role"
                        if factor_abundance else
                        "no compound moved a factor's own transcript the "
                        "needed way in enough contexts to report"))

        # ---- 8c-ter. which approved drugs move those factors?
        # The drug lookup can only return compounds somebody assayed at the
        # target, so a factor with no medicinal chemistry behind it comes
        # back empty however important it is. This asks the LINCS data
        # already on disk instead: whose own perturbation signature shifts
        # this factor's regulon? It needs no annotated binding, so it can
        # reach factors the other sources cannot.
        activity_drugs: dict = {}
        corpus_scans: dict = {}
        corpus_agreement: dict = {}
        if signed_regulons and any(c.available()
                                   for c in corpora_mod.CORPORA):
            # Capped: each factor is a full pass over every signature in
            # every corpus present, and the ranking below the top few is
            # not worth the wall clock.
            for r in ranked[:ACTIVITY_DRUG_FACTORS]:
                reg = signed_regulons.get(r.symbol.upper())
                if not reg:
                    continue
                scans = _optional(
                    f"perturbation corpora scan for {r.symbol}",
                    lambda g=reg, s=r.symbol: corpora_mod.scan(s, g), {})
                if not scans:
                    continue
                corpus_scans[r.symbol] = {k: v.to_dict()
                                          for k, v in scans.items()}
                agreed = corpora_mod.agreement(scans)
                if agreed:
                    corpus_agreement[r.symbol] = agreed
                # Pooled for the route tables, with the corpus stamped on
                # each hit so a rat liver result cannot be mistaken for a
                # human cell-line one.
                found = []
                for key, res in scans.items():
                    for h in res.hits:
                        h = dict(h)
                        h["corpus"] = key
                        found.append(h)
                if found:
                    found.sort(key=lambda h: (-h["n_significant"],
                                              -abs(h["median_z"])))
                    activity_drugs[r.symbol] = found
            # ok=True even when nothing was found. The scan RAN; finding
            # no compound is a negative result, not a failed lookup, and
            # conflating the two is the thing provenance.py exists to stop.
            ledger.record(
                "activity_drugs", requested=True, ok=True,
                yielded=sum(len(v) for v in activity_drugs.values()),
                detail=(f"{sum(len(v) for v in activity_drugs.values())} "
                        f"compound(s) across {len(activity_drugs)} factor(s) "
                        f"from "
                        f"{', '.join(c.key for c in corpora_mod.CORPORA if c.available())}"
                        f"; {sum(len(v) for v in corpus_agreement.values())} "
                        f"agree across more than one corpus"
                        if activity_drugs else
                        "the available corpora were scanned; no compound "
                        "moved any testable factor's regulon consistently "
                        "enough to report"))

        # ------------------ 8d. empirical: did any drug actually move this gene?
        lincs = {}
        # Supplying any LINCS input enables the step. `lincs_gmt` used
        # to be the exception: the files were loaded and the step
        # never ran, so a run with --lincs-gmt reported no LINCS row
        # at all - not even that it had been asked for.
        if (use_lincs or lincs_file or lincs_gmt) and not demo:
            sig_hits = _optional(
                "LINCS signature query",
                lambda: (sig_mod.load_signature_table(
                             lincs_file, prom.gene_symbol, direction=goal)
                         if lincs_file else
                         # A local Enrichr .gmt beats SigCom outright: it
                         # works, it is reproducible, and SigCom has been
                         # returning 500 for days at a time.
                         sig_mod.load_enrichr_gmt(prom.gene_symbol,
                                                  direction=goal)
                         if sig_mod.autoload_gmt_available(goal) else
                         sig_mod.query_gene(client, prom.gene_symbol,
                                            direction=goal)),
                [])
            sig_hits = sig_mod.flag_non_specific(sig_hits)
            ledger.record(
                "lincs", requested=True, ok=bool(sig_hits), yielded=len(sig_hits),
                detail=(f"{len(sig_hits)} signature hit(s) for this gene"
                        if sig_hits else
                        "LINCS returned no signature for this gene. Either no "
                        "approved compound has moved it in the profiled cell "
                        "lines, or the query failed — check the gene is "
                        "expressed in the L1000 panel before reading this as "
                        "a negative."))
            if signature_file:
                ref = sig_mod.load_signature_file(signature_file)
                lincs["reference_signature_genes"] = len(ref)
            lincs["hits"] = [h.to_dict() for h in sig_hits[:25]]
            lincs["summary"] = sig_mod.summarise(sig_hits, prom.gene_symbol)

        # ------------------------------- 8e. human genetic evidence
        genetic = {}
        if use_genetics and not demo:
            # A local GWAS catalogue first: on a blocked network the
            # Open Targets call is what empties this section.
            _local_gen = gen_mod.local_support(prom.gene_symbol)
            cis = _optional("Open Targets genetics",
                            lambda: gen_mod.target_consequences(
                                client, prom.gene_symbol),
                            gen_mod.GeneticEvidence(prom.gene_symbol))
            genetic["cis"] = cis.to_dict()
            genetic["summary"] = gen_mod.summarise(cis)
            _phen = (genetic["cis"] or {}).get("phenotypes") or []
            if _local_gen:
                genetic["local_gwas"] = _local_gen
            ledger.record("genetics", requested=True,
                          # A catalogue that was read and found nothing is
                          # assessed, not successful. `ok` means evidence
                          # arrived, so zero hits leaves it False and the
                          # detail below says which source was consulted.
                          ok=bool(_phen) or _n_local_gwas(_local_gen) > 0,
                          yielded=len(_phen) + _n_local_gwas(_local_gen),
                          detail=_genetics_detail(_phen, _local_gen))
            genetic["trans"] = [
                gen_mod.tf_link_support(client, r.symbol, prom.gene_symbol).to_dict()
                for r in ranked[:3]]

            # Actually run the trans-eQTL lookup the advice above describes.
            # The eQTLGen significant file is 2 MB; telling someone to go
            # check it by hand, once per factor, in identical wording, was
            # advice rather than an answer.
            def _trans_eqtl():
                any_hit = eqtl_mod.gene_has_any_trans_eqtl(
                    prom.gene_symbol, client=client)
                out = []
                for r in ranked[:6]:
                    locus = None
                    if any_hit:
                        # Only worth an Ensembl lookup when there is
                        # something for it to be near.
                        info = _optional(
                            f"locus of {r.symbol}",
                            lambda s=r.symbol: genes._lookup_symbol(
                                client, s, species), None)
                        if info and info.get("seq_region_name"):
                            locus = (info["seq_region_name"], info["start"],
                                     info["end"])
                    out.append(eqtl_mod.trans_support(
                        r.symbol, prom.gene_symbol, tf_locus=locus,
                        client=client))
                return out

            supports = _optional("eQTLGen trans-eQTL lookup", _trans_eqtl, [])
            if supports:
                genetic["trans_eqtl"] = [s.to_dict() for s in supports]
                genetic["trans_eqtl_summary"] = eqtl_mod.summarise(
                    prom.gene_symbol, supports)
                ledger.record(
                    "trans_eqtl", requested=True,
                    ok=any(s.checked for s in supports),
                    yielded=sum(1 for s in supports if s.supported),
                    detail=(f"{sum(1 for s in supports if s.supported)} of "
                            f"{len(supports)} factor(s) have trans-eQTL "
                            f"support" if any(s.checked for s in supports)
                            else "eQTLGen file not retrieved"))

        # ------- 8e0. approved drugs resembling our own best candidates
        # Asked of the report's strongest rows, and with everything already
        # annotated at the factor subtracted - otherwise a similarity search
        # rediscovers the table it was given and looks productive.
        analog_hits = {}
        if find_analogs and not demo and grades:
            analog_hits = _optional(
                "chemical-analog search",
                lambda: analog_mod.search(client, [g.to_dict() for g in grades]),
                {})

        # ------- 8e1. the routes a promoter-level analysis cannot see
        # Everything above acts on transcription. These act after it, around
        # it, or on a different gene entirely - and for some genes that is
        # where the whole problem is.
        beyond = {}
        if goal == "up":
            beyond["post_transcriptional"] = {
                "routes": [r.to_dict() for r in
                           ptx_mod.routes_for(goal, has_ptc=has_ptc)],
                "summary": ptx_mod.summarise(
                    ptx_mod.routes_for(goal, has_ptc=has_ptc),
                    has_ptc=has_ptc),
            }
            comp = comp_mod.assess(prom.gene_symbol, goal=goal,
                                   client=None if demo else client)
            comp["summary"] = comp_mod.summarise(comp)
            beyond["compensation"] = comp

            react = react_mod.assess(prom.gene_symbol, goal=goal)
            react["summary"] = react_mod.summarise(react)
            beyond["reactivation"] = react

            beyond["product_stability"] = stab_mod.product_stabilisation(
                prom.gene_symbol, goal=goal)

        # ---------------- 8e2. design a direct modality for this gene
        modality = {}
        if design_modality:
            cdna = (None if demo else
                    _optional("cDNA fetch",
                              lambda: genes.fetch_cdna(client,
                                                       prom.transcript_id),
                              None))
            modality = modal_mod.design(
                prom.gene_symbol, goal=goal, tissue=tissue,
                transcript=cdna, promoter=prom)

        # ---------------- 8f. design the experiment that settles it
        screen_plan = {}
        if design_screen:
            screen_plan = screens_mod.design(
                prom, [r.symbol for r in ranked[:20]],
                [g.to_dict() for g in grades],
                [t.to_dict() for t in direct_hits]).to_dict()

        # -------------------------- 9. what beats the current standard of care
        improvements = {}
        if find_improvements:
            improvements = discovery.discover(
                client, prom.gene_symbol, grades, ranked, direct_hits,
                outcomes=estimates, goal=goal, demo=demo)
            # CTD is empirical, not annotation: a curator recorded that
            # someone dosed cells and this gene's expression moved. It was
            # only ever used for the novelty check, so a run could report
            # "no measured evidence" while holding 152 curated expression
            # effects — which is what happened on the real CCR5 run.
            if not demo:
                n_ctd = improvements.get("ctd_links_found") or 0
                rl = (improvements.get("reported_links") or {}).get("counts") or {}
                ledger.record(
                    "ctd", requested=True, ok=bool(n_ctd), yielded=n_ctd,
                    detail=(f"{n_ctd} chemical(s) with curated expression "
                            f"effects on this gene; {rl.get('wanted', 0)} move "
                            f"it the way you want" if n_ctd else
                            "no curated chemical-gene expression effects "
                            "returned for this gene"))

        # Saturation check: when every surviving route comes from one or
        # two receptor families, that is the signature of the attrition
        # above rather than a finding about this gene.
        from . import paralogs as _par
        _fams = {}
        for g in grades:
            f = _par.family_of(g.tf_symbol or "")
            if f:
                _fams.setdefault(f.name, set()).add(g.tf_symbol)
        attrition["receptor_families"] = {k: sorted(v)
                                          for k, v in _fams.items()}
        attrition["family_saturated"] = bool(grades) and len(_fams) <= 2             and len(_fams) > 0

        block_for_catalogue = {
            "direct_therapies": [t.to_dict() for t in direct_hits],
            "drug_grades": [g.to_dict() for g in grades],
            "mechanisms": mechanisms,
            "stabilisation": stabilisation,
            "ligand_axes": ligand_axes,
            "upstream": upstream_hits,
            "lincs": lincs,
            "screens": screens,
            # Both of these were computed above and left out of the dict
            # the catalogue actually reads, so their routes never reached
            # the tables even though their report sections rendered.
            "activity_drugs": activity_drugs,
            "factor_abundance": factor_abundance,
            "undirected_factors": undirected_factors,
            "docking": list(docking_runs.values()),
            # The factors at this promoter and their roles. The docking
            # collector needs both: to tell a dock against one of them
            # from a dock against something else, and to say whether the
            # conformation docked is the one you want stabilised.
            # `annotation` is a dict subclass, so this is a key read.
            "tf_roles": {r.symbol: (r.annotation or {}).get("role")
                         for r in ranked},
            "sensitivity": sensitivity,
            "product_screen": product_screen,
            "promoter_dna": promoter_dna_scan,
            "gene_symbol": prom.gene_symbol,
        }
        all_routes = cat_mod.build(block_for_catalogue, goal=goal)

        # Where a factor has both a dock and a ligand screen, intersect
        # them. The two share no input - one reads pocket shape, the other
        # reads chemistry measured at the receptor - so agreement is worth
        # more than a better score from either. Both were computed and
        # never compared.
        dock_consensus = {}
        if docking_runs and screens:
            from . import docking as _dk
            for sym, sc in screens.items():
                run = docking_runs.get(sym)
                if not run:
                    continue
                hits = (sc.get("hits") if isinstance(sc, dict)
                        else getattr(sc, "hits", None)) or []
                both = _dk.consensus(run, hits)
                if both:
                    dock_consensus[sym] = both

        if screen_receptors and not demo:
            ledger.record("screen", requested=True, ok=bool(screens),
                          yielded=len(screens),
                          detail=(f"{len(screens)} receptor(s) screened"
                                  if screens else "no receptor screened"))

        outputs.append({
            "promoter": prom.to_dict(),
            "sequence": prom.sequence,
            "n_sites": len(hits),
            "scan_quality": quality,
            "salt_forms_collapsed": collapsed if grade_drugs else {},
            "attrition": attrition,
            "beyond_transcription": beyond,
            "analogs": analog_hits,
            "evidence_verdict": prov_mod.verdict(ledger, quality),
            "neighbours": neighbours,
            "n_tfs": len(by_tf),
            "tfs": [r.to_dict() for r in ranked],
            "direction_unknown": {
                sym: [d.to_dict() for d in ds]
                for sym, ds in (unknown_direction or {}).items()},
            "cell_cycle": cycle,
            "cell_cycle_summary": cycle_summary,
            "cell_cycle_blocked": sorted(blocked),
            "drug_grades": [g.to_dict() for g in grades],
            "outcome_estimates": estimates,
            "improvements": improvements,
            "all_routes": all_routes,
            "mechanisms": mechanisms,
            "upstream": upstream_hits,
            "stabilisation": stabilisation,
            "ligand_axes": ligand_axes,
            "cell_type": celltype,
            "tf_activity": activity,
            "activity_drugs": activity_drugs,
            "factor_abundance": factor_abundance,
            "undirected_factors": undirected_factors,
            "corpus_scans": corpus_scans,
            "corpus_agreement": corpus_agreement,
            "docking": list(docking_runs.values()),
            "docking_consensus": dock_consensus,
            "dependencies": dependencies,
            "sensitivity": sensitivity,
            "product_screen": product_screen,
            "promoter_dna": promoter_dna_scan,
            "modality_design": modality,
            "screen_plan": screen_plan,
            "lincs": lincs,
            "genetics": genetic,
            "screens": screens,
            "hits": [vars(h) for h in hits],
        })

    # How much the --all-tss windows duplicate each other. CCR5's two
    # overlapped by 4275 of 6000 bp, so they were largely one analysis
    # counted twice — which inflates every per-gene route count downstream.
    tss_overlap = nbr_mod.window_overlap(
        [o["promoter"] for o in outputs]) if len(outputs) > 1 else []

    return {
        "direct_therapies": [t.to_dict() for t in direct_hits],
        "tss_window_overlap": tss_overlap,
        "evidence_sources": ledger.to_dict(),
        "evidence_verdict": prov_mod.verdict(
            ledger, (outputs[0] or {}).get("scan_quality") if outputs else None),
        "query": {
            "gene": gene, "species": species, "upstream": upstream,
            "downstream": downstream, "pvalue": pvalue, "demo": demo,
            "mask_repeats": mask_repeats, "goal": goal,
            "tissue": tissue, "tf_weight": tf_weight,
            "proliferation": proliferation, "cell_type_check": cell_type_check,
            "use_lincs": use_lincs,
            "use_genetics": use_genetics,
            "regulon_file": regulon_file, "include_direct": include_direct,
            "screen_receptors": screen_receptors,
            "compare_mechanisms": compare_mechanisms,
            "upstream_network": upstream_network,
            "find_improvements": find_improvements,
            "chip_bed": chip_bed, "chip_api": chip_api,
            "approved_only": approved_only,
            "n_motifs": len(motif_set),
        },
        "results": outputs,
        "weights": WEIGHTS,
        "drug_weights": drug_scoring.DRUG_WEIGHTS,
    }


# ------------------------------------------------------------------ reports
def _fmt_offset(o) -> str:
    return f"{o:+d}" if isinstance(o, int) else "?"


def _dep_label(tf: dict) -> str:
    """The DepMap verdict, shortened for a table cell."""
    v = tf.get("dependency")
    if not v:
        return "not scored"
    return {
        "essential in nearly every line": "yes, nearly all",
        "essential in most lines": "yes, most",
        "essential in a substantial minority": "some lines",
        "essential in a minority of lines": "some lines",
        "essential in a handful of lines": "a few lines",
        "dispensable across the panel": "no",
    }.get(v, v)


def _corpus_note(res: dict) -> list:
    """Which corpora ran, which did not, and what agreed across them."""
    from . import corpora as cp

    scans = res.get("corpus_scans") or {}
    agree = res.get("corpus_agreement") or {}
    if not scans:
        return []

    present, absent = set(), set()
    for per_tf in scans.values():
        for key, r in per_tf.items():
            (present if r.get("available") else absent).add(key)

    L = []
    if len(present) > 1:
        L.append(f"**Scanned across {len(present)} independent corpora** "
                 f"({', '.join(sorted(present))}). They differ in species, "
                 f"tissue, platform and dose, so they share few systematic "
                 f"errors - which is the point of running more than one. "
                 f"Every signature in a single corpus shares that corpus's "
                 f"biases, so no amount of internal replication detects "
                 f"them.")
    else:
        only = ", ".join(sorted(present)) or "none"
        L.append(f"**Only one corpus was available ({only}).** Every "
                 f"signature in it shares its platform, species and dose "
                 f"regime, so internal replication cannot detect a "
                 f"systematic error common to all of them. A second corpus "
                 f"is what would.")
    if absent:
        L.append("")
        L.append(f"Not on disk, and therefore skipped rather than searched: "
                 f"{', '.join(sorted(absent))}. Absent is not the same as "
                 f"finding nothing. Run "
                 f"`python -m promoter_tf_drugs.corpora` for the downloads.")
    L.append("")

    flat = [a for per_tf in agree.values() for a in per_tf]
    if flat:
        consistent = [a for a in flat if a["consistent"]]
        conflicting = [a for a in flat if not a["consistent"]]
        if consistent:
            L.append("#### Agreeing across corpora")
            L.append("")
            L.append("The strongest empirical result this package can "
                     "produce. These compounds move the factor the same "
                     "way in systems that share almost no methodology.")
            L.append("")
            L.append("| Drug | Corpora | Direction |")
            L.append("|---|---|---|")
            for a in consistent[:15]:
                L.append(f"| {a['drug']} | "
                         f"{', '.join(sorted(a['corpora']))} | "
                         f"{a['direction']} |")
            L.append("")
        if conflicting:
            L.append(f"**{len(conflicting)} compound(s) appear in more than "
                     f"one corpus with opposite directions** "
                     f"({', '.join(a['drug'] for a in conflicting[:6])}). "
                     f"That usually means the systems are not comparable - "
                     f"rat liver in vivo against a human cell line, at "
                     f"different doses - rather than that either reading is "
                     f"wrong. It is shown because hiding it would make the "
                     f"agreeing rows look more uniform than the data is.")
            L.append("")
    elif len(present) > 1:
        L.append("No compound moved a regulon consistently in more than one "
                 "corpus. With this little overlap between their compound "
                 "sets that is the expected outcome, not a negative result.")
        L.append("")

    for key in sorted(present):
        c = cp.BY_KEY.get(key)
        if c and c.species != "human":
            L.append(f"*{c.name} is {c.species} ({c.system}). Gene symbols "
                     f"cross a species boundary and the regulon is not "
                     f"guaranteed to be conserved, so treat a hit found "
                     f"only there as a lead to check in a human system.*")
            L.append("")
    return L


def _docking_section(runs: list, consensus: dict) -> list:
    """Docking validation, and the compounds two independent methods both
    picked.

    The routes from a dock already appear in the tables. What does not
    appear anywhere is whether the pocket was validated at all, which is
    the only thing that separates a docking score from a random number.
    A screen that cannot separate known actives from decoys has not been
    shown to rank anything, and its hits are reported as unvalidated
    rather than withheld.
    """
    if not runs:
        return []
    from .docking import MIN_DOCKING_AUC, MIN_EF1

    L = ["### Structure-based docking", ""]
    L.append("Docking asks whether a molecule fits the pocket. It says "
             "nothing about direction: a compound that binds a factor may "
             "activate it, block it, or do neither to transcription of this "
             "gene.")
    L.append("")
    L.append("| Pocket | Docked | Failed | ROC AUC | EF 1% | Validated? |")
    L.append("|---|---|---|---|---|---|")
    for r in runs:
        auc = r.get("auc")
        ef1 = r.get("ef1")
        L.append(f"| {r.get('target') or '?'} | {r.get('n_docked', 0)} | "
                 f"{r.get('n_failed', 0)} | "
                 f"{f'{auc:.2f}' if auc is not None else 'not tested'} | "
                 f"{f'{ef1:.1f}x' if ef1 is not None else 'not tested'} | "
                 f"{'yes' if r.get('usable') else '**no**'} |")
    L.append("")
    L.append(f"Validated means the same scoring function separated known "
             f"actives from decoys in this pocket, at AUC "
             f"{MIN_DOCKING_AUC:.2f} or better and {MIN_EF1:.0f}x or better "
             f"enrichment in the top 1%. Below that the geometry may still "
             f"be right and the score is not evidence.")
    L.append("")

    failed = [r for r in runs if not r.get("usable")]
    if failed:
        n_f = len(failed)
        L.append(f"**{n_f} {'pocket' if n_f == 1 else 'pockets'} did not "
                 f"validate** "
                 f"({', '.join(str(r.get('target')) for r in failed[:5])}). "
                 f"Their hits are listed for the record and should not be "
                 f"ranked against hits from a pocket that did.")
        L.append("")
    for r in runs:
        for c in (r.get("caveats") or [])[:3]:
            L.append(f"- {r.get('target')}: {c}")
    if any(r.get("caveats") for r in runs):
        L.append("")

    if consensus:
        n = sum(len(v) for v in consensus.values())
        L.append(f"#### {n} {'compound' if n == 1 else 'compounds'} "
                 f"picked by both methods")
        L.append("")
        L.append("Docking reads pocket shape. The ligand screen reads "
                 "chemistry measured at this receptor. They share no input, "
                 "so agreement is worth more than a better score from "
                 "either - and a compound either method missed is not thereby "
                 "ruled out.")
        L.append("")
        L.append("| Factor | Drug | Dock score | Ligand screen quadrant |")
        L.append("|---|---|---|---|")
        for sym, hits in sorted(consensus.items()):
            for h in hits[:8]:
                sc = h.get("docking_score")
                L.append(f"| {sym} | **{h.get('drug')}** | "
                         f"{f'{sc:.1f}' if isinstance(sc, (int, float)) else EMDASH}"
                         f" | {h.get('ligand_quadrant') or EMDASH} |")
        L.append("")
        L.append("Still a binding hypothesis. Neither method addresses "
                 "direction, and neither was run in the tissue this gene "
                 "matters in.")
        L.append("")
    elif any(r.get("usable") for r in runs):
        L.append("No compound reached the top slice of both the dock and the "
                 "ligand screen. That is a real disagreement between two "
                 "independent arguments, not a missing result - treat the "
                 "single-method hits above accordingly.")
        L.append("")
    return L


def _sensitivity_section(sens: dict) -> list:
    """Drugs that behave as though they remove a factor.

    The only section here whose candidates needed no annotation of any
    kind, so it is the only one that can name a compound nobody has
    connected to the gene.
    """
    if not sens:
        return []
    from .sensitivity import MIN_ABS_R, MAX_Q, MIN_LINES

    assessed = {k: v for k, v in sens.items() if v.get("assessed")}
    if not assessed:
        reasons = {v.get("error") for v in sens.values() if v.get("error")}
        L = ["### Drugs that behave like losing the factor", ""]
        L.append("Nothing was scanned. " + (
            "; ".join(r for r in reasons if r)[:400] or "No reason recorded."))
        L.append("")
        return L

    L = ["### Drugs that behave like losing the factor", ""]
    L.append("Every other route in this report starts from something "
             "somebody wrote down - a curated ligand, a motif, a "
             "structure, a signed edge. None of them can surface a drug "
             "nobody has connected to the factor. This one reads a "
             "measured phenotype instead: if a compound hurts the same "
             "cell lines that knocking the factor out hurts, it behaves "
             "as though it removes the factor.")
    L.append("")

    rows = [(sym, h) for sym, v in sorted(assessed.items())
            for h in (v.get("hits") or [])]
    if not rows:
        n = sum(v.get("n_drugs_tested", 0) for v in assessed.values())
        L.append(f"{n} compound(s) were correlated across "
                 f"{len(assessed)} factor(s) and none cleared "
                 f"|r| >= {MIN_ABS_R} at q <= {MAX_Q}. That is the common "
                 f"outcome. It bounds what this method contributes here "
                 f"rather than saying no such drug exists: the panel is "
                 f"proliferating cancer lines, and a factor that matters "
                 f"elsewhere leaves no trace in it.")
        L.append("")
        return L

    L.append("| Factor | Drug | Partial r | Raw r | Lines | q | Behaves as |")
    L.append("|---|---|---|---|---|---|---|")
    for sym, h in sorted(rows, key=lambda t: -abs(t[1].get("r") or 0))[:20]:
        L.append(f"| {sym} | **{h.get('drug')}** | {h.get('r')} | "
                 f"{h.get('r_raw')} | {h.get('n_lines')} | "
                 f"{h.get('q_value')} | {h.get('direction')} |")
    L.append("")
    L.append(f"**Partial r** controls for how easily each cell line dies "
             f"under the whole panel. Some lines die under almost "
             f"anything and those same lines suffer most knockouts, which "
             f"puts a correlation between almost any drug and almost any "
             f"gene. The control moves a number either way - it removes a "
             f"shared component that can inflate a correlation or mask "
             f"one - so where the two columns differ a lot the raw value "
             f"was mostly about the cell lines. Estimates use at least "
             f"{MIN_LINES} shared lines.")
    L.append("")
    L.append("**This is not a target.** A correlated drug may act on the "
             "factor, on anything else the same lines need, or on a "
             "pathway the factor sits in, and the first is the least "
             "likely of the three for any one hit. Nothing here records "
             "any of them binding anything: the finding is a reason to "
             "run a binding assay, not a result from one.")
    L.append("")
    L.append("**A proliferation readout.** Both panels score growth of "
             "transformed cells in a dish, so a factor whose role is in a "
             "differentiated or post-mitotic cell looks silent here for "
             "reasons that have nothing to do with its importance.")
    L.append("")
    return L


def _depmap_section(deps: dict, res: dict) -> list:
    """What a cell loses without each factor, and what that does not mean."""
    from . import depmap as dm

    if not deps:
        return []
    n_scanned = len(res.get("tfs") or [])
    L = ["### What a cell loses without each factor (DepMap)", ""]
    L.append("The regulon column counts how many genes a factor controls, "
             "from literature curation, which measures how much anyone has "
             "published as much as what the factor does. This is the same "
             "collateral question answered by measurement: CRISPR knockout "
             "of each factor across roughly 1,100 cancer cell lines.")
    L.append("")
    L.append("It also reaches factors the curated tables cannot. A factor "
             "nobody has written a regulon for still has a knockout "
             "phenotype here.")
    L.append("")
    L.append(dm.summarise({k: _as_dep(v) for k, v in deps.items()},
                          requested=n_scanned))
    L.append("")
    L.append("| Factor | Cells need it? | Mean effect | Lines dependent | "
             "Most dependent lineage |")
    L.append("|---|---|---|---|---|")
    rows = sorted(deps.values(), key=lambda d: d.get("mean_effect") or 0.0)
    for d in rows[:25]:
        # Only when there is a dependency to localise. For a dispensable
        # factor the "most dependent" lineage is whichever one drifted
        # furthest below zero, which is noise dressed as a finding.
        lin = (d.get("top_lineages") or [{}])[0]
        lin_cell = "—"
        if d.get("n_dependent") and lin.get("lineage"):
            lin_cell = f"{lin['lineage']} ({lin['mean_effect']:+.2f})"
        L.append(f"| {d['gene']} | {d['verdict']} | "
                 f"{d['mean_effect']:+.2f} | "
                 f"{d['n_dependent']}/{d['n_lines']} "
                 f"({d['fraction_dependent']:.0%}) | {lin_cell} |")
    L.append("")
    L.append("**Mean effect** is on DepMap's own scale: 0 means knockout "
             "changed nothing, and -1 is the median of known "
             "common-essential genes. A line counts as dependent below "
             "-0.5.")
    L.append("")
    L.append("#### What this does not mean")
    L.append("")
    L.append("- **These are cancer cell lines dividing in a dish.** "
             "Essentiality here is a statement about proliferating "
             "transformed cells. A factor whose job is in a post-mitotic "
             "or differentiated cell looks dispensable for reasons that "
             "have nothing to do with its importance in your tissue.")
    L.append("- **Knockout is not inhibition.** The protein is removed "
             "entirely, including scaffolding roles no ligand would touch. "
             "Essential genes are drugged successfully at partial "
             "occupancy all the time, so a dependency is a reason for "
             "care rather than a veto.")
    L.append("- **Nothing here is about your gene.** This bounds the cost "
             "of acting on a factor. Whether acting on it moves your gene "
             "is the rest of the report.")
    L.append("- **A dispensable result is not permission.** It bounds the "
             "proliferative cost and nothing else.")
    L.append("")
    sel = [d for d in deps.values() if d.get("selective")]
    if sel:
        L.append(f"**The selective ones are the interesting ones.** "
                 f"{', '.join(d['gene'] for d in sel[:6])} matter in a "
                 f"minority of lines rather than all of them, so the "
                 f"requirement has a context. The lineage column is where "
                 f"to look for what that context is, and whether your "
                 f"tissue shares it is the question this panel cannot "
                 f"answer.")
        L.append("")
    return L


def _as_dep(d: dict):
    """A dict row back into the shape depmap.summarise reads."""
    from . import depmap as dm
    return dm.Dependency(**{k: v for k, v in d.items()
                            if k in dm.Dependency.__dataclass_fields__})


def _drug_glossary(routes: list, res: dict) -> list:
    """Every distinct drug in this block, with what it is and what it treats.

    One entry per drug rather than a column on each table. The tables are
    already wide, the same drug appears in several of them, and a sentence
    repeated five times is what made these reports long. Salts collapse
    onto one entry.
    """
    from . import indications as ind_mod

    mechanisms = {}
    for r in routes:
        name = r.get("drug")
        if name and name not in mechanisms:
            mechanisms[name] = (r.get("mechanism") or "").strip() or None

    described = ind_mod.describe_all(list(mechanisms), mechanisms=mechanisms)
    if not described:
        return []

    known = {k: v for k, v in described.items() if v.get("treats")}
    unknown = [k for k, v in described.items() if not v.get("treats")]

    L = ["### What these drugs are", ""]
    L.append("Every route in this report is off-label by construction: the "
             "drug was approved for one thing and is being considered for "
             "another. The approved use is the baseline that departure is "
             "measured from, so it is listed here once per drug rather "
             "than repeated in each table above.")
    L.append("")

    if known:
        L.append("| Drug | What it is | Approved to treat |")
        L.append("|---|---|---|")
        for name in sorted(known):
            d = known[name]
            L.append(f"| {name} | {d['what_it_is'] or '—'} | "
                     f"{d['treats']} |")
        L.append("")

    if unknown:
        L.append(f"**No description established for {len(unknown)}**: "
                 + ", ".join(sorted(unknown)[:25])
                 + ("…" if len(unknown) > 25 else "")
                 + ". That means this package could not find one offline, "
                 "not that the compound has no approved use. Several are "
                 "environmental agents and research chemicals from the CTD "
                 "tables, which have no indication because they are not "
                 "drugs.")
        L.append("")

    L.append(ind_mod.summarise(described))
    L.append("")
    return L


def _screen_scope_note(ran: list, skipped: list) -> list:
    """Why the screen covers the receptors it does, and no others.

    Written because the section looks nearly identical across unrelated
    genes, which invites the reading that these drugs are being selected
    for the gene. They are not. They are the receptors with enough
    measured ligands to train a model on.
    """
    from .agonism import PCHEMBL_CUTOFF
    from .screening import MIN_ACTIVES, MIN_BINDING_AUC

    out = ["## Virtual screening", ""]
    out.append(f"A screen needs a model, and a model needs data. This one "
               f"trains on measured binding at the receptor, so it runs "
               f"only where ChEMBL holds at least {MIN_ACTIVES} ligands at "
               f"pChEMBL {PCHEMBL_CUTOFF:.0f} or better, and only when the "
               f"trained model separates known actives from decoys at CV "
               f"AUC {MIN_BINDING_AUC} or above.")
    out.append("")
    if ran:
        out.append(f"**Screened:** {', '.join(sorted(ran))}.")
    if skipped:
        out.append(f"**Not screened:** {', '.join(sorted(skipped))} "
                   f"{chr(8212)} each one's reason is stated under its own "
                   f"heading below. A skip is always a limit of the "
                   f"available data, either too few measured ligands or a "
                   f"model that could not separate them from decoys, and "
                   f"never a judgement about the receptor.")
    out.append("")
    out.append("**This is why the same receptor families appear here "
               "whatever gene you run.** The filter selects for how much "
               "medicinal chemistry a receptor has attracted, not for "
               "anything about your promoter. Nuclear receptors, and the "
               "retinoid and steroid receptors above all, have ligand "
               "series running to thousands of compounds; most "
               "transcription factors have none at all. So a receptor "
               "reaching this section means it is well studied, and the "
               "drugs listed under it are that receptor's pharmacology, "
               "not a result about your gene.")
    out.append("")
    out.append("Read these rows as *candidate ligands for a factor that has "
               "a site in your promoter*. Whether moving that factor moves "
               "your gene is a separate question, and the activity and "
               "evidence sections above speak to it. A long list here is "
               "evidence about the receptor's chemistry; it is not evidence "
               "that the gene will respond.")
    out.append("")
    return out


def _regulon_column_note(any_regulon: bool, any_understudied: bool,
                         tfs: list) -> list:
    """What the regulon column means, and the two ways it misleads."""
    if not any_regulon:
        return ["*No regulon counts are available. Pass `--regulon-file` "
                "with a TRRUST-style table to fill the column; without one "
                "there is no basis for saying how much else a factor "
                "controls.*", ""]

    out = ["**Other genes it controls** is the collateral question. A drug "
           "acting through a factor moves everything that factor regulates, "
           "so this column says how wide the effect will be. Read it against "
           "the score: a strong route through a factor that controls two "
           "thousand genes is a worse starting point than a weaker route "
           "through one that controls thirty."]

    exact = [t for t in tfs if t.get("n_regulon") is not None
             and "tier" not in (t.get("regulon_provenance") or "")]
    tiered = [t for t in tfs if "tier" in (t.get("regulon_provenance") or "")]
    uncounted = [t for t in tfs if t.get("n_regulon") is None]

    if tiered:
        out.append("")
        out.append(f"A count written `~N` is a bundled estimate for a factor "
                   f"the regulon file does not list ({len(tiered)} here), not "
                   f"a count of curated edges. Treat it as an order of "
                   f"magnitude.")
    if uncounted:
        out.append("")
        out.append(f"`not counted` means the regulon table has no entry and "
                   f"there is no bundled estimate either ({len(uncounted)} "
                   f"here). It is not a small number. Most of the roughly "
                   f"1,600 human transcription factors have no curated "
                   f"regulon at all.")
    if any_understudied:
        out.append("")
        out.append("A ⚠ marks a count below the credible floor. These "
                   "are the trap on this column: a factor with eight curated "
                   "targets is far more likely to be understudied than "
                   "selective, and reading a small number as specificity "
                   "selects for obscurity. The same correction applies here "
                   "that `drug_scoring` applies to protein selectivity.")
    out.append("")
    out.append("Counts come from literature curation, so they measure what "
               "has been published as much as what the factor does. SP1 has "
               "thousands of recorded targets partly because it regulates "
               "many genes and partly because thousands of people assayed "
               "it.")
    out.append("")
    return out


def _genetics_detail(phenotypes: list, local: dict | None) -> str:
    """What the genetics step actually found, naming which source found it.

    Open Targets and a local GWAS catalogue are not the same evidence: the
    first fine-maps a locus to a gene with a trained model, the second
    reports the gene nearest the signal. Reporting one count would hide
    which of the two the row rests on.
    """
    parts = []
    if phenotypes:
        parts.append(f"{len(phenotypes)} associated phenotype(s) "
                     f"(Open Targets)")
    if local:
        n = _n_local_gwas(local)
        if n:
            parts.append(f"{n} genome-wide significant GWAS Catalog "
                         f"association(s), mapped or nearest gene, which is "
                         f"weaker than fine-mapped support")
        else:
            parts.append("a local GWAS Catalog was read and holds no "
                         "significant association for this gene")
    if not parts:
        return "no genetic association returned for this gene"
    return "; ".join(parts)


def _n_local_gwas(local: dict | None) -> int:
    """Significant associations only. `assessed` says the file was read."""
    if not local:
        return 0
    return int(local.get("n_significant") or 0)


def _sig_q(value) -> float:
    """A q-value for comparison, with None - not 0.0 - meaning missing.

    `value or 1` reads a q of exactly 0.0, the MOST significant result
    possible, as the least. `tfactivity._q` exists for this and this line
    reimplemented the broken version. Latent rather than live today only
    because `tfactivity._bh` floors its output at Q_FLOOR, an invariant
    held in another module.
    """
    return 1.0 if value is None else float(value)


def _route_cell(row: dict) -> str:
    """FDA route, compressed. "local only" is the one that rules a drug out."""
    r = row.get("route")
    if r == "local only":
        return "**local only**"
    if r in ("systemic", "unclear"):
        return r
    return "—"


def _bbb_footnote(counts: dict) -> list:
    """What the BBB column means, stated where it is read.

    `counts` maps the cells actually printed ("yes", "yes?", "no", "no?",
    "—") to how many rows carried them, so the note describes this table
    rather than warning about states that did not occur in it.
    """
    if not counts:
        return []
    total = sum(counts.values())
    n_unknown = counts.get("—", 0)
    n_yes = counts.get("yes", 0)
    n_borderline = counts.get("yes?", 0) + counts.get("no?", 0)

    out = ["**BBB?** is CNS MPO (Wager et al.) against the same scaled "
           "threshold the delivery section uses, so a drug cannot read "
           "\"yes\" here and \"unlikely\" there. It answers whether the "
           "molecule can reach brain at all: for a CNS target that decides "
           "whether a candidate is usable, and for a peripheral target it "
           "decides whether brain exposure is an unwanted side effect rather "
           "than the goal."]
    if n_yes or n_borderline:
        out.append("")
        out.append("A **yes** predicts *passive diffusion only*, and is the "
                   "weaker of the two answers. Active efflux by "
                   "P-glycoprotein at the endothelium pumps an otherwise "
                   "suitable molecule straight back out; known substrates "
                   "are checked here and a hit overrides a pass outright.")
    if n_borderline:
        out.append("")
        out.append(f"**A question mark means undecided** — the score landed "
                   f"within noise of the threshold, and {n_borderline} of "
                   f"{total} row(s) here did. This is not pedantry. "
                   f"**Atenolol clears a complete five-parameter MPO by 0.04 "
                   f"and is a textbook poor penetrant**, so a margin that "
                   f"thin carries no information. Read a question-marked "
                   f"cell as \"look this one up\", never as the answer it "
                   f"leans toward.")
    if n_unknown:
        out.append("")
        out.append(f"**— is not a no.** {n_unknown} of {total} row(s) "
                   f"returned no physicochemical properties — common for "
                   f"DGIdb- and Open Targets-sourced rows, which carry a "
                   f"name and a mechanism but no structure. Nothing was "
                   f"measured and nothing is being claimed; look the "
                   f"compound up before excluding it on this column.")
    out.append("")
    out.append("A **no** is the more reliable direction. A molecule failing "
               "on size, polarity or lipophilicity does not cross by "
               "diffusion, and no formulation fixes that short of bypassing "
               "the barrier — intrathecal delivery, convection-enhanced "
               "delivery, or focused ultrasound.")
    out.append("")
    return out


def markdown_report(result: dict, top_n: int = 25) -> str:
    q = result["query"]
    goal = q.get("goal", "up")
    verb = "Upregulation" if goal == "up" else "Downregulation"
    L: list[str] = []
    L.append(f"# {verb} strategy for {q['gene']} — promoter / TF / drug analysis")
    L.append("")
    L.append(f"Objective: **{'increase' if goal == 'up' else 'decrease'}** "
             f"transcription of {q['gene']}.")
    L.append("")
    if q["demo"]:
        L.append("> **DEMO MODE.** Synthetic sequence, illustrative matrices and a "
                 "hand-entered drug table. Not usable for real conclusions.")
        L.append("")

    # -------------------------------------------- can this run be read at all?
    # Printed before any ranking. A reader who sees only the tables below
    # cannot tell a measured negative from a lookup that broke, and every
    # table in this report is structurally identical in both cases.
    ev = result.get("evidence_verdict") or {}
    srcs = (result.get("evidence_sources") or {}).get("sources") or {}
    if ev or srcs:
        L.append("## Evidence available in this run")
        L.append("")
        badge = {"not interpretable": "⛔ NOT INTERPRETABLE",
                 "degraded": "⚠ DEGRADED",
                 "interpretable": "✓ INTERPRETABLE"}.get(ev.get("level"), "")
        if badge:
            L.append(f"**{badge}** — {ev.get('headline', '')}")
            L.append("")
        for reason in ev.get("reasons") or []:
            L.append(f"- {reason}")
        if ev.get("reasons"):
            L.append("")
        if srcs:
            L.append("| Source | State | What it returned |")
            L.append("|---|---|---|")
            for name, s in srcs.items():
                state = {"ok": "ok", "FAILED": "**FAILED**",
                         "not requested": "not requested"}.get(s["state"],
                                                               s["state"])
                L.append(f"| {name} | {state} | {s['detail']} |")
            L.append("")
        failed = [s for s in srcs.values() if s["state"] == "FAILED"]
        if failed:
            L.append("**A failed source is not a negative result.** These were "
                     "requested and did not deliver, so nothing below rests on "
                     "them:")
            L.append("")
            for s in failed:
                L.append(f"- **{s['name']}** — {s.get('impact', '')}")
            L.append("")
        missing = [s for s in srcs.values() if s["state"] == "not requested"]
        if missing:
            L.append("Not requested (enable to strengthen the run): "
                     + ", ".join(f"{s['name']} (`{s['enable_with']}`)"
                                 for s in missing if s.get("enable_with"))
                     + ".")
            L.append("")

    direct_hits = result.get("direct_therapies") or []
    if direct_hits:
        L.append(f"## Drugs targeting {q['gene']} directly")
        L.append("")
        L.append("These act on the gene, transcript or protein itself, skipping "
                 "the transcription-factor route entirely. Where one of these "
                 "exists it is almost always the better starting point: "
                 "sequence-directed agents affect your gene and essentially "
                 "nothing else.")
        L.append("")
        L.append("| Drug | Modality | Specificity | Raises protein? | Mechanism | "
                 "Approved | Source |")
        L.append("|---|---|---|---|---|---|---|")
        for t in direct_hits:
            flag = " ⚠withdrawn" if t["withdrawn"] else ""
            raises = "yes" if t["raises_protein"] else "activity only"
            L.append(f"| {t['drug']}{flag} | {t['modality']} | {t['specificity']} | "
                     f"{raises} | {t['mechanism']} | {t['approval_year'] or '—'} | "
                     f"{t['source']} |")
        L.append("")
        if any(not t["raises_protein"] for t in direct_hits):
            L.append("*\"activity only\" means the drug raises the function of protein "
                     "already present without increasing its abundance. If your readout "
                     "is a western blot rather than an activity assay, that distinction "
                     "decides whether the experiment appears to work.*")
            L.append("")
    elif result["query"].get("include_direct", True):
        L.append(f"## Drugs targeting {q['gene']} directly")
        L.append("")
        L.append("None found in the curated table or in ChEMBL/DGIdb for this "
                 "direction. The curated list is hand-maintained and goes out of "
                 "date; check FDA/EMA and clinicaltrials.gov before concluding "
                 "nothing exists.")
        L.append("")

    # ------------------------------------- the headline answer, up front
    # This is the question people actually bring to the tool. It used to be
    # answerable only by assembling eight tables yourself.
    if result["results"]:
        first = result["results"][0]
        burdens = {}
        for g in (first.get("drug_grades") or []):
            burdens[g["drug_name"].lower()] = sys_mod.clinical_burden(
                g["drug_name"], notes=" ".join(g.get("notes") or []),
                target=g.get("tf_symbol")).to_dict()
        preds = pred_mod.build(
            first, goal=goal,
            has_chip=bool(q.get("chip_bed") or q.get("chip_api")),
            burdens=burdens)
        L.extend(pred_mod.render(preds, q["gene"], goal,
                                 verdict=result.get("evidence_verdict")))

    ov = result.get("tss_window_overlap") or []
    if ov:
        L.append("> **⚠ These promoter windows overlap each other.** "
                 + " ".join(o["note"] for o in ov[:2])
                 + " Counts and route lists below are therefore not "
                 "independent across the windows.")
        L.append("")

    for res in result["results"]:
        p = res["promoter"]
        L.append(f"## Promoter of {p['gene_symbol']} ({p['transcript_id']})")
        L.append("")
        L.append(f"- Transcript chosen: **{p['transcript_selection']}**")
        L.append(f"- Assembly: {p['assembly']}  |  Strand: {'+' if p['strand'] == 1 else '-'}")
        L.append(f"- TSS: chr{p['chrom']}:{p['tss']}")
        L.append(f"- Window: chr{p['chrom']}:{p['start']}-{p['end']} "
                 f"(-{p['upstream']} / +{p['downstream']}, {p['sequence_length']} bp)")
        L.append(f"- Motifs scanned: {q['n_motifs']} at p ≤ {q['pvalue']:g} per position")
        L.append(f"- Non-redundant sites: {res['n_sites']} across {res['n_tfs']} TFs")
        L.append("")
        q2 = res.get("scan_quality") or {}
        if q2:
            if q2["verdict"] == "no signal":
                L.append(f"> **⚠ This scan found no signal.** {q2['advice']}")
                L.append(">")
                L.append("> Everything below is ranked from a list of chance "
                         "matches. Fix the scan before reading any of it.")
            elif q2["verdict"] == "weak":
                L.append(f"> **⚠ Weak scan.** {q2['advice']}")
            else:
                L.append(f"*Scan quality: {q2['advice']}*")
            L.append("")
            if q2.get("fdp_note"):
                # The enrichment ratio understates this badly: 1.2x reads as
                # mild signal but means ~5 of every 6 rows are chance.
                L.append(f"*{q2['fdp_note']}*")
                L.append("")

        nb = res.get("neighbours") or {}
        if nb.get("checked") and nb.get("genes"):
            L.append("#### What else is in this window")
            L.append("")
            L.append("Every other specificity check here asks whether a drug "
                     "or factor affects other genes. This asks the prior "
                     "question: is the sequence that was scanned uniquely "
                     "this gene's?")
            L.append("")
            L.append(nb["summary"])
            L.append("")
            L.append("| Gene | Biotype | Location | Strand | Relation |")
            L.append("|---|---|---|---|---|")
            for g in nb["genes"][:12]:
                L.append(f"| {g['gene']} | {g['biotype']} | "
                         f"chr{g['chrom']}:{g['start']}-{g['end']} | "
                         f"{'+' if g['strand'] == 1 else '-'} | "
                         f"{g['relation']} |")
            L.append("")
            for g in nb["genes"][:4]:
                if g["relation"] != "nearby":
                    L.append(f"- **{g['gene']}** — {g['note']}")
            if nb.get("n_antisense"):
                L.append("")
                L.append("*An antisense overlap matters more than it looks. "
                         "The direction logic in this report assumes a factor "
                         "acts on YOUR promoter; where the window lies inside "
                         "an antisense transcript, a site may drive that "
                         "instead, and antisense RNA can stabilise or repress "
                         "the sense gene — so the sign does not carry.*")
            L.append("")
            L.append("*This checks annotation, not sequence uniqueness. BLAT "
                     "the window if you need to know it appears nowhere else "
                     "in the genome.*")
            L.append("")

        nrs = [t for t in res["tfs"] if t.get("nuclear_receptor")]
        if nrs:
            L.append("### Ligand-activated (nuclear) receptors with sites here")
            L.append("")
            L.append("These are the directly druggable ones.")
            L.append("")
            L.append("| TF | Role | Sites | Best rel. score | Best site (vs TSS) | ChIP peaks | Score |")
            L.append("|---|---|---|---|---|---|---|")
            for t in nrs[:top_n]:
                L.append(f"| {t['symbol']} | {t['role']} | {t['n_sites']} | "
                         f"{t['best_rel_score']} | {t['best_site_sequence']} "
                         f"({_fmt_offset(t['best_site_tss_offset'])}) | "
                         f"{t['chip_peaks']} | {t['score']} |")
            L.append("")

        at = res.get("attrition") or {}
        if at.get("n_tfs_with_sites"):
            L.append("#### Why the drug tables below look the way they do")
            L.append("")
            L.append(f"**{at['n_tfs_with_sites']} factor(s) had predicted "
                     f"sites here. {at['n_with_drug']} of them have an "
                     f"approved drug"
                     + (f", and all {at['n_drugged_nuclear_receptor']} are "
                        f"nuclear receptors.**" if at.get("all_drugged_are_nr")
                        else f" ({at['n_drugged_nuclear_receptor']} of those "
                             f"are nuclear receptors).**"))
            L.append("")
            L.append("Roughly 48 of the ~1,600 human transcription factors "
                     "have a ligand-binding pocket, and the steroid, "
                     "retinoid, thyroid and vitamin-D families are most of "
                     "them. So the drug tables are not a ranking that "
                     "hormone drugs won — they are what is left after every "
                     "factor without a ligand has been dropped.")
            L.append("")
            lost = at.get("top25_without_drug") or []
            if lost:
                L.append(f"**{at.get('n_top25_without_drug')} of the top 25 "
                         f"factors by score have no drug and appear in no "
                         f"drug table:** " + ", ".join(lost[:12])
                         + (" …" if len(lost) > 12 else "") + ".")
                L.append("")
                L.append("That is where the gene-specific signal is. A "
                         "half-site like AGGTCA occurs every few kilobases "
                         "by chance, so nuclear receptors surface for almost "
                         "any promoter; these factors surfaced for *this* "
                         "one. The routes that can reach them are the "
                         "upstream-network, stabilisation and "
                         "protein-level sections, not the ligand tables.")
                L.append("")

        L.append("### Top-ranked transcription factors")
        L.append("")
        has_dep = any(t.get("dependency") for t in (res.get("tfs") or []))
        dep_col = " Cells need it? |" if has_dep else ""
        dep_sep = "---|" if has_dep else ""
        L.append("| Rank | TF | Role | NR | Other genes it controls |"
                 + dep_col +
                 " Sites | Best rel. score | Offset | ChIP | "
                 "Approved drugs | Score |")
        L.append("|---|---|---|---|---|" + dep_sep +
                 "---|---|---|---|---|---|")
        any_regulon = False
        any_understudied = False
        for i, t in enumerate(res["tfs"][:top_n], 1):
            n = t.get("n_regulon")
            if n is None:
                reg = "not counted"
            else:
                any_regulon = True
                reg = f"~{n}" if "tier" in (t.get("regulon_provenance") or "") \
                    else str(n)
                if t.get("regulon_understudied"):
                    any_understudied = True
                    reg += " \u26a0"
            dep_cell = ""
            if has_dep:
                dep_cell = f" {_dep_label(t)} |"
            L.append(f"| {i} | {t['symbol']} | {t['role']} | "
                     f"{'yes' if t['nuclear_receptor'] else ''} | {reg} |"
                     + dep_cell +
                     f" {t['n_sites']} | "
                     f"{t['best_rel_score']} | {_fmt_offset(t['best_site_tss_offset'])} | "
                     f"{t['chip_peaks']} | {len(t['drugs'])} | {t['score']} |")
        L.append("")
        L.extend(_regulon_column_note(any_regulon, any_understudied,
                                      res.get("tfs") or []))

        L.extend(_depmap_section(res.get("dependencies") or {}, res))

        L.extend(_sensitivity_section(res.get("sensitivity") or {}))

        L.extend(_docking_section(res.get("docking") or [],
                                  res.get("docking_consensus") or {}))

        act = res.get("tf_activity") or {}
        if act.get("n_testable"):
            L.append("### Is the factor actually active? (inferred from its "
                     "targets)")
            L.append("")
            L.append("Everything above this point is a *sequence* argument: a "
                     "PWM says a factor could bind. This is a different "
                     "question asked of different data - do the genes that "
                     "factor is known to regulate actually move, in the "
                     "direction its mode of regulation predicts? That is the "
                     "VIPER idea (Alvarez et al. 2016), and it is the closest "
                     "this package gets to measuring a factor rather than "
                     "predicting one.")
            L.append("")
            L.append("**Why it belongs next to the table above.** At p=1e-4 "
                     "over a 6 kb window most motif calls are noise. The two "
                     "lines of evidence are independent, so a factor with "
                     "both a site and a coherently perturbed regulon is "
                     "supported in a way that neither alone can match - and "
                     "a strong site with a flat regulon is the expected "
                     "failure mode of the scan, not a surprise.")
            L.append("")
            L.append(act.get("summary", ""))
            L.append("")

            scanned = act.get("scanned") or []
            if scanned:
                L.append(f"#### Factors with a site here that could be "
                         f"tested")
                L.append("")
                L.append(f"{len(scanned)} of the "
                         f"{len(res.get('tfs') or [])} factors with a motif "
                         f"hit in this window have enough curated targets to "
                         f"test. The rest are not absent from the result - "
                         f"they were never testable, which is a gap in the "
                         f"regulon table and not evidence about the factor.")
                L.append("")
                L.append("| TF | Activity | NES | q | Targets measured | "
                         "Coverage | Signed | Confounded with |")
                L.append("|---|---|---|---|---|---|---|---|")
                for a in scanned[:top_n]:
                    conf = a.get("confounded_with")
                    j = a.get("confounded_jaccard")
                    cstr = (f"{conf} (J={j})" if conf and (j or 0) >= 0.3
                            else "—")
                    L.append(
                        f"| {a['tf']} | {a['direction']} | "
                        f"{a['nes'] if a['nes'] is not None else '—'} | "
                        f"{a['q_value'] if a['q_value'] is not None else '—'} | "
                        f"{a['n_measured']}/{a['n_regulon']} | "
                        f"{a['coverage']:.0%} | "
                        f"{a['n_signed']} | {cstr} |")
                L.append("")
            else:
                L.append("*No factor with a site in this promoter had enough "
                         "curated targets to test. With TRRUST that is the "
                         "common case rather than a failure: most of the "
                         "~1,600 human transcription factors have too few "
                         "curated targets to support any activity estimate.*")
                L.append("")

            top = [a for a in (act.get("top") or [])
                   if a["tf"] not in {s["tf"] for s in scanned}]
            if top:
                L.append(f"#### The strongest calls overall, site or not")
                L.append("")
                L.append("These have no motif hit in your promoter window, so "
                         "they are not routes to your gene. They are here "
                         "because they say what this tissue's regulatory state "
                         "looks like, which is the context any route has to "
                         "work inside.")
                L.append("")
                L.append("| TF | Activity | NES | q | Targets measured |")
                L.append("|---|---|---|---|---|")
                for a in top[:10]:
                    L.append(
                        f"| {a['tf']} | {a['direction']} | "
                        f"{a['nes'] if a['nes'] is not None else '—'} | "
                        f"{a['q_value'] if a['q_value'] is not None else '—'} | "
                        f"{a['n_measured']}/{a['n_regulon']} |")
                L.append("")

            L.append("#### What this cannot tell you")
            L.append("")
            L.append(f"- **It is a reference tissue average, not your "
                     f"sample.** The signature is {act['tissue']} from a "
                     f"healthy-tissue expression table, so it says nothing "
                     f"about a patient, a cell line, or a disease state.")
            L.append("- **Activity is not attributable to one factor.** "
                     "Regulons overlap, and co-regulated targets cannot "
                     "distinguish a factor from its neighbour. Where the "
                     "overlap is large the table names the other factor; the "
                     "signal is real and which of the two produced it is "
                     "not decided.")
            L.append("- **Nuclear receptors are the hardest case here.** A "
                     "ligand-activated receptor can be transcribed "
                     "everywhere while being active only where its ligand "
                     "reaches, so a quiet regulon does not mean the receptor "
                     "is absent - it means nothing is switching it on in "
                     "healthy tissue. For those factors, read a flat result "
                     "as uninformative rather than negative.")
            L.append("- **Direction needs a signed edge, and half of them are "
                     "not.** TRRUST records no mode for 46% of its human "
                     "edges, so those support a two-tailed test (did the "
                     "regulon move) and no directional one. A row with a "
                     "dash under NES is that case.")
            src = (act.get("regulon_source") or "").lower()
            if "collectri" not in src:
                L.append(f"- **Coverage is the binding constraint, and it is "
                         f"fixable.** {act['n_regulons']} factors had a "
                         f"regulon here, against roughly 1,600 human "
                         f"transcription factors, so most were never "
                         f"testable. CollecTRI covers about three times as "
                         f"many and signs every edge, which would widen this "
                         f"section and let more of the dash rows resolve to "
                         f"a direction. Fetch it with "
                         f"`import decoupler as dc; "
                         f"dc.get_collectri(organism='human')`, save it to "
                         f"`data/collectri.tsv`, and it is picked up "
                         f"automatically.")
            if not act.get("permutations"):
                L.append("- **The p-values are analytic and "
                         "anticonservative.** They assume targets are "
                         "independent, and co-expression is exactly the "
                         "structure being measured. Re-run with "
                         "`--activity-permutations 2000` before relying on "
                         "one.")
            L.append("")

        ad = res.get("activity_drugs") or {}
        if ad:
            n_pairs = sum(len(v) for v in ad.values())
            L.append("### Drugs that move a factor, found without any "
                     "record that they bind it")
            L.append("")
            L.append("Every drug table above came from asking a database "
                     "which compounds were assayed at the target. That "
                     "question cannot return anything for a factor nobody "
                     "has done chemistry on, which is most of them, and it "
                     "is why those tables fill up with nuclear receptors.")
            L.append("")
            L.append("This asks the LINCS data instead. Each compound in "
                     "LINCS has its own transcriptional signature, so the "
                     "same regulon statistic used above can be run on it: "
                     "does this compound move the factor's known targets? "
                     "No annotated binding is needed, so a factor with no "
                     "pharmacology at all can still return candidates.")
            L.append("")
            L.append(f"**{n_pairs} compound-factor pair(s) across "
                     f"{len(ad)} factor(s).**")
            L.append("")
            L.append("| Factor | Drug | Corpus | Effect on factor | "
                     "Median z | Profiles agreeing | Significant |")
            L.append("|---|---|---|---|---|---|---|")
            for tf, hits in ad.items():
                for h in hits[:10]:
                    L.append(
                        f"| {tf} | {h['drug']} | "
                        f"{h.get('corpus', 'lincs')} | {h['direction']} | "
                        f"{h['median_z']:+.2f} | "
                        f"{h['n_agreeing']}/{h['n_signatures']} "
                        f"({h['consistency']:.0%}) | "
                        f"{h['n_significant']} at q<0.05 |")
            L.append("")
            L.extend(_corpus_note(res))
            L.append("**Repetition is the evidence here, not any single "
                     "profile.** The scan tests one regulon against roughly "
                     "33,000 signatures, so nominal significance is reached "
                     "about 1,650 times by chance. Significance is counted "
                     "after Benjamini-Hochberg across those signatures, and "
                     "a compound has to clear three further bars: profiled "
                     "at least "
                     f"{tfa_mod.MIN_SIGNATURES} times, moving the regulon "
                     f"the same way in at least "
                     f"{tfa_mod.MIN_CONSISTENCY:.0%} of them, with a median "
                     f"effect of at least z {tfa_mod.MIN_MEDIAN_Z}. A "
                     "compound profiled once is never reported.")
            L.append("")
            L.append("#### What a row here does not mean")
            L.append("")
            L.append("- **It is not a binding claim.** What was measured is "
                     "that the factor's targets moved when the compound was "
                     "given. The effect can run through anything upstream "
                     "of the factor, and the compound may never touch it.")
            L.append("- **The attribution to this factor is the weak "
                     "step.** Regulons overlap, so a compound that moves a "
                     "shared set of targets will score against every factor "
                     "that claims them. The activity section above names "
                     "the overlapping factor for each one.")
            L.append("- **Proliferation is the confound to rule out "
                     "first.** A compound that arrests the cell cycle moves "
                     "every regulon enriched for cell-cycle genes, so a "
                     "growth inhibitor can look like a specific effect on "
                     "any such factor. If the agreeing compounds share an "
                     "antiproliferative mechanism, suspect that before "
                     "crediting the factor.")
            L.append("- **LINCS cell lines are not your tissue.** These are "
                     "mostly cancer lines at micromolar doses for six to "
                     "twenty-four hours, which is a different regime from a "
                     "therapeutic exposure.")
            L.append("")
            L.append("These rows enter the route tables as **tier 3**, "
                     "alongside curated mechanisms rather than with the "
                     "tier-2 measured ones. The observation is real; the "
                     "attribution to the factor is an inference, and tier 3 "
                     "is what that combination is worth.")
            L.append("")

        ar = res.get("all_routes") or {}
        if ar.get("routes"):
            L.append("### Every route that acts on this gene")
            L.append("")
            L.append(f"{ar['n_routes']} routes across "
                     f"{len(ar['by_category'])} mechanism categories, pooled from "
                     "every method in this run and sorted by **evidence tier "
                     "first**, score second.")
            L.append("")
            L.append("A tier-1 route with a mediocre score beats a tier-5 route "
                     "with a great one. These are not the same kind of claim and "
                     "averaging them would destroy the information you need.")
            L.append("")
            L.append("**Why established binders sit at tier 2, just under the "
                     "direct agents.** A drug already measured to bind one of "
                     "your factors has half its chain settled: the "
                     "drug-to-receptor step has assay data behind it, so the "
                     "pharmacology is de-risked and only the receptor-to-gene "
                     "step is a hypothesis. A curated tier-3 ligand is "
                     "annotation on *both* halves. But they stay below tier 1, "
                     "because a direct agent acts on your gene itself and "
                     "skips the receptor-to-gene inference entirely — a "
                     "stronger claim however well-characterised the binder "
                     "is. The cost of the indirect route is in the row: "
                     "everything else that factor regulates moves with your "
                     "gene.")
            L.append("")
            for tier_label, n in ar["by_tier"].items():
                L.append(f"- {n} × {tier_label}")
            L.append("")
            L.append("| Tier | Route | Drug | Acts on | Steps | Mechanism | "
                     "Str | Spec | Eff | Score |")
            L.append("|---|---|---|---|---|---|---|---|---|---|")
            for r in ar["routes"][:max(top_n, 30)]:
                a = r["axes"]
                L.append(f"| {r['tier']} | {r['category']} | {r['drug']} | "
                         f"{r['acts_on']} | {r['directness']} | "
                         f"{r['mechanism'][:45]} | {a['strength']} | "
                         f"{a['specificity']} | {a['efficiency']} | {r['score']} |")
            L.append("")
            L.append("Best in each category:")
            L.append("")
            for cat, r in ar["best_per_category"].items():
                L.append(f"- **{cat}** — {r['drug']} (tier {r['tier']}, "
                         f"score {r['score']})")
            L.append("")
            if ar.get("categories_absent"):
                L.append("*No route found in: "
                         + ", ".join(ar["categories_absent"])
                         + ". An absent category is usually a real gap in the "
                         "pharmacology, but check you enabled the relevant flag "
                         "(--lincs, --screen, --upstream-network).*")
                L.append("")
            L.append("*Steps = how far the drug acts from your gene: 0 on the "
                     "gene or transcript, 1 on a TF that binds the promoter, 2+ "
                     "upstream of that TF. Each step costs specificity.*")
            L.append("")

            guide = cat_mod.CATEGORY_GUIDE
            present = [c for c in ar["by_category"] if c in guide]
            if present:
                L.append("#### What these categories mean")
                L.append("")
                L.append("Several of these sound alike and behave very "
                         "differently. The distinction that most often decides "
                         "an experiment is whether a route changes protein "
                         "**abundance** or protein **activity** — a western "
                         "blot only sees the first.")
                L.append("")

                # Which drugs are actually in each category, and what assay
                # would see the effect. Printing only "(21 found)" left the
                # reader unable to connect a category to their own candidates.
                by_cat: dict = {}
                for r in ar["routes"]:
                    by_cat.setdefault(r["category"], []).append(r["drug"])

                L.append("| Category | What it changes | Your drugs in it | "
                         "Assay that sees it |")
                L.append("|---|---|---|---|")
                for c in present:
                    g = guide[c]
                    drugs = by_cat.get(c) or []
                    seen, uniq = set(), []
                    for d in drugs:
                        if d.lower() not in seen:
                            seen.add(d.lower())
                            uniq.append(d)
                    shown = ", ".join(uniq[:4])
                    if len(uniq) > 4:
                        shown += f" … +{len(uniq) - 4}"
                    L.append(f"| {c} | **{g.get('changes', '?')}** | "
                             f"{shown or '—'} | {g.get('assay', '—')} |")
                L.append("")
                if "direct: gene product" in present:
                    L.append("**Only one of these changes activity rather "
                             "than abundance: *direct: gene product*.** "
                             "Everything else works through transcription, so "
                             "it moves how much of your gene's product exists "
                             "and an RT-qPCR will see it. A ligand or "
                             "inhibitor on the product itself changes what "
                             "the protein already present does while its "
                             "level barely moves — which is why ivacaftor "
                             "reads as a failure on a western blot and works "
                             "in patients.")
                else:
                    L.append("*Every category present here works through "
                             "transcription, so all of them move your gene's "
                             "product **abundance** and an RT-qPCR is the "
                             "right first readout. The one exception in the "
                             "taxonomy — a drug on the gene product itself, "
                             "which changes activity only — did not appear "
                             "in this run.*")
                L.append("")

                for c in present:
                    g = guide[c]
                    drugs = by_cat.get(c) or []
                    seen, uniq = set(), []
                    for d in drugs:
                        if d.lower() not in seen:
                            seen.add(d.lower())
                            uniq.append(d)
                    L.append(f"**{c}** ({ar['by_category'][c]} found)")
                    L.append("")
                    if uniq:
                        L.append(f"- *Yours:* {', '.join(uniq[:8])}"
                                 + (f" … +{len(uniq) - 8} more"
                                    if len(uniq) > 8 else ""))
                    L.append(f"- {g['what']}")
                    L.append(f"- *Why:* {g['why']}")
                    L.append(f"- *Cost:* {g['cost']}")
                    L.append(f"- *Differs:* {g['differs']}")
                    L.append("")
                absent = [c for c in (ar.get("categories_absent") or [])
                          if c in guide]
                if absent:
                    L.append("Not found here, and what you would be missing:")
                    L.append("")
                    for c in absent:
                        L.append(f"- **{c}** — {guide[c]['what']}")
                    L.append("")

        ar_g = res.get("all_routes") or {}
        if ar_g.get("routes"):
            L.extend(_drug_glossary(ar_g["routes"], res))

        if res.get("drug_grades"):
            L.append("### Graded candidates")
            L.append("")
            L.append("Combined = geometric mean of TF confidence and drug quality. "
                     "A weak TF prediction caps the grade no matter how good the drug.")
            L.append("")
            L.append("| Grade | Combined | Drug | BBB? | TF | TF conf. | "
                     "Drug qual. | Potency | Protein select. | Gene select. | "
                     "Evidence | Flags |")
            L.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
            bbb_seen: dict = {}
            for g in res["drug_grades"][:top_n]:
                c = g["components"]
                flags = "; ".join(g["notes"])[:70] or "—"
                v = dlv_mod.bbb_verdict(
                    g["drug_name"],
                    (g.get("evidence") or {}).get("properties"))
                bbb_seen[v["cell"]] = bbb_seen.get(v["cell"], 0) + 1
                L.append(f"| **{g['grade']}** | {g['combined']} | {g['drug_name']} | "
                         f"{v['cell']} | "
                         f"{g['tf_symbol']} | {g['tf_confidence']} | {g['drug_quality']} | "
                         f"{c['potency']} | {c['selectivity']} | "
                         f"{c.get('regulon', '—')} | {c.get('evidence_depth', '—')} | "
                         f"{flags} |")
            L.append("")
            L.extend(_bbb_footnote(bbb_seen))

        if res.get("cell_cycle_summary") or res.get("cell_cycle_blocked"):
            L.append("### Cell-cycle compatibility")
            L.append("")
            L.append(f"Tissue proliferation state: **{q.get('proliferation')}**")
            L.append("")
            if res.get("cell_cycle_summary"):
                L.append(res["cell_cycle_summary"])
                L.append("")
            if res.get("cell_cycle_blocked"):
                L.append("**Removed as incompatible:** "
                         + ", ".join(res["cell_cycle_blocked"]))
                L.append("")
                L.append("These are not ranked lower — they are removed. In "
                         "post-mitotic tissue, E2F/DREAM target promoters are "
                         "actively repressed rather than merely quiet, so a route "
                         "through a G1/S factor does not work at any dose.")
                L.append("")
            restricted = [(s, a) for s, a in (res.get("cell_cycle") or {}).items()
                          if a["duty"] < 1.0 and a["compatible"]]
            if restricted:
                L.append("| TF | Phase | Duty cycle | Effect multiplier |")
                L.append("|---|---|---|---|")
                for sym, a in restricted[:top_n]:
                    L.append(f"| {sym} | {a['phase']} | {a['duty']:.0%} | "
                             f"{a['multiplier']} |")
                L.append("")
                L.append("*A partial duty cycle dilutes the population-level "
                         "effect in asynchronous cells: per-cell efficacy can be "
                         "much higher than a bulk RT-qPCR suggests. Synchronise or "
                         "use a phase-resolved readout before calling it a "
                         "failure.*")
                L.append("")

        ups = res.get("upstream") or {}
        if ups:
            L.append("### Upstream routes to factors with no direct drug")
            L.append("")
            L.append("Most TFs have no approved ligand. These paths walk upstream "
                     "through signed causal edges — kinases, phosphatases, "
                     "receptors — to proteins that *are* druggable. Signs "
                     "multiply along the path, so an inhibitor of an inhibitor "
                     "counts as activating.")
            L.append("")
            L.append("Distance buys druggability and spends selectivity: a node "
                     "three hops up reaches your gene through one branch and "
                     "everything else through the others. Read the fan-out column.")
            L.append("")
            for sym, hits in list(ups.items())[:6]:
                L.append(f"**{sym}** (no direct pharmacology)")
                L.append("")
                L.append("| Upstream target | Hops | Fan-out | Path | Act on it | "
                         "Drugs | Score |")
                L.append("|---|---|---|---|---|---|---|")
                for h in hits:
                    bp = h["best_path"]
                    act = "activate" if h["net_sign"] == 1 else "inhibit"
                    names = ", ".join(d["name"] for d in h["drugs"][:3]) or "—"
                    L.append(f"| {h['node']} | {bp['hops']} | {bp['fan_out']} | "
                             f"{bp['description']} | {act} | {names} | "
                             f"{h['score']} |")
                L.append("")
                warned = [w for h in hits for w in h["warnings"]]
                if warned:
                    for w in sorted(set(warned))[:3]:
                        L.append(f"- ⚠ {w}")
                    L.append("")

        mechs = res.get("mechanisms") or {}
        if mechs:
            L.append("### Mechanism options: orthosteric vs allosteric")
            L.append("")
            L.append("Allosteric sites are not constrained to bind the endogenous "
                     "ligand, so they diverge across a receptor family where the "
                     "orthosteric pocket cannot. That usually buys subtype "
                     "selectivity — at the cost of a bounded maximal effect and a "
                     "hard dependence on endogenous tone.")
            L.append("")
            for sym, cmp in list(mechs.items())[:5]:
                L.append(f"**{sym}**")
                L.append("")
                L.append("| Mechanism | Drugs | Max effect | Subtype selectivity | "
                         "Outcompeted by ligand? | Needs endogenous ligand | "
                         "Examples |")
                L.append("|---|---|---|---|---|---|---|")
                for r in cmp["by_mechanism"]:
                    star = " ⭐" if r["mechanism"] in ("PAM", "NAM") else ""
                    L.append(f"| {r['mechanism']}{star} | {r['n_drugs']} | "
                             f"{r['max_effect']} | {r['subtype_selectivity']} | "
                             f"{'**yes**' if r.get('surmountable') else 'no'} | "
                             f"{'YES' if r['needs_endogenous_ligand'] else 'no'} | "
                             f"{', '.join(r['examples'])} |")
                L.append("")
                if cmp.get("surmountability"):
                    L.append(f"- {cmp['surmountability']}")
                    L.append("")
                if cmp.get("summary"):
                    L.append(f"*{cmp['summary']}*")
                    L.append("")
            L.append("**Before pursuing any allosteric option, verify that the "
                     "pathway has endogenous tone in your cell type.** A PAM with "
                     "no ligand to modulate does nothing at all, and that is a "
                     "measurement, not a prediction — nothing in this package can "
                     "check it for you.")
            L.append("")

        lax = res.get("ligand_axes") or {}
        if lax:
            L.append("### Changing the natural ligand instead of the receptor")
            L.append("")
            L.append("A nuclear receptor's activity is set by how much of its "
                     "endogenous ligand is present, and the enzymes that make "
                     "or destroy that ligand are often far more druggable than "
                     "the receptor. Anastrozole, finasteride and metyrapone all "
                     "work this way — none of them bind a receptor at all.")
            L.append("")
            L.append("| Strategy | Effect on the receptor |")
            L.append("|---|---|")
            L.append("| supplement ligand | up |")
            L.append("| inhibit catabolism | up |")
            L.append("| inhibit synthesis | down |")
            L.append("| sequester ligand | down |")
            L.append("")

            # Work out which notes every receptor repeats, and say them once.
            _allnotes = [set((a.get("notes") or [])[:4]) for a in lax.values()]
            _shared_ligand_notes = (set.intersection(*_allnotes)
                                    if len(_allnotes) > 1 else set())
            for n2 in sorted(_shared_ligand_notes):
                L.append(f"- {n2}")
            if _shared_ligand_notes:
                L.append("")

            for sym, a in list(lax.items())[:6]:
                L.append(f"#### {sym} — ligand: {a.get('ligand') or 'unknown'}")
                L.append("")
                if a.get("enzymes"):
                    L.append(f"Controlling enzymes: "
                             f"{', '.join(a['enzymes'])}")
                    L.append("")
                L.append("| Drug | Strategy | Acts on | Local? | Status |")
                L.append("|---|---|---|---|---|")
                for r in a["routes"][:top_n]:
                    yr = f" ({r['approval_year']})" if r["approval_year"] else ""
                    L.append(f"| {r['drug']}{yr} | {r['strategy']} | "
                             f"{r['acts_on']} | "
                             f"{'**yes**' if r['tissue_selective'] else 'no'} | "
                             f"{r['status']} |")
                L.append("")

                if a.get("chains"):
                    L.append("The chain, link by link:")
                    L.append("")
                    for c in a["chains"][:top_n]:
                        L.append(f"- {c}")
                    L.append("")

                per_drug = [(r["drug"], r["caveats"]) for r in a["routes"]
                            if r.get("caveats")]
                if per_drug:
                    L.append("Drug-specific caveats:")
                    L.append("")
                    for drug, cav in per_drug[:top_n]:
                        L.append(f"- **{drug}** — {cav[0]}")
                    L.append("")

                if a.get("opposite_direction"):
                    opp = ", ".join(
                        f"{o['drug']} ({o['strategy']})"
                        for o in a["opposite_direction"][:6])
                    L.append(f"*Same axis, opposite direction — do not confuse "
                             f"these with the list above: {opp}.*")
                    L.append("")

                # Notes shared by every receptor on this axis are
                # properties of the axis, not of the receptor, and printing
                # them per receptor is the same text three times.
                _own = [n2 for n2 in a["notes"][:4]
                        if n2 not in _shared_ligand_notes]
                for n2 in _own:
                    L.append(f"- {n2}")
                if _own:
                    L.append("")

            covered = lig_mod.covered_receptors()
            scanned = [t["symbol"] for t in (res.get("tfs") or [])][:25]
            without = [t for t in scanned if t not in covered]
            if without:
                L.append(f"**No ligand axis for:** {', '.join(without[:20])}"
                         + (" …" if len(without) > 20 else ""))
                L.append("")
                L.append("Most transcription factors are ligand-independent, so "
                         "this is usually a real absence rather than a gap in "
                         "the table. Receptors with a curated axis here: "
                         + ", ".join(f"{k} ({v})" for k, v in covered.items())
                         + ".")
                L.append("")

        stab = res.get("stabilisation") or {}
        if stab:
            L.append("### Locking the factor into an active state")
            L.append("")
            L.append("Nuclear receptors have a ligand-binding pocket; almost no "
                     "other TF does, which is why the graded tables keep "
                     "returning NR agonists. Stabilisation sidesteps the pocket: "
                     "every factor has protein-protein interfaces and "
                     "degradation machinery, and both are druggable.")
            L.append("")
            # Split by whether anything actually exists for the factor. The
            # previous version printed the same paragraph and the same four
            # exemplars under every pocket-less TF, which read as though a
            # route were available for each of them. For most, none is.
            actionable = {s: a for s, a in stab.items() if a.get("actionable")}
            nothing = [s for s, a in stab.items() if not a.get("actionable")]

            for sym, a in list(actionable.items())[:6]:
                pocket = ("has a ligand pocket" if a["has_ligand_pocket"]
                          else "**no ligand pocket**")
                L.append(f"**{sym}** — {pocket} · **a curated route exists**")
                L.append("")
                L.append("| Drug | Mechanism | Binds | Raises | Status |")
                L.append("|---|---|---|---|---|")
                for st in a["stabilisers"][:top_n]:
                    yr = f" ({st['approval_year']})" if st["approval_year"] else ""
                    L.append(f"| {st['drug']}{yr} | {st['mechanism']} | "
                             f"{st['target_node']} | {st['raises']} | "
                             f"{st['status']} |")
                L.append("")
                for n2 in a["notes"][:3]:
                    L.append(f"- {n2}")
                L.append("")

            if nothing:
                L.append(f"#### No stabilisation route for "
                         f"{len(nothing)} factor(s)")
                L.append("")
                L.append(", ".join(sorted(nothing)))
                L.append("")
                L.append("These have no ligand pocket, no direct drug, and "
                         "**nothing curated against them or their degradation "
                         "machinery**. For these factors stabilisation is a "
                         "research programme, not a route you can take: "
                         "making it real needs the E3 ligase that degrades "
                         "the factor *and* a ligand for that ligase, and "
                         "neither is published here.")
                L.append("")
                L.append("The mechanism class is real — it has shipped "
                         "approved drugs, listed once below — but those act "
                         "on **other proteins entirely**, and none of them "
                         "is evidence that any factor above is tractable.")
                L.append("")

            if any(a.get("exemplars") for a in stab.values()):
                ex_list = next(a["exemplars"] for a in stab.values()
                               if a.get("exemplars"))
                L.append("<details><summary>Approved proof the mechanism "
                         "class works — on unrelated targets</summary>")
                L.append("")
                for ex in ex_list[:4]:
                    L.append(f"- **{ex['drug']}** ({ex['approval_year']}) — "
                             f"{ex['mechanism']}: {ex['description']}")
                L.append("")
                L.append("</details>")
                L.append("")
            L.append("*Conformational stabilisers raise activity with protein "
                     "levels unchanged; degradation blockers raise levels. A "
                     "western blot cannot see the first kind — pick the assay to "
                     "match the mechanism.*")
            L.append("")

        ct = res.get("cell_type") or {}
        if ct:
            L.append("### Cell-type reach")
            L.append("")
            L.append("Three things must overlap: where the drug goes, where the "
                     "gene is expressed, and where the transcription factor is. "
                     "Any empty intersection kills the route regardless of motif "
                     "quality or potency.")
            L.append("")
            tfs = ct.get("tfs") or {}
            if tfs:
                L.append("| TF | Co-expressed with gene | Shared tissues |")
                L.append("|---|---|---|")
                for sym, ov in list(tfs.items())[:top_n]:
                    if not ov.get("known"):
                        verdict = "unknown"
                    elif ov.get("compatible"):
                        verdict = f"yes ({ov.get('n_shared')})"
                    else:
                        verdict = "**NO**"
                    L.append(f"| {sym} | {verdict} | "
                             f"{', '.join(ov.get('shared') or [])[:60] or '—'} |")
                L.append("")
                if any(o.get("compatible") is False for o in tfs.values()):
                    L.append("*Factors with no shared tissue were heavily "
                             "down-weighted: a hepatocyte factor cannot regulate "
                             "your gene in neurons, whatever its motif score.*")
                    L.append("")
            delivery_rows = ct.get("delivery") or []
            assessed = [d for d in delivery_rows if d.get("cns_mpo") is not None]
            unassessed = [d for d in delivery_rows if d.get("cns_mpo") is None]

            for d in assessed[:top_n]:
                verdict = ("likely" if d.get("likely_penetrant")
                           else "unlikely")
                L.append(f"- **{d['drug']}** — CNS MPO {d['cns_mpo']}/"
                         f"{d['mpo_parameters']}, {verdict} to reach CNS")
                # Notes belong to THIS drug. Printing them outside the
                # per-drug block dumped every failed assessment under
                # whichever drug happened to be printed last.
                seen_notes = []
                for n2 in d.get("notes", []):
                    if n2 not in seen_notes:
                        seen_notes.append(n2)
                for n2 in seen_notes[:2]:
                    L.append(f"  - {n2}")
            if assessed:
                L.append("")

            if unassessed:
                L.append(f"*{len(unassessed)} drug(s) could not be assessed for "
                         "CNS exposure — no structure or physicochemical "
                         "properties were returned for them, which is common "
                         "for entries sourced from DGIdb rather than ChEMBL. "
                         "They are neither passed nor failed on this axis: "
                         + ", ".join(d["drug"] for d in unassessed[:15])
                         + (" …" if len(unassessed) > 15 else "") + ".*")
                L.append("")

        lincs = res.get("lincs") or {}
        if lincs.get("hits"):
            L.append("### Measured expression changes (LINCS L1000)")
            L.append("")
            L.append("**This is the only empirical evidence in the report.** "
                     "Everything else predicts; this measured what happened to "
                     "the transcriptome when cells were dosed. Where it agrees "
                     "with a predicted route, you have a real candidate.")
            L.append("")
            L.append("| Drug | Effect | Weighted | Cell line | Dose | Time | "
                     "Evidence |")
            L.append("|---|---|---|---|---|---|---|")
            for h in lincs["hits"][:top_n]:
                kind = ("measured" if h["is_landmark"] is True
                        else "INFERRED" if h["is_landmark"] is False else "unknown")
                flag = " ⚠non-specific" if h["non_specific"] else ""
                L.append(f"| {h['drug']}{flag} | {h['score']} | "
                         f"{h['weighted_score']} | {h['cell_line'] or '—'} | "
                         f"{h['dose'] or '—'} | {h['timepoint'] or '—'} | {kind} |")
            L.append("")
            if lincs.get("summary"):
                L.append(lincs["summary"])
                L.append("")
            L.append("*L1000 measures ~978 landmark transcripts and infers the "
                     "other ~11,350. If your gene shows INFERRED, that number is a "
                     "model prediction, not a measurement — confirm by RT-qPCR "
                     "before relying on it. Check too that your gene is expressed "
                     "in the cell line listed, or absence of signal means "
                     "nothing.*")
            L.append("")

        gen = res.get("genetics") or {}
        if gen.get("summary") or (gen.get("cis") or {}).get("phenotypes"):
            L.append("### Human genetic evidence")
            L.append("")
            if gen.get("summary"):
                L.append(gen["summary"])
                L.append("")
            phen = (gen.get("cis") or {}).get("phenotypes") or []
            if phen:
                L.append("| Phenotype | Genetic score | Overall |")
                L.append("|---|---|---|")
                for pz in phen[:top_n]:
                    L.append(f"| {pz['phenotype']} | {pz['genetic_score']} | "
                             f"{pz['overall_score']} |")
                L.append("")
            # The actual trans-eQTL answer, not advice to go and find it.
            # This used to print the same "look for trans-eQTLs of X in
            # eQTLGen" paragraph once per factor — identical wording, no
            # result. The significant-results file is 2 MB.
            teq = gen.get("trans_eqtl") or []
            if teq:
                L.append("#### Does variation in the factor actually move "
                         "this gene? (trans-eQTL)")
                L.append("")
                L.append("The strongest evidence available anywhere in this "
                         "report. A motif says a factor *could* bind; a "
                         "trans-eQTL says natural genetic variation near it "
                         "**does** change your gene's expression, in real "
                         "human tissue, with exposure randomised by meiosis.")
                L.append("")
                if gen.get("trans_eqtl_summary"):
                    L.append(gen["trans_eqtl_summary"])
                    L.append("")
                supported = [t for t in teq if t.get("supported")]
                if supported:
                    L.append("| Factor | Variants | Best SNP | Z | p | n |")
                    L.append("|---|---|---|---|---|---|")
                    for t in supported[:top_n]:
                        h = (t.get("hits") or [{}])[0]
                        L.append(f"| **{t['tf']}** | {len(t.get('hits') or [])} | "
                                 f"{h.get('snp', '—')} | {h.get('z', 0):.1f} | "
                                 f"{h.get('p', 0):.2g} | {h.get('n', '—')} |")
                    L.append("")
                else:
                    tested = ", ".join(t["tf"] for t in teq[:8])
                    L.append(f"*No factor tested ({tested}) has trans-eQTL "
                             f"support for this gene.*")
                    L.append("")
                L.append("*eQTLGen is whole blood at n=31,684, FDR<0.05. A "
                         "hit is strong causal evidence. **Absence is weak** "
                         "— trans effects are small and detecting them needs "
                         "enormous power, so a null overwhelmingly means "
                         "underpowered rather than unregulated, and a "
                         "relationship operating outside blood would not "
                         "appear at all.*")
                L.append("")
            else:
                for t in (gen.get("trans") or [])[:3]:
                    if t.get("notes"):
                        L.append(f"- **{t['tf_symbol']} → {t['gene']}**: "
                                 f"{t['notes'][0]}")
                L.append("")

        an = res.get("analogs") or {}
        if an:
            L.append("### Approved drugs that resemble these candidates")
            L.append("")
            L.append(analog_mod.summarise(an))
            L.append("")
            L.append("Everything already annotated at the factor has been "
                     "subtracted. That matters: if a compound looks like a "
                     "known ligand and anyone assayed it there, it is "
                     "already in the tables above — so a raw similarity "
                     "search rediscovers its own input and looks productive.")
            L.append("")
            L.append("| Factor | Candidate | Resembles | Similarity |")
            L.append("|---|---|---|---|")
            for target, blk in an.items():
                for a2 in blk["analogs"][:6]:
                    L.append(f"| {target} | **{a2['drug']}** | "
                             f"{a2['similar_to']} | {a2['similarity']}% |")
            L.append("")
            L.append("*Similarity is not activity. Two molecules at this "
                     "Tanimoto can differ in the one substituent the pocket "
                     "cares about, and for a nuclear receptor a small change "
                     "can flip agonist to antagonist — so even if binding "
                     "holds, the direction is unknown. Each row is a binding "
                     "assay, not a finding.*")
            L.append("")

        bt = res.get("beyond_transcription") or {}
        if bt:
            L.append("### Routes a promoter analysis cannot see")
            L.append("")
            L.append("Everything above acts on transcription. These act "
                     "after it, around it, or on a different gene entirely. "
                     "For some genes that is where the whole problem is: if "
                     "the transcript is destroyed, the allele is silenced, "
                     "or the protein is degraded too fast, more "
                     "transcription does not help and may not be possible.")
            L.append("")

            comp = bt.get("compensation") or {}
            if comp.get("summary"):
                L.append("#### Could a relative do the job instead?")
                L.append("")
                L.append(comp["summary"])
                L.append("")
                for n2 in (comp.get("notes") or [])[:3]:
                    L.append(f"- {n2}")
                L.append("")
                k = comp.get("known_redundancy") or {}
                if k and k.get("paralog") and k["paralog"] != "(none functional)":
                    L.append(f"| | |")
                    L.append(f"|---|---|")
                    L.append(f"| paralog | **{k['paralog']}** |")
                    L.append(f"| why it substitutes | {k['why_it_substitutes']} |")
                    L.append(f"| why it is low | {k['silencing']} |")
                    L.append(f"| precedent | {k['precedent']} |")
                    L.append("")
                    for c in (k.get("caveats") or [])[:2]:
                        L.append(f"- {c}")
                    L.append("")

            react = bt.get("reactivation") or {}
            if react.get("applicable"):
                L.append("#### Is this a silenced allele rather than a low gene?")
                L.append("")
                for n2 in (react.get("notes") or [])[:2]:
                    L.append(f"- {n2}")
                L.append("")
                L.append("| Mechanism | Silencing | Status | Where selectivity comes from |")
                L.append("|---|---|---|---|")
                for r2 in (react.get("routes") or [])[:4]:
                    L.append(f"| {r2['mechanism']} | {r2['silencing_type']} | "
                             f"{r2['status']} | {r2['specificity'][:70]} |")
                L.append("")

            ptx = bt.get("post_transcriptional") or {}
            if ptx.get("routes"):
                L.append("#### After transcription")
                L.append("")
                if ptx.get("summary"):
                    L.append(ptx["summary"])
                    L.append("")
                L.append("| Mechanism | Raises | Status | Applies when |")
                L.append("|---|---|---|---|")
                for r2 in ptx["routes"][:6]:
                    L.append(f"| {r2['mechanism']} | {r2['raises']} | "
                             f"{r2['status']} | {r2['applies_when'][:58]} |")
                L.append("")
                lead = ptx["routes"][0]
                L.append(f"**{lead['mechanism'].capitalize()}** — "
                         f"{lead['what_it_does']}.")
                L.append("")
                for c in (lead.get("caveats") or [])[:2]:
                    L.append(f"- {c}")
                L.append("")

            ps = bt.get("product_stability") or {}
            if ps.get("routes"):
                L.append("#### Raising the protein instead of the transcript")
                L.append("")
                for n2 in (ps.get("notes") or [])[:2]:
                    L.append(f"- {n2}")
                L.append("")
                L.append("| Mechanism | Specificity | Status | Why |")
                L.append("|---|---|---|---|")
                for r2 in ps["routes"]:
                    L.append(f"| {r2['mechanism']} | {r2['specificity']} | "
                             f"{r2['status']} | {r2['specificity_note'][:62]} |")
                L.append("")
                L.append("*Three of these four raise most of the proteome, "
                         "which makes them excellent diagnostics and poor "
                         "therapies. The pharmacological chaperone is the "
                         "exception — and it needs a chaperone to exist for "
                         "your protein.*")
                L.append("")

        imp = res.get("improvements") or {}
        if imp:
            L.append("### Candidates that may beat the current standard")
            L.append("")
            env = imp.get("baseline_envelope")
            base = imp.get("baseline")
            if env:
                L.append(f"**The bar: strength {env['strength']}, specificity "
                         f"{env['specificity']}, efficiency {env['efficiency']}** "
                         f"({env['evidence']}) — the best established value on "
                         f"*each axis*, not one drug's scores. Representative "
                         f"incumbent: {base['name']} ({base['route']}).")
                L.append("")
                L.append(f"{imp['n_previously_linked']} of {imp['n_candidates']} "
                         f"candidates already have a reported link to this gene "
                         f"(CTD returned {imp['ctd_links_found']} chemicals with "
                         f"curated expression effects).")
            else:
                L.append("**No incumbent found.** Nothing in the curated tables or "
                         "CTD currently controls this gene's expression. There is "
                         "no bar to clear, so candidates below are ranked on their "
                         "own merits — check the gap is real rather than a hole in "
                         "the databases.")
            L.append("")
            acc = imp.get("accounting") or {}
            if acc:
                L.append(f"*Accounting for all {acc.get('total')} candidates: "
                         f"{acc.get('previously_linked')} already linked, "
                         f"{acc.get('novel_and_on_front')} novel and non-dominated, "
                         f"{acc.get('novel_but_dominated')} novel but beaten on "
                         f"every axis by another candidate, "
                         f"{acc.get('novel_no_advantage')} with no advantage. "
                         f"All of them are in results.json under "
                         f"`improvements.all_candidates`.*")
                L.append("")

            rl = imp.get("reported_links") or {}
            rc = rl.get("counts") or {}
            if rc.get("total"):
                L.append("#### Compounds already reported to change this gene")
                L.append("")
                L.append("These come from CTD's curated chemical–gene "
                         "interactions: someone measured an expression change "
                         "and a curator recorded it. That makes them the most "
                         "concrete evidence in this report — but the direction "
                         "is what decides whether each is a lead or a warning.")
                L.append("")
                L.append(f"{rc['total']} chemicals total: **{rc['wanted']}** "
                         f"move it the way you want, {rc['opposite']} move it "
                         f"the opposite way, {rc['unclear']} unclear.")
                L.append("")
                L.append("**Expect the 'predicted by this run?' column to say "
                         "no, and do not read that as a failed match.** CTD "
                         "indexes MeSH, which is a vocabulary for the "
                         "toxicology and pharmacology literature rather than "
                         "a drug list, so what comes back is three different "
                         "kinds of thing:")
                L.append("")
                L.append(f"- **chemical** ({rc.get('n_chemical', 0)}) "
                         f"{chr(8212)} a compound. The only kind that can "
                         f"match a candidate.")
                L.append(f"- **exposure** ({rc.get('n_exposure', 0)}) "
                         f"{chr(8212)} formaldehyde, asbestos, bisphenol, "
                         f"tobacco smoke. Real measured effects on your "
                         f"gene, and not things anyone can be given.")
                L.append(f"- **category** ({rc.get('n_category', 0)}) "
                         f"{chr(8212)} MeSH headings such as "
                         f"\"Anti-HIV Agents\" or \"Plant Extracts\". Not "
                         f"compounds at all, so the column is blank rather "
                         f"than no.")
                L.append("")
                matched = rc.get("matched_any_direction", 0)
                L.append(f"Across all three directions, **{matched} of "
                         f"{rc['total']} matched a candidate from this run**. "
                         + ("A low number here is the normal result: this "
                            "report's drug tables are approved pharmacology, "
                            "and most of CTD is not. It means the two "
                            "sources are looking at different parts of "
                            "chemical space, not that either is wrong."
                            if matched <= 5 else
                            "That overlap is the part of CTD this analysis "
                            "independently arrived at."))
                L.append("")

                want = rl.get("wanted_direction") or []
                if want:
                    L.append(f"**Reported to move {q['gene']} "
                             f"{'UP' if goal == 'up' else 'DOWN'}** — leads:")
                    L.append("")
                    L.append("| Chemical | Kind | Predicted by this run? | "
                             "CTD actions |")
                    L.append("|---|---|---|---|")
                    # Compounds first: they are the only rows where the
                    # match column carries information.
                    _order = {"chemical": 0, "exposure": 1, "category": 2}
                    for e in sorted(want,
                                    key=lambda x: _order.get(x.get("kind"), 3)
                                    )[:top_n]:
                        kind = e.get("kind", "chemical")
                        if kind == "category":
                            mark = f"{chr(8212)} (not a compound)"
                        elif e["in_candidates"]:
                            mark = "**yes**"
                        else:
                            mark = "no"
                        L.append(f"| {e['chemical']} | {kind} | {mark} | "
                                 f"{e['actions']} |")
                    if len(want) > top_n:
                        L.append(f"| … {len(want) - top_n} more | | |")
                    L.append("")

                only = rl.get("reported_but_not_predicted") or []
                if only:
                    L.append(f"**{len(only)} compounds have a measured "
                             f"effect on this gene and no mechanism from "
                             f"this analysis.** Exposures and MeSH headings "
                             f"are excluded from that count, because the "
                             f"scan was never going to propose formaldehyde "
                             f"and listing it as a miss says nothing.")
                    L.append("")
                    L.append("Most of what remains is research chemistry "
                             "rather than medicine: tool compounds, probes "
                             "and natural-product fractions that are "
                             "measured because somebody had them on a "
                             "shelf. Read them as evidence the gene IS "
                             "movable, and check the few that are approved "
                             "drugs individually.")
                    L.append("")
                    L.append(", ".join(e["chemical"] for e in only[:25])
                             + (" …" if len(only) > 25 else ""))
                    L.append("")

                opp = rl.get("opposite_direction") or []
                if opp:
                    L.append(f"**Reported to move it the WRONG way** — avoid, "
                             f"or use as negative controls:")
                    L.append("")
                    L.append(", ".join(e["chemical"] for e in opp[:25])
                             + (" …" if len(opp) > 25 else ""))
                    L.append("")
                L.append("*CTD aggregates across cell types, doses and assays. "
                         "A 'both' or contradictory annotation usually means "
                         "the effect is context-dependent rather than that the "
                         "data are wrong — check the underlying study before "
                         "relying on any single row.*")
                L.append("")

            novel = imp.get("novel_improvements") or []
            if novel:
                heading = ("#### Not previously linked, and better on at least "
                           "one axis" if imp.get("has_incumbent")
                           else "#### Not previously linked (no incumbent to beat)")
                L.append(heading)
                L.append("")
                L.append("| Candidate | Route | TF | Strength | Specificity | "
                         "Efficiency | Wins on | Trade-offs |")
                L.append("|---|---|---|---|---|---|---|---|")
                for c in novel[:top_n]:
                    a = c["axes"]
                    L.append(f"| {c['name']} | {c['route']} | {c['tf_symbol'] or '—'} | "
                             f"{a['strength']} | {a['specificity']} | {a['efficiency']} | "
                             f"{', '.join(c['wins_on'])} | "
                             f"{', '.join(c['loses_on']) or 'none'} |")
                L.append("")
                if imp.get("has_incumbent"):
                    L.append("*These are predicted scores compared against "
                             "measured ones, so a margin was required before "
                             "counting a win. They are hypotheses to test, not "
                             "improvements demonstrated.*")
                else:
                    L.append("*No incumbent exists, so these cleared no bar — they "
                             "are simply the candidates that scored above 0.5 on "
                             "at least one axis. Hypotheses to test.*")
            else:
                if imp.get("has_incumbent"):
                    L.append("#### No unlinked candidate beat the standard of care")
                    L.append("")
                    L.append("Every candidate that scored well is already a "
                             "reported regulator of this gene, or failed to clear "
                             "the per-axis bar by more than the prediction error. "
                             "This is the common result and it is informative: it "
                             "says the obvious pharmacology has been explored.")
                else:
                    L.append("#### No candidate scored above 0.5 on any axis")
                    L.append("")
                    L.append("There is no incumbent here, so nothing was filtered "
                             "out for failing to beat one — the candidates "
                             "themselves are simply weak. Widen the promoter "
                             "window or check whether this gene is regulated "
                             "through an enhancer outside it.")
            L.append("")

            dominated = imp.get("novel_dominated") or []
            if dominated:
                L.append("Novel but beaten on every axis by another candidate "
                         "(shown so nothing is hidden): "
                         + ", ".join(c["name"] for c in dominated[:12])
                         + (" ..." if len(dominated) > 12 else ""))
                L.append("")

            front = [c for c in (imp.get("pareto_front") or [])
                     if c.get("wins_on") or c.get("verdict", "").startswith("current")]
            if front:
                L.append("#### Non-dominated set (nothing beats these on every axis)")
                L.append("")
                for c in front[:top_n]:
                    linked = "known link" if c["previously_linked"] else "NEW link"
                    L.append(f"- **{c['name']}** ({linked}) — {c['verdict']}")
                L.append("")

        if res.get("outcome_estimates"):
            L.append("### Predicted transcriptional outcome")
            L.append("")
            L.append("P(gene moves) = P(engages) x P(direction) x P(activates here) "
                     "x P(responsive). Four uncertain links multiply, so these "
                     "numbers are small by construction — a chain at 0.7 per link "
                     "lands near 0.24. Anything reporting 0.9 here has hidden a link.")
            L.append("")
            # Anything named above (the discovery winner, the non-dominated
            # set, the incumbent) must appear here regardless of truncation.
            imp2 = res.get("improvements") or {}
            highlighted = {c["name"].lower() for c in
                           (imp2.get("novel_improvements") or [])
                           + (imp2.get("pareto_front") or [])}
            if imp2.get("baseline"):
                highlighted.add(imp2["baseline"]["name"].lower())

            ests = res["outcome_estimates"]
            pinned = [e for e in ests if e["drug"].lower() in highlighted]
            rest = [e for e in ests if e["drug"].lower() not in highlighted]
            shown = pinned + rest[:max(top_n - len(pinned), 0)]

            L.append("| Drug | TF | P(engage) | P(direction) | P(promoter) | "
                     "P(responsive) | **P(gene moves)** | Verdict | Limiting step |")
            L.append("|---|---|---|---|---|---|---|---|---|")
            for e in shown:
                serm = " ⚠tissue-dependent" if e["tissue_dependent"] else ""
                star = " ★" if e["drug"].lower() in highlighted else ""
                L.append(f"| {e['drug']}{star}{serm} | {e['tf_symbol']} | "
                         f"{e['p_engages']} | "
                         f"{e['p_direction']} | {e['p_activates']} | "
                         f"{e['p_responsive']} | **{e['p_point']}** "
                         f"({e['p_low']}–{e['p_high']}) | {e['verdict']} | "
                         f"{e['limiting_step']} |")
            L.append("")
            if highlighted:
                L.append("★ = named in the comparison against the standard of "
                         "care above, pinned here so the two sections cannot "
                         "disagree.")
                L.append("")
            spread = ([e["p_point"] for e in ests] or [0])
            if max(spread) - min(spread) < 0.05:
                L.append("**These numbers barely differ, and that is the "
                         "honest result.** Two of the four links are the same "
                         "for every candidate: direction defaults to 0.75 "
                         "without `--predict-direction`, and responsiveness is "
                         "fixed at 0.40 because nothing in this package can "
                         "measure it. The chain therefore conveys *magnitude* — "
                         "how uncertain this route is overall — not a ranking. "
                         "Use the grades and the axes above to order "
                         "candidates; use this to decide whether the route is "
                         "worth pursuing at all.")
                L.append("")
            top = shown[0] if shown else res["outcome_estimates"][0]
            L.append(f"**Most informative next experiment** (for {top['drug']}): "
                     f"{top['next_experiment']}")
            L.append("")
            if any(e["tissue_dependent"] for e in res["outcome_estimates"]):
                L.append("*Rows marked tissue-dependent sit in the SERM regime: "
                         "tamoxifen antagonises ERalpha in breast and activates it "
                         "in endometrium. For those, direction is a property of "
                         "your cell type, not of the molecule, and no structural "
                         "model can settle it — a coactivator recruitment assay can.*")
                L.append("")

        all_screens = res.get("screens") or {}
        if all_screens:
            ran = [s for s, v in all_screens.items() if not v.get("error")]
            skipped = [s for s, v in all_screens.items() if v.get("error")]
            L.extend(_screen_scope_note(ran, skipped))

        for sym, sc in all_screens.items():
            L.append(f"### Virtual screen of approved drugs against {sym}")
            L.append("")
            if sc.get("error"):
                L.append(sc["error"])
                L.append("")
                continue
            enr = sc.get("enrichment", {})
            L.append(f"Binding model, CV AUC {sc.get('model_auc')} — trained on "
                     f"{sc.get('n_actives')} measured actives vs "
                     f"{sc.get('n_decoys')} decoys, screened "
                     f"{sc.get('library_size')} approved drugs.")
            if sc.get("direction_note"):
                L.append("")
                L.append(f"*{sc['direction_note']}*")
            if not enr.get("usable", False):
                L.append("")
                L.append(f"**Screen rejected:** {enr.get('reason')} "
                         "No candidates are reported from it.")
                L.append("")
                continue
            L.append("")
            L.append("| Quadrant | Drug | What it is | Approved to treat | "
                     "Measured pChEMBL | P(binds), model | Novelty |")
            L.append("|---|---|---|---|---|---|---|")
            for h in sc["hits"][:top_n]:
                pch = h.get("measured_pchembl")
                prob = h.get("probability")
                # A model score on a molecule the model was trained on is an
                # in-sample fit. It is shown, because hiding it would hide
                # the enrichment diagnostic, but marked so it is not read as
                # a prediction and not mistaken for the ranking key.
                prob_cell = (f"{prob} (in-sample)"
                             if h.get("probability_in_sample") else prob)
                L.append(f"| {h['quadrant']} | {h['drug']} | "
                         f"{h.get('indication_hint') or '—'} | "
                         f"{h.get('treats') or 'not established here'} | "
                         f"{pch if pch is not None else '—'} | "
                         f"{prob_cell if prob is not None else '—'} | "
                         f"{h['novelty'] if h.get('novelty') is not None else '—'} |")
            L.append("")
            L.append("Two independent axes here. **ESTABLISHED LIGAND, "
                     "untested on this gene** is the most actionable class: "
                     "the drug reliably binds this receptor and has no "
                     "reported effect on your gene, so the pharmacology is "
                     "already de-risked and only the receptor-to-gene step "
                     "is a hypothesis. NOVEL HYPOTHESIS rows are hypotheses "
                     "on both counts, novel chemistry and an untested gene "
                     "link, so expect most to be wrong and treat each as "
                     "needing a binding assay first.")
            L.append("")
            L.append("**Read the indication column before the scores.** A "
                     "receptor only reaches this section when it has a "
                     "large ligand series, and such series are rarely "
                     "homogeneous: a topical gel, an oral oncology agent "
                     "and a withdrawn compound can sit on adjacent rows "
                     "with near-identical numbers. The scores rank binding "
                     "to the receptor and say nothing about whether the "
                     "drug could be given for your purpose.")
            L.append("")
            L.append("**The two score columns are not interchangeable, and "
                     "each quadrant has only one real number.** An "
                     "established ligand was a training molecule for this "
                     "model, so its model score is an in-sample fit and its "
                     "novelty is structurally zero, because the nearest "
                     "known ligand the lookup finds is itself. Those rows "
                     "are ranked by measured pChEMBL, which is data. "
                     "Predicted rows have no measurement, so they are ranked "
                     "by the model score, which is a prediction. Comparing a "
                     "pChEMBL to a probability compares a measurement to a "
                     "guess.")
            L.append("")

        du = res.get("direction_unknown") or {}
        if du:
            n = sum(len(v) for v in du.values())
            L.append("### Drug–gene links with no recorded direction")
            L.append("")
            L.append(f"{n} interaction(s) across {len(du)} factor(s) come from "
                     "curated sources that record a drug–gene relationship "
                     "**without saying which way it goes**. They are held back "
                     "from the recommendations above, because the entire "
                     "question here is direction and an undirected link cannot "
                     "answer it.")
            L.append("")
            L.append("They are still leads — someone curated a real "
                     "relationship — but they need a literature check, not a "
                     "prediction. Many are chemotherapeutics or broad stressors "
                     "whose association with a transcription factor reflects "
                     "cell-wide perturbation rather than a specific mechanism.")
            L.append("")
            L.append("| TF | Drug | Reported by | DGIdb score |")
            L.append("|---|---|---|---|")
            rows = [(sym, d) for sym, ds in du.items() for d in ds]
            rows.sort(key=lambda r: -(r[1].get("interaction_score") or 0))
            for sym, d in rows[:top_n]:
                score = d.get("interaction_score")
                L.append(f"| {sym} | {d['name']} | {d.get('sources') or '—'} | "
                         f"{round(score, 2) if score else '—'} |")
            if len(rows) > top_n:
                L.append(f"| … {len(rows) - top_n} more | | | |")
            L.append("")
            L.append("*Pass `--include-unknown-direction` to fold these into "
                     "the main tables anyway — but every such row will read "
                     "UNKNOWN / unclear, which is what they are.*")
            L.append("")

        L.append("### Approved drugs acting on these factors")
        L.append("")
        rows = []
        for t in res["tfs"]:
            for d in t["drugs"]:
                direction = drugs_mod.desired_direction(t["role"] or "unknown", goal)
                want = {"up": "increase", "down": "decrease"}.get(direction, "either")
                rows.append((t["score"], t["symbol"], t["role"], d, want))
        if not rows:
            kind = "activating" if goal == "up" else "blocking"
            L.append(f"No approved drug with {kind} pharmacology was found for any TF "
                     "matched in this promoter. That is the usual result outside the "
                     "nuclear-receptor family; see the README for indirect strategies.")
        else:
            L.append("| TF | TF role | Do what to this TF | Drug | Action | "
                     "Drug direction | First approval | Source |")
            L.append("|---|---|---|---|---|---|---|---|")
            # Untruncated, this ran to 850 rows per promoter — 2,550 lines of
            # a 3,943-line report, and unreadable at any length.
            ordered = sorted(rows, key=lambda r: -r[0])
            for _, sym, role, d, want in ordered[:top_n]:
                flag = " ⚠withdrawn" if d.get("withdrawn") else ""
                L.append(f"| {sym} | {role} | {want} | {d['name']}{flag} | "
                         f"{d['action_type']} | {d['direction']} | "
                         f"{d.get('first_approval') or '—'} | {d['source']} |")
            if len(ordered) > top_n:
                L.append(f"| … {len(ordered) - top_n} more, see "
                         f"tf_summary.tsv | | | | | | | |")
        L.append("")

    # Experimental design is deferred to the end of the report.
    # It is how you would TEST any of the above, and putting a
    # CRISPRa protocol between the reader and the drug list buried
    # the answer under the method.
    lab: list[str] = []
    # These two sections describe the GENE, not a particular promoter window,
    # so they are emitted once. Repeating them per TSS tripled the report.
    for res in result["results"][:1]:
        mo = res.get("modality_design") or {}
        if mo:
            lab.append(f"## Designing a direct modality for {mo['gene']}")
            lab.append("")
            d = mo["delivery"]
            flag = "" if mo["feasible_in_vivo"] else " ⚠"
            lab.append(f"**Delivery to {d['compartment']}: {d['status']}**{flag} — "
                     f"{d['route']}")
            lab.append("")
            lab.append(f"{d['note']}")
            lab.append("")
            lab.append(f"**{mo['verdict']}**")
            lab.append("")
            for opt in mo["options"]:
                lab.append(f"### {opt['modality']}")
                lab.append("")
                lab.append(opt["rationale"])
                lab.append("")
                if opt.get("note"):
                    lab.append(f"*{opt['note']}*")
                    lab.append("")
                cands = opt.get("candidates") or []
                # CRISPRi guides also carry a "sequence" key; the oligo branch
                # is distinguished by having an antisense strand.
                if cands and "antisense" in cands[0]:
                    lab.append("| Target site (5'→3') | Guide/antisense | Pos | "
                             "GC | Notes |")
                    lab.append("|---|---|---|---|---|")
                    for c in cands[:6]:
                        w = "; ".join(c["warnings"])[:45] or "clean"
                        lab.append(f"| `{c['sequence']}` | `{c['antisense']}` | "
                                 f"{c['position']} | {c['gc']:.0%} | {w} |")
                    lab.append("")
                elif cands:
                    lab.append("| Protospacer | PAM | Strand | vs TSS | Notes |")
                    lab.append("|---|---|---|---|---|")
                    for c in cands[:6]:
                        w = "; ".join(c["warnings"])[:40] or "clean"
                        lab.append(f"| `{c['sequence']}` | {c['pam']} | "
                                 f"{c['strand']} | {c['tss_offset']:+d} | {w} |")
                    lab.append("")
            lab.append(f"*{mo['off_target_warning']}*")
            lab.append("")

        sp = next((r.get("screen_plan") for r in result["results"]
                   if r.get("screen_plan")), None)
        if not sp:
            continue
        lab.append("## Screen design")
        lab.append("")
        if len(result["results"]) > 1:
            lab.append(f"*One plan for the gene. The reporter insert below uses "
                     f"the {result['results'][0]['promoter']['transcript_id']} "
                     f"window; clone a second construct for the other TSS if "
                     f"you want to compare them.*")
            lab.append("")
        lab.append("The outcome chain above is almost always limited by *gene "
                 "responsiveness in your cell type* — the link this package "
                 "cannot see. A screen measures it instead of assuming it. "
                 "Everything above prioritises this experiment; it does not "
                 "replace it.")
        lab.append("")
        lab.append("### Reporter construct")
        lab.append("")
        lab.append(f"- Insert: {sp['reporter']['insert']}")
        lab.append(f"- Vector: {sp['reporter']['vector']}")
        lab.append("- Controls:")
        for c in sp["reporter"]["controls"]:
            lab.append(f"  - {c}")
        lab.append(f"- {sp['reporter']['caveat']}")
        lab.append("")
        cr = sp["crispra"]
        lab.append("### CRISPRa screen")
        lab.append("")
        lab.append(f"- Approach: {cr['approach']}")
        lab.append(f"- Focused library: {cr['library']['focused']}")
        lab.append(f"- Unbiased: {cr['library']['genome_wide']}")
        lab.append(f"- Coverage: {cr['library']['coverage']}")
        lab.append(f"- Readout: {cr['readout']}")
        lab.append("")
        gl = cr["guides_at_target_promoter"]
        if gl:
            lab.append(f"Positive-control sgRNAs at your promoter "
                     f"({cr['n_usable_guides']} clean of {len(gl)} shown):")
            lab.append("")
            lab.append("| Protospacer (20 nt) | PAM | Strand | vs TSS | GC | Note |")
            lab.append("|---|---|---|---|---|---|")
            for g in gl[:10]:
                note = g["warnings"][0][:45] if g["warnings"] else "ok"
                lab.append(f"| `{g['sequence']}` | {g['pam']} | {g['strand']} | "
                         f"{g['tss_offset']:+d} | {g['gc']:.0%} | {note} |")
            lab.append("")
        lab.append(f"*{cr['note']}*")
        lab.append("")
        cp = sp["compounds"]
        lab.append("### Compound screen")
        lab.append("")
        for lib in cp["libraries"]:
            lab.append(f"- {lib}")
        lab.append("")
        lab.append(f"- Design: {cp['design']}")
        lab.append(f"- Counter-screen: {cp['counter_screen']}")
        lab.append("")
        if cp["cherry_pick"]:
            n_shown = len(cp["cherry_pick"][:24])
            lab.append("#### Cherry-pick list")
            lab.append("")
            lab.append("**What this is.** \"Cherry-picking\" means having the "
                     "vendor pull these named compounds into a custom plate "
                     "rather than buying a whole library. It is a cost and "
                     "time decision: ~48 wells you can run in a week, against "
                     "~1,520 compounds for Prestwick or ~3,000 for the "
                     "Selleck set.")
            lab.append("")
            lab.append("**Where these came from.** Every compound below was "
                     "proposed earlier in *this run* — nothing external was "
                     "consulted to build the list:")
            lab.append("")
            for src in cp.get("pulled_from", []):
                lab.append(f"- {src}")
            lab.append("")
            if cp.get("priority_key"):
                lab.append("Priority means:")
                lab.append("")
                for p, meaning in sorted(cp["priority_key"].items()):
                    lab.append(f"- **{p}** — {meaning}")
                lab.append("")
            if cp.get("cannot_do"):
                lab.append(f"**What it cannot do.** {cp['cannot_do']}")
                lab.append("")
            total = len(cp["cherry_pick"])
            lab.append(f"Showing {n_shown} of {total} selected "
                     f"(from {cp.get('n_available', total)} candidates in "
                     f"this run), in priority order:")
            lab.append("")
            lab.append("| Priority | Compound | Why it is on the plate | Pulled from |")
            lab.append("|---|---|---|---|")
            for c in cp["cherry_pick"][:24]:
                lab.append(f"| {c.get('priority', '—')} | {c['compound']} | "
                         f"{c['rationale']} | {c.get('source', '—')} |")
            lab.append("")
        lab.append(f"*{cp['note']}*")
        lab.append("")
        lab.append("### What to expect")
        lab.append("")
        for k, v in sp["expectations"].items():
            lab.append(f"- **{k.replace('_', ' ')}**: {v}")
        lab.append("")

    L.append("---")
    L.append("")
    L.append("### How to read this")
    L.append("")
    if goal == "down":
        L.append("**For downregulation, check the direct modalities first.** "
                 "Approved siRNA (patisiran, givosiran, inclisiran, vutrisiran) and "
                 "antisense oligonucleotides (inotersen, nusinersen, tofersen) knock "
                 "down a named transcript with sequence-level specificity. A TF-level "
                 "intervention hits every other gene that factor controls. Use this "
                 "report when a direct modality is not an option, or to understand "
                 "*why* a drug you already have lowers your gene.")
        L.append("")
    L.append("A PWM hit means the sequence *resembles* a binding site. It is a "
             "hypothesis generator, not evidence of binding, and not evidence that "
             "binding would increase transcription of this gene in your cell type. "
             "Rank order here is a heuristic (weights printed in the JSON output), "
             "not a validated predictor. Every row needs wet-lab confirmation: "
             "promoter-reporter assay for the element, agonist dose-response with "
             "RT-qPCR/protein readout for the drug, and a site-directed mutant of the "
             "element to show the drug effect runs through it.")
    L.append("")

    # ------------------------------------- experimental design, at the end
    # Deferred deliberately. These sections are how you TEST the analysis,
    # and they are long: a CRISPRa protocol and a plate layout sitting
    # between the reader and the drug list buried the answer under the
    # method. Everything above is what the model concluded; everything
    # below is how you would find out whether it is right.
    if lab:
        L.append("---")
        L.append("")
        L.append("# Testing any of this in the lab")
        L.append("")
        L.append("Everything above is prediction. This part is experiment — "
                 "how you would find out whether any of it holds, and it is "
                 "placed last because it is long, not because it is "
                 "optional.")
        L.append("")
        L.append("**If you run one thing, run the reporter assay with the "
                 "scrambled-motif control.** A drug that induces the "
                 "wild-type promoter but not the mutant has an effect "
                 "running through the element this report predicted. A drug "
                 "that induces both means the prediction was wrong. That "
                 "single control turns a hit into a mechanism.")
        L.append("")
        L.extend(lab)
    return "\n".join(L)


# ===================================================================== #
# The short version.                                                    #
# ===================================================================== #
# report.md is long on purpose: every number in it is qualified, because
# an unqualified number from this pipeline would be read as stronger than
# it is. That makes it the wrong document to open when you want to know
# what to do next.
#
# This is the other document. It carries the shortlist, how much the
# evidence is worth, and what would change the answer - and it points at
# report.md for everything else rather than restating it in shorter and
# therefore less qualified words. The two rules it follows: say nothing
# here that report.md does not say, and never drop a caveat to make a row
# fit on one line.

SUMMARY_TOP = 10


def _clip(text: str, limit: int = 58) -> str:
    """Trim to a word boundary. A mid-word cut ("renal fai") reads as a
    rendering bug and costs the reader the one fact the column carries."""
    text = (text or "").strip()
    if len(text) <= limit:
        return text
    cut = text[:limit].rsplit(" ", 1)[0].rstrip(",;:")
    return (cut or text[:limit]) + "…"


def _verdict_line(verdict: dict) -> str:
    """How much the ranking is worth, in one sentence."""
    level = (verdict or {}).get("level", "unknown")
    return {
        "strong": "**Evidence: strong.** Measured evidence supports this "
                  "ranking from more than one independent source.",
        "usable": "**Evidence: usable.** At least one measured source "
                  "backs the ranking; read the tiers before acting.",
        "degraded": "**Evidence: degraded.** Read the ranking as a reading "
                    "order, not as evidence. A pillar of the analysis is "
                    "missing.",
        "not interpretable": "**Evidence: NOT INTERPRETABLE.** Too much is "
                             "missing for the ranking to mean anything. "
                             "Fix the inputs before reading the rows.",
    }.get(level, f"**Evidence: {level}.**")


# ===================================================================== #
# One row per drug: which methods found it, and what it should do.       #
# ===================================================================== #
# The route tables are one row per ROUTE, so a drug reached three ways
# appears three times and a reader cannot see that it was reached three
# ways. This inverts that: one row per compound, with the methods that
# found it counted.
#
# That count is the most useful number in the report. These methods share
# very little - a curated ligand annotation, a measured expression
# signature, a regulon statistic, a pocket shape, a chemical similarity -
# so they fail for unrelated reasons. A compound found by three is
# supported in a way no single method's score can express, and the route
# tables were hiding it.
#
# What it is NOT: a vote. Methods differ enormously in strength, so two
# weak agreements do not equal one measured result. The evidence tier
# still decides, and the method count breaks ties within a tier.

# Which methods are genuinely independent of each other. Two rows from
# sources sharing an input are not corroboration, so they collapse to one
# method here: from_grades and from_ligands both read the same ChEMBL
# annotation, and the two signature routes both read the same corpora.
# Keys are the `source` strings catalogue.py actually emits, read from it
# rather than guessed - the first version invented plausible names and
# matched none of them, so every row reported zero methods.
METHOD_OF_SOURCE = {
    "direct": "direct therapy for the gene",
    "drug_scoring": "annotated ligand (ChEMBL/DGIdb/Open Targets)",
    "ligands": "endogenous ligand availability",
    "signatures": "measured signature on your gene",
    "tfactivity": "regulon shift in perturbation corpora",
    "abundance": "measured change in the factor's own transcript",
    "screening": "ligand-based model (trained on this receptor)",
    "docking": "structure-based docking",
    "network": "upstream signalling network",
    "stabilizers": "protein stability / degradation",
    "catalogue": "chromatin-level agent",
}

# How much a method's agreement is worth, for ordering only. A measured
# result on your own gene is not the same kind of thing as a similarity
# search, and a flat count would say it was.
METHOD_WEIGHT = {
    "direct therapy for the gene": 1.0,
    "measured signature on your gene": 0.9,
    "annotated ligand (ChEMBL/DGIdb/Open Targets)": 0.7,
    "regulon shift in perturbation corpora": 0.6,
    "upstream signalling network": 0.5,
    "measured change in the factor's own transcript": 0.5,
    "endogenous ligand availability": 0.5,
    "ligand-based model (trained on this receptor)": 0.4,
    "protein stability / degradation": 0.4,
    "structure-based docking": 0.3,
    "chromatin-level agent": 0.1,
}


def unmapped_sources(result: dict) -> set:
    """Route sources with no method name, so a new source cannot be added
    without appearing in this table. Pinned by a test."""
    seen = set()
    for res in result.get("results", []):
        for r in (res.get("all_routes") or {}).get("routes") or []:
            src = r.get("source")
            if src and src not in METHOD_OF_SOURCE:
                seen.add(src)
    return seen


def _expected_change(routes: list, goal: str, outcomes: dict,
                     measured: dict) -> dict:
    """What should happen to the target gene, with every number it has.

    Three quantities are available and they answer different questions.

    `p_move` is the probability the gene moves at all, from outcome.py's
    four-step chain, with the interval it carries. It is a real computed
    number. It is also nearly constant across candidates, because two of
    the four steps are fixed - direction defaults to 0.75 without
    --predict-direction and responsiveness is pinned at 0.40 since nothing
    here can measure it - so it conveys how uncertain this CLASS of route
    is, not which candidate is better.

    `effect` is a measured signed effect size on your gene, present only
    where a corpus supplied one. A set-membership hit from a .gmt carries
    no magnitude at all and is reported as such rather than as its nominal
    placeholder score.

    What is deliberately absent is a predicted fold-change. Deriving one
    would need the share of the gene's transcription driven by that factor,
    the factor's occupancy change at the achieved dose, and the dose
    reached in your tissue. None of the three is available, so a number
    would be fabricated.
    """
    want = "increase" if goal == "up" else "decrease"
    drug = (routes[0].get("drug") or "").strip().lower()

    out = {"direction": want, "p_move": None, "p_low": None, "p_high": None,
           "effect": None, "effect_basis": None, "chain": None}

    o = outcomes.get(drug)
    if o:
        out["p_move"] = o.get("p_point")
        out["p_low"] = o.get("p_low")
        out["p_high"] = o.get("p_high")
        out["limiting_step"] = o.get("limiting_step")
        out["chain"] = {
            "engages the target": o.get("p_engages"),
            "right direction": o.get("p_direction"),
            "acts on the promoter": o.get("p_activates"),
            "gene is responsive": o.get("p_responsive"),
        }

    m = measured.get(drug)
    if m is not None:
        out["effect"] = m["score"]
        out["effect_basis"] = m["basis"]

    if m is not None:
        out["basis"] = ("MEASURED: this compound moved your gene in a "
                        "cell-based assay, in the direction you want")
    elif any(r.get("source") == "signatures" for r in routes):
        out["basis"] = ("MEASURED: this compound moved your gene in a "
                        "cell-based assay, in the direction you want")
    elif any(r.get("source") == "direct" for r in routes):
        out["basis"] = ("established therapy acting on the gene or its "
                        "product")
    else:
        out["basis"] = ("PREDICTED from the factor's role and the drug's "
                        "action; no measurement of your gene supports it")
    return out


def _measured_effects(result: dict) -> dict:
    """{drug: signed effect on the gene} from whatever corpus supplied one.

    A `.gmt` hit is set membership and carries no magnitude - its score is
    a nominal placeholder - so it is reported as a direction without a
    number rather than as that placeholder.
    """
    from . import signatures as sig
    from .congeners import salt_free_name

    out: dict = {}
    for res in result.get("results", []):
        for h in ((res.get("lincs") or {}).get("hits") or []):
            name = salt_free_name((h.get("drug") or "").strip().lower())
            if not name:
                continue
            score = h.get("score")
            if score is None:
                continue
            nominal = abs(float(score)) == sig.GMT_SET_SCORE
            if nominal:
                out.setdefault(name, {
                    "score": None,
                    "basis": ("set membership in a consensus signature: the "
                              "compound moved this gene, with no magnitude "
                              "attached")})
            else:
                prev = out.get(name)
                if prev is None or prev["score"] is None or \
                        abs(score) > abs(prev["score"]):
                    out[name] = {
                        "score": round(float(score), 2),
                        "basis": ("measured signed effect on this gene in a "
                                  "perturbation corpus (z-like, not a "
                                  "fold-change)")}
    return out

def corroboration(result: dict, top_n: int = 20) -> list:
    """One row per drug across the top routes, with the methods that found it."""
    from .congeners import salt_free_name
    from . import indications as ind_mod, systemic as sys_mod

    goal = result["query"].get("goal", "up")
    routes = []
    outcomes: dict = {}
    for res in result["results"]:
        routes.extend((res.get("all_routes") or {}).get("routes") or [])
        for o in (res.get("outcome_estimates") or []):
            key = (o.get("drug") or "").strip().lower()
            # Keep the most favourable route's estimate for a drug reached
            # more than one way; the row reports the best available path.
            if key and (key not in outcomes
                        or (o.get("p_point") or 0)
                        > (outcomes[key].get("p_point") or 0)):
                outcomes[key] = o
    if not routes:
        return []
    measured = _measured_effects(result)

    # Properties for the BBB call, keyed salt-free so a salt form in the
    # route table finds the parent's structure.
    props: dict = {}
    for res in result["results"]:
        for gr in (res.get("drug_grades") or []):
            key = salt_free_name((gr.get("drug_name") or "").lower())
            pr = (gr.get("evidence") or {}).get("properties")
            if key and pr and key not in props:
                props[key] = pr

    # Rank by tier then score, and take the top routes to summarise.
    routes.sort(key=lambda r: (r.get("tier", 9), -(r.get("score") or 0)))

    grouped: dict = {}
    order: list = []
    for r in routes:
        stem = salt_free_name((r.get("drug") or "").lower())
        if not stem:
            continue
        if stem not in grouped:
            if len(order) >= top_n:
                continue          # the table is full; later routes only
                                  # add methods to drugs already in it
            order.append(stem)
            grouped[stem] = {"drug": r.get("drug"), "routes": []}
        grouped[stem]["routes"].append(r)

    out = []
    for stem in order:
        g = grouped[stem]
        rs = g["routes"]
        best = min(rs, key=lambda r: r.get("tier", 9))

        methods = []
        for r in rs:
            m = METHOD_OF_SOURCE.get(r.get("source"))
            if m and m not in methods:
                methods.append(m)

        desc = ind_mod.describe(
            g["drug"], mechanism=(best.get("mechanism") or "").strip() or None)
        cls = sys_mod.classify(g["drug"])
        _pr = props.get(salt_free_name((g["drug"] or "").lower()))
        bbb = dlv_mod.bbb_verdict(g["drug"], _pr)
        # Assay-interference structures undermine the EVIDENCE behind a
        # predicted route, which is what liabilities.py was written for
        # and was never connected to anything. It never demotes an
        # approved drug for being approved - only the confidence in a
        # prediction about it.
        liab = liab_mod.assess(g["drug"], (_pr or {}).get("canonical_smiles"))
        exp = _expected_change(rs, goal, outcomes, measured)

        out.append({
            "drug": g["drug"],
            "n_methods": len(methods),
            "methods": methods,
            "weight": round(sum(METHOD_WEIGHT.get(m, 0.2)
                                for m in methods), 2),
            "best_tier": best.get("tier"),
            "acts_on": sorted({r.get("acts_on") for r in rs
                               if r.get("acts_on")})[:3],
            "categories": sorted({r.get("category") for r in rs}),
            "mechanism": (best.get("mechanism") or "").strip(),
            "what_it_is": desc["what_it_is"],
            "treats": desc["treats"],
            "systems": (cls.name if cls else None),
            "systems_why": (cls.why if cls else None),
            "bbb": bbb["cell"],
            "bbb_why": bbb["why"],
            "route": desc.get("route"),
            "route_note": desc.get("route_note"),
            "generic_available": desc.get("generic_available"),
            "pains": bool(liab.pains),
            "evidence_multiplier": liab.evidence_multiplier,
            "expected_change": exp["direction"],
            "expected_basis": exp["basis"],
            "p_move": exp["p_move"],
            "p_low": exp["p_low"],
            "p_high": exp["p_high"],
            "effect_size": exp["effect"],
            "effect_basis": exp["effect_basis"],
            "chain": exp.get("chain"),
            "limiting_step": exp.get("limiting_step"),
            "n_routes": len(rs),
        })

    out.sort(key=lambda r: (r["best_tier"], -r["weight"], -r["n_methods"]))
    return out


def corroboration_table(result: dict, top_n: int = 20) -> list:
    """The cross-method table, as markdown lines."""
    rows = corroboration(result, top_n=top_n)
    if not rows:
        return []
    goal = result["query"].get("goal", "up")
    gene = result["query"]["gene"]

    L = [f"## All {len(rows)} candidates, and which methods found each", ""]
    L.append("One row per compound rather than per route, so a drug reached "
             "several ways is visible as one candidate with several "
             "supports.")
    L.append("")
    L.append(f"| Drug | Methods | Tier | Acts on | What it is | Systems it "
             f"acts on | Reach | BBB? | {gene} change | P(moves) | "
             f"Measured effect |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|")
    any_p, any_effect = False, False
    for r in rows:
        sys_cell = r["systems"] or "no systemic class flagged"
        what = r["what_it_is"] or r["mechanism"] or "\u2014"
        tag = ("measured" if "MEASURED" in r["expected_basis"]
               else "established" if "established" in r["expected_basis"]
               else "predicted")

        p = r.get("p_move")
        if p is None:
            p_cell = "not estimated"
        else:
            any_p = True
            lo, hi = r.get("p_low"), r.get("p_high")
            p_cell = (f"{p:.2f} [{lo:.2f}\u2013{hi:.2f}]"
                      if lo is not None and hi is not None else f"{p:.2f}")

        eff = r.get("effect_size")
        if eff is not None:
            any_effect = True
            eff_cell = f"{eff:+.2f}"
        elif r.get("effect_basis"):
            eff_cell = "direction only"
        else:
            eff_cell = "\u2014"

        L.append(f"| **{r['drug']}** | {r['n_methods']} | "
                 f"{r['best_tier']} | "
                 f"{', '.join(r['acts_on']) or EMDASH} | "
                 f"{_clip(what, 34)} | {_clip(sys_cell, 28)} | "
                 f"{_route_cell(r)} | "
                 f"{r.get('bbb') or '—'} | "
                 f"{r['expected_change']} ({tag}) | {p_cell} | "
                 f"{eff_cell} |")
    L.append("")

    # What the two numeric columns are, and what neither of them is.
    L.append(f"**There is no predicted fold-change here, and that is not an "
             f"omission.** Deriving one would need the share of {gene}'s "
             f"transcription driven by that factor, the factor's occupancy "
             f"change at the dose achieved, and the dose reached in your "
             f"tissue. None of the three is available to this package, so "
             f"any number in a fold-change column would be fabricated. The "
             f"two numbers below are the ones it can actually compute.")
    L.append("")
    if any_p:
        lo_p = min(r["p_move"] for r in rows if r.get("p_move") is not None)
        hi_p = max(r["p_move"] for r in rows if r.get("p_move") is not None)
        L.append(f"**P(moves)** is the probability {gene} moves at all, "
                 f"with its interval - the product of four steps: the drug "
                 f"engaging the target, acting in the right direction, "
                 f"acting through the promoter, and the gene being "
                 f"responsive in your cells.")
        L.append("")
        L.append(f"Read it as a magnitude of uncertainty, **not as a "
                 f"ranking**. Every value here sits between {lo_p:.2f} and "
                 f"{hi_p:.2f}, and that narrow band is not a coincidence: "
                 f"two of the four steps are the same constant for every "
                 f"row. Direction defaults to 0.75 without "
                 f"`--predict-direction`, and responsiveness is pinned at "
                 f"0.40 because nothing in this package can measure it. So "
                 f"the number tells you how uncertain this CLASS of route "
                 f"is; order on tier and method count instead.")
        L.append("")
        limits = {}
        for r in rows:
            st = r.get("limiting_step")
            if st:
                limits[st] = limits.get(st, 0) + 1
        if limits:
            top = max(limits.items(), key=lambda kv: kv[1])
            L.append(f"For {top[1]} of {len(rows)} rows the limiting step is "
                     f"the same: **{top[0]}**. That is where the probability "
                     f"is lost, and the experiment that raises it is the one "
                     f"worth doing first.")
            L.append("")
    if any_effect:
        L.append("**Measured effect** is a signed effect size on your gene "
                 "from a perturbation corpus, z-like rather than a "
                 "fold-change. `direction only` means the corpus recorded "
                 "that the compound moved the gene without a magnitude, "
                 "which is what a consensus gene-set membership carries.")
        L.append("")

    multi = [r for r in rows if r["n_methods"] > 1]
    if multi:
        L.append(f"**{len(multi)} of {len(rows)} were found by more than one "
                 f"method.** These methods share little - a curated ligand "
                 f"annotation, a measured expression signature, a regulon "
                 f"statistic, a pocket shape - so they fail for unrelated "
                 f"reasons, and agreement between them is worth more than a "
                 f"better score from any one.")
        L.append("")
        for r in multi[:8]:
            L.append(f"- **{r['drug']}**: {'; '.join(r['methods'])}")
        L.append("")
    else:
        L.append(f"**No compound was found by more than one method.** With "
                 f"{len({m for r in rows for m in r['methods']})} method(s) "
                 f"contributing, that is expected rather than informative: "
                 f"corroboration needs two methods to have run and both to "
                 f"reach the same compound. `report.md` lists the flags that "
                 f"turn the others on.")
        L.append("")

    L.append("**The method count is not a vote.** These methods differ "
             "enormously in strength, so two weak agreements do not equal "
             "one measured result. Tier still decides the order and the "
             "count breaks ties within a tier.")
    L.append("")
    bbb_counts: dict = {}
    for r in rows:
        c = r.get("bbb") or "—"
        bbb_counts[c] = bbb_counts.get(c, 0) + 1
    L.extend(_bbb_footnote(bbb_counts))

    local = [r for r in rows if r.get("route") == "local only"]
    if local:
        n_l = len(local)
        L.append(f"**Reach** is the FDA route of administration, from the "
                 f"Orange Book. **{n_l} "
                 f"{'candidate here is' if n_l == 1 else 'candidates here are'}"
                 f" approved only by local routes** "
                 f"({', '.join(r['drug'] for r in local[:6])}) - topical, "
                 f"ophthalmic or inhaled. Those cannot reach a target "
                 f"outside the tissue they are applied to, however they "
                 f"score above, so for an internal gene they are not "
                 f"usable without reformulation.")
        L.append("")
    elif any(r.get("route") for r in rows):
        L.append("**Reach** is the FDA route of administration. Every "
                 "candidate here has a systemic route, so none is ruled "
                 "out for being topical-only.")
        L.append("")
    else:
        L.append("**Reach** is blank for every row because no Orange Book "
                 "is on disk. That is unknown route, not systemic - a "
                 "topical-only drug looks identical. See "
                 "`orangebook.SOURCE` for the download.")
        L.append("")

    pains = [r for r in rows if r.get("pains")]
    if pains:
        n_p = len(pains)
        L.append(f"**{n_p} "
                 f"{'candidate matches' if n_p == 1 else 'candidates match'} "
                 f"a pan-assay interference pattern** "
                 f"({', '.join(r['drug'] for r in pains[:5])}). They are "
                 f"approved drugs and this says nothing about their "
                 f"safety. What it undermines is the EVIDENCE: these "
                 f"structures read active across unrelated assays through "
                 f"aggregation, redox cycling or fluorescence "
                 f"interference, so a measured or predicted activity for "
                 f"them is likelier to be an artefact than the tier "
                 f"suggests.")
        L.append("")

    L.append(f"**The {gene} change column is a direction.** `measured` "
             f"means a compound moved {gene} in a cell-based assay; "
             f"`established` means it is an approved therapy acting on the "
             f"gene or its product; `predicted` means it follows from the "
             f"factor's role and the drug's action, with no measurement of "
             f"{gene} behind it.")
    L.append("")
    return L


def summary_report(result: dict, top_n: int = SUMMARY_TOP) -> str:
    """A one-page summary. Everything here is also in report.md."""
    q = result["query"]
    goal = q.get("goal", "up")
    verb = "raise" if goal == "up" else "lower"
    L: list[str] = []

    L.append(f"# {q['gene']}: how to {verb} it - summary")
    L.append("")
    L.append(f"Objective: **{verb} {q['gene']}**. "
             f"{'Demo data - not a real analysis. ' if q.get('demo') else ''}"
             f"Full detail, with every qualification, in `report.md`.")
    L.append("")

    verdict = result.get("evidence_verdict") or {}
    L.append(_verdict_line(verdict))
    for r in (verdict.get("reasons") or [])[:3]:
        L.append(f"- {r}")
    # Only when the verdict's own reasons have not already said it. The
    # first draft printed the same 78% twice, two lines apart.
    fdp = verdict.get("false_discovery_proportion")
    said = any("chance" in r or "expected to be" in r
               for r in (verdict.get("reasons") or []))
    if fdp is not None and not said:
        L.append(f"- Of the motif sites called, about **{fdp:.0%} are "
                 f"expected to be chance matches**.")
    L.append("")

    # ---- the shortlist
    preds = []
    for res in result["results"]:
        preds.extend((res.get("all_routes") or {}).get("routes") or [])
    seen, rows = set(), []
    for r in sorted(preds, key=lambda x: (x.get("tier", 9),
                                          -(x.get("score") or 0))):
        key = (r.get("drug", "").lower(), r.get("acts_on"))
        if key in seen:
            continue
        seen.add(key)
        rows.append(r)
        if len(rows) >= top_n:
            break

    if rows:
        from . import indications as ind_mod
        L.append(f"## Top {len(rows)} to look at first")
        L.append("")
        L.append("Ordered by evidence tier, then score. A tier-1 route with "
                 "a mediocre score beats a tier-5 route with a great one.")
        L.append("")
        L.append("| # | Drug | Approved to treat | Acts on | Tier | Why it "
                 "is here |")
        L.append("|---|---|---|---|---|---|")
        for i, r in enumerate(rows, 1):
            d = ind_mod.describe(
                r.get("drug", ""),
                mechanism=(r.get("mechanism") or "").strip() or None)
            treats = d["treats"] or "not established here"
            L.append(f"| {i} | **{r.get('drug')}** | {_clip(treats)} | "
                     f"{r.get('acts_on') or EMDASH} | {r.get('tier')} | "
                     f"{r.get('category')} |")
        L.append("")
        L.append("**Every row is off-label.** The third column is what the "
                 "drug was approved for; using it to move your gene is a "
                 "different purpose, and the distance between the two is "
                 "part of what you are judging.")
        L.append("")
    else:
        L.append("## No route was found")
        L.append("")
        L.append("Nothing cleared the thresholds. `report.md` has the "
                 "per-stage accounting for where candidates were lost.")
        L.append("")

    L.extend(corroboration_table(result, top_n=20))

    # ---- what ran and what did not: the thing that decides how to read it
    sources = ((result.get("evidence_sources") or {}).get("sources") or {})
    if sources:
        ok = [k for k, v in sources.items() if v.get("state") == "ok"]
        failed = [k for k, v in sources.items()
                  if v.get("state") not in ("ok", None)
                  and v.get("requested")]
        if ok or failed:
            L.append("## What the answer rests on")
            L.append("")
            if ok:
                L.append(f"- **Ran:** {', '.join(sorted(ok))}")
            if failed:
                L.append(f"- **Did not run:** {', '.join(sorted(failed))}. "
                         f"A source that failed is not a negative result; "
                         f"`report.md` names the flag that would enable "
                         f"each.")
            L.append("")
        elif not q.get("demo"):
            L.append("## What the answer rests on")
            L.append("")
            L.append("No optional evidence source was requested, so this "
                     "ranking is motif prediction and drug annotation "
                     "only. `report.md` lists the flags that add measured "
                     "evidence.")
            L.append("")

    # ---- the things most likely to mislead, from the block itself
    cautions = []
    for res in result["results"]:
        at = res.get("attrition") or {}
        if at.get("n_drugged_nuclear_receptor"):
            cautions.append(
                f"The drug tables are dominated by nuclear receptors "
                f"because {at.get('n_with_drug', 0)} of "
                f"{at.get('n_tfs_with_sites', 0)} factors with a site here "
                f"have any ligand at all. That is attrition, not a "
                f"ranking.")
        if at.get("n_top25_without_drug"):
            names = ", ".join((at.get("top25_without_drug") or [])[:4])
            cautions.append(
                f"{at['n_top25_without_drug']} of the best-scoring factors "
                f"have no approved ligand ({names}) and appear in no drug "
                f"table. The gene-specific signal is partly in them.")
        act = res.get("tf_activity") or {}
        if act.get("n_testable") and not [
                a for a in (act.get("scanned") or [])
                if _sig_q(a.get("q_value")) < 0.05]:
            cautions.append(
                "No factor with a site in this promoter has a detectably "
                "perturbed regulon in the tissue checked, so these routes "
                "rest on sequence alone.")
        ov = result.get("tss_window_overlap") or []
        if ov:
            cautions.append(
                "The promoter windows overlap each other, so counts across "
                "them are not independent.")
        break

    if cautions:
        L.append("## Read these before the table")
        L.append("")
        for c in cautions[:4]:
            L.append(f"- {c}")
        L.append("")

    # ---- what to do next
    L.append("## Next step")
    L.append("")
    L.append(f"Every row above is a hypothesis. The experiment that settles "
             f"one is the same in each case: a promoter-reporter assay for "
             f"the element, a dose-response with RT-qPCR or protein readout "
             f"for the drug, and a site-directed mutant of the element to "
             f"show the drug effect runs through it. Without the third, a "
             f"change in {q['gene']} does not establish the mechanism this "
             f"report proposes.")
    L.append("")
    L.append("*Not medical advice. Hypotheses for laboratory testing.*")
    L.append("")
    return "\n".join(L)


def write_outputs(result: dict, outdir: str | Path, *, top_n: int = 25) -> list[Path]:
    """Write every output; top_n limits report tables, not exported data."""
    if isinstance(top_n, bool) or not isinstance(top_n, int) or top_n < 1:
        raise ValueError("top_n must be a positive integer")
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    written = []

    p = outdir / "report.md"
    p.write_text(markdown_report(result, top_n=top_n), encoding="utf-8")
    written.append(p)

    # The short version. report.md is long because every number in it is
    # qualified; this is the page to open first.
    p = outdir / "summary.md"
    p.write_text(summary_report(result), encoding="utf-8")
    written.append(p)

    p = outdir / "results.json"
    p.write_text(json.dumps(result, indent=2, default=str), encoding="utf-8")
    written.append(p)

    p = outdir / "tf_summary.tsv"
    cols = ["symbol", "name", "role", "nuclear_receptor", "n_sites", "best_rel_score",
            "best_score_bits", "best_site_tss_offset", "best_site_sequence",
            "matrix_id", "chip_peaks", "score", "n_drugs", "drugs"]
    lines = ["\t".join(cols)]
    for res in result["results"]:
        for t in res["tfs"]:
            names = "; ".join(d["name"] for d in t["drugs"])
            row = [t.get(c) for c in cols[:-2]] + [len(t["drugs"]), names]
            lines.append("\t".join("" if v is None else str(v) for v in row))
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    written.append(p)

    # BED of predicted sites, for genome-browser inspection.
    p = outdir / "sites.bed"
    lines = ['track name="predicted_TFBS" description="PWM hits" useScore=1']
    for res in result["results"]:
        prom = res["promoter"]
        for h in res["hits"]:
            strand = h["strand"]
            if prom["strand"] == 1:
                gs = prom["start"] - 1 + h["start"]
                ge = prom["start"] - 1 + h["end"]
            else:
                gs = prom["end"] - h["end"]
                ge = prom["end"] - h["start"]
                strand = {"+": "-", "-": "+"}.get(strand, strand)
            score = max(0, min(1000, int(h["rel_score"] * 1000)))
            lines.append("\t".join([
                f"chr{prom['chrom']}", str(gs), str(ge),
                f"{h['motif_name']}|{h['matrix_id']}", str(score), strand,
            ]))
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    written.append(p)

    p = outdir / "promoter.fa"
    fa = []
    for res in result["results"]:
        prom = res["promoter"]
        fa.append(f">{prom['gene_symbol']}|{prom['transcript_id']}|"
                  f"chr{prom['chrom']}:{prom['start']}-{prom['end']}|"
                  f"strand={prom['strand']}|-{prom['upstream']}/+{prom['downstream']}")
        seq = res["sequence"]
        fa += [seq[i:i + 60] for i in range(0, len(seq), 60)]
    p.write_text("\n".join(fa) + "\n", encoding="utf-8")
    written.append(p)

    return written
