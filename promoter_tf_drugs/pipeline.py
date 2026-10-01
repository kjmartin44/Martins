"""Orchestration: gene -> promoter -> motifs -> activators -> approved drugs."""

from __future__ import annotations

import json
import logging
import math
from dataclasses import dataclass, field
from pathlib import Path

from . import (annotate, catalogue as cat_mod, cellcycle, congeners as cong_mod,
               delivery as dlv_mod,
               demo_data, direct, discovery, drug_scoring, eqtl as eqtl_mod,
               modalities as modal_mod,
               drugs as drugs_mod, evidence, genes, motifs as motifs_mod,
               genetics as gen_mod, ligands as lig_mod, modulators as mod_mod,
               neighbours as nbr_mod, network as net_mod,
               predictions as pred_mod, systemic as sys_mod,
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
            "drugs": [d.to_dict() for d in self.drugs],
        }


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
        lincs_gmt: tuple | None = None, expression_file: str | None = None,
        promoter_fasta: str | None = None,
        gtf: str | None = None,
        use_genetics: bool = False,
        upstream_network: bool = False,
        upstream_depth: int = 3,
        approved_only: bool = True, use_dgidb: bool = True,
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
    if lincs_gmt:
        sig_mod.DEFAULT_GMT_PATHS['up'] = (lincs_gmt[0],)
        sig_mod.DEFAULT_GMT_PATHS['down'] = (lincs_gmt[1],)
    if gtf:
        # Removes the neighbour check's network dependency.
        nbr_mod.load_gtf(gtf)
    if promoter_fasta:
        # Removes the one stage that had no offline route.
        genes.load_promoter_fasta(promoter_fasta)

    # ------------------------------------- 0. drugs that skip the TF route
    # Done first because if one of these exists it outranks everything below.
    direct_hits = []
    if include_direct:
        uniprot = None
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
                use_dgidb=use_dgidb, goal=goal,
                allow_unknown_direction=allow_unknown_direction)
        for s, ds in drug_map.items():
            if s in by_tf:
                by_tf[s].drugs = ds

        # Only *usable* ChIP counts as measured binding. Asking for ChIP and
        # getting a broken response previously scored every factor a hard zero
        # on the 20-point evidence axis — presenting a failed lookup as though
        # binding had been measured and ruled out. A failed source falls back
        # to the unverified discount instead.
        for r in by_tf.values():
            r.score = score_tf(r, has_chip=chip_src.usable, goal=goal)

        ranked = sorted(by_tf.values(), key=lambda r: -r.score)

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
            for r in ranked[:6]:
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

        # ------------------ 8d. empirical: did any drug actually move this gene?
        lincs = {}
        if (use_lincs or lincs_file) and not demo:
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
            cis = _optional("Open Targets genetics",
                            lambda: gen_mod.target_consequences(
                                client, prom.gene_symbol),
                            gen_mod.GeneticEvidence(prom.gene_symbol))
            genetic["cis"] = cis.to_dict()
            genetic["summary"] = gen_mod.summarise(cis)
            _phen = (genetic["cis"] or {}).get("phenotypes") or []
            ledger.record("genetics", requested=True, ok=bool(_phen),
                          yielded=len(_phen),
                          detail=(f"{len(_phen)} associated phenotype(s)"
                                  if _phen else
                                  "no genetic association returned for this "
                                  "gene"))
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

        block_for_catalogue = {
            "direct_therapies": [t.to_dict() for t in direct_hits],
            "drug_grades": [g.to_dict() for g in grades],
            "mechanisms": mechanisms,
            "stabilisation": stabilisation,
            "ligand_axes": ligand_axes,
            "upstream": upstream_hits,
            "lincs": lincs,
            "screens": screens,
        }
        all_routes = cat_mod.build(block_for_catalogue, goal=goal)

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

        L.append("### Top-ranked transcription factors")
        L.append("")
        L.append("| Rank | TF | Role | NR | Sites | Best rel. score | Offset | ChIP | Approved drugs | Score |")
        L.append("|---|---|---|---|---|---|---|---|---|---|")
        for i, t in enumerate(res["tfs"][:top_n], 1):
            L.append(f"| {i} | {t['symbol']} | {t['role']} | "
                     f"{'yes' if t['nuclear_receptor'] else ''} | {t['n_sites']} | "
                     f"{t['best_rel_score']} | {_fmt_offset(t['best_site_tss_offset'])} | "
                     f"{t['chip_peaks']} | {len(t['drugs'])} | {t['score']} |")
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

        if res.get("drug_grades"):
            L.append("### Graded candidates")
            L.append("")
            L.append("Combined = geometric mean of TF confidence and drug quality. "
                     "A weak TF prediction caps the grade no matter how good the drug.")
            L.append("")
            L.append("| Grade | Combined | Drug | TF | TF conf. | Drug qual. | "
                     "Potency | Protein select. | Gene select. | Evidence | Flags |")
            L.append("|---|---|---|---|---|---|---|---|---|---|---|")
            for g in res["drug_grades"][:top_n]:
                c = g["components"]
                flags = "; ".join(g["notes"])[:70] or "—"
                L.append(f"| **{g['grade']}** | {g['combined']} | {g['drug_name']} | "
                         f"{g['tf_symbol']} | {g['tf_confidence']} | {g['drug_quality']} | "
                         f"{c['potency']} | {c['selectivity']} | "
                         f"{c.get('regulon', '—')} | {c.get('evidence_depth', '—')} | "
                         f"{flags} |")
            L.append("")

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

                for n2 in a["notes"][:4]:
                    L.append(f"- {n2}")
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

                want = rl.get("wanted_direction") or []
                if want:
                    L.append(f"**Reported to move {q['gene']} "
                             f"{'UP' if goal == 'up' else 'DOWN'}** — leads:")
                    L.append("")
                    L.append("| Chemical | Predicted by this run? | CTD actions |")
                    L.append("|---|---|---|")
                    for e in want[:top_n]:
                        L.append(f"| {e['chemical']} | "
                                 f"{'yes' if e['in_candidates'] else '**no**'} | "
                                 f"{e['actions']} |")
                    if len(want) > top_n:
                        L.append(f"| … {len(want) - top_n} more | | |")
                    L.append("")

                only = rl.get("reported_but_not_predicted") or []
                if only:
                    L.append(f"**{len(only)} of these were NOT predicted by this "
                             f"run at all.** They have a measured effect on your "
                             "gene and no mechanism from this analysis — which "
                             "makes them worth reading closely, since an "
                             "empirical hit the model missed is either a route "
                             "outside its scope or a clue that the scan is "
                             "wrong.")
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

        for sym, sc in (res.get("screens") or {}).items():
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
            L.append("| Quadrant | Drug | P(wanted direction) | Novelty | Note |")
            L.append("|---|---|---|---|---|")
            for h in sc["hits"][:top_n]:
                L.append(f"| {h['quadrant']} | {h['drug']} | {h['probability']} | "
                         f"{h['novelty']} | {(h['notes'] or [''])[0][:60]} |")
            L.append("")
            L.append("*Two independent axes here. **ESTABLISHED LIGAND, "
                     "untested on this gene** is the most actionable class: "
                     "the drug reliably binds this receptor and has no "
                     "reported effect on your gene, so the pharmacology is "
                     "already de-risked and only the receptor-to-gene step is "
                     "a hypothesis. NOVEL HYPOTHESIS rows are hypotheses on "
                     "both counts — novel chemistry AND an untested gene link "
                     "— so expect most to be wrong and treat each as requiring "
                     "a binding assay first.*")
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


def write_outputs(result: dict, outdir: str | Path, *, top_n: int = 25) -> list[Path]:
    """Write all five outputs; top_n limits report tables, not exported data."""
    if isinstance(top_n, bool) or not isinstance(top_n, int) or top_n < 1:
        raise ValueError("top_n must be a positive integer")
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    written = []

    p = outdir / "report.md"
    p.write_text(markdown_report(result, top_n=top_n), encoding="utf-8")
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
