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
from . import validation as _validation
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


# ===================================================================== #
# What run() accepts.                                                   #
# ===================================================================== #
# fifty-eight arguments, of which about fifteen used to be enforced by a
# block of hand-written checks below. The annotations covered the rest
# and Python never read them, so `upstream="1000"` reached the window
# arithmetic and failed somewhere that looked unrelated.
#
# `why` carries what a bad value DOES, which is the part that is not
# obvious from the rule. See validation.py for why this is a table here
# rather than a pydantic model.
_V = _validation.Spec
_FLAG = _V("bool")

RUN_SPECS: dict = {
    "gene": _V("str", why="a symbol like CCR5 or an Ensembl gene ID."),
    "upstream": _V("int", minimum=0,
                   why="bases before the TSS; with downstream it sets the "
                       "scanned window."),
    "downstream": _V("int", minimum=0,
                     why="bases after the TSS."),
    "pvalue": _V("float", minimum=0, maximum=1, exclusive_min=True,
                 exclusive_max=True,
                 why="per-base motif threshold. At 1.0 every position "
                     "calls a site."),
    "species": _V("str", why="an Ensembl species name."),
    "pfm_file": _V("path", optional=True, must_exist=True,
                   why="JASPAR matrices. Without it the motif set is the "
                       "small bundled one."),
    "chip_bed": _V("path", optional=True, must_exist=True,
                   why="measured occupancy. Its absence is not evidence "
                       "that a factor does not bind."),
    "matrix_ids": _V("list", optional=True),
    "all_tss": _FLAG,
    "mask_repeats": _FLAG,
    "goal": _V("str", choices=("up", "down"),
               why="which way you want the gene to move. It decides "
                   "whether an activator is wanted on or off, so a wrong "
                   "value inverts every recommendation."),
    "grade_drugs": _FLAG,
    "tissue": _V("str", optional=True, choices=("cns",),
                 why="down-weights drugs unlikely to reach the "
                     "compartment."),
    "tf_weight": _V("float", minimum=0, maximum=1,
                    why="how much factor-level evidence counts against "
                        "drug-level evidence."),
    "predict_direction": _FLAG,
    "regulon_file": _V("path", optional=True, must_exist=True,
                       why="TF to target mapping."),
    "include_direct": _FLAG,
    "estimate_outcome": _FLAG,
    "screen_receptors": _FLAG,
    "find_improvements": _FLAG,
    "chip_api": _FLAG,
    "compare_mechanisms": _FLAG,
    "proliferation": _V("str",
                        choices=("cycling", "post-mitotic", "mixed",
                                 "unknown"),
                        why="whether the target cells divide, which "
                            "decides if an antiproliferative is a "
                            "mechanism or a side effect."),
    "allow_unknown_direction": _FLAG,
    "cell_type_check": _FLAG,
    "include_stabilizers": _FLAG,
    "design_screen": _FLAG,
    "design_modality": _FLAG,
    "use_lincs": _FLAG,
    "signature_file": _V("path", optional=True, must_exist=True),
    "lincs_file": _V("path", optional=True, must_exist=True),
    "lincs_geneinfo": _V("any", optional=True),
    "tf_activity": _FLAG,
    "depmap_file": _V("path", optional=True),
    "prism_file": _V("path", optional=True),
    "sensitivity_scan": _FLAG,
    "screen_product": _FLAG,
    "promoter_dna": _FLAG,
    "docking_results": _V("path", optional=True),
    "activity_tissue": _V("str", optional=True),
    "activity_permutations": _V("int", minimum=0,
                                why="null draws for the activity score. "
                                    "0 means no empirical p-value is "
                                    "computed, not that it is 1."),
    "lincs_gmt": _V("tuple", optional=True),
    "expression_file": _V("path", optional=True, must_exist=True),
    "promoter_fasta": _V("path", optional=True),
    "has_ptc": _V("bool", optional=True,
                  why="None means unknown, which is different from False. "
                      "Whether the transcript carries a premature stop "
                      "decides if raising it helps or makes more truncated "
                      "protein."),
    "gtf": _V("path", optional=True),
    "genome": _V("path", optional=True),
    "use_genetics": _FLAG,
    "upstream_network": _FLAG,
    "upstream_depth": _V("int", minimum=1,
                         why="how many edges upstream to walk."),
    "approved_only": _FLAG,
    "use_dgidb": _FLAG,
    "use_opentargets": _FLAG,
    "find_analogs": _FLAG,
    "demo": _FLAG,
    "offline": _FLAG,
    "cache_dir": _V("str", optional=True),
    "min_ic": _V("float", minimum=0,
                 why="information content floor for a motif."),
}


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
    # locals() here is exactly the arguments, before anything else is
    # bound. Every one is checked against RUN_SPECS, so an annotation
    # stops being a comment.
    _validation.check_all(RUN_SPECS, dict(locals()))
    gene = gene.strip()
    # The one rule that spans two arguments, so it has no place in a
    # per-argument table.
    if upstream + downstream == 0:
        raise ValueError(
            "upstream and downstream must define a non-empty window")
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
    # Say which KIND of regulon arrived. A ChIP table and a curated one
    # support different claims, and the ledger is the line that is
    # supposed to make that visible.
    _reg_kind, _reg_phrase = (regulon_mod.describe_source(regulon_file)
                              if regulon_file and regulons
                              else ("none", ""))
    ledger.record("regulon", requested=bool(regulon_file),
                  ok=bool(regulons), yielded=len(regulons or {}),
                  detail=(f"{len(regulons or {})} {_reg_phrase}"
                          if regulons else
                          "no regulon file; coarse bundled tiers in use"))
    # Established from this run's arguments, every run. The table is
    # module state, so without the clear a second analysis in the same
    # process inherits the first one's and reports it as its own.
    dlv_mod.reset_expression()
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


# Report formatting lives in report.py. Re-exported because these two
# are public and callers import them from here.
from .report import (                                   # noqa: F401
    markdown_report, summary_report, write_outputs,
    corroboration, corroboration_table, unmapped_sources)
