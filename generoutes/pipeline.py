"""Orchestration: run every route for one gene and rank what they return."""

from __future__ import annotations

import json
import logging
import math
from dataclasses import dataclass, field
from pathlib import Path

from . import crosscheck
from . import (annotate, catalogue as cat_mod, cellcycle, congeners as cong_mod,
               delivery as dlv_mod, fda as fda_mod, knockout as knockout_mod,
               liabilities as liab_mod,
               sensitivity as sens_mod,
               promoterdna as dna_mod,
               tfactivity as tfa_mod,
               regulators as reg_mod,
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
from . import progress as _progress
from . import magnitude as mag_mod
from . import constraint as con_mod
from .http import Client

from . import fame as fame_mod
from . import isoforms
from . import poisonexon
from . import proteomics as prot_mod
from . import utr as utr_mod
from . import proliferation as prolif_mod
from . import standard_of_care as soc_mod

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
    # "bound", "not bound" or "not assayed". The last two look identical
    # in a peak count of zero and mean opposite things, so the state is
    # kept rather than re-derived from the count.
    chip_state: str = "not assayed"
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
            "chip_state": self.chip_state,
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


def _activity_regulons(regulon_file) -> tuple:
    """`(signed regulons, the file they came from)` for `--tf-activity`.

    The activity test asks whether a factor's targets moved the way the
    factor's sign says they should. That needs a table that signs its edges.

    A ChIP occupancy table (ReMap) says a factor was found near a gene and
    never whether it raises or lowers it, so every edge is unsigned and the
    test can only ask whether the targets moved at all. The launcher passes
    one as `--regulon-file` for the regulon breadth count, which used to make
    it the activity table as well and leave a signed table (CollecTRI, which
    signs every edge) sitting unused in data/. So when the file given is
    occupancy and a signed table can be found, the signed one is used for
    activity and the file given is still used for breadth. A curated or
    CollecTRI file given on purpose is the user's choice and wins. With no
    file given, whatever is in data/ is used, CollecTRI first.
    """
    if regulon_file and regulon_mod.describe_source(regulon_file)[0] == "chip":
        auto, auto_path = _optional(
            "signed regulon autoload", tfa_mod.autoload_regulons, ({}, None))
        if auto:
            log.info("--tf-activity uses %s, which signs its edges, instead "
                     "of %s, which is unsigned occupancy; %s is still used "
                     "for the regulon breadth count",
                     auto_path, regulon_file, regulon_file)
            return auto, auto_path
        log.warning("--tf-activity was given %s, which is unsigned ChIP "
                    "occupancy, and no signed table was found, so it can "
                    "only test whether targets moved, not which way. Put "
                    "CollecTRI at %s. %s", regulon_file,
                    " or ".join(tfa_mod.DEFAULT_REGULON_PATHS),
                    tfa_mod.COLLECTRI_SOURCE)
    if regulon_file:
        parsed = _optional(
            "signed regulon parse",
            lambda: tfa_mod.load_signed_regulons(regulon_file), None)
        return parsed, regulon_file
    found, source = _optional(
        "signed regulon autoload", tfa_mod.autoload_regulons, ({}, None))
    if not found:
        log.warning(
            "--tf-activity found no regulon table. Pass --regulon-file, or "
            "put one at %s. %s",
            " or ".join(tfa_mod.DEFAULT_REGULON_PATHS),
            tfa_mod.COLLECTRI_SOURCE)
    return found, source


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


# What a missing peak is worth when the factor was never assayed: the
# same unverified discount a run with no ChIP table at all applies.
UNVERIFIED_CHIP = 0.35


def chip_state(res: TFResult, has_chip: bool, assayed=None) -> str:
    """`bound`, `not bound` or `not assayed` for one factor.

    `assayed` is the set of symbols the ChIP table covers anywhere in the
    genome, or None when that is unknown. Only a factor the table
    demonstrably covers can be said not to be bound; for the rest the
    table is silent, and silence is not absence.
    """
    if not has_chip:
        return "not assayed"
    if res.chip_peaks:
        return "bound"
    if assayed and res.symbol.upper() in assayed:
        return "not bound"
    return "not assayed"


def score_tf(res: TFResult, has_chip: bool, goal: str = "up",
             assayed=None) -> float:
    b = res.best_hit
    if b is None:
        return 0.0
    w = WEIGHTS
    s = w["motif_quality"] * b.rel_score
    s += w["site_count"] * min(len(res.hits), 3) / 3
    s += w["proximity"] * _proximity(b.tss_offset)

    role = res.annotation.get("role", "unknown")
    s += w["activator_role"] * ROLE_CREDIT[goal].get(role, 0.0)

    state = chip_state(res, has_chip, assayed)
    if state == "bound":
        s += w["chip_evidence"]
    elif state == "not bound":
        pass                      # covered by the table and not bound here
    else:
        s += w["chip_evidence"] * UNVERIFIED_CHIP * min(b.rel_score, 1.0)

    if res.drugs:
        s += w["druggable"]
    return s


# The window a no-signal scan is retried at (bases before and after the start).
WIDE_UPSTREAM = 5000
WIDE_DOWNSTREAM = 1000

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
    "tissue": _V("str", optional=True, choices=("cns", "skin"),
                 why="cns down-weights drugs unlikely to reach the "
                     "compartment; skin keeps drugs approved only for "
                     "topical use in the reports."),
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
    "knockout_gmt": _V("tuple", optional=True),
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
    "widen_no_signal": _FLAG,
    "ctd_route": _FLAG,
    "null_shuffles": _V("int", minimum=0,
                        why="shuffled copies of each window used to estimate "
                            "how many motif sites are chance; 0 uses the "
                            "formula instead."),
    "upstream_factors": _V("int", minimum=1,
                           why="how many ranked factors the upstream search "
                               "covers."),
    "upstream_drugs_per_node": _V("int", minimum=1,
                                  why="how many drugs are listed for each "
                                      "upstream node."),
    "regulator_route": _FLAG,
    "regulator_measured": _FLAG,
    "regulator_steps": _V("int", minimum=1, maximum=2,
                          why="1 lists drugs at the gene's own recorded "
                              "regulators; 2 also lists drugs at the factors "
                              "that regulate those."),
    "approved_only": _FLAG,
    "fda_only": _FLAG,
    "include_local_only": _FLAG,
    "use_dgidb": _FLAG,
    "use_opentargets": _FLAG,
    "find_analogs": _FLAG,
    "demo": _FLAG,
    "offline": _FLAG,
    "cache_dir": _V("str", optional=True),
    "max_requests": _V("int", optional=True, minimum=1),
    "min_ic": _V("float", minimum=0,
                 why="information content floor for a motif."),
    "scenario": _V("str", optional=True,
                   choices=("haploinsufficiency", "hypomorph",
                            "biallelic_hypomorph",
                            "toxic_gain", "null", "silenced_allele",
                            "modifier"),
                   why="the shape of the dosage problem. It fixes what "
                       "counts as enough, and for two of them says that "
                       "changing transcription cannot help at all. "
                       "Declared rather than inferred: nothing here can "
                       "read a variant."),
    "broken_allele_makes": _V("str",
                              choices=("null", "truncated", "unknown"),
                              why="what the broken copy produces. "
                                  "Raising transcription raises both "
                                  "copies, so a truncated product is "
                                  "amplified alongside the working one."),
    "target_fraction": _V("float", minimum=0, maximum=1,
                          exclusive_min=True,
                          why="for a lowering plan, the share of the "
                              "current level you want left. 0.5 is a "
                              "halving."),
    "condition": _V("str", optional=True,
                    why="the condition being treated, for a gene with no "
                        "mutation in it."),
    "avoid_genotoxic": _V("bool", optional=True,
                          why="the patient already has chromosomal "
                              "instability or is cancer-predisposed, so a "
                              "compound that works by damaging DNA or the "
                              "spindle causes the defect being treated."),
}


# The arguments that name a file. Recorded in the result's provenance
# with size, modification time and - when the file is small enough to be
# worth hashing - its sha256, so a saved result says which bytes it read
# and not only which path it was pointed at.
RUN_FILE_ARGS = tuple(
    name for name, spec in RUN_SPECS.items() if spec.kind == "path")


def query_block(args: dict, **derived) -> dict:
    """What the run was asked to do, and what it worked out from that.

    Built from `RUN_SPECS` rather than written out by hand. The
    hand-written version recorded twenty-four of fifty-eight arguments,
    so a result could not be audited against the command that produced
    it - and when the panel turned out to be dropping arguments on the
    way in, the saved result could not show it, because the dropped ones
    were not among those recorded.

    `derived` is for values the run computed rather than received, like
    the number of motifs actually loaded. They are not arguments and are
    kept distinguishable by being passed separately.

    Paths are recorded as given, not resolved: the record should say
    what the run was pointed at, not where that happened to be today.
    """
    out: dict = {}
    for name in RUN_SPECS:
        value = args.get(name)
        # Path objects are not JSON-serialisable, and a run started from
        # Python may well pass one.
        out[name] = str(value) if isinstance(value, Path) else value
    out.update(derived)
    return out


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
        activity_permutations: int = 0, lincs_gmt: tuple | None = None,
        knockout_gmt: tuple | None = None, expression_file: str | None = None,
        promoter_fasta: str | None = None,
        has_ptc: bool | None = None,
        gtf: str | None = None,
        genome: str | None = None,
        use_genetics: bool = False,
        upstream_network: bool = False,
        upstream_depth: int = 3,
        null_shuffles: int = motifs_mod.NULL_SHUFFLES,
        widen_no_signal: bool = True,
        ctd_route: bool = True,
        upstream_factors: int = net_mod.UPSTREAM_FACTORS,
        upstream_drugs_per_node: int = net_mod.UPSTREAM_DRUGS_PER_NODE,
        regulator_route: bool = False,
        regulator_steps: int = 1,
        regulator_measured: bool = False,
        approved_only: bool = True, fda_only: bool = True,
        include_local_only: bool = False, use_dgidb: bool = True,
        use_opentargets: bool = True,
        find_analogs: bool = False,
        demo: bool = False, offline: bool = False,
        cache_dir: str | None = None, min_ic: float = 6.0,
        max_requests: int | None = None,
        scenario: str | None = None,
        broken_allele_makes: str = "unknown",
        target_fraction: float = 0.5,
        condition: str | None = None,
        avoid_genotoxic: bool = False) -> dict:
    """Run the full pipeline and return a JSON-serialisable result dict."""
    # locals() here is exactly the arguments, before anything else is
    # bound. Every one is checked against RUN_SPECS, so an annotation
    # stops being a comment, and the same snapshot is what the result's
    # query block records - so what a run was asked to do cannot drift
    # from what it reports it was asked to do.
    _args = dict(locals())
    _validation.check_all(RUN_SPECS, _args)
    gene = gene.strip()
    # `gene` is the only argument reassigned below, so the snapshot
    # needs refreshing once and then matches throughout.
    _args["gene"] = gene
    # The one rule that spans two arguments, so it has no place in a
    # per-argument table.
    if upstream + downstream == 0:
        raise ValueError(
            "upstream and downstream must define a non-empty window")
    # Likewise a scenario and a direction that point opposite ways. Both are
    # declared, neither inferred, so only a check can notice.
    _conflict = mag_mod.goal_conflict(scenario, goal)
    if _conflict:
        raise ValueError(_conflict)
    _ckw = {"offline": offline or demo, "max_requests": max_requests}
    client = Client(cache_dir=cache_dir, **_ckw) if cache_dir \
        else Client(**_ckw)
    _progress.start(f"{gene}: {'raising' if goal == 'up' else 'lowering'}"
                    + (", demo data" if demo else "")
                    + (", offline" if offline and not demo else ""))

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

    signed_regulons, regulon_source = (
        _activity_regulons(regulon_file) if tf_activity else (None, None))
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
    _uniprot_failed = False
    if include_direct:
        if not demo:
            try:
                ann = annotate.annotate_tfs(client, [gene])
                uniprot = (ann.get(gene.upper()) or {}).get("uniprot")
            except Exception as exc:                       # noqa: BLE001
                # This was `pass`. The search below still runs without the
                # accession, but a lookup that failed and a gene with no
                # accession look identical afterwards, and the progress
                # line then reads "none annotated" as a finding about the
                # gene. Say what happened.
                log.warning("UniProt lookup for %s failed (%s: %s); the "
                            "direct-drug search continues without an "
                            "accession, so an empty result here is weaker "
                            "evidence than usual", gene,
                            type(exc).__name__, str(exc)[:100])
                _uniprot_failed = True
        direct_hits = direct.find_direct(
            client, gene, uniprot, goal=goal,
            approved_only=approved_only, use_dgidb=use_dgidb, demo=demo)

    _progress.step(0, "drugs acting on the gene itself",
                   f"{len(direct_hits)} found" if direct_hits
                   else ("none annotated (the UniProt lookup failed, so "
                         "this is not a finding)" if _uniprot_failed
                         else "none annotated"))

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

    _progress.step(1, "promoter sequence",
                   f"{len(promoters)} transcript(s)")

    # ------------------------------------------------------------ 2. motifs
    motif_set = (demo_data.demo_motifs() if demo
                 else motifs_mod.load_motifs(client, pfm_file, matrix_ids))
    log.info("scanning with %d motifs at p <= %g", len(motif_set), pvalue)
    _progress.step(2, "motifs loaded",
                   f"{len(motif_set)} at p <= {pvalue:g}; scanning")

    def _scan_window(prom):
        """`(hits, quality, scan_info)` for one promoter window."""
        scan_info: dict = {}
        hits = motifs_mod.scan_promoter(
            prom.sequence, motif_set, pvalue=pvalue, min_ic=min_ic,
            offset_fn=prom.offset_to_tss, skip_repeats=mask_repeats and not demo,
            info=scan_info)
        hits = motifs_mod.collapse_overlapping(hits)
        scannable = motifs_mod.scannable_length(
            prom.sequence, skip_repeats=mask_repeats and not demo)
        # How many sites chance gives is measured on shuffled copies of this
        # window, since the formula that assumes every motif matches at its
        # stated p-value undercounted on a masked promoter. 0 shuffles, or a
        # window with too little scannable sequence, uses the formula.
        null = None
        if null_shuffles >= 10 and scannable >= 100:
            null = motifs_mod.empirical_null(
                prom.sequence, motif_set, pvalue=pvalue, min_ic=min_ic,
                skip_repeats=mask_repeats and not demo,
                n_shuffles=null_shuffles)
        quality = motifs_mod.scan_quality(
            len(hits), scan_info.get("n_scanned", len(motif_set)), scannable,
            pvalue, null=null,
            observed_scores=[h.pvalue for h in hits
                             if h.pvalue is not None])
        quality["motifs_loaded"] = len(motif_set)
        quality["motifs_scanned"] = scan_info.get("n_scanned")
        quality["motifs_unreachable"] = scan_info.get("n_unreachable")
        return hits, quality, scan_info

    outputs = []

    def _analyse_promoter(prom) -> dict:
        """Everything this run has to say about one promoter window.

        A function rather than a loop body so its names are its own.
        It reads the run-level values above it - the client, the
        ledger, the shared screen cache - and rebinds none of them.
        """
        hits, quality, scan_info = _scan_window(prom)
        # A window with no signal says nothing about the gene, and real sites
        # are often further out. Scanning 5 genes at four window sizes found
        # 20 to 48 sites above chance at -5000/+1000 where the core window gave
        # none (ABCD2, CYP24A1), though the chance share did not fall. So when
        # the requested window shows no signal, the same transcript is
        # rescanned wider, and the wider window is kept only if it does show
        # some. What was tried is recorded either way.
        quality["window_widened"] = None
        if (widen_no_signal and not demo
                and quality["verdict"] == "no signal"
                and (prom.upstream < WIDE_UPSTREAM
                     or prom.downstream < WIDE_DOWNSTREAM)):
            _wide = _optional(
                "wider promoter window",
                lambda: genes.get_promoters(
                    client, gene, species=species,
                    upstream=max(prom.upstream, WIDE_UPSTREAM),
                    downstream=max(prom.downstream, WIDE_DOWNSTREAM),
                    all_tss=all_tss, mask="soft" if mask_repeats else None),
                None)
            _match = next((w for w in _wide or []
                           if w.transcript_id == prom.transcript_id), None)
            if _match is not None:
                _h2, _q2, _i2 = _scan_window(_match)
                _try = {"from": [prom.upstream, prom.downstream],
                        "to": [_match.upstream, _match.downstream],
                        "sites_before": quality["observed"],
                        "sites_after": _q2["observed"],
                        "chance_after": _q2["expected_by_chance"],
                        "verdict_after": _q2["verdict"]}
                if _q2["verdict"] != "no signal":
                    _try["kept"] = True
                    prom, hits, quality, scan_info = _match, _h2, _q2, _i2
                    quality["window_widened"] = _try
                else:
                    _try["kept"] = False
                    quality["window_widened"] = _try
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

        _progress.step(3, "factors annotated",
                       f"{len(symbols)} with a site in this window")

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

        _progress.step(
            4, "measured occupancy",
            f"peaks for {sum(1 for r in by_tf.values() if r.chip_peaks)} "
            f"of {len(by_tf)} factor(s)" if chip_status.get("ok")
            else "none; factor ranking rests on sequence alone")

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

        # Which factors the ChIP table covers at all, from the regulon
        # built out of the same BED genome-wide. Without one, coverage is
        # unknown and a missing peak goes back to meaning nothing.
        _assayed = (set(regulons) if regulons and _reg_kind == "chip"
                    else None)
        for r in by_tf.values():
            r.chip_state = chip_state(r, chip_src.usable, _assayed)
            r.score = score_tf(r, has_chip=chip_src.usable, goal=goal,
                               assayed=_assayed)

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

        _progress.step(5, "factors worth pursuing",
                       f"{len(actionable)} of {len(by_tf)} have a role or "
                       f"are nuclear receptors")

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
                # "ok" used to mean only that nothing crashed. A model that
                # trained for no receptor looked the same as one that worked,
                # so say how many grades actually got a call.
                n_called = sum(
                    1 for g in grades
                    if getattr(getattr(g, "_prediction", None), "predicted",
                               "undetermined") != "undetermined")
                ledger.record(
                    "direction_model", requested=True,
                    ok=broke == 0 and bool(grades),
                    yielded=n_called,
                    detail=(f"{broke} of {len(grades)} grade(s) could not be "
                            f"direction-checked" if broke else
                            f"{n_called} of {len(grades)} grade(s) got a "
                            f"model call" + ("" if n_called else
                                             "; too few labelled agonists and "
                                             "antagonists to train")))

            # Collapse salt forms. One pharmacology was occupying five rows
            # (Liothyronine / Liothyronine Sodium / Levothyroxine Sodium ...),
            # pushing distinct chemistry off the bottom of every table and
            # inflating the candidate count the Pareto accounting reports.
            # Done here so every downstream consumer — catalogue, outcome
            # chain, discovery — inherits one row per active moiety.
            grades, collapsed = cong_mod.collapse_grades(grades)
            for g in grades:
                alts = collapsed.get((g.drug_name, g.tf_symbol))
                if alts:
                    g.notes.append(
                        "same active moiety as " + ", ".join(alts[:4])
                        + (f" (+{len(alts) - 4} more)" if len(alts) > 4 else ""))

        _progress.step(6, "drug/factor pairs graded",
                       f"{len(grades)} pair(s)" if grade_drugs
                       else "skipped, --no-grade")

        # ---- 6b. does the drug reach the cells expressing the gene?
        # Runs HERE, straight after grading, because it is the last
        # thing that changes a factor's score. Further down - where it
        # used to be - the stabilisation, ligand-axis, mechanism and
        # upstream sections had already taken their `ranked[:N]` from
        # the unpenalised order, so one report showed two different
        # top lists and a factor absent from the gene's tissue stayed
        # in half of them.
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


        # ------------------------------ 7. propagate to a transcriptional call
        estimates = []
        if estimate_outcome:
            by_symbol = {r.symbol: r for r in ranked}
            for g in grades:
                est = outcome_mod.estimate(
                    g.drug_name, g.tf_symbol, prom.gene_symbol, goal=goal,
                    grade=g, tf_result=by_symbol.get(g.tf_symbol),
                    prediction=getattr(g, "_prediction", None))
                d = est.to_dict()
                d["verdict"] = est.verdict
                d["next_experiment"] = outcome_mod.narrowing_experiment(est)
                estimates.append(d)

            # Order by the chain's own result. Inheriting grade order meant the
            # discovery winner could fall outside the printed rows while worse
            # candidates were shown — two rankings disagreeing silently.
            estimates.sort(key=lambda e: -e["p_point"])

        # ---- does the gene change when the factor is lost?
        # Read from perturbation experiments (see knockout.py). It replaces a
        # fixed placeholder for "responsiveness"; with no files on disk the
        # estimates carry no `knockout` entry and the report says so.
        knockout_results: dict = {}
        if not demo:
            _ko_paths = knockout_mod.resolve_paths(knockout_gmt)
            if _ko_paths:
                knockout_results = _optional(
                    "knockout evidence",
                    lambda: knockout_mod.assess_many(
                        prom.gene_symbol, [r.symbol for r in ranked],
                        _ko_paths), {})
                for d in estimates:
                    d["knockout"] = knockout_results.get(
                        d["tf_symbol"].upper())
                n_meas = sum(1 for v in knockout_results.values()
                             if v["class"] == knockout_mod.MEASURED)
                n_test = sum(1 for v in knockout_results.values()
                             if v["n_experiments"])
                ledger.record(
                    "knockout", requested=True, ok=bool(knockout_results),
                    yielded=n_meas,
                    detail=(f"{n_test} of {len(knockout_results)} factor(s) "
                            f"were perturbed in the knockout data; the gene "
                            f"moved for {n_meas}"
                            if knockout_results else
                            "the knockout files could not be read"))

        _progress.step(7, "transcriptional call",
                       f"{len(estimates)} route(s) estimated"
                       if estimate_outcome else "skipped, --no-outcome")

        # ------------------------------------- 8. optional library screen
        screens = {}
        # The receptors this screen could be asked about: factors with a drug
        # that names its ChEMBL target. Counted whether or not the screen
        # runs, so the ledger can tell "nothing to screen" from "could not".
        _screen_eligible = {getattr(d, "target_chembl_id", None)
                            for r in ranked for d in r.drugs} - {None}
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
                # ok whether or not anything was found: this is a local scan
                # of a sequence in hand, so it either ran or it did not, and
                # "no feature beat the shuffled null" is its answer. Recording
                # that as a failure kept runs at "degraded" for a finished,
                # negative analysis.
                ledger.record(
                    "promoter_dna", requested=True, ok=True,
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

        # One batched read of what every approved drug targets, built the
        # first time a stage needs it. None when no local drug table is on
        # disk, only approved drugs were not asked for, or the build failed.
        _index_box: list = []

        def _approved_index():
            if not _index_box:
                idx = None
                if approved_only:
                    try:
                        from . import localdrugs as _ld
                        _lib = _ld.library()
                        if _lib:
                            idx = drugs_mod.approved_target_index(
                                client, _lib)
                    except Exception as exc:                  # noqa: BLE001
                        log.warning("could not build the approved-drug index "
                                    "(%s); asking about upstream nodes one "
                                    "at a time, capped at %d per factor",
                                    str(exc)[:100],
                                    net_mod.UPSTREAM_MAX_CHECKED)
                _index_box.append(idx)
            return _index_box[0]

        # ------------------- 8c. multi-hop upstream routes to undruggable TFs
        upstream_hits = {}
        upstream_search = {}
        if upstream_network and not demo:
            edges = net_mod.fetch_network(client)
            # The high-scoring factors with no ligand are precisely what this
            # search exists for, so they go in first rather than depending on
            # where they happen to fall in the ranking.
            _priority = [r for r in ranked if r.symbol in
                         set(attrition.get("top25_without_drug") or [])]
            # Two fewer than the factors searched, which is the six it always
            # was at the default of eight.
            _rest = [r for r in ranked[:max(6, upstream_factors - 2)]
                     if r not in _priority]
            _up = {"factors_searched": 0, "factors_errored": 0,
                   "nodes_reachable": 0, "nodes_checked": 0,
                   "lookups_failed": 0, "factors_truncated": 0,
                   "max_checked": net_mod.UPSTREAM_MAX_CHECKED}
            _drug_cache: dict = {}
            # One batched read of what every approved drug targets, when the
            # local drug table is on disk, so each upstream node is a lookup
            # and none has to be skipped. Without it, or if it fails, each
            # node is asked about on its own and the walk is capped.
            _drug_index = _approved_index() if edges else None
            _up["drug_index"] = _drug_index is not None
            for r in (_priority + _rest)[:upstream_factors] if edges else []:
                # The point of this step is TFs with no direct pharmacology.
                if r.drugs:
                    continue
                try:
                    srch = net_mod.find_druggable_upstream(
                        client, r.symbol, goal=goal,
                        tf_role=r.annotation.get("role", "activator"),
                        max_depth=upstream_depth, edges=edges,
                        approved_only=approved_only, drug_cache=_drug_cache,
                        drug_index=_drug_index,
                        drugs_per_node=upstream_drugs_per_node)
                except Exception as exc:                      # noqa: BLE001
                    log.warning("upstream search failed for %s: %s", r.symbol, exc)
                    _up["factors_errored"] += 1
                    continue
                _up["factors_searched"] += 1
                _up["nodes_reachable"] += srch.n_reachable
                _up["nodes_checked"] += srch.n_checked
                _up["lookups_failed"] += srch.n_lookup_failed
                _up["factors_truncated"] += int(srch.truncated)
                if srch.stopped:
                    log.warning("upstream search for %s stopped early: %s",
                                r.symbol, srch.stopped)
                if srch.hits:
                    upstream_hits[r.symbol] = [h.to_dict() for h in srch.hits[:10]]
            upstream_search = _up
            _scope = (f"{_up['nodes_checked']} of {_up['nodes_reachable']} "
                      f"reachable upstream node(s) checked across "
                      f"{_up['factors_searched']} factor(s)")
            if not edges:
                ledger.record("upstream", requested=True, ok=False,
                              yielded=0,
                              detail="the OmniPath network could not be "
                                     "loaded, so no upstream search was "
                                     "possible")
            elif not _up["factors_searched"] and not _up["factors_errored"]:
                ledger.record("upstream", requested=True, ok=True, yielded=0,
                              detail="every leading factor already has a "
                                     "drug acting on it directly, so no "
                                     "upstream search was needed")
            elif upstream_hits:
                ledger.record("upstream", requested=True, ok=True,
                              yielded=len(upstream_hits),
                              detail=f"druggable routes for "
                                     f"{len(upstream_hits)} factor(s); "
                                     f"{_scope}")
            elif (not _up["factors_searched"]
                  or _up["lookups_failed"] * 2 > _up["nodes_checked"]):
                ledger.record(
                    "upstream", requested=True, ok=False, yielded=0,
                    detail=f"{_up['lookups_failed']} of "
                           f"{_up['nodes_checked']} drug lookup(s) failed, so "
                           f"this says nothing about whether a druggable "
                           f"route exists; {_scope}")
            else:
                _more = (f"; {_up['factors_truncated']} factor(s) reached "
                         f"the limit of {_up['max_checked']} nodes with more "
                         f"still reachable, so a route further out is not "
                         f"ruled out" if _up["factors_truncated"] else "")
                ledger.record(
                    "upstream", requested=True, ok=True, yielded=0,
                    detail=f"no node checked has an approved drug that moves "
                           f"it the needed way; {_scope}{_more}")

        # --------------- 8c2. drugs at the gene's known regulators
        # The reverse of the promoter scan: take the factors CollecTRI and
        # the knockout data record as regulating this gene, skip those that
        # already have a drug route, and ask what approved drug acts on each.
        regulator_hits: dict = {}
        regulators_undrugged: list = []
        if regulator_route and not demo:
            _reg_tables = signed_regulons
            if not _reg_tables:
                _reg_tables, _ = _optional(
                    "signed regulon autoload", tfa_mod.autoload_regulons,
                    ({}, None))
            _reg_index = _approved_index()
            if not _reg_tables:
                ledger.record(
                    "regulators", requested=True, ok=False, yielded=0,
                    detail="no signed regulon table (CollecTRI) was found, "
                           "so the known regulators of the gene could not "
                           "be listed")
            elif _reg_index is None:
                ledger.record(
                    "regulators", requested=True, ok=False, yielded=0,
                    detail=("the regulator route reads the approved-drug "
                            "index, which is not built when investigational "
                            "drugs are included" if not approved_only else
                            "the approved-drug index could not be built "
                            "(it needs the local drug table), so no drug "
                            "could be matched to a regulator"))
            else:
                _regs = _optional(
                    "known regulators",
                    lambda: reg_mod.candidates(
                        prom.gene_symbol, _reg_tables,
                        knockout_paths=knockout_gmt,
                        steps=regulator_steps,
                        skip=[r.symbol for r in ranked if r.drugs]), None)
                if _regs is None:
                    ledger.record("regulators", requested=True, ok=False,
                                  yielded=0,
                                  detail="the regulator lookup failed")
                else:
                    _binders, _binder_failed = None, []
                    if regulator_measured:
                        # Asked only for a regulator the mechanism index has
                        # nothing for, one lookup each, cached by the client.
                        from . import localdrugs as _ld2
                        _blib, _bcache = _ld2.library(), {}

                        def _binders(sym):
                            if sym not in _bcache:
                                try:
                                    _bcache[sym] = drugs_mod.measured_binders(
                                        client, sym, _blib)
                                except Exception as exc:      # noqa: BLE001
                                    log.debug("measured binders for %s: %s",
                                              sym, exc)
                                    _binder_failed.append(sym)
                                    _bcache[sym] = []
                            return _bcache[sym]
                    regulator_hits = reg_mod.find_drugs(
                        _regs, goal, _reg_index, binders=_binders)
                    # The other half of the same walk. These produce no
                    # route - there is no molecule to give - but they are
                    # what is actually recorded as running the gene, and
                    # dropping them left the section listing whatever
                    # happened to be druggable.
                    regulators_undrugged = reg_mod.without_drugs(
                        _regs, goal, _reg_index, binders=_binders)
                    _n_drugs = sum(len(v["drugs"])
                                   for v in regulator_hits.values())
                    ledger.record(
                        "regulators", requested=True, ok=True,
                        yielded=len(regulator_hits),
                        detail=(f"{len(_regs)} recorded regulator(s) of the "
                                f"gene had no drug route yet; "
                                f"{len(regulator_hits)} have an approved "
                                f"drug that moves them the needed way "
                                f"({_n_drugs} route(s)); "
                                f"{len(regulators_undrugged)} have none, and "
                                f"are listed rather than dropped"
                                + (f"; {len(_binder_failed)} measured-activity "
                                   f"lookup(s) failed, which is not a finding "
                                   f"that none exist" if _binder_failed
                                   else "")))

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
        # Demo mode uses synthetic factors and must not score them against a
        # real DepMap release that happens to be sitting in data/. A file
        # named with --depmap-file is still read.
        _dm_path = (None if demo and not depmap_file
                    else depmap_mod.resolve(depmap_file))
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
            # NOT `by_tf`: that name already holds this window's
            # factors, and reusing it here made the reported `n_tfs`
            # the size of the regulon table rather than the number of
            # factors with a site in this promoter.
            activity_by_tf = {a.tf.upper(): a for a in acts}

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
                a = activity_by_tf.get(r.symbol.upper())
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

        # ---- the transcript's own UTRs: translation and decay
        # The coverage audit marked these two steps "listed only", and
        # both can be asked of sequence this package already fetches. It
        # annotates a lever; it does not search for a drug that pulls it.
        utr_scan: dict = {}
        if not demo:
            _five, _three = _optional(
                "transcript UTRs",
                lambda: genes.fetch_utrs(client, prom.transcript_id),
                (None, None))
            if _five or _three:
                utr_scan = _optional(
                    "UTR feature scan",
                    lambda: utr_mod.scan(_five or "", _three or ""), {}) or {}
            _lever = sum(
                1 for k in ("translation", "stability")
                for v in [utr_scan.get(k) or {}]
                if (v.get("enriched")
                    or any((v.get(sub) or {}).get("enriched")
                           for sub in ("pentamer", "nonamer",
                                       "g_quadruplex", "m6a_drach",
                                       "pumilio"))
                    # The TOP motif is positional, so it is present or
                    # absent rather than enriched, and it is as real a
                    # lever as any of the counted ones.
                    or (v.get("top_motif") or {}).get("present")))
            ledger.record(
                "utr", requested=True, ok=bool(utr_scan), yielded=_lever,
                detail=(f"5'UTR {utr_scan.get('five_prime_length', 0)} bp, "
                        f"3'UTR {utr_scan.get('three_prime_length', 0)} bp; "
                        f"{_lever} feature(s) beat a shuffle of the same "
                        f"sequence" if utr_scan else
                        "the transcript's UTRs could not be derived, so "
                        "translation and mRNA stability were not examined"))

        # ---- which drugs change how much of the PROTEIN exists?
        # Asked of the gene's own product, which is what every route
        # above is a proxy for. Coverage is the limit, not the question:
        # public drug-proteomics reaches a small fraction of the approved
        # library, so an absent compound is unexamined and not negative.
        protein_abundance: list = []
        if not demo:
            _prot = _optional(
                "protein abundance table",
                lambda: prot_mod.load(prom.gene_symbol,
                                      want="up" if goal == "up" else "down"),
                [])
            protein_abundance = [h.to_dict() for h in (_prot or [])]
            ledger.record(
                # Not on disk means not asked, not failed.
                "protein_abundance", requested=prot_mod.available(),
                ok=prot_mod.available(),
                yielded=len(protein_abundance),
                detail=(f"{len(protein_abundance)} compound(s) measured to "
                        f"move {prom.gene_symbol} protein the needed way"
                        if prot_mod.available() else prot_mod.SOURCE))

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
            # The recorded regulators that no drug reaches are exactly what
            # this step exists for, and most of them never appear in the
            # motif ranking, so they were never asked. Direct ones only: a
            # two-step regulator's sign is a product of two edges, and the
            # list runs to hundreds. The sign gives the role - a factor
            # that RAISES the gene has to go up to raise it.
            _asked = {(r.get("symbol") or "").upper() for r in rows}
            _from_regulators = 0
            for _reg in regulators_undrugged:
                _sym = (_reg.get("symbol") or "").strip()
                if _reg.get("via") or not _sym or _sym.upper() in _asked:
                    continue
                _asked.add(_sym.upper())
                _from_regulators += 1
                rows.append({
                    "symbol": _sym,
                    "role": ("activator" if _reg.get("sign") == 1
                             else "repressor")})
            found = _optional(
                "factor abundance scan",
                lambda: abund_mod.find(rows, goal=goal), [])
            factor_abundance = [h.to_dict() for h in (found or [])]
            undirected_factors = abund_mod.undirected_factors(rows, goal=goal)
            _scope = (f"{len(rows)} factor(s) asked, {_from_regulators} of "
                      f"them recorded regulators the motif scan did not "
                      f"rank")
            ledger.record(
                "factor_abundance", requested=True, ok=True,
                yielded=len(factor_abundance),
                detail=(f"{len(factor_abundance)} compound(s) move a "
                        f"factor's own transcript the needed way; {_scope}; "
                        f"{len(undirected_factors)} factor(s) skipped for "
                        f"an unrecorded role"
                        if factor_abundance else
                        f"no compound moved a factor's own transcript the "
                        f"needed way in enough contexts to report "
                        f"({_scope})"))

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
                          detail=_genetics_detail(
                              _phen, _local_gen,
                              error=(genetic["cis"] or {}).get("error")))
            genetic["trans"] = [
                gen_mod.tf_link_support(client, r.symbol, prom.gene_symbol).to_dict()
                for r in ranked[:3]]

            # Actually run the trans-eQTL lookup the advice above describes.
            # The eQTLGen significant file is 2 MB; telling someone to go
            # check it by hand, once per factor, in identical wording, was
            # advice rather than an answer.
            def _trans_eqtl():
                # A table that cannot be read must not take the whole
                # section with it. `trans_support` reports an unreadable
                # table per factor as NOT answered, and the ledger reads
                # `checked`; letting this raise would skip both and leave
                # the report silent about a source it was asked for.
                any_hit = _optional(
                    "eQTLGen index",
                    lambda: eqtl_mod.gene_has_any_trans_eqtl(
                        prom.gene_symbol, client=client), False)
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

            # The approved-drug index turns "E3 ligase inhibition" from a
            # generic mechanism into a question with an answer: is there a
            # drug at the ligase that destroys THIS protein?
            beyond["product_stability"] = stab_mod.product_stabilisation(
                prom.gene_symbol, goal=goal,
                drug_index=None if demo else _optional(
                    'approved-drug index for the degradation route',
                    _approved_index, None))

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
            # Which endpoint, decided by what the caller DECLARED about
            # the alleles and never guessed. A promoter reporter reads
            # transcription from a cloned fragment; where the problem is
            # the protein, it cannot see a stabiliser, a splice modifier
            # or a readthrough agent, so it would screen for the wrong
            # thing with full confidence.
            _why = []
            if scenario == "biallelic_hypomorph":
                _why.append("the only working copy is itself impaired")
            elif scenario == "hypomorph":
                _why.append("both copies are impaired")
            if has_ptc:
                _why.append("the transcript carries a premature stop")
            if broken_allele_makes == "truncated":
                _why.append("the broken copy makes a truncated product")
            screen_plan = screens_mod.design(
                prom, [r.symbol for r in ranked[:20]],
                [g.to_dict() for g in grades],
                [t.to_dict() for t in direct_hits],
                protein_endpoint=bool(_why),
                why=", and ".join(_why)).to_dict()

        # -------------------------- 9. what beats the current standard of care
        _progress.step(8, "library screen",
                       "run" if (screen_receptors and not demo)
                       else "not requested")
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
                # Three outcomes, two of which used to share a
                # sentence. The wording lives in discovery so it can be
                # tested without a network round trip.
                assessed = improvements.get("ctd_assessed", True)
                ctd_detail = discovery.ctd_ledger_detail(
                    n_ctd, rl.get("wanted", 0), assessed)
                ledger.record(
                    "ctd", requested=True, ok=bool(n_ctd) or assessed,
                    yielded=n_ctd, detail=ctd_detail)

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

        # Approved HDAC, DNMT and EZH2 inhibitors beyond the curated list.
        _chromatin_extra: list = []
        if not demo:
            _idx = _optional("approved-drug index for chromatin agents",
                             _approved_index, None)
            if _idx:
                _chromatin_extra = cat_mod.chromatin_agents_from_index(_idx)

        # Where CTD and LINCS disagree, and which ranked drugs CTD records
        # moving the gene the other way. Changes no tier or score.
        _crosscheck = None
        if not demo and improvements:
            _crosscheck = _optional(
                "CTD and LINCS cross-check",
                lambda: crosscheck.build(
                    prom.gene_symbol, goal, improvements,
                    use_lincs=bool(use_lincs or lincs_file or lincs_gmt)),
                None)

        block_for_catalogue = {
            "crosscheck": _crosscheck,
            "chromatin_extra": _chromatin_extra,
            # What CTD records for this gene, for the literature route.
            "literature": ({
                "gene": prom.gene_symbol,
                "reported_links": improvements.get("reported_links"),
                "papers": improvements.get("ctd_papers"),
                "protein_links": improvements.get("reported_protein_links"),
                "protein_papers": improvements.get("ctd_protein_papers")}
                if (ctd_route and improvements) else None),
            "direct_therapies": [t.to_dict() for t in direct_hits],
            "drug_grades": [g.to_dict() for g in grades],
            "mechanisms": mechanisms,
            "stabilisation": stabilisation,
            "ligand_axes": ligand_axes,
            "upstream": upstream_hits,
            "upstream_search": upstream_search,
            "regulators": regulator_hits,
            "lincs": lincs,
            "screens": screens,
            # Both of these were computed above and left out of the dict
            # the catalogue actually reads, so their routes never reached
            # the tables even though their report sections rendered.
            "activity_drugs": activity_drugs,
            "factor_abundance": factor_abundance,
            "protein_abundance": protein_abundance,
            "undirected_factors": undirected_factors,
            "docking": list(docking_runs.values()),
            # The factors at this promoter and their roles. The docking
            # collector needs both: to tell a dock against one of them
            # from a dock against something else, and to say whether the
            # conformation docked is the one you want stabilised.
            # `annotation` is a dict subclass, so this is a key read.
            # The recorded regulators go in too, with their sign as the
            # role: most never appear in `ranked`, so without them a
            # route through MYC or FOXM1 has no recorded direction and
            # the cancer-risk check can make no claim about it.
            "tf_roles": {
                **{r["symbol"]: ("activator" if r.get("sign") == 1
                                 else "repressor")
                   for r in regulators_undrugged if not r.get("via")},
                **{r.symbol: (r.annotation or {}).get("role")
                   for r in ranked}},
            "sensitivity": sensitivity,
            "product_screen": product_screen,
            "promoter_dna": promoter_dna_scan,
            "gene_symbol": prom.gene_symbol,
            # Each factor's confidence, role and the way each of its drugs
            # moves it, so a ligand found by a screen is scored like a graded
            # one. Read by `catalogue.from_screen`.
            "tf_info": cat_mod.tf_info_from_tfs(
                [r.to_dict() for r in ranked]),
        }
        _mitotic_gene = (cellcycle.is_mitotic_gene(gene)
                         or cellcycle.is_mitotic_gene(prom.gene_symbol))
        _fdp = quality.get("false_discovery_proportion")
        all_routes = cat_mod.build(block_for_catalogue, goal=goal,
                                   mitotic_gene=_mitotic_gene,
                                   scan_fdp=_fdp,
                                   avoid_genotoxic=avoid_genotoxic)
        # Compounds that move a large share of the transcriptome lose score
        # within their tier. The routes are built once to learn which
        # compounds are candidates, and again only if one of them is that
        # broad, so a run with no LINCS files or no broad candidate pays
        # nothing and scores as it always did.
        # Did the compound move the gene, or move the cell cycle? Asked
        # only for a gene whose transcript tracks the cycle, because for
        # any other gene it is not a confound. Measured on the compound's
        # own signatures; see proliferation.py.
        _prolif: dict = {}
        if _mitotic_gene and not demo:
            _cands = [r.get("drug") for r in
                      (all_routes.get("routes") or []) if r.get("drug")]
            _prolif = _optional(
                "proliferation confound check",
                lambda: prolif_mod.assess(prom.gene_symbol, _cands,
                                          direction=goal), {}) or {}
            _n_tracks = sum(1 for v in _prolif.values()
                            if v.get("tracks_proliferation"))
            ledger.record(
                "proliferation",
                requested=sig_mod.autoload_gmt_available(goal),
                ok=bool(_prolif),
                yielded=_n_tracks,
                detail=(f"{len(_prolif)} candidate(s) checked against the "
                        f"mitotic programme; {_n_tracks} move it along with "
                        f"the gene" if _prolif else
                        "no LINCS .gmt on disk, so no candidate could be "
                        "checked for moving the cell cycle rather than this "
                        "gene"))

        # How many OTHER genes the literature records each candidate as
        # moving. The retrospective found this is the control that matters
        # most for reading the order: a baseline of drug popularity beat
        # every route here. See fame.py.
        _fame: dict = {}
        if not demo:
            _ftable = _optional("fame table", fame_mod.load, None)
            _names = [r.get("drug") for r in
                      (all_routes.get("routes") or []) if r.get("drug")]
            _fame = fame_mod.assess(_names, _ftable)
            _n_fame = sum(1 for v in _fame.values() if v.get("promiscuous"))
            ledger.record(
                "fame", requested=_ftable is not None,
                ok=bool(_fame), yielded=_n_fame,
                detail=(f"{len(_fame)} candidate(s) found in the fame "
                        f"table over {_ftable['n_genes']} gene(s); "
                        f"{_n_fame} are recorded for most of them"
                        if _fame and _ftable else
                        "no drug-fame table, so no candidate could be "
                        "checked for being a compound the literature "
                        "records for every gene. " + fame_mod.SOURCE))

        if not demo:
            _b = _optional("LINCS compound breadth",
                           lambda: _compound_breadth([{"all_routes": all_routes}]),
                           {})
            _very = {k: v for k, v in (_b.get("compounds") or {}).items()
                     if v.get("very_broad")}
            if (_very or any(v.get("tracks_proliferation")
                             for v in _prolif.values())
                    or any(v.get("promiscuous") for v in _fame.values())):
                all_routes = cat_mod.build(block_for_catalogue, goal=goal,
                                           breadth=_b.get("compounds"),
                                           mitotic_gene=_mitotic_gene,
                                           scan_fdp=_fdp,
                                           proliferation=_prolif,
                                           avoid_genotoxic=avoid_genotoxic,
                                           fame=_fame)

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
            _ok, _detail = _screen_outcome(screens, _screen_eligible)
            ledger.record("screen", requested=True, ok=_ok,
                          yielded=len(screens), detail=_detail)

        # Which of these factors a shuffle could not have produced. Done
        # here rather than in scan_quality because the measured peaks are
        # per factor and that function sees only sequence.
        if quality:
            quality["occupancy"] = motifs_mod.occupancy_partition(
                [r.to_dict() for r in ranked], quality)

        return {
            "promoter": prom.to_dict(),
            "sequence": prom.sequence,
            "n_sites": len(hits),
            "scan_quality": quality,
            # `collapsed` is keyed by (drug, factor), which json cannot
            # write, so the saved field joins the pair into one string.
            "salt_forms_collapsed": (
                {f"{drug} @ {tf}": alts
                 for (drug, tf), alts in collapsed.items()}
                if grade_drugs else {}),
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
            "knockout": knockout_results,
            "improvements": improvements,
            "crosscheck": _crosscheck,
            "all_routes": all_routes,
            "mechanisms": mechanisms,
            "upstream": upstream_hits,
            "upstream_search": upstream_search,
            "regulators": regulator_hits,
            "regulators_without_drugs": regulators_undrugged,
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
            "proliferation": _prolif,
            "protein_abundance": protein_abundance,
            "utr": utr_scan,
            "fame": _fame,
            "modality_design": modality,
            "screen_plan": screen_plan,
            "lincs": lincs,
            "genetics": genetic,
            "screens": screens,
            "hits": [vars(h) for h in hits],
        }

    for prom in promoters:
        outputs.append(_analyse_promoter(prom))

    # How much the --all-tss windows duplicate each other. CCR5's two
    # overlapped by 4275 of 6000 bp, so they were largely one analysis
    # counted twice — which inflates every per-gene route count downstream.
    tss_overlap = nbr_mod.window_overlap(
        [o["promoter"] for o in outputs]) if len(outputs) > 1 else []

    # The requirement, when a mechanism was declared, and what the
    # population says about this gene's tolerance of a dose change.
    _requirement = None
    if scenario:
        try:
            _ask = mag_mod.ask_for(
                scenario, gene, broken_allele_makes=broken_allele_makes,
                target_fraction=target_fraction, goal=goal,
                condition=condition or "")
            _routes = [r for o in outputs
                       for r in ((o.get("all_routes") or {}).get("routes")
                                 or [])]
            _effects = mag_mod.measured_effects(_routes, gene)
            _requirement = {
                **_ask.to_dict(),
                "effects": _effects,
                "summary": mag_mod.summarise(_ask, _effects),
            }
        except Exception as exc:                            # noqa: BLE001
            log.warning("requirement not computed: %s", str(exc)[:120])

    def _assess_dosage() -> dict:
        c = con_mod.assess(gene)
        return {**c.to_dict(),
                "summary": con_mod.summarise(
                    c, goal,
                    one_allele_lost=scenario == "haploinsufficiency",
                    broken_allele_makes=broken_allele_makes)}

    _dosage = _optional("dosage constraint", _assess_dosage, None)

    # The declared scenario against the measured constraint. Both are now
    # in hand, which is the earliest this can be asked: `goal_conflict`
    # runs before the analysis because it needs only two arguments, and
    # this needs a table that is read during it.
    if _requirement is not None and _dosage:
        _clash = mag_mod.constraint_conflict(
            scenario, _dosage.get("pli"), _dosage.get("loeuf"))
        if _clash:
            _requirement["constraint_conflict"] = _clash
            log.warning("declared scenario %r sits badly with pLI %s",
                        scenario, _dosage.get("pli"))

    _compound_breadth_info = ({} if demo else _optional(
        "LINCS compound breadth", lambda: _compound_breadth(outputs), {}))

    _brake = {} if demo else _brake_exons(client, gene, species, goal)
    _mitotic_gene_flag = (cellcycle.is_mitotic_gene(gene) or any(
        cellcycle.is_mitotic_gene(p.gene_symbol) for p in promoters))

    _result = {
        "direct_therapies": [t.to_dict() for t in direct_hits],
        # Exons only decay-prone isoforms use: candidates for an antisense
        # oligonucleotide that raises output. Annotation, not a design.
        "brake_exons": _brake,
        # True for a gene whose job is mitosis; see cellcycle.MITOTIC_GENES.
        "mitotic_gene": _mitotic_gene_flag,
        # How many genes each candidate compound moves in LINCS, as a caution
        # about specificity. Annotation only; see signatures.signature_breadth.
        "compound_breadth": _compound_breadth_info,
        # What is already approved for the condition this gene causes. Not a
        # route to the gene and not part of any score; see standard_of_care.
        "standard_of_care": soc_mod.for_gene(gene),
        "tss_window_overlap": tss_overlap,
        "evidence_sources": ledger.to_dict(),
        "evidence_verdict": prov_mod.verdict(
            ledger, (outputs[0] or {}).get("scan_quality") if outputs else None),
        "query": query_block(_args, n_motifs=len(motif_set)),
        # Whether the reports show only FDA approved drugs. The full set
        # stays in this file; see fda.py.
        "fda_filter": fda_mod.status(
            demo=demo, requested=fda_only, tissue=tissue,
            activity_tissue=activity_tissue,
            include_local_only=include_local_only),
        # How much the gene has to move, and whether its dose is under
        # selection. Both are about the GENE rather than about any
        # route, so they sit beside the query rather than in a block.
        "requirement": _requirement,
        "dosage": _dosage,
        # What produced this, as against what it was asked for. A run
        # from six months ago against a since-replaced ReMap file is
        # otherwise indistinguishable from one against the current one.
        # What the run cost, and whether a ceiling bound it. `refused`
        # above zero means evidence is missing because the caller
        # capped it, not because a source had nothing - and every table
        # below looks the same either way.
        "spend": client.spend(),
        "provenance": prov_mod.run_provenance({
            name: _args.get(name) for name in RUN_FILE_ARGS}),
        "results": outputs,
        "weights": WEIGHTS,
        "drug_weights": drug_scoring.DRUG_WEIGHTS,
    }
    _progress.step(9, "report assembled",
                   f"{sum(len(o.get('all_routes') or {}) for o in outputs)} "
                   f"route group(s) over {len(outputs)} promoter(s)")
    _progress.finish()
    return _result




def _brake_exons(client, gene: str, species: str, goal: str) -> dict:
    """Candidate brake exons for an up goal, from the Ensembl isoform list.

    Always asks Ensembl, even when a local GTF supplies the gene models: the
    GTF records carry no exon coordinates, and the basic set leaves out most
    decay isoforms. Returns `{}` when the goal is not to raise the gene, and
    a dict with an `error` when the isoforms could not be read, so a failed
    lookup is not reported as "no brake exon".
    """
    if goal != "up":
        return {}
    try:
        record = genes._lookup_symbol(client, gene, species)
        found = poisonexon.find(record)
        # A candidate exon is only worth an oligonucleotide if the isoforms
        # carrying it carry message. Weighed here rather than inside
        # poisonexon so that module stays a pure read of the annotation.
        found["isoform_weight"] = isoforms.weigh(
            record, isoforms.fetch(client, gene))
        return found
    except Exception as exc:                                # noqa: BLE001
        log.warning("brake-exon search for %s failed: %s", gene,
                    str(exc)[:100])
        return {"gene": gene, "candidates": [], "n_nmd": 0,
                "error": f"the isoform list could not be read "
                         f"({str(exc)[:80]})",
                "note": ""}


def _compound_breadth(outputs: list) -> dict:
    """How broadly each candidate compound acts, from the local LINCS files.

    Empty when the files are not on disk. Added to the result as an
    annotation; nothing is scored or filtered on it.
    """
    breadth = sig_mod.signature_breadth()
    if breadth is None:
        return {}
    names = []
    for block in outputs or []:
        for r in ((block.get("all_routes") or {}).get("routes") or []):
            if r.get("drug"):
                names.append(r["drug"])
    found = sig_mod.breadth_for(names, breadth)
    return {"broad_percentile": sig_mod.BROAD_PERCENTILE,
            "min_signatures": sig_mod.MIN_SIGNATURES_FOR_BREADTH,
            "n_signatures": breadth["n_signatures"],
            "median_genes_all": breadth["median_genes_all"],
            "compounds": found}


def _screen_outcome(screens: dict, eligible: set) -> tuple:
    """(ok, detail) for the receptor screen's ledger row.

    The screen runs on receptors that already have an approved drug naming
    their ChEMBL target. A promoter whose factors have none (HTT's: 49
    factors, no approved drug on any) leaves nothing to screen, which is a
    complete answer and not a failure. Recording it as a failure kept the
    run's verdict at "degraded" for a reason that was not a fault.
    """
    if screens:
        return True, f"{len(screens)} receptor(s) screened"
    if not eligible:
        return True, ("no factor has an approved drug that names its target, "
                      "so there was no receptor to screen")
    return False, (f"{len(eligible)} receptor(s) had drugs to build a screen "
                   f"from, and none could be screened")


def _genetics_detail(phenotypes: list, local: dict | None,
                     error: str | None = None) -> str:
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
        if error:
            # The query did not complete. "No association returned" would
            # read as a finding about the gene, and it is not one.
            return (f"the Open Targets query failed ({error[:110]}), so "
                    f"this says nothing about whether the gene has "
                    f"associations")
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
