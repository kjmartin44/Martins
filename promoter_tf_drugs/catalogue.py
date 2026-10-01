"""Every route that acts on your gene, in one table, on one scale.

The problem this fixes
----------------------
Each module in this package answers its own question and prints its own
section: direct agents, graded TF drugs, allosteric modulators, stabilisers,
upstream network routes, measured LINCS hits, virtual-screen hypotheses. They
were never compared, because each carried a different score on a different
scale. You could not answer "what are ALL the ways to move this gene, and which
is best?" — only "what did each module find?"

This collects them into one ranking.

Why it is not one number
------------------------
Collapsing everything into a single score would be the wrong fix, because these
claims are not the same kind of claim. "Approved, with clinical data on this
gene" and "a docking program liked this pose" differ by orders of magnitude in
how much they should move you, and averaging them destroys exactly the
information you need.

So every route carries an **evidence tier**, and the table sorts by tier first
and score second. A tier-1 route with a mediocre score beats a tier-5 route
with a great one, every time. The tiers:

  1 clinical      approved for this gene; the effect is documented in humans
  2 measured      observed experimentally — LINCS signature, ChIP-backed
                  binding at this promoter
  3 curated       a human curator annotated the mechanism (ChEMBL, GO); real
                  evidence, but about the drug-protein pair, not about your gene
  4 predicted     a model says so — motif scan, ligand-based direction
                  classifier, signed network path
  5 hypothesis    generated, untested — virtual screen hit, docking pose

Anything at tier 4 or 5 is an experiment you have not done yet. The point of
showing them together is to see how far you are reaching, not to hide it.

Directness
----------
Also tracked, because it is nearly independent of evidence: how many steps from
the drug to your gene. 0 = the drug acts on the gene or its transcript.
1 = acts on a TF that binds your promoter. 2+ = acts upstream of that TF. Each
step costs specificity, and the table shows the cost rather than burying it.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, asdict, field

from .discovery import Axes

log = logging.getLogger(__name__)

TIERS = {
    1: "clinical — approved for this gene, documented in humans",
    2: "measured — the drug's target engagement is observed (LINCS, ChIP, or a measured receptor binder)",
    3: "curated — mechanism annotated by a human curator",
    4: "predicted — a model says so",
    5: "hypothesis — generated, untested",
}

# Categories of thing that can act on a gene's expression.
CATEGORIES = [
    "direct: gene or transcript",
    "direct: gene product",
    "TF: orthosteric ligand",
    "TF: allosteric modulator",
    "TF: stabiliser or glue",
    "TF: endogenous ligand",
    "TF: degradation blocker",
    "upstream signalling",
    "empirical signature",
    "chromatin-level",
    "virtual screen hypothesis",
]

# Approved chromatin-level drugs. They will move almost any promoter, which is
# both why they belong here and why their specificity score is dismal.
CHROMATIN_DRUGS = [
    ("Vorinostat", "HDAC inhibitor", "up", 2006),
    ("Romidepsin", "HDAC inhibitor", "up", 2009),
    ("Belinostat", "HDAC inhibitor", "up", 2014),
    ("Panobinostat", "HDAC inhibitor", "up", 2015),
    ("Valproic acid", "weak HDAC inhibitor", "up", 1978),
    ("Azacitidine", "DNMT inhibitor", "up", 2004),
    ("Decitabine", "DNMT inhibitor", "up", 2006),
    ("Tazemetostat", "EZH2 inhibitor", "up", 2020),
]


# What each category actually is, why you would pick it, and what it costs.
# The `differs` field is the important one: several of these sound alike and
# behave very differently.
CATEGORY_GUIDE = {
    "direct: gene or transcript": {
        "changes": "abundance",
        "assay": "RT-qPCR then western/ELISA",
        "what": "Acts on your gene's own DNA or RNA — siRNA, antisense, "
                "splice-switching oligos, gene therapy.",
        "why": "Sequence-level specificity. Affects your gene and essentially "
               "nothing else, which no protein-targeted route can match.",
        "cost": "Delivery. Solved for liver (GalNAc) and CNS (intrathecal); "
                "not solved for immune cells, muscle or most tumours.",
        "differs": "Unlike 'direct: gene product', this changes how much "
                   "protein is made rather than what existing protein does.",
    },
    "direct: gene product": {
        "changes": "ACTIVITY ONLY",
        "assay": "a functional/activity assay — a western blot will show nothing",
        "what": "Binds the protein your gene encodes — an antibody, enzyme "
                "inhibitor or receptor ligand acting on the product itself.",
        "why": "Immediate, well-understood pharmacology; the protein is right "
               "there rather than three inferential steps away.",
        "cost": "Usually changes the protein's *activity*, not its abundance. "
                "If your readout is a western blot, it can look like nothing "
                "happened.",
        "differs": "Does not touch expression at all — it is here because it "
                   "changes functional output, which is often what you "
                   "actually want.",
    },
    "TF: orthosteric ligand": {
        "changes": "abundance",
        "assay": "RT-qPCR then western/ELISA",
        "what": "Binds the transcription factor's natural ligand pocket.",
        "why": "The best-characterised route, and where nearly all approved "
               "TF pharmacology lives.",
        "cost": "Only ~48 human TFs (the nuclear receptors) have such a "
                "pocket, and the factor's whole regulon moves with your gene.",
        "differs": "Clamps the receptor on or off regardless of endogenous "
                   "signalling — contrast the allosteric and ligand-axis "
                   "routes, which scale what is already happening.",
    },
    "TF: allosteric modulator": {
        "changes": "abundance",
        "assay": "RT-qPCR then western/ELISA",
        "what": "Binds a site away from the ligand pocket, amplifying or "
                "damping the endogenous ligand's effect (PAM/NAM).",
        "why": "Allosteric sites are not conserved across a receptor family, "
               "so subtype selectivity is achievable where orthosteric "
               "chemistry has failed. Circadian and tissue patterning survive.",
        "cost": "Bounded maximal effect, and it does nothing without "
                "endogenous ligand present.",
        "differs": "Same receptor as the orthosteric route, opposite "
                   "trade-off: less maximal effect, more selectivity.",
    },
    "TF: endogenous ligand": {
        "changes": "abundance",
        "assay": "RT-qPCR then western/ELISA; also measure the endogenous ligand, which is what actually moved",
        "what": "Changes how much natural ligand exists, by drugging the "
                "enzyme that makes or destroys it.",
        "why": "The enzyme is often far more druggable than the receptor, and "
               "some controlling enzymes are tissue-restricted, giving "
               "selectivity the receptor cannot.",
        "cost": "The ligand acts on every receptor it binds, everywhere — a "
                "systemic hormone change has a systemic phenotype. Homeostatic "
                "feedback also partly opposes you.",
        "differs": "Never touches the receptor. Anastrozole and finasteride "
                   "work this way and appear in no other category.",
    },
    "TF: stabiliser or glue": {
        "changes": "abundance",
        "assay": "RT-qPCR for your gene; note the DRUG acts on the factor's activity or level, not your gene's",
        "what": "Locks the factor in an active conformation, or glues it to a "
                "partner protein.",
        "why": "Needs no ligand pocket, so it reaches factors that orthosteric "
               "chemistry cannot touch at all.",
        "cost": "Mostly preclinical for TF targets; the approved examples are "
                "on other protein classes.",
        "differs": "Raises activity per molecule with protein levels "
                   "unchanged — the opposite of a degradation blocker.",
    },
    "TF: degradation blocker": {
        "changes": "abundance",
        "assay": "RT-qPCR for your gene; western for the FACTOR to confirm it accumulated",
        "what": "Stops the factor being destroyed, so more of it accumulates.",
        "why": "The class that actually shipped: NRF2 via KEAP1, HIF via PHD "
               "inhibitors.",
        "cost": "Effect is capped by the synthesis rate, and stabilising a "
                "factor stabilises it everywhere.",
        "differs": "Raises protein *abundance*, visible on a western blot — "
                   "unlike a conformational stabiliser, which does not.",
    },
    "upstream signalling": {
        "changes": "abundance",
        "assay": "RT-qPCR then western/ELISA; confirm the pathway node moved too, or you cannot attribute the effect",
        "what": "Drugs a kinase, phosphatase or receptor upstream of the "
                "transcription factor.",
        "why": "Turns 'this factor has no drugs' into a real list. Most TFs "
               "are only reachable this way.",
        "cost": "Specificity falls with every step: a node three hops up "
                "reaches your gene through one branch and thousands of other "
                "genes through the rest.",
        "differs": "The only category where the drug never touches your gene "
                   "or its TF — read the printed path and judge it yourself.",
    },
    "empirical signature": {
        "changes": "abundance (measured)",
        "assay": "already measured as mRNA — confirm by RT-qPCR in your own cells",
        "what": "A compound observed to move your gene in a cell-based "
                "transcriptomic assay, mechanism unknown.",
        "why": "The only evidence here that is measured rather than inferred. "
               "Where it agrees with a predicted route, you have a real "
               "candidate.",
        "cost": "No mechanism, so no way to predict whether it transfers to "
                "your cell type. Check the cell line and whether the value was "
                "measured or model-inferred.",
        "differs": "Every other category starts from mechanism and predicts an "
                   "effect. This starts from an observed effect and has no "
                   "mechanism.",
    },
    "chromatin-level": {
        "changes": "abundance, genome-wide",
        "assay": "RT-qPCR plus a housekeeping panel; if everything moved, nothing specific happened",
        "what": "HDAC, DNMT and EZH2 inhibitors — genome-wide transcriptional "
                "effects.",
        "why": "Reliable positive control that your reporter *can* be induced.",
        "cost": "Specificity ~0.05. These move thousands of promoters. Almost "
                "never an acceptable therapeutic route on that basis alone.",
        "differs": "Not targeted at your gene in any sense; it is in the table "
                   "so you can see what a non-specific option scores.",
    },
    "virtual screen hypothesis": {
        "changes": "unknown — predicted binding only",
        "assay": "a binding assay first; expression is two steps further on",
        "what": "An approved drug predicted to bind one of your TFs, never "
                "tested against it.",
        "why": "Real repurposing space — most approved drugs have been assayed "
               "against a tiny fraction of the proteome.",
        "cost": "An extrapolation beyond the training chemistry. Expect most "
                "to be wrong; each is a binding assay, not a finding.",
        "differs": "The only category with no prior evidence of the "
                   "drug–target pair existing at all.",
    },
}


@dataclass
class Route:
    drug: str
    category: str
    mechanism: str
    acts_on: str                 # the protein or nucleic acid the drug binds
    tier: int
    directness: int              # steps from drug to your gene
    axes: Axes
    approved: bool
    source: str                  # which module produced this
    score: float = 0.0
    notes: list = field(default_factory=list)

    def to_dict(self) -> dict:
        d = asdict(self)
        d["axes"] = self.axes.to_dict()
        d["tier_label"] = TIERS[self.tier]
        return d


def _score(r: Route) -> float:
    """Within a tier, rank on the axes and penalise distance from the gene."""
    base = 100.0 * (0.4 * r.axes.strength + 0.35 * r.axes.specificity
                    + 0.25 * r.axes.efficiency)
    return round(base / (1.0 + 0.25 * r.directness), 1)


# ------------------------------------------------------------- collectors
def from_direct(direct_hits: list[dict]) -> list[Route]:
    out = []
    for t in direct_hits:
        mech_text = (t.get("mechanism") or "").lower()
        # Small molecules can act on the transcript too - risdiplam is a
        # splicing modifier, not a drug against the SMN protein.
        nucleic = t["modality"] in (
            "siRNA", "antisense oligonucleotide",
            "splice-switching oligonucleotide", "gene therapy"
        ) or any(k in mech_text for k in ("splic", "mrna", "pre-mrna",
                                          "transcript", "exon"))
        r = Route(
            drug=t["drug"],
            category="direct: gene or transcript" if nucleic
                     else "direct: gene product",
            mechanism=t["mechanism"], acts_on=t["gene"],
            tier=1 if t["approved"] else 3, directness=0,
            axes=Axes(0.85 if t["raises_protein"] else 0.55,
                      t["specificity"], 0.75, "measured"),
            approved=t["approved"], source="direct")
        if not t["raises_protein"]:
            r.notes.append("raises activity, not abundance — a western blot "
                           "will show nothing")
        if t.get("withdrawn"):
            r.notes.append("WITHDRAWN")
            r.tier = 3
        out.append(r)
    return out


def from_grades(grades: list[dict], mechanisms: dict | None = None) -> list[Route]:
    """Graded TF drugs, split by orthosteric vs allosteric where known."""
    allosteric: dict[str, str] = {}
    for sym, cmp in (mechanisms or {}).items():
        for h in cmp.get("hits", []):
            if h.get("is_allosteric"):
                allosteric[h["drug"].lower()] = h["mechanism"]

    out = []
    for g in grades:
        c = g["components"]
        mech = allosteric.get(g["drug_name"].lower())
        out.append(Route(
            drug=g["drug_name"],
            category=("TF: allosteric modulator" if mech
                      else "TF: orthosteric ligand"),
            mechanism=mech or "ligand at the transcription factor",
            acts_on=g["tf_symbol"], tier=3, directness=1,
            axes=Axes(round(0.45 * c.get("potency", .5)
                            + 0.3 * c.get("action", .5)
                            + 0.25 * g["tf_confidence"], 3),
                      round(c.get("selectivity", .5) * c.get("regulon", .55), 3),
                      round(0.4 * c.get("potency", .5)
                            + 0.25 * c.get("maturity", .5)
                            + 0.2 * c.get("tissue", 1.0)
                            + 0.15 * c.get("evidence_depth", .5), 3),
                      "predicted"),
            approved=True, source="drug_scoring",
            notes=list(g.get("notes", []))[:2]))
    return out


def from_stabilisers(stabilisation: dict) -> list[Route]:
    out = []
    for sym, a in (stabilisation or {}).items():
        for st in a.get("stabilisers", []):
            approved = st["status"] == "approved"
            cat = ("TF: degradation blocker"
                   if st["mechanism"] == "degradation blocker"
                   else "TF: stabiliser or glue")
            r = Route(
                drug=st["drug"], category=cat, mechanism=st["mechanism"],
                acts_on=st["target_node"],
                tier=3 if approved else 5, directness=1,
                axes=Axes(0.75 if approved else 0.6, 0.65,
                          0.7 if approved else 0.4,
                          "curated" if approved else "predicted"),
                approved=approved, source="stabilizers")
            r.notes.append(f"raises {st['raises']}")
            if not approved:
                r.notes.append(f"NOT approved ({st['status']})")
            out.append(r)
    return out


def from_ligands(ligand_axes: dict) -> list[Route]:
    """Drugs that change the receptor's natural ligand rather than binding it."""
    out = []
    for sym, a in (ligand_axes or {}).items():
        for r in a.get("routes", []):
            approved = r["status"] == "approved"
            # One step: the drug binds an enzyme, the enzyme sets the ligand,
            # the ligand sets the receptor. Tissue-selective enzymes score
            # better because the ligand change stays local.
            spec = 0.6 if r["tissue_selective"] else 0.3
            rt = Route(
                drug=r["drug"], category="TF: endogenous ligand",
                mechanism=f"{r['strategy']} ({r['ligand']}) via {r['acts_on']}",
                acts_on=r["acts_on"],
                tier=3 if approved else 5, directness=2,
                axes=Axes(0.7 if approved else 0.5, spec,
                          0.65 if approved else 0.35,
                          "curated" if approved else "predicted"),
                approved=approved, source="ligands")
            rt.notes.extend(r.get("caveats", [])[:2])
            if r["tissue_selective"]:
                rt.notes.append("controlling enzyme is not expressed "
                                "everywhere: the ligand change is local")
            out.append(rt)
    return out


def from_upstream(upstream: dict) -> list[Route]:
    out = []
    for sym, hits in (upstream or {}).items():
        for h in hits:
            bp = h["best_path"]
            for d in h.get("drugs", []):
                out.append(Route(
                    drug=d["name"], category="upstream signalling",
                    mechanism=f"{'activate' if h['net_sign'] == 1 else 'inhibit'} "
                              f"{h['node']} → {sym} ({bp['hops']} hops)",
                    acts_on=h["node"], tier=4,
                    directness=1 + bp["hops"],
                    axes=Axes(round(bp["confidence"], 3),
                              round(bp["specificity"], 3), 0.6, "predicted"),
                    approved=bool(d.get("approved")), source="network",
                    notes=[bp["description"]]))
    return out


def from_lincs(lincs: dict) -> list[Route]:
    out = []
    for h in (lincs or {}).get("hits", []):
        if h.get("non_specific"):
            continue
        # Three levels, because CLUE records three and a binary collapsed
        # the middle one - which demoted every hit on a non-landmark gene
        # (CCR5, VDR, HBB are all "best inferred") to tier 4, BELOW a
        # curated ligand annotation. That understated it: the signature was
        # measured, and only the value for this gene was imputed, whereas a
        # tier-3 ligand route is annotation on two separate steps.
        fs = h.get("feature_space")
        measured = h.get("is_landmark") is True or fs == "landmark"
        if measured:
            tier, basis = 2, "measured"
        elif fs == "best inferred":
            tier, basis = 3, "imputed from measured landmarks"
        else:
            tier, basis = 4, "predicted"

        extra = []
        if measured:
            pass
        elif fs == "best inferred":
            extra.append("value IMPUTED from the 978 measured landmarks, not "
                         "measured for this gene - high confidence among "
                         "inferred genes, still not a measurement")
        else:
            extra.append("value INFERRED, not measured")

        out.append(Route(
            drug=h["drug"], category="empirical signature",
            mechanism="moved this gene in a cell-based transcriptomic assay",
            acts_on="(unknown / whole cell)",
            tier=tier, directness=0,
            axes=Axes(min(1.0, abs(h["score"]) / 5.0), 0.40, 0.6, basis),
            approved=True, source="signatures",
            notes=([f"cell line {h['cell_line']}"] if h.get("cell_line") else [])
                  + extra))
    return out


def from_screen(screens: dict) -> list[Route]:
    out = []
    for sym, sc in (screens or {}).items():
        if not (sc.get("enrichment") or {}).get("usable"):
            continue
        for h in sc.get("hits", []):
            quad = h["quadrant"]

            if quad == "ESTABLISHED LIGAND, untested on this gene":
                # Tier 2, not 3. What makes these different from every other
                # curated TF ligand is that the drug-receptor step is
                # MEASURED - that is what "established ligand" means - and
                # only the receptor-to-gene step is a hypothesis. A curated
                # tier-3 route is annotation on both halves.
                #
                # They stay BELOW tier 1: a direct agent acts on the gene
                # itself and skips the receptor-to-gene inference entirely,
                # so it is a stronger claim however well-characterised the
                # binder is. Below the direct drugs, above everything
                # curated, is exactly where the evidence puts them.
                out.append(Route(
                    drug=h["drug"], category="TF: orthosteric ligand",
                    mechanism=(f"binds {sym} (measured), and {sym} is "
                               f"predicted to act on this promoter - so the "
                               f"drug would reach your gene THROUGH {sym}, "
                               f"not by touching it directly"),
                    acts_on=sym, tier=2, directness=1,
                    axes=Axes(0.7, 0.5, 0.65, "measured binding"),
                    approved=True, source="screening",
                    notes=[f"half-measured: the {h['drug']}-{sym} interaction "
                           f"has assay data behind it, so the pharmacology is "
                           f"de-risked. What is untested is whether moving "
                           f"{sym} moves your gene - the one step a binding "
                           f"assay cannot tell you",
                           f"everything else {sym} regulates moves too, which "
                           f"a direct agent on your gene would avoid"]))
                continue

            if quad != "NOVEL HYPOTHESIS":
                continue
            out.append(Route(
                drug=h["drug"], category="virtual screen hypothesis",
                mechanism=f"predicted ligand at {sym}, never tested there",
                acts_on=sym, tier=5, directness=1,
                axes=Axes(round(h["probability"], 3), 0.5, 0.5, "predicted"),
                approved=True, source="screening",
                notes=["extrapolation beyond training chemistry; requires a "
                       "binding assay before it means anything"]))
    return out


def chromatin_routes(goal: str = "up") -> list[Route]:
    """Approved chromatin drugs. They move almost any promoter — which is the
    point, and the problem."""
    if goal != "up":
        return []
    out = []
    for name, mech, direction, year in CHROMATIN_DRUGS:
        r = Route(
            drug=name, category="chromatin-level", mechanism=mech,
            # Not "directness 0": acting genome-wide is not acting on your gene.
            # Scored as remote, because that is what it is.
            acts_on="chromatin (genome-wide)", tier=3, directness=3,
            axes=Axes(0.70, 0.05, 0.45, "curated"),
            approved=True, source="catalogue")
        r.notes.append(
            "acts genome-wide: will move thousands of promoters including "
            "yours. Useful as a positive control that your reporter CAN be "
            "induced; almost never an acceptable therapeutic route on "
            "specificity grounds")
        out.append(r)
    return out


# ------------------------------------------------------------------ driver
def build(result_block: dict, *, goal: str = "up",
          include_chromatin: bool = True) -> dict:
    """Assemble every route from one promoter's result block."""
    routes: list[Route] = []
    routes += from_direct(result_block.get("direct_therapies") or [])
    routes += from_grades(result_block.get("drug_grades") or [],
                          result_block.get("mechanisms"))
    routes += from_stabilisers(result_block.get("stabilisation") or {})
    routes += from_ligands(result_block.get("ligand_axes") or {})
    routes += from_upstream(result_block.get("upstream") or {})
    routes += from_lincs(result_block.get("lincs") or {})
    routes += from_screen(result_block.get("screens") or {})
    if include_chromatin:
        routes += chromatin_routes(goal)

    for r in routes:
        r.score = _score(r)

    # Tier first, always. A weak tier-1 route beats a brilliant tier-5 one.
    # Within a tier, score alone: directness is already priced into the score,
    # and sorting on it again would double-count it (which floated genome-wide
    # chromatin drugs above far better TF ligands).
    routes.sort(key=lambda r: (r.tier, -r.score))

    by_category: dict[str, int] = {}
    by_tier: dict[int, int] = {}
    for r in routes:
        by_category[r.category] = by_category.get(r.category, 0) + 1
        by_tier[r.tier] = by_tier.get(r.tier, 0) + 1

    best_per_category = {}
    for r in routes:
        if r.category not in best_per_category:
            best_per_category[r.category] = r.to_dict()

    return {
        "routes": [r.to_dict() for r in routes],
        "n_routes": len(routes),
        "by_category": by_category,
        "by_tier": {TIERS[k]: v for k, v in sorted(by_tier.items())},
        "best_per_category": best_per_category,
        "categories_absent": [c for c in CATEGORIES if c not in by_category],
    }
