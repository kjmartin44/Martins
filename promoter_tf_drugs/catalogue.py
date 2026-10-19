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
    "direct: predicted product binder",
    "promoter DNA binder",
    "TF: orthosteric ligand",
    "TF: inferred activity change",
    "TF: abundance change",
    "TF: predicted pocket binder",
    "TF: allosteric modulator",
    "TF: stabiliser or glue",
    "TF: endogenous ligand",
    "TF: degradation blocker",
    "upstream signalling",
    "off-pathway pocket binder",
    "measured phenocopy",
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
    "TF: predicted pocket binder": {
        "changes": "abundance",
        "assay": "a binding assay first, then a reporter for the factor, "
                 "then RT-qPCR for your gene",
        "what": "A compound predicted by docking to fit the factor's "
                "ligand pocket, from a structure rather than from any "
                "measured activity.",
        "why": "It needs no prior ligands for the receptor, so it can "
               "propose a scaffold the ligand-based model could not "
               "recognise. A compound picked by both is supported by two "
               "arguments sharing no input.",
        "cost": "Weakest evidence in the report. The pocket was shown to "
                "rank its known binders above decoys, which is not the "
                "same as the score being right for any one compound, and "
                "docking says nothing about agonist versus antagonist.",
        "differs": "Unlike 'TF: orthosteric ligand', nothing here was "
                   "measured at the receptor - the argument is entirely "
                   "from the shape of the pocket.",
    },
    "off-pathway pocket binder": {
        "changes": "abundance",
        "assay": "establish the link to your gene first; a binding assay "
                 "for this protein tells you nothing about the promoter "
                 "until that path exists",
        "what": "A docking run against a protein that is not one of the "
                "factors found at this promoter. Usually an upstream node "
                "someone chose to dock by hand.",
        "why": "It may still be the right target - an indirect route can "
               "be more druggable than the factor itself - but this run "
               "is evidence about a pocket, not about a path.",
        "cost": "The route to your gene is asserted by whoever launched "
                "the dock, not demonstrated here. The network section is "
                "where a path is evidenced, and this compound does not "
                "appear there.",
        "differs": "Unlike 'TF: predicted pocket binder', the protein "
                   "docked is not among this promoter's factors, so even "
                   "a perfect pose leaves the connection to your gene "
                   "unestablished.",
    },
    "measured phenocopy": {
        # Abundance, not activity. This acts on the FACTOR, so your gene's
        # transcription moves and RT-qPCR sees it. "activity" is reserved
        # for a drug on your gene's own product, which changes what the
        # existing protein does without moving how much there is.
        "changes": "abundance",
        "assay": "a binding assay against the factor, then a reporter, "
                 "then RT-qPCR for your gene - the correlation names a "
                 "candidate and establishes none of those steps",
        "what": "A compound whose effect on cell viability across the "
                "DepMap panel tracks what happens when this factor is "
                "knocked out, so it behaves as though it removes the "
                "factor.",
        "why": "It needs no annotation of any kind. Every other route "
               "here starts from something a curator, a motif or a "
               "structure already asserted, so none of them can surface a "
               "drug nobody has connected to the factor. This one can, "
               "because it reads a measured phenotype.",
        "cost": "Correlation across cell lines, not target engagement. A "
                "hit may act on the factor, on anything else the same "
                "lines need, or on a pathway the factor sits in, and the "
                "first is the least likely of the three. Both panels "
                "score proliferation of transformed cells, so a factor "
                "that matters in a differentiated cell leaves no trace.",
        "differs": "Unlike 'empirical signature', nothing was measured "
                   "about your gene or about transcription at all - the "
                   "readout is whether cells live.",
    },
    "TF: abundance change": {
        "changes": "abundance",
        "assay": "RT-qPCR and a western for the FACTOR first, then "
                 "RT-qPCR for your gene",
        "what": "A compound measured to change how much of the factor's "
                "own transcript exists, rather than what the factor does.",
        "why": "The only route here that needs the factor to be expressed "
               "rather than druggable, so it reaches factors with no "
               "ligand, no structure and no curated upstream kinase - "
               "which is most of them.",
        "cost": "More factor is not more factor ACTIVITY. For a "
                "ligand-activated receptor the two are nearly unrelated, "
                "and transcript is not protein. Broad agents that shut "
                "down transcription move every factor and are flagged.",
        "differs": "Unlike 'empirical signature', which is a compound "
                   "measured against YOUR GENE, this was measured against "
                   "the FACTOR - one further inferential step from what "
                   "you care about.",
    },
    "TF: inferred activity change": {
        "changes": "abundance",
        "assay": "RT-qPCR for your gene, plus a reporter or ChIP for the "
                 "factor to show the effect runs through it",
        "what": "A compound whose own LINCS signature moves this factor's "
                "known targets, found without any record that it binds the "
                "factor.",
        "why": "It is the only category here that can reach a factor with "
               "no medicinal chemistry behind it. Every other route needs "
               "somebody to have assayed a compound at the target first, "
               "which is why the other tables fill up with nuclear "
               "receptors.",
        "cost": "The attribution is weak. The regulon moved; which of the "
                "factors that share those targets moved it is not settled, "
                "and a compound that simply arrests the cell cycle shifts "
                "every proliferation-linked regulon.",
        "differs": "Unlike 'empirical signature', which is a drug measured "
                   "to move YOUR GENE, this is a drug measured to move the "
                   "FACTOR'S TARGETS - one further inferential step from "
                   "the thing you care about.",
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
    "direct: predicted product binder": {
        "changes": "activity",
        "assay": "a binding assay against the protein first, then a "
                 "functional assay for direction - RT-qPCR will not see "
                 "this route at all, and a western blot reads it as no "
                 "effect",
        "what": "A compound predicted to bind this gene's own protein "
                "and not recorded as acting on it - either because its "
                "chemistry resembles the protein's known ligands, or "
                "because it docks into the protein's pocket.",
        "why": "It is the only route here aimed at the gene's own "
               "product rather than at a factor upstream of it, so it "
               "skips every inference about whether moving a factor moves "
               "your gene. And unlike the curated product drugs it can "
               "propose a compound no one has annotated.",
        "cost": "Predicted binding and nothing more. The model is "
                "refused below a cross-validated AUC floor, which bounds "
                "how badly it can be ranking, but a hit is a hypothesis "
                "for a binding assay. Direction is separate again: "
                "whether binding inhibits or activates the protein "
                "decides whether it serves your goal.",
        "differs": "Unlike 'direct: gene product', nothing was measured "
                   "or curated for this compound at this protein - that "
                   "category is annotation, this is a model. Unlike 'TF: "
                   "predicted pocket binder', the pocket is on the gene's "
                   "own product, so no inference about a factor "
                   "controlling the gene is involved.",
    },
    "promoter DNA binder": {
        "changes": "abundance",
        "assay": "an EMSA or footprint on this promoter first, then a "
                 "reporter, then RT-qPCR - and a counter-screen on an "
                 "unrelated promoter carrying the same feature, because "
                 "selectivity is the whole question here",
        "what": "A small molecule with a documented preference for a "
                "sequence feature this promoter carries - a quadruplex, a "
                "GC-box, an AT tract - proposed because the feature is "
                "enriched against shuffles of this same sequence.",
        "why": "It is the only route that acts on the promoter itself "
               "rather than on a protein, so it needs no druggable factor "
               "and no annotation. Where a promoter depends on one "
               "feature, occupying it is a direct intervention.",
        "cost": "Sequence preference is not selectivity. These agents "
                "bind the feature wherever it occurs and most of the "
                "genome's copies are not in your gene, so the dose that "
                "engages your promoter engages thousands of others. "
                "Nearly every approved member is cytotoxic, which is how "
                "it works rather than a side effect. And a motif is not a "
                "folded structure in a cell.",
        "differs": "Unlike every TF route, no protein is targeted - which "
                   "removes the question of whether the factor controls "
                   "your gene and replaces it with whether the molecule "
                   "can tell your promoter from any other.",
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
    # `approved` above is what the SOURCE claims. Ten collectors hardcode
    # it to True, so on its own it means "this came from a drug table",
    # not "this is an approved drug".
    #
    # This is the checked version: True when the name is in a source of
    # approved drugs on this machine, False when such a source was read
    # and the name is absent, None when there was nothing to check
    # against. None is not False - with no drug table on disk every route
    # is unverified, and reading that as "not approved" would empty the
    # report rather than correct it.
    approval_verified: bool | None = None
    # The multiplier applied to `score` because the compound moves a large
    # share of the transcriptome in LINCS (see signatures.signature_breadth).
    # None when no penalty applied, which is the usual case and also the case
    # for a compound the LINCS files do not contain.
    breadth_penalty: float | None = None
    # The multiplier applied because the compound disturbs mitosis and the
    # gene being treated is a mitotic one (see systemic.MITOTIC_DISRUPTORS).
    mitotic_penalty: float | None = None

    def to_dict(self) -> dict:
        d = asdict(self)
        d["axes"] = self.axes.to_dict()
        d["tier_label"] = TIERS[self.tier]
        return d


# A route through a compound in the top decile of LINCS breadth keeps its
# tier and has its score halved, so within a tier it sits below narrower
# compounds. Not applied to tier: a broad compound with the best evidence is
# still the best evidence, and the flag beside it says the cost. The factor is
# a judgement, not a calibrated figure.
BROAD_SCORE_FACTOR = 0.5
# For a mitotic gene, a compound that disturbs mitosis loses score the same
# way and for the same kind of reason. Applied once: a compound that is both
# very broad and a mitotic disruptor is halved, not quartered.
MITOTIC_SCORE_FACTOR = 0.5


def _score(r: Route) -> float:
    """Within a tier, rank on the axes and penalise distance from the gene."""
    base = 100.0 * (0.4 * r.axes.strength + 0.35 * r.axes.specificity
                    + 0.25 * r.axes.efficiency)
    return round(base / (1.0 + 0.25 * r.directness), 1)


# ------------------------------------------------------------- collectors
_APPROVED_NAMES: set | None = None
# Whether a real drug LIBRARY was read, as opposed to only the curated
# tables this package ships. It decides whether a miss means anything:
# the curated tables hold a few hundred names, so a drug absent from them
# is usually just uncurated. Vemurafenib is approved and is not in them.
_APPROVED_FROM_LIBRARY = False


def approved_name_set() -> set | None:
    """Lowercased salt-free names of drugs known to be approved, or None.

    Built once from whatever is on disk, most to least complete: a merged
    local drug table, then the curated indication and drug-class tables
    this package ships. None when none of them is available, which must
    stay distinguishable from an empty set.

    It reads the DEFAULT library locations (see `localdrugs.resolve`), not
    a file passed with `--library`, and it is built once per process: a
    library that becomes readable later is not picked up until the next
    run. A read failure is logged when it happens.
    """
    global _APPROVED_NAMES, _APPROVED_FROM_LIBRARY
    if _APPROVED_NAMES is not None:
        return _APPROVED_NAMES or None

    from .congeners import salt_free_name
    names: set = set()
    try:
        from . import localdrugs as _ld
        rows = _ld.library()
        for row in rows:
            n = salt_free_name((row.get("name") or "").lower())
            if n:
                names.add(n)
        _APPROVED_FROM_LIBRARY = bool(rows)
    except Exception as exc:                                # noqa: BLE001
        # Said once, here. This used to be `pass`, so a locked, corrupt or
        # half-written library silently left approval checks resting on a
        # few hundred curated names - the effect (a vaguer "unknown" for
        # every drug not curated) is graceful, but nothing said WHY it had
        # changed. The set is built once per process, so this is the only
        # place it can be reported.
        log.warning("the local approved-drug library could not be read "
                    "(%s: %s); approval checks fall back to the %s curated "
                    "names and answer 'unknown' for anything else",
                    type(exc).__name__, str(exc)[:100], "few hundred")
    try:
        from . import indications as _ind
        names.update(salt_free_name(k.lower()) for k in _ind.CURATED)
    except Exception as exc:                                # noqa: BLE001
        log.warning("the curated indication table could not be read "
                    "(%s: %s)", type(exc).__name__, str(exc)[:100])
    # NOT systemic.CLASSES. That is a burden table, and it deliberately
    # contains investigational and tool compounds - CGP-60474, AT-7519,
    # ALW-II-38-3 - because their PHARMACOLOGY is what the burden claim is
    # about. Reading it as a list of approved drugs made every one of them
    # verify as approved, which is the opposite of what this function is
    # for and exactly the confusion it was added to end.
    names.discard("")
    _APPROVED_NAMES = names
    return names or None


def verify_approval(drug: str) -> bool | None:
    """Is this name in a source of approved drugs on this machine?

    None when nothing was available to check against. A tool compound
    like CGP-60474 and an approved drug nobody curated look identical to
    this function, so a False means "not in the sources read", never
    "proven unapproved".
    """
    pool = approved_name_set()
    if pool is None:
        return None
    from .congeners import salt_free_name
    if salt_free_name((drug or "").lower()) in pool:
        return True
    # A miss is only informative against a real library. The curated
    # tables alone are a few hundred names, so a miss there says the drug
    # is uncurated and nothing about its approval.
    return False if _APPROVED_FROM_LIBRARY else None


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
        # Tier 2 is "measured receptor binder" (see TIERS). A graded drug with
        # a measured potency of at least 1 uM against this factor and a
        # recorded action is that, exactly as a screen-found ligand is; it
        # used to be tier 3 only because the factor was not among the few
        # receptors screened, so the tier followed which screens ran.
        ev = g.get("evidence") or {}
        p = ev.get("pchembl_median")
        measured = (p is not None and p >= MEASURED_BINDING_PCHEMBL
                    and c.get("action", 0) >= 0.7)
        tier = 2 if measured else 3
        out.append(Route(
            drug=g["drug_name"],
            category=("TF: allosteric modulator" if mech
                      else "TF: orthosteric ligand"),
            mechanism=mech or "ligand at the transcription factor",
            acts_on=g["tf_symbol"], tier=tier, directness=1,
            axes=Axes(round(0.45 * c.get("potency", .5)
                            + 0.3 * c.get("action", .5)
                            + 0.25 * g["tf_confidence"], 3),
                      round(c.get("selectivity", .5) * c.get("regulon", .55), 3),
                      round(0.4 * c.get("potency", .5)
                            + 0.25 * c.get("maturity", .5)
                            + 0.2 * c.get("tissue", 1.0)
                            + 0.15 * c.get("evidence_depth", .5), 3),
                      "measured binding" if measured else "predicted"),
            approved=True, source="drug_scoring",
            notes=(([f"measured: pChEMBL {p} against {g['tf_symbol']} "
                     f"(n={ev.get('pchembl_n')}), with a recorded action"]
                    if measured else [])
                   + list(g.get("notes", []))[:2])))
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
            # Which way the DRUG must move the node:
            #     n = want_tf_sign * net_sign
            # This read net_sign alone, the same inverted expression
            # `network.attach_drugs` had, so the report said "inhibit" for
            # every repressor factor and every downregulation goal while
            # the pipeline selected an activator. The two disagreed on
            # screen.
            want_tf = h.get("want_tf_sign")
            if want_tf is None:
                # A results.json written before the sign was recorded.
                # Saying so beats printing a direction that may be wrong.
                verb = "move"
            else:
                verb = ("activate" if want_tf * h["net_sign"] == 1
                        else "inhibit")
            hops = bp["hops"]
            for d in h.get("drugs", []):
                out.append(Route(
                    drug=d["name"], category="upstream signalling",
                    mechanism=f"{verb} {h['node']} → {sym} "
                              f"({hops} hop{'' if hops == 1 else 's'})",
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


def from_docking(docking: dict, *, goal: str = "up",
                 tf_roles: dict | None = None,
                 gene: str | None = None) -> list[Route]:
    """Hits from a docking run performed somewhere else.

    Tier 4. The compound is predicted to fit a pocket, which is two
    inferences away from your gene: whether it binds, and then whether
    binding moves the factor the way you need. A docking score is also a
    poor estimator of affinity even when the pocket validated, so the
    axes carry the validation rather than the score.

    `docking` is one run's result dict, as written by
    `docking.dock_from_bundle`. A run that did not validate contributes
    nothing, because its scores carry no information.

    `tf_roles` maps this promoter's factors to their roles. It decides two
    things that were previously not asked:

    * whether the protein docked is one of this promoter's factors at all.
      A dock against an upstream node used to come back labelled "TF:
      predicted pocket binder" with directness 1, claiming a direct
      relationship with a gene it has no established path to.
    * how the pocket's conformation relates to the direction you want.
      Docking still says nothing about whether a compound agonises or
      antagonises; what it can say is that the pocket it was docked into
      is the conformation you want stabilised, or the opposite one.
    """
    runs = docking if isinstance(docking, list) else [docking]
    out = []
    for run in runs:
        if not run or not run.get("usable"):
            continue
        tf = run.get("target") or "?"
        auc, ef1 = run.get("auc"), run.get("ef1")
        pocket = run.get("pocket") or {}
        state = pocket.get("state")

        roles = {k.upper(): v for k, v in (tf_roles or {}).items()}
        # The gene's own product is the shortest route there is, and the
        # off-pathway test below would have called it the longest: it is
        # not one of the promoter's factors, so it was labelled as
        # establishing no route to the gene it IS.
        is_product = bool(gene) and tf.upper() == gene.strip().upper()
        on_pathway = is_product or (not roles) or tf.upper() in roles
        # For the product, moving it IS the goal - no factor sits between
        # it and the gene - so the conformation wanted follows the goal
        # directly rather than through a factor's role.
        role = "activator" if is_product else roles.get(tf.upper())

        # How the docked conformation relates to the direction wanted.
        # "supports" and "opposes" are about the POCKET, never about the
        # compound: a molecule that fits the active conformation may still
        # be an antagonist of it.
        fit, fit_note = "unaddressed", None
        if state not in ("active", "inactive"):
            # Asked and unanswerable, not never asked. The conformation
            # is read from bound coactivator or corepressor peptides,
            # which only nuclear receptors carry, so for a GPCR, a
            # kinase or an enzyme it cannot be determined at all.
            fit_note = (
                f"the conformation of this structure could not be read, "
                f"so nothing here says whether the pocket is the one you "
                f"want stabilised. That is a limit of the method - the "
                f"conformation is judged by a bound coactivator or "
                f"corepressor peptide, which only nuclear receptors "
                f"carry - and not a finding about {tf}")
        if state in ("active", "inactive") and role:
            from .drugs import desired_direction
            want = desired_direction(role, goal)
            want_state = ("active" if want == "up" else
                          "inactive" if want == "down" else None)
            if want_state:
                fit = "supports" if state == want_state else "opposes"
                fit_note = (
                    f"to move your gene {goal} through a {role} you want "
                    f"{tf} {'more' if want == 'up' else 'less'} active, so "
                    f"the conformation to stabilise is {want_state}. This "
                    f"pocket is the {state} one"
                    + ("." if fit == "supports" else
                       ". A compound that fits it is as likely to hold the "
                       "factor in the state you do NOT want. That is a "
                       "property of the structure chosen, not a measurement "
                       "of this compound - an antagonist can bind the "
                       "active conformation."))

        for h in run.get("hits", []):
            score = h.get("score")
            if score is None:
                continue
            notes = [
                f"structure {pocket.get('pdb_id')}"
                + (f", ligand {pocket.get('ligand_resname')}"
                   if pocket.get("ligand_resname") else ""),
                "a docking score says the molecule can occupy the "
                "pocket. It does not say whether it activates or "
                "inhibits, and that decides the direction your gene "
                "moves",
            ]
            # Unless something else did say. The direction comes from a
            # classifier over this receptor's own ligands of known
            # action, not from the docking, so it is reported as the
            # separate claim it is - with the confidence attached,
            # because a low-confidence call is a different thing from a
            # measurement.
            if h.get("direction") in ("agonist", "antagonist"):
                conf = h.get("direction_confidence") or "unstated"
                notes.append(
                    f"predicted {h['direction']} at {tf}, {conf} "
                    f"confidence, from structural similarity to ligands "
                    f"of known action rather than from the docking. A "
                    f"prediction about the receptor, not a measurement "
                    f"of this compound")
            if fit_note:
                notes.append(fit_note)
            if is_product:
                notes.append(
                    f"this pocket is on {tf} itself, not on a factor "
                    f"upstream of it. A compound binding here changes what "
                    f"the existing protein DOES rather than how much of it "
                    f"there is, so RT-qPCR will not see it and a western "
                    f"blot reads it as no effect")
            if not on_pathway:
                notes.append(
                    f"{tf} is not one of the factors found at this "
                    f"promoter, so this run establishes a pocket and not "
                    f"a route to your gene. Whatever connects them was "
                    f"asserted when the dock was launched, and is not "
                    f"evidenced here.")

            out.append(Route(
                drug=h.get("drug") or "?",
                category=("direct: predicted product binder" if is_product
                          else "TF: predicted pocket binder" if on_pathway
                          else "off-pathway pocket binder"),
                mechanism=(f"docks into the {tf} pocket at "
                           f"{score:.1f} kcal/mol"
                           + (f" ({state} conformation, {fit} the direction "
                              f"wanted)" if fit != "unaddressed" else
                              f" ({state} conformation)" if state and
                              state != "unknown" else "")),
                # Directness 1 means "acts on a factor at this promoter".
                # A protein that is not one of them is at least one more
                # step away, and how many more is unknown - 2 is the
                # smallest honest answer, not a measured hop count.
                acts_on=tf, tier=4,
                directness=0 if is_product else 1 if on_pathway else 2,
                # Strength is the POCKET's validated enrichment, not this
                # compound's score: the score separates binders from
                # non-binders and ranks binders badly, so reading it as a
                # per-compound strength would overstate what was done.
                axes=Axes(min(1.0, (auc or 0.0)), 0.3, 0.5,
                          f"pocket validated at AUC {auc}"
                          + (f", {ef1}-fold top-1% enrichment" if ef1
                             else "")),
                approved=True, source="docking",
                notes=notes))
    return out


def from_sensitivity(sens: dict, *, goal: str = "up",
                     tf_roles: dict | None = None) -> list[Route]:
    """Routes from a drug-sensitivity correlation across cell lines.

    Tier 4. The drug's effect on cells is measured, but the claim that it
    acts on THIS factor is inferred from a correlation, so it sits with
    the predicted routes rather than the measured ones - a curated ligand
    at tier 3 at least asserts a direct interaction, and this asserts
    none.

    Unlike docking, this readout carries a direction. A compound that
    phenocopies the knockout behaves as though it REMOVES the factor, so
    it serves a goal that needs the factor lowered and works against one
    that needs it raised. Routes pointing the wrong way are dropped here,
    which is the same treatment annotated ligands already get - the
    direction is part of what the correlation measured, not a guess about
    a conformation.
    """
    from .drugs import desired_direction

    out = []
    roles = {k.upper(): v for k, v in (tf_roles or {}).items()}
    for sym, res in (sens or {}).items():
        if not res or not res.get("assessed"):
            continue
        if res.get("readout") != "phenocopy":
            continue
        role = roles.get(sym.upper())
        want = desired_direction(role or "unknown", goal)
        n_lines_total = res.get("n_lines") or 0

        for h in res.get("hits", []):
            moves = h.get("direction")
            # "phenocopies" lowers the factor; "opposes" raises it.
            gives = ("down" if moves == "phenocopies" else
                     "up" if moves == "opposes" else None)
            if gives is None:
                continue
            if want in ("up", "down") and gives != want:
                continue

            r, raw = h.get("r"), h.get("r_raw")
            q = h.get("q_value")
            verb = ("behaves as though it removes" if gives == "down"
                    else "behaves as though it raises")
            out.append(Route(
                drug=h.get("drug") or "?",
                category="measured phenocopy",
                mechanism=(f"{verb} {sym}: its effect on cell viability "
                           f"tracks {sym} knockout across "
                           f"{h.get('n_lines')} line(s), partial r = {r}"),
                acts_on=sym, tier=4, directness=1,
                # Strength is the partial correlation, which is the
                # measured quantity. The raw value is in the notes and is
                # deliberately not used: it is mostly cell-line fragility.
                axes=Axes(min(1.0, abs(float(r or 0.0))), 0.4, 0.5,
                          f"partial r = {r} at q = {q}"),
                approved=True, source="sensitivity",
                notes=[
                    f"correlation across {h.get('n_lines')} of "
                    f"{n_lines_total} shared cell line(s), q = {q}",
                    f"raw correlation was {raw}; the reported {r} controls "
                    f"for how easily each line dies under the panel as a "
                    f"whole. That control moves the number either way - it "
                    f"removes a shared component that can inflate a "
                    f"correlation or mask one - so a large gap means the "
                    f"raw value was mostly about the cell lines",
                    "no record of this compound binding the factor. The "
                    "argument is that it behaves as though it does, which "
                    "is a reason to run a binding assay rather than a "
                    "substitute for one",
                    "a proliferation readout in cancer lines: it says "
                    "nothing about a differentiated tissue",
                ]))
    return out


def from_product_screen(screen: dict, gene: str) -> list[Route]:
    """Predicted binders of the gene's own protein product.

    Tier 5. A model trained on the product's known ligands says the
    chemistry fits; nothing was measured for this compound at this
    protein.

    Only the NOVEL HYPOTHESIS rows become routes. The established-ligand
    rows are the drugs already annotated against the product, which
    `from_direct` reports at tier 1 with far better standing - emitting
    them here too would double-count the same drug and bury the one thing
    this screen adds, which is the compounds nobody has annotated.
    """
    if not screen or not screen.get("is_gene_product"):
        return []
    if not (screen.get("enrichment") or {}).get("usable"):
        return []

    sym = (gene or screen.get("gene") or "?").upper()
    auc = screen.get("model_auc")
    out = []
    for h in screen.get("hits", []):
        if h.get("quadrant") != "NOVEL HYPOTHESIS":
            continue
        out.append(Route(
            drug=h.get("drug") or "?",
            category="direct: predicted product binder",
            mechanism=(f"predicted to bind the {sym} protein itself, by a "
                       f"model trained on its known ligands"),
            # Directness 0: it acts on the gene's own product, with no
            # factor in between. That is the shortest route in the report
            # and the tier is what says how little is known about it.
            acts_on=sym, tier=5, directness=0,
            axes=Axes(round(float(h.get("probability") or 0.0), 3), 0.6,
                      0.5, f"binding model, cross-validated AUC {auc}"),
            approved=True, source="product_screen",
            notes=list(h.get("notes") or [])[:3]))
    return out


def from_promoter_dna(dna: dict) -> list[Route]:
    """Agents preferring a sequence feature this promoter is enriched for.

    Tier 5. The feature is measured - it is in the sequence, and it beat a
    shuffled null - but that an agent preferring it would move this gene
    is entirely untested. Directness 0 because it acts at the promoter
    with no protein in between, which makes it the shortest route in the
    report and among the least evidenced, and the tier is what carries
    that.

    Specificity is pinned low on purpose. These bind their feature
    wherever it occurs, so the axis that asks how little else is touched
    is the one they fail.
    """
    hits = (dna or {}).get("candidates") or []
    out = []
    for h in hits:
        approved = bool(h.get("approved"))
        out.append(Route(
            drug=h.get("drug") or "?",
            category="promoter DNA binder",
            mechanism=(f"{h.get('mechanism')}; this promoter carries "
                       f"{h.get('n_sites')} {h.get('feature')} "
                       f"(expected {h.get('expected')} by composition, "
                       f"p={h.get('p_value')})"),
            acts_on="promoter DNA", tier=5, directness=0,
            axes=Axes(0.4, 0.15, 0.4,
                      f"sequence feature enriched at p={h.get('p_value')}"),
            approved=approved, source="promoter_dna",
            notes=[
                h.get("note") or "",
                "binds this feature wherever it occurs in the genome, not "
                "only in your promoter - selectivity is the question to "
                "ask first, and a counter-screen on an unrelated promoter "
                "carrying the same feature is the experiment",
                "a motif in sequence is not a folded structure in a cell",
            ] + ([] if approved else
                 ["not an approved drug; included because it is the "
                  "reference compound for this feature"])))
    return out


def from_abundance(hits: list) -> list[Route]:
    """Compounds that change how much of the factor exists.

    Tier 3, the same standing as a curated mechanism and one step below
    `from_lincs`. There the compound was measured against YOUR GENE; here
    it was measured against the factor, and that the factor controls your
    gene is the inference this rests on.

    A promiscuous compound - one moving most of the factors asked about -
    is kept and demoted rather than dropped, which is how plicamycin and
    the other broad agents are already handled.
    """
    out = []
    for h in (hits or []):
        ns = h.get("non_specific")
        nr = h.get("nuclear_receptor")
        # Specificity is the axis that carries the two things wrong with
        # this route: a promiscuous compound, and a receptor whose
        # abundance says little about its activity.
        spec = 0.45
        if ns:
            spec = 0.1
        elif nr:
            spec = 0.25
        out.append(Route(
            drug=h.get("drug") or "?",
            category="TF: abundance change",
            mechanism=(f"{'raises' if h.get('moved_factor') == 'up' else 'lowers'} "
                       f"{h.get('tf')} transcript in "
                       f"{h.get('n_signatures')} profile(s)"),
            acts_on=h.get("tf"),
            tier=int(h.get("tier") or 4), directness=1,
            axes=Axes(min(1.0, (h.get("n_signatures") or 1) / 10.0),
                      spec, 0.5,
                      f"measured on {h.get('tf')}'s own transcript, "
                      f"{h.get('n_signatures')} profile(s)"),
            approved=True, source="abundance",
            notes=list(h.get("caveats") or [])[:3]))
    return out


def from_activity_drugs(activity_drugs: dict) -> list[Route]:
    """Compounds whose LINCS signature moves a factor's regulon.

    Tier 3, not 2, and the distinction is the whole point. Tier 2 is "the
    drug's target engagement is observed". Here nothing about the drug's
    engagement with the factor is observed at all: what was measured is
    that the factor's TARGETS moved when the compound was given. The
    observation is real and the attribution to this particular factor is
    an inference, which is the same standing as a curated mechanism.

    Specificity is scored low on purpose. A compound reaching this table
    did so by moving a whole regulon, which is the opposite of acting
    narrowly, and the effect can run through anything upstream of the
    factor or through general toxicity.
    """
    out = []
    for tf, hits in (activity_drugs or {}).items():
        for h in hits:
            n_sig, n_tot = h["n_significant"], h["n_signatures"]
            out.append(Route(
                drug=h["drug"],
                category="TF: inferred activity change",
                mechanism=(f"{h['direction']} this factor's regulon across "
                           f"{n_sig} of {n_tot} LINCS profiles "
                           f"(median z {h['median_z']:+.2f})"),
                acts_on=tf,
                tier=3, directness=1,
                # Strength from how reproducible it was, not from one
                # profile: the repetition is the evidence here.
                axes=Axes(min(1.0, abs(h["median_z"]) / 4.0),
                          0.25, 0.5,
                          f"regulon shift in {n_sig}/{n_tot} profiles"),
                approved=True, source="tfactivity",
                notes=([f"consistent in {h['consistency']:.0%} of profiles, "
                        f"{len(h['cell_lines'])} cell line(s)"]
                       + list(h.get("caveats") or [])[:2])))
    return out


# A measured binder is "measured" for the engagement it shows. pChEMBL 6 is
# 1 uM, below which the binding is too weak to call the drug engaged at a
# dose anyone could give.
MEASURED_BINDING_PCHEMBL = 6.0
# The weight a drug that is a known binder with no recorded action gets for
# "action" in a strength: the same 0.40 `drug_scoring` gives an unknown action.
UNKNOWN_ACTION = 0.40


def _potency01(pchembl) -> float:
    """pChEMBL 5 -> 0 and 9 -> 1, as `drug_scoring._potency` reads it."""
    if pchembl is None:
        return 0.5
    return max(0.0, min(1.0, (float(pchembl) - 5.0) / 4.0))


def tf_info_from_tfs(tfs: list) -> dict:
    """`{symbol: {confidence, role, directions}}` from a result's `tfs` list.

    `confidence` is the factor's own score on 0 to 1, `role` its role at this
    promoter, and `directions` maps each of its drugs (by normalised name) to
    the way the drug moves the factor. The pipeline builds the same thing from
    its ranked factors; this reads it back from a saved result so routes can
    be rebuilt from one.
    """
    from .signatures import breadth_key
    out = {}
    for t in tfs or []:
        out[t["symbol"]] = {
            "confidence": (t.get("score") or 0.0) / 100.0,
            "role": t.get("role"),
            "directions": {breadth_key(d["name"]): d.get("direction")
                           for d in t.get("drugs", [])}}
    return out


def from_screen(screens: dict, tf_info: dict | None = None,
                goal: str = "up") -> list[Route]:
    """Routes from the receptor ligand screens.

    `tf_info` (see `tf_info_from_tfs`) lets an established ligand be scored
    the way a graded one is: potency from its own measured pChEMBL, the
    factor's confidence, and whether the drug is known to move the factor the
    way the goal needs. Without it the routes keep their old constant axes.
    """
    out = []
    for sym, sc in (screens or {}).items():
        if not (sc.get("enrichment") or {}).get("usable"):
            continue

        # A family model was trained on this receptor's RELATIVES and never
        # on the receptor. `transfer_auc` is the measured answer to whether
        # that chemistry carries over, checked against ligands the model
        # never saw - so it, not the model's own probability, is what a
        # route's strength may rest on.
        family = sc.get("family") if sc.get("family_model") else None
        t_auc = sc.get("transfer_auc")
        relatives = ", ".join(c["symbol"] for c in
                              (sc.get("contributing") or [])[:4])

        for h in sc.get("hits", []):
            quad = h["quadrant"]

            if (quad == "ESTABLISHED LIGAND, untested on this gene"
                    and tf_info is not None):
                route = _established_ligand_route(h, sym, tf_info, goal,
                                                  bool(family))
                if route is not None:
                    out.append(route)
                continue

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
                           f"a direct agent on your gene would avoid"]
                    + ([f"surfaced by a {family} family model, but the tier "
                        f"rests on the measured {h['drug']}-{sym} binding "
                        f"rather than on the model - which ranked it, and "
                        f"did not establish it"] if family else [])))
                continue

            if quad != "NOVEL HYPOTHESIS":
                continue
            if family:
                # Strength is the MEASURED transfer AUC, not the model's
                # own probability. That probability comes from a model
                # that never saw this receptor, so reading it as a
                # strength at this receptor is the overclaim the transfer
                # check exists to prevent.
                strength = round(float(t_auc), 3) if t_auc else 0.3
                out.append(Route(
                    drug=h["drug"], category="virtual screen hypothesis",
                    mechanism=(f"predicted ligand at {sym} by a "
                               f"{family} family model - trained on "
                               f"{relatives or 'relatives'}, not on {sym}"),
                    acts_on=sym, tier=5, directness=1,
                    axes=Axes(strength, 0.4, 0.5,
                              f"family model, transfer AUC {t_auc}"),
                    approved=True, source="screening",
                    notes=[
                        f"the model was never trained on {sym}. It ranks "
                        f"{sym}'s own known ligands at AUC {t_auc}, which "
                        f"is why it is reported at all - but a family "
                        f"model predicts FAMILY binding, and subtype "
                        f"selectivity is exactly what it cannot see",
                        "extrapolation beyond training chemistry; requires "
                        "a binding assay before it means anything",
                    ]))
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


def _established_ligand_route(h: dict, sym: str, tf_info: dict, goal: str,
                              family: bool) -> Route | None:
    """A tier-2 route for a receptor's measured ligand, scored like a graded one.

    The strength is the one `from_grades` uses (0.45 potency, 0.30 action, 0.25
    factor confidence), so the same drug at the same factor scores the same
    whichever collector found it. Three things were missing before: every such
    route had strength 0.7 whatever its potency, so dozens tied; a drug known
    to move the factor the wrong way ranked beside one known to move it the
    right way; and one with no recorded action ranked as high as a recorded
    agonist. Returns None for a drug known to move the factor against the
    goal, which would move the gene the wrong way.
    """
    from .drugs import desired_direction
    from .signatures import breadth_key
    info = tf_info.get(sym) or {}
    want = desired_direction(info.get("role") or "unknown", goal)
    direction = (info.get("directions") or {}).get(breadth_key(h["drug"]))
    known = direction in ("up", "down")
    if known and want != "either" and direction != want:
        return None
    action = 1.0 if known else UNKNOWN_ACTION
    conf = info.get("confidence")
    conf = 0.5 if conf is None else conf
    pchembl = h.get("measured_pchembl")
    strength = round(0.45 * _potency01(pchembl) + 0.30 * action
                     + 0.25 * conf, 3)
    notes = [f"half-measured: the {h['drug']}-{sym} interaction has assay "
             f"data behind it"
             + (f" (pChEMBL {pchembl})" if pchembl is not None else "")
             + ", so the pharmacology is de-risked. What is untested is "
             f"whether moving {sym} moves your gene - the one step a binding "
             f"assay cannot tell you",
             f"everything else {sym} regulates moves too, which a direct "
             f"agent on your gene would avoid"]
    if not known:
        notes.append(f"which way {h['drug']} moves {sym} is not recorded, so "
                     f"it may push the gene the wrong way; scored with an "
                     f"unknown action")
    if family:
        notes.append("surfaced by a family model, but the tier rests on the "
                     "measured binding and not on the model")
    return Route(
        drug=h["drug"], category="TF: orthosteric ligand",
        mechanism=(f"binds {sym} (measured), and {sym} is predicted to act on "
                   f"this promoter - so the drug would reach your gene "
                   f"THROUGH {sym}, not by touching it directly"),
        acts_on=sym, tier=2, directness=1,
        axes=Axes(strength, 0.5, 0.65, "measured binding"),
        approved=True, source="screening", notes=notes)


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
          include_chromatin: bool = True, breadth: dict | None = None,
          mitotic_gene: bool = False) -> dict:
    """Assemble every route from one promoter's result block.

    `breadth` is `signatures.breadth_for`'s output, keyed by
    `signatures.breadth_key`. A route through a compound marked `very_broad`
    has its score multiplied by `BROAD_SCORE_FACTOR` and says so in its notes.
    Without it nothing is penalised, so a run with no LINCS files, and demo
    mode, score exactly as before.

    `mitotic_gene` says the gene being treated is one whose job is mitosis. A
    route through a compound that disturbs mitosis then has its score
    multiplied by `MITOTIC_SCORE_FACTOR` and says so. It never changes a tier.
    """
    routes: list[Route] = []
    routes += from_direct(result_block.get("direct_therapies") or [])
    routes += from_grades(result_block.get("drug_grades") or [],
                          result_block.get("mechanisms"))
    routes += from_stabilisers(result_block.get("stabilisation") or {})
    routes += from_ligands(result_block.get("ligand_axes") or {})
    routes += from_upstream(result_block.get("upstream") or {})
    routes += from_lincs(result_block.get("lincs") or {})
    routes += from_screen(
        result_block.get("screens") or {},
        tf_info=result_block.get("tf_info") or None, goal=goal)
    routes += from_activity_drugs(
        result_block.get("activity_drugs") or {})
    routes += from_docking(result_block.get("docking") or {}, goal=goal,
                           tf_roles=result_block.get("tf_roles"),
                           gene=result_block.get("gene_symbol"))
    routes += from_abundance(
        result_block.get("factor_abundance") or [])
    routes += from_promoter_dna(result_block.get("promoter_dna") or {})
    routes += from_product_screen(
        result_block.get("product_screen") or {},
        result_block.get("gene_symbol") or "")
    routes += from_sensitivity(result_block.get("sensitivity") or {},
                               goal=goal,
                               tf_roles=result_block.get("tf_roles"))
    if include_chromatin:
        routes += chromatin_routes(goal)

    if breadth:
        from .signatures import breadth_key
    if mitotic_gene:
        from .systemic import is_mitotic_disruptor
    for r in routes:
        r.score = _score(r)
        scaled = False
        if mitotic_gene and is_mitotic_disruptor(r.drug):
            r.score = round(r.score * MITOTIC_SCORE_FACTOR, 1)
            r.mitotic_penalty = MITOTIC_SCORE_FACTOR
            scaled = True
            r.notes.append(
                f"score multiplied by {MITOTIC_SCORE_FACTOR}: this compound "
                f"disturbs mitosis directly, and the gene is one whose job "
                f"is mitosis, so it moves the gene by changing the cell's "
                f"cycle, not by treating the cause")
        rec = breadth.get(breadth_key(r.drug)) if breadth else None
        if rec and rec.get("very_broad"):
            if not scaled:
                r.score = round(r.score * BROAD_SCORE_FACTOR, 1)
            r.breadth_penalty = BROAD_SCORE_FACTOR
            r.notes.append(
                f"score multiplied by {BROAD_SCORE_FACTOR}: this compound "
                f"moves a median of {rec['median_genes']:.0f} genes across "
                f"{rec['n_signatures']} LINCS signatures "
                f"({rec['percentile']:.0f}th percentile), so lowering or "
                f"raising this gene is one of many changes it makes")
        # One place rather than eleven, so a new collector cannot forget.
        if r.approval_verified is None:
            r.approval_verified = verify_approval(r.drug)

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
