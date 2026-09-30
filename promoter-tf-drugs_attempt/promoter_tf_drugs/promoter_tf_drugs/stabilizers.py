"""Locking a transcription factor into an active state — without a ligand pocket.

Why this is the interesting route
---------------------------------
Everything else in this package eventually finds a nuclear receptor, because
nuclear receptors have a hydrophobic ligand-binding pocket and almost no other
transcription factor does. That single structural fact is why ~48 of ~1,600
human TFs carry essentially all the approved direct pharmacology, and why the
graded tables keep returning dexamethasone and calcitriol.

Stabilisation sidesteps the pocket entirely. Every TF has **protein-protein
interfaces** and every TF has **degradation machinery**, and both are
druggable. This is the mechanism class that reaches the factors nothing else
can.

Four mechanisms, routinely conflated, doing different things
------------------------------------------------------------
**Conformational stabiliser.** Binds the folded protein and locks it in an
active, high-affinity shape — the thing you described. Raises *activity per
molecule*. Pharmacological chaperones (migalastat for alpha-galactosidase A)
and tetramer stabilisers (tafamidis for transthyretin) are the clean approved
examples; eprenetapopt refolding mutant p53 is the TF-directed one, still
clinical.

**PPI stabiliser / molecular glue.** Does not bind either partner usefully on
its own; glues two proteins together by completing a composite interface.
Rapamycin (FKBP12-mTOR), tacrolimus and ciclosporin (FKBP12/cyclophilin-
calcineurin) are approved proof the class works. Fusicoccin stabilising
14-3-3 complexes with FOXO and TAZ is the TF-directed version.

**Degradation blocker.** Raises *abundance* rather than activity, by stopping
the protein being destroyed. The two approved successes are both here:
PHD inhibitors stabilising HIF (daprodustat, roxadustat) and KEAP1 modifiers
stabilising NRF2 (dimethyl fumarate, omaveloxolone).

**PPI disruptor that releases the TF.** Breaking an inhibitory complex frees
the factor — MDM2 inhibitors releasing p53, KEAP1-NRF2 disruptors. Mechanically
the opposite of a glue, same net direction.

And the mirror image, for downregulation: a glue can *destroy* a TF. Thalidomide
and lenalidomide recruit IKZF1/IKZF3 to cereblon for degradation. Those are
approved, and they are the most convincing demonstration that transcription
factors are druggable when you stop looking for a pocket.

Activity vs abundance — the distinction that decides your assay
---------------------------------------------------------------
A conformational stabiliser raises activity with protein levels unchanged; a
degradation blocker raises levels. If your readout is a western blot, the first
looks like a total failure. Every entry here is flagged, the same way
`direct.py` flags activity-only drugs.

Honest status
-------------
TF-directed stabilisers are mostly **preclinical or clinical, not approved**.
The approved wins are degradation blockers (NRF2, HIF) and glue-degraders
(IKZF). If a curated entry here is not approved, it says so — this is a
direction to pursue, not a prescription to write.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, asdict, field

from .http import Client, HttpError

log = logging.getLogger(__name__)

MECHANISMS = {
    "conformational stabiliser": {
        "raises": "activity",
        "needs_pocket": False,
        "note": "locks the folded protein in an active shape; protein levels "
                "unchanged, so a western blot will show nothing",
    },
    "molecular glue": {
        "raises": "activity",
        "needs_pocket": False,
        "note": "completes a composite interface between two proteins; requires "
                "the partner to be present in your cell type",
    },
    "degradation blocker": {
        "raises": "abundance",
        "needs_pocket": False,
        "note": "raises protein levels by preventing destruction; visible on a "
                "western blot, and the effect ceiling is set by synthesis rate",
    },
    "PPI disruptor (releasing)": {
        "raises": "activity",
        "needs_pocket": False,
        "note": "frees the factor from an inhibitory complex; effect depends on "
                "how much sequestered protein exists to release",
    },
    "glue-degrader": {
        "raises": "none (removes protein)",
        "needs_pocket": False,
        "note": "recruits the factor to an E3 ligase for destruction — a "
                "downregulation mechanism",
    },
}


@dataclass
class Stabiliser:
    tf_symbol: str
    drug: str
    mechanism: str
    direction: str                 # "up" | "down"
    raises: str                    # "activity" | "abundance" | "none (removes protein)"
    target_node: str               # what the drug actually binds
    status: str                    # "approved" | "clinical" | "preclinical"
    approval_year: int | None
    description: str
    notes: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


def _s(tf, drug, mech, direction, node, status, year, desc):
    m = MECHANISMS[mech]
    st = Stabiliser(tf, drug, mech, direction, m["raises"], node, status, year, desc)
    st.notes.append(m["note"])
    if status != "approved":
        st.notes.append(f"NOT approved ({status}) — a direction to pursue, not "
                        "an available treatment")
    return st


# --- Curated. Verify anything you intend to rely on; approvals move. --------
CURATED: dict[str, list[Stabiliser]] = {}
for _e in [
    # ---- approved: degradation blockers, the class that actually shipped ----
    ("NFE2L2", "Dimethyl fumarate", "degradation blocker", "up", "KEAP1",
     "approved", 2013,
     "Modifies KEAP1 cysteines, blocking NRF2 ubiquitination; NRF2 accumulates"),
    ("NFE2L2", "Omaveloxolone", "degradation blocker", "up", "KEAP1",
     "approved", 2023, "KEAP1 inhibitor; NRF2 stabilisation in Friedreich ataxia"),
    ("NFE2L2", "Monomethyl fumarate", "degradation blocker", "up", "KEAP1",
     "approved", 2020, "Active metabolite of dimethyl fumarate"),
    ("HIF1A", "Daprodustat", "degradation blocker", "up", "EGLN1/PHD",
     "approved", 2023,
     "Prolyl hydroxylase inhibitor; blocks VHL-mediated HIF degradation"),
    ("HIF1A", "Roxadustat", "degradation blocker", "up", "EGLN1/PHD",
     "approved", 2018, "PHD inhibitor; HIF stabilisation (EU/Japan/China)"),
    ("HIF1A", "Vadadustat", "degradation blocker", "up", "EGLN1/PHD",
     "approved", 2020, "PHD inhibitor"),
    ("EPAS1", "Roxadustat", "degradation blocker", "up", "EGLN1/PHD",
     "approved", 2018, "PHD inhibition stabilises HIF-2alpha as well as HIF-1alpha"),

    # ---- approved glue-degraders: TFs destroyed, not activated -------------
    ("IKZF1", "Lenalidomide", "glue-degrader", "down", "CRBN", "approved", 2005,
     "Cereblon molecular glue; recruits IKZF1 for ubiquitination and degradation"),
    ("IKZF3", "Lenalidomide", "glue-degrader", "down", "CRBN", "approved", 2005,
     "Same glue interface; IKZF3 is the second neosubstrate"),
    ("IKZF1", "Pomalidomide", "glue-degrader", "down", "CRBN", "approved", 2013,
     "Cereblon glue with a distinct neosubstrate profile"),
    ("IKZF1", "Thalidomide", "glue-degrader", "down", "CRBN", "approved", 1998,
     "The original cereblon glue"),
    ("GLI2", "Arsenic trioxide", "conformational stabiliser", "down", "GLI2",
     "approved", 2000,
     "Binds GLI2 directly and inhibits Hedgehog transcriptional output"),

    # ---- clinical / preclinical, TF-directed -------------------------------
    ("TP53", "Eprenetapopt", "conformational stabiliser", "up", "TP53",
     "clinical", None,
     "APR-246; refolds mutant p53 into a transcriptionally active conformation"),
    ("TP53", "Idasanutlin", "PPI disruptor (releasing)", "up", "MDM2",
     "clinical", None, "Blocks MDM2-p53, releasing p53 from degradation"),
    ("TP53", "Navtemadlin", "PPI disruptor (releasing)", "up", "MDM2",
     "clinical", None, "MDM2 antagonist"),
    ("NR4A2", "Amodiaquine", "conformational stabiliser", "up", "NR4A2",
     "approved", 1948,
     "Binds the Nurr1 ligand-binding domain and enhances its transcriptional "
     "activity; approved as an antimalarial, repurposing interest in Parkinson's"),
    ("NR4A2", "Hydroxychloroquine", "conformational stabiliser", "up", "NR4A2",
     "approved", 1955, "Related Nurr1 binder"),
    ("FOXO3", "Fusicoccin A", "molecular glue", "up", "YWHAB/14-3-3",
     "preclinical", None,
     "Stabilises the 14-3-3 complex; the canonical TF-directed glue"),
    ("RBM39", "Indisulam", "glue-degrader", "down", "DCAF15", "clinical", None,
     "DCAF15 glue recruiting RBM39 for degradation"),
    ("PPARG", "MRL-24", "conformational stabiliser", "up", "PPARG",
     "preclinical", None,
     "Partial agonist that blocks CDK5-mediated Ser273 phosphorylation without "
     "full agonism — activity gained without the classic TZD side effects"),
]:
    tf, drug, mech, direction, node, status, year, desc = _e
    CURATED.setdefault(tf, []).append(
        _s(tf, drug, mech, direction, node, status, year, desc))

# Non-TF exemplars: proof the mechanism class works clinically. Shown when a TF
# has nothing of its own, so the route is not dismissed as theoretical.
EXEMPLARS = [
    _s("—", "Tafamidis", "conformational stabiliser", "up", "TTR", "approved",
       2019, "Binds the transthyretin tetramer interface and prevents "
             "dissociation — the cleanest approved example of locking a protein "
             "into a stable shape"),
    _s("—", "Migalastat", "conformational stabiliser", "up", "GLA", "approved",
       2016, "Pharmacological chaperone; stabilises mutant alpha-galactosidase A "
             "into a foldable, traffickable form"),
    _s("—", "Sirolimus", "molecular glue", "up", "FKBP12-mTOR", "approved", 1999,
       "Glues FKBP12 to mTOR; neither interaction occurs without the drug"),
    _s("—", "Tacrolimus", "molecular glue", "up", "FKBP12-calcineurin",
       "approved", 1994, "Glues FKBP12 to calcineurin"),
]


def curated_for(tf: str, goal: str = "up") -> list[Stabiliser]:
    return [s for s in CURATED.get(tf.upper(), []) if s.direction == goal]


def has_ligand_pocket(tf: str, annotation: dict | None = None) -> bool:
    """Nuclear receptors have a pocket; almost nothing else does."""
    from .annotate import NUCLEAR_RECEPTORS
    if annotation and annotation.get("nuclear_receptor"):
        return True
    return tf.upper() in NUCLEAR_RECEPTORS


def chembl_stabilisers(client: Client, symbol: str, uniprot: str | None = None,
                       *, goal: str = "up",
                       approved_only: bool = True) -> list[Stabiliser]:
    """ChEMBL mechanisms annotated as stabilising rather than agonising."""
    from .drugs import chembl_drugs

    keys = ("STABILISER", "STABILIZER", "CHAPERONE")
    try:
        found = chembl_drugs(client, symbol, uniprot, approved_only)
    except HttpError:
        return []
    out = []
    for d in found:
        action = (d.action_type or "").upper()
        if not any(k in action for k in keys):
            continue
        out.append(Stabiliser(
            tf_symbol=symbol.upper(), drug=d.name,
            mechanism="conformational stabiliser", direction=goal,
            raises="activity", target_node=symbol.upper(),
            status="approved" if d.approved else "clinical",
            approval_year=d.first_approval,
            description=d.mechanism or action))
    return out


def assess(tf: str, *, goal: str = "up", annotation: dict | None = None,
           has_direct_drug: bool = False) -> dict:
    """Is stabilisation the right route for this factor, and what exists?"""
    hits = curated_for(tf, goal)
    pocket = has_ligand_pocket(tf, annotation)

    result = {
        "tf": tf.upper(),
        "has_ligand_pocket": pocket,
        "stabilisers": [s.to_dict() for s in hits],
        "n_approved": sum(1 for s in hits if s.status == "approved"),
        "recommended": False,
        # Is there an actual route for THIS factor, or only a mechanism class
        # that has worked on other proteins? The distinction decides whether
        # this is something you can do or something someone could research.
        "actionable": bool(hits),
        "status": "no route known",
        "notes": [],
    }

    if not pocket and not has_direct_drug:
        if hits:
            # A real, named agent against this factor or its degradation node.
            result["recommended"] = True
            result["status"] = "curated route exists"
            result["notes"].append(
                f"{tf.upper()} has no ligand pocket and no direct drug, but "
                f"there IS a curated stabilisation route for it — see the "
                f"table. This is the case the mechanism exists for.")
        else:
            # The common case, and previously over-claimed. Saying
            # "stabilisation recommended" and then listing four approved
            # drugs against unrelated proteins reads as though something
            # were available. Nothing is.
            result["status"] = "no route known"
            result["notes"].append(
                f"{tf.upper()} has no ligand pocket, no direct drug, and "
                f"**nothing curated against it or its degradation "
                f"machinery**. Stabilisation is not a route you can take "
                f"here — it is a research programme. Making it real needs "
                f"two things nobody has published for this factor: the E3 "
                f"ligase that degrades it, and a ligand for that ligase.")
    elif pocket:
        result["notes"].append(
            f"{tf.upper()} has a ligand-binding pocket, so conventional agonists "
            "are available and are usually the simpler route. Stabilisation is "
            "worth considering where you need activity without the full agonist "
            "phenotype — the PPARG Ser273 story is the example.")

    if any(s.raises == "abundance" for s in hits):
        result["notes"].append(
            "Degradation blockers raise protein LEVELS; conformational "
            "stabilisers raise activity with levels unchanged. Pick your assay "
            "accordingly — a western blot cannot see the second kind.")
    return result
