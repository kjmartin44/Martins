"""How many close relatives share the pocket you are trying to hit?

The axis this adds
------------------
`drug_scoring.selectivity` counts how many proteins a drug binds.
`regulon` counts how many genes that protein controls. `systemic` counts how
many body systems move. None of them asks the question that decides whether a
*selective* ligand can exist at all:

    how many close paralogues share the binding site?

An orthosteric pocket is under selection to bind the same endogenous ligand
across a receptor family, so it is **conserved** and hard to discriminate.
Thyroid hormone binds THRA and THRB because both evolved to bind thyroid
hormone. That is not an off-target defect in any particular drug; it is a
property of the family, and no amount of medicinal chemistry makes it go away
cheaply.

This is visible in a real CCR5 run's own output:

  * THRA **and** THRB both surfaced, and dextrothyroxine, liothyronine and
    levothyroxine each scored against *both*
  * RARA **and** RARG both surfaced, with tretinoin, trifarotene and
    alitretinoin crossing between them

Why not just count ligands?
---------------------------
Because ligand count measures **study effort**, not selectivity — the same
trap `regulon.regulon_evidence` and `drug_scoring.evidence_depth` already
correct for. A factor with two known ligands usually means nobody did the
chemistry. Worse, the causation runs backwards: nuclear receptors have many
ligands *because* they have a pocket, and the pocket is why they are
druggable at all. Few ligands correlates with *no pocket*, which is the
hardest case, not the cleanest.

Family size is the variable that actually predicts cross-reactivity.

What this does NOT mean
-----------------------
**A large family is a risk, not a verdict.** Subtype-selective agents exist
for several of these families and were often developed precisely because the
pan-agonist was unusable — resmetirom for THRB, trifarotene for RARG. Where
such an agent exists it is named, because "this family is promiscuous" is the
wrong conclusion to draw about a drug specifically designed not to be.

**A small family is not selectivity either.** A receptor with no paralogue
can still be hit by a drug that binds forty unrelated proteins. That is the
`selectivity` axis, and this one says nothing about it.

**Whole-protein similarity is not pocket similarity.** These groupings are by
curated receptor family, which tracks the ligand-binding domain far better
than a global sequence identity would — but two members can share a family
and still differ in the residues that matter. A binding panel is the test;
this is a prior.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Family:
    name: str
    members: frozenset
    ligand: str                  # what the shared pocket evolved to bind
    note: str                    # the cross-reactivity reality
    selective: tuple = ()        # agents that DID achieve subtype selectivity


# Curated because the relationship that matters is the ligand-binding domain,
# which a sequence-identity cutoff approximates poorly. Ensembl paralogy can
# refine this when its homology endpoint is reachable; it is not required.
FAMILIES: tuple = (
    Family("thyroid hormone receptor", frozenset({"THRA", "THRB"}),
           "triiodothyronine (T3)",
           "T3 binds both; the classic pan-agonist problem, and the reason "
           "cardiac effects limited thyromimetics for decades",
           ("resmetirom (THRB-selective, liver-directed)",)),
    Family("retinoic acid receptor", frozenset({"RARA", "RARB", "RARG"}),
           "all-trans retinoic acid",
           "tretinoin and alitretinoin are pan-RAR; subtype selectivity is "
           "achievable but most approved retinoids do not have it",
           ("trifarotene (RARG-selective)", "tamibarotene (RARA-selective)")),
    Family("retinoid X receptor", frozenset({"RXRA", "RXRB", "RXRG"}),
           "9-cis retinoic acid",
           "RXR is also the obligate heterodimer partner for RAR, PPAR, VDR, "
           "LXR and FXR — so an RXR ligand perturbs all of those axes, which "
           "is a broader problem than family cross-reactivity alone",
           ()),
    Family("peroxisome proliferator-activated receptor",
           frozenset({"PPARA", "PPARD", "PPARG"}),
           "fatty acids and eicosanoids",
           "fibrates and glitazones have meaningful but incomplete subtype "
           "preference; dual and pan agonists are common",
           ("seladelpar (PPARD)", "pioglitazone (PPARG-preferring)")),
    Family("estrogen receptor", frozenset({"ESR1", "ESR2"}),
           "estradiol",
           "most approved ER drugs bind both; ERbeta-selective agents have "
           "largely failed to reach the clinic",
           ()),
    Family("estrogen-related receptor", frozenset({"ESRRA", "ESRRB", "ESRRG"}),
           "no confirmed endogenous ligand (orphan)",
           "constitutively active orphans — ligands are synthetic and "
           "cross-reactivity within the family is poorly characterised",
           ()),
    Family("3-ketosteroid receptor",
           frozenset({"NR3C1", "NR3C2", "PGR", "AR"}),
           "3-ketosteroids (cortisol, aldosterone, progesterone, DHT)",
           "the most clinically consequential cross-reactivity in the NR "
           "superfamily: spironolactone hits AR and PGR, progestins hit GR "
           "and MR, mifepristone hits PGR and GR",
           ("eplerenone (MR-selective vs spironolactone)",)),
    Family("xenobiotic sensor", frozenset({"NR1I2", "NR1I3", "VDR"}),
           "structurally diverse lipophiles",
           "PXR and CAR have famously promiscuous pockets evolved to detect "
           "many xenobiotics; VDR is in the same subfamily but far more "
           "discriminating, so calcitriol is not a PXR/CAR problem",
           ("calcitriol (VDR, does not meaningfully engage PXR/CAR)",)),
    Family("liver X receptor", frozenset({"NR1H2", "NR1H3"}),
           "oxysterols",
           "LXRalpha and LXRbeta share ligands; hepatic lipogenesis via "
           "LXRalpha is the dose-limiting effect",
           ()),
    Family("farnesoid X receptor", frozenset({"NR1H4"}), "bile acids",
           "no close paralogue — obeticholic acid is not fighting a "
           "conserved pocket, which is part of why FXR was tractable",
           ()),
    Family("pregnane-adjacent singleton", frozenset({"NR1H5P"}), "orphan",
           "pseudogene in human; listed so it is not silently treated as a "
           "druggable relative of FXR", ()),
    Family("REV-ERB", frozenset({"NR1D1", "NR1D2"}), "heme",
           "paralogous circadian repressors with shared ligand chemistry",
           ()),
    Family("photoreceptor / TLX orphan", frozenset({"NR2E1", "NR2E3"}),
           "orphan", "distant orphans with no approved pharmacology", ()),
    Family("atypical NR0B", frozenset({"NR0B1", "NR0B2"}),
           "no DNA-binding domain",
           "DAX-1 and SHP lack a DBD and act as corepressors, so promoter "
           "logic applies to them differently from the rest of the family",
           ()),
    Family("germ cell nuclear factor", frozenset({"NR6A1"}), "orphan",
           "no close paralogue, but also no ligand pharmacology", ()),
    Family("RAR-related orphan receptor",
           frozenset({"RORA", "RORB", "RORC"}), "sterols",
           "RORgamma-t inverse agonists are the developed class; "
           "cross-reactivity with RORalpha is a recurring liability",
           ()),
    Family("COUP-TF", frozenset({"NR2F1", "NR2F2", "NR2F6"}),
           "orphan", "orphan repressors with no established pharmacology",
           ()),
    Family("HNF4", frozenset({"HNF4A", "HNF4G"}), "fatty acids",
           "structurally similar, and effectively undrugged", ()),
    Family("NR4A orphan", frozenset({"NR4A1", "NR4A2", "NR4A3"}),
           "no classical ligand pocket",
           "the pocket is filled with bulky side chains; ligands act at an "
           "alternate surface, so family logic applies weakly here",
           ()),
    Family("NR5A", frozenset({"NR5A1", "NR5A2"}), "phospholipids",
           "SF-1 and LRH-1 share phospholipid binding", ()),
    Family("TR2/TR4", frozenset({"NR2C1", "NR2C2"}), "orphan",
           "orphan receptors, no approved pharmacology", ()),
)

_INDEX: dict = {m: f for f in FAMILIES for m in f.members}


def family_of(symbol: str) -> Family | None:
    return _INDEX.get((symbol or "").strip().upper())


def paralog_specificity(symbol: str) -> dict:
    """Can a ligand at this receptor plausibly be subtype-selective?

    Returns a 0..1 score where **higher = fewer competing relatives = a
    selective ligand is more achievable**, matching the direction of
    `regulon.transcriptional_specificity` and `delivery.expression_breadth`
    so there is one convention across every specificity axis.

    The curve is a heuristic, not a fit: 1 / (1 + 0.5 * n_relatives), giving
    1.0 alone, 0.67 with one relative, 0.5 with two, 0.4 with three.
    """
    sym = (symbol or "").strip().upper()
    fam = family_of(sym)
    if not fam:
        return {"score": None, "family": None, "n_relatives": None,
                "relatives": [], "assessed": False,
                "note": (f"{sym} is not in the curated receptor-family table. "
                         f"Most transcription factors have no ligand pocket "
                         f"at all, so this is usually a real absence rather "
                         f"than a gap — but it is UNASSESSED, not selective.")}

    others = sorted(fam.members - {sym})
    score = round(1.0 / (1.0 + 0.5 * len(others)), 3)
    note = (f"{sym} shares the {fam.name} pocket ({fam.ligand}) with "
            f"{len(others)} relative(s): {', '.join(others)}. {fam.note}."
            if others else
            f"{sym} has no close paralogue in the {fam.name} family, so a "
            f"selective ligand is not fighting a conserved pocket.")
    return {"score": score, "family": fam.name, "n_relatives": len(others),
            "relatives": others, "assessed": True, "note": note,
            "ligand": fam.ligand,
            "selective_agents": list(fam.selective)}


def cross_reactive_pairs(symbols) -> list:
    """Which of these factors are relatives of each other?

    Two family members both scoring well is not two independent routes — one
    ligand will usually hit both, so they are one bet. In a real CCR5 run
    THRA/THRB and RARA/RARG each surfaced as separate rows.
    """
    present = [s.upper() for s in symbols if family_of(s)]
    seen, out = set(), []
    for a in present:
        fam = family_of(a)
        group = sorted(fam.members & set(present))
        key = (fam.name, tuple(group))
        if len(group) > 1 and key not in seen:
            seen.add(key)
            out.append({
                "family": fam.name, "members": group, "ligand": fam.ligand,
                "note": (f"{', '.join(group)} are relatives sharing the "
                         f"{fam.ligand} pocket — a ligand at one will "
                         f"usually engage the others, so these are one bet, "
                         f"not {len(group)} independent routes"),
                "selective_agents": list(fam.selective),
            })
    return out
