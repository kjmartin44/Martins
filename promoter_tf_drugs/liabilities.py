"""Structural liabilities: what makes a predicted hit not worth believing.

What this is for, and what it is not for
----------------------------------------
The obvious use of a structural alert is to filter out toxic molecules.
That is the wrong use here, and applying it would throw away the best
candidates in the report.

Every compound in these tables is **already FDA-approved**. It has been
through toxicology, through trials, and into people. A rule-of-five
violation or a reactive group on an approved drug is not news: it is a
property of a molecule that was approved anyway, and plenty of the most
useful drugs in medicine fail these filters outright. Doxorubicin is a
quinone. Cisplatin is a reactive metal complex. Both are alerts and both
are approved for good reasons.

So the approval itself answers most of the safety question, and the
package already reads the parts that approval does not answer - withdrawal
and black-box warnings from the drug annotation, and delivery from
`delivery.py`.

The question these alerts DO answer
-----------------------------------
Whether a *prediction* about a molecule is believable.

PAINS - pan-assay interference compounds - are structures that come up
active in assay after assay for reasons that have nothing to do with the
target: they aggregate, they are redox-cyclers, they quench fluorescence,
they covalently modify protein. A PAINS match does not mean the compound
is dangerous. It means the measured activity that a model was trained on,
or the activity a model is now predicting, may be an artefact of the
assay. That is a reason to distrust the *evidence*, which is exactly what
this report is made of.

So a flag here demotes confidence in a predicted route. It never demotes
an approved drug for being approved.

Where each catalogue comes from
-------------------------------
PAINS (Baell and Holloway 2010) is the interference set and the one that
matters most here. Brenk is a set of groups unwanted in lead discovery,
which is a different and much more aggressive question than "is this
approved drug usable" - it is reported and deliberately not acted on. QED
and the Lipinski and Veber rules are oral-drug-likeness heuristics, useful
for a predicted novel compound and close to meaningless for a drug already
being dosed orally.

RDKit ships all of the catalogues, so none of this needs the network.

What the PAINS filter misses
----------------------------
It is a set of SMARTS patterns, not a judgement, and the set is not
exhaustive. Quercetin matches. **Curcumin does not**, despite being the
most cited interference compound in the literature - its reputation rests
on catechol and Michael-acceptor chemistry rather than on a pattern in
Baell's published set. Verified here against RDKit's catalogue.

So a clean PAINS result is weak evidence of a clean compound, while a
match is strong evidence of a dirty one. Read it one way only. The Brenk
catalogue does flag curcumin, on the Michael acceptor, which is why both
are reported even though only PAINS changes anything.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

log = logging.getLogger(__name__)

_CATALOG = None
_CATALOG_FAILED = False


@dataclass
class Liability:
    """Structural flags for one molecule."""
    drug: str
    pains: list = field(default_factory=list)
    brenk: list = field(default_factory=list)
    qed: float | None = None
    lipinski_violations: int | None = None
    veber_pass: bool | None = None
    # Confidence multiplier for a PREDICTED route, never for an approved
    # drug's safety. 1.0 means nothing here undermines the prediction.
    evidence_multiplier: float = 1.0
    notes: list = field(default_factory=list)
    assessed: bool = False

    def to_dict(self) -> dict:
        return dict(self.__dict__)


def _catalog():
    """The combined PAINS and Brenk catalogue, built once."""
    global _CATALOG, _CATALOG_FAILED
    if _CATALOG is not None or _CATALOG_FAILED:
        return _CATALOG
    try:
        from rdkit.Chem import FilterCatalog
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams
                          .FilterCatalogs.PAINS)
        params.AddCatalog(FilterCatalog.FilterCatalogParams
                          .FilterCatalogs.BRENK)
        _CATALOG = FilterCatalog.FilterCatalog(params)
    except Exception as exc:                                # noqa: BLE001
        log.info("rdkit filter catalogues unavailable (%s); structural "
                 "liabilities are not assessed", str(exc)[:80])
        _CATALOG_FAILED = True
    return _CATALOG


def assess(drug_name: str, smiles: str | None) -> Liability:
    """Structural flags for one molecule. `assessed` is False without a
    structure, which is not the same as a clean result."""
    liab = Liability(drug=drug_name)
    if not smiles:
        liab.notes.append(
            "no structure available, so nothing was checked. This is not a "
            "clean result")
        return liab

    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import Descriptors, QED
        RDLogger.DisableLog("rdApp.*")
        mol = Chem.MolFromSmiles(smiles)
    except ImportError:
        liab.notes.append("rdkit not installed; structural checks skipped")
        return liab
    if mol is None:
        liab.notes.append("structure could not be parsed")
        return liab

    liab.assessed = True

    cat = _catalog()
    if cat is not None:
        for match in cat.GetMatches(mol):
            desc = match.GetDescription() or ""
            props = list(match.GetPropList())
            where = match.GetProp("FilterSet") if "FilterSet" in props else ""
            if "pains" in desc.lower() or "pains" in where.lower():
                liab.pains.append(desc)
            else:
                liab.brenk.append(desc)

    try:
        liab.qed = round(QED.qed(mol), 3)
    except Exception as exc:                                # noqa: BLE001
        log.debug("QED failed for %s: %s", drug_name, type(exc).__name__)
        liab.notes.append("QED could not be computed for this structure")

    # Lipinski, counted rather than pass/fail: an approved drug breaking one
    # rule is unremarkable and the count is what carries information.
    try:
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        liab.lipinski_violations = sum(
            [mw > 500, logp > 5, hbd > 5, hba > 10])
        tpsa = Descriptors.TPSA(mol)
        rot = Descriptors.NumRotatableBonds(mol)
        liab.veber_pass = bool(tpsa <= 140 and rot <= 10)
    except Exception as exc:                                # noqa: BLE001
        # The fields stay None, which means "not computed", not "no
        # violations". Say so on the record rather than only in a log.
        log.debug("drug-likeness descriptors failed for %s: %s", drug_name,
                  type(exc).__name__)
        liab.notes.append("Lipinski and Veber descriptors could not be "
                          "computed for this structure, so their absence "
                          "is not a clean result")

    # Only PAINS moves the multiplier, and only because it undermines the
    # EVIDENCE. Brenk and the drug-likeness rules are reported and not
    # acted on: they ask whether this would make a good lead to start
    # from, which is a question already settled by the drug being approved.
    if liab.pains:
        liab.evidence_multiplier = 0.6
        liab.notes.append(
            f"matches {len(liab.pains)} pan-assay interference "
            f"substructure(s) ({'; '.join(liab.pains[:2])}). PAINS compounds "
            f"read as active across unrelated assays through aggregation, "
            f"redox cycling or fluorescence interference, so a measured or "
            f"predicted activity for this molecule is more likely to be an "
            f"artefact. This lowers confidence in the EVIDENCE and says "
            f"nothing about the drug's safety")
    if liab.brenk:
        liab.notes.append(
            f"matches {len(liab.brenk)} substructure(s) unwanted in lead "
            f"discovery. For an approved drug this is reported and not "
            f"acted on - the question it asks was settled when the compound "
            f"was approved")
    if liab.lipinski_violations:
        liab.notes.append(
            f"{liab.lipinski_violations} Lipinski violation(s). Noted "
            f"because it predicts oral absorption poorly for a NEW compound; "
            f"for one already dosed in people it is a description, not a "
            f"problem")
    return liab


def assess_many(drugs: dict) -> dict:
    """`{name: smiles}` to `{name: Liability}`."""
    return {name: assess(name, smi) for name, smi in (drugs or {}).items()}


def summarise(liabs: dict) -> str:
    if not liabs:
        return ""
    assessed = [x for x in liabs.values() if x.assessed]
    if not assessed:
        return ("No structures were available, so no compound was checked "
                "for assay-interference substructures. That is a gap in the "
                "data, not a clean result.")
    pains = [x for x in assessed if x.pains]
    unchecked = len(liabs) - len(assessed)

    bits = [f"{len(assessed)} compound(s) checked for structural "
            f"liabilities."]
    if pains:
        bits.append(
            f"**{len(pains)} match a pan-assay interference pattern** "
            f"({', '.join(x.drug for x in pains[:5])}). These are not unsafe "
            f"- they are approved drugs - but activity measured for them "
            f"reads positive across unrelated assays for reasons unconnected "
            f"to any target, so the evidence behind their routes is weaker "
            f"than the tier alone suggests.")
    else:
        bits.append("None matched a pan-assay interference pattern, so no "
                    "route here is weakened on that ground.")
    if unchecked:
        bits.append(f"{unchecked} had no structure and were not checked.")
    return " ".join(bits)
