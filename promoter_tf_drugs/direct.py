"""Drugs that act on the target gene directly, bypassing transcription factors.

Why this belongs at the top of the report
-----------------------------------------
The TF route in the rest of this package is indirect and lossy: predict a
binding site, hope the factor binds it in your cell type, hope it activates
rather than represses, hope the drug moves the factor enough, and accept that
everything else the factor controls moves too.

A drug that targets your gene's own product or transcript skips all of it. If
one exists, it is almost always the better starting point, and it is the first
thing a reviewer will ask whether you checked.

Two categories are collected here:

1. **Protein-level.** Approved drugs whose ChEMBL/DGIdb target *is* your gene
   product. Direction matters as everywhere else: to raise functional output
   you want agonists, activators, chaperones or stabilisers; to lower it,
   inhibitors, antagonists or degraders.

   Note a real distinction the database cannot make for you: an *enzyme
   activator* raises the activity of protein already present; it does not
   raise protein abundance. If your readout is a western blot, that drug will
   look like a failure even when it works. Check what you are actually
   measuring.

2. **Nucleic-acid and gene-level.** Antisense oligonucleotides, siRNA,
   splice-switching oligos and gene therapies act on the transcript or the gene
   itself with sequence-level specificity - far better specificity than any TF
   intervention can offer.

   The upregulating members of this class are the interesting ones for the
   original question, because there is no general "siRNA for turning genes up":
   nusinersen and risdiplam raise SMN protein by correcting SMN2 splicing, the
   DMD exon-skippers restore a partly functional dystrophin, and the AAV gene
   therapies supply a working copy outright.

The bundled table below is curated by hand and **will go out of date**. Approvals
happen constantly. Treat it as a prompt to go and check, not as an authority,
and verify anything you intend to rely on against FDA/EMA labelling.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, asdict

from .drugs import (ACTIVATING_ACTIONS, INHIBITING_ACTIONS, Drug, chembl_drugs,
                    dgidb_drugs)
from .http import Client

log = logging.getLogger(__name__)

# Modality -> how specific it is to your gene (0..1). Sequence-directed agents
# score highest: they act on your transcript, not on a protein family.
MODALITY_SPECIFICITY = {
    "gene therapy": 1.00,
    "siRNA": 0.98,
    "antisense oligonucleotide": 0.97,
    "splice-switching oligonucleotide": 0.97,
    "monoclonal antibody": 0.85,
    "enzyme replacement": 0.90,
    "small molecule": 0.55,
}


@dataclass
class DirectTherapy:
    gene: str
    drug: str
    modality: str
    direction: str            # "up" | "down"
    mechanism: str
    approval_year: int | None
    approved: bool
    withdrawn: bool
    raises_protein: bool      # True = more protein; False = more activity only
    specificity: float
    source: str

    def to_dict(self) -> dict:
        return asdict(self)


def _t(gene, drug, modality, direction, mech, year, *, raises_protein=True,
       withdrawn=False):
    return DirectTherapy(
        gene=gene, drug=drug, modality=modality, direction=direction,
        mechanism=mech, approval_year=year, approved=True, withdrawn=withdrawn,
        raises_protein=raises_protein,
        specificity=MODALITY_SPECIFICITY.get(modality, 0.6),
        source="curated",
    )


# --- Curated approved gene-directed therapies -------------------------------
# Verify before relying on any row. Approval years are first approval in a
# major jurisdiction and may differ between FDA and EMA.
CURATED: dict[str, list[DirectTherapy]] = {}
for _e in [
    # ---- raising functional protein output (the hard direction) ----
    ("SMN2", "Nusinersen", "splice-switching oligonucleotide", "up",
     "Promotes SMN2 exon 7 inclusion, raising full-length SMN protein", 2016),
    ("SMN2", "Risdiplam", "small molecule", "up",
     "SMN2 pre-mRNA splicing modifier; oral, raises full-length SMN", 2020),
    ("SMN1", "Onasemnogene abeparvovec", "gene therapy", "up",
     "AAV9 delivery of a functional SMN1 copy", 2019),
    ("DMD", "Eteplirsen", "splice-switching oligonucleotide", "up",
     "Exon 51 skipping restores the dystrophin reading frame", 2016),
    ("DMD", "Golodirsen", "splice-switching oligonucleotide", "up",
     "Exon 53 skipping", 2019),
    ("DMD", "Viltolarsen", "splice-switching oligonucleotide", "up",
     "Exon 53 skipping", 2020),
    ("DMD", "Casimersen", "splice-switching oligonucleotide", "up",
     "Exon 45 skipping", 2021),
    ("DMD", "Delandistrogene moxeparvovec", "gene therapy", "up",
     "AAV delivery of a micro-dystrophin transgene", 2023),
    ("RPE65", "Voretigene neparvovec", "gene therapy", "up",
     "AAV2 delivery of functional RPE65", 2017),
    ("F9", "Etranacogene dezaparvovec", "gene therapy", "up",
     "AAV5 delivery of factor IX Padua variant", 2022),
    ("F8", "Valoctocogene roxaparvovec", "gene therapy", "up",
     "AAV5 delivery of B-domain-deleted factor VIII", 2023),
    ("HBB", "Betibeglogene autotemcel", "gene therapy", "up",
     "Lentiviral beta-globin gene addition to autologous HSCs", 2022),
    ("HBG1", "Exagamglogene autotemcel", "gene therapy", "up",
     "CRISPR disruption of the BCL11A erythroid enhancer de-represses HbF", 2023),
    ("HBG1", "Hydroxyurea", "small molecule", "up",
     "Induces fetal haemoglobin; mechanism incompletely defined", 1998),
    ("CFTR", "Ivacaftor", "small molecule", "up",
     "CFTR potentiator: raises channel open probability, not protein level",
     2012, False),
    ("CFTR", "Elexacaftor/tezacaftor/ivacaftor", "small molecule", "up",
     "Corrector plus potentiator: improves CFTR folding, trafficking and gating",
     2019),

    # ---- lowering output ----
    ("TTR", "Patisiran", "siRNA", "down", "Hepatic TTR mRNA knockdown", 2018),
    ("TTR", "Vutrisiran", "siRNA", "down", "Hepatic TTR mRNA knockdown", 2022),
    ("TTR", "Inotersen", "antisense oligonucleotide", "down", "TTR mRNA knockdown", 2018),
    ("TTR", "Eplontersen", "antisense oligonucleotide", "down", "TTR mRNA knockdown", 2023),
    ("PCSK9", "Inclisiran", "siRNA", "down", "Hepatic PCSK9 mRNA knockdown", 2021),
    ("PCSK9", "Evolocumab", "monoclonal antibody", "down",
     "Neutralises circulating PCSK9 protein", 2015, False),
    ("PCSK9", "Alirocumab", "monoclonal antibody", "down",
     "Neutralises circulating PCSK9 protein", 2015, False),
    ("ALAS1", "Givosiran", "siRNA", "down", "Hepatic ALAS1 mRNA knockdown", 2019),
    ("HAO1", "Lumasiran", "siRNA", "down", "Hepatic HAO1 mRNA knockdown", 2020),
    ("LDHA", "Nedosiran", "siRNA", "down", "Hepatic LDHA mRNA knockdown", 2023),
    ("SOD1", "Tofersen", "antisense oligonucleotide", "down",
     "Intrathecal SOD1 mRNA knockdown", 2023),
    ("APOC3", "Volanesorsen", "antisense oligonucleotide", "down",
     "APOC3 mRNA knockdown (EU)", 2019),
]:
    gene, drug, modality, direction, mech, year = _e[:6]
    raises = _e[6] if len(_e) > 6 else True
    CURATED.setdefault(gene, []).append(
        _t(gene, drug, modality, direction, mech, year, raises_protein=raises))

# A couple of notable withdrawals, kept visible rather than deleted.
CURATED.setdefault("APOB", []).append(
    _t("APOB", "Mipomersen", "antisense oligonucleotide", "down",
       "APOB-100 mRNA knockdown; withdrawn over hepatotoxicity", 2013,
       withdrawn=True))


# --------------------------------------------------------------- collection
def curated_for_gene(gene: str, goal: str = "up") -> list[DirectTherapy]:
    want = goal
    return [t for t in CURATED.get(gene.upper(), []) if t.direction == want]


def protein_level_drugs(client: Client, gene: str, uniprot: str | None = None, *,
                        goal: str = "up", approved_only: bool = True,
                        use_dgidb: bool = True) -> list[DirectTherapy]:
    """Approved drugs whose ChEMBL/DGIdb target is this gene's own product."""
    want = ACTIVATING_ACTIONS if goal == "up" else INHIBITING_ACTIONS
    direction = goal

    found: list[Drug] = chembl_drugs(client, gene, uniprot, approved_only)
    if use_dgidb:
        found += dgidb_drugs(client, [gene], approved_only=approved_only)

    out, seen = [], set()
    for d in found:
        action = (d.action_type or "").upper()
        if not any(a in action for a in want):
            continue
        key = d.name.lower()
        if key in seen:
            continue
        seen.add(key)

        # An activator raises activity of existing protein, not its abundance.
        raises_protein = not any(
            w in action for w in ("AGONIST", "ACTIVATOR", "OPENER", "MODULATOR"))
        out.append(DirectTherapy(
            gene=gene.upper(), drug=d.name, modality="small molecule",
            direction=direction, mechanism=d.mechanism or action,
            approval_year=d.first_approval, approved=d.approved,
            withdrawn=d.withdrawn, raises_protein=raises_protein,
            specificity=MODALITY_SPECIFICITY["small molecule"], source=d.source,
        ))
    return out


def find_direct(client: Client, gene: str, uniprot: str | None = None, *,
                goal: str = "up", approved_only: bool = True,
                use_dgidb: bool = True, demo: bool = False) -> list[DirectTherapy]:
    """Everything acting directly on this gene, curated plus database."""
    results = curated_for_gene(gene, goal)
    if not demo:
        try:
            results += protein_level_drugs(
                client, gene, uniprot, goal=goal,
                approved_only=approved_only, use_dgidb=use_dgidb)
        except Exception as exc:                        # noqa: BLE001
            log.warning("direct drug lookup failed for %s: %s", gene, exc)

    seen, dedup = set(), []
    for t in results:
        if t.drug.lower() in seen:
            continue
        seen.add(t.drug.lower())
        dedup.append(t)

    # Sequence-directed first, then protein-raising, then recent.
    dedup.sort(key=lambda t: (t.withdrawn, -t.specificity, not t.raises_protein,
                              -(t.approval_year or 0)))
    return dedup
