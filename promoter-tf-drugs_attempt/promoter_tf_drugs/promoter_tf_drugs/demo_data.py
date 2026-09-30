"""Self-contained fixtures so the pipeline runs with zero network access.

Everything here is clearly labelled DEMO_*. The count matrices are built from
published consensus strings, not from JASPAR's actual site counts, and the drug
table is a small hand-entered subset. Use `--demo` to smoke-test the code and
read the report format; use the live pipeline for anything you intend to act
on.
"""

from __future__ import annotations

import numpy as np

from .annotate import TFAnnotation
from .drugs import Drug
from .motifs import Motif

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "W": "AT", "S": "GC", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


def consensus_to_pfm(consensus: str, n_sites: int = 100, strength: float = 0.88) -> np.ndarray:
    """Build an L x 4 count matrix from an IUPAC consensus string."""
    rows = []
    for ch in consensus.upper():
        allowed = IUPAC[ch]
        row = np.zeros(4)
        for i, b in enumerate("ACGT"):
            if b in allowed:
                row[i] = strength / len(allowed)
            else:
                row[i] = (1 - strength) / (4 - len(allowed))
        if len(allowed) == 4:
            row[:] = 0.25
        rows.append(row * n_sites)
    return np.array(rows)


# consensus, JASPAR-style name  (name drives TF symbol extraction)
DEMO_CONSENSUS = [
    ("DEMO0001.1", "NR3C1",        "AGAACANNNTGTTCT"),   # GRE, inverted repeat 3
    ("DEMO0002.1", "PPARG::RXRA",  "AGGTCAAAGGTCA"),     # DR1
    ("DEMO0003.1", "RARA::RXRA",   "AGGTCANNNNNAGGTCA"), # DR5
    ("DEMO0004.1", "VDR::RXRA",    "AGGTCANNNAGGTCA"),   # DR3
    ("DEMO0005.1", "HNF4A",        "GGGGTCAAAGGTCA"),
    ("DEMO0006.1", "SP1",          "GGGGCGGGGC"),
    ("DEMO0007.1", "CREB1",        "TGACGTCA"),
    ("DEMO0008.1", "NFKB1",        "GGGACTTTCC"),
    ("DEMO0009.1", "STAT3",        "TTCCCGGAA"),
    ("DEMO0010.1", "REST",         "TTCAGCACCACGGACAGCGCC"),  # NRSE, repressor
]


def demo_motifs() -> list[Motif]:
    return [Motif(mid, name, consensus_to_pfm(cons))
            for mid, name, cons in DEMO_CONSENSUS]


# Sites planted at these offsets relative to the TSS, to give the scanner
# something real to find.
PLANTED = [
    ("AGAACAGTGTGTTCT", -850),    # GRE
    ("AGGTCAAAGGTCA", -610),      # PPARG/RXRA DR1
    ("AGGTCAGGCTTAGGTCA", -430),  # RARA/RXRA DR5
    ("GGGGCGGGGC", -95),          # SP1
    ("GGGGCGGGGC", -62),
    ("TGACGTCA", -140),
    ("GGGACTTTCC", -310),
    ("TTCAGCACCACGGACAGCGCC", -720),  # REST/NRSF
    ("TATAAAAG", -30),            # TATA box
]


def demo_promoter_sequence(upstream: int = 1000, downstream: int = 200,
                           gc: float = 0.60, seed: int = 20260929) -> str:
    rng = np.random.default_rng(seed)
    probs = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    seq = list("".join(rng.choice(list("ACGT"), size=upstream + downstream, p=probs)))
    for site, offset in PLANTED:
        idx = upstream + offset          # offset is negative -> index < upstream
        if 0 <= idx <= len(seq) - len(site):
            seq[idx:idx + len(site)] = list(site)
    return "".join(seq)


DEMO_ANNOTATIONS_RAW = {
    # symbol: (role, nuclear_receptor, uniprot, full name)
    "NR3C1":  ("dual",      True,  "P04150", "glucocorticoid receptor"),
    "PPARG":  ("activator", True,  "P37231", "peroxisome proliferator activated receptor gamma"),
    "RXRA":   ("activator", True,  "P19793", "retinoid X receptor alpha"),
    "RARA":   ("activator", True,  "P10276", "retinoic acid receptor alpha"),
    "VDR":    ("activator", True,  "P11473", "vitamin D receptor"),
    "HNF4A":  ("activator", True,  "P41235", "hepatocyte nuclear factor 4 alpha"),
    "SP1":    ("activator", False, "P08047", "Sp1 transcription factor"),
    "CREB1":  ("activator", False, "P16220", "cAMP responsive element binding protein 1"),
    "NFKB1":  ("dual",      False, "P19838", "nuclear factor kappa B subunit 1"),
    "STAT3":  ("activator", False, "P40763", "signal transducer and activator of transcription 3"),
    "REST":   ("repressor", False, "Q13127", "RE1 silencing transcription factor"),
}


def demo_annotations(symbols: list[str]) -> dict[str, TFAnnotation]:
    out = {}
    for s in symbols:
        s = s.upper()
        role, nr, uni, name = DEMO_ANNOTATIONS_RAW.get(s, ("unknown", False, None, None))
        out[s] = TFAnnotation(
            symbol=s, name=name, entrez=None, uniprot=uni, role=role,
            go_activator=["GO:0001228"] if role in ("activator", "dual") else [],
            go_repressor=["GO:0001227"] if role in ("repressor", "dual") else [],
            nuclear_receptor=nr, ligand_activated=nr, source="DEMO",
        )
    return out


def _d(name, sym, action, mech, direction, year, oral=True):
    return Drug(name=name, chembl_id=None, target_symbol=sym, action_type=action,
                mechanism=mech, approved=True, max_phase=4.0, first_approval=year,
                withdrawn=False, oral=oral, direction=direction, source="DEMO")


DEMO_DRUGS = {
    "NR3C1": [
        _d("Dexamethasone", "NR3C1", "AGONIST", "Glucocorticoid receptor agonist", "up", 1958),
        _d("Prednisolone", "NR3C1", "AGONIST", "Glucocorticoid receptor agonist", "up", 1955),
        _d("Hydrocortisone", "NR3C1", "AGONIST", "Glucocorticoid receptor agonist", "up", 1952),
    ],
    "PPARG": [
        _d("Pioglitazone", "PPARG", "AGONIST", "PPAR gamma agonist", "up", 1999),
        _d("Rosiglitazone", "PPARG", "AGONIST", "PPAR gamma agonist", "up", 1999),
    ],
    "RARA": [
        _d("Tretinoin", "RARA", "AGONIST", "Retinoic acid receptor agonist", "up", 1971),
        _d("Alitretinoin", "RARA", "AGONIST", "Retinoid receptor agonist", "up", 1999),
    ],
    "RXRA": [
        _d("Bexarotene", "RXRA", "AGONIST", "Retinoid X receptor agonist", "up", 1999),
    ],
    "VDR": [
        _d("Calcitriol", "VDR", "AGONIST", "Vitamin D receptor agonist", "up", 1978),
        _d("Paricalcitol", "VDR", "AGONIST", "Vitamin D receptor agonist", "up", 1998, oral=False),
    ],
    # REST is a repressor: the useful direction is *down*.
    "REST": [
        _d("Valproic Acid", "REST", "INHIBITOR",
           "Indirect: HDAC inhibition relieves REST/CoREST-mediated silencing",
           "down", 1978),
    ],
    # Blocking pharmacology. Note how much easier this direction is to fill:
    # antagonists exist where agonists do not.
    "_DOWN_": [],
}

# Merge antagonists in. Kept separate above only for readability.
for _sym, _entries in {
    "NR3C1": [("Mifepristone", "ANTAGONIST", "Glucocorticoid receptor antagonist", 1988)],
    "SP1":   [("Plicamycin", "INHIBITOR",
               "Binds GC-rich DNA and displaces Sp1", 1970)],
    "ESR1":  [("Tamoxifen", "ANTAGONIST", "Selective estrogen receptor modulator", 1977),
              ("Fulvestrant", "DEGRADER", "Estrogen receptor degrader", 2002)],
    "AR":    [("Enzalutamide", "ANTAGONIST", "Androgen receptor antagonist", 2012)],
}.items():
    DEMO_DRUGS.setdefault(_sym, []).extend(
        _d(n, _sym, a, m, "down", y) for n, a, m, y in _entries)
DEMO_DRUGS.pop("_DOWN_", None)

# Allosteric entries, so the orthosteric-vs-allosteric comparison has something
# to show offline. Illustrative only.
for _sym, _entries in {
    "PPARG": [("Demo-Pam-1", "POSITIVE ALLOSTERIC MODULATOR",
               "Binds the PPARgamma alternate site; amplifies endogenous ligand",
               "up", 2015)],
    "VDR":   [("Demo-Pam-2", "POSITIVE ALLOSTERIC MODULATOR",
               "Allosteric potentiator of calcitriol signalling", "up", 2018)],
}.items():
    DEMO_DRUGS.setdefault(_sym, []).extend(
        _d(n, _sym, a, m, dr, y) for n, a, m, dr, y in _entries)


def demo_drugs(annotations: dict, goal: str = "up") -> dict[str, list[Drug]]:
    from .drugs import filter_drugs_by_goal
    raw = {s: list(DEMO_DRUGS[s]) for s in annotations if s in DEMO_DRUGS}
    actionable, _unknown = filter_drugs_by_goal(raw, annotations, goal)
    return actionable


# --- demo pharmacology evidence -------------------------------------------
# (median pChEMBL, distinct annotated targets, withdrawn, black-box)
DEMO_EVIDENCE = {
    "Dexamethasone":  (8.0, 2, False, False),
    "Prednisolone":   (7.4, 2, False, False),
    "Hydrocortisone": (7.1, 3, False, False),
    "Pioglitazone":   (6.7, 3, False, True),    # CHF black-box
    "Rosiglitazone":  (7.8, 2, False, True),
    "Tretinoin":      (8.3, 6, False, True),
    "Alitretinoin":   (7.9, 6, False, True),
    "Bexarotene":     (7.7, 3, False, True),
    "Calcitriol":     (9.3, 1, False, False),
    "Paricalcitol":   (8.6, 1, False, False),
    "Valproic Acid":  (4.2, 8, False, True),    # weak, promiscuous
    "Mifepristone":   (8.6, 3, False, True),
    "Plicamycin":     (6.5, 4, True, True),     # withdrawn
    "Tamoxifen":      (8.1, 4, False, True),
    "Fulvestrant":    (8.8, 1, False, False),
    "Enzalutamide":   (7.6, 2, False, True),
}

# Rough physicochemistry, for the --tissue cns demo.
DEMO_PROPS = {
    "Dexamethasone": {"full_mwt": 392.5, "psa": 94.8, "cx_logp": 1.9},
    "Calcitriol":    {"full_mwt": 416.6, "psa": 60.7, "cx_logp": 5.1},
    "Valproic Acid": {"full_mwt": 144.2, "psa": 37.3, "cx_logp": 2.8},
    "Plicamycin":    {"full_mwt": 1085.2, "psa": 331.0, "cx_logp": -0.4},
}


def demo_evidence(drug):
    """evidence_lookup callable for drug_scoring.grade_all in demo mode."""
    from .drug_scoring import DrugEvidence
    pch, n_t, withdrawn, box = DEMO_EVIDENCE.get(drug.name, (None, None, False, False))
    return DrugEvidence(
        pchembl_median=pch, pchembl_n=6 if pch else 0, n_targets=n_t,
        withdrawn=withdrawn or bool(drug.withdrawn), black_box=box,
        properties=DEMO_PROPS.get(drug.name, {}), source="DEMO",
    )
