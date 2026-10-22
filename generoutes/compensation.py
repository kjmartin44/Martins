"""You may not need to raise your gene at all.

The question this asks
----------------------
Every other module here takes the target gene as fixed and looks for a way to
raise it. But the clinical goal is almost never "more of this transcript" —
it is "restore this function". If a paralog does the same job and is easier
to reach, raising that instead solves the same problem.

This is not a theoretical move. It is how several of the most successful
programmes in the field actually worked:

    SMA         SMN1 is deleted. Nobody raises SMN1. Risdiplam and
                nusinersen redirect splicing of SMN2, the paralog, which
                normally makes mostly non-functional protein.
    Sickle cell HBB carries the sickling mutation. Hydroxyurea does not touch
                HBB — it raises HBG1/HBG2, the fetal paralog, which was
                silenced after birth.
    DMD         Dystrophin is absent. Utrophin is its autosomal paralog,
                expressed in fetal muscle and downregulated in adults; raising
                it has been a serious programme for twenty years.

In each case the target gene was undruggable and the paralog was not. A tool
that only ever asks about the gene you named cannot find any of them.

Why this is a different search
------------------------------
`paralogs.py` also talks about paralogs, but asks the opposite question: how
many relatives share a *receptor's binding pocket*, as a cross-reactivity
risk. This asks whether a relative shares your gene's *function*, as an
opportunity. Same word, inverted use.

What decides whether it works
-----------------------------
Three things, and all three have to hold:

**Functional redundancy.** The paralog must do the job, not merely resemble
the sequence. Utrophin substitutes for dystrophin structurally; many close
paralogs do not substitute at all.

**Silencing that can be reversed.** The useful cases are paralogs that are
*switched off* rather than absent — a fetal gene after the developmental
switch, or a second copy that splices unproductively. Reactivating something
the cell already knows how to make is far easier than inventing expression.

**A tractable route of its own.** The paralog has to be more druggable than
the original, or nothing is gained.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

log = logging.getLogger(__name__)


@dataclass
class Redundancy:
    gene: str
    paralog: str
    relationship: str
    why_it_substitutes: str
    silencing: str               # why the paralog is low, if it is
    precedent: str
    status: str                  # "approved" | "clinical" | "preclinical"
    caveats: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return dict(self.__dict__)


# --- Curated. Short and certainly incomplete; see `summarise` on that. ------
PAIRS: tuple = (
    Redundancy(
        gene="SMN1", paralog="SMN2",
        relationship="near-identical duplicate, 5q13",
        why_it_substitutes=(
            "identical protein when exon 7 is included; the paralog differs "
            "by a single silent C>T that disrupts an exonic splice enhancer"),
        silencing="not silenced — mis-spliced, skipping exon 7 ~90% of the time",
        precedent="risdiplam, nusinersen (both FDA-approved) redirect SMN2 "
                  "splicing rather than touching SMN1",
        status="approved",
        caveats=["copy number of the paralog varies between patients and "
                 "predicts both severity and response"]),
    Redundancy(
        gene="HBB", paralog="HBG1/HBG2",
        relationship="fetal paralog in the beta-globin locus",
        why_it_substitutes=(
            "fetal haemoglobin carries oxygen and, in sickle cell, directly "
            "interferes with HbS polymerisation"),
        silencing="developmentally switched off after birth, largely by BCL11A",
        precedent="hydroxyurea (FDA-approved) raises HbF; BCL11A-directed "
                  "gene therapies (exa-cel) target the switch itself",
        status="approved",
        caveats=["hydroxyurea's HbF mechanism is still incompletely defined",
                 "the response is highly variable between patients"]),
    Redundancy(
        gene="DMD", paralog="UTRN",
        relationship="autosomal paralog of dystrophin",
        why_it_substitutes=(
            "binds the same dystroglycan complex and stabilises the "
            "sarcolemma; overexpression rescues the mdx mouse"),
        silencing="expressed in fetal and regenerating muscle, downregulated "
                  "in mature fibres",
        precedent="ezutromid reached clinical trials and failed on efficacy; "
                  "the approach remains active",
        status="clinical",
        caveats=["no approved utrophin upregulator",
                 "utrophin is not a perfect structural substitute, and the "
                 "level needed for benefit is unsettled"]),
    Redundancy(
        gene="SERPINA1", paralog="(none functional)",
        relationship="no redundant paralog",
        why_it_substitutes="—",
        silencing="—",
        precedent="augmentation therapy supplies the protein instead",
        status="approved",
        caveats=["listed as the contrasting case: where no paralog exists, "
                 "protein replacement rather than paralog induction is the "
                 "route"]),
)

_BY_GENE = {r.gene.upper(): r for r in PAIRS}


def known_redundancy(gene: str) -> Redundancy | None:
    return _BY_GENE.get((gene or "").strip().upper())


def find_paralogs(client, gene: str, *, species: str = "human",
                  min_identity: float = 30.0, limit: int = 12) -> list:
    """Within-species paralogs from Ensembl, as candidates for redundancy.

    Sequence paralogy is necessary and nowhere near sufficient: it says two
    genes share ancestry, not that one does the other's job. Every hit needs
    a functional argument before it means anything, which is why the curated
    table above carries `why_it_substitutes` as a required field.
    """
    url = (f"https://rest.ensembl.org/homology/symbol/{species}/{gene}"
           f"?type=paralogues;format=condensed;content-type=application/json")
    # A failure propagates. It used to return [], which reads as "Ensembl
    # knows no paralogs"; the caller says which of the two happened.
    res = client.get(url)

    homologies = []
    for block in (res or {}).get("data", []):
        homologies.extend(block.get("homologies", []))

    out = []
    for h in homologies:
        target = h.get("target") or {}
        ident = target.get("perc_id")
        if ident is None or ident < min_identity:
            continue
        out.append({
            "paralog": target.get("id"),
            "identity": round(float(ident), 1),
            "type": h.get("type"),
            "note": (f"{round(float(ident),1)}% identity. Shared ancestry "
                     f"only — whether it does this gene's job is a separate "
                     f"question this cannot answer."),
        })
    out.sort(key=lambda r: -r["identity"])
    return out[:limit]


def assess(gene: str, *, goal: str = "up", client=None) -> dict:
    """Is raising a relative a better route than raising this gene?"""
    if goal != "up":
        return {}

    known = known_redundancy(gene)
    result = {
        "gene": gene.upper(),
        "known_redundancy": known.to_dict() if known else None,
        "paralogs": [],
        "notes": [],
    }

    if known and known.paralog != "(none functional)":
        result["notes"].append(
            f"**{known.gene} has a documented functional paralog: "
            f"{known.paralog}.** {known.why_it_substitutes.capitalize()}. "
            f"{known.precedent}.")
        result["notes"].append(
            f"The paralog is low because it is {known.silencing} — which is "
            f"the useful case: reactivating something the cell already knows "
            f"how to make is far easier than inventing expression.")
        result["recommended"] = True
    elif known:
        result["notes"].append(
            f"No functional paralog for {known.gene}. {known.precedent}.")
        result["recommended"] = False
    else:
        result["recommended"] = None
        result["notes"].append(
            f"{gene.upper()} is not in the curated redundancy table, which is "
            f"short and certainly incomplete. Absence here is not evidence "
            f"that no paralog substitutes — it means nobody has written one "
            f"down in this file.")

    if client is not None:
        try:
            result["paralogs"] = find_paralogs(client, gene)
        except Exception as exc:                            # noqa: BLE001
            result["paralogs"] = []
            log.warning("Ensembl paralog lookup for %s failed: %s", gene,
                        str(exc)[:120])
            result["notes"].append(
                f"The Ensembl paralog lookup failed ({str(exc)[:100]}), so "
                f"the absence of sequence paralogs here is not a finding.")
        if result["paralogs"]:
            result["notes"].append(
                f"{len(result['paralogs'])} sequence paralog(s) found in "
                f"Ensembl. Sequence paralogy says these share ancestry, NOT "
                f"that any of them does this gene's job — each needs a "
                f"functional argument before it is a route.")

    return result


def summarise(result: dict) -> str:
    if not result:
        return ""
    if result.get("recommended"):
        k = result["known_redundancy"]
        return (f"**Consider raising {k['paralog']} instead of "
                f"{result['gene']}.** This is how the most successful "
                f"programmes in this space actually worked — the target gene "
                f"was undruggable and the paralog was not.")
    if result.get("recommended") is False:
        return (f"No functional paralog for {result['gene']}: raising a "
                f"relative is not available here, and supplying the protein "
                f"may be the better question.")
    return ("Whether a paralog could substitute for this gene is unanswered. "
            "It is worth asking explicitly, because the question 'restore "
            "this function' has more answers than 'raise this transcript'.")
