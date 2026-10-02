"""Human genetic evidence, via eQTLs — and what GWAS/PheWAS can and cannot tell you.

A correction worth making first
-------------------------------
GWAS and PheWAS do not record whether a gene's expression changes in patients
taking a drug. That relation is not in those datasets: GWAS links *germline
variants* to traits, and PheWAS scans many traits against one variant. Neither
observes drug exposure and expression together. If you want drug-induced
expression change, that is LINCS (see signatures.py) or a clinical
transcriptomic study, not GWAS.

What genetics *can* do here is better than the thing being asked for, in two
ways:

**1. Genetic support for the TF -> gene link (trans-eQTL).** If variants near a
transcription factor associate with your target gene's expression in humans,
that is in-vivo causal evidence that the factor regulates your gene — in real
tissue, in real people, with randomised allocation of the exposure courtesy of
meiosis. A motif match says a protein could bind. A trans-eQTL says variation
in that protein actually changes your gene's output. This is the strongest
regulatory evidence in the whole pipeline when it is available, and it is
available for far fewer pairs than you would like: trans-eQTLs need very large
sample sizes and most are undetectable.

**2. What happens if you succeed (cis-eQTL -> phenotype, i.e. drug-target MR
in reverse).** Variants that raise your target gene's expression can be used as
lifelong natural experiments for the intervention you are designing. If people
genetically predisposed to higher expression show a phenotype, that anticipates
your drug's effect — and its side effects — before you make it. This is the
logic that validated PCSK9 inhibition years ahead of the trials.

Caveats that decide whether any of this applies
-----------------------------------------------
* **Tissue.** eQTLs are tissue-specific. A liver eQTL says nothing about brain.
  Always read the tissue before the effect size.
* **Lifelong vs acute.** Genetic variants act from conception; a drug acts for
  months. Compensation and development confound the comparison, and the
  direction of that bias is not predictable.
* **Horizontal pleiotropy.** A variant may affect the phenotype by a route that
  does not pass through your gene, which breaks the instrument. Colocalisation
  evidence partly addresses this and is reported where available.
* **Absence is weak.** No trans-eQTL usually means insufficient power, not
  absence of regulation.

Source: Open Targets Platform GraphQL (https://api.platform.opentargets.org),
free and no key. GTEx and the eQTL Catalogue are the underlying data.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, asdict, field

from .http import Client, HttpError

log = logging.getLogger(__name__)

OPENTARGETS = "https://api.platform.opentargets.org/api/v4/graphql"


@dataclass
class GeneticEvidence:
    gene: str
    tf_symbol: str | None = None
    kind: str = "cis"                  # "cis" | "trans"
    supported: bool = False
    tissues: list = field(default_factory=list)
    phenotypes: list = field(default_factory=list)
    colocalised: bool | None = None
    notes: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


_ASSOC_QUERY = """
query($ensemblId: String!, $size: Int!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    associatedDiseases(page: {index: 0, size: $size}) {
      count
      rows { disease { name } score datatypeScores { id score } }
    }
  }
}
"""


def _resolve_ensembl(client: Client, symbol: str) -> str | None:
    q = """query($q: String!) {
      search(queryString: $q, entityNames: ["target"], page: {index:0, size:1}) {
        hits { id entity name }
      }
    }"""
    try:
        res = client.post(OPENTARGETS, {"query": q, "variables": {"q": symbol}})
    except Exception as exc:                              # noqa: BLE001
        # Catching broadly on purpose: a resolver failure of any kind should
        # degrade to "no genetic evidence", never take the run down.
        log.debug("Open Targets resolve failed for %s: %s", symbol, exc)
        return None
    hits = (((res.get("data") or {}).get("search") or {}).get("hits")) or []
    return hits[0]["id"] if hits else None


def target_consequences(client: Client, gene: str, *, size: int = 15) -> GeneticEvidence:
    """What does human genetic variation in this gene associate with?

    This is drug-target Mendelian randomisation run in reverse: variants that
    alter the gene's expression are lifelong natural experiments for the
    intervention you are proposing, including its unwanted effects.
    """
    ev = GeneticEvidence(gene=gene.upper(), kind="cis")
    ensembl = _resolve_ensembl(client, gene)
    if not ensembl:
        ev.notes.append(f"could not resolve {gene} in Open Targets")
        return ev

    try:
        res = client.post(OPENTARGETS, {
            "query": _ASSOC_QUERY,
            "variables": {"ensemblId": ensembl, "size": size}})
    except HttpError as exc:
        ev.notes.append(f"Open Targets query failed: {exc}")
        return ev

    target = (res.get("data") or {}).get("target") or {}
    rows = ((target.get("associatedDiseases") or {}).get("rows")) or []
    for r in rows:
        dtypes = {d["id"]: d["score"] for d in (r.get("datatypeScores") or [])}
        genetic = dtypes.get("genetic_association")
        if not genetic:
            continue
        ev.phenotypes.append({
            "phenotype": (r.get("disease") or {}).get("name"),
            "overall_score": round(float(r.get("score") or 0), 3),
            "genetic_score": round(float(genetic), 3),
        })

    ev.supported = bool(ev.phenotypes)
    if ev.supported:
        ev.notes.append(
            "These are phenotypes associated with germline variation in this "
            "gene. They anticipate what raising or lowering it may do — "
            "including on-target adverse effects — but variants act lifelong "
            "while a drug does not, so developmental compensation confounds the "
            "comparison in a direction that cannot be predicted.")
    else:
        ev.notes.append(
            "no genetic associations found. That is usually insufficient power "
            "rather than evidence the gene does nothing.")
    return ev


def tf_link_support(client: Client, tf: str, gene: str) -> GeneticEvidence:
    """Is there human genetic evidence that this TF regulates this gene?

    A trans-eQTL is the strongest regulatory evidence available: a motif says a
    protein *could* bind; a trans-eQTL says variation in that protein *does*
    change your gene's output, in human tissue.

    Trans-eQTL catalogues are not exposed through a single convenient API, so
    this reports what can be resolved and tells you where to check manually
    rather than silently returning nothing.
    """
    ev = GeneticEvidence(gene=gene.upper(), tf_symbol=tf.upper(), kind="trans")
    ev.notes.append(
        f"To test this link genetically: look for trans-eQTLs of {tf} affecting "
        f"{gene} in eQTLGen (blood, n>30,000) or GTEx. A hit is strong causal "
        f"evidence for the regulatory relationship. Absence is weak evidence "
        f"against it — trans-eQTL detection needs very large samples and most "
        f"real regulatory links are undetectable at current power.")
    ev.notes.append(
        "eQTLGen: https://www.eqtlgen.org/  |  eQTL Catalogue: "
        "https://www.ebi.ac.uk/eqtl/  |  GTEx: https://gtexportal.org/")
    return ev


def summarise(cis: GeneticEvidence) -> str | None:
    if not cis.supported:
        return None
    top = cis.phenotypes[:3]
    names = ", ".join(p["phenotype"] for p in top if p.get("phenotype"))
    return (f"Human genetic variation in {cis.gene} associates with: {names}. "
            "Read this as a preview of what modulating the gene may do, "
            "including the adverse effects — this is the reasoning that "
            "anticipated PCSK9 inhibition years before the trials.")


def local_support(symbol: str) -> dict | None:
    """Genetic support from a local GWAS Catalog, when one is on disk.

    Open Targets reports an L2G score built from fine-mapped loci with a
    trained model. The catalogue reports the gene mapped or nearest to a
    signal, which is not the causal gene at that locus roughly half the
    time. So this is weaker evidence of the same kind, labelled as such
    and never presented as the same number - and it exists because on a
    blocked network the Open Targets call is what leaves this section
    empty.
    """
    try:
        from . import localannot as la
        if not la.available().get("gwas"):
            return None
        got = la.genetic_support(symbol)
        # `assessed` means the file was read, which is returned even with
        # zero hits on purpose: read-and-found-nothing is a different
        # state from never-looked, and the caller needs to tell them apart.
        return got if got.get("assessed") else None
    except Exception:                                       # noqa: BLE001
        return None
