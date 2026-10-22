"""Decide which of the motif-matched TFs are plausible *activators*.

The request was for "transcriptional activator receptors", which collapses two
distinct categories. Both are annotated here, separately, because they behave
very differently when you try to drug them:

  activator      GO:0001228 / GO:0001216 - the factor has annotated
                 transcription-activating activity. Most of these are not
                 directly ligand-controllable.

  nuclear receptor (GO:0004879 / GO:0098531, plus a curated NR list) - a
                 ligand-activated TF. This is the subset with real small-
                 molecule agonists in the clinic, and almost certainly what
                 "activator receptor" is reaching for.

Caveats that the code cannot fix for you:
  * Activator/repressor is context-dependent. NR3C1 activates at GREs and
    represses at nGREs and via tethering to AP-1/NF-kB on the same cell. Many
    factors carry both GO terms; they are reported as "dual".
  * GO annotations are incomplete. Absence of GO:0001228 is weak evidence.
"""

from __future__ import annotations

import logging

from .http import Client, HttpError

log = logging.getLogger(__name__)

MYGENE = "https://mygene.info/v3"

ACTIVATOR_GO = {
    "GO:0001228",  # DNA-binding transcription activator activity, RNA pol II-specific
    "GO:0001216",  # DNA-binding transcription activator activity
    "GO:0003713",  # transcription coactivator activity
}
REPRESSOR_GO = {
    "GO:0001227",  # DNA-binding transcription repressor activity, RNA pol II-specific
    "GO:0001217",  # DNA-binding transcription repressor activity
    "GO:0003714",  # transcription corepressor activity
}
LIGAND_ACTIVATED_GO = {
    "GO:0004879",  # nuclear receptor activity
    "GO:0098531",  # ligand-activated transcription factor activity
    "GO:0003707",  # nuclear steroid receptor activity (legacy)
}

# The 48 human nuclear receptors. Bundled so NR classification never depends on
# a network call or on GO coverage.
NUCLEAR_RECEPTORS = {
    "NR0B1", "NR0B2",
    "NR1D1", "NR1D2",
    "NR1H2", "NR1H3", "NR1H4",
    "NR1I2", "NR1I3",
    "NR2C1", "NR2C2", "NR2E1", "NR2E3", "NR2F1", "NR2F2", "NR2F6",
    "NR3C1", "NR3C2",
    "NR4A1", "NR4A2", "NR4A3",
    "NR5A1", "NR5A2",
    "NR6A1",
    "AR", "ESR1", "ESR2", "PGR",
    "ESRRA", "ESRRB", "ESRRG",
    "HNF4A", "HNF4G",
    "PPARA", "PPARD", "PPARG",
    "RARA", "RARB", "RARG",
    "RXRA", "RXRB", "RXRG",
    "RORA", "RORB", "RORC",
    "THRA", "THRB", "VDR",
}


class TFAnnotation(dict):
    """dict with convenience accessors; keeps JSON serialisation trivial."""

    @property
    def role(self) -> str:
        return self.get("role", "unknown")

    @property
    def is_activator(self) -> bool:
        return self.role in ("activator", "dual")

    @property
    def is_nuclear_receptor(self) -> bool:
        return bool(self.get("nuclear_receptor"))


def _blank(symbol: str) -> TFAnnotation:
    return TFAnnotation(
        symbol=symbol, name=None, entrez=None, uniprot=None,
        role="unknown", go_activator=[], go_repressor=[],
        nuclear_receptor=symbol in NUCLEAR_RECEPTORS,
        ligand_activated=symbol in NUCLEAR_RECEPTORS,
        source="curated-NR-list" if symbol in NUCLEAR_RECEPTORS else "none",
    )


def _go_ids(hit: dict) -> list[dict]:
    go = hit.get("go") or {}
    mf = go.get("MF") or []
    if isinstance(mf, dict):
        mf = [mf]
    bp = go.get("BP") or []
    if isinstance(bp, dict):
        bp = [bp]
    return list(mf) + list(bp)


def _from_hit(symbol: str, hit: dict) -> TFAnnotation:
    ann = _blank(symbol)
    terms = _go_ids(hit)
    ids = {t.get("id") for t in terms if isinstance(t, dict)}

    act = sorted(ids & ACTIVATOR_GO)
    rep = sorted(ids & REPRESSOR_GO)
    lig = sorted(ids & LIGAND_ACTIVATED_GO)

    if act and rep:
        role = "dual"
    elif act:
        role = "activator"
    elif rep:
        role = "repressor"
    else:
        role = "unknown"

    uni = hit.get("uniprot") or {}
    ann.update(
        name=hit.get("name"),
        entrez=hit.get("entrezgene"),
        uniprot=uni.get("Swiss-Prot") if isinstance(uni, dict) else None,
        role=role,
        go_activator=act,
        go_repressor=rep,
        ligand_activated=bool(lig) or ann["nuclear_receptor"],
        source="mygene/GO",
    )
    return ann


def annotate_tfs(client: Client, symbols: list[str]) -> dict[str, TFAnnotation]:
    """Batch-annotate TF gene symbols. Never raises; degrades to the NR list."""
    symbols = sorted({s.upper() for s in symbols if s})
    out = {s: _blank(s) for s in symbols}

    # Local files first. On a network that blocks MyGene the role falls back
    # to "unknown" for every factor, and the role is what decides which way
    # a drug has to push - abundance.py skips any factor without one. HGNC
    # and a GOA file answer the same question from disk, using the same GO
    # term sets, so the two paths cannot disagree.
    try:
        from . import localannot as _la
        if any(_la.available()[k] for k in ("hgnc", "goa")):
            local = _la.annotate_local(symbols)
            resolved = 0
            for sym, ann in local.items():
                if ann.get("role") != "unknown" or ann.get("uniprot"):
                    out[sym] = ann
                    resolved += 1
            if resolved == len(symbols):
                log.info("annotated %d factor(s) from local HGNC/GOA; "
                         "MyGene not called", resolved)
                return out
            if resolved:
                log.info("annotated %d of %d factor(s) locally; asking "
                         "MyGene for the rest", resolved, len(symbols))
                symbols = [s for s in symbols if s not in local
                           or (local[s].get("role") == "unknown"
                               and not local[s].get("uniprot"))]
    except Exception as exc:                                # noqa: BLE001
        log.debug("local annotation unavailable: %s", str(exc)[:80])

    for i in range(0, len(symbols), 100):
        chunk = symbols[i:i + 100]
        try:
            res = client.post(f"{MYGENE}/query", {
                "q": chunk,
                "scopes": "symbol,alias",
                "species": "human",
                "fields": "symbol,name,entrezgene,uniprot.Swiss-Prot,go.MF,go.BP",
            })
        except HttpError as exc:
            log.warning("MyGene lookup failed (%s); falling back to curated NR list", exc)
            continue

        if isinstance(res, dict):
            res = [res]
        for hit in res:
            if hit.get("notfound"):
                continue
            sym = (hit.get("query") or hit.get("symbol") or "").upper()
            if sym not in out:
                continue
            # Several records can come back per query; keep the most informative.
            cand = _from_hit(sym, hit)
            if out[sym]["role"] == "unknown" or cand["role"] != "unknown":
                out[sym] = cand

    return out
