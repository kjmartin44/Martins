"""Gene annotation and genetic evidence from files, not web APIs.

What this replaces
------------------
`annotate.py` asks MyGene for a factor's name, UniProt accession and GO
terms, and derives the activator/repressor role from those terms. `genetics`
asks Open Targets whether human genetics supports the gene. Both are web
APIs, and on a network that blocks them the role falls back to "unknown" for
every factor - which matters more than it sounds, because the role decides
which direction a drug has to push and `abundance.py` skips any factor
without one.

Three files cover it, and all three are plain tables:

  HGNC complete set    symbol, name, Entrez id, UniProt accession
  GOA human (.gaf)     symbol to GO term, with the aspect
  GWAS Catalog         gene to trait associations with p-values

The role logic is NOT reimplemented here. It imports ACTIVATOR_GO,
REPRESSOR_GO and LIGAND_ACTIVATED_GO from `annotate`, so an offline run and
an online run cannot disagree about whether a factor is an activator. That
was the thing most worth getting right: two code paths deciding direction by
different rules would be worse than having only one of them.

What the offline path loses
---------------------------
MyGene resolves aliases and withdrawn symbols; HGNC carries its own alias
and previous-symbol columns, which are read, but the coverage is not
identical. A symbol absent from HGNC gets a blank annotation rather than a
guessed one.

Open Targets genetics is an L2G score built from fine-mapped GWAS loci with
a trained model. The GWAS Catalog is the underlying association list and no
more: a reported gene for a trait at a p-value. That is weaker, and it is
labelled as a different kind of evidence rather than presented as the same
number. A catalogue association says somebody found a signal near this gene;
it does not say the gene is the causal one at that locus.
"""

from __future__ import annotations

import csv
import gzip
import io
import logging
from pathlib import Path

from .annotate import (ACTIVATOR_GO, LIGAND_ACTIVATED_GO, NUCLEAR_RECEPTORS,
                       REPRESSOR_GO, TFAnnotation)

log = logging.getLogger(__name__)

DEFAULT_HGNC_PATHS = ("data/hgnc_complete_set.txt",
                      "data/hgnc_complete_set.tsv", "data/hgnc.tsv")
DEFAULT_GOA_PATHS = ("data/goa_human.gaf.gz", "data/goa_human.gaf",
                     "data/goa_human.gz")
DEFAULT_GWAS_PATHS = ("data/gwas_catalog_associations.tsv",
                      "data/gwas_catalog.tsv")

SOURCE = """\
HGNC   https://www.genenames.org/download/archive/  (hgnc_complete_set.txt,
       about 15 MB) -> data/hgnc_complete_set.txt
GOA    https://current.geneontology.org/annotations/goa_human.gaf.gz
       (about 10 MB) -> data/goa_human.gaf.gz
GWAS   https://www.ebi.ac.uk/gwas/docs/file-downloads  (all associations,
       about 100 MB) -> data/gwas_catalog_associations.tsv

All optional. Each one that is present replaces a web call; each one that is
absent leaves that field unresolved rather than guessed.
"""

_HGNC: dict = {}
_GOA: dict = {}
_GWAS: dict = {}
_LOADED = {"hgnc": None, "goa": None, "gwas": None}


def _resolve(explicit, candidates):
    for c in ([explicit] if explicit else list(candidates)):
        if c and Path(c).exists():
            return str(c)
    return None


def _open(path):
    """Text handle for a plain or gzipped file, sniffed from the bytes."""
    p = Path(path)
    with open(p, "rb") as probe:
        gz = probe.read(2) == b"\x1f\x8b"
    if gz:
        return io.TextIOWrapper(gzip.open(p, "rb"), encoding="utf-8",
                                errors="replace")
    return open(p, "r", encoding="utf-8", errors="replace")


# --------------------------------------------------------------- HGNC
def load_hgnc(path=None) -> int:
    """symbol -> {name, uniprot, entrez}, plus aliases and old symbols.

    Idempotent: the file is ~15 MB and this is called per factor.
    """
    global _HGNC
    if _HGNC:
        return len(_HGNC)
    p = _resolve(path, DEFAULT_HGNC_PATHS)
    if not p:
        return 0

    with _open(p) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = {c.lower(): c for c in (reader.fieldnames or [])}

        def col(*names):
            for n in names:
                if n in cols:
                    return cols[n]
            return None

        c_sym = col("symbol", "approved symbol")
        c_name = col("name", "approved name")
        c_uni = col("uniprot_ids", "uniprot ids", "uniprot")
        c_ent = col("entrez_id", "ncbi gene id", "entrez")
        c_alias = col("alias_symbol", "alias symbols")
        c_prev = col("prev_symbol", "previous symbols")
        if not c_sym:
            log.warning("%s has no symbol column; HGNC not loaded", p)
            return 0

        for row in reader:
            sym = (row.get(c_sym) or "").strip().upper()
            if not sym:
                continue
            uni = (row.get(c_uni) or "").strip() if c_uni else ""
            rec = {
                "symbol": sym,
                "name": (row.get(c_name) or "").strip() or None
                        if c_name else None,
                # HGNC can list several accessions; the first is canonical.
                "uniprot": (uni.split("|")[0].strip() or None) if uni
                           else None,
                "entrez": (row.get(c_ent) or "").strip() or None
                          if c_ent else None,
            }
            _HGNC[sym] = rec
            # Aliases and retired symbols point at the record, so a motif
            # table using an old name still resolves.
            for c in (c_alias, c_prev):
                if not c:
                    continue
                for a in (row.get(c) or "").split("|"):
                    a = a.strip().upper()
                    if a and a not in _HGNC:
                        _HGNC[a] = rec

    _LOADED["hgnc"] = p
    log.info("HGNC: %d symbol(s) including aliases from %s",
             len(_HGNC), Path(p).name)
    return len(_HGNC)


# ---------------------------------------------------------------- GOA
def load_goa(path=None) -> int:
    """symbol -> set of GO ids, from a GAF. Idempotent."""
    global _GOA
    if _GOA:
        return len(_GOA)
    p = _resolve(path, DEFAULT_GOA_PATHS)
    if not p:
        return 0

    n = 0
    with _open(p) as fh:
        for line in fh:
            if not line or line[0] == "!":
                continue
            f = line.rstrip("\n").split("\t")
            # GAF 2.2: 2=symbol, 3=qualifier, 4=GO id, 6=evidence, 8=aspect
            if len(f) < 9:
                continue
            qualifier = f[3]
            if "NOT" in qualifier:
                # An explicit NOT annotation is evidence AGAINST the term,
                # so counting it would invert the role it contributes to.
                continue
            sym = f[2].strip().upper()
            go = f[4].strip()
            if sym and go.startswith("GO:"):
                _GOA.setdefault(sym, set()).add(go)
                n += 1

    _LOADED["goa"] = p
    log.info("GOA: %d annotation(s) over %d symbol(s) from %s",
             n, len(_GOA), Path(p).name)
    return len(_GOA)


def annotate_local(symbols) -> dict:
    """`{symbol: TFAnnotation}` from HGNC and GOA, same shape as MyGene.

    The role is decided by the GO sets imported from `annotate`, so this
    path and the network path cannot disagree about direction.
    """
    load_hgnc()
    load_goa()
    out: dict = {}
    for sym in symbols:
        s = (sym or "").strip().upper()
        if not s:
            continue
        rec = _HGNC.get(s) or {}
        # GOA is keyed on the approved symbol, so an alias found in HGNC has
        # to be followed to its canonical name or it loses its GO terms -
        # and with them its role, which is what decides direction.
        canonical = (rec.get("symbol") or s).upper()
        ids = _GOA.get(s) or _GOA.get(canonical) or set()
        act = sorted(ids & ACTIVATOR_GO)
        rep = sorted(ids & REPRESSOR_GO)
        lig = sorted(ids & LIGAND_ACTIVATED_GO)
        role = ("dual" if act and rep else "activator" if act
                else "repressor" if rep else "unknown")

        is_nr = s in NUCLEAR_RECEPTORS or canonical in NUCLEAR_RECEPTORS
        ann = TFAnnotation(
            symbol=sym, name=rec.get("name"), entrez=rec.get("entrez"),
            uniprot=rec.get("uniprot"), role=role,
            nuclear_receptor=is_nr,
            go_activator=act, go_repressor=rep,
            ligand_activated=bool(lig) or is_nr,
            source="local HGNC/GOA",
        )
        if not ids:
            ann["note"] = (
                "no GO annotation on disk for this symbol, so the role is "
                "unresolved rather than absent")
        out[sym] = ann
    return out


def uniprot_for(symbol: str) -> str | None:
    """Gene symbol to UniProt accession from the local HGNC table.

    One accessor, because `TFAnnotation` subclasses dict and three callers
    read the accession with `getattr`, which returns None for every gene
    and fails silently - the docking CLI refused every `--gene` and the
    offline potency lookup never ran.

    Aliases and retired symbols resolve, since `load_hgnc` indexes them.
    Returns None when no HGNC table is on disk, which means unknown.
    """
    if not load_hgnc():
        return None
    s = (symbol or "").strip().upper()
    return (_HGNC.get(s) or {}).get("uniprot")


def available() -> dict:
    """Which local annotation files are present."""
    return {
        "hgnc": _resolve(None, DEFAULT_HGNC_PATHS),
        "goa": _resolve(None, DEFAULT_GOA_PATHS),
        "gwas": _resolve(None, DEFAULT_GWAS_PATHS),
    }


# --------------------------------------------------------------- GWAS
def load_gwas(path=None) -> int:
    """gene symbol -> [{trait, p_value, pubmed}] from the GWAS Catalog.

    Reported genes only. The catalogue names the gene nearest or mapped to
    a signal, which is not the same as the causal gene at that locus, and
    `genetic_support` says so rather than presenting it as Open Targets'
    L2G score.
    """
    global _GWAS
    if _GWAS:
        return len(_GWAS)
    p = _resolve(path, DEFAULT_GWAS_PATHS)
    if not p:
        return 0

    with _open(p) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = {c.lower(): c for c in (reader.fieldnames or [])}

        def col(*names):
            for n in names:
                if n in cols:
                    return cols[n]
            return None

        c_gene = col("mapped_gene", "reported gene(s)", "mapped gene")
        c_trait = col("disease/trait", "mapped_trait", "trait")
        c_p = col("p-value", "p_value", "pvalue")
        c_pub = col("pubmedid", "pubmed_id")
        if not (c_gene and c_trait):
            log.warning("%s lacks a gene or trait column; GWAS not loaded", p)
            return 0

        for row in reader:
            raw = (row.get(c_gene) or "").strip()
            if not raw:
                continue
            try:
                pval = float((row.get(c_p) or "").strip()) if c_p else None
            except ValueError:
                pval = None
            trait = (row.get(c_trait) or "").strip()
            for g in raw.replace(";", ",").replace(" - ", ",").split(","):
                g = g.strip().upper()
                if not g or g == "NR":
                    continue
                _GWAS.setdefault(g, []).append({
                    "trait": trait, "p_value": pval,
                    "pubmed": (row.get(c_pub) or "").strip()
                              if c_pub else None})

    _LOADED["gwas"] = p
    log.info("GWAS Catalog: %d gene(s) with associations from %s",
             len(_GWAS), Path(p).name)
    return len(_GWAS)


# Genome-wide significance. An association above this is not a finding.
GWAS_SIGNIFICANT = 5e-8


def genetic_support(symbol: str, path=None) -> dict:
    """Catalogue associations for one gene, with what they do not establish."""
    load_gwas(path)
    s = (symbol or "").strip().upper()
    if not _GWAS:
        return {"assessed": False, "n": 0, "source": None,
                "note": ("no GWAS Catalog on disk, so human genetic support "
                         "was not checked. That is different from there "
                         "being none.")}

    hits = _GWAS.get(s, [])
    sig = [h for h in hits
           if h["p_value"] is not None and h["p_value"] <= GWAS_SIGNIFICANT]
    traits = sorted({h["trait"] for h in sig})[:8]
    out = {
        "assessed": True, "source": "GWAS Catalog (local)",
        "n": len(hits), "n_significant": len(sig), "traits": traits,
        "caveats": [
            "the catalogue names the gene mapped or nearest to a signal, "
            "which is not the same as the causal gene at that locus - the "
            "nearest gene is wrong perhaps half the time",
            "this is an association list, not Open Targets' L2G score: "
            "weaker evidence of the same kind, and labelled as such",
        ],
    }
    if not hits:
        out["note"] = (
            f"{s} has no association in the catalogue. For a gene nobody has "
            f"powered a study for, that is absence of evidence rather than "
            f"evidence of absence.")
    return out


def summarise() -> str:
    avail = available()
    have = [k for k, v in avail.items() if v]
    missing = [k for k, v in avail.items() if not v]
    if not have:
        return ("No local annotation files are present, so gene names, "
                "UniProt accessions, GO-derived roles and genetic support "
                "all need the network.\n" + SOURCE)
    bits = [f"Local annotation in use: {', '.join(have)}."]
    if "goa" in have:
        bits.append(
            "GO terms come from the file, and the activator/repressor "
            "decision uses the same term sets as the online path, so the "
            "two cannot disagree.")
    if missing:
        bits.append(
            f"Still needing the network: {', '.join(missing)}. Each absent "
            f"file leaves its field unresolved rather than guessed.")
    return " ".join(bits)
