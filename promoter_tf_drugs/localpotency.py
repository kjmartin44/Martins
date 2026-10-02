"""Measured binding potency from a file, not from the ChEMBL API.

The gap this closes
-------------------
`localdrugs` made the drug library and the known-ligand sets available
offline, and said plainly what it could not supply: potency. The
Repurposing Hub records that a compound acts on a target, not how hard it
binds. So on an offline run the potency axis in `drug_scoring` had nothing
to read, the established-ligand ranking in `screening` lost its order, and
`analogs` had no measured anchor.

Two files can supply it, and both are keyed on UniProt, which is what
`localannot` now resolves from HGNC.

  BindingDB   one TSV, free, Ki/Kd/IC50/EC50 per target-ligand pair with
              UniProt accessions and SMILES. About 1 GB uncompressed.
              https://www.bindingdb.org/bind/chemsearch/marvin/Download.jsp

  ChEMBL      the SQLite release. Authoritative and much larger (tens of
              GB unpacked), and queryable with the standard library's
              sqlite3, so it needs no extra dependency.
              https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/

Either produces the same shape `agonism.fetch_actives` returns -
`{chembl_id, smiles, pchembl}` - so nothing downstream has to know which
one answered.

How a potency becomes a pChEMBL
-------------------------------
pChEMBL is -log10 of the molar affinity, so 1 nM is 9 and 1 uM is 6, and
larger is more potent. BindingDB reports nM, sometimes with a `>` or `<`
prefix meaning the assay could not bracket the value.

**A censored measurement is dropped, not converted.** ">10000 nM" means the
compound was inactive up to that concentration, and recording it as
pChEMBL 5 would turn a failed experiment into a weak binding constant. That
one decision matters more than anything else in this module, because those
rows are common and treating them as data biases every potency downward.

Where several measurements exist for a pair, the median is taken rather
than the best. Taking the best would systematically favour compounds that
have been assayed often, which is the study-effort confound the rest of
this package works to avoid.
"""

from __future__ import annotations

import csv
import gzip
import io
import logging
import math
import statistics
from pathlib import Path

log = logging.getLogger(__name__)

DEFAULT_BINDINGDB_PATHS = (
    "data/BindingDB_All.tsv", "data/bindingdb_all.tsv",
    "data/BindingDB_All.tsv.gz", "data/bindingdb.tsv",
)
DEFAULT_CHEMBL_SQLITE_PATHS = (
    "data/chembl.db", "data/chembl_35.db", "data/chembl_35_sqlite/chembl_35.db",
)

SOURCE = """\
BindingDB  https://www.bindingdb.org/bind/chemsearch/marvin/Download.jsp
           "BindingDB_All.tsv" -> data/BindingDB_All.tsv  (about 1 GB)

ChEMBL     https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/
           the SQLite release -> data/chembl.db  (tens of GB unpacked)

BindingDB is the lighter of the two and enough for a potency axis. Neither
is required: without one, potency stays unmeasured rather than guessed.
"""

# Affinity columns BindingDB uses, in the order to prefer them. Ki and Kd
# are equilibrium constants; IC50 and EC50 depend on assay conditions, so
# they are used only when nothing better exists for the pair.
_AFFINITY_COLS = (
    ("Ki (nM)", "Ki"),
    ("Kd (nM)", "Kd"),
    ("IC50 (nM)", "IC50"),
    ("EC50 (nM)", "EC50"),
)
_UNIPROT_COLS = (
    "UniProt (SwissProt) Primary ID of Target Chain",
    "UniProt (SwissProt) Primary ID of Target Chain 1",
    "UniProt (TrEMBL) Primary ID of Target Chain",
)
_SMILES_COLS = ("Ligand SMILES", "SMILES")
_ID_COLS = ("ChEMBL ID of Ligand", "BindingDB MonomerID",
            "BindingDB Reactant_set_id")


def resolve(explicit=None, candidates=DEFAULT_BINDINGDB_PATHS) -> str | None:
    for c in ([explicit] if explicit else list(candidates)):
        if c and Path(c).exists():
            return str(c)
    return None


def _open(path):
    p = Path(path)
    with open(p, "rb") as probe:
        gz = probe.read(2) == b"\x1f\x8b"
    if gz:
        return io.TextIOWrapper(gzip.open(p, "rb"), encoding="utf-8",
                                errors="replace")
    return open(p, "r", encoding="utf-8", errors="replace")


def to_pchembl(value: str) -> tuple:
    """(pchembl, censored). nM in, -log10(M) out.

    `censored` is True for a `>` or `<` bound, and the caller must drop
    those: ">10000 nM" is a compound that did nothing up to 10 uM, and
    converting it to pChEMBL 5 would record a failed experiment as a weak
    binding constant.
    """
    v = (value or "").strip()
    if not v:
        return None, False
    censored = v[0] in "><~"
    v = v.lstrip("><~= ").strip()
    try:
        nm = float(v)
    except ValueError:
        return None, censored
    if nm <= 0:
        return None, censored
    return round(-math.log10(nm * 1e-9), 2), censored


def load_bindingdb(uniprots, path=None, *, min_pchembl: float = 5.0) -> dict:
    """`{uniprot: [{chembl_id, smiles, pchembl, measure, n}]}` in one pass.

    The file is about 1 GB, so every target wanted in a run is extracted
    together rather than re-reading it per target.
    """
    want = {u.strip().upper() for u in uniprots if u}
    if not want:
        return {}
    p = resolve(path)
    if not p:
        log.info("no BindingDB file on disk; potency stays unmeasured. %s",
                 SOURCE.splitlines()[0])
        return {}

    # {(uniprot, smiles): {measure: [pchembl, ...]}}
    pairs: dict = {}
    n_rows = n_censored = 0

    with _open(p) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = reader.fieldnames or []
        uni_c = next((c for c in _UNIPROT_COLS if c in cols), None)
        smi_c = next((c for c in _SMILES_COLS if c in cols), None)
        id_c = next((c for c in _ID_COLS if c in cols), None)
        aff = [(c, label) for c, label in _AFFINITY_COLS if c in cols]
        if not (uni_c and smi_c and aff):
            log.warning("%s does not look like BindingDB (missing UniProt, "
                        "SMILES or affinity columns); potency not loaded", p)
            return {}

        for row in reader:
            n_rows += 1
            uni = (row.get(uni_c) or "").strip().upper()
            if uni not in want:
                continue
            smi = (row.get(smi_c) or "").strip()
            if not smi:
                continue
            for col, label in aff:
                pch, censored = to_pchembl(row.get(col))
                if censored:
                    n_censored += 1
                    continue
                if pch is None:
                    continue
                rec = pairs.setdefault((uni, smi), {})
                rec.setdefault(label, []).append(pch)
                rec["_id"] = (row.get(id_c) or "").strip() if id_c else None

    out: dict = {}
    for (uni, smi), measures in pairs.items():
        # Ki and Kd before IC50 and EC50: the first two are equilibrium
        # constants and the others depend on how the assay was run.
        chosen = next((label for _, label in _AFFINITY_COLS
                       if label in measures), None)
        if chosen is None:
            continue
        vals = measures[chosen]
        out.setdefault(uni, []).append({
            "chembl_id": measures.get("_id") or None,
            "smiles": smi,
            # Median, not best: taking the best would favour compounds
            # assayed many times, which is study effort rather than
            # potency.
            "pchembl": round(statistics.median(vals), 2),
            "measure": chosen,
            "n_measurements": len(vals),
        })

    for uni in out:
        out[uni] = [r for r in out[uni] if r["pchembl"] >= min_pchembl]
        out[uni].sort(key=lambda r: -r["pchembl"])

    log.info("BindingDB: %d target(s) of %d requested, %d ligand(s) over "
             "pChEMBL %.1f, from %d row(s); %d censored value(s) dropped",
             len(out), len(want), sum(len(v) for v in out.values()),
             min_pchembl, n_rows, n_censored)
    return out


def load_chembl_sqlite(uniprots, path=None, *,
                       min_pchembl: float = 5.0) -> dict:
    """The same shape from a ChEMBL SQLite release.

    Uses the standard library's sqlite3, so the authoritative source needs
    no extra dependency - only the disk space.
    """
    import sqlite3

    want = [u.strip().upper() for u in uniprots if u]
    if not want:
        return {}
    p = resolve(path, DEFAULT_CHEMBL_SQLITE_PATHS)
    if not p:
        return {}

    sql = """
        SELECT cs.accession, m.chembl_id, st.canonical_smiles,
               act.pchembl_value, act.standard_type
        FROM component_sequences cs
        JOIN target_components tc ON tc.component_id = cs.component_id
        JOIN assays a            ON a.tid = tc.tid
        JOIN activities act      ON act.assay_id = a.assay_id
        JOIN molecule_dictionary m ON m.molregno = act.molregno
        JOIN compound_structures st ON st.molregno = act.molregno
        WHERE cs.accession IN ({})
          AND act.pchembl_value IS NOT NULL
          AND act.standard_type IN ('Ki','Kd','IC50','EC50')
    """.format(",".join("?" * len(want)))

    out: dict = {}
    try:
        con = sqlite3.connect(f"file:{p}?mode=ro", uri=True)
        try:
            grouped: dict = {}
            for acc, cid, smi, pch, typ in con.execute(sql, want):
                if not smi or pch is None:
                    continue
                grouped.setdefault((acc.upper(), smi),
                                   {"id": cid, "vals": [], "type": typ}
                                   )["vals"].append(float(pch))
            for (acc, smi), g in grouped.items():
                out.setdefault(acc, []).append({
                    "chembl_id": g["id"], "smiles": smi,
                    "pchembl": round(statistics.median(g["vals"]), 2),
                    "measure": g["type"],
                    "n_measurements": len(g["vals"]),
                })
        finally:
            con.close()
    except Exception as exc:                                # noqa: BLE001
        log.warning("ChEMBL SQLite query failed (%s); potency not loaded "
                    "from it", str(exc)[:120])
        return {}

    for acc in out:
        out[acc] = [r for r in out[acc] if r["pchembl"] >= min_pchembl]
        out[acc].sort(key=lambda r: -r["pchembl"])
    log.info("ChEMBL SQLite: %d target(s), %d ligand(s) over pChEMBL %.1f",
             len(out), sum(len(v) for v in out.values()), min_pchembl)
    return out


def actives_for(uniprot: str, *, bindingdb=None, chembl_db=None,
                min_pchembl: float = 5.0) -> list:
    """Measured ligands of one target, from whichever file is present.

    ChEMBL is preferred when both are there: it is the source the online
    path uses, so an offline run agrees with an online one.
    """
    u = (uniprot or "").strip().upper()
    if not u:
        return []
    if resolve(chembl_db, DEFAULT_CHEMBL_SQLITE_PATHS):
        got = load_chembl_sqlite([u], chembl_db, min_pchembl=min_pchembl)
        if got.get(u):
            return got[u]
    return load_bindingdb([u], bindingdb, min_pchembl=min_pchembl).get(u, [])


def available() -> dict:
    return {"bindingdb": resolve(None, DEFAULT_BINDINGDB_PATHS),
            "chembl_sqlite": resolve(None, DEFAULT_CHEMBL_SQLITE_PATHS)}


def summarise(actives: list, symbol: str | None = None) -> str:
    avail = available()
    if not any(avail.values()):
        return ("No local potency file is present, so measured affinity is "
                "unavailable and anything ranking on it is unranked rather "
                "than scored zero.\n" + SOURCE)
    src = "ChEMBL SQLite" if avail["chembl_sqlite"] else "BindingDB"
    if not actives:
        return (f"{src} is loaded and holds no measured ligand above the "
                f"potency floor for {symbol or 'this target'}. For a target "
                f"nobody has assayed, that is a gap in the literature rather "
                f"than evidence the target is undruggable.")
    best = max(a["pchembl"] for a in actives)
    kinds = sorted({a["measure"] for a in actives})
    return (f"{len(actives)} measured ligand(s) for {symbol or 'this target'} "
            f"from {src}, best pChEMBL {best:.1f}, from {', '.join(kinds)}. "
            f"Censored values such as '>10000 nM' were dropped rather than "
            f"converted, since a failed experiment is not a weak binding "
            f"constant, and repeated measurements are summarised by their "
            f"median rather than their best.")
