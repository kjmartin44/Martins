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
    # The distilled projection first: same columns, a fraction of the
    # size, and reading the 9 GB original once per target is minutes of
    # re-reading the same file.
    "data/bindingdb_distilled.tsv",
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


BINDINGDB_REST = "https://bindingdb.org/rest/getLigandsByUniprot"

# BindingDB prefixes every field in the JSON response, and names the
# response object after a typo of its own method. Both are theirs.
_BDB_ROOT = "getLindsByUniprotResponse"


def fetch_bindingdb(client, uniprots, path=None) -> dict:
    """Add targets to the local potency table, one request each.

    The full BindingDB release is about 1.5 GB zipped, and almost all of
    it is targets this gene has nothing to do with. One accession at a
    time through the cached client is a few hundred kilobytes, re-reads
    from cache afterwards, and honours `--offline` like every other call
    here.

    Merges into whatever is already on disk rather than replacing it, so
    a table built from the cache keeps its rows and gains these.
    """
    import csv as _csv

    want = [u.strip().upper() for u in uniprots if u and u.strip()]
    if not want:
        return {"added": 0, "targets": {}, "error": "no accession given"}

    # Default to the file ALREADY IN USE, not the first candidate path.
    # `resolve` returns the first that exists, so writing to the first
    # candidate created a second file that shadowed it - fetching CCR5
    # took ESR1 from 2,020 actives to zero, because the new file won the
    # lookup and held only CCR5.
    out_path = Path(path or resolve() or DEFAULT_BINDINGDB_PATHS[0])
    if not out_path.suffix:
        out_path = out_path.with_suffix(".tsv")

    # Keyed on (accession, smiles, assay type) so a re-fetch updates in
    # place instead of duplicating.
    rows: dict = {}
    header = ["UniProt (SwissProt) Primary ID of Target Chain",
              "Ligand SMILES", "BindingDB MonomerID",
              "Ki (nM)", "Kd (nM)", "IC50 (nM)", "EC50 (nM)"]
    if out_path.exists():
        with _open(out_path) as fh:
            reader = _csv.DictReader(fh, delimiter="\t")
            # Carry over any column the existing file has that this one
            # does not, including whichever identifier column it uses, so
            # a merge never drops data that was already there.
            for c in (reader.fieldnames or []):
                if c not in header:
                    header.append(c)
            for row in reader:
                uni = (row.get(header[0]) or "").strip().upper()
                smi = (row.get(header[1]) or "").strip()
                if uni and smi:
                    rows[(uni, smi)] = dict(row)

    before = len(rows)
    per_target: dict = {}
    for uni in want:
        try:
            res = client.get(BINDINGDB_REST,
                             {"uniprot": uni,
                              "response": "application/json"})
        except Exception as exc:                            # noqa: BLE001
            per_target[uni] = f"failed: {type(exc).__name__}"
            log.warning("BindingDB fetch failed for %s: %s", uni,
                        str(exc)[:120])
            continue
        root = (res or {}).get(_BDB_ROOT) or {}
        recs = root.get("bdb.affinities") or []
        kept = 0
        for a in recs:
            smi = (a.get("bdb.smile") or "").strip()
            val = (a.get("bdb.affinity") or "").strip()
            kind = (a.get("bdb.affinity_type") or "").strip().upper()
            col = {"KI": "Ki (nM)", "KD": "Kd (nM)", "IC50": "IC50 (nM)",
                   "EC50": "EC50 (nM)"}.get(kind)
            if not (smi and val and col):
                continue
            rec = rows.setdefault((uni, smi), {
                header[0]: uni, header[1]: smi,
                header[2]: str(a.get("bdb.monomerid") or "")})
            # Keep the relation: BindingDB writes ">10000" for a compound
            # that did nothing, and `to_pchembl` is what decides that is
            # not a binding constant.
            prev = rec.get(col)
            if not prev or _looser(prev, val):
                rec[col] = val
            kept += 1
        per_target[uni] = kept
        log.info("BindingDB %s: %d affinity record(s)", uni, kept)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(header) + "\n")
        for (uni, smi), rec in sorted(rows.items()):
            fh.write("\t".join(str(rec.get(c, "") or "")
                                for c in header) + "\n")

    return {"path": str(out_path), "targets": per_target,
            "rows_before": before, "rows_after": len(rows),
            "added": len(rows) - before}


def _looser(prev: str, new: str) -> bool:
    """Prefer an exact measurement, then the more potent one."""
    def rank(t):
        t = (t or "").strip()
        try:
            return (1 if t[:1] in "><~" else 0,
                    float(t.lstrip("><~= ").strip()))
        except (ValueError, IndexError):
            return (2, float("inf"))
    return rank(new) < rank(prev)


DISTILLED = "data/bindingdb_distilled.tsv"


def distil_bindingdb(src=None, out=DISTILLED, *,
                     min_pchembl: float = 0.0) -> dict:
    """Project the full release down to the columns this package reads.

    BindingDB_All is 9 GB across 640 columns and seven of them are used
    here. A pass takes under three minutes, which is tolerable once and
    not once per target - so this writes a small file with the same
    columns, and everything afterwards reads that instead.

    `min_pchembl` defaults to 0: every row with a target, a structure
    and a parseable affinity is carried.

    That is not timidity. `load_bindingdb` collects ALL of a
    target-ligand pair's measurements and then prefers Ki, then Kd, then
    IC50, then EC50, because the first two are equilibrium constants.
    Dropping a weak Ki ROW removes it from that pair, so a stronger IC50
    wins the preference instead and the pair clears a threshold its best
    evidence says it should fail. Measured: at 4.0, CCR5 gained three
    ligands and PPARG seventy-nine against the same read of the full
    file, while the distilled file held strictly FEWER rows.

    Any row-level threshold does this. The projection must change the
    columns and nothing else.
    """
    import csv as _csv

    src = resolve(src)
    if not src:
        return {"error": "no BindingDB file on disk to distil. " + SOURCE}
    out_path = Path(out)
    if Path(src).resolve() == out_path.resolve():
        return {"error": f"{out} is the source; nothing to do"}

    keep = ["UniProt (SwissProt) Primary ID of Target Chain 1",
            "Ligand SMILES", "ChEMBL ID of Ligand",
            "Ki (nM)", "Kd (nM)", "IC50 (nM)", "EC50 (nM)"]
    n_in = n_out = 0
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with _open(src) as fh, open(out_path, "w", encoding="utf-8",
                                newline="") as dst:
        reader = _csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return {"error": f"{src} is empty"}

        # Column INDICES, not DictReader: building a 640-key dict per row
        # is most of the cost of reading this file.
        idx = {}
        for name in keep:
            if name in header:
                idx[name] = header.index(name)
        uni_i = next((idx[c] for c in _UNIPROT_COLS if c in idx), None)
        smi_i = next((idx[c] for c in _SMILES_COLS if c in idx), None)
        if uni_i is None or smi_i is None:
            return {"error": f"{src} has no UniProt or SMILES column; "
                             f"this does not look like BindingDB"}
        aff_i = [(idx[c], c) for c, _ in _AFFINITY_COLS if c in idx]
        id_i = next((idx[c] for c in _ID_COLS if c in idx), None)

        dst.write("\t".join(
            ["UniProt (SwissProt) Primary ID of Target Chain",
             "Ligand SMILES", "ChEMBL ID of Ligand"]
            + [c for _, c in aff_i]) + "\n")

        width = max([uni_i, smi_i] + [i for i, _ in aff_i]
                    + ([id_i] if id_i is not None else [])) + 1
        for row in reader:
            n_in += 1
            if len(row) < width:
                continue
            uni = row[uni_i].strip().upper()
            smi = row[smi_i].strip()
            if not uni or not smi:
                continue
            # Every value on a kept row goes through unchanged. Blanking
            # the weak ones changed which assay type wins the Ki > Kd >
            # IC50 > EC50 preference downstream, so a pair whose best
            # evidence was a weak Ki came back as a potent IC50 binder.
            # The threshold decides which ROWS are worth carrying, never
            # which measurement within one is visible.
            vals = [row[i].strip() for i, _ in aff_i]
            usable = False
            for raw in vals:
                pch, _censored = to_pchembl(raw)
                # A censored bound counts as a measurement worth carrying:
                # the lookup is what decides ">10000 nM" is not a binding
                # constant, and deciding it here would hide it.
                if pch is not None and pch >= min_pchembl:
                    usable = True
                    break
                if raw.strip()[:1] in "><~" and raw.strip()[1:].strip():
                    usable = True
                    break
            if not usable:
                continue
            dst.write("\t".join(
                [uni, smi, (row[id_i].strip() if id_i is not None else "")]
                + vals) + "\n")
            n_out += 1

    size = out_path.stat().st_size
    log.info("distilled %s: %d of %d row(s) kept, %.0f MB",
             out_path.name, n_out, n_in, size / 1e6)
    return {"path": str(out_path), "rows_in": n_in, "rows_out": n_out,
            "bytes": size, "min_pchembl": min_pchembl}


# Per-process memo of what has been read, keyed by (file, accession).
# `actives_for` is called once per target and each call is a full pass, so
# without this a run screening seven targets reads the file seven times.
_LOADED: dict = {}


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
    memo_key = (resolve(bindingdb), u, min_pchembl)
    if memo_key in _LOADED:
        return _LOADED[memo_key]
    if resolve(chembl_db, DEFAULT_CHEMBL_SQLITE_PATHS):
        got = load_chembl_sqlite([u], chembl_db, min_pchembl=min_pchembl)
        if got.get(u):
            _LOADED[memo_key] = got[u]
            return got[u]
    got = load_bindingdb([u], bindingdb, min_pchembl=min_pchembl).get(u, [])
    _LOADED[memo_key] = got
    return got


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
