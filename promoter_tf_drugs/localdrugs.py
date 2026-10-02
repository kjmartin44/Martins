"""Approved drugs and their known targets, from a file instead of an API.

Why this exists
---------------
The docking and screening paths need three things from ChEMBL: the gene's
ChEMBL target id, the measured ligands at that target for the validation
gate, and the approved-drug library to screen. On a network that blocks or
breaks those calls, none of it runs - and a corporate network is the normal
case, not the exception.

The circuit breaker already stops a dead host costing more than one
request's worth of waiting. That makes a broken endpoint cheap; it does not
make the analysis possible. For that the data has to be on disk.

The one file that covers all three
----------------------------------
The Broad Drug Repurposing Hub publishes a single annotated table: around
6,800 compounds with SMILES, clinical phase, mechanism of action and target
gene symbols. From that one file you get the screening library, the clinical
status to filter it by, and the known ligands per target that the docking
validation gate needs - which is the whole ChEMBL dependency, offline, in a
few megabytes.

  https://repurposing.broadinstitute.org/download
  the file named repurposing_drugs_*.txt plus repurposing_samples_*.txt,
  or the merged export; any of them is accepted.

What a local table cannot give you
----------------------------------
Potency. The Hub records that a compound acts on a target, not how hard it
binds, so an `actives` set from it has no pChEMBL values. Everything that
reads potency - the established-ligand ranking in `screening`, the
`drug_scoring` potency axis - degrades to unranked rather than silently
scoring zero. The docking gate does not need potency, only the identity of
the known binders, so docking is unaffected.

It is also a curated snapshot. A compound approved after the file was
published is absent, and absence here is not evidence that nothing exists.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

log = logging.getLogger(__name__)

DEFAULT_LIBRARY_PATHS = (
    "data/repurposing_drugs.txt", "data/repurposing_hub.tsv",
    "data/approved_library.tsv", "data/approved_library.csv",
)

SOURCE = (
    "Broad Drug Repurposing Hub: https://repurposing.broadinstitute.org/"
    "download - save the drugs table as data/repurposing_drugs.txt. "
    "Any two-column name/SMILES table also works."
)

# Column names these tables actually use, in the order to try them.
_NAME_COLS = ("pert_iname", "name", "drug", "drug_name", "compound",
              "preferred_name", "title")
_SMILES_COLS = ("smiles", "canonical_smiles", "SMILES", "structure")
_TARGET_COLS = ("target", "targets", "gene", "gene_symbol", "target_symbol")
_PHASE_COLS = ("clinical_phase", "phase", "max_phase", "max_clinical_phase")

# Values in a clinical-phase column that mean "in people, approved".
_APPROVED = {"launched", "approved", "4", "4.0", "phase 4", "phase_4"}


def _sniff(path: Path) -> str:
    """Delimiter, from the header rather than the file extension."""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        head = fh.readline()
    counts = {d: head.count(d) for d in ("\t", ",", ";", "|")}
    best = max(counts, key=counts.get)
    return best if counts[best] else "\t"


def _pick(fieldnames, candidates) -> str | None:
    lowered = {(f or "").strip().lower(): f for f in (fieldnames or [])}
    for c in candidates:
        if c.lower() in lowered:
            return lowered[c.lower()]
    return None


def resolve(explicit=None, candidates=DEFAULT_LIBRARY_PATHS) -> str | None:
    for c in ([explicit] if explicit else list(candidates)):
        if c and Path(c).exists():
            return str(c)
    return None


def read_table(path, *, approved_only: bool = True) -> list:
    """`[{name, smiles, targets, phase}]` from a local drug table.

    Rows with no structure are dropped: a name alone cannot be docked or
    fingerprinted, so keeping them would pad the library with compounds
    that silently fail later.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{path} not found. {SOURCE}")

    delim = _sniff(path)
    out: list = []
    n_rows = n_no_smiles = n_filtered = 0

    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        name_c = _pick(reader.fieldnames, _NAME_COLS)
        smi_c = _pick(reader.fieldnames, _SMILES_COLS)
        tgt_c = _pick(reader.fieldnames, _TARGET_COLS)
        ph_c = _pick(reader.fieldnames, _PHASE_COLS)

        if smi_c is None:
            raise ValueError(
                f"{path} has no SMILES column. Looked for "
                f"{', '.join(_SMILES_COLS)}; found "
                f"{', '.join(reader.fieldnames or [])}. A table without "
                f"structures cannot be screened or docked.")

        for row in reader:
            n_rows += 1
            smi = (row.get(smi_c) or "").strip()
            if not smi:
                n_no_smiles += 1
                continue
            phase = (row.get(ph_c) or "").strip() if ph_c else ""
            if approved_only and phase:
                if phase.strip().lower() not in _APPROVED:
                    n_filtered += 1
                    continue
            targets = []
            if tgt_c:
                raw = (row.get(tgt_c) or "").strip()
                targets = [t.strip().upper() for t in
                           raw.replace(";", "|").replace(",", "|").split("|")
                           if t.strip()]
            name = (row.get(name_c) or "").strip() if name_c else ""
            out.append({"name": name or smi[:24], "smiles": smi,
                        "targets": targets, "phase": phase or None,
                        "chembl_id": None})

    log.info("local drug table %s: %d usable of %d row(s)"
             "%s%s", path.name, len(out), n_rows,
             f", {n_no_smiles} without a structure" if n_no_smiles else "",
             f", {n_filtered} not approved" if n_filtered else "")
    if approved_only and not any(r["phase"] for r in out):
        log.info("no clinical-phase column was found in %s, so nothing was "
                 "filtered on approval status - the whole table is being "
                 "treated as the library", path.name)
    return out


def library(path=None, *, approved_only: bool = True) -> list:
    """The screening library from a local table, or [] when none is on disk."""
    p = resolve(path)
    if not p:
        log.info("no local drug table found; the library needs the network. "
                 "%s", SOURCE)
        return []
    return read_table(p, approved_only=approved_only)


def actives_for(symbol: str, rows: list, *, uniprot: str | None = None,
                with_potency: bool = True) -> list:
    """Known ligands of one target, from the table's own annotation.

    This is what replaces `agonism.fetch_actives` offline. It carries no
    pChEMBL, because the table records that a compound acts on a target and
    not how hard it binds - so `pchembl` is None rather than a fabricated
    number, and anything ranking on potency degrades to unranked.
    """
    want = (symbol or "").strip().upper()
    if not want:
        return []

    # A local potency file, when there is one, answers this better: it
    # carries measured affinity rather than only the fact of an
    # interaction, so the potency axis has something to read.
    if with_potency and uniprot:
        try:
            from . import localpotency as lpot
            if any(lpot.available().values()):
                measured = lpot.actives_for(uniprot)
                if measured:
                    log.info("local potency file: %d measured ligand(s) for "
                             "%s (%s)", len(measured), want, uniprot)
                    return measured
        except Exception as exc:                            # noqa: BLE001
            log.debug("local potency unavailable: %s", str(exc)[:80])

    out = []
    for r in rows:
        if want in (r.get("targets") or []):
            out.append({"chembl_id": r.get("chembl_id"),
                        "name": r.get("name"),
                        "smiles": r["smiles"], "pchembl": None})
    log.info("local table: %d known ligand(s) for %s, without potency",
             len(out), want)
    return out


def summarise(rows: list, symbol: str | None = None) -> str:
    if not rows:
        return ("No local drug table is on disk, so the library and the "
                "known-ligand set both need the network. " + SOURCE)
    n_t = len({t for r in rows for t in (r.get("targets") or [])})
    bits = [f"{len(rows)} compound(s) with structures from a local table, "
            f"annotated against {n_t} target(s)."]
    if symbol:
        n = len(actives_for(symbol, rows))
        bits.append(
            f"{n} of them are recorded as acting on {symbol}"
            + (", which is the validation set the docking gate needs."
               if n else " - too few to validate a pocket against, so "
                         "docking will refuse rather than report."))
    bits.append(
        "A local table carries no potency, so anything that ranks on it is "
        "unranked rather than scored zero.")
    return " ".join(bits)
