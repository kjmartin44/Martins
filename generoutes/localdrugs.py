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
from dataclasses import dataclass
from pathlib import Path

from .congeners import salt_free_name

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
# Read off the real files rather than guessed. Getting these wrong does
# not raise: the name silently falls back to the SMILES string, which is
# what the report then prints and what the dedup fallback keys on.
_NAME_COLS = (
    "pert_iname",            # Broad Repurposing Hub
    "inn",                   # DrugCentral (structures.smiles)
    "drugname", "drug_name", "drug name",   # TTD, SuperDRUG2
    "name", "drug", "compound", "preferred_name", "title",
    "generic_name", "common_name", "ingredient", "active_ingredient",
    "molecule_name", "pref_name",           # ChEMBL exports
)
_SMILES_COLS = (
    "smiles", "canonical_smiles", "SMILES", "structure",
    "canonical smiles", "isomeric_smiles", "openbabel_smiles",
)
_TARGET_COLS = ("target", "targets", "gene", "gene_symbol", "target_symbol",
                "gene_name", "uniprot", "accession")
# Identifier columns. Dropping these cost the offline screening paths
# more than it looked: every screen keys its known ligands, its potency
# lookups and - worst - its decoy EXCLUSION on chembl_id, so a library
# where every id is None let training actives be sampled as decoys.
_ID_COLS = ("chembl_id", "molecule_chembl_id", "chembl", "drug_id",
            "broad_id", "pert_id", "id")

_PHASE_COLS = ("clinical_phase", "phase", "max_phase", "max_clinical_phase",
               "highest_clinical_phase", "status", "development_status")

# Values in a clinical-phase column that mean "in people, approved".
_APPROVED = {"launched", "approved", "4", "4.0", "phase 4", "phase_4"}


# ===================================================================== #
# Where else approved molecules come from.                              #
# ===================================================================== #
# No single source has all of them, and they disagree in useful ways: one
# carries target annotations, another carries approval dates per agency,
# another is already curated for docking. Merging on the parent structure
# is what makes using several of them safe.
#
# The counts are distinct molecules, not products, and they overlap
# heavily. Adding a second source typically contributes a few hundred new
# moieties rather than doubling anything - `merge` reports exactly how
# many, per source, so the decision is not guesswork.

@dataclass(frozen=True)
class DrugSource:
    key: str
    name: str
    approx_molecules: int
    url: str
    filename: str
    gives: str
    caveat: str


SOURCES: tuple = (
    DrugSource(
        "repurposing_hub", "Broad Drug Repurposing Hub", 6800,
        "https://repurposing.broadinstitute.org/download",
        "data/repurposing_drugs.txt",
        "SMILES, clinical phase, mechanism of action, target gene symbols",
        "the only one here that gives target symbols AND phase in one "
        "file, which is why it is the default; about 2,400 of its entries "
        "are launched and the rest are earlier phases"),
    DrugSource(
        "drugcentral", "DrugCentral", 4900,
        "https://drugcentral.org/download",
        "data/drugcentral.smiles",
        "SMILES, INN, targets, indications, and approval dates per agency",
        "the best metadata of any free source, and the only one that says "
        "WHICH agency approved a drug - so it can separate FDA-approved "
        "from approved-elsewhere, which ChEMBL's max_phase cannot"),
    DrugSource(
        "chembl_approved", "ChEMBL, max_phase 4", 4000,
        "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/",
        "data/chembl_approved.tsv",
        "SMILES and clinical phase for drugs approved anywhere",
        "counts salts as separate molecules and includes drugs never "
        "approved in the US; the merge collapses the salts, and "
        "DrugCentral is what distinguishes the agencies"),
    DrugSource(
        "zinc_fda", "ZINC curated FDA subset", 1600,
        "https://zinc.docking.org/substances/subsets/fda/",
        "data/zinc_fda.smi",
        "SMILES only, already curated to the approved set",
        "no names, no targets, no phase - just structures, identified by "
        "ZINC id. Worth having because the subset exists to be docked, so "
        "its curation is doing the work the metadata would"),
    DrugSource(
        "ttd", "Therapeutic Target Database", 2800,
        "https://db.idrblab.net/ttd/full-data-download",
        "data/ttd_drugs.txt",
        "SMILES, drug names, and the targets each drug acts on",
        "target coverage is good and the clinical status vocabulary is its "
        "own, so phases may not map onto the standard values"),
)

BY_KEY = {s.key: s for s in SOURCES}


def source_status() -> list:
    """Which of the known sources are on disk."""
    return [{"key": s.key, "name": s.name, "file": s.filename,
             "present": Path(s.filename).exists(),
             "approx_molecules": s.approx_molecules}
            for s in SOURCES]


def describe_sources() -> str:
    """What to download, and what each one adds that the others do not."""
    L = ["Approved-drug sources this package can merge", ""]
    for s in SOURCES:
        mark = "PRESENT" if Path(s.filename).exists() else "missing"
        L.append(f"  [{mark:7s}] {s.name}  (~{s.approx_molecules:,} "
                 f"molecules)")
        L.append(f"              {s.url}")
        L.append(f"              -> {s.filename}")
        L.append(f"              gives: {s.gives}")
        L.append(f"              note:  {s.caveat}")
        L.append("")
    L.append("They overlap heavily, so a second source usually adds a few "
             "hundred new moieties rather than doubling the library. Merge "
             "them and the per-source accounting says how many:")
    L.append("")
    L.append("  python -m generoutes.localdrugs --merge "
             + " ".join(s.filename for s in SOURCES[:3])
             + " --out data/library.tsv")
    L.append("")
    L.append("If you only take two: the Repurposing Hub for target "
             "annotations and phase, and DrugCentral for per-agency "
             "approval - which is the only way to say FDA-approved rather "
             "than approved somewhere.")
    return "\n".join(L)


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


# Extensions that mean "no header, SMILES first, optional id second".
_SMI_SUFFIXES = (".smi", ".ism", ".can")


def _looks_headerless(path: Path, delim: str) -> bool:
    """True when the first line is data rather than column names.

    ZINC's curated subsets, and most .smi files, have no header at all:
    the first line is already a molecule. Passing that to DictReader
    consumes a real compound as the header and then fails to find a SMILES
    column, which is how `zinc_fda.smi` was rejected outright.
    """
    if path.suffix.lower() in _SMI_SUFFIXES:
        return True
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        first = fh.readline().strip()
    if not first:
        return False
    fields = [f for f in first.replace(delim, " ").split() if f]
    if not fields:
        return False
    known = set(_NAME_COLS) | {c.lower() for c in _SMILES_COLS}
    if any(f.lower() in known for f in fields):
        return False
    # A header field is a word; a SMILES is not. Ring-closure digits and
    # bond characters are the giveaway.
    return any(ch in fields[0] for ch in "()=#[]@") or fields[0].isupper()


def read_smiles_file(path, *, approved_only: bool = True) -> list:
    """A headerless SMILES file: structure first, optional identifier next.

    This is what ZINC's `fda` and `world` subsets ship as, and they are
    worth supporting for a reason beyond format: that subset exists to be
    docked, so it is already curated to the approved set.

    `approved_only` cannot be applied here and is not silently accepted.
    The format carries no phase column, so there is nothing to filter on:
    approval is whatever the file's own curation asserts. A caller asking
    for it gets a warning rather than the quiet impression that a filter
    ran, which is how a `world` subset could be read as approved-only.
    """
    path = Path(path)
    if approved_only:
        log.warning(
            "%s is a plain SMILES file with no phase column, so "
            "approved_only could not be applied. Every structure in the "
            "file is returned and approval rests on the file's own "
            "curation - ZINC's `fda` subset is approved, `world` is not.",
            path.name)
    out: list = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            fields = line.split()
            if not fields:
                continue
            smi = fields[0]
            if smi.lower() in {"smiles", "structure"}:
                continue            # a one-word header, if there is one
            name = fields[1] if len(fields) > 1 else smi[:24]
            out.append({"name": name, "smiles": smi, "targets": [],
                        # No phase column exists, so the file's own
                        # curation is the only claim about approval - which
                        # is why the subset matters more than the format.
                        "phase": None, "chembl_id": None})
    log.info("SMILES file %s: %d structure(s)", path.name, len(out))
    return out


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
    if _looks_headerless(path, delim):
        return read_smiles_file(path, approved_only=approved_only)

    out: list = []
    n_rows = n_no_smiles = n_filtered = 0

    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        name_c = _pick(reader.fieldnames, _NAME_COLS)
        smi_c = _pick(reader.fieldnames, _SMILES_COLS)
        tgt_c = _pick(reader.fieldnames, _TARGET_COLS)
        ph_c = _pick(reader.fieldnames, _PHASE_COLS)
        id_c = _pick(reader.fieldnames, _ID_COLS)

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
            # Gated on the COLUMN, not on the cell. Gating on the cell let
            # a blank pass while "Phase 1" was dropped, so a row with a gap
            # in its phase column entered the library as an approved drug.
            # The "no phase column at all, so nothing to filter on" case is
            # the `if ph_c else ""` above, and the log line below says so.
            if approved_only and ph_c:
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
            cid = (row.get(id_c) or "").strip() if id_c else ""
            out.append({"name": name or smi[:24], "smiles": smi,
                        "targets": targets, "phase": phase or None,
                        # Read, not discarded. None only when the table
                        # genuinely has no identifier column.
                        "chembl_id": cid or None,
                        # Which column the name came from, so resolve_name
                        # can prefer an INN over a brand name.
                        "_name_column": name_c})

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


# ===================================================================== #
# Merging several sources into one library.                             #
# ===================================================================== #
# "There are 20,000 FDA-approved drugs and I only have 3,000" is the usual
# reaction to a library this size, and the premise is a counting difference
# rather than a gap.
#
# The FDA's often-quoted ~20,000 counts approved PRODUCTS: every strength,
# every formulation, every manufacturer of the same molecule listed
# separately. Metformin alone accounts for hundreds of entries. Distinct
# active moieties with a defined structure is a far smaller number -
# roughly 1,500 to 2,000 small molecules ever approved in the US, and about
# 4,000 in ChEMBL at max_phase 4, which also counts drugs approved
# elsewhere and counts salts separately.
#
# So a library of ~3,000 is not a seventh of the drugs. It is most of the
# distinct molecules, and the way to make it genuinely larger is to WIDEN
# THE DEFINITION deliberately - approved anywhere rather than in the US,
# compounds in trials, withdrawn drugs with human safety data - not to hunt
# for missing FDA approvals.
#
# This merges sources and reports what each contributed, so widening is a
# decision with a number attached.

# A withdrawn drug is not a lower phase: it REACHED approval and was then
# pulled. Putting it at -1 sorted it below preclinical, so asking for
# withdrawn drugs dragged in every lab compound. It sits at 4 with a
# separate flag, and `include_withdrawn` decides it.
PHASE_ORDER = {"launched": 4, "approved": 4, "4": 4, "4.0": 4,
               "withdrawn": 4,
               "phase 3": 3, "3": 3, "3.0": 3,
               "phase 2": 2, "2": 2, "2.0": 2,
               "phase 1": 1, "1": 1, "1.0": 1,
               "preclinical": 0, "": 0}
WITHDRAWN = {"withdrawn", "withdrawn from market", "suspended"}


def parent_and_key(smiles: str) -> tuple:
    """`(parent SMILES, dedup key)`, or `(None, None)`.

    Both come out of one standardisation because both need it and it is
    the expensive part. `merge` kept the key and threw the structure
    away, storing whichever row created the group - so a molecule
    collapsed "on the parent structure", as `explain_counts` tells the
    reader, could carry a salt form's SMILES. 98 of 2,119 in the shipped
    library did, alendronic acid with its sodium among them, and that
    SMILES is what docking and the fingerprint screens then use.

    An InChIKey alone does NOT merge salts: the counterion is part of the
    molecule, so metformin and metformin hydrochloride get entirely
    different keys. Measured, not assumed. So the parent is extracted
    first - largest fragment, uncharged, tautomer-canonicalised - and the
    key is the connectivity layer of that molecule's InChIKey.

    Residual limitation, also measured: a compound drawn in an explicit
    protonation state can still miss. Metformin written with a charged
    guanidinium does not merge with neutral metformin, because a
    delocalised charge cannot be canonically neutralised. Databases
    normally write the neutral parent or a simple salt, so the common cases
    are covered and this one is rare.
    """
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem.MolStandardize import rdMolStandardize
    except ImportError:
        return None, None
    RDLogger.DisableLog("rdApp.*")
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, None
        mol = rdMolStandardize.ChargeParent(mol)
        mol = rdMolStandardize.TautomerEnumerator().Canonicalize(mol)
        key = Chem.MolToInchiKey(mol)
        return (Chem.MolToSmiles(mol), key[:14]) if key else (None, None)
    except Exception as exc:                                # noqa: BLE001
        # One structure RDKit cannot standardise gets no key, so it is
        # matched by name only. Per molecule, so debug rather than warning.
        log.debug("no structure key for %s: %s", smiles[:40],
                  type(exc).__name__)
        return None, None


# ===================================================================== #
# Conflicts between sources, and how each is resolved.                  #
# ===================================================================== #
# Sources disagree, and the first version of merge() resolved every
# disagreement silently. One of those silences was dangerous.
#
# Rofecoxib: the Repurposing Hub called it launched and ChEMBL called it
# withdrawn. Because a withdrawn row was FILTERED OUT before anything was
# recorded from it, the withdrawal information went with it, and the drug
# survived in the library through a source that did not know it had been
# pulled - reading `launched`, `withdrawn=False`. That is Vioxx, removed
# for cardiovascular deaths, presented as an approved option.
#
# So withdrawal is now collected from every row before any filtering, and
# it propagates on contact: if ANY source says a molecule was withdrawn,
# the merged molecule is withdrawn. That is deliberately not the rule used
# for clinical phase, where the maximum wins. Phase is a claim about how
# far something got and the best-informed source is usually right;
# withdrawal is a safety signal and the cost of ignoring one is not
# symmetric with the cost of carrying one too many.

# How much a source's name is to be trusted, highest first. An INN is an
# international non-proprietary name, which is what the report should
# print; a brand name and a ZINC id are both worse, for different reasons.
NAME_PREFERENCE = ("inn", "pert_iname", "drugname", "drug_name",
                   "drug name", "pref_name", "name")


def _name_rank(column: str | None) -> int:
    col = (column or "").strip().lower()
    return (NAME_PREFERENCE.index(col) if col in NAME_PREFERENCE
            else len(NAME_PREFERENCE))


def _looks_like_an_id(name: str) -> bool:
    """ZINC ids, ChEMBL ids and the like are identifiers, not names."""
    n = (name or "").strip()
    if not n:
        return True
    upper = n.upper()
    if upper.startswith(("ZINC", "CHEMBL", "DB0", "DB1", "CID")):
        return True
    # mostly digits is an accession rather than a drug
    return sum(c.isdigit() for c in n) > len(n) * 0.5


def resolve_name(candidates: list) -> tuple:
    """(chosen name, conflict or None) from [(name, source, column)].

    Prefers a real name over an identifier, then an INN over a brand name,
    then lowercase over shouted. A disagreement between two genuine names
    is reported rather than resolved away, because a brand name silently
    replacing a generic one changes what the reader searches for.
    """
    real = [(n, s, c) for n, s, c in candidates if not _looks_like_an_id(n)]
    pool = real or candidates
    if not pool:
        return None, None

    # Agreement between sources beats column provenance. A column named
    # INN is supposed to hold an international non-proprietary name, but a
    # curation slip can put a brand name there - and when two of three
    # sources say "metformin" and one says "Glucophage", the count is the
    # better signal than the header.
    agree: dict = {}
    for n, _s, _c in pool:
        agree[n.strip().lower()] = agree.get(n.strip().lower(), 0) + 1

    def key(item):
        n, _src, col = item
        return (-agree[n.strip().lower()], _name_rank(col),
                n.isupper(), len(n), n.lower())

    pool = sorted(pool, key=key)
    chosen = pool[0][0]

    distinct = {n.strip().lower() for n, _, _ in real}
    conflict = None
    if len(distinct) > 1:
        conflict = {
            "kind": "name",
            "chosen": chosen,
            "alternatives": sorted(distinct - {chosen.strip().lower()}),
            "by_source": {s: n for n, s, _ in real},
            "resolution": ("kept the name most sources agree on, then "
                           "best column provenance - an INN over a brand "
                           "name over an identifier. Agreement comes first "
                           "because a curation slip can put a brand name "
                           "in an INN column. The alternatives are listed "
                           "because a brand name replacing a generic one "
                           "changes what you would search for"),
        }
    return chosen, conflict


def _phase_conflict(entries: list) -> dict | None:
    """entries: [(raw_phase, rank, source)]."""
    ranks = {r for _, r, _ in entries}
    if len(ranks) <= 1:
        return None
    best = max(entries, key=lambda e: e[1])
    return {
        "kind": "phase",
        "chosen": best[0],
        "by_source": {s: p for p, _, s in entries},
        "resolution": ("took the highest phase any source claims. A source "
                       "recording a lower phase is usually an older "
                       "snapshot, but the one claiming approval is the one "
                       "to check if this drug matters to you"),
    }


def _target_conflict(by_source: dict) -> dict | None:
    """Disjoint target sets usually mean the sources are annotating
    different things about the same molecule, not the same thing twice."""
    sets = {s: set(t) for s, t in by_source.items() if t}
    if len(sets) < 2:
        return None
    all_t = set().union(*sets.values())
    common = set.intersection(*sets.values())
    if common:
        return None
    return {
        "kind": "targets",
        "chosen": sorted(all_t),
        "by_source": {s: sorted(t) for s, t in sets.items()},
        "resolution": ("unioned, and flagged because the sources name NO "
                       "target in common. That often means they are "
                       "annotating different relationships - an efficacy "
                       "target in one and an uptake transporter in the "
                       "other, as with metformin's PRKAB1 against SLC22A1 "
                       "- so treat a route through either as unconfirmed "
                       "by the other source"),
    }


def parent_key(smiles: str) -> str | None:
    """The dedup key alone, for callers that do not need the structure."""
    return parent_and_key(smiles)[1]


def merge(paths, *, min_phase: int = 4, include_withdrawn: bool = False,
          fallback_name_key: bool = True) -> dict:
    """Combine several drug tables into one deduplicated library.

    `min_phase` is the widening dial: 4 keeps approved or launched only, 3
    adds phase 3, 1 adds everything given to a human. `include_withdrawn`
    is separate, because a withdrawn drug reached approval and was then
    pulled - it carries human safety data and a documented reason for
    removal, so it is an opt-in rather than a point on the phase scale.

    Returns the library plus an accounting of what each source contributed
    and what was collapsed, because more rows is only an improvement when
    the extra rows are distinct molecules.
    """
    merged: dict = {}
    per_source: dict = {}
    n_rows = n_no_key = n_below_phase = 0

    # Pass one: read everything, key it, and record every source's claim.
    # Nothing is filtered yet. The first version filtered withdrawn rows
    # here, which threw away the withdrawal along with the row - so
    # rofecoxib survived through a source that did not know it had been
    # pulled and read back as `launched`, `withdrawn=False`.
    for path in paths:
        src = Path(path).name
        try:
            rows = read_table(path, approved_only=False)
        except Exception as exc:                            # noqa: BLE001
            per_source[src] = {
                "error": type(exc).__name__ + ": "
                         + str(exc).replace("\n", " ")[-220:]}
            continue
        added = already = 0
        for r in rows:
            n_rows += 1
            raw_phase = (r.get("phase") or "").strip().lower()
            phase = PHASE_ORDER.get(raw_phase, 0)
            parent, key = parent_and_key(r["smiles"])
            unparsed = False
            if key is None:
                if not fallback_name_key:
                    n_no_key += 1
                    continue
                stem = salt_free_name((r.get("name") or "").lower())
                if not stem:
                    n_no_key += 1
                    continue
                key = "name:" + stem
                unparsed = True

            if key in merged:
                already += 1
                rec = merged[key]
            else:
                added += 1
                rec = merged[key] = {
                    # The standardised parent, not whichever salt form
                    # happened to create the group. `explain_counts` tells
                    # the reader the forms were "collapsed on the parent
                    # structure", and this is what makes that true - it is
                    # also the structure docking and the fingerprint
                    # screens use for the molecule.
                    "key": key, "smiles": parent or r["smiles"],
                    "structure_unparsed": unparsed,
                    "sources": [], "n_forms": 0, "chembl_id": None,
                    "_names": [], "_phases": [], "_targets": {},
                }
            # The first identifier any form of this molecule carries. The
            # screens key their known ligands, their potency lookups and -
            # worst - their decoy EXCLUSION on it, so a merged library with
            # None everywhere let training actives be sampled as decoys.
            if rec["chembl_id"] is None and r.get("chembl_id"):
                rec["chembl_id"] = r["chembl_id"]
            if src not in rec["sources"]:
                rec["sources"].append(src)
            rec["n_forms"] += 1
            name = (r.get("name") or "").strip()
            if name:
                rec["_names"].append((name, src, r.get("_name_column")))
            rec["_phases"].append((r.get("phase") or "", phase, src))
            if r.get("targets"):
                rec["_targets"].setdefault(src, set()).update(r["targets"])
        per_source[src] = {"rows": len(rows), "new": added,
                           "already_present": already}

    # Pass two: resolve each molecule's conflicts, then filter.
    conflicts: list = []
    library: list = []
    n_withdrawn = 0
    for rec in merged.values():
        name, name_conflict = resolve_name(rec["_names"])
        phases = rec["_phases"]
        best_phase = max(phases, key=lambda e: e[1]) if phases else ("", 0, "")
        # Withdrawal propagates on contact. Any source saying a molecule
        # was pulled makes it withdrawn, which is deliberately NOT the
        # max-phase rule used above: phase is a claim about how far
        # something got and the best-informed source is usually right,
        # while withdrawal is a safety signal whose costs are asymmetric.
        withdrawn_by = sorted({src for raw, _, src in phases
                               if raw.strip().lower() in WITHDRAWN})
        is_withdrawn = bool(withdrawn_by)

        if name_conflict:
            name_conflict["drug"] = name
            conflicts.append(name_conflict)
        pc = _phase_conflict(phases)
        if pc:
            pc["drug"] = name
            conflicts.append(pc)
        tc = _target_conflict(rec["_targets"])
        if tc:
            tc["drug"] = name
            conflicts.append(tc)
        if is_withdrawn and len(withdrawn_by) < len(rec["sources"]):
            conflicts.append({
                "kind": "withdrawn",
                "drug": name,
                "chosen": "withdrawn",
                "by_source": {src: ("withdrawn" if src in withdrawn_by
                                    else raw or "listed")
                              for raw, _, src in phases},
                "resolution": (
                    f"{', '.join(withdrawn_by)} records this as withdrawn "
                    f"and the other source(s) do not. Withdrawal wins: a "
                    f"drug that was pulled must not be presented as an "
                    f"approved option because one snapshot predates the "
                    f"withdrawal"),
                "severity": "safety",
            })

        if is_withdrawn and not include_withdrawn:
            n_withdrawn += 1
            continue
        if best_phase[1] < min_phase:
            n_below_phase += 1
            continue

        targets = sorted(set().union(*rec["_targets"].values())
                         if rec["_targets"] else set())
        library.append({
            "name": name, "smiles": rec["smiles"], "targets": targets,
            "phase": ("withdrawn" if is_withdrawn else best_phase[0] or None),
            "withdrawn": is_withdrawn,
            "withdrawn_by": withdrawn_by,
            "chembl_id": rec["chembl_id"], "sources": rec["sources"],
            "n_forms": rec["n_forms"], "key": rec["key"],
            "structure_unparsed": rec["structure_unparsed"],
        })

    library = sorted(library, key=lambda r: (r["name"] or ""))
    collapsed = sum(r["n_forms"] - 1 for r in library)
    return {
        "library": library,
        "n_distinct": len(library),
        "n_rows_read": n_rows,
        "n_salt_forms_collapsed": collapsed,
        "n_below_phase": n_below_phase,
        "n_unkeyable": n_no_key,
        "n_withdrawn_excluded": n_withdrawn,
        "per_source": per_source,
        "min_phase": min_phase,
        "include_withdrawn": include_withdrawn,
        "conflicts": conflicts,
        "n_conflicts": len(conflicts),
        "n_safety_conflicts": sum(1 for c in conflicts
                                  if c.get("severity") == "safety"),
    }


def explain_counts(result: dict) -> str:
    """Why the number is what it is."""
    L = [f"{result['n_distinct']} distinct molecules from "
         f"{result['n_rows_read']} rows across "
         f"{len(result['per_source'])} source(s)."]
    if result["n_salt_forms_collapsed"]:
        L.append(
            f"{result['n_salt_forms_collapsed']} row(s) were salt, hydrate "
            f"or duplicate forms of a molecule already present, collapsed "
            f"on the parent structure - screening the hydrochloride and the "
            f"free base as two candidates doubles the compute and finds the "
            f"same thing twice.")
    if result["n_below_phase"]:
        L.append(f"{result['n_below_phase']} row(s) were below the "
                 f"requested clinical phase.")
    if result["n_unkeyable"]:
        L.append(f"{result['n_unkeyable']} row(s) had neither a usable "
                 f"structure nor a name.")
    if result.get("n_safety_conflicts"):
        L.append(
            f"**{result['n_safety_conflicts']} safety conflict(s):** a "
            f"source records a drug as withdrawn while another still lists "
            f"it. Withdrawal wins in every case, because a drug that was "
            f"pulled must not be offered as approved just because one "
            f"snapshot predates the withdrawal.")
    if result.get("n_conflicts"):
        kinds: dict = {}
        for c in result["conflicts"]:
            kinds[c["kind"]] = kinds.get(c["kind"], 0) + 1
        L.append(f"{result['n_conflicts']} conflict(s) between sources: "
                 + ", ".join(f"{n} {k}" for k, n in sorted(kinds.items()))
                 + ". Each is recorded with what every source said and how "
                   "it was resolved, under `conflicts`.")
    if result.get("n_withdrawn_excluded"):
        L.append(f"{result['n_withdrawn_excluded']} withdrawn drug(s) were "
                 f"excluded. They reached approval and were then pulled, so "
                 f"they carry human safety data and a documented reason for "
                 f"removal; include_withdrawn keeps them.")
    L.append(
        "On the counting: the FDA's often-quoted ~20,000 is approved "
        "PRODUCTS - every strength, formulation and manufacturer of the "
        "same molecule counted separately. Distinct active moieties with a "
        "structure is roughly 1,500 to 2,000 for the US, and about 4,000 in "
        "ChEMBL at max_phase 4, which also counts drugs approved elsewhere. "
        "A library of a few thousand is most of the molecules rather than a "
        "fraction of them.")
    if result["min_phase"] >= 4:
        L.append(
            "To make it genuinely larger, widen the definition on purpose: "
            "min_phase 3 adds phase 3 and min_phase 1 adds everything "
            "given to a human, while include_withdrawn keeps drugs that "
            "were approved and later pulled. Each step changes what a hit "
            "means, so it is a decision rather than a default.")
    return " ".join(L)


def write_table(library: list, path) -> Path:
    """Write a merged library back out as a table this package can read."""
    import csv as _csv

    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = _csv.writer(fh, delimiter="\t")
        # `chembl_id` is a column here because `read_table` reads one
        # (_ID_COLS), so a library written out and read back keeps the
        # identifier the screens exclude their own actives by.
        w.writerow(["name", "clinical_phase", "chembl_id", "target",
                    "smiles", "sources", "n_forms"])
        for r in library:
            w.writerow([r.get("name") or "", r.get("phase") or "",
                        r.get("chembl_id") or "",
                        "|".join(r.get("targets") or []), r["smiles"],
                        "|".join(r.get("sources") or []), r.get("n_forms", 1)])
    return path


def main(argv=None) -> int:
    """`python -m generoutes.localdrugs` - inspect or merge tables."""
    import argparse

    p = argparse.ArgumentParser(
        prog="python -m generoutes.localdrugs",
        description="Inspect a local drug table, or merge several into one "
                    "deduplicated library.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=SOURCE)
    p.add_argument("--sources", action="store_true",
                   help="list the approved-drug sources this can merge, "
                        "which are present, and what each adds that the "
                        "others do not")
    p.add_argument("--merge", nargs="+", metavar="FILE",
                   help="tables to combine. Salt, hydrate and duplicate "
                        "forms collapse onto the parent structure")
    p.add_argument("--out", metavar="FILE",
                   help="write the merged library here, in a format this "
                        "package reads back")
    p.add_argument("--min-phase", type=int, default=4, metavar="N",
                   help="4 approved or launched only, 3 adds phase 3, 1 "
                        "adds anything given to a human. This is the dial "
                        "that genuinely enlarges the library, because there "
                        "are no missing FDA approvals to find (default 4)")
    p.add_argument("--include-withdrawn", action="store_true",
                   help="keep drugs that were approved and later pulled. "
                        "They carry human safety data and a documented "
                        "reason for removal, which is why this is separate "
                        "from --min-phase")
    p.add_argument("--orange-book", nargs="?", const=True, metavar="FILE",
                   help="annotate each molecule with its FDA route of "
                        "administration from the Orange Book, and report how "
                        "many of the library are approved only by a local "
                        "route. Reads products.txt or EOBZIP.zip from disk; "
                        "no network. Pass a path, or the flag alone to look "
                        "in data/")
    p.add_argument("--target", metavar="SYMBOL",
                   help="also report how many of the merged compounds are "
                        "annotated against this target")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)s %(name)s: %(message)s")

    if args.sources:
        print(describe_sources())
        return 0

    if not args.merge:
        found = resolve()
        if not found:
            print("No local drug table found.")
            print()
            print(SOURCE)
            return 1
        rows = read_table(found)
        print(f"{found}: {len(rows)} approved compound(s) with structures")
        print()
        print(summarise(rows, args.target))
        return 0

    res = merge(args.merge, min_phase=args.min_phase,
                include_withdrawn=args.include_withdrawn)
    print(f"Merged library: {res['n_distinct']} distinct molecules")
    print()
    for src, info in res["per_source"].items():
        if "error" in info:
            print(f"  {src:34s} FAILED: {info['error']}")
        else:
            print(f"  {src:34s} {info['rows']:6d} rows -> "
                  f"{info['new']:5d} new, "
                  f"{info['already_present']:5d} already present")
    print()
    print(explain_counts(res))
    if args.orange_book:
        # Route of administration is the one thing the FDA's 20,000-product
        # file has that no structure source does, and the reason to read it
        # at all: a topical-only approval is not a route to an internal
        # gene however well the molecule binds.
        from . import orangebook as ob
        acc = ob.enrich(res["library"],
                        None if args.orange_book is True else args.orange_book)
        print()
        print(ob.summarise(acc))

    if args.target:
        n = sum(1 for r in res["library"]
                if args.target.strip().upper() in (r.get("targets") or []))
        print()
        print(f"{n} of them are annotated against {args.target.upper()}.")
    if args.out:
        write_table(res["library"], args.out)
        print()
        print(f"Wrote {args.out}. Use it with --library on the docking and "
              f"screening paths.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
