"""Drug to PROTEIN abundance, where the transcript side has LINCS.

Why this is a separate module from `signatures`
-----------------------------------------------
Every abundance claim elsewhere in this package is about mRNA.
`abundance.py` asks whether a compound moved a factor's own TRANSCRIPT;
`signatures.py` asks the same of the target gene. For a disease where the
pathogenic allele makes an impaired or unstable PROTEIN, that gap is the
whole question: a stabiliser, a chaperone or a degradation blocker raises
protein with the transcript flat, and nothing in this package could see
it.

Which datasets answer this, and which do not
--------------------------------------------
This matters more than the reader, because three well-known resources are
routinely named for this question and two of them do not answer it.

- **LINCS P100** is a reduced-representation PHOSPHOproteomic assay of
  about 90 phosphopeptides. It reports signalling state under drug, not
  the abundance of a protein you choose, and it will not contain most
  proteins. Its sibling **LINCS GCP** measures about 60 histone marks.
  Neither is a drug-to-protein-abundance corpus.
- **CPTAC** is proteogenomics of human tumours. It answers "how much of
  this protein is in this tissue" well, and "does this drug raise it" not
  at all, because it has no drug perturbation arm.
- **Thermal proteome profiling** (TPP, Meltome) measures a melting-point
  shift, which is target ENGAGEMENT rather than abundance. A drug can
  bind a protein without changing how much of it exists.

The resources that do answer it are dose-resolved proteomics under drug:
**decryptE** (around 144 drugs, dose-resolved, whole proteome) and
**ProTargetMiner** (around 50 anticancer compounds). Both are smaller
than LINCS by orders of magnitude, which is the honest limit: the protein
side of this question has a small fraction of the transcript side's
coverage, and a compound absent from the table is unexamined rather than
negative.

So this reads a TABLE rather than one named source. Give it whichever you
can obtain, or one from your own experiment.
"""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass, field
from pathlib import Path

log = logging.getLogger(__name__)

DEFAULT_PATHS = ("data/drug_protein_abundance.tsv",)

# A change this size, in log2, before a row is worth reporting. Proteomics
# is noisier per measurement than RNA-seq, and a 1.2x shift on one run is
# not a finding.
MIN_ABS_LOG2FC = 0.58                      # about 1.5x
# And this many independent observations of it, for the same reason.
MIN_OBSERVATIONS = 2

SOURCE = (
    "No drug-to-protein-abundance table on disk. decryptE (dose-resolved "
    "proteome under ~144 drugs) and ProTargetMiner (~50 anticancer "
    "compounds) answer this question; LINCS P100 (~90 phosphopeptides) "
    "and CPTAC (tumour proteogenomics, no drug arm) do not. Put a table "
    "at data/drug_protein_abundance.tsv with columns "
    "drug, protein, log2fc [, cell_line, dose, time, source].")


@dataclass
class ProteinHit:
    """One compound measured to change one protein's abundance."""
    drug: str
    protein: str
    direction: str                  # "up" | "down"
    median_log2fc: float
    n_observations: int
    cell_lines: list = field(default_factory=list)
    sources: list = field(default_factory=list)
    notes: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


def available(path=None) -> bool:
    cands = [path] if path else list(DEFAULT_PATHS)
    return any(c and Path(c).exists() for c in cands)


def _resolve(path=None) -> Path | None:
    cands = [path] if path else list(DEFAULT_PATHS)
    return next((Path(c) for c in cands if c and Path(c).exists()), None)


def load(protein: str, *, path=None, want: str = "up") -> list[ProteinHit]:
    """Compounds measured to move `protein` the way `want` asks.

    Columns, header required, order free, tab or comma:

        drug, protein, log2fc [, cell_line] [, dose] [, time] [, source]

    `log2fc` is signed, positive meaning the drug RAISED the protein.
    `protein` is matched on the gene symbol, case-insensitively, so a
    UniProt accession column is ignored if present.

    Rows for the same compound are pooled: the median change is reported
    with how many observations it rests on, because one measurement of a
    1.5x shift in one cell line is not evidence and the table has no way
    of saying so itself.
    """
    import csv
    import gzip
    import io
    import statistics

    found = _resolve(path)
    if found is None:
        return []
    target = (protein or "").strip().upper()
    if not target:
        return []

    fh = (io.TextIOWrapper(gzip.open(found, "rb"), encoding="utf-8",
                           errors="replace") if found.suffix == ".gz"
          else open(found, encoding="utf-8", errors="replace"))
    per_drug: dict = {}
    try:
        with fh:
            sample = fh.read(8192)
            fh.seek(0)
            delim = "\t" if sample.count("\t") >= sample.count(",") else ","
            for row in csv.DictReader(fh, delimiter=delim):
                row = {(k or "").strip().lower(): v for k, v in row.items()}
                if (row.get("protein") or row.get("gene") or
                        "").strip().upper() != target:
                    continue
                raw = row.get("log2fc") or row.get("logfc") or row.get("score")
                try:
                    fc = float(raw)
                except (TypeError, ValueError):
                    continue
                name = (row.get("drug") or "").strip()
                if not name:
                    continue
                rec = per_drug.setdefault(_key(name), {
                    "name": name, "fc": [], "cells": set(), "srcs": set()})
                rec["fc"].append(fc)
                if row.get("cell_line"):
                    rec["cells"].add(row["cell_line"].strip())
                if row.get("source"):
                    rec["srcs"].add(row["source"].strip())
    except OSError as exc:
        log.warning("could not read the protein-abundance table %s (%s)",
                    found, exc)
        return []

    out: list[ProteinHit] = []
    for rec in per_drug.values():
        n = len(rec["fc"])
        med = statistics.median(rec["fc"])
        if n < MIN_OBSERVATIONS or abs(med) < MIN_ABS_LOG2FC:
            continue
        moved = "up" if med > 0 else "down"
        if moved != want:
            continue
        hit = ProteinHit(
            drug=rec["name"], protein=target, direction=moved,
            median_log2fc=round(med, 3), n_observations=n,
            cell_lines=sorted(rec["cells"])[:6],
            sources=sorted(rec["srcs"])[:3])
        hit.notes.append(
            f"measured on the PROTEIN, not the transcript: median "
            f"{2 ** med:.1f}x over {n} observation(s)"
            + (f" in {', '.join(hit.cell_lines[:3])}" if hit.cell_lines
               else ""))
        if n < 3:
            hit.notes.append(
                "two observations is the floor this module accepts; treat "
                "it as a lead to repeat rather than a measurement")
        out.append(hit)
    out.sort(key=lambda h: (-h.n_observations, -abs(h.median_log2fc)))
    log.info("protein abundance for %s: %d compound(s) from %s",
             target, len(out), found.name)
    return out


def _key(name: str) -> str:
    from .signatures import breadth_key
    return breadth_key(name)


def summarise(hits: list, protein: str, *, path=None) -> str:
    """One paragraph for the report."""
    if _resolve(path) is None:
        # One line. The dataset comparison - which resources answer this
        # question and which only look as though they do - belongs in
        # docs/methods.md, not in the middle of a result.
        return (f"**No compound was checked against {protein} PROTEIN**, so "
                f"every abundance claim in this report is about the "
                f"transcript. A drug-to-protein table at "
                f"`data/drug_protein_abundance.tsv` would change that; see "
                f"[methods](methods.md) for which public datasets actually "
                f"answer this and which do not.")
    if not hits:
        return (f"**A drug-to-protein table was read and holds no compound "
                f"that raises {protein}.** That is a measured negative over "
                f"whatever the table covers, which on current public data "
                f"is at most a few hundred drugs - so it is not a statement "
                f"about the approved library as a whole.")
    return (
        f"**{len(hits)} compound(s) are measured to move {protein} protein "
        f"the needed way.** This is the only evidence in the report about "
        f"the protein rather than the transcript, which matters because a "
        f"stabiliser or a degradation blocker raises protein with the "
        f"transcript flat - invisible to every other route here. Coverage "
        f"is the limit: public drug-proteomics covers a small fraction of "
        f"the approved library, so a compound absent from it is unexamined "
        f"and not negative.")
