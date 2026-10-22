"""Does your gene change when a factor is lost? Read from perturbation data.

Why this exists
---------------
The route from a drug to your gene runs through a factor, and the step the
package could never see was whether your gene responds when that factor
changes in a cell. `outcome.py` filled the gap with a fixed 0.40, the same for
every factor, which told a reader nothing. This module replaces the placeholder
with what has been measured: experiments in which the factor was knocked out,
knocked down or overexpressed and the genes that changed were recorded.

Two Enrichr libraries are read, both from local files:

* `LINCS_L1000_CRISPR_KO_Consensus_Sigs`: for each gene knocked out by CRISPR
  in the L1000 cell lines, the consensus lists of genes that went up and down.
* `TF_Perturbations_Followed_by_Expression`: knockouts, knockdowns and
  overexpression of transcription factors from GEO, in human and mouse, with
  the genes up and down.

What a result says
------------------
For one gene and one factor, across every experiment on that factor:

* `measured response`: the gene is among the changed genes in at least one
  experiment. The report says how many, and which way. A gene that falls when
  the factor is lost is one the factor activates, and one that rises is one it
  represses. Overexpression is read the other way round. Experiments that
  disagree are reported as disagreeing.
* `tested, no response`: the factor was perturbed and the gene is in none of
  the lists. This is weak. Each list holds only the few hundred genes that
  changed most, so a gene can respond and not be on it.
* `not tested`: nothing in these files perturbs the factor. This is the usual
  case, and it is not evidence about the factor.

What it does not say
--------------------
It is not a probability and it does not replace one. The lists are set
membership with no magnitude. The cells are cancer lines, embryonic cells and
mouse tissue, not necessarily yours. A knockout is the loss of a factor, and a
drug changes its activity by some other amount. Mutation, truncation and fusion
experiments are not read, since they can raise or lower activity. A response
here says the factor can move the gene somewhere. It does not say it will in
your cells, or that a drug acting on the factor will do it.
"""

from __future__ import annotations

import logging
from pathlib import Path

log = logging.getLogger(__name__)

DEFAULT_KNOCKOUT_PATHS = ("data/knockout_lincs_crispr_ko.gmt",
                          "data/knockout_tf_perturbations.gmt")

SOURCE = ("Enrichr libraries LINCS_L1000_CRISPR_KO_Consensus_Sigs and "
          "TF_Perturbations_Followed_by_Expression, saved as "
          "data/knockout_lincs_crispr_ko.gmt and "
          "data/knockout_tf_perturbations.gmt. The launcher offers to fetch "
          "them: https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text"
          "&libraryName=<library>")

MEASURED = "measured response"
NO_RESPONSE = "tested, no response"
UNTESTED = "not tested"

# Words in a TF_Perturbations term that say the factor was lost or raised.
_LOSS = frozenset({"KO", "KD", "SIRNA", "SHRNA", "DEFICIENCY", "DEFICIENT",
                   "DEPLETION", "DEPLET", "INACTIVATION", "ABLATION",
                   "INHIBITION", "DELETION", "SILENCING", "KNOCKOUT", "DN",
                   "DOMNEG", "NULL"})
_GAIN = frozenset({"OE", "ACTIVATION", "ACTVATION", "INDUCTION", "INDUCED",
                   "HYPERACTIVE"})


def resolve_paths(paths=None) -> list:
    """The files that exist among `paths`, or among the default names."""
    return [str(p) for p in (paths or DEFAULT_KNOCKOUT_PATHS)
            if Path(p).is_file()]


def available(paths=None) -> bool:
    return bool(resolve_paths(paths))


def _parse(term: str) -> tuple | None:
    """`(factor, experiment id, kind, direction, source)` for a library term.

    `kind` is "loss" or "gain" and `direction` is "up" or "down" as the
    library states it. Returns None for a term that is not read: a mutation,
    a fusion, or anything that is not clearly a loss or a gain.
    """
    toks = term.split()
    if len(toks) >= 2 and toks[-1].upper() in ("UP", "DOWN"):
        direction = toks[-1].lower()
        exp = " ".join(toks[:-1])
        if len(toks) == 2:
            # "GENE Up": a CRISPR knockout of the gene, consensus over lines.
            return toks[0].upper(), exp, "loss", direction, "LINCS CRISPR KO"
        kind_word = toks[1].upper()
        if kind_word in _LOSS:
            kind = "loss"
        elif kind_word in _GAIN:
            kind = "gain"
        else:
            return None
        return toks[0].upper(), exp, kind, direction, "GEO TF perturbation"
    return None


def _vote(kind: str, direction: str) -> str:
    """"activates" when the gene falls on loss or rises on gain, else
    "represses"."""
    lost_down = kind == "loss" and direction == "down"
    gained_up = kind == "gain" and direction == "up"
    return "activates" if (lost_down or gained_up) else "represses"


def assess_many(gene: str, factors, paths=None) -> dict:
    """`{factor: result}` for one gene against several factors.

    Streams the files once and keeps nothing but counts, since the libraries
    hold about three million gene memberships. Each result has `class`,
    `n_experiments` (distinct perturbations of the factor), `n_responded`,
    `activates`, `represses`, `sources`, `species` and `experiments` (the first
    few that responded, for reading).
    """
    gene = (gene or "").strip().upper()
    want = {f.strip().upper() for f in factors if f}
    seen: dict = {f: {} for f in want}
    for path in resolve_paths(paths):
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                term, _, rest = line.rstrip("\r\n").partition("\t")
                parsed = _parse(term)
                if parsed is None or parsed[0] not in want:
                    continue
                factor, exp, kind, direction, source = parsed
                rec = seen[factor].setdefault(
                    exp, {"kind": kind, "source": source, "votes": set()})
                # rest begins with the empty description field.
                if gene in rest.split("\t"):
                    rec["votes"].add(_vote(kind, direction))
    return {f: _summarise(seen[f]) for f in want}


def regulators_of(gene: str, paths=None) -> dict:
    """`{factor: {activates, represses, sign}}` for every factor whose
    perturbation changed `gene`.

    The reverse of `assess_many`: one pass over the files, keeping the
    experiments whose lists name the gene. `activates` and `represses` count
    experiments, not list memberships, and `sign` is +1 when more experiments
    say the factor raises the gene, -1 when more say it lowers it, and 0 on a
    tie. Every gene that was perturbed appears as a factor here, so a caller
    that wants transcription factors only has to filter.
    """
    gene = (gene or "").strip().upper()
    seen: dict = {}
    if not gene:
        return {}
    for path in resolve_paths(paths):
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                term, _, rest = line.rstrip("\r\n").partition("\t")
                parsed = _parse(term)
                if parsed is None:
                    continue
                factor, exp, kind, direction, _source = parsed
                if gene not in rest.split("\t"):
                    continue
                seen.setdefault(factor, {}).setdefault(exp, set()).add(
                    _vote(kind, direction))
    out = {}
    for factor, exps in seen.items():
        act = sum(1 for v in exps.values() if "activates" in v)
        rep = sum(1 for v in exps.values() if "represses" in v)
        out[factor] = {"activates": act, "represses": rep,
                       "sign": (1 if act > rep else -1 if rep > act else 0)}
    return out


def _species(exp: str, source: str) -> str:
    if source == "LINCS CRISPR KO":
        return "human cell lines"
    up = exp.upper().split()
    return "mouse" if "MOUSE" in up else "human" if "HUMAN" in up else "other"


def _summarise(experiments: dict) -> dict:
    n = len(experiments)
    out = {"n_experiments": n, "n_responded": 0, "activates": 0,
           "represses": 0, "sources": sorted({e["source"]
                                              for e in experiments.values()}),
           "species": sorted({_species(k, v["source"])
                              for k, v in experiments.items()}),
           "experiments": []}
    for exp, rec in sorted(experiments.items()):
        if not rec["votes"]:
            continue
        out["n_responded"] += 1
        for v in rec["votes"]:
            out[v] += 1
        if len(out["experiments"]) < 5:
            out["experiments"].append(
                {"experiment": exp, "reads": "/".join(sorted(rec["votes"]))})
    if n == 0:
        out["class"] = UNTESTED
    elif out["n_responded"]:
        out["class"] = MEASURED
    else:
        out["class"] = NO_RESPONSE
    return out


def direction_word(res: dict) -> str | None:
    """Which way the factor moves the gene, as measured, or None."""
    a, r = res.get("activates", 0), res.get("represses", 0)
    if a and not r:
        return "activates"
    if r and not a:
        return "represses"
    if a and r:
        return "disagrees"
    return None


def cell(res: dict | None) -> str:
    """A short phrase for a table cell."""
    if not res:
        return "not checked"
    cls = res.get("class")
    n, k = res.get("n_experiments", 0), res.get("n_responded", 0)
    if cls == MEASURED:
        word = direction_word(res)
        what = {"activates": "falls when the factor is lost",
                "represses": "rises when the factor is lost",
                "disagrees": "moves both ways"}[word]
        return f"gene {what}, {k} of {n}"
    if cls == NO_RESPONSE:
        return f"not among top changes, 0 of {n}"
    return "—"


def role_conflict(res: dict | None, role: str | None) -> str | None:
    """A sentence when the measured direction contradicts the factor's role.

    The role at the promoter comes from annotation. When knockout data say the
    opposite, one of them is wrong for this gene, and the reader should know
    before choosing whether to activate or block the factor.
    """
    word = direction_word(res or {})
    if word in (None, "disagrees") or not role:
        return None
    r = role.lower()
    if word == "activates" and r.startswith("repress"):
        return ("annotated as a repressor, yet the gene falls when it is "
                "lost in the knockout data")
    if word == "represses" and r.startswith("activ"):
        return ("annotated as an activator, yet the gene rises when it is "
                "lost in the knockout data")
    return None
