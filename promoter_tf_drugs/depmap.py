"""What happens to a cell when a factor is removed, measured not curated.

The question this answers
-------------------------
`regulon.py` counts how many genes a factor controls, from literature
curation, and that count measures how much anyone has published as much as
what the factor does. For most factors it is absent entirely: a real CCR5
run left 83 of them at "not counted".

DepMap asks a different question with a measurement. CRISPR knockout of
every gene in roughly 1,100 cancer cell lines, scored by how much the cells
suffer. It does not say what the factor regulates. It says what the cell
cannot do without it, which is the collateral question the regulon count is
standing in for, answered directly.

How to read a Chronos score
---------------------------
Zero means knocking the gene out changed nothing. Minus one is the median
of known common-essential genes, so it is the natural unit: a score near -1
means losing this gene is about as bad as losing a ribosomal protein. A
line counts as dependent below -0.5.

Three readouts, and the third is the useful one:

  mean effect        how bad losing it is, averaged over cell lines
  fraction dependent in how many lines it matters at all
  selectivity        whether that is a few lines or all of them

A factor essential everywhere is a poor target for anything chronic. A
factor essential nowhere is safer to perturb and tells you little else. A
factor essential in a narrow set of lines is the interesting case, because
the dependency has a context and your tissue may or may not be in it.

What this is not
----------------
**Cancer cell lines dividing in a dish.** Essentiality there is a statement
about proliferating transformed cells, not about toxicity in differentiated
tissue in a patient. A factor that matters only in a post-mitotic neuron
will look dispensable here, and that is a property of the panel, not of the
factor.

**Knockout, not inhibition.** The protein is gone, including whatever
scaffolding it does that no ligand would touch, and gone completely rather
than partially. Plenty of essential genes are drugged successfully at
partial occupancy.

**Silent about your gene.** Nothing here says the factor regulates your
gene, or which direction a drug would move it. It bounds the cost of acting
on the factor; the benefit has to come from somewhere else in the report.

So a dependency finding is a reason to be careful, and the absence of one is
not permission.
"""

from __future__ import annotations

import csv
import logging
import statistics
from dataclasses import dataclass, field
from pathlib import Path

log = logging.getLogger(__name__)

# Picked up automatically, in order. The DepMap release names its files
# consistently, so the common spellings are all accepted.
DEFAULT_EFFECT_PATHS = (
    "data/CRISPRGeneEffect.csv", "data/crispr_gene_effect.csv",
    "data/depmap_gene_effect.csv",
)
DEFAULT_MODEL_PATHS = ("data/Model.csv", "data/model.csv")
DEFAULT_ESSENTIAL_PATHS = (
    "data/CRISPRInferredCommonEssentials.csv",
    "data/common_essentials.csv",
)

SOURCE = ("DepMap: https://depmap.org/portal/data_page/ - download "
          "CRISPRGeneEffect.csv from the current release into data/. "
          "Model.csv adds lineage, and is optional.")

# DepMap's own conventions, not invented here.
DEPENDENT_BELOW = -0.5      # a line counts as dependent below this
COMMON_ESSENTIAL_MEDIAN = -1.0
PAN_ESSENTIAL_FRACTION = 0.90
# Below this fraction a dependency has a context worth identifying rather
# than being a general requirement for staying alive.
SELECTIVE_BELOW = 0.25


@dataclass
class Dependency:
    """One factor's CRISPR knockout phenotype across the DepMap panel."""
    gene: str
    mean_effect: float | None = None
    median_effect: float | None = None
    min_effect: float | None = None
    n_lines: int = 0
    n_dependent: int = 0
    fraction_dependent: float = 0.0
    spread: float | None = None             # SD across lines
    verdict: str = "not measured"           # see _verdict
    selective: bool = False
    top_lineages: list = field(default_factory=list)
    burden: float | None = None             # 0..1, HIGHER = more collateral
    # DepMap's own common-essential call, when that file is on disk. None
    # means it was not consulted, which is not the same as False: a gene
    # absent from a list nobody loaded has not been cleared.
    depmap_common_essential: bool | None = None
    caveats: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return dict(self.__dict__)


def _verdict(mean: float, frac: float) -> str:
    if frac >= PAN_ESSENTIAL_FRACTION:
        return "essential in nearly every line"
    if frac >= 0.5:
        return "essential in most lines"
    if frac >= SELECTIVE_BELOW:
        return "essential in a substantial minority"
    if frac >= 0.05:
        return "essential in a minority of lines"
    if frac > 0.0:
        return "essential in a handful of lines"
    return "dispensable across the panel"


# ------------------------------------------------------------ file reading
def _symbol(column: str) -> str:
    """`A1BG (1)` -> `A1BG`. DepMap suffixes the Entrez id onto every
    column name, so a plain symbol lookup finds nothing without this."""
    return column.split("(")[0].strip().upper()


def resolve(explicit=None, candidates=DEFAULT_EFFECT_PATHS):
    """First existing path, or None. Never raises - a missing DepMap file
    costs you this section and nothing else."""
    for c in ([explicit] if explicit else list(candidates)):
        if c and Path(c).exists():
            return str(c)
    return None


def load_gene_effect(path, genes) -> dict:
    """`{SYMBOL: [effect per cell line]}` for the requested genes only.

    CRISPRGeneEffect.csv is cell lines by genes and about 500 MB, so it is
    streamed once and only the requested columns are kept. Loading the whole
    matrix to read eight columns is the obvious way to write this and the
    reason it would be unusable.

    Also returns the cell-line ids, in row order, under the key
    `__models__`, so a lineage lookup can line them up afterwards.
    """
    want = {g.strip().upper() for g in genes if g}
    if not want:
        return {}

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{path} not found. {SOURCE}")

    out: dict = {}
    models: list = []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh)
        try:
            header = next(reader)
        except StopIteration:
            return {}

        idx: dict = {}
        for i, col in enumerate(header[1:], start=1):
            sym = _symbol(col)
            if sym in want and sym not in idx:
                idx[sym] = i
        if not idx:
            log.info("none of the %d requested gene(s) appear in %s",
                     len(want), path.name)
            return {}

        for sym in idx:
            out[sym] = []
        for row in reader:
            if not row:
                continue
            models.append(row[0])
            for sym, i in idx.items():
                if i >= len(row):
                    out[sym].append(None)
                    continue
                v = row[i].strip()
                try:
                    out[sym].append(float(v))
                except ValueError:
                    out[sym].append(None)     # DepMap leaves blanks

    out["__models__"] = models
    log.info("DepMap: %d of %d requested gene(s) found across %d cell line(s)",
             len(idx), len(want), len(models))
    return out


def load_lineages(path=None, *, beside=None) -> dict:
    """`{ModelID: lineage}` from Model.csv. Empty when absent.

    `beside` is the gene-effect file. DepMap ships Model.csv in the same
    download, so a user who points --depmap-file at their download folder
    expects the lineages to come with it; looking only in data/ silently
    dropped the whole lineage column.
    """
    cands = list(DEFAULT_MODEL_PATHS)
    if beside:
        here = Path(beside).parent
        cands = [str(here / n) for n in ("Model.csv", "model.csv")] + cands
    p = resolve(path, cands)
    if not p:
        return {}
    out: dict = {}
    with open(p, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            mid = (row.get("ModelID") or row.get("DepMap_ID") or "").strip()
            lin = (row.get("OncotreeLineage") or row.get("lineage")
                   or row.get("primary_disease") or "").strip()
            if mid and lin:
                out[mid] = lin
    log.info("DepMap lineages: %d cell line(s)", len(out))
    return out


def load_common_essentials(path=None) -> set:
    """DepMap's own inferred common-essential list, if it is on disk."""
    p = resolve(path, DEFAULT_ESSENTIAL_PATHS)
    if not p:
        return set()
    # The header is a row like any other to csv.reader, and the release
    # has used "Essentials", "Gene" and "gene" across versions - so
    # discarding one literal let the others through as a gene called
    # ESSENTIALS. Skip the first row when it is not a gene entry.
    _HEADERS = {"GENE", "GENES", "ESSENTIALS", "ESSENTIAL",
                "COMMON_ESSENTIALS", "SYMBOL"}
    out: set = set()
    with open(p, newline="", encoding="utf-8", errors="replace") as fh:
        for i, row in enumerate(csv.reader(fh)):
            if not row or not row[0].strip():
                continue
            sym = _symbol(row[0])
            if i == 0 and sym in _HEADERS:
                continue
            out.add(sym)
    out -= _HEADERS
    log.info("DepMap common essentials: %d gene(s)", len(out))
    return out


# -------------------------------------------------------------- assessment
def assess(gene: str, effects: list, *, lineages: list | None = None,
           common_essentials: set | None = None) -> Dependency:
    """Summarise one gene's knockout phenotype.

    `effects` is the per-cell-line Chronos score, with None for lines where
    the gene was not scored. `lineages` is the matching lineage per line
    when Model.csv was available.

    `common_essentials` is DepMap's own inferred list. It is a second
    opinion rather than a duplicate: that call is made over the whole
    release, while the fraction computed here sees only the lines in the
    matrix handed to this function, which may be a subset. Where the two
    disagree, both are reported - a disagreement is information about how
    much the local panel can support, not a number to pick between.
    """
    d = Dependency(gene=gene.upper())
    vals = [v for v in (effects or []) if v is not None]
    if common_essentials is not None:
        d.depmap_common_essential = gene.upper() in common_essentials
    if len(vals) < 10:
        d.caveats.append(
            f"only {len(vals)} cell line(s) carry a score for this gene, "
            f"which is too few to summarise - treat this as not measured")
        if d.depmap_common_essential:
            d.caveats.append(
                "DepMap's own common-essential list does include this gene. "
                "Too few lines here to confirm it, and that list is the "
                "better-powered call, so treat the factor as essential")
        return d

    d.n_lines = len(vals)
    d.mean_effect = round(statistics.fmean(vals), 3)
    d.median_effect = round(statistics.median(vals), 3)
    d.min_effect = round(min(vals), 3)
    d.spread = round(statistics.pstdev(vals), 3)
    dependent = [v for v in vals if v < DEPENDENT_BELOW]
    d.n_dependent = len(dependent)
    d.fraction_dependent = round(len(dependent) / float(len(vals)), 3)
    d.verdict = _verdict(d.mean_effect, d.fraction_dependent)
    d.selective = 0.0 < d.fraction_dependent < SELECTIVE_BELOW

    # Burden on the same convention systemic.clinical_burden uses: higher
    # means more collateral. Anchored on DepMap's own scale, where -1 is
    # the median common essential.
    d.burden = round(min(1.0, max(0.0, -d.mean_effect / 1.0)), 3)

    if lineages:
        per: dict = {}
        for v, lin in zip(effects or [], lineages):
            if v is None or not lin:
                continue
            per.setdefault(lin, []).append(v)
        ranked = sorted(
            ((lin, statistics.fmean(vs), len(vs))
             for lin, vs in per.items() if len(vs) >= 5),
            key=lambda t: t[1])
        d.top_lineages = [{"lineage": lin, "mean_effect": round(m, 3),
                           "n_lines": n} for lin, m, n in ranked[:5]]

    # Caveats that depend on what was found.
    if d.fraction_dependent >= PAN_ESSENTIAL_FRACTION:
        d.caveats.append(
            "removing this factor harms nearly every line tested, so a drug "
            "that strongly inhibits it is likely cytotoxic rather than "
            "selectively regulatory. That is a reason for caution, not a "
            "veto: essential genes are drugged successfully at partial "
            "occupancy, and knockout is not inhibition")
    elif d.selective:
        d.caveats.append(
            f"the dependency is confined to {d.n_dependent} of {d.n_lines} "
            f"lines, so it has a context. Whether your tissue is in that "
            f"context is the question, and this panel cannot answer it")
    elif d.fraction_dependent == 0.0:
        d.caveats.append(
            "no line in the panel needed this factor to grow. That bounds "
            "the proliferative cost and says nothing about a differentiated "
            "tissue, where the factor's actual job may lie")

    if d.depmap_common_essential is True and \
            d.fraction_dependent < PAN_ESSENTIAL_FRACTION:
        d.caveats.append(
            f"DepMap calls this gene common-essential, while only "
            f"{d.fraction_dependent:.0%} of the {d.n_lines} line(s) read here "
            f"fall below the dependency threshold. That list is inferred "
            f"over the whole release, so where the two differ it is the "
            f"better-powered call - the local fraction is reported, not "
            f"preferred")
    elif d.depmap_common_essential is False and \
            d.fraction_dependent >= PAN_ESSENTIAL_FRACTION:
        d.caveats.append(
            f"{d.fraction_dependent:.0%} of lines here depend on this gene, "
            f"yet DepMap does not list it as common-essential. Its method "
            f"requires consistency across the full panel, so this pattern "
            f"usually means the lines read here are not representative")

    d.caveats.append(
        "DepMap is cancer cell lines dividing in vitro. A factor whose role "
        "is in a post-mitotic or differentiated cell looks dispensable here "
        "for reasons that have nothing to do with its importance")
    return d


def assess_many(genes, *, effect_path=None, model_path=None,
                lineages: dict | None = None, essential_path=None) -> dict:
    """`{SYMBOL: Dependency}` for every gene DepMap scored.

    Genes absent from the file are omitted rather than returned as
    dispensable, because "not in the release" and "knockout did nothing"
    are different statements and only one of them is reassuring.
    """
    p = resolve(effect_path)
    if not p:
        log.info("no DepMap file found; dependency is not assessed. %s",
                 SOURCE)
        return {}

    data = load_gene_effect(p, genes)
    if not data:
        return {}
    models = data.pop("__models__", [])

    lin_map = (lineages if lineages is not None
               else load_lineages(model_path, beside=p))
    per_line = [lin_map.get(m) for m in models] if lin_map else None

    # Loaded once for the whole batch rather than per gene: it is a file
    # read, and `load_common_essentials` returns an empty set when there is
    # none. Empty must stay None here, because an empty set would make
    # every gene read as "DepMap says not essential".
    ess = load_common_essentials(essential_path) or None

    out = {}
    for sym, effects in data.items():
        d = assess(sym, effects, lineages=per_line, common_essentials=ess)
        if d.n_lines:
            out[sym] = d
    return out


def summarise(deps: dict, *, requested: int = 0) -> str:
    """Prose for the report."""
    if not deps:
        return ("No DepMap data was loaded, so the cost of acting on these "
                "factors is unmeasured. " + SOURCE)

    pan = [d for d in deps.values()
           if d.fraction_dependent >= PAN_ESSENTIAL_FRACTION]
    sel = [d for d in deps.values() if d.selective]
    free = [d for d in deps.values() if d.fraction_dependent == 0.0]

    def _n(items, singular, plural):
        names = ", ".join(d.gene for d in items[:5])
        return (f"{len(items)} ({names}) {singular}" if len(items) == 1
                else f"{len(items)} ({names}) {plural}")

    bits = [f"{len(deps)} of {requested or len(deps)} factor(s) were scored "
            f"in DepMap."]
    if pan:
        bits.append(
            _n(pan, "is needed", "are needed")
            + " by nearly every cell line tested, so strong inhibition is "
              "likely cytotoxic.")
    if sel:
        bits.append(
            _n(sel, "shows", "show")
            + " a dependency confined to a minority of lines, which is the "
              "informative case: the requirement has a context, and whether "
              "your tissue shares it is the question to ask next.")
    if free:
        bits.append(
            _n(free, "was", "were")
            + " dispensable throughout. That bounds the proliferative cost "
              "of removing them and says nothing about a differentiated "
              "tissue.")
    return " ".join(bits)
