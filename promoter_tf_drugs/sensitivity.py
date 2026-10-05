"""Which approved drugs behave as though they act on your factor.

The question this answers
-------------------------
Every other route in this package starts from something somebody wrote
down: a curated ligand annotation, a motif, a structure, a signed edge in
a network. That shared dependence is a blind spot. A drug nobody has ever
linked to your factor cannot be found by any of them, however obvious the
link would be if anyone looked.

This asks a question that needs no annotation at all. DepMap measures what
happens to roughly 1,100 cancer cell lines when each gene is knocked out.
PRISM measures what happens to the same lines under roughly 6,000
compounds, most of them approved. Both are matrices over the same cell
lines, so they can be compared column by column.

If knocking out your factor hurts the same cell lines that a drug hurts,
the drug is behaving as though it removes the factor. That is a measured
phenotypic argument, and it is made without any record of the compound
binding anything.

Two readouts
------------
**Phenocopy.** Correlate a drug's viability profile against the factor's
CRISPR knockout profile. Both scales run negative for harm - Chronos is
negative where knockout hurts, PRISM log fold change is negative where the
drug kills - so a POSITIVE correlation means the drug and the knockout
hurt the same lines.

**Expression dependence.** Correlate a drug's viability profile against
the gene's expression across the same lines. Here the scales run opposite
ways: expression is higher where there is more transcript, viability is
more negative where the drug kills more, so sensitivity rising with
expression appears as a NEGATIVE correlation. The two readouts therefore
have opposite expected signs, which is a property of the source scales and
not of the biology.

The confounder that dominates both
-----------------------------------
Cell lines are not equally fragile. Some die easily under almost anything,
and those same lines suffer most knockouts. That shared sensitivity puts a
large positive correlation between almost any drug and almost any gene,
and reading it as a specific relationship is the main way this analysis
goes wrong.

So the number reported is a partial correlation, controlling for each
line's mean response across the whole drug panel. The raw correlation is
reported beside it, and the gap between them is usually large. Ranking on
the raw value would mostly rank cell-line fragility.

What this does not establish
----------------------------
**Not a target.** A drug whose profile matches your factor's knockout may
act on the factor, or anywhere else in a pathway the factor is required
for, or on something that happens to be needed by the same lines. The
correlation is agnostic between those and the first is the least likely of
the three for any given hit.

**A proliferation readout.** Both panels score growth of transformed cells
in a dish. A factor whose job is in a differentiated or post-mitotic cell
looks silent here, for reasons that have nothing to do with its
importance, and so does every drug that acts on it.

**Nothing about your gene.** This is about the FACTOR, not the promoter it
binds. Whether moving the factor moves your gene is the same open question
every other TF route in this report leaves open.

**Correlation across lines, not a dose response.** A strong correlation
built on a narrow range of either profile is fragile, so the number of
lines behind each estimate is reported and small overlaps are refused
rather than scored.
"""

from __future__ import annotations

import csv
import logging
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

log = logging.getLogger(__name__)

# PRISM ships under several names across releases; all the common
# spellings are accepted so a download needs no renaming.
DEFAULT_PRISM_PATHS = (
    "data/Repurposing_Public_24Q2_LFC_COLLAPSED.csv",
    "data/Repurposing_Public_LFC.csv",
    "data/prism_lfc.csv",
    "data/PRISM_Repurposing_Public.csv",
)

DEFAULT_COMPOUND_PATHS = (
    "data/Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv",
    "data/Repurposing_Public_Compound_List.csv",
    "data/prism_compounds.csv",
)

DEFAULT_EXPRESSION_PATHS = (
    "data/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
    "data/CCLE_expression.csv",
    "data/ccle_expression.csv",
)

SOURCE = (
    "PRISM drug sensitivity: https://depmap.org/portal/data_page/ - "
    "download the Repurposing Public LFC matrix into data/, and the "
    "compound list beside it to get drug names instead of BRD ids. "
    "Expression is optional and comes from the same page.")

# Below this many shared cell lines a correlation is not estimated. Ten
# points can produce r = 0.9 from noise, and this analysis exists to find
# a handful of hits among thousands of tests.
MIN_LINES = 25

# Reported as a hit. Deliberately modest: after partialling out line
# fragility the surviving correlations are small, and the multiple-testing
# correction is what decides whether one is real.
MIN_ABS_R = 0.20
MAX_Q = 0.10

# A q-value floor, so the most significant result cannot round to zero and
# read as missing through a falsy test somewhere downstream.
Q_FLOOR = 1e-6


@dataclass
class DrugCorrelation:
    """One drug's relationship to one gene's profile across cell lines."""
    drug: str
    column: str                      # the matrix column it came from
    r: float | None = None           # partial, controlling line fragility
    r_raw: float | None = None       # before that control
    n_lines: int = 0
    p_value: float | None = None
    q_value: float | None = None
    readout: str = "phenocopy"       # or "expression"
    direction: str = "unclear"       # "phenocopies" | "opposes" | "unclear"

    def to_dict(self) -> dict:
        return dict(self.__dict__)


@dataclass
class SensitivityResult:
    """Everything one gene's scan produced, including why it found nothing."""
    gene: str
    assessed: bool = False
    readout: str = "phenocopy"
    n_drugs_tested: int = 0
    n_lines: int = 0
    hits: list = field(default_factory=list)
    caveats: list = field(default_factory=list)
    error: str | None = None

    def to_dict(self) -> dict:
        d = dict(self.__dict__)
        d["hits"] = [h.to_dict() if isinstance(h, DrugCorrelation) else h
                     for h in self.hits]
        return d


def resolve(explicit=None, candidates=DEFAULT_PRISM_PATHS):
    """First existing path, or None. Never raises: a missing PRISM file
    costs this section and nothing else in the run."""
    for c in ([explicit] if explicit else list(candidates)):
        if c and Path(c).exists():
            return str(c)
    return None


def available() -> dict:
    """Which of the three files this module can use are on disk."""
    return {
        "prism": resolve(None, DEFAULT_PRISM_PATHS),
        "compounds": resolve(None, DEFAULT_COMPOUND_PATHS),
        "expression": resolve(None, DEFAULT_EXPRESSION_PATHS),
    }


def load_compound_names(path=None) -> dict:
    """`{column id: drug name}` from the PRISM compound list.

    Without it the matrix columns are BRD identifiers, which are correct
    and unreadable. A missing file leaves the identifier in place rather
    than dropping the drug.
    """
    p = resolve(path, DEFAULT_COMPOUND_PATHS)
    if not p:
        return {}
    out: dict = {}
    with open(p, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh)
        cols = {(c or "").strip().lower(): c for c in (reader.fieldnames or [])}

        def col(*names):
            for n in names:
                if n in cols:
                    return cols[n]
            return None

        c_id = col("iduplicate", "broad_id", "column_name", "sample_id",
                   "id", "drug_id")
        c_name = col("name", "drug_name", "pert_iname", "compound")
        if not (c_id and c_name):
            log.warning("%s has no id/name column pair; PRISM columns stay "
                        "as identifiers", Path(p).name)
            return {}
        for row in reader:
            key = (row.get(c_id) or "").strip()
            name = (row.get(c_name) or "").strip()
            if key and name:
                out.setdefault(key, name)
    log.info("PRISM compound names: %d entries from %s", len(out),
             Path(p).name)
    return out


def load_prism(path=None) -> dict:
    """The PRISM matrix as `{models, columns, values}`.

    `values` is a float array, cell lines by compounds, with NaN where a
    compound was not screened in a line. The whole matrix is kept because
    every column is a candidate: unlike the gene-effect loader, there is
    no short list of columns to select.

    Returns `{}` when no file is on disk, which the caller must treat as
    not measured rather than as no drug correlating.
    """
    p = resolve(path, DEFAULT_PRISM_PATHS)
    if not p:
        log.info("no PRISM file on disk; drug sensitivity is not assessed. "
                 "%s", SOURCE.splitlines()[0])
        return {}

    models: list = []
    rows: list = []
    with open(p, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh)
        try:
            header = next(reader)
        except StopIteration:
            return {}
        columns = [c.strip() for c in header[1:]]
        if not columns:
            return {}
        for row in reader:
            if not row:
                continue
            models.append(row[0].strip())
            vals = np.full(len(columns), np.nan)
            for i, raw in enumerate(row[1:len(columns) + 1]):
                raw = raw.strip()
                if not raw:
                    continue
                try:
                    vals[i] = float(raw)
                except ValueError:
                    continue
            rows.append(vals)

    if not rows:
        return {}
    values = np.vstack(rows)
    log.info("PRISM: %d cell line(s) x %d compound(s) from %s",
             values.shape[0], values.shape[1], Path(p).name)
    return {"models": models, "columns": columns, "values": values}


def load_expression(path, gene: str) -> dict:
    """`{cell line: expression}` for one gene, streamed column-selected.

    The expression matrix is cell lines by genes and large, and one gene
    is wanted, so it is read the way `depmap.load_gene_effect` reads its
    own: one pass, one column kept.
    """
    p = resolve(path, DEFAULT_EXPRESSION_PATHS)
    if not p:
        return {}
    want = (gene or "").strip().upper()
    if not want:
        return {}

    out: dict = {}
    with open(p, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh)
        try:
            header = next(reader)
        except StopIteration:
            return {}
        idx = None
        for i, col in enumerate(header[1:], start=1):
            # The expression matrix uses the same `SYMBOL (entrez)` form.
            if col.split("(")[0].strip().upper() == want:
                idx = i
                break
        if idx is None:
            log.info("%s is not a column in %s", want, Path(p).name)
            return {}
        for row in reader:
            if not row or idx >= len(row):
                continue
            raw = row[idx].strip()
            if not raw:
                continue
            try:
                out[row[0].strip()] = float(raw)
            except ValueError:
                continue
    log.info("expression for %s: %d cell line(s)", want, len(out))
    return out


def _bh(pvals: list) -> list:
    """Benjamini-Hochberg q-values, order preserved, floored."""
    vals = [(i, p) for i, p in enumerate(pvals) if p is not None]
    if not vals:
        return [None] * len(pvals)
    vals.sort(key=lambda t: t[1])
    m = len(vals)
    out = [None] * len(pvals)
    running = 1.0
    for rank in range(m, 0, -1):
        i, p = vals[rank - 1]
        running = min(running, p * m / rank)
        out[i] = max(Q_FLOOR, round(min(1.0, running), 6))
    return out


def _p_from_r(r: float, n: int, n_controls: int = 0) -> float | None:
    """Two-sided p for a correlation, by the usual t transform.

    `n_controls` is how many variables were partialled out, which costs a
    degree of freedom each. Ignoring that would make a partial correlation
    look very slightly more significant than it is.
    """
    df = n - 2 - n_controls
    if df <= 0:
        return None
    r = max(-0.999999, min(0.999999, float(r)))
    t = abs(r) * math.sqrt(df / (1.0 - r * r))
    # Student t survival function without scipy, by the incomplete beta
    # identity: P(|T| > t) = I_{df/(df+t^2)}(df/2, 1/2)
    x = df / (df + t * t)
    return max(0.0, min(1.0, _betainc_half(x, df / 2.0)))


def _betainc_half(x: float, a: float) -> float:
    """Regularised incomplete beta I_x(a, 1/2), by continued fraction.

    Only this one shape is needed, for the t distribution's tail. NumPy
    has no incomplete beta and scipy is not a dependency here, so it is
    computed directly rather than adding one.
    """
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    b = 0.5
    lbeta = (math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b))
    front = math.exp(a * math.log(x) + b * math.log(1.0 - x) - lbeta) / a
    # Lentz's algorithm for the continued fraction.
    f, c, d = 1.0, 1.0, 0.0
    for i in range(0, 300):
        m = i // 2
        if i == 0:
            num = 1.0
        elif i % 2 == 0:
            num = (m * (b - m) * x) / ((a + 2 * m - 1) * (a + 2 * m))
        else:
            num = -((a + m) * (a + b + m) * x) / ((a + 2 * m) * (a + 2 * m + 1))
        d = 1.0 + num * d
        if abs(d) < 1e-30:
            d = 1e-30
        d = 1.0 / d
        c = 1.0 + num / c
        if abs(c) < 1e-30:
            c = 1e-30
        f *= c * d
        if abs(1.0 - c * d) < 1e-10:
            break
    return min(1.0, max(0.0, front * (f - 1.0)))


def _pairwise_r(matrix: np.ndarray, target: np.ndarray) -> tuple:
    """Pearson r of every column against `target`, pairwise complete.

    Returns (r, n) arrays. Columns with fewer than two shared points, or
    no variance on either side, come back as NaN rather than as zero - a
    correlation that could not be computed is not a correlation of none.
    """
    n_cols = matrix.shape[1]
    out_r = np.full(n_cols, np.nan)
    out_n = np.zeros(n_cols, dtype=int)

    t_ok = np.isfinite(target)
    for j in range(n_cols):
        col = matrix[:, j]
        ok = t_ok & np.isfinite(col)
        n = int(ok.sum())
        out_n[j] = n
        if n < 3:
            continue
        a = col[ok]
        b = target[ok]
        sa, sb = a.std(), b.std()
        if sa == 0.0 or sb == 0.0:
            continue
        out_r[j] = float(((a - a.mean()) * (b - b.mean())).mean() / (sa * sb))
    return out_r, out_n


def _partial(r_xy: np.ndarray, r_xz: np.ndarray,
             r_yz: np.ndarray | float) -> np.ndarray:
    """Partial correlation of x and y controlling for z.

        r_xy.z = (r_xy - r_xz r_yz) / sqrt((1 - r_xz^2)(1 - r_yz^2))

    Used to take out the one thing that otherwise dominates: some cell
    lines die under almost anything, and those same lines suffer most
    knockouts, so almost every drug correlates with almost every gene.
    """
    denom = np.sqrt(np.clip(1.0 - r_xz ** 2, 1e-12, None)
                    * np.clip(1.0 - np.asarray(r_yz) ** 2, 1e-12, None))
    return (r_xy - r_xz * r_yz) / denom


def correlate(prism: dict, target: dict, *, readout: str = "phenocopy",
              gene: str = "?", names: dict | None = None,
              min_lines: int = MIN_LINES, min_abs_r: float = MIN_ABS_R,
              max_q: float = MAX_Q, top_n: int = 25) -> SensitivityResult:
    """Correlate every drug's viability profile against one gene profile.

    `target` maps cell line id to the gene's value: a Chronos score for
    the phenocopy readout, an expression value for the other. Lines absent
    from either side are dropped pairwise, so one sparsely screened drug
    does not shrink the rest.
    """
    res = SensitivityResult(gene=(gene or "?").upper(), readout=readout)
    if not prism or not prism.get("columns"):
        res.error = ("no PRISM matrix loaded, so no drug was tested. That "
                     "is a missing file, not an absence of correlated "
                     "drugs. " + SOURCE.splitlines()[0])
        return res
    if not target:
        res.error = (f"{res.gene} has no profile to correlate against - it "
                     f"is not a column in the file this readout needs. Not "
                     f"assessed.")
        return res

    models = prism["models"]
    values = prism["values"]
    vec = np.array([target.get(m, np.nan) for m in models], dtype=float)
    shared = int((np.isfinite(vec) & np.isfinite(values).any(axis=1)).sum())
    res.n_lines = shared
    if shared < min_lines:
        res.error = (
            f"only {shared} cell line(s) are in both the PRISM matrix and "
            f"the {readout} profile for {res.gene}, below the {min_lines} "
            f"needed. A correlation over fewer points is noise with a "
            f"decimal place, so none is reported.")
        return res

    # The confounder: how hard each line is to kill, averaged over the
    # whole panel. Taken out of every correlation below.
    with np.errstate(invalid="ignore"):
        fragility = np.nanmean(values, axis=1)

    r_xy, n_xy = _pairwise_r(values, vec)
    r_xz, _ = _pairwise_r(values, fragility)
    ok_z = np.isfinite(fragility) & np.isfinite(vec)
    if int(ok_z.sum()) >= 3 and fragility[ok_z].std() > 0 \
            and vec[ok_z].std() > 0:
        a, b = fragility[ok_z], vec[ok_z]
        r_yz = float(((a - a.mean()) * (b - b.mean())).mean()
                     / (a.std() * b.std()))
    else:
        r_yz = 0.0

    r_part = _partial(r_xy, np.nan_to_num(r_xz), r_yz)

    cands, pvals = [], []
    for j, col in enumerate(prism["columns"]):
        if n_xy[j] < min_lines or not np.isfinite(r_part[j]):
            continue
        p = _p_from_r(r_part[j], int(n_xy[j]), n_controls=1)
        if p is None:
            continue
        cands.append(DrugCorrelation(
            drug=(names or {}).get(col, col), column=col,
            r=round(float(r_part[j]), 4),
            r_raw=round(float(r_xy[j]), 4) if np.isfinite(r_xy[j]) else None,
            n_lines=int(n_xy[j]), p_value=p, readout=readout))
        pvals.append(p)

    res.n_drugs_tested = len(cands)
    if not cands:
        res.error = (
            f"no compound was screened in at least {min_lines} of the lines "
            f"that also carry a {readout} value for {res.gene}.")
        return res

    for c, q in zip(cands, _bh(pvals)):
        c.q_value = q
        # Sign conventions differ between the two readouts, and this is
        # the one place that has to be right. Chronos and PRISM both run
        # negative for harm, so a drug that hurts the same lines as the
        # knockout correlates POSITIVELY. Expression runs the other way:
        # more transcript is a higher number while more killing is a lower
        # one, so sensitivity rising with expression is NEGATIVE.
        if readout == "phenocopy":
            c.direction = ("phenocopies" if (c.r or 0) > 0 else "opposes")
        else:
            c.direction = ("tracks expression" if (c.r or 0) < 0
                           else "inverse to expression")

    res.assessed = True
    hits = [c for c in cands
            if abs(c.r or 0.0) >= min_abs_r
            and (c.q_value if c.q_value is not None else 1.0) <= max_q]
    hits.sort(key=lambda c: -abs(c.r or 0.0))
    res.hits = hits[:top_n]

    res.caveats.append(
        f"{len(cands)} compound(s) tested against {res.gene}; "
        f"{len(hits)} cleared |r| >= {min_abs_r} at q <= {max_q}. The "
        f"correlation is partial, controlling for how easily each line "
        f"dies under the panel as a whole - without that control most of "
        f"these numbers measure cell-line fragility.")
    if hits:
        worst_gap = max(abs((h.r_raw or 0) - (h.r or 0)) for h in hits)
        res.caveats.append(
            f"the control moved a reported correlation by up to "
            f"{worst_gap:.2f}. It can move either way: a shared fragility "
            f"component inflates a correlation when the drug and the "
            f"knockout both track it, and masks one when they do not, so "
            f"a large gap means the raw number was mostly about the cell "
            f"lines rather than about this factor.")
    res.caveats.append(
        "a correlated drug may act on this factor, on anything else the "
        "same lines need, or on a pathway the factor sits in. This ranks "
        "candidates for a binding assay; it does not identify a target.")
    res.caveats.append(
        "both panels score proliferation of transformed cells in a dish, "
        "so a factor whose role is in a differentiated cell looks silent "
        "here for reasons unrelated to its importance")
    return res


def phenocopy(gene: str, *, prism=None, effect_path=None, names=None,
              **kw) -> SensitivityResult:
    """Drugs whose viability profile matches this gene's knockout profile."""
    from . import depmap as dm

    prism = prism if prism is not None else load_prism()
    if not prism:
        return correlate({}, {}, gene=gene, readout="phenocopy")

    p = dm.resolve(effect_path)
    if not p:
        res = SensitivityResult(gene=gene.upper(), readout="phenocopy")
        res.error = ("no DepMap gene-effect file on disk, so there is no "
                     "knockout profile to match drugs against. " + dm.SOURCE)
        return res

    data = dm.load_gene_effect(p, [gene])
    models = data.pop("__models__", [])
    effects = data.get(gene.strip().upper())
    if not effects:
        res = SensitivityResult(gene=gene.upper(), readout="phenocopy")
        res.error = (f"{gene.upper()} is not scored in the DepMap release "
                     f"on disk, so no knockout profile exists to match. "
                     f"Absent from the file is not the same as a gene "
                     f"whose knockout does nothing.")
        return res

    target = {m: v for m, v in zip(models, effects) if v is not None}
    return correlate(prism, target, readout="phenocopy", gene=gene,
                     names=names if names is not None
                     else load_compound_names(), **kw)


def expression_dependence(gene: str, *, prism=None, expression_path=None,
                          names=None, **kw) -> SensitivityResult:
    """Drugs whose effect tracks how much of this gene a line expresses."""
    prism = prism if prism is not None else load_prism()
    if not prism:
        return correlate({}, {}, gene=gene, readout="expression")
    target = load_expression(expression_path, gene)
    return correlate(prism, target, readout="expression", gene=gene,
                     names=names if names is not None
                     else load_compound_names(), **kw)


def summarise(res: SensitivityResult) -> str:
    """Prose for the report, including the case where nothing was found."""
    if res is None:
        return ""
    if not res.assessed:
        return res.error or "drug sensitivity was not assessed."

    if not res.hits:
        return (
            f"{res.n_drugs_tested} compound(s) were correlated against "
            f"{res.gene}'s {res.readout} profile across {res.n_lines} cell "
            f"line(s) and none cleared the threshold. That is the common "
            f"outcome and it is informative: after controlling for how "
            f"easily each line dies, no approved compound in the panel "
            f"behaves as though it acts on this factor. It bounds what "
            f"this method can contribute rather than saying no such drug "
            f"exists - the panel is proliferating cancer lines, and a "
            f"factor that matters elsewhere leaves no trace in it.")

    lead = res.hits[0]
    if res.readout == "phenocopy":
        what = (f"**{lead.drug}** hurts the same cell lines that losing "
                f"{res.gene} hurts (partial r = {lead.r}, raw "
                f"{lead.r_raw}, {lead.n_lines} lines, q = {lead.q_value})")
    else:
        what = (f"**{lead.drug}** is most effective in the lines "
                f"expressing {res.gene} most (partial r = {lead.r}, raw "
                f"{lead.r_raw}, {lead.n_lines} lines, q = {lead.q_value})")

    return (
        f"{len(res.hits)} of {res.n_drugs_tested} compound(s) correlate "
        f"with {res.gene}'s {res.readout} profile. {what}. Nothing here "
        f"records any of them binding {res.gene}: the argument is that "
        f"they behave as though they do, which is a reason to run a "
        f"binding assay and not a substitute for one.")
