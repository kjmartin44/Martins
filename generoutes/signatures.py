"""LINCS L1000 transcriptomic signatures — the only empirical link in the package.

Everything else here is inference: a motif resembles a site, a factor is
annotated as an activator, a drug is annotated as an agonist. L1000 measured
what actually happened to the transcriptome when cells were dosed. That makes
it the most valuable evidence available and, for the same reason, the evidence
most often over-read.

Two distinct uses
-----------------
**Direct query.** Which approved drugs moved *your gene* in the direction you
want? Skips every inferential hop.

**Connectivity query (the one you asked about).** Find drugs whose whole
signature resembles that of a known activator of your gene, *regardless of
annotated target*. This is what CMap was designed for, and it is how you find a
compound whose official mechanism has nothing to do with your pathway but whose
cellular effect does. It is also where most false positives come from.

Four things that decide whether a hit is real
---------------------------------------------
**1. Landmark vs inferred.** L1000 physically measures ~978 landmark
transcripts and *computationally infers* the other ~11,350. If your gene is not
a landmark, its value is a model prediction from the landmarks, not a
measurement. Inferred values are much weaker evidence, and the distinction is
almost never stated in papers that use this data. Hits are weighted
accordingly, and the report always says which kind you have.

**2. Cell line.** Signatures exist in a fixed panel (MCF7, PC3, A375, HA1E,
VCAP, HEPG2, HT29, A549 and others). If your gene is not expressed in the line
that produced the signature, there is nothing to move and absence of signal
means nothing. Always read the cell line before the score.

**3. Non-specific stress dominates.** A large fraction of strong L1000
signatures are some version of "this cell is dying": HSP induction, ATF4/CHOP,
cell-cycle collapse. Cytotoxic and promiscuous compounds therefore connect to
almost everything, and a naive connectivity search returns them for any query.
Compounds whose signatures are strong across many unrelated queries are flagged
as non-specific rather than reported as hits.

**4. Dose and time.** 6 h and 24 h signatures of the same compound can differ
substantially, and high-dose signatures drift toward the stress response above.

Sources
-------
  SigCom LINCS  https://maayanlab.cloud/sigcom-lincs  (free API)
  CLUE / CMap   https://clue.io  (registration; the canonical resource)
  iLINCS        http://www.ilincs.org

You can also export a signature from CLUE by hand and pass it with
--signature-file; the scoring logic here works on a local file with no network.
"""

from __future__ import annotations

import logging
import math
import re
from dataclasses import dataclass, asdict, field
from pathlib import Path

from .http import Client, HttpError

log = logging.getLogger(__name__)

SIGCOM = "https://maayanlab.cloud/sigcom-lincs/metadata-api"
SIGCOM_DATA = "https://maayanlab.cloud/sigcom-lincs/data-api/api/v1"

# Evidence weighting for the landmark/inferred distinction.
LANDMARK_WEIGHT = 1.0
INFERRED_WEIGHT = 0.40

# A compound connecting to this fraction of unrelated queries is behaving like
# a general stress signature rather than a specific mechanism.
PROMISCUITY_FLAG = 0.25

# A small, well-known subset of the 978 landmarks, for offline checks. The full
# list comes from the API; this is a fallback so the landmark distinction is
# never silently skipped.
KNOWN_LANDMARKS = {
    "HSPA8", "HMOX1", "DDIT3", "ATF3", "EGR1", "JUN", "FOS", "MYC", "CDKN1A",
    "GADD45A", "TP53", "RELB", "NFKBIA", "SOD1", "CAT", "TXNRD1", "NQO1",
    "HSPA1A", "DNAJB1", "PPARG", "NR3C1", "VDR", "ESR1", "AR", "STAT3",
    "CEBPB", "SREBF1", "HIF1A", "FOXO3", "E2F1", "CCND1", "MKI67",
}


@dataclass
class SignatureHit:
    drug: str
    gene: str
    direction: str              # "up" | "down"
    score: float                # signed effect (z-like) or connectivity tau
    cell_line: str | None = None
    dose: str | None = None
    timepoint: str | None = None
    is_landmark: bool | None = None
    evidence_weight: float = INFERRED_WEIGHT
    weighted_score: float = 0.0
    query_type: str = "direct"  # "direct" | "connectivity"
    reference_drug: str | None = None
    non_specific: bool = False
    feature_space: str | None = None   # landmark | best inferred | inferred
    notes: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


# ------------------------------------------------------------ landmark check
def is_landmark(gene: str, client: Client | None = None) -> bool | None:
    """True/False if known, None if it cannot be determined.

    None matters: it means "we do not know whether this was measured or
    inferred", which is different from knowing it was inferred.
    """
    g = gene.upper()
    if g in KNOWN_LANDMARKS:
        return True
    if client is None:
        return None
    try:
        res = client.get(f"{SIGCOM}/entities/find", {
            "filter": f'{{"where":{{"meta.symbol":"{g}"}}}}'})
    except HttpError as exc:
        log.warning("SigCom lookup for %s failed (%s); it is treated as "
                    "not found, which is not a finding", g, str(exc)[:80])
        return None
    rows = res if isinstance(res, list) else res.get("entities", []) or []
    for r in rows:
        meta = (r or {}).get("meta") or {}
        for key in ("pr_is_lm", "is_landmark", "landmark"):
            if key in meta:
                return bool(meta[key])
    return None


def weight_for(landmark: bool | None) -> tuple[float, str | None]:
    if landmark is True:
        return LANDMARK_WEIGHT, None
    if landmark is False:
        return INFERRED_WEIGHT, (
            "gene is NOT an L1000 landmark: this value was computationally "
            "inferred from the 978 measured transcripts, not measured. Treat it "
            "as substantially weaker evidence and confirm by RT-qPCR")
    return INFERRED_WEIGHT, (
        "could not determine whether this gene is an L1000 landmark; scored as "
        "inferred to stay conservative")


# --------------------------------------------------------- connectivity math
def connectivity(sig_a: dict[str, float], sig_b: dict[str, float]) -> float:
    """Signed cosine over shared genes. -1 opposed, +1 concordant.

    CMap's published tau uses a weighted Kolmogorov-Smirnov statistic against a
    reference population; a signed cosine over the overlap is a transparent
    stand-in that needs no background distribution. It ranks similarly for
    strong signatures and does not pretend to be a calibrated tau.
    """
    shared = set(sig_a) & set(sig_b)
    if len(shared) < 10:
        return 0.0
    a = [sig_a[g] for g in shared]
    b = [sig_b[g] for g in shared]
    na = math.sqrt(sum(x * x for x in a))
    nb = math.sqrt(sum(x * x for x in b))
    if na == 0 or nb == 0:
        return 0.0
    return sum(x * y for x, y in zip(a, b)) / (na * nb)


def load_signature_file(path: str | Path) -> dict[str, float]:
    """Read a two-column gene,score CSV/TSV exported from CLUE or iLINCS."""
    path = Path(path)
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    sig: dict[str, float] = {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            parts = [p.strip().strip('"') for p in line.rstrip("\n").split(sep)]
            if len(parts) < 2:
                continue
            try:
                sig[parts[0].upper()] = float(parts[1])
            except ValueError:
                if i == 0:
                    continue        # header
    log.info("loaded signature with %d genes from %s", len(sig), path.name)
    return sig


def flag_non_specific(hits: list[SignatureHit],
                      threshold: float = PROMISCUITY_FLAG) -> list[SignatureHit]:
    """Mark compounds that connect to an implausible share of queries.

    Cytotoxic and promiscuous compounds produce a generic stress signature that
    resembles everything. Without this filter, a connectivity search returns
    the same two dozen compounds regardless of what you asked.
    """
    if not hits:
        return hits
    counts: dict[str, int] = {}
    for h in hits:
        counts[h.drug.lower()] = counts.get(h.drug.lower(), 0) + 1
    n_queries = len({h.reference_drug or h.gene for h in hits}) or 1

    for h in hits:
        n_seen = counts[h.drug.lower()]
        share = n_seen / n_queries
        # A compound appearing in a single query cannot be promiscuous, however
        # small the query set. Requiring >1 avoids flagging specific hits when
        # few queries were run.
        if n_queries >= 4 and n_seen > 1 and share >= threshold:
            h.non_specific = True
            h.notes.append(
                f"connects to {share:.0%} of queries — this is the signature of a "
                "cytotoxic or promiscuous compound (generic stress response), "
                "not evidence of a specific mechanism")
            h.weighted_score *= 0.25
    return hits


# ----------------------------------------------------------------- querying
def query_gene(client: Client, gene: str, *, direction: str = "up",
               limit: int = 50) -> list[SignatureHit]:
    """Approved drugs whose L1000 signature moves this gene the wanted way."""
    lm = is_landmark(gene, client)
    weight, note = weight_for(lm)

    try:
        payload = {("up_entities" if direction == "up" else "down_entities"): [gene.upper()],
                   "limit": limit, "database": "l1000_cp"}
        res = client.post(f"{SIGCOM_DATA}/enrich/ranktwosided", payload)
    except HttpError as exc:
        log.warning("LINCS query failed for %s: %s", gene, exc)
        return []

    rows = (res.get("results") if isinstance(res, dict) else res) or []
    hits: list[SignatureHit] = []
    for r in rows:
        if not isinstance(r, dict):
            continue
        meta = r.get("meta") or r
        name = (meta.get("pert_name") or meta.get("drug") or meta.get("id") or "?")
        # `a or b` would read a z-score of exactly 0.0 - a gene that did
        # not move - as missing and silently substitute a different
        # column. These arrive parsed from JSON, so 0.0 is a real float
        # here, unlike the CSV path below where every value is a string
        # and "0" is truthy.
        _z = r.get("zscore")
        score = float(_z if _z is not None else (r.get("score") or 0.0))
        h = SignatureHit(
            drug=str(name).title(), gene=gene.upper(), direction=direction,
            score=round(score, 3), cell_line=meta.get("cell_line"),
            dose=meta.get("pert_dose"), timepoint=meta.get("pert_time"),
            is_landmark=lm, evidence_weight=weight,
            weighted_score=round(abs(score) * weight, 3), query_type="direct")
        if note:
            h.notes.append(note)
        if h.cell_line:
            h.notes.append(
                f"measured in {h.cell_line}: confirm your gene is expressed in "
                "that line, or absence of signal means nothing")
        hits.append(h)

    hits.sort(key=lambda h: -h.weighted_score)
    return hits


def connected_to_reference(reference: dict[str, float],
                           candidates: dict[str, dict[str, float]], *,
                           min_score: float = 0.3,
                           top_n: int = 25) -> list[SignatureHit]:
    """Rank candidate signatures by connectivity to a reference signature.

    `reference` is the signature of a known activator of your gene. Drugs whose
    whole-transcriptome response resembles it are candidates *whatever their
    annotated target*, which is the entire point of the approach.
    """
    hits: list[SignatureHit] = []
    for drug, sig in candidates.items():
        c = connectivity(reference, sig)
        if abs(c) < min_score:
            continue
        h = SignatureHit(
            drug=drug, gene="(signature)", direction="up" if c > 0 else "down",
            score=round(c, 3), query_type="connectivity",
            evidence_weight=LANDMARK_WEIGHT,
            weighted_score=round(abs(c), 3))
        h.notes.append(
            "matched on whole-signature similarity, not on target annotation. "
            "Similar signatures can arise from shared cytotoxicity rather than "
            "shared mechanism — check the compound is not simply killing cells")
        hits.append(h)
    hits.sort(key=lambda h: -h.weighted_score)
    return hits[:top_n]


def summarise(hits: list[SignatureHit], gene: str) -> str:
    if not hits:
        return (f"No LINCS signature found moving {gene}. Check whether {gene} is "
                "expressed in the L1000 cell panel at all — if it is not, this "
                "says nothing about your drug.")
    lm = hits[0].is_landmark
    kind = ("measured directly (landmark transcript)" if lm is True
            else "computationally inferred, not measured" if lm is False
            else "landmark status unknown")
    specific = [h for h in hits if not h.non_specific]
    return (f"{len(hits)} signature hits for {gene}; values are {kind}. "
            f"{len(hits) - len(specific)} flagged as non-specific stress "
            f"signatures. This is the only empirical evidence in the report — "
            f"where it agrees with the predicted routes, you have a real "
            f"candidate.")


# ------------------------------------------------- offline signature table
def load_signature_table(path, gene: str, *, direction: str = "up",
                         min_abs_score: float = 1.5) -> list[SignatureHit]:
    """Read drug/gene expression effects from a local table.

    For when the LINCS API is unreachable, or you want a pinned, reproducible
    input. Expected columns (header required, order free, comma or tab):

        drug, gene, score  [, cell_line] [, dose] [, time] [, is_landmark]

    `score` is signed — positive means the drug raised the gene. Export one
    from CLUE (https://clue.io), iLINCS (http://www.ilincs.org) or SigCom
    LINCS, or assemble it from your own RNA-seq.

    The same landmark/inferred weighting applies: if the table carries an
    `is_landmark` column it is honoured, otherwise the value is scored as
    inferred, which is the conservative reading.
    """
    import csv
    import gzip
    import io
    from pathlib import Path

    path = Path(path)
    fh = (io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8",
                           errors="replace") if path.suffix == ".gz"
          else open(path, encoding="utf-8", errors="replace"))
    target = gene.upper()
    hits: list[SignatureHit] = []

    with fh:
        sample = fh.read(8192)
        fh.seek(0)
        delim = "\t" if sample.count("\t") >= sample.count(",") else ","
        for row in csv.DictReader(fh, delimiter=delim):
            row = {(k or "").strip().lower(): v for k, v in row.items()}
            if (row.get("gene") or "").strip().upper() != target:
                continue
            try:
                score = float(row.get("score") or row.get("zscore") or 0)
            except ValueError:
                continue
            if abs(score) < min_abs_score:
                continue
            if (direction == "up") != (score > 0):
                continue

            lm_raw = (row.get("is_landmark") or "").strip().lower()
            lm = (True if lm_raw in ("1", "true", "yes", "y")
                  else False if lm_raw in ("0", "false", "no", "n") else None)
            weight, note = weight_for(lm)

            h = SignatureHit(
                drug=(row.get("drug") or row.get("pert_name") or "?").title(),
                gene=target, direction=direction, score=round(score, 3),
                cell_line=row.get("cell_line") or row.get("cell"),
                dose=row.get("dose"), timepoint=row.get("time"),
                is_landmark=lm, evidence_weight=weight,
                weighted_score=round(abs(score) * weight, 3),
                query_type="direct")
            h.notes.append("from a local signature table, not the LINCS API")
            if note:
                h.notes.append(note)
            hits.append(h)

    hits.sort(key=lambda h: -h.weighted_score)
    log.info("loaded %d local signature hits for %s", len(hits), target)
    return hits

# ------------------------------------------------- offline LINCS (Enrichr)
# SigCom is the richest LINCS interface and also the least reliable: it was
# returning HTTP 500 for an entire working session, which costs ~34s of
# patient retrying per endpoint and leaves the report with no empirical
# evidence at all.
#
# Enrichr publishes the same L1000 chemical-perturbation data as two plain
# .gmt files, 16 MB and 14 MB, and they download in seconds:
#
#   https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text
#       &libraryName=LINCS_L1000_Chem_Pert_up
#       &libraryName=LINCS_L1000_Chem_Pert_down
#
# Fetch them once into data/ and the direct query runs offline forever.
#
# What is gained and lost against the API. Gained: it works, it is pinned and
# reproducible, and the signature name carries the cell line, timepoint and
# dose that the report already treats as decisive. Lost: the API returns a
# continuous score per gene, whereas a .gmt is set membership - a gene is in
# the signature or it is not. So a hit here means "this compound moved your
# gene in a consensus signature", with no magnitude. Treated as a fixed
# effect size and labelled as set membership, never as a measured z-score.
GMT_SET_SCORE = 2.0          # nominal; membership carries no magnitude

DEFAULT_GMT_PATHS = {
    "up": ("data/lincs_chem_pert_up.gmt", "lincs_chem_pert_up.gmt"),
    "down": ("data/lincs_chem_pert_down.gmt", "lincs_chem_pert_down.gmt"),
}


def _parse_sig_name(name: str) -> dict:
    """Pull cell line, timepoint, drug and dose out of a signature label.

    Enrichr labels look like `LJP006 HA1E 24H-AT-7519-10`: experiment and
    cell line and duration, then the drug and dose joined by hyphens. The
    drug itself frequently contains hyphens (AT-7519, DCC-2036), so the dose
    is taken as the trailing numeric field and the drug is whatever remains.
    """
    head, _, tail = name.partition("-")
    parts = head.split()
    cell = parts[1] if len(parts) > 1 else None
    time = parts[2] if len(parts) > 2 else None

    drug, dose = tail, None
    if tail:
        bits = tail.rsplit("-", 1)
        if len(bits) == 2:
            try:
                float(bits[1])
                drug, dose = bits[0], bits[1]
            except ValueError:
                pass
    return {"drug": (drug or name).strip(), "cell_line": cell,
            "timepoint": time, "dose": dose}


# --------------------------------------------------- how broadly a compound acts
# A .gmt signature lists the genes that passed a cutoff in one compound,
# cell line, dose and time. The lists differ a lot in length (median 119 genes
# up and down together, longest 648), so the length says how many genes that
# condition moved. A compound whose signatures are long moves a large part of
# the transcriptome whichever gene you asked about, and a hit on your gene from
# it is weak evidence of specificity.
#
# This is a caution and not a verdict. The files hold set membership with no
# scores, so they cannot say whether YOUR gene stands out among the genes a
# compound moved. The length rises with dose, depends on the cell line, and
# the cutoff is whatever Enrichr used. Nothing here changes a score.
BROAD_PERCENTILE = 75            # top quartile of all signatures
# The tier of breadth that costs a route score, not just a flag: the top
# decile. Chosen to catch the compounds that move a large part of the
# transcriptome whatever gene is asked about (vorinostat, trichostatin A,
# emetine) and to leave the merely broad ones (the CDK inhibitors, about
# the 77th percentile) flagged but not scored down.
VERY_BROAD_PERCENTILE = 90
MIN_SIGNATURES_FOR_BREADTH = 3   # fewer than this is too few to judge

_BREADTH_CACHE: dict = {}


def breadth_key(name: str) -> str:
    """A compound's name reduced so spellings of it meet: case, hyphens and
    spaces, and trailing salt words ("Vorinostat", "trichostatin-a",
    "Emetine Dihydrochloride Hydrate")."""
    from .congeners import salt_free_name
    spaced = re.sub(r"[-_]+", " ", (name or "").lower())
    return re.sub(r"[^a-z0-9]", "", salt_free_name(spaced))


def _count_genes(path) -> dict:
    """`{signature name: number of genes}` for one .gmt file."""
    counts: dict = {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2 or not f[0]:
                continue
            counts[f[0]] = len([g for g in f[2:] if g.strip()])
    return counts


def signature_breadth(up_path=None, down_path=None) -> dict | None:
    """How many genes each compound moves, from the local .gmt files.

    Returns `{"n_signatures", "median_genes_all", "compounds": {key: {
    "median_genes", "n_signatures", "percentile"}}}`, or None when either file
    is missing. `percentile` is the share of ALL signatures that moved fewer
    genes than this compound's median. Read once per process and remembered.
    """
    import bisect
    import statistics
    from pathlib import Path as _P

    def find(direction, explicit):
        cands = [explicit] if explicit else list(DEFAULT_GMT_PATHS[direction])
        return next((_P(c) for c in cands if c and _P(c).exists()), None)

    up, down = find("up", up_path), find("down", down_path)
    if up is None or down is None:
        return None
    key = (str(up.resolve()), up.stat().st_size, up.stat().st_mtime_ns,
           str(down.resolve()), down.stat().st_size, down.stat().st_mtime_ns)
    if key in _BREADTH_CACHE:
        return _BREADTH_CACHE[key]

    n_up, n_down = _count_genes(up), _count_genes(down)
    per_signature = {n: n_up.get(n, 0) + n_down.get(n, 0)
                     for n in set(n_up) | set(n_down)}
    everything = sorted(per_signature.values())
    by_compound: dict = {}
    for name, n in per_signature.items():
        by_compound.setdefault(breadth_key(_parse_sig_name(name)["drug"]),
                               []).append(n)
    compounds = {}
    for k, counts in by_compound.items():
        if not k:
            continue
        med = statistics.median(counts)
        compounds[k] = {
            "median_genes": med, "n_signatures": len(counts),
            "percentile": round(100.0 * bisect.bisect_left(everything, med)
                                / len(everything), 1)}
    out = {"n_signatures": len(everything),
           "median_genes_all": statistics.median(everything),
           "compounds": compounds}
    _BREADTH_CACHE[key] = out
    return out


def breadth_for(names, breadth: dict | None) -> dict:
    """Breadth records for the compounds in `names`, keyed by `breadth_key`.

    Each record adds `broad` (top quartile, with enough signatures to judge)
    and `thin` (too few signatures to say either way). Compounds the LINCS
    files do not contain are left out, which is not the same as "narrow".
    """
    if not breadth:
        return {}
    out = {}
    for name in names:
        k = breadth_key(name)
        rec = breadth["compounds"].get(k)
        if rec is None or k in out:
            continue
        thin = rec["n_signatures"] < MIN_SIGNATURES_FOR_BREADTH
        out[k] = dict(
            rec, thin=thin,
            broad=(not thin and rec["percentile"] >= BROAD_PERCENTILE),
            very_broad=(not thin
                        and rec["percentile"] >= VERY_BROAD_PERCENTILE))
    return out


def load_enrichr_gmt(gene: str, *, direction: str = "up",
                     up_path=None, down_path=None) -> list:
    """Compounds whose L1000 consensus signature contains `gene`.

    Reads the direction you asked for: `up` returns compounds that RAISED the
    gene, `down` those that lowered it.
    """
    from pathlib import Path as _P

    want = direction if direction in ("up", "down") else "up"
    explicit = {"up": up_path, "down": down_path}[want]
    candidates = ([explicit] if explicit else list(DEFAULT_GMT_PATHS[want]))

    path = next((_P(c) for c in candidates if c and _P(c).exists()), None)
    if path is None:
        raise FileNotFoundError(
            f"no local LINCS .gmt for direction {want!r}. Download it once:\n"
            f"  curl -L -o data/lincs_chem_pert_{want}.gmt "
            f'"https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text'
            f'&libraryName=LINCS_L1000_Chem_Pert_{want}"')

    target = gene.strip().upper()
    hits: list = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            if target not in {g.strip().upper() for g in fields[2:] if g.strip()}:
                continue
            meta = _parse_sig_name(fields[0])
            h = SignatureHit(
                drug=meta["drug"], gene=target, direction=want,
                score=GMT_SET_SCORE if want == "up" else -GMT_SET_SCORE,
                cell_line=meta["cell_line"], dose=meta["dose"],
                timepoint=meta["timepoint"],
                # Real confidence, not a guess: CLUE says whether this
                # gene was measured on the array or imputed from the 978
                # landmarks. Defaulting to None demoted every hit on a
                # non-landmark gene below curated annotation.
                is_landmark=(feature_space(target) == "landmark"
                             if feature_space(target) else None),
                query_type="direct")
            h.feature_space = feature_space(target)
            h.notes.append(confidence_note(target))
            h.notes.append(
                "set membership in a LINCS consensus signature (Enrichr "
                ".gmt): the compound moved this gene, with no magnitude "
                "attached - confirm by RT-qPCR before relying on it")
            hits.append(h)

    for h in hits:
        h.evidence_weight = INFERRED_WEIGHT
        h.weighted_score = round(h.score * h.evidence_weight, 3)

    log.info("LINCS .gmt (%s): %d compound signature(s) contain %s",
             path.name, len(hits), target)
    return hits


def autoload_gmt_available(direction: str = "up") -> bool:
    from pathlib import Path as _P
    return any(_P(c).exists() for c in DEFAULT_GMT_PATHS.get(direction, ()))

# ------------------------------------------- L1000 measurement confidence
# The package carried a binary is_landmark, and `catalogue.from_lincs` read
# it as "tier 2 if measured, else tier 4" - so an empirical hit on a
# non-landmark gene ranked BELOW a curated ligand annotation. That is too
# harsh for most genes, because CLUE distinguishes three levels and the
# binary collapses the middle one:
#
#   landmark        978  directly measured on the array
#   best inferred  9196  imputed from the landmarks, high confidence
#   inferred       2154  imputed, lower confidence
#
# CCR5, VDR and HBB are all "best inferred". Treating that as equivalent to
# "a model says so" understates it: the signature itself was measured, and
# the value for your gene was imputed by a model trained on 978 genes that
# were. Tiering it alongside curated annotation is the honest middle.
#
#   https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/geneinfo_beta.txt
_FEATURE_SPACE: dict = {}
DEFAULT_GENEINFO_PATHS = ("data/lincs_geneinfo.txt", "lincs_geneinfo.txt")


def load_geneinfo(path=None) -> int:
    """Index each gene's L1000 measurement confidence."""
    from pathlib import Path as _P
    TAB = chr(9)
    if _FEATURE_SPACE:
        return len(_FEATURE_SPACE)
    for cand in ([path] if path else list(DEFAULT_GENEINFO_PATHS)):
        if not cand or not _P(cand).exists():
            continue
        text = _P(cand).read_text(encoding="utf-8", errors="replace")
        lines = text.splitlines()
        if not lines:
            continue
        hdr = lines[0].split(TAB)
        try:
            i_sym = hdr.index("gene_symbol")
            i_fs = hdr.index("feature_space")
        except ValueError:
            continue
        for ln in lines[1:]:
            f = ln.split(TAB)
            if len(f) > max(i_sym, i_fs) and f[i_sym]:
                _FEATURE_SPACE[f[i_sym].strip().upper()] = f[i_fs].strip()
        log.info("L1000 measurement confidence indexed for %d genes from %s",
                 len(_FEATURE_SPACE), _P(cand).name)
        return len(_FEATURE_SPACE)
    return 0


def feature_space(gene: str) -> str | None:
    """"landmark" | "best inferred" | "inferred", or None if unknown."""
    load_geneinfo()
    return _FEATURE_SPACE.get((gene or "").strip().upper())


def confidence_note(gene: str) -> str:
    fs = feature_space(gene)
    if fs == "landmark":
        return (f"{gene.upper()} is one of the 978 L1000 LANDMARK genes, so "
                f"its value was physically measured on the array. This is "
                f"the strongest form this evidence takes.")
    if fs == "best inferred":
        return (f"{gene.upper()} is 'best inferred': the signature was "
                f"measured, but the value for this gene was imputed from the "
                f"978 landmarks rather than measured directly. High "
                f"confidence among inferred genes, and still not a "
                f"measurement of your gene - confirm by RT-qPCR.")
    if fs == "inferred":
        return (f"{gene.upper()} is 'inferred' at the lower confidence level. "
                f"Its value is a model output two steps from any "
                f"measurement; treat it as a lead, not evidence.")
    return (f"{gene.upper()} is not in the L1000 gene space, or the "
            f"confidence table is not loaded, so whether any value for it "
            f"was measured or imputed is unknown - treated as inferred.")
