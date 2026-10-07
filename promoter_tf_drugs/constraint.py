"""How much dosage this gene tolerates, in either direction.

Every route in this package ends with "and then your gene moves". That
assumes two things nothing else here checks: that there is room to move
it, and that moving it is safe. For a rare loss-of-function condition
those are the whole question - you are not looking for any change, you
are looking for enough change, without causing a second problem.

Population constraint answers part of it. gnomAD counts how many
loss-of-function variants a gene carries across ~800,000 people and
compares that with how many would be expected by chance. A gene with far
fewer than expected has been under selection: losing a copy mattered to
the people who carried it.

**What this does and does not tell you**

It tells you whether dosage matters. A gene where losing one copy is
selected against is a gene whose output level is being watched by
something, which makes a partial loss plausible as a disease mechanism
and makes overshooting a real risk.

It does not tell you that a drug can raise the gene, how far, or that
raising it helps. Nothing in a population database speaks to that.

**The recessive trap, which matters here**

A gene can be entirely LoF-tolerant in the heterozygous state - high
LOEUF, pLI near zero - and still cause severe disease when BOTH copies
are hit. Mosaic variegated aneuploidy is biallelic BUB1B; carriers are
healthy, so selection never saw the phenotype. Reading "tolerant" as
"this gene does not matter" would be exactly backwards for the
conditions this module exists to serve, and the summary says so rather
than leaving it to be inferred.

**Sources, both offline**

  gnomAD constraint metrics (one row per transcript):
      https://gnomad.broadinstitute.org/downloads
      -> data/gnomad_constraint.tsv

  ClinGen dosage sensitivity curation (haplo- and triplosensitivity,
  curated by a panel rather than computed):
      https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv
      -> data/clingen_dosage.tsv

Neither is required. Without them this reports that dosage was not
assessed, which is not the same as a gene being tolerant.
"""

from __future__ import annotations

import csv
import logging
import re
from dataclasses import dataclass, field, asdict
from pathlib import Path

log = logging.getLogger(__name__)

DEFAULT_GNOMAD_PATHS = (
    "data/gnomad_constraint.tsv",
    "data/gnomad.v4.1.constraint_metrics.tsv",
    "data/gnomad.v2.1.1.lof_metrics.by_gene.txt",
)
DEFAULT_CLINGEN_PATHS = (
    "data/clingen_dosage.tsv",
    "data/ClinGen_gene_curation_list_GRCh38.tsv",
)

# pLI is the probability the gene is intolerant of losing one copy. The
# 0.9 line is gnomAD's own convention and is what the literature cites.
PLI_SENSITIVE = 0.9

# LOEUF is the upper bound of the observed/expected ratio for LoF
# variants: lower means more constrained. gnomAD recommends ranking on
# it rather than on pLI, and treats below 0.35 as the constrained end.
LOEUF_SENSITIVE = 0.35
# Above this the gene demonstrably tolerates losing a copy.
LOEUF_TOLERANT = 1.0

_CACHE: dict = {}


@dataclass(frozen=True)
class Source:
    """Which file answered, and which release of it.

    gnomAD v2.1.1 and v4.1 disagree about the same gene - different
    cohort, different transcript set, LOEUF recomputed - so a pLI
    without a release is a number without a provenance. ClinGen is
    re-dated weekly and genes move from "no evidence available" to
    curated between releases.

    `release` is None when the file does not say, which is honest:
    guessing a release is worse than recording none, because the point
    of the field is to make two runs comparable.
    """
    name: str
    path: str
    release: str | None = None
    release_from: str | None = None     # where the release was read


@dataclass
class Constraint:
    """What is known about this gene's tolerance of a dosage change."""
    symbol: str
    assessed: bool = False
    pli: float | None = None
    loeuf: float | None = None
    haploinsufficient: bool | None = None     # ClinGen curation
    triplosensitive: bool | None = None       # ClinGen curation
    dosage: str = "unknown"                   # sensitive|tolerant|unknown
    notes: list = field(default_factory=list)
    # The files this came from. Empty when nothing was on disk, which
    # is what `assessed is False` already says.
    sources: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


def resolve(explicit=None, candidates=DEFAULT_GNOMAD_PATHS) -> str | None:
    """First of `candidates` that exists, or `explicit` when given."""
    if explicit:
        return str(explicit) if Path(explicit).is_file() else None
    for c in candidates:
        if Path(c).is_file():
            return c
    return None


def _column(header: list, *names) -> int | None:
    low = [h.strip().lower() for h in header]
    for n in names:
        if n in low:
            return low.index(n)
    return None


# Where each loaded file came from, keyed as _CACHE is. Populated by
# the loaders so `assess` can report it without re-reading anything.
_SOURCES: dict = {}


def _gnomad_release(path: str, header: list) -> tuple:
    """(release, where it was read) for a gnomAD constraint table."""
    m = re.search(r"v?(\d+\.\d+(?:\.\d+)?)", Path(path).name)
    if m:
        return f"v{m.group(1)}", "filename"
    # No version in the name. The column spelling still separates the
    # generations: v4 nests its fields, v2 flattens them.
    low = {h.strip().lower() for h in header}
    if "lof.oe_ci.upper" in low or "lof.pli" in low:
        return "v4 series", "column names"
    if "oe_lof_upper" in low:
        return "v2 series", "column names"
    return None, None


def _clingen_release(lines: list) -> tuple:
    """(release, where it was read) for a ClinGen curation list.

    The file dates itself on a comment line near the top, in the form
    "#07 Oct,2026". Only comment lines are searched, so a gene whose
    name happens to look like a date cannot be mistaken for one.
    """
    for line in lines[:10]:
        if not line.startswith("#"):
            break
        m = re.match(r"#\s*(\d{1,2}\s+[A-Za-z]{3,}\s*,?\s*\d{4})\s*$",
                     line)
        if m:
            return re.sub(r"\s+", " ", m.group(1)).strip(), "file header"
    return None, None


def _flag(row: list, i: int | None) -> bool:
    """A gnomAD boolean column, which is spelled several ways."""
    if i is None or len(row) <= i:
        return False
    return str(row[i]).strip().lower() in ("true", "t", "1", "yes")


def load_gnomad(path=None) -> dict:
    """{SYMBOL: (pli, loeuf)} from a gnomAD constraint table."""
    p = resolve(path)
    if not p:
        return {}
    key = ("gnomad", p)
    if key in _CACHE:
        return _CACHE[key]

    out: dict = {}
    try:
        with open(p, "r", encoding="utf-8", errors="replace") as fh:
            rdr = csv.reader(fh, delimiter="\t")
            header = next(rdr, None)
            if not header:
                return {}
            i_sym = _column(header, "gene", "gene_symbol", "symbol")
            i_pli = _column(header, "pli", "lof.pli", "lof_pli")
            i_loeuf = _column(header, "loeuf", "oe_lof_upper",
                              "lof.oe_ci.upper")
            i_mane = _column(header, "mane_select", "mane")
            i_canon = _column(header, "canonical")
            if i_sym is None:
                log.warning("%s has no gene column; constraint not loaded",
                            Path(p).name)
                return {}

            # Which transcript speaks for the gene. v4.1 carries 211,523
            # transcript rows for 18,204 genes, so this is a choice made
            # about twelve times per gene rather than an occasional
            # tie-break.
            #
            # Keeping the lowest LOEUF looks neutral and is not: a
            # minimum over a dozen draws lands below the typical one, so
            # every multi-transcript gene drifts toward "constrained".
            # On the real file 9,207 of 18,204 genes differ from their
            # MANE Select value that way, 922 of them by more than 0.2,
            # which is wider than the gap between LOEUF_SENSITIVE and
            # LOEUF_TOLERANT. gnomAD applies its own LOEUF rank and
            # decile to MANE Select only, falling back to canonical, so
            # that is the order used here. The lowest-LOEUF rule
            # survives as the last resort, for the 644 genes carrying
            # neither flag.
            best: dict = {}
            for row in rdr:
                if len(row) <= i_sym:
                    continue
                sym = row[i_sym].strip().upper()
                if not sym:
                    continue

                def num(i):
                    if i is None or len(row) <= i:
                        return None
                    try:
                        return float(row[i])
                    except ValueError:
                        return None

                pli, loeuf = num(i_pli), num(i_loeuf)
                if pli is None and loeuf is None:
                    continue
                rank = (2 if _flag(row, i_mane)
                        else 1 if _flag(row, i_canon) else 0)
                prev = best.get(sym)
                if prev is None or rank > prev[0] or (
                        rank == prev[0] and loeuf is not None
                        and (prev[2] is None or loeuf < prev[2])):
                    best[sym] = (rank, pli, loeuf)
            out = {sym: (pli, loeuf) for sym, (_, pli, loeuf) in best.items()}
    except OSError as exc:
        log.warning("could not read %s: %s", p, exc)
        return {}

    rel, how = _gnomad_release(p, header)
    _SOURCES[key] = Source("gnomAD constraint", p, rel, how)
    log.info("gnomAD constraint: %d gene(s) from %s (release %s)",
             len(out), Path(p).name, rel or "unrecorded")
    _CACHE[key] = out
    return out


def load_clingen(path=None) -> dict:
    """{SYMBOL: (haploinsufficient, triplosensitive)} from ClinGen."""
    p = resolve(path, DEFAULT_CLINGEN_PATHS)
    if not p:
        return {}
    key = ("clingen", p)
    if key in _CACHE:
        return _CACHE[key]

    out: dict = {}
    try:
        # ClinGen's header is itself a '#' line - five comment lines,
        # then "#Gene Symbol<TAB>Gene ID<TAB>...". Dropping every '#'
        # line before looking for a header therefore threw the header
        # away and read the first gene as one, no column matched, and
        # this returned {} for all 1,700 curated genes without a word.
        # A module whose job is to separate "not assessed" from
        # "tolerant" was reporting the first as the second.
        #
        # So: the last '#' line before the data is the header.
        with open(p, "r", encoding="utf-8", errors="replace") as fh:
            lines = fh.read().splitlines()
        header_line, data = None, []
        for line in lines:
            if line.startswith("#"):
                if not data:
                    header_line = line.lstrip("#")
                continue
            if line.strip():
                data.append(line)
        if header_line is None and data:
            header_line, data = data[0], data[1:]
        if header_line is None:
            return {}
        header = next(csv.reader([header_line], delimiter="\t"))
        rdr = csv.reader(data, delimiter="\t")
        i_sym = _column(header, "gene symbol", "gene_symbol", "symbol",
                        "#gene symbol")
        i_hi = _column(header, "haploinsufficiency score",
                       "haploinsufficiency_score")
        i_ts = _column(header, "triplosensitivity score",
                       "triplosensitivity_score")
        if i_sym is None:
            # Loudly. Silence is what let the bug above live: a caller
            # cannot tell an empty return from a parse failure.
            log.warning("%s: no gene column in header %r; "
                        "ClinGen dosage not loaded",
                        Path(p).name, header[:4])
            return {}
        for row in rdr:
            if len(row) <= i_sym:
                continue
            sym = row[i_sym].strip().upper()
            if not sym:
                continue

            def score(i):
                """True, False, or None - never a number treated as a scale.

                ClinGen's field mixes an evidence grade with two
                category codes that sit above the grade's range, so a
                `>=` test reads both of them as the strongest possible
                evidence:

                  3   sufficient evidence for dosage pathogenicity
                  2   some evidence
                  1   little evidence
                  0   NO EVIDENCE AVAILABLE  -> not assessed
                  30  autosomal recessive    -> not a dosage statement
                  40  DOSAGE SENSITIVITY UNLIKELY -> assessed, negative
                  "Not yet evaluated"        -> not assessed

                0 is the one that matters most here. It is the file
                saying nobody has looked, and returning False for it
                manufactures a measured negative out of an absence -
                for 1,306 of the genes in the current release.
                """
                if i is None or len(row) <= i:
                    return None
                try:
                    v = int(str(row[i]).strip())
                except ValueError:
                    return None          # "Not yet evaluated"
                if v == 30:              # recessive: says nothing about dose
                    return None
                if v == 40:              # curated as insensitive
                    return False
                if v <= 0:               # no evidence available
                    return None
                return v >= 2
            out[sym] = (score(i_hi), score(i_ts))
    except OSError as exc:
        log.warning("could not read %s: %s", p, exc)
        return {}

    rel, how = _clingen_release(lines)
    _SOURCES[key] = Source("ClinGen dosage curation", p, rel, how)
    log.info("ClinGen dosage: %d gene(s) from %s (release %s)",
             len(out), Path(p).name, rel or "unrecorded")
    _CACHE[key] = out
    return out


def assess(symbol: str, *, gnomad_path=None, clingen_path=None
           ) -> Constraint:
    """What is known about dosage tolerance for one gene."""
    sym = (symbol or "").strip().upper()
    c = Constraint(symbol=sym)
    if not sym:
        return c

    pli, loeuf = (load_gnomad(gnomad_path).get(sym) or (None, None))
    hi, ts = (load_clingen(clingen_path).get(sym) or (None, None))
    gp, cp = resolve(gnomad_path), resolve(clingen_path, DEFAULT_CLINGEN_PATHS)
    c.sources = [asdict(src) for src in
                 (_SOURCES.get(("gnomad", gp)), _SOURCES.get(("clingen", cp)))
                 if src is not None]
    c.pli, c.loeuf = pli, loeuf
    c.haploinsufficient, c.triplosensitive = hi, ts
    c.assessed = any(v is not None for v in (pli, loeuf, hi, ts))

    if not c.assessed:
        c.notes.append(
            "dosage was not assessed: no gnomAD constraint or ClinGen "
            "curation on disk for this gene. That is a missing file, not "
            "a finding that the gene tolerates any dose.")
        return c

    if hi:
        c.dosage = "sensitive"
        c.notes.append(
            "ClinGen curates this gene as haploinsufficient: one working "
            "copy is not enough, which is direct evidence that its level "
            "matters.")
        # The population may say the opposite, and a reader shown one
        # verdict cannot tell. On the current releases 44 genes are
        # curated haploinsufficient while carrying more loss-of-function
        # variants than chance predicts. PCSK9 is one of them, and the
        # reason is instructive: its loss has a phenotype, which is what
        # ClinGen records, but the phenotype is low LDL, so selection
        # never removed the carriers.
        if loeuf is not None and loeuf > LOEUF_TOLERANT:
            c.notes.append(
                f"The two sources disagree. LOEUF {loeuf:.2f} says the "
                f"population carries about as many loss-of-function "
                f"variants as chance predicts. A curated phenotype that "
                f"selection did not act on is usually late onset, "
                f"incompletely penetrant, or not harmful - so read what "
                f"ClinGen curated, not just that it did.")
    elif pli is not None and pli >= PLI_SENSITIVE:
        c.dosage = "sensitive"
        c.notes.append(
            f"pLI {pli:.2f}: losing one copy is selected against in the "
            f"population, so the gene's level is under control and a "
            f"partial loss is a plausible mechanism.")
        # 181 genes are constrained in gnomAD while ClinGen records no
        # evidence for haploinsufficiency. Those are different claims:
        # ClinGen asks whether a deletion has been reported to cause
        # disease, gnomAD asks whether carriers are depleted. Absence of
        # the first is routine for a gene nobody has deleted.
        if hi is False:
            c.notes.append(
                "ClinGen reviewed this gene and did not find evidence "
                "for haploinsufficiency. That is a weaker statement "
                "than it sounds beside the constraint above: it is "
                "about whether a deletion has been reported, not about "
                "whether one would be tolerated.")
    elif loeuf is not None and loeuf < LOEUF_SENSITIVE:
        c.dosage = "sensitive"
        c.notes.append(
            f"LOEUF {loeuf:.2f}: among the most constrained genes, which "
            f"means its dose matters.")
    elif loeuf is not None and loeuf > LOEUF_TOLERANT:
        c.dosage = "tolerant"
        c.notes.append(
            f"LOEUF {loeuf:.2f}: the population carries about as many "
            f"loss-of-function variants as chance predicts, so losing ONE "
            f"copy is tolerated.")
        c.notes.append(
            "That is not evidence the gene is unimportant. A recessive "
            "condition is invisible to this measure entirely - carriers "
            "are healthy, so selection never saw the phenotype. Mosaic "
            "variegated aneuploidy works this way: biallelic BUB1B, with "
            "unaffected heterozygous parents.")

    if ts:
        c.notes.append(
            "ClinGen curates this gene as triplosensitive: too much is "
            "itself harmful. Any plan to raise it needs a ceiling, and "
            "the package cannot predict where that ceiling is.")
    return c


# What the broken allele produces, which decides whether raising both
# copies helps. Not a property of the gene, so it has to be told.
ALLELE_PRODUCTS = ("null", "truncated", "unknown")


def headroom_note(c: Constraint, goal: str, *,
                  one_allele_lost: bool = False,
                  broken_allele_makes: str = "unknown") -> str | None:
    """What this implies for a plan to move the gene in `goal`.

    `one_allele_lost` is the common rare-disease case: one copy is
    broken and the plan is to raise the other to compensate.
    `broken_allele_makes` is what the broken copy produces - "null" when
    it makes nothing (a deletion, or a transcript cleared by
    nonsense-mediated decay), "truncated" when it makes a shortened or
    altered protein that survives.

    Deliberately short of a recommendation. Constraint says what
    selection did about losing the gene; it cannot say whether a drug
    will raise it, by how much, or whether that helps.
    """
    if not c.assessed and not one_allele_lost:
        return None

    if goal == "up" and one_allele_lost:
        bits = [
            "**One copy broken, raising the other.** The target is "
            "roughly a doubling of what the working allele makes, to "
            "restore the normal total. That is at the achievable end as "
            "these things go - far easier than asking for a level the "
            "gene has never had."]

        if c.assessed and c.dosage == "sensitive":
            bits.append(
                "The gene is dosage-sensitive, which supports the "
                "premise: its level is under selection, so a halving "
                "plausibly matters and restoring it plausibly helps.")
        elif c.assessed and c.dosage == "tolerant":
            bits.append(
                "Note the tension: the population tolerates losing one "
                "copy of this gene, so whatever makes your case "
                "symptomatic is not explained by the heterozygous state "
                "alone. The working allele may itself be hypomorphic, "
                "or the broken one may be doing harm rather than merely "
                "being absent.")

        # The hazard that decides whether the strategy works at all.
        bits.append(
            "**Transcriptional upregulation is not allele-specific.** "
            "Nothing in this package acts on one copy: a drug at the "
            "promoter, or at a factor that binds it, raises both.")
        if broken_allele_makes == "null":
            bits.append(
                "With a null allele that is harmless - the broken copy "
                "makes nothing to amplify, so the extra transcription "
                "lands on functional protein.")
        elif broken_allele_makes == "truncated":
            bits.append(
                "That is the problem here: the broken copy makes a "
                "truncated protein, so raising transcription raises the "
                "bad product alongside the good. If it competes with "
                "the working protein the net effect can be worse than "
                "no treatment. An allele-selective approach, or one "
                "aimed at the protein rather than the transcript, "
                "avoids this; raising the promoter does not.")
        else:
            bits.append(
                "Settle this before acting on anything below: does the "
                "broken copy make nothing, or does it make a shortened "
                "protein that survives? If nonsense-mediated decay "
                "clears that transcript the answer is nothing, and "
                "raising both copies is safe. If it does not - a "
                "last-exon stop, an in-frame deletion, a missense - "
                "you would be raising the bad product too. Pass "
                "`--has-ptc` once you know.")

        if c.triplosensitive:
            bits.append(
                "ClinGen curates this gene as triplosensitive, so do "
                "not overshoot. Restoring the normal total is a dose "
                "every healthy carrier demonstrates is tolerated; "
                "exceeding it is not.")
        return " ".join(bits)

    if goal == "up":
        if c.triplosensitive:
            return ("Raising this gene carries a known hazard: ClinGen "
                    "curates it as triplosensitive, so an overshoot is "
                    "not a neutral outcome. Nothing here predicts the "
                    "dose at which that starts.")
        if c.dosage == "sensitive":
            return ("The gene's level is under selection, so a partial "
                    "loss is a plausible mechanism and a modest increase "
                    "is a coherent goal. How much is enough is not "
                    "something this package can tell you.")
        if c.dosage == "tolerant":
            return ("One copy of this gene is dispensable in the "
                    "population, so if your condition is recessive the "
                    "relevant question is how much function the variant "
                    "alleles retain - not what the population tolerates.")
    elif goal == "down":
        if c.dosage == "sensitive":
            return ("This gene is dosage-sensitive, so lowering it is "
                    "not obviously safe: the population shows that "
                    "having less of it is selected against.")
    return None


def available() -> dict:
    return {"gnomad": resolve(None, DEFAULT_GNOMAD_PATHS),
            "clingen": resolve(None, DEFAULT_CLINGEN_PATHS)}


def summarise(c: Constraint, goal: str = "up", *,
              one_allele_lost: bool = False,
              broken_allele_makes: str = "unknown") -> str:
    """One paragraph for the report."""
    if not c.assessed and one_allele_lost:
        return headroom_note(c, goal, one_allele_lost=True,
                             broken_allele_makes=broken_allele_makes) or ""
    if not c.assessed:
        return (f"Dosage tolerance for {c.symbol} was not assessed - no "
                f"gnomAD constraint or ClinGen curation on disk. "
                f"See `constraint` in the documentation for the two "
                f"files, both of which work offline.")
    bits = [f"**{c.symbol} dosage: {c.dosage}.**"]
    bits.extend(c.notes)
    note = headroom_note(c, goal, one_allele_lost=one_allele_lost,
                         broken_allele_makes=broken_allele_makes)
    if note:
        bits.append(note)
    return " ".join(bits)


def reset_cache() -> None:
    _CACHE.clear()
    _SOURCES.clear()
