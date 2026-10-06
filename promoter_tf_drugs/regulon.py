"""How many genes does this TF control? The specificity that actually matters.

The drug-level `selectivity` score asks how many *proteins* a compound binds.
That is the wrong denominator for this pipeline. Suppose you find a compound
that hits exactly one target, with 1 nM affinity, no off-targets at all - a
perfect score on molecular selectivity. If that target is SP1, the drug will
move tens of thousands of promoters. Your gene is one of them.

So specificity here has two independent layers, and a drug needs both:

  molecular       how many proteins does the drug bind?      (drug_scoring)
  transcriptional how many genes does that protein control?  (this module)

A narrow TF is worth a great deal. Getting your gene via a factor with 40
targets is a far cleaner experiment - and a far more plausible therapy - than
getting it via a global regulator, even when the motif score is worse.

Data sources, best first
------------------------
  TRRUST v2   curated human TF-target pairs, one small TSV.
              https://www.grnpedia.org/trrust/  -> trrust_rawdata.human.tsv
              Pass with --regulon-file. This is the recommended route.
  DoRothEA    confidence-tiered regulons (A-E); A/B only is a good filter.
  ChIP-Atlas  measured target genes per TF, larger and noisier.

Without a file, a bundled qualitative tier table is used. Those tiers are
deliberately coarse - "broad" vs "narrow", not invented counts - because a
precise-looking fake number is worse than an honest bucket.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

log = logging.getLogger(__name__)

# Coarse tiers for well-characterised factors. Used only when no regulon file
# is supplied. Approximate order of magnitude, not counts.
TIER_TARGETS = {"global": 5000, "broad": 1200, "moderate": 300, "narrow": 60}

GLOBAL_REGULATORS = {
    "SP1", "SP3", "MYC", "MAX", "TP53", "NFKB1", "RELA", "CREB1", "E2F1", "E2F4",
    "JUN", "JUNB", "FOS", "EGR1", "YY1", "CTCF", "TBP", "NFYA", "NFYB", "ELK1",
    "ETS1", "USF1", "USF2", "SRF", "ATF1", "ATF2", "GABPA", "KLF4", "KLF5",
}
BROAD_REGULATORS = {
    "NR3C1", "NR3C2", "ESR1", "ESR2", "AR", "PGR", "PPARG", "PPARA", "RARA",
    "RXRA", "HIF1A", "EPAS1", "STAT1", "STAT3", "STAT5A", "SMAD2", "SMAD3",
    "SREBF1", "SREBF2", "FOXO1", "FOXO3", "NFE2L2", "AHR", "CEBPB", "TCF7L2",
    "RUNX1", "SOX2", "POU5F1", "HSF1",
}
MODERATE_REGULATORS = {
    "VDR", "HNF4A", "HNF1A", "THRA", "THRB", "NR1H3", "NR1H4", "NR1I2", "NR1I3",
    "ESRRA", "REST", "GATA1", "GATA3", "GATA4", "MEF2C", "MYOD1", "PAX6",
    "NEUROD1", "TFEB", "NR4A1", "RORC", "PPARD", "NR1D1",
}


def tier_for(symbol: str) -> str:
    s = symbol.upper()
    if s in GLOBAL_REGULATORS:
        return "global"
    if s in BROAD_REGULATORS:
        return "broad"
    if s in MODERATE_REGULATORS:
        return "moderate"
    return "unknown"


# ------------------------------------------------------------- file loading
def describe_source(path, *, sample: int = 200) -> tuple:
    """`(kind, phrase)` for a regulon file, read from its evidence column.

    The two tables this accepts are not the same kind of claim. TRRUST
    names a PubMed id: a person read that paper and asserted the factor
    regulates the gene. ReMap names cell lines: the factor was found
    BOUND near the gene, which is a different and much weaker statement,
    because most occupancy has no measured effect on transcription.

    Calling both "curated" promotes one of them a tier, in the line of
    the report whose job is to say where the evidence came from.
    """
    path = Path(path)
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    pmid = chip = other = 0
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for i, line in enumerate(fh):
                if i >= sample:
                    break
                parts = line.rstrip("\n").split(sep)
                if len(parts) < 4:
                    other += 1
                    continue
                ev = parts[3].strip().strip('"')
                if ev.upper().startswith("REMAP"):
                    chip += 1
                elif ev.replace(";", "").replace(",", "").strip().isdigit():
                    pmid += 1
                else:
                    other += 1
    except OSError:
        return "unknown", "targets from a regulon file"

    if chip > pmid and chip > other:
        return "chip", ("TFs with targets called from ChIP occupancy. "
                        "A peak near a gene is binding, not regulation")
    if pmid > chip and pmid > other:
        return "curated", "TFs with literature-curated targets"
    return "unknown", ("TFs with targets from a regulon file of "
                       "unrecorded provenance")


def load_regulons(path: str | Path) -> dict[str, set]:
    """Parse TRRUST-style TSV: TF <tab> target <tab> mode <tab> pmid.

    Also accepts a plain two-column TF/target TSV, and DoRothEA CSV exports
    with tf/target columns.
    """
    path = Path(path)
    regulons: dict[str, set] = defaultdict(set)
    sep = "," if path.suffix.lower() == ".csv" else "\t"

    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        first = fh.readline()
        cols = [c.strip().strip('"').lower() for c in first.split(sep)]
        has_header = "tf" in cols or "source" in cols
        if has_header:
            tf_i = cols.index("tf") if "tf" in cols else cols.index("source")
            tg_i = (cols.index("target") if "target" in cols
                    else cols.index("gene") if "gene" in cols else 1)
        else:
            tf_i, tg_i = 0, 1
            fh.seek(0)

        for line in fh:
            parts = [p.strip().strip('"') for p in line.rstrip("\n").split(sep)]
            if len(parts) <= max(tf_i, tg_i):
                continue
            tf, target = parts[tf_i].upper(), parts[tg_i].upper()
            if tf and target and tf != target:
                regulons[tf].add(target)

    log.info("loaded regulons for %d TFs from %s", len(regulons), path.name)
    return dict(regulons)


def regulon_size(symbol: str, regulons: dict | None = None) -> tuple[int | None, str]:
    """(n_target_genes, provenance). None when genuinely unknown."""
    s = symbol.upper()
    if regulons and s in regulons:
        return len(regulons[s]), "regulon file"
    tier = tier_for(s)
    if tier != "unknown":
        return TIER_TARGETS[tier], f"bundled tier: {tier}"
    return None, "unknown"


# -------------------------------------------------------------- the scoring
def transcriptional_specificity(symbol: str, regulons: dict | None = None
                                ) -> tuple[float, int | None, str | None]:
    """0..1. A TF controlling ~50 genes scores ~1.0; ~5000 scores ~0.2.

    Log-scaled, because the difference between 50 and 500 targets matters far
    more than between 5000 and 5500.
    """
    import math

    n, provenance = regulon_size(symbol, regulons)
    if n is None:
        return 0.55, None, "regulon size unknown; scored neutral"
    if n <= 50:
        return 1.0, n, None

    score = max(0.15, 1.0 - 0.42 * math.log10(n / 50.0))
    note = None
    if n >= 2000:
        note = (f"{symbol} is a global regulator (~{n} target genes); any drug "
                "acting through it will move much more than your gene")
    elif n >= 600:
        note = f"{symbol} controls ~{n} genes; expect broad transcriptional effects"
    return round(score, 3), n, note


# Below this many curated targets, the count says more about how much anyone
# has studied the factor than about how many genes it controls. TRRUST is
# literature-derived: a TF with 8 recorded targets is far likelier to be
# under-assayed than to be genuinely narrow.
MIN_CREDIBLE_REGULON = 15


@dataclass
class RegulonEvidence:
    """How much the regulon count can be *trusted*, as distinct from its value.

    `transcriptional_specificity` answers "how broad is this factor". That is
    the wrong question to ask of a number derived from curation effort. A drug
    acting through a factor nobody has studied is not thereby specific — it is
    unmeasured, and the two must not score alike.

    This is the same correction `drug_scoring.evidence_depth` already applies
    to protein selectivity, which the README states plainly: "A low target
    count can mean selective, or can mean nobody looked."
    """
    symbol: str
    n: int | None
    provenance: str
    understudied: bool = False
    trustworthy: bool = True
    note: str | None = None

    def to_dict(self) -> dict:
        return {"symbol": self.symbol, "n": self.n,
                "provenance": self.provenance,
                "understudied": self.understudied,
                "trustworthy": self.trustworthy, "note": self.note}


def regulon_evidence(symbol: str, regulons: dict | None = None
                     ) -> RegulonEvidence:
    """Whether this factor's regulon count is usable as specificity evidence.

    Surfaces the provenance that `transcriptional_specificity` computes and
    throws away — "regulon file", "bundled tier: broad" and "unknown" are
    three different epistemic states that a single float collapses into one.
    """
    n, provenance = regulon_size(symbol, regulons)

    if n is None:
        return RegulonEvidence(
            symbol, None, provenance, trustworthy=False,
            note=("no regulon information for this factor; its breadth is "
                  "unmeasured, which is not the same as narrow"))

    if provenance.startswith("bundled tier"):
        return RegulonEvidence(
            symbol, n, provenance, trustworthy=False,
            note=(f"{n} is an order-of-magnitude bucket, not a count — supply "
                  f"--regulon-file for a real number"))

    if n < MIN_CREDIBLE_REGULON:
        return RegulonEvidence(
            symbol, n, provenance, understudied=True, trustworthy=False,
            note=(f"only {n} curated target(s) for {symbol}. Below "
                  f"{MIN_CREDIBLE_REGULON} this measures curation depth rather "
                  f"than biology — the factor is UNDERSTUDIED, which must not "
                  f"be read as specific"))

    return RegulonEvidence(symbol, n, provenance, trustworthy=True)
