"""Did the compound move the gene, or just move the cell cycle?

The confound this exists to measure
-----------------------------------
A mitotic gene's transcript rises and falls with the cell-cycle phase. BUB1B,
CCNB1, PLK1 and the rest are low in G1 and high at G2/M, so the share of cells
in G2/M sets the measured mRNA level of all of them together. Anything that
pushes a culture into cycle raises every one; anything that arrests it lowers
every one. Neither has touched the promoter.

That is not a hypothetical. In the CTD records for BUB1B, the chemicals that
LOWER it are doxorubicin, etoposide, methotrexate, fluorouracil and dasatinib
- antiproliferatives - and the chemicals that RAISE it are tumour promoters,
mitogens and carcinogens: TPA, croton oil, phenobarbital, estradiol,
methylcholanthrene, dimethylhydrazine. The dominant determinant of BUB1B mRNA
across the whole published record is the proliferation index, not regulation
at its promoter.

`cellcycle.MITOTIC_DISRUPTORS` already halves the score of compounds on a
hand-written list. This measures the same thing instead of listing it: a
compound that moves this gene and also moves the rest of the mitotic
programme the SAME way is, on the evidence available, moving the programme.

What it can and cannot say
--------------------------
It reads LINCS consensus set membership (the Enrichr .gmt files), which
records that a compound moved a gene and not by how much. So this is a
co-occurrence test, not a regression: it says the gene travelled with its
programme, which is what a confound looks like, and cannot say what share of
the change the programme accounts for. A compound that genuinely acts at the
promoter of a cell-cycle gene AND arrests the cycle would also be flagged;
that is the right error to make, because the two are not separable in this
data and only a phase-resolved experiment separates them.

A compound absent from the files gets no verdict, which is not a pass.
"""

from __future__ import annotations

import logging
from pathlib import Path

log = logging.getLogger(__name__)

# The G2/M transcriptional programme, as a readout of what share of the
# culture is dividing. Deliberately the well-measured core rather than an
# exhaustive list: these are L1000 landmarks or close, so a consensus
# signature that moved the programme almost always names several of them.
MARKERS = frozenset("""
MKI67 CCNB1 CCNB2 CCNA2 CDK1 PLK1 TOP2A AURKA AURKB BUB1 BUB1B TTK CDC20
UBE2C CENPF KIF11 KIF20A NUSAP1 TPX2 PTTG1 BIRC5 CKS2 RRM2 TYMS MCM2 PCNA
""".split())

# A signature has to name at least this many markers before it counts as
# having moved the programme. One is noise; three is a pattern.
MIN_MARKERS = 3
# A compound is called out when the share of ITS signatures that move the
# programme sits in this percentile of all compounds. Measured against the
# background rather than an absolute share, because the background is not
# zero: across 33,133 LINCS signatures, 1.9% move three or more of these
# markers on their own. An absolute threshold would have called that
# signal; a percentile calls what stands out from it.
PROGRAMME_PERCENTILE = 95
# Below this many signatures a share is one or two observations.
MIN_SIGNATURES = 5
# And a floor, so a quiet file cannot make an unremarkable compound the
# 95th percentile of nothing.
MIN_SHARE = 0.05


def _drug_of(name: str) -> str:
    from .signatures import _parse_sig_name, breadth_key
    return breadth_key(_parse_sig_name(name)["drug"])


def scan_gmt(path, *, gene: str, wanted: set | None = None) -> dict:
    """One pass: per compound, how often it moves the mitotic programme.

    `gene` is excluded from the marker set, so a compound is never counted
    as confounded on the strength of the very gene being asked about.
    `wanted`, when given, limits the result to those compound keys.
    """
    target = (gene or "").strip().upper()
    markers = MARKERS - {target}
    per_drug: dict = {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            key = _drug_of(fields[0])
            if wanted is not None and key not in wanted:
                continue
            genes = {g.strip().upper() for g in fields[2:] if g.strip()}
            rec = per_drug.setdefault(key, {"n_signatures": 0,
                                            "n_with_programme": 0,
                                            "markers": set()})
            rec["n_signatures"] += 1
            hit = genes & markers
            if len(hit) >= MIN_MARKERS:
                rec["n_with_programme"] += 1
                rec["markers"] |= hit
    return per_drug


def _cutoff(per_drug: dict) -> tuple:
    """`(share to beat, background share)` from the file's own distribution.

    Self-calibrating, the same way `signatures.BROAD_PERCENTILE` is: the
    question "does this compound move the programme unusually often" has
    no answer until "how often does anything" is known.
    """
    shares = sorted(r["n_with_programme"] / r["n_signatures"]
                    for r in per_drug.values()
                    if r["n_signatures"] >= MIN_SIGNATURES)
    total = sum(r["n_signatures"] for r in per_drug.values())
    moved = sum(r["n_with_programme"] for r in per_drug.values())
    background = (moved / total) if total else 0.0
    if not shares:
        return 1.0, background
    i = min(len(shares) - 1,
            int(len(shares) * PROGRAMME_PERCENTILE / 100))
    return max(shares[i], MIN_SHARE), background


def assess(gene: str, drugs, *, direction: str = "up",
           up_path=None, down_path=None) -> dict:
    """`{drug key: verdict}` for the compounds named in `drugs`.

    `direction` is the way the run wants the gene to move; the same
    direction's signature file is read, because the confound is the gene
    and its programme moving TOGETHER.

    The whole file is scanned, not only the compounds asked about, because
    the verdict is a comparison against what every other compound does.
    """
    from .signatures import DEFAULT_GMT_PATHS, breadth_key

    want = direction if direction in ("up", "down") else "up"
    explicit = {"up": up_path, "down": down_path}[want]
    cands = [explicit] if explicit else list(DEFAULT_GMT_PATHS[want])
    path = next((Path(c) for c in cands if c and Path(c).exists()), None)
    if path is None:
        return {}

    wanted = {breadth_key(d) for d in drugs if d}
    if not wanted:
        return {}
    everything = scan_gmt(path, gene=gene)
    cut, background = _cutoff(everything)

    out: dict = {}
    for key in wanted:
        rec = everything.get(key)
        if not rec or not rec["n_signatures"]:
            continue
        n = rec["n_signatures"]
        share = rec["n_with_programme"] / n
        out[key] = {
            "n_signatures": n,
            "n_with_programme": rec["n_with_programme"],
            "share": round(share, 3),
            "background_share": round(background, 4),
            "cutoff": round(cut, 3),
            "percentile": PROGRAMME_PERCENTILE,
            "markers": sorted(rec["markers"])[:8],
            "tracks_proliferation": (n >= MIN_SIGNATURES and share >= cut),
            "direction": want,
        }
    log.info("proliferation check on %s: %d of %d compound(s) move the "
             "mitotic programme with it (cut %.3f, background %.4f)", gene,
             sum(1 for v in out.values() if v["tracks_proliferation"]),
             len(out), cut, background)
    return out


def note_for(rec: dict, gene: str) -> str:
    """The sentence that goes on a route whose compound tracks the cycle."""
    return (
        f"moves the mitotic programme with it: in "
        f"{rec['n_with_programme']} of {rec['n_signatures']} LINCS "
        f"signature(s) ({rec['share']:.0%}, against "
        f"{rec['background_share']:.1%} for signatures at large) this "
        f"compound also moves {', '.join(rec['markers'][:4])} the same way. "
        f"{gene} is a cell-cycle gene, so this is what a change in the "
        f"share of dividing cells looks like, and not evidence of anything "
        f"at the promoter. Read any follow-up from phase-resolved or "
        f"synchronised cells.")


def summarise(assessed: dict, gene: str) -> str:
    """One paragraph for the report."""
    if not assessed:
        return ""
    n = sum(1 for v in assessed.values() if v["tracks_proliferation"])
    any_rec = next(iter(assessed.values()))
    if not n:
        return (f"**No candidate compound moved the mitotic programme with "
                f"{gene} unusually often.** Each was checked against "
                f"{len(MARKERS) - 1} other cell-cycle transcripts across its "
                f"LINCS signatures, and compared with the "
                f"{any_rec['background_share']:.1%} of all signatures that "
                f"move the programme anyway. None reached the "
                f"{any_rec['percentile']}th percentile, so none of these "
                f"rows stands out as a change in how many cells are "
                f"dividing. That is not a clean bill: the test needs the "
                f"compound to be in the files at all.")
    return (
        f"**{n} of {len(assessed)} checked compound(s) move the mitotic "
        f"programme along with {gene}.** {gene}'s transcript tracks the "
        f"cell-cycle phase, so a compound that moves the rest of the "
        f"programme the same way is most simply read as changing how many "
        f"cells are dividing rather than as acting on this gene. Their "
        f"scores are halved and each says so. The bar is the "
        f"{any_rec['percentile']}th percentile of all compounds, against a "
        f"background of {any_rec['background_share']:.1%} of signatures "
        f"moving the programme anyway. This is a measurement over LINCS set "
        f"membership, not a list of known antimitotics; a compound the "
        f"files do not cover gets no verdict, which is not a pass.")
