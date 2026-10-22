"""Translation and mRNA stability, read off the transcript's own UTRs.

Why this exists
---------------
The coverage audit marked two of the six steps to a protein as "listed
only": the report named the mechanism and looked nothing up for the gene.
Both can be asked of sequence this package already fetches, with no new
data source:

- **translation** — an upstream open reading frame in the 5'UTR. A uORF
  captures scanning ribosomes and represses translation of the main ORF
  downstream, so a transcript that carries one has a translational lever
  that a transcript without one does not.
- **mRNA stability** — AU-rich elements in the 3'UTR. The AUUUA
  pentamer, and the UUAUUUAUU nonamer that is the functional unit,
  recruit the decay machinery.

What this is not
----------------
It is an annotation of the transcript, not a drug search. Finding a uORF
says a lever exists; it does not say any approved compound pulls it, and
no public corpus of compound-by-translation or compound-by-half-life
effects is read here. The audit reflects that: the step becomes
*searched* because something gene-specific ran and returned a number, and
the gap column still says no drug was matched to it.

Everything is judged against a dinucleotide shuffle of the same sequence,
for the reason `promoterdna` does it: AUUUA occurs constantly in an
AU-rich 3'UTR by composition alone, and reporting that as a finding would
be reporting the base content.
"""

from __future__ import annotations

import logging
import random
import re

log = logging.getLogger(__name__)

SHUFFLES = 50
MIN_ENRICHMENT_P = 0.05
# A uORF has to clear this to be worth reporting. Very short ones are
# common and weakly repressive.
MIN_UORF_CODONS = 2

ARE_PENTAMER = re.compile(r"(?=ATTTA)")
# The nonamer is the functional element; a lone pentamer is much weaker.
ARE_NONAMER = re.compile(r"(?=TTATTTATT)")

# Four G-runs close enough to fold back on each other. A G-quadruplex in
# the 5'UTR stalls the scanning ribosome, so it represses translation and
# resolving it raises output. The loop bound is the usual 1-7.
G_QUADRUPLEX = re.compile(r"(?=G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,})")

# DRACH, the m6A consensus: D=A/G/T, R=A/G, H=A/C/T. Methylated adenosine
# in the 3'UTR recruits YTHDF readers and shortens half-life. The motif is
# common enough that only a clear excess over shuffle means anything, and
# even then it says where to look rather than where the mark is.
DRACH = re.compile(r"(?=[AGT][AG]AC[ACT])")

# The Pumilio response element. Much more specific than DRACH, and PUM1/2
# binding it shortens half-life, so an excess is a named decay factor
# rather than an anonymous one.
PUMILIO = re.compile(r"(?=TGTA\w ATA)".replace(" ", ""))

# A 5'TOP motif: a C at the cap followed by a pyrimidine run. It is
# positional rather than countable, so it is reported as present or not
# and never tested against a shuffle.
TOP_MIN_PYRIMIDINES = 4
TOP_MAX_PYRIMIDINES = 15


def _shuffle(seq: str, rng) -> str:
    from .promoterdna import _dinucleotide_shuffle
    return _dinucleotide_shuffle(seq, rng)


def find_uorfs(five_prime: str) -> list[dict]:
    """Upstream ORFs: an ATG with an in-frame stop inside the 5'UTR.

    A uORF that runs past the end of the 5'UTR is out-of-frame overlap
    with the main ORF rather than a discrete uORF, and is reported
    separately because it represses differently.
    """
    seq = (five_prime or "").upper().replace("U", "T")
    out = []
    for m in re.finditer(r"(?=ATG)", seq):
        start = m.start()
        for i in range(start + 3, len(seq) - 2, 3):
            codon = seq[i:i + 3]
            if codon in ("TAA", "TAG", "TGA"):
                codons = (i - start) // 3
                if codons >= MIN_UORF_CODONS:
                    out.append({"start": start, "stop": i + 3,
                                "codons": codons, "overlaps_cds": False,
                                "kozak": _kozak(seq, start)})
                break
        else:
            out.append({"start": start, "stop": None,
                        "codons": (len(seq) - start) // 3,
                        "overlaps_cds": True, "kozak": _kozak(seq, start)})
    return out


def _kozak(seq: str, atg: int) -> str:
    """Strong, adequate or weak, by the -3 and +4 positions."""
    minus3 = seq[atg - 3] if atg >= 3 else ""
    plus4 = seq[atg + 3] if len(seq) > atg + 3 else ""
    strong_minus3 = minus3 in ("A", "G")
    if strong_minus3 and plus4 == "G":
        return "strong"
    if strong_minus3 or plus4 == "G":
        return "adequate"
    return "weak"


def _top_motif(five: str) -> dict:
    """A 5'TOP motif at the cap, which puts the transcript under mTOR.

    Only meaningful at the real transcription start. The sequence here is
    the annotated 5'UTR, so a hit is as good as that annotation and no
    better; a transcript whose start site is wrong by a few bases loses
    the motif entirely.
    """
    if not five or five[0] != "C":
        return {"present": False,
                "why": "the annotated 5' end is not a C, which a TOP motif "
                       "requires"}
    run = 0
    for ch in five[1:TOP_MAX_PYRIMIDINES + 1]:
        if ch not in "CT":
            break
        run += 1
    if run < TOP_MIN_PYRIMIDINES:
        return {"present": False,
                "why": f"the 5' end is C followed by only {run} "
                       f"pyrimidine(s); a TOP motif needs "
                       f"{TOP_MIN_PYRIMIDINES}"}
    return {"present": True, "n_pyrimidines": run,
            "motif": five[:run + 1],
            "why": "translation of a TOP transcript is set by mTORC1, so "
                   "its output moves with mTOR signalling rather than with "
                   "its own promoter"}


def _count(pat, seq: str) -> int:
    return sum(1 for _ in pat.finditer(seq))


def _enrichment(observed: int, seq: str, counter, *, shuffles: int,
                seed: int) -> dict:
    rng = random.Random(seed)
    null = [counter(_shuffle(seq, rng)) for _ in range(shuffles)]
    exp = sum(null) / float(len(null)) if null else 0.0
    p = (sum(1 for m in null if m >= observed) + 1) / float(len(null) + 1)
    return {"n_observed": observed, "expected": round(exp, 2),
            "p_value": round(p, 4),
            "enriched": bool(observed > exp and p <= MIN_ENRICHMENT_P)}


def scan(five_prime: str, three_prime: str, *, shuffles: int = SHUFFLES,
         seed: int = 0) -> dict:
    """Both UTRs, each feature against a shuffle of its own sequence."""
    five = (five_prime or "").upper().replace("U", "T")
    three = (three_prime or "").upper().replace("U", "T")
    out: dict = {"five_prime_length": len(five),
                 "three_prime_length": len(three),
                 "shuffles": shuffles}

    if len(five) >= 30:
        uorfs = find_uorfs(five)
        discrete = [u for u in uorfs if not u["overlaps_cds"]]
        out["translation"] = {
            **_enrichment(len(discrete), five,
                          lambda s: len([u for u in find_uorfs(s)
                                         if not u["overlaps_cds"]]),
                          shuffles=shuffles, seed=seed),
            "uorfs": discrete[:12],
            "n_overlapping": sum(1 for u in uorfs if u["overlaps_cds"]),
            "gc": round((five.count("G") + five.count("C")) / len(five), 3),
            "g_quadruplex": _enrichment(
                _count(G_QUADRUPLEX, five), five,
                lambda s: _count(G_QUADRUPLEX, s),
                shuffles=shuffles, seed=seed),
            "top_motif": _top_motif(five),
        }
    else:
        out["translation"] = {"error": f"5'UTR is {len(five)} bp, too short "
                                       f"to scan for upstream ORFs"}

    if len(three) >= 50:
        out["stability"] = {
            "pentamer": _enrichment(_count(ARE_PENTAMER, three), three,
                                    lambda s: _count(ARE_PENTAMER, s),
                                    shuffles=shuffles, seed=seed),
            "nonamer": _enrichment(_count(ARE_NONAMER, three), three,
                                   lambda s: _count(ARE_NONAMER, s),
                                   shuffles=shuffles, seed=seed),
            "au_content": round(
                (three.count("A") + three.count("T")) / len(three), 3),
            "m6a_drach": _enrichment(_count(DRACH, three), three,
                                     lambda s: _count(DRACH, s),
                                     shuffles=shuffles, seed=seed),
            "pumilio": _enrichment(_count(PUMILIO, three), three,
                                   lambda s: _count(PUMILIO, s),
                                   shuffles=shuffles, seed=seed),
        }
    else:
        out["stability"] = {"error": f"3'UTR is {len(three)} bp, too short "
                                     f"to scan for AU-rich elements"}
    return out


DEFAULT_SITES_PATHS = ("data/transcript_sites.tsv",)

# How the public sources map onto the one table. Printed when the table
# is missing, so the reader is told what to fetch rather than that
# nothing is known.
SITES_SOURCES = (
    ("TargetScan", "predicted microRNA sites, conserved families",
     "https://www.targetscan.org/ -> Predicted_Targets_Info.txt"),
    ("miRTarBase", "microRNA targets with experimental support",
     "https://mirtarbase.cuhk.edu.cn/ -> hsa_MTI.xlsx"),
    ("ENCODE eCLIP", "measured RNA-binding-protein footprints, ~150 RBPs",
     "https://www.encodeproject.org/ -> eCLIP IDR peaks, BED"),
    ("m6A atlases", "measured methylation sites rather than DRACH motifs",
     "REPIC or m6A-Atlas"),
)


def annotated_sites(gene: str, *, path=None) -> dict:
    """Named microRNA, RBP or m6A sites on this transcript, from a table.

    Columns, header required, order free, tab or comma:

        gene, factor, kind [, region] [, evidence] [, source]

    `kind` is free text - "microRNA", "RBP", "m6A" - and is used only to
    group the output. `evidence` separates measured sites from predicted
    ones; a row without it is reported as unspecified rather than
    promoted to measured.
    """
    import csv
    from pathlib import Path

    found = next((Path(c) for c in ([path] if path else list(DEFAULT_SITES_PATHS))
                  if c and Path(c).exists()), None)
    if found is None:
        return {}
    want = (gene or "").strip().upper()
    rows: list = []
    try:
        with open(found, encoding="utf-8", errors="replace") as fh:
            sample = fh.read(8192)
            fh.seek(0)
            delim = "\t" if sample.count("\t") >= sample.count(",") else ","
            for row in csv.DictReader(fh, delimiter=delim):
                row = {(k or "").strip().lower(): (v or "").strip()
                       for k, v in row.items()}
                if (row.get("gene") or "").upper() != want:
                    continue
                if not row.get("factor"):
                    continue
                rows.append({"factor": row["factor"],
                             "kind": row.get("kind") or "unspecified",
                             "region": row.get("region") or "",
                             "evidence": row.get("evidence") or "unspecified",
                             "source": row.get("source") or found.name})
    except OSError:
        return {}
    by_kind: dict = {}
    for r in rows:
        by_kind.setdefault(r["kind"], []).append(r)
    return {"gene": want, "path": str(found), "n": len(rows),
            "by_kind": by_kind}


def summarise_sites(sites: dict, gene: str) -> str:
    """Who acts on the transcript, or what to fetch so the run can say."""
    if not sites:
        lines = [f"**No named regulator of the {gene} transcript was read.** "
                 f"The scans above say an element is in excess; they cannot "
                 f"say which factor uses it. Any of these, reshaped into "
                 f"`data/transcript_sites.tsv` with columns "
                 f"gene, factor, kind [, region, evidence, source], would:"]
        for name, what, where in SITES_SOURCES:
            lines.append(f"- **{name}** - {what} ({where})")
        return "\n".join(lines)
    parts = [f"**{sites['n']} named regulator(s) of the {gene} transcript** "
             f"from `{sites['path']}`."]
    for kind, rows in sorted(sites["by_kind"].items()):
        measured = [r for r in rows
                    if r["evidence"].lower() not in ("predicted",
                                                     "unspecified")]
        names = ", ".join(sorted({r["factor"] for r in rows})[:10])
        parts.append(f"{kind}: {names}"
                     + (f" ({len(measured)} of {len(rows)} measured rather "
                        f"than predicted)" if rows else ""))
    return " ".join(parts)


def summarise(s: dict, gene: str = "this gene") -> str:
    """Two sentences a reader can act on, or the reason there are none."""
    if not s:
        return ""
    parts = []
    t = s.get("translation") or {}
    if t.get("error"):
        parts.append(t["error"] + ".")
    else:
        n, exp = t.get("n_observed", 0), t.get("expected", 0)
        if n and t.get("enriched"):
            strong = sum(1 for u in t.get("uorfs") or []
                         if u.get("kozak") == "strong")
            parts.append(
                f"The {s['five_prime_length']} bp 5'UTR carries **{n} "
                f"upstream ORF(s)** against {exp} expected on shuffled "
                f"sequence (p {t['p_value']})"
                + (f", {strong} with a strong Kozak" if strong else "")
                + f". A uORF captures scanning ribosomes, so {gene} has a "
                  f"translational lever - though nothing here searched for "
                  f"a compound that pulls it.")
        elif n:
            parts.append(
                f"The 5'UTR carries {n} upstream ORF(s), which is about "
                f"what shuffled sequence of the same composition gives "
                f"({exp}), so this is not a feature of {gene} in "
                f"particular.")
        else:
            parts.append(f"No upstream ORF in the {s['five_prime_length']} "
                         f"bp 5'UTR, so translational repression by that "
                         f"route is not available here.")
        g4 = t.get("g_quadruplex") or {}
        if g4.get("enriched"):
            parts.append(
                f"It does carry {g4['n_observed']} G-quadruplex motif(s) "
                f"against {g4['expected']} expected (p {g4['p_value']}); a "
                f"quadruplex stalls the scanning ribosome, so resolving it "
                f"is a second translational lever.")
        top = t.get("top_motif") or {}
        if top.get("present"):
            parts.append(
                f"The 5' end is a TOP motif (C plus "
                f"{top['n_pyrimidines']} pyrimidines), so this "
                f"transcript's translation follows mTORC1 - which makes "
                f"mTOR signalling a route to its output, and makes any "
                f"mTOR inhibitor in the ranking below a likely way to "
                f"LOWER it.")
    st = s.get("stability") or {}
    if st.get("error"):
        parts.append(st["error"] + ".")
    else:
        non = st.get("nonamer") or {}
        pen = st.get("pentamer") or {}
        if non.get("enriched") or pen.get("enriched"):
            which = "nonamer" if non.get("enriched") else "pentamer"
            rec = non if non.get("enriched") else pen
            parts.append(
                f"The {s['three_prime_length']} bp 3'UTR is enriched for "
                f"the AU-rich {which}: {rec['n_observed']} against "
                f"{rec['expected']} expected (p {rec['p_value']}). That "
                f"marks the transcript for rapid decay, so its half-life "
                f"is a lever - and again, no compound was matched to it.")
        else:
            parts.append(
                f"The 3'UTR is {st.get('au_content', 0):.0%} AU but holds "
                f"no more AU-rich elements than shuffled sequence of the "
                f"same composition, so there is no specific decay element "
                f"to aim at.")
        pum = st.get("pumilio") or {}
        if pum.get("enriched"):
            parts.append(
                f"It is enriched for the Pumilio response element "
                f"({pum['n_observed']} against {pum['expected']}, p "
                f"{pum['p_value']}), which names PUM1/PUM2 as a decay "
                f"factor acting on it.")
        m6a = st.get("m6a_drach") or {}
        if m6a.get("enriched"):
            parts.append(
                f"DRACH, the m6A consensus, is also in excess "
                f"({m6a['n_observed']} against {m6a['expected']}). DRACH "
                f"is a weak signal on its own - most occurrences are never "
                f"methylated - so this points at an m6A atlas rather than "
                f"standing as a finding.")
    return " ".join(parts)
