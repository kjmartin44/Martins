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
        }
    else:
        out["stability"] = {"error": f"3'UTR is {len(three)} bp, too short "
                                     f"to scan for AU-rich elements"}
    return out


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
    return " ".join(parts)
