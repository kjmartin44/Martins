"""Small molecules aimed at the promoter DNA rather than at a protein.

The route this covers
---------------------
Every other route in this package acts on a protein: a transcription
factor, something upstream of one, or the gene's own product. This one
acts on the promoter itself. Some small molecules bind DNA with a
sequence preference, and where a promoter carries the feature they prefer
they can block or stabilise what happens there. Mithramycin occupying
GC-boxes and displacing SP1 is the worked example.

What is read off the sequence
-----------------------------
Four features, all computed from the promoter this run already has, with
no network and no structure:

  G-quadruplex   four runs of three or more Gs with short loops, on
                 either strand. A G4 in a promoter can hold transcription
                 open or shut depending on where it sits, and a ligand
                 that stabilises it freezes that state.
  i-motif        the C-rich complement of the same thing.
  GC-box         the SP1 consensus, which is what mithramycin occupies.
  AT tract       a run of A/T wide enough for a minor-groove binder.

Why a count alone is not a finding
----------------------------------
G-quadruplex motifs are common in human promoters, so finding one says
little on its own. Every count here is therefore compared against the
same promoter dinucleotide-shuffled, which preserves base and
dinucleotide composition while destroying arrangement. A GC-rich promoter
will produce G4 motifs by composition alone and the shuffle says so.

What this cannot tell you
-------------------------
**A motif is not a structure.** The sequence can fold into a quadruplex;
whether it does in your cells, at your locus, in that chromatin state, is
not something a regex knows.

**Sequence preference is not selectivity.** These agents prefer a
feature, not your promoter. Any other promoter carrying the same feature
is bound too, and most of the genome's G4s are not in your gene.

**Nearly every approved DNA binder is cytotoxic.** That is not incidental
- binding everyone's DNA is how they work as chemotherapy. The agents
below are flagged with their systemic class, and for a chronic target
that flag is usually the end of the conversation.

**Direction is not addressed.** Stabilising a G4 in one promoter
represses it and in another relieves repression, depending on the strand
and the position relative to the TSS. Nothing here predicts which.
"""

from __future__ import annotations

import logging
import random
import re
from dataclasses import dataclass, field

log = logging.getLogger(__name__)

# Four G-runs of >=3, separated by loops of 1-7. The standard permissive
# pattern; a stricter one misses the long-loop quadruplexes that fold
# perfectly well.
G4 = re.compile(r"(?=(G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,}))")
IMOTIF = re.compile(r"(?=(C{3,}\w{1,7}C{3,}\w{1,7}C{3,}\w{1,7}C{3,}))")
# SP1's GC-box. The classical consensus, which is what mithramycin sits in.
GC_BOX = re.compile(r"(?=(GGGGCGGGG|GGGCGG|CCGCCC))")
AT_TRACT = re.compile(r"(?=([AT]{6,}))")

SHUFFLES = 200          # dinucleotide shuffles for the null
MIN_ENRICHMENT_P = 0.05


@dataclass
class DNAAgent:
    """A small molecule with a documented DNA sequence preference."""
    name: str
    prefers: str            # which feature key it needs
    mechanism: str
    approved: bool
    note: str


# Hand-curated and deliberately short. Membership is a claim that the
# agent has a documented sequence preference, not that it is a sensible
# therapy - most of these are chemotherapy and the burden table says so.
AGENTS: tuple[DNAAgent, ...] = (
    DNAAgent(
        "Mithramycin", "gc_box",
        "occupies GC-rich minor groove and displaces SP1 from GC-boxes",
        True,
        "the worked example of a promoter-directed small molecule, and "
        "also a cytotoxic antibiotic withdrawn from general oncology use "
        "for toxicity"),
    DNAAgent(
        "Plicamycin", "gc_box",
        "the same agent under its other name; displaces SP1 from GC-boxes",
        True,
        "identical molecule to mithramycin, listed because sources name "
        "it either way"),
    DNAAgent(
        "Actinomycin D", "gc_box",
        "intercalates at GpC steps and blocks elongation",
        True,
        "sequence preference is real and weak; it stops transcription "
        "broadly and is dosed as chemotherapy"),
    DNAAgent(
        "Mitoxantrone", "g4",
        "intercalates and stabilises G-quadruplex structures",
        True,
        "a topoisomerase II poison first; the G4 affinity is documented "
        "but is not why it is given"),
    DNAAgent(
        "Doxorubicin", "g4",
        "intercalates, with reported quadruplex affinity",
        True,
        "an anthracycline whose clinical effect is DNA damage and "
        "topoisomerase poisoning, not promoter selectivity"),
    DNAAgent(
        "Berberine", "g4",
        "binds and stabilises G-quadruplexes",
        False,
        "a plant alkaloid sold as a supplement, not an approved drug; "
        "included because it is the G4 binder people most often have "
        "access to"),
    DNAAgent(
        "Pyridostatin", "g4",
        "selective G-quadruplex stabiliser",
        False,
        "a tool compound and the reference G4 ligand in the literature; "
        "no clinical use"),
    DNAAgent(
        "Quarfloxin", "g4",
        "disrupts a quadruplex-nucleolin complex",
        False,
        "reached phase 2 in oncology and was not approved; the only G4 "
        "ligand to get that far"),
    DNAAgent(
        "Telomestatin", "g4",
        "macrocyclic G-quadruplex stabiliser",
        False,
        "a natural product and a research reagent"),
    DNAAgent(
        "Pentamidine", "at_tract",
        "binds AT-rich minor groove",
        True,
        "approved as an antiparasitic; the DNA binding is real and is "
        "not what it is prescribed for"),
    DNAAgent(
        "Netropsin", "at_tract",
        "classical AT-selective minor-groove binder",
        False,
        "a reference compound, never a drug"),
    DNAAgent(
        "Distamycin", "at_tract",
        "classical AT-selective minor-groove binder",
        False,
        "a reference compound, never a drug"),
)


@dataclass
class Feature:
    key: str
    label: str
    n_observed: int = 0
    positions: list = field(default_factory=list)
    expected: float = 0.0
    p_value: float | None = None
    enriched: bool = False

    def to_dict(self) -> dict:
        return {"key": self.key, "label": self.label,
                "n_observed": self.n_observed,
                "positions": self.positions[:20],
                "expected": round(self.expected, 2),
                "p_value": self.p_value, "enriched": self.enriched}


def _revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def _dinucleotide_shuffle(seq: str, rng: random.Random) -> str:
    """Shuffle preserving dinucleotide composition, by an Eulerian walk.

    A plain shuffle preserves only base composition, and G-runs are a
    dinucleotide property - GG is what makes a quadruplex possible. A
    mononucleotide null would therefore call almost every GC-rich
    promoter enriched.
    """
    if len(seq) < 4:
        return seq
    edges: dict = {}
    for a, b in zip(seq, seq[1:]):
        edges.setdefault(a, []).append(b)
    for v in edges.values():
        rng.shuffle(v)
    out = [seq[0]]
    cur = seq[0]
    for _ in range(len(seq) - 1):
        nxt = edges.get(cur)
        if not nxt:
            # Walk stranded: fall back to a remaining base so the length
            # is preserved rather than silently returning a short string.
            left = [k for k, v in edges.items() if v]
            if not left:
                break
            cur = rng.choice(left)
            out.append(cur)
            continue
        cur = nxt.pop()
        out.append(cur)
    return "".join(out)


def _count(pattern, seq: str) -> tuple:
    hits = [(m.start(), m.group(1)) for m in pattern.finditer(seq)]
    # Overlapping lookahead matches would count one quadruplex several
    # times; keep the first match of each non-overlapping block.
    kept, last_end = [], -1
    for start, text in hits:
        if start >= last_end:
            kept.append(start)
            last_end = start + len(text)
    return len(kept), kept


def _scannable_runs(sequence: str) -> list:
    """`(offset, run)` for each stretch of unmasked, unambiguous sequence.

    Soft-masked (lowercase) bases are repeats, and the motif scan already
    drops them - `pwm.encode` treats lowercase as N under `skip_repeats`,
    because Alu and simple repeats are a major source of spurious hits for
    GC-rich and homopolymer-like patterns. Those are exactly the patterns
    here, so uppercasing the window put them all back, and the shuffled
    null cannot catch it: the null is built from the same string.

    Split rather than deleted, so no feature is called across a junction
    that does not exist in the genome. `[ACGT]` excludes lowercase and N
    at once, which is why nothing is uppercased on the way in.
    """
    return [(m.start(), m.group(0))
            for m in re.finditer("[ACGT]+", sequence or "")]


def _count_runs(pattern, runs) -> tuple:
    """`_count` over several runs, with positions in the original frame."""
    n, pos = 0, []
    for off, run in runs:
        c, p = _count(pattern, run)
        n += c
        pos.extend(off + i for i in p)
    return n, pos


def scan_features(sequence: str, *, shuffles: int = SHUFFLES,
                  seed: int = 0) -> dict:
    """Sequence features of one promoter, each against a shuffled null."""
    runs = _scannable_runs(sequence)
    total = sum(len(r) for _, r in runs)
    out: dict = {}
    if total < 50:
        return {"error": f"promoter sequence has {total} bp of unmasked, "
                         f"unambiguous sequence, too short to scan for "
                         f"DNA-binding features"}

    rc_runs = [(0, _revcomp(r)) for _, r in runs]
    gc = sum(r.count("G") + r.count("C")
             for _, r in runs) / float(total)

    specs = [
        ("g4", "G-quadruplex", G4, True),
        ("imotif", "i-motif (C-rich)", IMOTIF, False),
        ("gc_box", "SP1 GC-box", GC_BOX, False),
        ("at_tract", "AT tract (>=6)", AT_TRACT, False),
    ]
    rng = random.Random(seed)
    # Each run shuffled on its own, so the null has the same run structure
    # as the observation and the two are comparable.
    shuffled = [[(0, _dinucleotide_shuffle(r, rng)) for _, r in runs]
                for _ in range(shuffles)]

    feats: dict = {}
    for key, label, pat, both_strands in specs:
        n, pos = _count_runs(pat, runs)
        if both_strands:
            n_rc, _ = _count_runs(pat, rc_runs)
            n += n_rc
        null = []
        for sh in shuffled:
            m, _ = _count_runs(pat, sh)
            if both_strands:
                m2, _ = _count_runs(pat, [(0, _revcomp(r)) for _, r in sh])
                m += m2
            null.append(m)
        exp = sum(null) / float(len(null)) if null else 0.0
        # One-sided empirical p with the usual +1 correction, so a feature
        # never comes back at p = 0 from a finite null.
        p = (sum(1 for m in null if m >= n) + 1) / float(len(null) + 1)
        feats[key] = Feature(key=key, label=label, n_observed=n,
                             positions=pos, expected=exp, p_value=round(p, 4),
                             enriched=bool(n > exp and p <= MIN_ENRICHMENT_P))

    out["gc_content"] = round(gc, 3)
    out["length"] = total
    out["features"] = {k: f.to_dict() for k, f in feats.items()}
    out["shuffles"] = shuffles
    return out


def candidates(scan: dict, *, require_enriched: bool = True) -> list:
    """Agents whose preferred feature this promoter actually carries.

    `require_enriched` keeps only features that beat the shuffled null.
    Without it a GC-rich promoter matches the GC-box agents purely by
    composition, which is the finding this module exists to avoid
    reporting.
    """
    if not scan or scan.get("error"):
        return []
    feats = scan.get("features") or {}
    out = []
    for a in AGENTS:
        f = feats.get(a.prefers)
        if not f or not f["n_observed"]:
            continue
        if require_enriched and not f["enriched"]:
            continue
        out.append({
            "drug": a.name, "prefers": a.prefers,
            "feature": f["label"], "n_sites": f["n_observed"],
            "expected": f["expected"], "p_value": f["p_value"],
            "mechanism": a.mechanism, "approved": a.approved,
            "note": a.note,
        })
    return out


def summarise(scan: dict, hits: list) -> str:
    """Prose, including the case where the promoter carries nothing."""
    if not scan:
        return "the promoter sequence was not scanned for DNA features."
    if scan.get("error"):
        return scan["error"]

    feats = scan.get("features") or {}
    enriched = [f["label"] for f in feats.values() if f["enriched"]]
    present = [f"{f['label']} x{f['n_observed']}"
               for f in feats.values() if f["n_observed"]]

    head = (f"{scan['length']} bp at {scan['gc_content']:.0%} GC. "
            f"Present: {', '.join(present) or 'none of the four'}.")
    if not enriched:
        return (head + " None of them is more common than in "
                f"{scan['shuffles']} dinucleotide shuffles of this same "
                f"sequence, so what is here is explained by base "
                f"composition. No DNA-directed agent is proposed: a "
                f"feature a promoter has only by composition is a feature "
                f"most other promoters have too.")

    body = (head + f" Enriched against {scan['shuffles']} shuffles of the "
            f"same sequence: {', '.join(enriched)}.")
    if not hits:
        return body + " No curated agent prefers those features."
    approved = [h["drug"] for h in hits if h["approved"]]
    return (
        body + f" {len(hits)} agent(s) have a documented preference for "
        f"them" + (f", of which {len(approved)} are approved "
                   f"({', '.join(approved)})" if approved else
                   ", none of them approved") +
        ". Sequence preference is not selectivity: these bind the feature "
        "wherever it occurs, and most of the genome's copies are not in "
        "your gene. Nearly every approved member is cytotoxic, which is "
        "how it works rather than a side effect.")
