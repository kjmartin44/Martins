"""Other perturbation corpora, and what agreement between them is worth.

Why more than LINCS
-------------------
`tfactivity.drugs_moving_factor` asks whether a compound's own expression
signature moves a factor's regulon. Run against LINCS alone that is one
platform, one set of mostly cancer cell lines, one dose regime and one
laboratory pipeline. A compound can look convincing there for reasons
peculiar to that setup: L1000 measures 978 genes and infers the rest, the
lines are transformed and proliferating, and the exposures are hours at
micromolar.

Every systematic error of that kind is shared by all 33,000 LINCS
signatures, so no amount of internal replication detects it. The only
thing that does is a corpus built differently. These are those corpora.

What each one is, and what is wrong with it
-------------------------------------------
**CREEDS** is crowd-curated from GEO: real published experiments,
extracted by hand. Different platforms, tissues and labs, which is the
point, and heterogeneous quality with few replicates per signature, which
is the cost.

**DrugMatrix** is toxicogenomics - rat tissue dosed in vivo, mostly liver.
In vivo and whole-organ is a genuinely different measurement from a cell
line, and it is **rat**, so every gene symbol crosses a species boundary
and liver dominates.

**TG-GATEs** is the same idea from a different consortium, with rat in
vivo and both rat and human primary hepatocytes. The human hepatocyte arm
is the most comparable to a human therapeutic exposure of anything here.
Still liver, still mostly toxicology doses.

**ARCHS4** is not a perturbation corpus and is deliberately not treated as
one. It is uniformly reprocessed RNA-seq across GEO, and its useful output
here is co-expression: the genes that move with a factor across thousands
of samples. That is a regulon, not a drug effect, so it is exposed as an
alternative input to the activity statistic rather than as a drug source.

What agreement is worth, and what it is not
-------------------------------------------
A compound that shifts a factor's regulon in LINCS and in DrugMatrix has
done so in a human cell line and in rat liver, on different platforms, in
different labs. The shared systematic errors are few, so agreement is
worth far more than a stronger result in either alone.

Disagreement is the harder case and usually means nothing. These corpora
differ in species, tissue, dose, duration and platform all at once, so a
compound active in one and not another may simply not have been tested
comparably. Absence in a corpus is only evidence when the compound was
actually profiled in it, which is why coverage is reported per corpus
rather than inferred.

None of this fixes the attribution problem. Every corpus here measures the
regulon, so every one of them is subject to the same confound: a factor
sharing targets with another cannot be separated from it, and a compound
that arrests the cell cycle moves every proliferation-linked regulon. More
corpora make the observation more robust and do not make it more specific.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

log = logging.getLogger(__name__)


@dataclass(frozen=True)
class Corpus:
    """One perturbation corpus and the caveats that travel with it."""
    key: str
    name: str
    species: str
    system: str                 # what was actually dosed
    paths: tuple                # (up, down) or (single,) for signed-in-name
    enrichr: tuple = ()         # Enrichr library names, for the download
    split: bool = True          # separate up/down files
    caveats: tuple = ()

    def available(self) -> bool:
        from pathlib import Path
        return all(any(Path(c).exists() for c in group)
                   for group in self.paths)

    def resolved(self) -> tuple:
        from pathlib import Path
        out = []
        for group in self.paths:
            hit = next((c for c in group if Path(c).exists()), None)
            if hit is None:
                return ()
            out.append(hit)
        return tuple(out)


def _p(*names):
    """Accept the file in data/ or beside the working directory."""
    return tuple(f"data/{n}" for n in names) + names


CORPORA: tuple = (
    Corpus(
        key="lincs", name="LINCS L1000 chemical perturbations",
        species="human", system="cancer cell lines, hours at micromolar",
        paths=(_p("lincs_chem_pert_up.gmt"), _p("lincs_chem_pert_down.gmt")),
        enrichr=("LINCS_L1000_Chem_Pert_up", "LINCS_L1000_Chem_Pert_down"),
        caveats=(
            "978 genes are measured and the rest are inferred, so a result "
            "on a non-landmark gene is an imputation",
            "transformed proliferating lines at doses above therapeutic",
        )),
    Corpus(
        key="creeds", name="CREEDS single drug perturbations",
        species="human and mouse", system="published GEO experiments",
        paths=(_p("creeds_drug_up.gmt"), _p("creeds_drug_down.gmt")),
        enrichr=("Drug_Perturbations_from_GEO_up",
                 "Drug_Perturbations_from_GEO_down"),
        caveats=(
            "crowd-extracted from heterogeneous published experiments, so "
            "platform, tissue and quality vary between signatures",
            "few replicates per signature; several are a single contrast",
            "mixes human and mouse studies, and the symbols are not "
            "harmonised here",
        )),
    Corpus(
        key="drugmatrix", name="DrugMatrix toxicogenomics",
        species="rat", system="tissue dosed in vivo, mostly liver",
        paths=(_p("drugmatrix_up.gmt"), _p("drugmatrix_down.gmt")),
        enrichr=("DrugMatrix",),
        caveats=(
            "RAT in vivo. Every symbol crosses a species boundary, and rat "
            "and human regulons are not interchangeable",
            "liver dominates the panel, so a factor whose role is "
            "elsewhere is poorly covered",
            "toxicology doses, chosen to produce an effect rather than to "
            "model a therapeutic exposure",
        )),
    Corpus(
        key="tggates", name="Open TG-GATEs",
        species="rat and human", system="rat in vivo and primary hepatocytes",
        paths=(_p("tggates_up.gmt"), _p("tggates_down.gmt")),
        enrichr=(),
        caveats=(
            "the human arm is primary hepatocytes, which is the closest "
            "thing here to a human therapeutic exposure; the rat arm is "
            "not",
            "liver only",
            "about 170 compounds, so most of the approved library is "
            "simply absent rather than negative",
        )),
)

BY_KEY = {c.key: c for c in CORPORA}

DOWNLOAD = """\
Each corpus is a pair of Enrichr gene-set libraries, in the same .gmt
format the LINCS files already use. Fetch the ones you want into data/:

  CREEDS
    curl -L -o data/creeds_drug_up.gmt "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Drug_Perturbations_from_GEO_up"
    curl -L -o data/creeds_drug_down.gmt "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Drug_Perturbations_from_GEO_down"

  DrugMatrix  (one library; up and down are distinguished in the term name)
    curl -L -o data/drugmatrix_raw.gmt "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=DrugMatrix"
    python -m promoter_tf_drugs.corpora --split data/drugmatrix_raw.gmt --as drugmatrix

  TG-GATEs is not on Enrichr. Open TG-GATEs publishes the raw arrays at
  https://dbarchive.biosciencedbc.jp/en/open-tggates/ ; a derived up/down
  .gmt placed at data/tggates_{up,down}.gmt is picked up automatically.

Nothing is required. A corpus that is absent is skipped and said to be
absent, which is different from a corpus that ran and found nothing.
"""


# ------------------------------------------------------------- splitting
def split_signed_gmt(path, out_prefix, *, up_markers=("-up", "_up", " up"),
                     down_markers=("-dn", "-down", "_dn", "_down",
                                   " down")) -> dict:
    """Split a single .gmt whose direction is in the term name.

    DrugMatrix ships as one library where each term ends in something like
    `-up` or `-dn`. The scan needs the two directions in separate files,
    the way the LINCS libraries already come.
    """
    from pathlib import Path

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{path} not found")
    up_lines, down_lines, skipped = [], [], 0

    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip():
                continue
            term = line.split("\t", 1)[0].strip().lower()
            if any(term.endswith(m) or m in term for m in down_markers):
                down_lines.append(line)
            elif any(term.endswith(m) or m in term for m in up_markers):
                up_lines.append(line)
            else:
                skipped += 1

    up_path = Path(f"{out_prefix}_up.gmt")
    down_path = Path(f"{out_prefix}_down.gmt")
    up_path.write_text("".join(up_lines), encoding="utf-8")
    down_path.write_text("".join(down_lines), encoding="utf-8")
    log.info("split %s: %d up, %d down, %d unlabelled and dropped",
             path.name, len(up_lines), len(down_lines), skipped)
    return {"up": str(up_path), "down": str(down_path),
            "n_up": len(up_lines), "n_down": len(down_lines),
            "n_unlabelled": skipped}


# ------------------------------------------------------------- the scan
@dataclass
class CorpusResult:
    corpus: str
    available: bool
    hits: list = field(default_factory=list)
    n_signatures: int = 0
    error: str | None = None
    caveats: tuple = ()

    def to_dict(self) -> dict:
        return dict(self.__dict__)


def scan(tf: str, regulon: dict, *, keys=None, **kw) -> dict:
    """Run the regulon statistic over every available corpus.

    Returns `{corpus key: CorpusResult}`. A corpus with no files present
    is reported as unavailable rather than as finding nothing, because
    those are different statements and only one of them is evidence.
    """
    from . import tfactivity as ta

    want = [BY_KEY[k] for k in keys] if keys else list(CORPORA)
    out: dict = {}
    for c in want:
        resolved = c.resolved()
        if not resolved:
            out[c.key] = CorpusResult(
                corpus=c.key, available=False,
                error=(f"{c.name} is not on disk. It is skipped, which is "
                       f"not the same as finding nothing in it."),
                caveats=c.caveats)
            continue
        up, down = resolved[0], resolved[-1]
        try:
            hits = ta.drugs_moving_factor(tf, regulon, up_path=up,
                                          down_path=down, **kw)
        except Exception as exc:                            # noqa: BLE001
            out[c.key] = CorpusResult(corpus=c.key, available=True,
                                      error=str(exc)[:160],
                                      caveats=c.caveats)
            continue
        out[c.key] = CorpusResult(
            corpus=c.key, available=True,
            hits=[h.to_dict() for h in hits], n_signatures=len(hits),
            caveats=c.caveats)
    return out


def agreement(results: dict, *, min_corpora: int = 2) -> list:
    """Compounds that move the factor the same way in several corpora.

    The payoff for running more than one. These corpora differ in species,
    tissue, platform and dose, so they share few systematic errors and a
    compound appearing in two of them is supported in a way no amount of
    internal replication within one corpus can match.

    Matching is on the salt-free name, because the corpora do not share a
    compound vocabulary.
    """
    from .congeners import salt_free_name

    seen: dict = {}
    for key, res in (results or {}).items():
        if not res.available if hasattr(res, "available") else not res.get(
                "available"):
            continue
        hits = res.hits if hasattr(res, "hits") else res.get("hits") or []
        for h in hits:
            stem = salt_free_name((h.get("drug") or "").lower())
            if not stem:
                continue
            rec = seen.setdefault(stem, {"drug": h.get("drug"),
                                         "corpora": {}, "directions": set()})
            rec["corpora"][key] = {
                "direction": h.get("direction"),
                "median_z": h.get("median_z"),
                "n_signatures": h.get("n_signatures"),
            }
            rec["directions"].add(h.get("direction"))

    out = []
    for stem, rec in seen.items():
        if len(rec["corpora"]) < min_corpora:
            continue
        consistent = len(rec["directions"]) == 1
        out.append({
            "drug": rec["drug"],
            "n_corpora": len(rec["corpora"]),
            "corpora": rec["corpora"],
            "consistent": consistent,
            "direction": (next(iter(rec["directions"])) if consistent
                          else "conflicting"),
            "note": (
                f"moves this factor's regulon in {len(rec['corpora'])} "
                f"independent corpora "
                f"({', '.join(sorted(rec['corpora']))}), which differ in "
                f"species, tissue, platform and dose. They share few "
                f"systematic errors, so agreement across them is worth "
                f"more than a stronger result in any one"
                if consistent else
                f"appears in {len(rec['corpora'])} corpora with "
                f"CONFLICTING directions. That usually means the corpora "
                f"are not comparable - different species, tissue or dose - "
                f"rather than that one is wrong"),
        })
    out.sort(key=lambda r: (-r["n_corpora"], not r["consistent"]))
    return out


def summarise(results: dict, agreed: list | None = None) -> str:
    present = [k for k, r in (results or {}).items()
               if (r.available if hasattr(r, "available")
                   else r.get("available"))]
    absent = [k for k in (results or {}) if k not in present]
    if not present:
        return ("No perturbation corpus is on disk beyond what is already "
                "loaded, so this scan had nothing independent to check "
                "against. " + DOWNLOAD.split("\n")[0])

    counts = []
    for k in present:
        r = results[k]
        n = r.n_signatures if hasattr(r, "n_signatures") else r.get(
            "n_signatures", 0)
        counts.append(f"{k} ({n})")
    bits = [f"Scanned {len(present)} corpus/corpora: {', '.join(counts)}."]
    if absent:
        bits.append(
            f"{len(absent)} absent and therefore skipped "
            f"({', '.join(sorted(absent))}) - absent is not the same as "
            f"finding nothing.")
    if agreed:
        consistent = [a for a in agreed if a["consistent"]]
        bits.append(
            f"**{len(consistent)} compound(s) agree across corpora**, which "
            f"is the result worth acting on: these differ in species, "
            f"tissue, platform and dose, so they share few systematic "
            f"errors.")
        conflicting = [a for a in agreed if not a["consistent"]]
        if conflicting:
            bits.append(
                f"{len(conflicting)} appear in more than one corpus with "
                f"conflicting directions, which usually reflects "
                f"incomparable systems rather than a wrong answer.")
    elif len(present) > 1:
        bits.append(
            "No compound moved the regulon consistently in more than one "
            "corpus. With this few overlapping compounds that is the "
            "expected outcome rather than a negative result.")
    return " ".join(bits)


def main(argv=None) -> int:
    """`python -m promoter_tf_drugs.corpora` - status and splitting."""
    import argparse

    p = argparse.ArgumentParser(
        prog="python -m promoter_tf_drugs.corpora",
        description="Which perturbation corpora are on disk, and how to "
                    "prepare the ones that need it.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=DOWNLOAD)
    p.add_argument("--split", metavar="GMT",
                   help="split a single .gmt whose direction is in the term "
                        "name into up and down files")
    p.add_argument("--as", dest="prefix", metavar="NAME",
                   help="output prefix for --split, e.g. drugmatrix -> "
                        "data/drugmatrix_{up,down}.gmt")
    args = p.parse_args(argv)

    if args.split:
        if not args.prefix:
            p.error("--split needs --as NAME")
        prefix = args.prefix
        if "/" not in prefix and "\\\\" not in prefix:
            prefix = f"data/{prefix}"
        res = split_signed_gmt(args.split, prefix)
        print(f"wrote {res['up']} ({res['n_up']} terms)")
        print(f"wrote {res['down']} ({res['n_down']} terms)")
        if res["n_unlabelled"]:
            print(f"{res['n_unlabelled']} term(s) carried no direction "
                  f"marker and were dropped")
        return 0

    print("Perturbation corpora")
    print("--------------------")
    for c in CORPORA:
        state = "present" if c.available() else "absent"
        print(f"  {c.key:12s} {state:8s} {c.species:18s} {c.system}")
    print()
    if not all(c.available() for c in CORPORA):
        print(DOWNLOAD)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
