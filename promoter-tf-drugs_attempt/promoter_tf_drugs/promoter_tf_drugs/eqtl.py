"""Trans-eQTL evidence: does variation in the factor actually move your gene?

Why this is the strongest link in the package
---------------------------------------------
Everything else here is inference. A motif says a protein *could* bind. A GO
term says a factor *is annotated* as an activator. A trans-eQTL says that
natural genetic variation near the transcription factor **does** change your
gene's expression, in real human tissue, with exposure randomised by meiosis.
That is causal evidence of the regulatory relationship, and nothing else in
this pipeline comes close to it.

The report used to end each factor with "look for trans-eQTLs of X affecting
your gene in eQTLGen" — advice, not an answer, repeated verbatim per factor.
The significant-results file is 2 MB. There is no reason not to just look.

What the lookup can and cannot conclude
---------------------------------------
**A hit is strong.** Variation near the factor associates with your gene's
expression across tens of thousands of people.

**Absence is weak, and this matters more than it sounds.** eQTLGen's
significant trans-eQTL set covers 6,297 genes. Most genes are simply not in
it, because trans effects are small and detecting them needs enormous power —
so "no trans-eQTL" overwhelmingly means *underpowered*, not *unregulated*.
CCR5, for one, has none at all, which says nothing about whether CCR5 is
trans-regulated.

**It is whole blood.** eQTLGen is a blood eQTL resource. A regulatory
relationship that operates in hepatocytes or neurons will not appear, however
real it is.

Source
------
  eQTLGen Consortium, trans-eQTLs FDR<0.05, n=31,684 whole blood
  https://www.eqtlgen.org/
"""

from __future__ import annotations

import gzip
import io
import logging
from dataclasses import dataclass, field
from pathlib import Path

log = logging.getLogger(__name__)

TRANS_URL = ("https://molgenis26.gcc.rug.nl/downloads/eqtlgen/trans-eqtl/"
             "2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-"
             "BonferroniAdded.txt.gz")

# How near a SNP must be to the factor's locus to count as "variation in this
# factor". Regulatory variants sit outside the gene body, so a window rather
# than an overlap; 1 Mb is the conventional cis window.
CIS_WINDOW = 1_000_000

_CACHE: dict = {"rows": None, "targets": None}


@dataclass
class TransHit:
    snp: str
    snp_chr: str
    snp_pos: int
    gene: str
    z: float
    p: float
    fdr: float
    n: int

    def to_dict(self) -> dict:
        return dict(self.__dict__)


@dataclass
class TransSupport:
    tf: str
    gene: str
    checked: bool = False
    supported: bool = False
    hits: list = field(default_factory=list)
    note: str = ""

    def to_dict(self) -> dict:
        return {"tf": self.tf, "gene": self.gene, "checked": self.checked,
                "supported": self.supported,
                "hits": [h.to_dict() for h in self.hits], "note": self.note}


def _fetch(client) -> bytes:
    """Download the table once, to a file beside the HTTP cache.

    It cannot go through `Client.request`: that caches by writing the payload
    as gzipped JSON *text* and decodes responses with `errors="replace"`,
    either of which destroys a gzip stream. A 2 MB bulk file belongs on disk
    anyway, the way `--pfm-file` and `--regulon-file` already work.
    """
    import urllib.request

    cache_dir = Path(getattr(client, "cache_dir", ".")) / "eqtlgen"
    cache_dir.mkdir(parents=True, exist_ok=True)
    local = cache_dir / "trans-eqtl-fdr0.05.txt.gz"

    if local.exists() and local.stat().st_size > 100_000:
        return local.read_bytes()

    if getattr(client, "offline", False):
        raise RuntimeError(
            f"offline and {local} not present; fetch it once with network "
            f"access, or pass a local copy")

    req = urllib.request.Request(
        TRANS_URL, headers={"User-Agent": "promoter-tf-drugs"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        raw = resp.read()
    if len(raw) < 100_000:
        raise RuntimeError(
            f"eQTLGen download was only {len(raw)} bytes — expected ~2 MB, "
            f"so this is probably an error page")
    tmp = local.with_suffix(".tmp")
    tmp.write_bytes(raw)
    tmp.replace(local)
    log.info("cached eQTLGen trans-eQTLs (%.1f MB) at %s",
             len(raw) / 1e6, local)
    return raw


def load_trans(client=None, path: str | None = None) -> list:
    """Parse the eQTLGen significant trans-eQTL table (cached in memory).

    2 MB gzipped, ~59,800 rows. `path` reads a local copy; otherwise the file
    is downloaded once into the client's cache directory, so a second run
    costs nothing.
    """
    if _CACHE["rows"] is not None:
        return _CACHE["rows"]

    if path:
        raw = Path(path).read_bytes()
    elif client is not None:
        raw = _fetch(client)
    else:
        raise ValueError("need a client or a local path")

    if isinstance(raw, str):
        raw = raw.encode("utf-8", "replace")
    text = gzip.open(io.BytesIO(raw), "rt", encoding="utf-8",
                     errors="replace") if raw[:2] == b"\x1f\x8b" else \
        io.StringIO(raw.decode("utf-8", "replace"))

    rows, targets = [], set()
    with text as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(hdr):
                continue
            try:
                rows.append(TransHit(
                    snp=f[idx["SNP"]], snp_chr=f[idx["SNPChr"]],
                    snp_pos=int(f[idx["SNPPos"]]),
                    gene=f[idx["GeneSymbol"]].upper(),
                    z=float(f[idx["Zscore"]]), p=float(f[idx["Pvalue"]]),
                    fdr=float(f[idx["FDR"]]), n=int(f[idx["NrSamples"]])))
            except (ValueError, KeyError, IndexError):
                continue
            targets.add(rows[-1].gene)

    _CACHE["rows"], _CACHE["targets"] = rows, targets
    log.info("eQTLGen: %d significant trans-eQTLs over %d target genes",
             len(rows), len(targets))
    return rows


def gene_has_any_trans_eqtl(gene: str, client=None, path=None) -> bool:
    load_trans(client, path)
    return gene.upper() in (_CACHE["targets"] or set())


def trans_support(tf: str, gene: str, *, tf_locus: tuple | None = None,
                  client=None, path=None) -> TransSupport:
    """Does variation near `tf` associate with `gene` expression in trans?

    `tf_locus` is (chrom, start, end) for the factor. Without it the lookup
    can still answer the prior question — whether the target gene has any
    trans-eQTL at all — which for most genes is the whole answer.
    """
    res = TransSupport(tf=tf.upper(), gene=gene.upper())
    try:
        rows = load_trans(client, path)
    except Exception as exc:                                # noqa: BLE001
        res.note = (f"eQTLGen lookup unavailable ({str(exc)[:80]}); the "
                    f"trans-eQTL question was NOT answered")
        return res

    res.checked = True
    for_gene = [h for h in rows if h.gene == res.gene]
    if not for_gene:
        res.note = (
            f"{res.gene} has NO significant trans-eQTL in eQTLGen "
            f"(n=31,684 whole blood, FDR<0.05) — it is not among the 6,297 "
            f"genes with one. So this test returns nothing for {res.tf} or "
            f"for any other factor, and there is no point checking them one "
            f"by one. Absence here is weak evidence: trans effects are small "
            f"and detecting them needs enormous power, so this overwhelmingly "
            f"means underpowered rather than unregulated — and eQTLGen is "
            f"blood, so a relationship operating elsewhere would not show.")
        return res

    if not tf_locus:
        res.note = (
            f"{res.gene} has {len(for_gene)} significant trans-eQTL(s) in "
            f"eQTLGen, but the locus of {res.tf} was not supplied, so whether "
            f"any of them sits near {res.tf} was not tested.")
        return res

    chrom, start, end = tf_locus
    chrom = str(chrom).replace("chr", "")
    near = [h for h in for_gene
            if h.snp_chr == chrom
            and (start - CIS_WINDOW) <= h.snp_pos <= (end + CIS_WINDOW)]
    res.hits = sorted(near, key=lambda h: h.p)[:10]
    res.supported = bool(near)
    if near:
        best = res.hits[0]
        res.note = (
            f"**Genetic support.** {len(near)} variant(s) within 1 Mb of "
            f"{res.tf} associate with {res.gene} expression in trans "
            f"(best {best.snp}, Z={best.z:.1f}, p={best.p:.2g}, "
            f"n={best.n}). Direction: variation raising {res.tf} activity "
            f"tracks {'higher' if best.z > 0 else 'lower'} {res.gene}. This "
            f"is in-vivo causal evidence for the regulatory link — the "
            f"strongest available anywhere in this report.")
    else:
        res.note = (
            f"{res.gene} does have {len(for_gene)} trans-eQTL(s), but none "
            f"within 1 Mb of {res.tf}. Weak evidence against this particular "
            f"link, not against trans-regulation of {res.gene}.")
    return res


def summarise(gene: str, supports: list) -> str:
    """One statement for the whole factor list, instead of N identical ones."""
    if not supports:
        return ""
    if not any(s.checked for s in supports):
        return (f"Trans-eQTL evidence was not retrieved, so the strongest "
                f"available test of these regulatory links was not run.")
    if all(s.checked and not s.hits and "NO significant trans-eQTL" in s.note
           for s in supports):
        return (f"**No factor can be supported this way for {gene.upper()}.** "
                f"{gene.upper()} has no significant trans-eQTL in eQTLGen at "
                f"all, so the test returns nothing for every factor here. "
                f"That is a property of the dataset's power and of the tissue "
                f"(whole blood), not evidence that {gene.upper()} is "
                f"unregulated.")
    hit = [s for s in supports if s.supported]
    if hit:
        return (f"**{len(hit)} of {len(supports)} factors have human genetic "
                f"support for regulating {gene.upper()}**: "
                + ", ".join(s.tf for s in hit)
                + ". This is causal evidence in real tissue, and it "
                  "outranks every motif and annotation in this report.")
    return (f"{gene.upper()} has trans-eQTLs, but none near any of the "
            f"factors tested. Weak evidence against these particular links.")
