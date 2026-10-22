"""Fetch every reference file the analysis reads from data/, in one go.

    python scripts/setup_data.py              # show what is missing, ask, fetch
    python scripts/setup_data.py --yes        # fetch without asking
    python scripts/setup_data.py --list       # show the inventory and stop
    python scripts/setup_data.py --small      # skip the files over about 500 MB
    python scripts/setup_data.py --only remap,collectri
    python scripts/setup_data.py --with-optional   # also the extra LINCS corpora

`scripts\\setup_data.bat` runs this with the checkout's Python.

These files are git-ignored: large, versioned upstream, and re-fetchable. Each
item is skipped when a file the loaders accept is already in data/. A download
goes to a `.part` file and is kept only if it passes a check on its content (a
gzip that opens, a zip that holds the expected member, a table with the expected
header, a published checksum), so an error page or a dropped connection never
leaves something that looks like data.

Nothing here validates a drug or a method. It only fetches inputs, from the
sources listed in `docs/datasets.md`. Release names change, so a URL that stops
working is reported as a failure with the source to check by hand.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.request
import zipfile
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

AGENT = "Mozilla/5.0 generoutes-setup"
FIGSHARE = "https://api.figshare.com/v2/articles/{}/files?page_size=1000"
DEPMAP_ARTICLE = 25880521      # DepMap Public 24Q2
PRISM_ARTICLE = 25917643       # Repurposing Public 24Q2
# The ids the launcher has always used, kept as a fallback if the listing
# cannot be read. Their checksums are what the launcher checks.
FIGSHARE_FALLBACK = {
    "CRISPRGeneEffect.csv": (46489063, "d10d73692672380b3792231549e97d5a"),
    "Model.csv": (46489732, None),
    "CRISPRInferredCommonEssentials.csv": (46489597, None),
    "Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv":
        (46630984, "042b15d3ebcc138a80a571b9286451ae"),
    "Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv":
        (46630981, None),
    "OmicsExpressionProteinCodingGenesTPMLogp1.csv": (46490878, None),
}
ENRICHR = ("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text"
           "&libraryName={}")
OMNIPATH_COLLECTRI = (
    "https://omnipathdb.org/interactions?datasets=collectri&genesymbols=yes"
    "&organisms=9606&fields=sources,references&license=academic&format=tsv")
GENCODE_DIR = ("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
               "latest_release/")
BINDINGDB_ZIP = ("https://www.bindingdb.org/rwd/bind/downloads/"
                 "BindingDB_All_202610_tsv.zip")


class SetupError(Exception):
    """One item could not be fetched or did not pass its check."""


# ------------------------------------------------------------- network
def _open(url: str, timeout: int = 60):
    req = urllib.request.Request(url, headers={"User-Agent": AGENT})
    return urllib.request.urlopen(req, timeout=timeout)


def get_text(url: str) -> str:
    with _open(url) as r:
        return r.read().decode("utf-8", "replace")


def get_json(url: str):
    return json.loads(get_text(url))


def download(url: str, dest, *, check=None, retries: int = 3) -> Path:
    """Fetch `url` to `dest` through a `.part` file, kept only if `check` passes.

    `check(path)` returns None when the file is good and a reason when it is
    not. Connection errors are retried; a failed check is not, since the same
    URL would answer the same way.
    """
    dest = Path(dest)
    part = dest.with_name(dest.name + ".part")
    dest.parent.mkdir(parents=True, exist_ok=True)
    last = None
    for attempt in range(1, retries + 1):
        try:
            with _open(url, timeout=120) as r, open(part, "wb") as fh:
                total = int(r.headers.get("Content-Length") or 0)
                done, mark = 0, 0
                while True:
                    block = r.read(1 << 20)
                    if not block:
                        break
                    fh.write(block)
                    done += len(block)
                    if total and done * 10 // total > mark:
                        mark = done * 10 // total
                        print(f"    {mark * 10}% of {total / 1e6:.0f} MB",
                              flush=True)
            if total and done != total:
                raise OSError(f"got {done} of {total} bytes")
            last = None
            break
        except OSError as exc:
            last = exc
            if part.exists():
                part.unlink()
            if attempt < retries:
                time.sleep(2 * attempt)
    if last is not None:
        raise SetupError(f"could not fetch {url}: {last}")
    reason = check(part) if check else None
    if reason:
        part.unlink(missing_ok=True)
        raise SetupError(f"{dest.name}: {reason}")
    os.replace(part, dest)
    return dest


# ------------------------------------------------------------- checks
def contains(token: str):
    def check(path):
        with open(path, "rb") as fh:
            head = fh.read(2_000_000)
        return None if token.encode() in head else (
            f"does not contain {token!r}, so it is not the expected file")
    return check


def first_line_has(*tokens: str):
    def check(path):
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            line = fh.readline()
        return None if all(t in line for t in tokens) else (
            f"header lacks {tokens}, so it is not the expected table")
    return check


def gzip_opens(path):
    try:
        with gzip.open(path, "rb") as fh:
            fh.read(1024)
    except (OSError, EOFError):
        return "is not a readable gzip file"
    return None


def zip_holds(suffix: str):
    def check(path):
        try:
            names = zipfile.ZipFile(path).namelist()
        except zipfile.BadZipFile:
            return "is not a zip file"
        return None if any(n.lower().endswith(suffix) for n in names) else (
            f"holds no {suffix} member")
    return check


def md5_is(want: str | None):
    def check(path):
        if not want:
            return None
        h = hashlib.md5()
        with open(path, "rb") as fh:
            for block in iter(lambda: fh.read(1 << 20), b""):
                h.update(block)
        return None if h.hexdigest() == want.lower() else (
            f"checksum {h.hexdigest()} does not match the published {want}")
    return check


def gmt_ok(min_lines: int):
    def check(path):
        import check_gmt
        ok, n, good = check_gmt.check(path, min_lines)
        return None if ok else (
            f"is not a gene set library ({good} of {n} lines well formed)")
    return check


def all_of(*checks):
    def check(path):
        for c in checks:
            reason = c(path)
            if reason:
                return reason
        return None
    return check


# ------------------------------------------------------------- figshare
def figshare_file(article: int, name: str) -> tuple:
    """`(download url, md5 or None)` for a named file of a figshare article."""
    try:
        for f in get_json(FIGSHARE.format(article)):
            if f.get("name") == name:
                return f["download_url"], (f.get("computed_md5") or None)
    except (OSError, ValueError, KeyError):
        pass
    if name in FIGSHARE_FALLBACK:
        fid, md5 = FIGSHARE_FALLBACK[name]
        return f"https://ndownloader.figshare.com/files/{fid}", md5
    raise SetupError(f"{name} is not in figshare article {article}")


def fetch_figshare(data: Path, article: int, name: str, *, check=None):
    url, md5 = figshare_file(article, name)
    download(url, data / name,
             check=all_of(md5_is(md5), check) if check else md5_is(md5))


# ------------------------------------------------------------- items
@dataclass
class Item:
    key: str
    title: str
    size: str                      # as fetched, approximate
    group: str                     # "core", "big" or "optional"
    have: tuple                    # per file, the names that satisfy it
    fetch: object                  # (data_dir) -> None
    loses: str = ""
    needs: tuple = ()              # keys that must be present first
    manual: str = ""               # where to get it by hand if this fails


def _reg(key: str) -> tuple:
    """Every filename the loaders accept for a registered dataset."""
    from generoutes import datasets
    for ds in datasets.registry():
        if ds.key == key:
            return tuple(ds.paths)
    raise KeyError(key)


def present(item: Item, data: Path) -> bool:
    for alts in item.have:
        if not any(_exists(data, a) for a in alts):
            return False
    return True


def _exists(data: Path, name: str) -> bool:
    p = Path(name)
    if p.parts and p.parts[0] == "data":
        p = Path(*p.parts[1:])
    f = data / p
    return f.is_file() and f.stat().st_size > 0


def _fetch_simple(url, name, check):
    def run(data):
        download(url, data / name, check=check)
    return run


def _fetch_collectri(data: Path):
    """CollecTRI as source, target, weight. decoupler's table when installed,
    otherwise OmniPath's CollecTRI dataset converted to the same columns."""
    dest = data / "collectri.tsv"
    try:
        import decoupler as dc
        getter = (getattr(dc, "get_collectri", None)
                  or getattr(getattr(dc, "op", None), "collectri", None))
        if getter is not None:
            df = getter(organism="human")
            if len(df) > 20000 and {"source", "target", "weight"} <= set(
                    df.columns):
                df.to_csv(dest, sep="\t", index=False)
                print("    written from decoupler")
                return
    except Exception as exc:                                  # noqa: BLE001
        print(f"    decoupler not used ({str(exc)[:60]})")
    raw = data / "collectri.omnipath.tsv"
    download(OMNIPATH_COLLECTRI, raw,
             check=first_line_has("source_genesymbol", "is_stimulation"))
    try:
        rows = convert_omnipath(raw.read_text(encoding="utf-8"))
    finally:
        raw.unlink(missing_ok=True)
    if len(rows) < 20000:
        raise SetupError(f"OmniPath gave only {len(rows)} signed edges")
    with open(dest, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["source", "target", "weight"])
        w.writerows(rows)
    print(f"    written from OmniPath ({len(rows)} signed edges). This is "
          f"OmniPath's current CollecTRI dataset, which can differ from "
          f"decoupler's table; install decoupler to match it.")


def convert_omnipath(text: str) -> list:
    """`[(source, target, +1 or -1)]` from OmniPath interaction rows.

    An edge that both stimulates and inhibits, or does neither, has no single
    sign and is left out, as are edges with no gene symbol.
    """
    out = []
    for r in csv.DictReader(io.StringIO(text), delimiter="\t"):
        s, t = r.get("source_genesymbol"), r.get("target_genesymbol")
        if not s or not t:
            continue
        stim = r.get("is_stimulation") == "True"
        inh = r.get("is_inhibition") == "True"
        if stim != inh:
            out.append((s, t, 1 if stim else -1))
    return out


def _fetch_gencode(data: Path):
    html = get_text(GENCODE_DIR)
    names = sorted(set(re.findall(
        r"gencode\.v(\d+)\.basic\.annotation\.gtf\.gz", html)), key=int)
    if not names:
        raise SetupError("the GENCODE listing names no basic annotation GTF")
    name = f"gencode.v{names[-1]}.basic.annotation.gtf.gz"
    download(GENCODE_DIR + name, data / "gencode.basic.gtf.gz",
             check=gzip_opens)
    print(f"    GENCODE release {names[-1]}")


def _fetch_hpa(data: Path):
    z = data / "rna_tissue_consensus.zip"
    download("https://www.proteinatlas.org/download/tsv/"
             "rna_tissue_consensus.tsv.zip", z,
             check=zip_holds("rna_tissue_consensus.tsv"))
    try:
        with zipfile.ZipFile(z) as zf:
            member = next(n for n in zf.namelist()
                          if n.lower().endswith("rna_tissue_consensus.tsv"))
            with zf.open(member) as src, open(
                    data / "rna_tissue_consensus.tsv", "wb") as out:
                shutil.copyfileobj(src, out)
    finally:
        z.unlink(missing_ok=True)


def _fetch_lincs(data: Path):
    for lib, name in (("LINCS_L1000_Chem_Pert_up", "lincs_chem_pert_up.gmt"),
                      ("LINCS_L1000_Chem_Pert_down",
                       "lincs_chem_pert_down.gmt")):
        if not _exists(data, name):
            download(ENRICHR.format(lib), data / name, check=gmt_ok(1000))
    if not _exists(data, "lincs_geneinfo.txt"):
        download("https://s3.amazonaws.com/macchiato.clue.io/builds/"
                 "LINCS2020/geneinfo_beta.txt", data / "lincs_geneinfo.txt",
                 check=first_line_has("gene_symbol", "feature_space"))


def _fetch_knockout(data: Path):
    for lib, name, n in (
            ("LINCS_L1000_CRISPR_KO_Consensus_Sigs",
             "knockout_lincs_crispr_ko.gmt", 1000),
            ("TF_Perturbations_Followed_by_Expression",
             "knockout_tf_perturbations.gmt", 500)):
        if not _exists(data, name):
            download(ENRICHR.format(lib), data / name, check=gmt_ok(n))


def _fetch_orangebook(data: Path):
    download("https://www.fda.gov/media/76860/download?attachment",
             data / "EOBZIP.zip", check=zip_holds("products.txt"))


def _fetch_gwas(data: Path):
    import unpack_gwas
    z = data / "gwas_catalog_associations.zip"
    download("https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
             "gwas-catalog-associations-full.zip", z,
             check=zip_holds(".tsv"))
    try:
        unpack_gwas.unpack(z, data / "gwas_catalog_associations.tsv")
    except (OSError, ValueError, zipfile.BadZipFile) as exc:
        raise SetupError(f"the GWAS zip could not be unpacked: {exc}")
    finally:
        z.unlink(missing_ok=True)


def _fetch_depmap(data: Path):
    """Each of the three files only when no accepted name for it is present,
    so a 420 MB effect table already on disk is not fetched again."""
    for key, name in (("depmap_effect", "CRISPRGeneEffect.csv"),
                      ("depmap_model", "Model.csv"),
                      ("depmap_essential",
                       "CRISPRInferredCommonEssentials.csv")):
        if not any(_exists(data, a) for a in _reg(key)):
            fetch_figshare(data, DEPMAP_ARTICLE, name)


def _fetch_prism(data: Path):
    for name in ("Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv",
                 "Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv"):
        if not _exists(data, name):
            fetch_figshare(data, PRISM_ARTICLE, name)


def _fetch_ccle(data: Path):
    fetch_figshare(data, DEPMAP_ARTICLE,
                   "OmicsExpressionProteinCodingGenesTPMLogp1.csv")


def _fetch_remap(data: Path):
    download("https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/"
             "remap2022_nr_macs2_hg38_v1_0.bed.gz", data / "remap2022.bed.gz",
             check=gzip_opens)


def _run(cmd: list):
    print("    " + " ".join(str(c) for c in cmd))
    done = subprocess.run([str(c) for c in cmd], cwd=ROOT)
    if done.returncode:
        raise SetupError(f"{cmd[2] if len(cmd) > 2 else cmd[0]} failed")


def _build_chip_regulon(data: Path):
    _run([sys.executable, "-m", "generoutes.chipregulon",
          "--bed", data / "remap2022.bed.gz",
          "--gtf", data / "gencode.basic.gtf.gz",
          "--out", data / "chip_regulon.tsv"])


def build_library(data: Path, client=None):
    """data/approved_library.tsv from ChEMBL: approved small molecules with a
    structure, with the withdrawn flag and oral route when ChEMBL has them."""
    from generoutes.http import Client
    client = client or Client()
    chembl = "https://www.ebi.ac.uk/chembl/api/data"
    rows, offset = [], 0
    while True:
        res = client.get(f"{chembl}/molecule", {
            "max_phase": 4, "molecule_type": "Small molecule",
            "limit": 200, "offset": offset, "format": "json"})
        mols = res.get("molecules", [])
        for m in mols:
            smi = (m.get("molecule_structures") or {}).get("canonical_smiles")
            if not smi:
                continue
            rows.append([(m.get("pref_name") or m["molecule_chembl_id"]),
                         smi, m["molecule_chembl_id"],
                         int(float(m.get("max_phase") or 0)),
                         1 if m.get("withdrawn_flag") else "",
                         1 if m.get("oral") else ""])
        offset += 200
        if not mols or not (res.get("page_meta") or {}).get("next"):
            break
    if len(rows) < 2500:
        raise SetupError(f"ChEMBL returned only {len(rows)} approved "
                         f"molecules with a structure; expected about 3,300")
    rows.sort(key=lambda r: str(r[0]).upper())
    with open(data / "approved_library.tsv", "w", encoding="utf-8",
              newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["name", "smiles", "chembl_id", "max_phase", "withdrawn",
                    "oral"])
        w.writerows(rows)
    print(f"    {len(rows)} approved molecules")


def _fetch_library(data: Path):
    build_library(data)


def _fetch_bindingdb(data: Path):
    z = data / "BindingDB_All.zip"
    raw = data / "BindingDB_All.tsv"
    download(BINDINGDB_ZIP, z, check=zip_holds("bindingdb_all.tsv"))
    try:
        with zipfile.ZipFile(z) as zf:
            member = next(n for n in zf.namelist()
                          if n.lower().endswith("bindingdb_all.tsv"))
            print("    unpacking about 9 GB; it is distilled and deleted next")
            with zf.open(member) as src, open(raw, "wb") as out:
                shutil.copyfileobj(src, out, 1 << 22)
        _run([sys.executable, "-m", "generoutes.localpotency",
              "--distil", raw])
    finally:
        z.unlink(missing_ok=True)
        raw.unlink(missing_ok=True)


def _fetch_corpora(data: Path):
    for lib, name in (("Drug_Perturbations_from_GEO_up", "creeds_drug_up.gmt"),
                      ("Drug_Perturbations_from_GEO_down",
                       "creeds_drug_down.gmt")):
        if not _exists(data, name):
            download(ENRICHR.format(lib), data / name, check=gmt_ok(100))
    if not (_exists(data, "drugmatrix_up.gmt")
            and _exists(data, "drugmatrix_down.gmt")):
        raw = data / "drugmatrix_raw.gmt"
        download(ENRICHR.format("DrugMatrix"), raw, check=gmt_ok(100))
        _run([sys.executable, "-m", "generoutes.corpora",
              "--split", raw, "--as", "drugmatrix"])


def items() -> list:
    """Every file group, in the order a fetch should run."""
    return [
        Item("jaspar", "JASPAR 2024 motifs", "0.3 MB", "core",
             (_reg("motifs"),),
             _fetch_simple("https://jaspar.elixir.no/download/data/2024/CORE/"
                           "JASPAR2024_CORE_vertebrates_non-redundant_pfms_"
                           "jaspar.txt", "jaspar2024.txt", contains(">MA")),
             "the scan runs on a handful of bundled matrices"),
        Item("trrust", "TRRUST regulons", "0.3 MB", "core",
             (("data/trrust.tsv",),),
             _fetch_simple("https://www.grnpedia.org/trrust/data/"
                           "trrust_rawdata.human.tsv", "trrust.tsv",
                           contains("\t")),
             "regulon breadth falls back to bundled tiers"),
        Item("collectri", "CollecTRI signed regulons", "2.5 MB", "core",
             (("data/collectri.tsv", "data/collectri.csv"),),
             _fetch_collectri,
             "--regulator-route and --tf-activity have no signed table"),
        Item("clingen", "ClinGen dosage curation", "0.3 MB", "core",
             (_reg("clingen"),),
             _fetch_simple("https://ftp.clinicalgenome.org/"
                           "ClinGen_gene_curation_list_GRCh38.tsv",
                           "ClinGen_gene_curation_list_GRCh38.tsv",
                           contains("Haploinsufficiency Score")),
             "dosage rests on population constraint alone"),
        Item("gnomad", "gnomAD v4.1 constraint", "91 MB", "core",
             (_reg("gnomad"),),
             _fetch_simple("https://storage.googleapis.com/gcp-public-data--"
                           "gnomad/release/4.1/constraint/"
                           "gnomad.v4.1.constraint_metrics.tsv",
                           "gnomad.v4.1.constraint_metrics.tsv",
                           contains("lof.oe_ci.upper")),
             "the report says dosage was not assessed"),
        Item("orangebook", "FDA Orange Book", "1 MB", "core",
             (_reg("orangebook"),), _fetch_orangebook,
             "the FDA only filter cannot tell approved from not, and route "
             "of administration is unknown"),
        Item("gencode", "GENCODE basic annotation", "35 MB", "core",
             (_reg("gtf"),), _fetch_gencode,
             "gene lookup and the neighbour check need the network"),
        Item("hpa", "Human Protein Atlas tissue RNA", "5 MB zipped", "core",
             (_reg("expression"),), _fetch_hpa,
             "--cell-type-check cannot run and the tissue list is empty"),
        Item("hgnc", "HGNC complete set", "17 MB", "core",
             (_reg("hgnc"),),
             _fetch_simple("https://storage.googleapis.com/public-download-"
                           "files/hgnc/tsv/tsv/hgnc_complete_set.txt",
                           "hgnc_complete_set.txt",
                           first_line_has("hgnc_id", "symbol")),
             "symbol resolution needs Ensembl"),
        Item("goa", "GO annotations, human", "15 MB", "core",
             (_reg("goa"),),
             _fetch_simple("https://current.geneontology.org/annotations/"
                           "goa_human.gaf.gz", "goa_human.gaf.gz",
                           gzip_opens),
             "a factor's activator or repressor role is less often known"),
        Item("lincs", "LINCS signatures and gene info", "31 MB", "core",
             (("data/lincs_chem_pert_up.gmt",),
              ("data/lincs_chem_pert_down.gmt",),
              _reg("lincs_geneinfo")), _fetch_lincs,
             "no measured drug to expression link"),
        Item("knockout", "Enrichr knockout libraries", "19 MB", "core",
             (("data/knockout_lincs_crispr_ko.gmt",),
              ("data/knockout_tf_perturbations.gmt",)), _fetch_knockout,
             "responsiveness is not measured for any factor"),
        Item("library", "Approved drug library (ChEMBL)", "0.3 MB", "core",
             (_reg("library"),), _fetch_library,
             "the docking screen and the upstream drug index need ChEMBL "
             "queries one at a time"),
        Item("depmap", "DepMap 24Q2 CRISPR effect, models, essentials",
             "420 MB", "core",
             (_reg("depmap_effect"), _reg("depmap_model"),
              _reg("depmap_essential")), _fetch_depmap,
             "what a cell loses without a factor is not measured"),
        Item("prism", "PRISM 24Q2 repurposing matrix", "73 MB", "core",
             (_reg("prism"), _reg("prism_compounds")), _fetch_prism,
             "the measured phenocopy route cannot run"),
        Item("gwas", "GWAS Catalog associations", "70 MB zipped, 600 MB",
             "core", (_reg("gwas"),), _fetch_gwas,
             "genetics rests on Open Targets alone"),
        Item("remap", "ReMap 2022 ChIP-seq peaks", "1.4 GB", "big",
             (_reg("chip"),), _fetch_remap,
             "factor ranking rests on sequence alone"),
        Item("chip_regulon", "ChIP-derived regulon (built)", "27 MB built",
             "core", (_reg("regulon_chip"),), _build_chip_regulon,
             "factors curation never reached have no target count",
             needs=("remap", "gencode")),
        Item("bindingdb", "BindingDB measured affinities", "600 MB, 9 GB "
             "while unpacking", "big", (_reg("potency"),), _fetch_bindingdb,
             "the potency axis has nothing to read"),
        Item("ccle", "DepMap expression (CCLE)", "about 500 MB", "big",
             (_reg("ccle_expression"),), _fetch_ccle,
             "the expression half of the sensitivity scan is skipped"),
        Item("corpora", "CREEDS and DrugMatrix signatures", "about 60 MB",
             "optional",
             (("data/creeds_drug_up.gmt",), ("data/creeds_drug_down.gmt",),
              ("data/drugmatrix_up.gmt",), ("data/drugmatrix_down.gmt",)),
             _fetch_corpora, "agreement across corpora is not checked"),
    ]


# ------------------------------------------------------------- planning
def plan(data: Path, *, only=None, skip=(), small=False,
         with_optional=False) -> tuple:
    """`(to_fetch, present, left_out)` as lists of Items."""
    chosen, have, left = [], [], []
    for it in items():
        if only and it.key not in only:
            continue
        if it.key in skip:
            left.append((it, "skipped by --skip"))
        elif present(it, data):
            have.append(it)
        elif it.group == "optional" and not (with_optional or
                                             (only and it.key in only)):
            left.append((it, "optional; add --with-optional"))
        elif it.group == "big" and small and not (only and it.key in only):
            left.append((it, "over 500 MB; left out by --small"))
        else:
            chosen.append(it)
    return chosen, have, left


def fetch_all(data: Path, chosen: list) -> dict:
    """Run each fetch. `{key: None or the reason it failed}`."""
    data.mkdir(parents=True, exist_ok=True)
    result: dict = {}
    by_key = {it.key: it for it in items()}
    for it in chosen:
        blocked = [k for k in it.needs
                   if not present(by_key[k], data)]
        if blocked:
            result[it.key] = "needs " + ", ".join(blocked) + " first"
            print(f"[{it.key}] skipped: {result[it.key]}")
            continue
        print(f"[{it.key}] {it.title} ({it.size})", flush=True)
        try:
            it.fetch(data)
            if not present(it, data):
                raise SetupError("finished, but the expected file is "
                                 "not in data/")
            result[it.key] = None
            print("    done")
        except SetupError as exc:
            result[it.key] = str(exc)
            print(f"    FAILED: {exc}")
        except Exception as exc:                              # noqa: BLE001
            result[it.key] = f"{type(exc).__name__}: {str(exc)[:120]}"
            print(f"    FAILED: {result[it.key]}")
    return result


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--data-dir", default=str(ROOT / "data"))
    p.add_argument("--yes", "-y", action="store_true",
                   help="fetch without asking")
    p.add_argument("--list", action="store_true",
                   help="show what is present and missing, then stop")
    p.add_argument("--small", action="store_true",
                   help="skip the files over about 500 MB")
    p.add_argument("--with-optional", action="store_true",
                   help="also fetch the extra LINCS corpora")
    p.add_argument("--only", default="", help="comma separated item keys")
    p.add_argument("--skip", default="", help="comma separated item keys")
    args = p.parse_args(argv)

    known = {it.key for it in items()}
    only = {k for k in args.only.split(",") if k}
    skip = {k for k in args.skip.split(",") if k}
    bad = (only | skip) - known
    if bad:
        p.error(f"unknown item(s) {sorted(bad)}; known: {sorted(known)}")
    data = Path(args.data_dir)
    chosen, have, left = plan(data, only=only or None, skip=skip,
                              small=args.small,
                              with_optional=args.with_optional)

    print(f"Reference files in {data}")
    for it in have:
        print(f"  have     {it.key:13} {it.title}")
    for it in chosen:
        print(f"  MISSING  {it.key:13} {it.title} ({it.size})")
    for it, why in left:
        print(f"  left out {it.key:13} {it.title} - {why}")
    if args.list or not chosen:
        if not chosen:
            print("Nothing to fetch.")
        return 0
    if not args.yes:
        if not sys.stdin.isatty():
            print("Not a terminal, so nothing was fetched. Pass --yes.")
            return 1
        ans = input(f"Fetch {len(chosen)} item(s)? Y/n: ").strip().lower()
        if ans == "n":
            return 0
    result = fetch_all(data, chosen)
    failed = {k: v for k, v in result.items() if v}
    print()
    print(f"Fetched {len(result) - len(failed)} of {len(result)} item(s).")
    by_key = {it.key: it for it in items()}
    for k, why in failed.items():
        print(f"  {k}: {why}")
        print(f"      without it: {by_key[k].loses}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
