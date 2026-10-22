"""The setup script that fetches the git-ignored reference files.

Nothing here touches the network. `download`, `get_text` and `get_json` are
replaced, so these tests check the decisions the script makes: what counts as
present, what it refuses to keep, what it skips, and that each fetch leaves the
file the loaders look for.
"""

from __future__ import annotations

import gzip
import importlib.util
import io
import re
import shlex
import zipfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load(name):
    path = ROOT / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"_{name}", path)
    mod = importlib.util.module_from_spec(spec)
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    sys.modules[spec.name] = mod        # dataclasses looks the module up
    spec.loader.exec_module(mod)
    return mod


sd = _load("setup_data")


class _Resp:
    def __init__(self, body: bytes, length=True):
        self._b = io.BytesIO(body)
        self.headers = {"Content-Length": str(len(body))} if length else {}

    def read(self, n=-1):
        return self._b.read(n)

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


# --------------------------------------------------------------- download
def test_a_download_is_kept_only_when_it_passes_its_check(tmp_path, monkeypatch):
    monkeypatch.setattr(sd, "_open", lambda url, timeout=60: _Resp(b"<html>"))
    with pytest.raises(sd.SetupError, match="not the expected file"):
        sd.download("u", tmp_path / "m.txt", check=sd.contains(">MA"))
    assert not (tmp_path / "m.txt").exists()
    assert not (tmp_path / "m.txt.part").exists()

    monkeypatch.setattr(sd, "_open",
                        lambda url, timeout=60: _Resp(b">MA0001 X\nA [1]"))
    got = sd.download("u", tmp_path / "m.txt", check=sd.contains(">MA"))
    assert got.read_bytes().startswith(b">MA")
    assert not (tmp_path / "m.txt.part").exists()


def test_a_dropped_connection_is_retried_then_reported(tmp_path, monkeypatch):
    calls = []

    def flaky(url, timeout=60):
        calls.append(url)
        raise OSError("reset")
    monkeypatch.setattr(sd, "_open", flaky)
    monkeypatch.setattr(sd.time, "sleep", lambda s: None)
    with pytest.raises(sd.SetupError, match="could not fetch"):
        sd.download("http://x", tmp_path / "f", retries=3)
    assert len(calls) == 3 and not (tmp_path / "f").exists()


def test_a_short_read_is_not_kept(tmp_path, monkeypatch):
    class Short(_Resp):
        def __init__(self):
            super().__init__(b"abc")
            self.headers = {"Content-Length": "10"}
    monkeypatch.setattr(sd, "_open", lambda url, timeout=60: Short())
    monkeypatch.setattr(sd.time, "sleep", lambda s: None)
    with pytest.raises(sd.SetupError):
        sd.download("u", tmp_path / "f", retries=1)
    assert not (tmp_path / "f").exists()


# ----------------------------------------------------------------- checks
def test_the_content_checks_refuse_what_is_not_the_file(tmp_path):
    gz = tmp_path / "a.gz"
    gz.write_bytes(gzip.compress(b"x" * 5000))
    assert sd.gzip_opens(gz) is None
    page = tmp_path / "page.gz"
    page.write_bytes(b"<html>error</html>")
    assert "gzip" in sd.gzip_opens(page)

    z = tmp_path / "a.zip"
    with zipfile.ZipFile(z, "w") as zf:
        zf.writestr("Products.TXT", "x")
    assert sd.zip_holds("products.txt")(z) is None
    assert "holds no" in sd.zip_holds("other.tsv")(z)
    assert "not a zip" in sd.zip_holds("products.txt")(page)

    assert sd.md5_is(None)(page) is None
    assert "does not match" in sd.md5_is("0" * 32)(page)

    hdr = tmp_path / "h.tsv"
    hdr.write_text("hgnc_id\tsymbol\n1\tA\n")
    assert sd.first_line_has("hgnc_id", "symbol")(hdr) is None
    assert sd.first_line_has("gene_symbol")(hdr) is not None


# ------------------------------------------------------------- conversions
def test_omnipath_edges_keep_one_sign_and_drop_the_ambiguous():
    head = "source_genesymbol\ttarget_genesymbol\tis_stimulation\tis_inhibition"
    text = "\n".join([head,
                      "MYC\tTERT\tTrue\tFalse",
                      "RELA\tCYP3A4\tFalse\tTrue",
                      "A\tB\tTrue\tTrue",
                      "C\tD\tFalse\tFalse",
                      "\tE\tTrue\tFalse"])
    assert sd.convert_omnipath(text) == [("MYC", "TERT", 1),
                                         ("RELA", "CYP3A4", -1)]


def test_a_collectri_file_the_package_reads_comes_out(tmp_path, monkeypatch):
    head = "source_genesymbol\ttarget_genesymbol\tis_stimulation\tis_inhibition"
    body = "\n".join(f"TF{i % 300}\tG{i}\tTrue\tFalse" for i in range(20500))

    def fake(url, dest, check=None, retries=3):
        Path(dest).write_text(head + "\n" + body + "\n")
        return dest
    monkeypatch.setattr(sd, "download", fake)
    monkeypatch.setitem(__import__("sys").modules, "decoupler", None)
    sd._fetch_collectri(tmp_path)
    from generoutes import tfactivity
    regs = tfactivity.load_signed_regulons(tmp_path / "collectri.tsv",
                                           min_targets=5)
    assert regs["TF0"]["G0"] == 1
    assert not (tmp_path / "collectri.omnipath.tsv").exists()


def test_a_collectri_with_too_few_edges_is_refused(tmp_path, monkeypatch):
    def fake(url, dest, check=None, retries=3):
        Path(dest).write_text("source_genesymbol\ttarget_genesymbol\t"
                              "is_stimulation\tis_inhibition\nA\tB\tTrue\t"
                              "False\n")
    monkeypatch.setattr(sd, "download", fake)
    monkeypatch.setitem(__import__("sys").modules, "decoupler", None)
    with pytest.raises(sd.SetupError, match="only 1 signed"):
        sd._fetch_collectri(tmp_path)
    assert not (tmp_path / "collectri.tsv").exists()


# ---------------------------------------------------------------- figshare
def test_a_figshare_file_is_found_by_name_with_its_checksum(monkeypatch):
    monkeypatch.setattr(sd, "get_json", lambda url: [
        {"name": "Other.csv", "download_url": "u0", "computed_md5": "a"},
        {"name": "Model.csv", "download_url": "u1", "computed_md5": "b" * 32}])
    assert sd.figshare_file(1, "Model.csv") == ("u1", "b" * 32)


def test_figshare_falls_back_to_the_known_ids_when_the_listing_fails(
        monkeypatch):
    def down(url):
        raise OSError("down")
    monkeypatch.setattr(sd, "get_json", down)
    url, md5 = sd.figshare_file(1, "CRISPRGeneEffect.csv")
    assert url.endswith("/46489063") and md5 == "d10d73692672380b3792231549e97d5a"
    with pytest.raises(sd.SetupError):
        sd.figshare_file(1, "Unknown.csv")


# ---------------------------------------------------------------- planning
def _touch(data, *names):
    for n in names:
        p = data / n
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("x")


def test_a_present_file_is_left_alone_under_any_name_the_loader_accepts(
        tmp_path):
    _touch(tmp_path, "crispr_gene_effect.csv", "Model.csv")
    chosen, have, _ = sd.plan(tmp_path)
    assert "depmap" in {i.key for i in chosen}, "essentials file is missing"
    by = {i.key: i for i in sd.items()}
    assert sd.present(by["trrust"], tmp_path) is False
    _touch(tmp_path, "trrust.tsv")
    assert sd.present(by["trrust"], tmp_path) is True
    # an empty file is not a file
    (tmp_path / "jaspar2024.txt").write_text("")
    assert sd.present(by["jaspar"], tmp_path) is False


def test_the_depmap_fetch_asks_only_for_the_parts_that_are_missing(
        tmp_path, monkeypatch):
    _touch(tmp_path, "CRISPRGeneEffect.csv", "Model.csv")
    asked = []
    monkeypatch.setattr(sd, "fetch_figshare",
                        lambda data, art, name, **k: asked.append(name))
    sd._fetch_depmap(tmp_path)
    assert asked == ["CRISPRInferredCommonEssentials.csv"]


def test_small_leaves_out_the_big_files_and_optional_needs_asking(tmp_path):
    chosen, _, left = sd.plan(tmp_path, small=True)
    keys = {i.key for i in chosen}
    assert "remap" not in keys and "bindingdb" not in keys
    assert "corpora" not in keys
    assert {i.key for i, _ in left} >= {"remap", "bindingdb", "corpora"}
    chosen, _, _ = sd.plan(tmp_path, with_optional=True)
    assert "corpora" in {i.key for i in chosen} and "remap" in {
        i.key for i in chosen}
    chosen, _, _ = sd.plan(tmp_path, only={"remap"}, small=True)
    assert [i.key for i in chosen] == ["remap"], "naming an item overrides"


# ----------------------------------------------------------------- running
def test_a_failure_is_recorded_and_the_rest_still_run(tmp_path, monkeypatch):
    by = {i.key: i for i in sd.items()}

    def boom(data):
        raise sd.SetupError("no route to host")

    def writes(data):
        (data / "jaspar2024.txt").write_text(">MA0001 X")
    monkeypatch.setattr(by["trrust"], "fetch", boom)
    monkeypatch.setattr(by["jaspar"], "fetch", writes)
    monkeypatch.setattr(sd, "items", lambda: list(by.values()))
    out = sd.fetch_all(tmp_path, [by["trrust"], by["jaspar"]])
    assert out == {"trrust": "no route to host", "jaspar": None}


def test_a_fetch_that_writes_nothing_is_reported_not_trusted(tmp_path,
                                                             monkeypatch):
    by = {i.key: i for i in sd.items()}
    monkeypatch.setattr(by["jaspar"], "fetch", lambda data: None)
    monkeypatch.setattr(sd, "items", lambda: list(by.values()))
    out = sd.fetch_all(tmp_path, [by["jaspar"]])
    assert "not in data/" in out["jaspar"]


def test_the_chip_regulon_waits_for_the_peaks_and_the_annotation(tmp_path):
    by = {i.key: i for i in sd.items()}
    out = sd.fetch_all(tmp_path, [by["chip_regulon"]])
    assert out["chip_regulon"].startswith("needs remap, gencode")


# ------------------------------------------- each fetch leaves what it should
def _fake_download(monkeypatch):
    def fake(url, dest, check=None, retries=3):
        dest = Path(dest)
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.suffix == ".zip":
            with zipfile.ZipFile(dest, "w") as zf:
                zf.writestr("rna_tissue_consensus.tsv", "Gene\tTissue\n")
                zf.writestr("products.txt", "x")
        else:
            dest.write_bytes(b"data")
        return dest
    monkeypatch.setattr(sd, "download", fake)
    monkeypatch.setattr(sd, "get_text", lambda url: (
        "gencode.v49.basic.annotation.gtf.gz gencode.v50.basic.annotation"
        ".gtf.gz"))
    monkeypatch.setattr(sd, "get_json", lambda url: [])


@pytest.mark.parametrize("key", [
    "jaspar", "trrust", "clingen", "gnomad", "orangebook", "gencode", "hpa",
    "hgnc", "goa", "lincs", "knockout", "depmap", "prism", "ccle", "remap"])
def test_each_fetch_leaves_a_file_the_loaders_accept(tmp_path, monkeypatch,
                                                     key):
    _fake_download(monkeypatch)
    by = {i.key: i for i in sd.items()}
    assert not sd.present(by[key], tmp_path)
    out = sd.fetch_all(tmp_path, [by[key]])
    assert out == {key: None}, out
    assert sd.present(by[key], tmp_path)


def test_gencode_takes_the_newest_release_in_the_listing(tmp_path,
                                                         monkeypatch):
    got = []
    _fake_download(monkeypatch)
    monkeypatch.setattr(
        sd, "download", lambda url, dest, check=None, retries=3:
        (got.append(url), Path(dest).write_bytes(b"x")))
    monkeypatch.setattr(sd, "get_text", lambda url: (
        "gencode.v9.basic.annotation.gtf.gz gencode.v50.basic.annotation"
        ".gtf.gz gencode.v100.basic.annotation.gtf.gz"))
    sd._fetch_gencode(tmp_path)
    assert got[0].endswith("gencode.v100.basic.annotation.gtf.gz")


def test_the_library_builder_writes_the_columns_the_loader_reads(tmp_path):
    from generoutes import localdrugs

    class Chembl:
        def get(self, url, params=None, **kw):
            off = params["offset"]
            mols = [{"molecule_chembl_id": f"CHEMBL{off + i}",
                     "pref_name": f"DRUG {off + i}",
                     "max_phase": "4.0",   # ChEMBL sends a string
                     "withdrawn_flag": i == 0, "oral": i % 2 == 0,
                     "molecule_structures": {"canonical_smiles": "CCO"}}
                    for i in range(200)]
            return {"molecules": mols,
                    "page_meta": {"next": "x" if off < 2800 else None}}
    sd.build_library(tmp_path, client=Chembl())
    rows = localdrugs.library(str(tmp_path / "approved_library.tsv"))
    assert len(rows) == 3000
    first = tmp_path.joinpath("approved_library.tsv").read_text().splitlines()
    assert first[0].split("\t") == ["name", "smiles", "chembl_id",
                                    "max_phase", "withdrawn", "oral"]


def test_the_library_builder_refuses_a_short_answer(tmp_path):
    class Chembl:
        def get(self, url, params=None, **kw):
            return {"molecules": [], "page_meta": {}}
    with pytest.raises(sd.SetupError, match="only 0"):
        sd.build_library(tmp_path, client=Chembl())
    assert not (tmp_path / "approved_library.tsv").exists()


# ----------------------------------------------------------------- command
def test_the_command_line_checks_names_and_will_not_fetch_unasked(
        tmp_path, capsys, monkeypatch):
    with pytest.raises(SystemExit):
        sd.main(["--data-dir", str(tmp_path), "--only", "nonsense"])
    monkeypatch.setattr(sd.sys.stdin, "isatty", lambda: False)
    fetched = []
    monkeypatch.setattr(sd, "fetch_all",
                        lambda data, chosen: fetched.append(chosen) or {})
    assert sd.main(["--data-dir", str(tmp_path), "--only", "jaspar"]) == 1
    assert not fetched, "no terminal and no --yes means nothing is fetched"
    assert sd.main(["--data-dir", str(tmp_path), "--only", "jaspar",
                    "--list"]) == 0
    assert not fetched
    out = capsys.readouterr().out
    assert "MISSING  jaspar" in out


def test_the_exit_status_says_whether_everything_arrived(tmp_path,
                                                         monkeypatch):
    monkeypatch.setattr(sd, "fetch_all",
                        lambda data, chosen: {"jaspar": "down"})
    assert sd.main(["--data-dir", str(tmp_path), "--only", "jaspar",
                    "--yes"]) == 1
    monkeypatch.setattr(sd, "fetch_all", lambda data, chosen: {"jaspar": None})
    assert sd.main(["--data-dir", str(tmp_path), "--only", "jaspar",
                    "--yes"]) == 0


# ------------------------------------------------- the inventory and wrappers
def test_every_table_the_registry_names_is_fetched_or_explained():
    """A table added to datasets.py must either be fetched here or be on the
    short list of ones that cannot be."""
    from generoutes import datasets
    covered = {a for it in sd.items() for alts in it.have for a in alts}
    uncovered = {ds.key for ds in datasets.registry()
                 if not any(p in covered for p in ds.paths)}
    # chembl: tens of GB and optional; genome: 3 GB and optional;
    # fda_label_cmax: built from a finished run by build_fda_cmax.py;
    # ctd_fame: counted from your own saved runs by build_fame_table.py,
    # so there is nothing to fetch - it does not exist until you have run
    # some genes.
    # drug_protein_abundance: no public file to fetch - decryptE and
    # ProTargetMiner are downloads a person makes and places by hand.
    assert uncovered == {"chembl", "genome", "fda_label_cmax",
                         "ctd_fame", "drug_protein_abundance"}, uncovered


def test_the_batch_wrapper_and_the_launcher_point_at_the_script():
    bat = (ROOT / "scripts" / "setup_data.bat").read_text(encoding="utf-8")
    assert r"scripts\setup_data.py %*" in bat
    assert "SETUP_NOPAUSE" in bat
    launcher = (ROOT / "scripts" / "run_everything.bat").read_text(
        encoding="utf-8")
    assert r"call scripts\setup_data.bat --yes" in launcher
    for name in re.findall(r"data\\([\w.]+)", launcher.split(
            "NEEDSETUP=1", 1)[0].split("set \"NEEDSETUP=\"")[-1]):
        assert name  # the names it checks are real tokens


# ----------------------------------------- the launcher's own commands parse
def _commands(text: str) -> list:
    joined = re.sub(r"\^\r?\n\s*", " ", text)
    return [ln.strip() for ln in joined.splitlines()]


SUBS = {
    "!GENE!": "BUB1B", "!DIR!": "up", "!OUT!": "out", "!PROLIF!": "cycling",
    "!REGULON!": "--regulon-file data/chip_regulon.tsv",
    "!CHIP!": "--chip-bed data/remap2022.bed.gz",
    "!WINDOW!": "--upstream 1000 --downstream 200 --pvalue 1e-5",
    "!TFW!": "--tf-weight 0.7", "!TISS!": "--tissue cns",
    "!EXPR!": "--expression-file data/rna_tissue_consensus.tsv",
    "!GTF!": "--gtf data/gencode.basic.gtf.gz",
    "!LINCS!": "--lincs --lincs-gmt data/up.gmt data/down.gmt "
               "--lincs-geneinfo data/g.txt",
    "!ACTIVITY!": '--tf-activity --activity-tissue "liver" '
                  "--activity-permutations 2000",
    "!SENS!": "--sensitivity-scan", "!DOCKING!": "--docking-results d.json",
    "!UPFACT!": "--upstream-factors 17",
    "!DOSAGE!": "--scenario haploinsufficiency --broken-allele-makes "
                "truncated", "!LIBRARY!": "",
}


def _expand(cmd: str) -> list:
    for k, v in SUBS.items():
        cmd = cmd.replace(k, v)
    cmd = cmd.replace('"%PY%"', "python").replace("\\", "/")
    return shlex.split(cmd, posix=True)


def test_the_launcher_analysis_and_funnel_commands_parse_with_the_real_parsers():
    from generoutes import cli, panel
    text = (ROOT / "scripts" / "run_everything.bat").read_text(
        encoding="utf-8")
    cmds = _commands(text)
    analysis = [c for c in cmds if "-m generoutes --gene" in c
                and "--all-tss" in c and "--no-grade" not in c]
    funnel = [c for c in cmds if "-m generoutes.panel" in c]
    assert len(analysis) == 1 and len(funnel) == 1
    for cmd, parser in ((analysis[0], cli.build_parser()),
                        (funnel[0], panel.build_parser())):
        argv = _expand(cmd)
        argv = argv[argv.index("-m") + 2:]
        args = parser.parse_args(argv)
        kw = cli.run_kwargs_from_args(args)
        assert kw["regulator_route"] is True
        assert kw["upstream_network"] is True
        assert kw["use_genetics"] is True
        assert kw["upstream_factors"] == 17, "the analysis and the funnel agree"


def test_the_launcher_passes_every_analysis_switch_the_cli_offers():
    """Which CLI switches the launcher leaves out is a choice made here, so a
    new analysis flag cannot be left off unnoticed."""
    from generoutes import cli
    text = (ROOT / "scripts" / "run_everything.bat").read_text(
        encoding="utf-8")
    left_out = {
        # a different input source, or found by name under data/
        "--species", "--min-ic", "--matrix-id", "--genome", "--promoter-fasta",
        "--knockout-gmt", "--lincs-file", "--signature-file", "--prism-file",
        # limits whose defaults are what the launcher wants
        "--upstream-drugs-per-node", "--null-shuffles",
        # off on purpose: leads with no direction, see docs/checks.md
        "--regulator-measured",
        # the opposite of a default the launcher wants
        "--no-widen",
        "--no-ctd-route",
        "--depmap-file", "--cache-dir", "--max-requests",
        # opposite of a default, or narrows the output
        "--no-mask-repeats", "--include-unknown-direction",
        "--include-local-only", "--include-investigational",
        "--include-non-fda", "--no-dgidb", "--no-opentargets", "--quiet",
        # about the run's mode, or about a particular gene or patient
        "--demo", "--offline", "--has-ptc", "--help",
        # a declaration about the patient, like --has-ptc: the launcher
        # cannot know it, and defaulting it either way would be a claim
        "--avoid-genotoxic",
    }
    missing = []
    for a in cli.build_parser()._actions:
        opts = [o for o in a.option_strings if o.startswith("--")]
        if opts and opts[0] not in text and opts[0] not in left_out:
            missing.append(opts[0])
    assert not missing, (
        f"the launcher never passes {missing}; add it, or list it as left "
        f"out on purpose")
