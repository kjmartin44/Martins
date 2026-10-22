"""Regression tests for input, export, CLI, and HTTP behavior."""

import gzip
import io
import json
import urllib.error
from concurrent.futures import ThreadPoolExecutor

import pytest

from generoutes import cli, pipeline, shortlist, motifs
from generoutes.http import Client, HttpError, _retry_delay


@pytest.mark.parametrize("kwargs", [
    {"gene": " "}, {"upstream": -1}, {"downstream": -1},
    {"upstream": 0, "downstream": 0}, {"upstream": 1.5},
    {"pvalue": 0}, {"pvalue": 1}, {"pvalue": float("nan")},
    {"tf_weight": 1.1}, {"tf_weight": float("inf")},
    {"min_ic": -1}, {"goal": "sideways"}, {"tissue": "other"},
    {"proliferation": "other"}, {"upstream_depth": 0},
])
def test_invalid_inputs_fail_before_client_creation(monkeypatch, kwargs):
    def unexpected(*a, **k):
        pytest.fail("invalid inputs reached the HTTP client")
    monkeypatch.setattr(pipeline, "Client", unexpected)
    with pytest.raises(ValueError):
        pipeline.run(**({"gene": "SMN1"} | kwargs))


def test_missing_file_fails_before_network(tmp_path, monkeypatch):
    monkeypatch.setattr(pipeline, "Client", lambda **kw: pytest.fail("client created"))
    with pytest.raises(FileNotFoundError, match="pfm_file"):
        pipeline.run("SMN1", pfm_file=str(tmp_path / "absent.pfm"))


@pytest.fixture
def demo(tmp_path):
    return pipeline.run("DEMOGENE", demo=True, cache_dir=str(tmp_path / "cache"))


def test_saved_and_printed_report_use_same_limit(demo, tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(cli, "run", lambda *a, **kw: demo)
    assert cli.main(["--gene", "DEMOGENE", "--demo", "--top", "2",
                     "--out", str(tmp_path / "out")]) == 0
    saved = (tmp_path / "out" / "report.md").read_text(encoding="utf-8")
    assert capsys.readouterr().out == saved + "\n"
    assert saved == pipeline.markdown_report(demo, top_n=2)
    assert json.loads((tmp_path / "out" / "results.json").read_text()) == demo


@pytest.mark.parametrize("strand, expected_start, expected_end, expected_strands", [
    (1, 102, 107, ["+", "-"]), (-1, 113, 118, ["-", "+"]),
])
def test_bed_coordinates_and_genomic_strands(
        demo, tmp_path, strand, expected_start, expected_end, expected_strands):
    block = demo["results"][0]
    block["promoter"].update(start=101, end=120, strand=strand)
    hit = dict(block["hits"][0], start=2, end=7)
    block["hits"] = [dict(hit, strand="+"), dict(hit, strand="-")]
    pipeline.write_outputs(demo, tmp_path)
    rows = [line.split("\t") for line in
            (tmp_path / "sites.bed").read_text().splitlines()[1:]]
    assert [(int(r[1]), int(r[2])) for r in rows] == [
        (expected_start, expected_end)] * 2
    assert [r[5] for r in rows] == expected_strands


def test_cli_error_has_no_traceback(tmp_path, capsys):
    assert cli.main(["--gene", "SMN1", "--pvalue", "nan",
                     "--out", str(tmp_path / "out")]) == 1
    captured = capsys.readouterr()
    assert "pvalue" in captured.err
    assert "Traceback" not in captured.err
    assert not (tmp_path / "out").exists()


def test_shortlist_zero_overrides_preset(tmp_path, monkeypatch):
    source = tmp_path / "results.json"
    source.write_text('{"results": []}', encoding="utf-8")
    observed = {}
    def capture(data, **criteria):
        observed.update(criteria)
        return []
    monkeypatch.setattr(shortlist, "filter_routes", capture)
    assert shortlist.main([str(source), "--preset", "specific",
                           "--max-steps", "0", "--min-specificity", "0"]) == 0
    assert observed["max_steps"] == 0
    assert observed["min_specificity"] == 0


def test_cache_separates_method_format_and_headers(tmp_path, monkeypatch):
    client = Client(tmp_path, min_interval=0)
    calls = []
    def respond(req, **kw):
        calls.append(req)
        response = io.BytesIO(b'{"value": 1}')
        response.headers = {}
        return response
    monkeypatch.setattr("urllib.request.urlopen", respond)
    url = "https://example.org/data"
    assert client.get(url) == {"value": 1}
    assert client.get(url, parse_json=False) == '{"value": 1}'
    assert client.request(url, method="POST") == {"value": 1}
    assert client.get(url, headers={"Accept": "text/plain"}) == {"value": 1}
    client.offline = True
    assert client.get(url) == {"value": 1}
    assert len(calls) == 4


def test_cache_writers_do_not_collide(tmp_path):
    client = Client(tmp_path)
    path = client._path("https://example.org/shared", None)
    with ThreadPoolExecutor(max_workers=6) as pool:
        list(pool.map(lambda i: client._write_cache(path, {"payload": i}), range(30)))
    with gzip.open(path, "rt") as fh:
        assert json.load(fh)["payload"] in range(30)
    assert not list(path.parent.glob("*.tmp"))


def test_malformed_cache_is_a_miss(tmp_path):
    client = Client(tmp_path, offline=True)
    path = client._path("https://example.org/data", None)
    with gzip.open(path, "wt") as fh:
        json.dump({"unexpected": "shape"}, fh)
    with pytest.raises(HttpError, match="offline"):
        client.get("https://example.org/data")


def test_query_parameters_preserve_existing_query_and_fragment(tmp_path, monkeypatch):
    client = Client(tmp_path)
    monkeypatch.setattr(client, "request", lambda url, **kw: url)
    assert client.get("https://example.org/?old=1#anchor", {"new": [2, 3]}) == (
        "https://example.org/?old=1&new=2&new=3#anchor")


@pytest.mark.parametrize("value, expected", [("3", 3), ("-1", 0),
    ("invalid", 4), ("nan", 4), ("inf", 4), (None, 4)])
def test_retry_after_values(value, expected):
    assert _retry_delay(value, 2) == expected


def test_retry_after_http_date(monkeypatch):
    monkeypatch.setattr("time.time", lambda: 0)
    assert _retry_delay("Thu, 01 Jan 1970 00:01:00 GMT", 0) == 60


def test_redirect_cycle_is_bounded(tmp_path, monkeypatch):
    calls = []
    def redirect(req, **kw):
        calls.append(req.full_url)
        target = "/b" if req.full_url.endswith("/a") else "/a"
        raise urllib.error.HTTPError(req.full_url, 307, "redirect",
                                     {"Location": target}, None)
    monkeypatch.setattr("urllib.request.urlopen", redirect)
    client = Client(tmp_path, min_interval=0)
    with pytest.raises(HttpError, match="too many redirects"):
        client.post("https://example.org/a", {}, max_redirects=2)
    assert len(calls) == 3


def test_final_failed_attempt_does_not_sleep(tmp_path, monkeypatch):
    def unavailable(req, **kw):
        raise urllib.error.HTTPError(req.full_url, 503, "unavailable", {}, None)
    monkeypatch.setattr("urllib.request.urlopen", unavailable)
    monkeypatch.setattr("time.sleep", lambda t: pytest.fail("unnecessary sleep"))
    with pytest.raises(HttpError, match="after 1 attempts"):
        Client(tmp_path, min_interval=0).get("https://example.org", retries=1)


def test_matrix_filter_applies_to_local_file(tmp_path):
    pfm = tmp_path / "motifs.pfm"
    rows = "A [ 8 1 ]\nC [ 1 8 ]\nG [ 1 1 ]\nT [ 1 1 ]\n"
    pfm.write_text(">MA1.1 TF1\n" + rows + ">MA2.1 TF2\n" + rows, encoding="utf-8")
    client = Client(tmp_path / "cache", offline=True)
    loaded = motifs.load_motifs(client, pfm, ["MA2.1"])
    assert [m.matrix_id for m in loaded] == ["MA2.1"]
    with pytest.raises(ValueError, match="MA3.1"):
        motifs.load_motifs(client, pfm, ["MA2.1", "MA3.1"])


def test_demo_blocks_network_even_with_chip_api(tmp_path, monkeypatch):
    monkeypatch.setattr("urllib.request.urlopen",
                        lambda *a, **kw: pytest.fail("demo used the network"))
    result = pipeline.run("DEMOGENE", demo=True, chip_api=True,
                          cache_dir=str(tmp_path / "cache"))
    assert result["query"]["demo"] is True
    assert result["results"][0]["promoter"]["assembly"] == "DEMO"


def test_failed_cache_write_keeps_successful_response(tmp_path, monkeypatch, caplog):
    calls = []
    def respond(req, **kw):
        calls.append(req)
        response = io.BytesIO(b'{"value": 1}')
        response.headers = {}
        return response
    def denied(*a, **kw):
        raise PermissionError("file is in use")
    monkeypatch.setattr("urllib.request.urlopen", respond)
    monkeypatch.setattr("pathlib.Path.replace", denied)
    client = Client(tmp_path, min_interval=0)
    assert client.get("https://example.org/data") == {"value": 1}
    assert len(calls) == 1
    assert "could not cache response" in caplog.text
    assert not list(tmp_path.rglob("*.tmp"))
# ===================================================================== #
# Printing a report on a console that cannot spell it.                  #
# ===================================================================== #
def _cp1252_stream():
    """A stdout like a default Windows console: cp1252, reconfigurable."""
    return io.TextIOWrapper(io.BytesIO(), encoding="cp1252",
                            errors="strict", line_buffering=False)


def test_explain_prints_a_report_holding_an_arrow_on_a_cp1252_console(
        tmp_path, monkeypatch, capsysbinary):
    """`generoutes-explain` died with a UnicodeEncodeError on any result
    whose notes held an arrow, which the route chains build routinely. It
    read as a bug in the analysis rather than in how it was printed."""
    from generoutes import explain, provenance

    result = {"query": {"gene": "G", "goal": "up"},
              "results": [{"all_routes": {"routes": [{
                  "drug": "Bexarotene", "category": "TF: orthosteric ligand",
                  "mechanism": "agonist", "acts_on": "RXRA", "tier": 2,
                  "directness": 1, "score": 50.0, "approved": True,
                  "tier_label": "measured",
                  "source": "drug_scoring",
                  "axes": {"strength": 0.5, "specificity": 0.5,
                           "efficiency": 0.5},
                  "notes": ["Bexarotene \u2192 RXRA \u2192 G"]}]}}]}
    path = tmp_path / "results.json"
    path.write_text(json.dumps(result), encoding="utf-8")

    monkeypatch.setattr("sys.stdout", _cp1252_stream())
    monkeypatch.setattr("sys.stderr", _cp1252_stream())
    assert explain.main([str(path), "bexarotene"]) == 0


def test_shortlist_does_the_same(tmp_path, monkeypatch):
    from generoutes import shortlist as sl

    result = {"query": {"gene": "G", "goal": "up"},
              "results": [{"all_routes": {"routes": [{
                  "drug": "Bexarotene", "category": "TF: orthosteric ligand",
                  "mechanism": "agonist", "acts_on": "RXRA", "tier": 2,
                  "directness": 1, "score": 50.0, "approved": True,
                  "tier_label": "measured",
                  "source": "drug_scoring",
                  "axes": {"strength": 0.5, "specificity": 0.5,
                           "efficiency": 0.5},
                  "notes": ["A \u2192 B"]}]}}]}
    path = tmp_path / "results.json"
    path.write_text(json.dumps(result), encoding="utf-8")

    monkeypatch.setattr("sys.stdout", _cp1252_stream())
    monkeypatch.setattr("sys.stderr", _cp1252_stream())
    assert sl.main([str(path)]) == 0


def test_the_console_guard_leaves_a_stream_it_cannot_change_alone():
    """A pipe opened in binary has no reconfigure; that is not an error."""
    from generoutes.provenance import use_utf8_console

    class Fixed:
        encoding = "ascii"

    import sys
    before = (sys.stdout, sys.stderr)
    try:
        sys.stdout = sys.stderr = Fixed()
        use_utf8_console()          # must not raise
    finally:
        sys.stdout, sys.stderr = before
