"""Every field of a result has to survive `json.dumps`.

`collapse_grades` keys its alternates by `(drug, factor)`, because the same
drug at two factors is two different routes. A tuple is not a JSON key, so
a run where any salt form collapsed used to finish the whole analysis and
then die in `write_outputs`, losing it.
"""

import json

import pytest

from generoutes import congeners as cong_mod
from generoutes import pipeline, report


@pytest.fixture
def _one_salt_collapsed(monkeypatch):
    """Make the run look like it folded a salt form into its parent."""
    real = cong_mod.collapse_grades

    def fake(grades):
        kept, _info = real(grades)
        return kept, {("Liothyronine", "THRB"): ["Liothyronine Sodium"]}

    monkeypatch.setattr(cong_mod, "collapse_grades", fake)


def test_a_run_that_collapsed_a_salt_form_still_writes(_one_salt_collapsed,
                                                       tmp_path):
    result = pipeline.run("DEMOGENE", demo=True, offline=True)
    blocks = [b for b in (result.get("results") or [])
              if b.get("salt_forms_collapsed")]
    assert blocks, "the fixture should have folded one salt form"
    collapsed = blocks[0]["salt_forms_collapsed"]
    assert all(isinstance(k, str) for k in collapsed)

    written = report.write_outputs(result, str(tmp_path / "out"))
    assert written
    saved = json.loads((tmp_path / "out" / "results.json")
                       .read_text(encoding="utf-8"))
    assert saved["results"][0]["salt_forms_collapsed"] == collapsed
