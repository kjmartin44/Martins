import importlib.util
from pathlib import Path

from generoutes import localdrugs

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "build_hub_library.py"
spec = importlib.util.spec_from_file_location("build_hub_library", SCRIPT)
hub = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hub)

DRUGS = (
    "!Source\tThe Drug Repurposing Hub\t\t\t\t\n"
    "pert_iname\tclinical_phase\tmoa\ttarget\tdisease_area\tindication\n"
    "aspirin\tLaunched\tCOX inhibitor\tPTGS1|PTGS2\t\t\n"
    "mystery\tPhase 1\tunknown\tABC1\t\t\n"
)
SAMPLES = (
    "!Source\tThe Drug Repurposing Hub\n"
    "broad_id\tpert_iname\tsmiles\n"
    "BRD-1\tAspirin\tCC(=O)Oc1ccccc1C(=O)O\n"
)


def test_joins_structures_and_drops_drugs_without_one(tmp_path):
    drugs, samples, out = (tmp_path / n for n in ("d.txt", "s.txt", "o.tsv"))
    drugs.write_text(DRUGS)
    samples.write_text(SAMPLES)
    assert hub.main([str(drugs), str(samples), "--out", str(out)]) == 0
    rows = localdrugs.library(str(out))
    assert [r["name"] for r in rows] == ["aspirin"]
    assert rows[0]["targets"] == ["PTGS1", "PTGS2"]
    assert rows[0]["smiles"].startswith("CC(=O)O")
