"""The audit's checks against public lists, with the services replaced."""

import pytest

from generoutes import audit
from generoutes.annotate import NUCLEAR_RECEPTORS

HEADER = "HGNC ID\tApproved symbol\tLocus type\n"


def _group(symbols, *, pseudo=()):
    rows = [f"HGNC:{i}\t{s}\tgene with protein product\n"
            for i, s in enumerate(sorted(symbols))]
    rows += [f"HGNC:9{i}\t{s}\tpseudogene\n" for i, s in enumerate(pseudo)]
    return HEADER + "".join(rows)


def test_nuclear_receptor_check_passes_when_the_lists_agree():
    text = _group(NUCLEAR_RECEPTORS, pseudo=["NR1H5P"])
    r = audit.check_nuclear_receptors(text)
    assert r["missing_from_ours"] == []
    assert r["not_in_reference"] == []
    assert r["reference"] == len(NUCLEAR_RECEPTORS)    # the pseudogene is skipped


def test_nuclear_receptor_check_names_a_receptor_we_lack_and_one_we_add():
    text = _group((NUCLEAR_RECEPTORS - {"RARA"}) | {"NR9Z9"})
    r = audit.check_nuclear_receptors(text)
    assert r["missing_from_ours"] == ["NR9Z9"]
    assert r["not_in_reference"] == ["RARA"]


def test_nuclear_receptor_check_refuses_a_page_that_is_not_the_table():
    """A portal page parses to no rows. That must not read as 'all extra'."""
    with pytest.raises(ValueError):
        audit.check_nuclear_receptors("<html>Please sign in</html>")


class _Dictionary:
    """A stand-in for RCSB: knows a few codes and their types."""

    TYPES = {"MSE": "L-PEPTIDE LINKING", "HOH": "NON-POLYMER",
             "GOL": "NON-POLYMER", "GLC": "D-saccharide, alpha linking"}

    def post(self, url, data):
        ids = data["variables"]["ids"]
        return {"data": {"chem_comps": [
            {"chem_comp": {"id": i, "type": self.TYPES[i]}}
            if i in self.TYPES else None for i in ids]}}


def test_residue_check_reports_unknown_and_mistyped_codes(monkeypatch):
    from generoutes import docking, exploratory

    monkeypatch.setattr(exploratory, "SUGARS", frozenset({"GLC", "MSE", "MAL"}))
    monkeypatch.setattr(docking, "MODIFIED_RESIDUES",
                        frozenset({"MSE", "GOL", "ZZZ"}))
    monkeypatch.setattr(docking, "SOLVENTS", frozenset({"HOH", "WAT"}))
    monkeypatch.setattr(docking, "CRYSTALLISATION_ADDITIVES",
                        frozenset({"GOL", "MSE"}))
    r = audit.check_residue_codes(_Dictionary())
    assert r["MODIFIED_RESIDUES"]["not_found"] == ["ZZZ"]
    assert r["MODIFIED_RESIDUES"]["wrong_type"] == ["GOL"]    # a ligand
    assert r["SOLVENTS"] == {"codes": 2, "not_found": [], "wrong_type": []}
    assert r["CRYSTALLISATION_ADDITIVES"]["wrong_type"] == ["MSE"]
    # MAL is a retired disaccharide code, so absence is not reported
    assert r["exploratory.SUGARS"] == {"codes": 3, "not_found": [],
                                       "wrong_type": ["MSE"]}


def test_the_shipped_residue_tables_hold_only_upper_case_three_character_codes():
    """The offline half of the PDB check: a code that can never match."""
    from generoutes import docking, exploratory

    for table in (docking.MODIFIED_RESIDUES, docking.SOLVENTS,
                  docking.CRYSTALLISATION_ADDITIVES, exploratory.SUGARS):
        assert all(c == c.upper() and 1 <= len(c) <= 3 for c in table)


def test_the_audit_names_every_drug_table_it_claims_to_cover():
    names = audit.all_curated_names()
    for table in ("systemic.CLASSES", "systemic.MITOTIC_DISRUPTORS",
                  "systemic.GENOTOXIC", "indications.CURATED",
                  "standard_of_care.CONDITIONS", "promoterdna.AGENTS",
                  "stabilizers.CURATED", "delivery.PGP_SUBSTRATES",
                  "ligands.ROUTES", "direct.CURATED",
                  "catalogue.CHROMATIN_DRUGS",
                  "retrospective.POSITIVE_CONTROLS"):
        assert names.get(table), table


# ---------------------------------------------------------------- gene symbols
class _Hgnc:
    """A stand-in for the HGNC search: knows which symbols are approved."""

    APPROVED = {"MYC", "WWTR1", "TP53"}

    def __init__(self):
        self.queries = []

    def get(self, url, params=None, headers=None):
        self.queries.append(url)
        assert headers and headers["Accept"] == "application/json"
        asked = url.rsplit("symbol:", 1)[1].strip("()").replace("%20", " ")
        docs = [{"symbol": s} for s in asked.split(" OR ")
                if s in self.APPROVED]
        return {"response": {"docs": docs}}


def test_symbol_check_names_a_stale_symbol_and_the_table_that_holds_it():
    tables = {"systemic.ONCOGENES": {"MYC", "TAZ"},
              "reactivation.IMPRINTED": {"TP53", "SLC22A18"}}
    bad = audit.check_gene_symbols(_Hgnc(), tables)
    assert bad == {"systemic.ONCOGENES": ["TAZ"],
                   "reactivation.IMPRINTED": ["SLC22A18"]}


def test_symbol_check_asks_in_chunks_and_leaves_known_non_symbols_alone():
    client = _Hgnc()
    tables = {"a": {f"G{i}" for i in range(95)} | {"MYC", "SNORD116"}}
    bad = audit.check_gene_symbols(client, tables, chunk=40)
    assert len(client.queries) == 3                    # 97 symbols, 40 a time
    assert "SNORD116" not in bad["a"] and "MYC" not in bad["a"]


def test_the_shipped_gene_lists_use_the_symbols_the_audit_found_current():
    from generoutes import reactivation, systemic

    assert "TAZ" not in systemic.ONCOGENES        # TAFAZZIN's old symbol
    assert "WWTR1" in systemic.ONCOGENES
    assert "SLC67A1" in reactivation.IMPRINTED
    tables = audit.all_curated_symbols()
    assert len(tables) >= 15 and all(tables.values())


# ------------------------------------------------------------ tables on disk
def test_local_checks_report_absent_files_as_not_checked(monkeypatch):
    from generoutes import orangebook, signatures

    monkeypatch.setattr(signatures, "load_geneinfo", lambda *a, **k: 0)
    monkeypatch.setattr(orangebook, "route_vocabulary", lambda *a, **k: set())
    assert audit.check_local_references() == {"landmarks": None,
                                              "orangebook_routes": None}


def test_local_checks_name_a_landmark_the_clue_table_does_not_list(monkeypatch):
    from generoutes import orangebook, signatures

    monkeypatch.setattr(signatures, "KNOWN_LANDMARKS", {"JUN", "VDR"})
    monkeypatch.setattr(signatures, "load_geneinfo", lambda *a, **k: 2)
    monkeypatch.setattr(signatures, "feature_space",
                        lambda g: "landmark" if g == "JUN" else "best inferred")
    monkeypatch.setattr(orangebook, "route_vocabulary",
                        lambda *a, **k: {"ORAL", "SPINAL"})
    r = audit.check_local_references()
    assert r["landmarks"] == ["VDR"]
    assert r["orangebook_routes"]["unclassified"] == ["SPINAL"]
    assert "PARENTERAL" in r["orangebook_routes"]["not_in_this_release"]


def test_the_landmark_fallback_no_longer_lists_genes_clue_calls_inferred():
    from generoutes import signatures

    for gene in ("VDR", "ESR1", "AR", "E2F1", "MKI67", "CEBPB", "SREBF1",
                 "DDIT3", "ATF3", "SOD1", "NQO1"):
        assert gene not in signatures.KNOWN_LANDMARKS
        assert signatures.is_landmark(gene) is None    # unknown, not "yes"


# ----------------------------------------------------------- Orange Book routes
def test_oral_contraceptive_pack_routes_read_as_oral():
    from generoutes import orangebook

    assert orangebook.split_routes("ORAL-28") == ["ORAL"]
    assert orangebook.split_routes("oral-21, rectal") == ["ORAL", "RECTAL"]
    assert orangebook.split_routes("INTRAMUSCULAR, INTRAVENOUS") == [
        "INTRAMUSCULAR", "INTRAVENOUS"]
    assert orangebook.split_routes("") == []


def test_a_drug_sold_only_as_a_28_day_pack_reaches_the_body(
        tmp_path, monkeypatch):
    from generoutes import orangebook

    monkeypatch.setattr(orangebook, "_PRODUCTS", {})
    monkeypatch.setattr(orangebook, "_LOADED", None)
    header = ("Ingredient~DF;Route~Trade_Name~Applicant~Strength~Appl_Type~"
              "Appl_No~Product_No~TE_Code~Approval_Date~RLD~RS~Type~"
              "Applicant_Full_Name")
    row = ("NORETHINDRONE~TABLET;ORAL-28~CAMILA~MAYNE~0.35MG~A~076177~001~"
           "AB~Jan 01, 2003~No~No~RX~Mayne")
    path = tmp_path / "products.txt"
    path.write_text("\n".join([header, row]) + "\n", encoding="utf-8")
    got = orangebook.lookup("norethindrone", path=str(path))
    assert got["reach"] == "systemic" and got["unclassified_routes"] == []
    assert orangebook.route_vocabulary(str(path)) == {"ORAL"}
