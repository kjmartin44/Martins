"""Reading peak plasma concentrations from FDA label text.

The sentences below are invented in the shapes real labels use, so the tests
check the rules without copying a label. What the rules refuse matters as much
as what they read: a figure from an interaction study, a special population or a
metabolite is not the blood level of the approved dose.
"""

from __future__ import annotations

import pytest

from generoutes import fdalabel as fl

PLUS = "±"
BAD = "�"                 # what a lost plus-minus sign arrives as


def one(text, heading="Absorption"):
    return fl.candidates(f"{heading} {text}")


def vals(acc):
    return [(c["value"], c["unit"], c["dose"]) for c in acc]


# ----------------------------------------------------------------- reading
def test_a_cmax_with_a_dose_and_unit_is_read():
    acc, rej = one("Following a single oral dose of 20 mg, the mean Cmax was "
                   "40 ng/mL (range 35 to 45 ng/mL).")
    assert vals(acc) == [(40.0, "ng/mL", "20 mg")] and not rej
    assert acc[0]["kind"] == "single dose"


def test_the_spaced_c_max_and_peak_wordings_are_read():
    for text in ("The maximum plasma concentration (C max ) is 792 (49) ng/mL.",
                 "Peak plasma concentrations of 792 ng/mL were reached.",
                 "The mean peak tretinoin concentration was 792 ng/mL.",
                 "The highest observed concentration was 792 ng/mL."):
        acc, _ = one(text)
        assert [c["value"] for c in acc] == [792.0], text


def test_a_list_of_figures_shares_the_unit_and_matches_the_doses():
    acc, _ = one("Peak plasma levels are approximately 0.6, 1.1, 1.4 and "
                 "1.8 mcg/mL for 500, 1,000, 1,500, and 2,000 mg once-daily "
                 "doses, respectively.")
    assert vals(acc) == [(0.6, "ug/mL", "500 mg"), (1.1, "ug/mL", "1,000 mg"),
                         (1.4, "ug/mL", "1,500 mg"), (1.8, "ug/mL", "2,000 mg")]


def test_a_figure_without_its_own_unit_takes_the_later_one():
    acc, _ = one(f"the mean {PLUS} SD peak concentrations after the first "
                 f"dose was 394 {PLUS} 89 and after 1 week of continuous "
                 f"treatment was 138 {PLUS} 139 ng/mL.")
    assert [(c["value"], c["inherited"]) for c in acc] == [(394.0, True),
                                                           (138.0, False)]


def test_a_lost_plus_minus_sign_does_not_break_the_reading():
    acc, _ = one(f"the mean {BAD} SD peak concentration was 138 {BAD} 139 ng/mL")
    assert [c["value"] for c in acc] == [138.0]


def test_the_top_of_a_range_is_not_a_second_peak():
    acc, _ = one("an average peak plasma concentration of 40 ng/mL (range 35 "
                 "to 45 ng/mL) occurred.")
    assert [c["value"] for c in acc] == [40.0]


def test_a_number_that_is_a_dose_or_a_time_is_not_a_concentration():
    acc, _ = one("Peak concentration occurred at 2 hours after 100 mg, and "
                 "was 12 ng/mL.")
    assert [c["value"] for c in acc] == [12.0]
    assert fl.candidates("Absorption Time to reach peak concentration was "
                         "between 1 and 2 hours.")[0] == []


def test_an_auc_unit_is_not_a_dose():
    acc, _ = one("The Cmax is 792 ng/mL and exposure is 2.9 mcg h/mL at 25 mg.")
    assert acc[0]["dose"] == "25 mg"


def test_units_are_normalised_the_way_exposure_reads_them():
    assert fl._unit("mcg/mL") == "ug/mL" and fl._unit(BAD + "g/mL") == "ug/mL"
    assert fl._unit("µM") == "uM" and fl._unit("nmol/L") == "nM"
    assert fl._unit("μg / mL") == "ug/mL"


# ---------------------------------------------------------------- refusing
@pytest.mark.parametrize("sentence,word", [
    ("In patients with renal impairment the Cmax was 900 ng/mL.", "renal"),
    ("In elderly subjects the mean peak concentration was 900 ng/mL.",
     "elderly"),
    ("Coadministration with ketoconazole raised the Cmax to 900 ng/mL.",
     "Coadministration"),
    ("With a high-fat meal the Cmax was 900 ng/mL.", "high-fat"),
    ("The Cmax of the active metabolite was 900 ng/mL.", "metabolite"),
    ("In children the Cmax was 900 ng/mL.", "children"),
    ("In female subjects the peak concentration was 900 ng/mL.", "female"),
])
def test_a_special_population_or_condition_is_set_aside_with_its_reason(
        sentence, word):
    acc, rej = one(sentence)
    assert not acc and len(rej) == 1 and word in rej[0]["why"]


def test_a_value_under_a_population_or_interaction_heading_is_set_aside():
    text = ("Absorption The mean Cmax was 40 ng/mL after 20 mg. "
            "Renal Impairment The mean Cmax was 90 ng/mL after 20 mg. "
            "Drug Interaction Studies The Cmax was 70 ng/mL after 20 mg.")
    acc, rej = fl.candidates(text)
    assert [c["value"] for c in acc] == [40.0]
    assert len(rej) == 2 and all("heading" in r["why"] for r in rej)


def test_a_sentence_without_a_peak_phrase_is_not_read():
    assert fl.candidates("Absorption The half-life is 12 ng/mL.")[0] == []
    assert fl.candidates("")[0] == []


# --------------------------------------------------------------- the choice
class _Fda:
    """A stand-in for the HTTP client serving one openFDA response."""

    def __init__(self, labels):
        self.labels, self.calls = labels, []

    def get(self, url, params=None, **kw):
        self.calls.append(params)
        return {"results": self.labels}


def _label(text, generic="Tretinoin", route="ORAL", extra=None):
    lab = {"openfda": {"generic_name": [generic], "substance_name": [generic],
                       "brand_name": ["Brand"], "route": [route]},
           "pharmacokinetics": [text], "effective_time": "20240131",
           "set_id": "abc"}
    lab.update(extra or {})
    return lab


def test_a_drug_gets_a_row_with_its_source_dose_and_sentence():
    client = _Fda([_label("Absorption Following a single oral dose of 20 mg "
                          "the mean Cmax was 40 ng/mL. Distribution The drug "
                          "is 99% bound to plasma proteins.")])
    got = fl.read_drug(client, "Tretinoin", 300.0)
    assert got["status"] == "row"
    row = got["row"]
    assert row["cmax_nM"] == pytest.approx(40 * 1e3 / 300.0, rel=1e-3)
    assert row["dose"] == "20 mg" and row["fu"] == 0.01
    assert "FDA label Brand (oral), 2024-01-31" in row["source"]
    assert "Cmax was 40 ng/mL" in row["note"] and "99% bound" in row["note"]


def test_labels_that_repeat_a_figure_count_it_once():
    text = "Absorption The mean Cmax was 40 ng/mL after 20 mg."
    got = fl.read_drug(_Fda([_label(text), _label(text)]), "Tretinoin", 300.0)
    assert got["status"] == "row" and len(got["accepted"]) == 1
    assert got["accepted"][0]["n_labels"] == 2


def test_figures_that_disagree_widely_go_to_review_not_to_the_table():
    text = ("Absorption After 10 mg the mean Cmax was 4 ng/mL. "
            "After 40 mg the mean Cmax was 400 ng/mL.")
    got = fl.read_drug(_Fda([_label(text)]), "Tretinoin", 300.0)
    assert got["status"] == "review" and got["row"] is None
    assert "more than 30-fold" in got["reason"]


def test_a_mass_unit_without_a_weight_is_not_guessed():
    got = fl.read_drug(_Fda([_label("Absorption The mean Cmax was 40 ng/mL "
                                    "after 20 mg.")]), "Tretinoin", None)
    assert got["status"] == "review" and "molecular weight" in got["reason"]
    got = fl.read_drug(_Fda([_label("Absorption The mean Cmax was 2 uM after "
                                    "20 mg.")]), "Tretinoin", None)
    assert got["row"]["cmax_nM"] == 2000


def test_a_label_with_no_peak_sentence_gives_nothing_to_read():
    got = fl.read_drug(_Fda([_label("Absorption The drug is absorbed with a "
                                    "Tmax of 2 hours.")]), "Tretinoin", 300.0)
    assert got["status"] == "none" and "no peak concentration" in got["reason"]


def test_everything_set_aside_is_review_and_says_why():
    got = fl.read_drug(_Fda([_label("Absorption In renal impairment the Cmax "
                                    "was 90 ng/mL.")]), "Tretinoin", 300.0)
    assert got["status"] == "review" and "renal" in got["reason"]
    assert got["rejected"][0]["label"].startswith("FDA label")


def test_combination_labels_and_other_drugs_are_not_used():
    combo = _label("Absorption The Cmax was 40 ng/mL.",
                   generic="amlodipine and atorvastatin")
    other = _label("Absorption The Cmax was 40 ng/mL.", generic="Isotretinoin")
    got = fl.read_drug(_Fda([combo, other]), "Tretinoin", 300.0)
    assert got["labels"] == 0 and got["status"] == "none"
    assert "no single ingredient label" in got["reason"]


def test_a_bare_name_finds_a_label_spelled_with_its_salt():
    lab = _label("Absorption The Cmax was 40 ng/mL after 20 mg.",
                 generic="Lapatinib Ditosylate")
    assert fl.read_drug(_Fda([lab]), "Lapatinib", 500.0)["status"] == "row"


def test_a_failed_search_is_reported_not_raised():
    class Down(_Fda):
        def get(self, url, params=None, **kw):
            raise RuntimeError("503 unavailable")
    got = fl.read_drug(Down([]), "Tretinoin", 300.0)
    assert got["status"] == "none" and got["labels"] == 0


def test_clinical_pharmacology_is_used_when_there_is_no_pk_section():
    lab = _label("", extra={"clinical_pharmacology": [
        "Absorption The mean Cmax was 40 ng/mL after 20 mg."]})
    lab.pop("pharmacokinetics")
    assert fl.read_drug(_Fda([lab]), "Tretinoin", 300.0)["status"] == "row"


# ---------------------------------------------------------- protein binding
@pytest.mark.parametrize("text,fu,upper", [
    ("It is approximately 98.4% bound to plasma proteins.", 0.016, False),
    ("Protein binding is greater than 95%.", 0.05, True),
    ("The drug is >99% protein bound.", 0.01, True),
])
def test_an_unbound_fraction_is_read_from_a_binding_statement(text, fu, upper):
    got, note = fl.unbound_fraction(text)
    assert got == pytest.approx(fu)
    assert ("at most" in note) is upper


def test_conflicting_or_absent_binding_gives_no_fraction():
    assert fl.unbound_fraction("It is 90% bound. Elsewhere 95% bound.")[0] is None
    assert fl.unbound_fraction("Nothing about binding.")[0] is None
    assert fl.unbound_fraction("It is 6% to 12% bound to plasma proteins.")[0] \
        is None, "a range is not one fraction"
    assert fl.unbound_fraction("The metabolite is 99% bound.")[0] is None
