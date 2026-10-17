from popoff.citations import verify_citation


def test_year_and_negated_evo2_claim_are_allowed():
    record = {"integers": [1, 0], "floats": []}
    claim = "CFD is the Doench 2016 product. Evo 2 is not a cutting probability."
    verdict = verify_citation(claim, record)
    assert verdict["accepted"]


def test_clinical_and_evo2_as_truth_are_rejected():
    record = {"integers": [], "floats": []}
    assert verify_citation("Evo 2 confirms a cutting probability.", record)["accepted"] is False
    assert verify_citation("This guide is safe for patients.", record)["accepted"] is False
