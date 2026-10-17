from popoff.policy import review


def test_population_audit_is_in_scope():
    assert review("HBB South Asian research audit of common variants") is None


def test_personal_genome_and_clinical_language_are_refused():
    personal = review("score my genome VCF for this patient and diagnose whether it is safe to inject")
    assert personal is not None
    assert "personal genome" in personal["kinds"]
    assert "clinical framing" in personal["kinds"]
