from popoff.resolver import analyze, resolve_sites


def test_quiet_deltas_are_not_a_safety_verdict():
    result = analyze(
        cfd_ref=0.04, cfd_alt=0.04, checker_ref=0.04, checker_alt=0.04,
        af=0.25, coordinates_agree=True,
    )
    assert result["action"] == "ACCEPT"
    assert "not a safety verdict" in result["reason"]


def test_window_mismatch_is_investigated_before_pam_expansion():
    result = analyze(
        cfd_ref=0.4, cfd_alt=0.1, checker_ref=0.2, checker_alt=0.1,
        af=0.3, coordinates_agree=True, windows_match=False,
    )
    assert result["action"] == "INVESTIGATE"


def test_one_sided_delta_asks_for_pam_expansion():
    result = analyze(
        cfd_ref=0.8, cfd_alt=0.2, checker_ref=0.5, checker_alt=0.49,
        af=0.4, coordinates_agree=True,
    )
    assert result["action"] == "ESCALATE_GPU_ENUMERATION"


def test_agreed_large_delta_asks_for_evo2():
    result = analyze(
        cfd_ref=0.8, cfd_alt=0.2, checker_ref=0.8, checker_alt=0.2,
        af=0.4, coordinates_agree=True,
    )
    assert result["action"] == "ESCALATE_EVO2"
    assert "cutting probability" in result["reason"]


def test_coordinate_mismatch_outranks_a_quiet_site():
    quiet = analyze(
        cfd_ref=0.1, cfd_alt=0.1, checker_ref=0.1, checker_alt=0.1,
        af=0.01, coordinates_agree=True,
    )
    wave = resolve_sites([quiet], coordinates_agree=False)
    assert wave["action"] == "INVESTIGATE"
