from popoff.cfd import cfd_score


def test_published_poly_g_example():
    # Doench calculator doctest: three PAM-proximal A mismatches on a poly-G guide.
    score = cfd_score("G" * 20, "G" * 17 + "AAA" + "GGG")
    assert abs(score - 0.4635989007074176) < 1e-9


def test_perfect_match_is_one():
    assert cfd_score("ATCGATCGATCGATCGATCG", "ATCGATCGATCGATCGATCGAGG") == 1.0
