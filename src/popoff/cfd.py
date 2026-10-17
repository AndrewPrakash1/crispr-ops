"""Doench 2016 CFD score.

The mismatch and PAM weights are the published coefficients (Supplementary
Table 19 of Doench et al., Nature Biotechnology 2016, doi:10.1038/nbt.3437).
Position 1 is PAM-distal. The score is the product of one weight per mismatch
and one PAM-dinucleotide weight. A perfect spacer with a GG PAM scores 1.

This is a reimplementation of that product, checked against the published
calculator's known value for a 3-mismatch poly-G guide (0.4635989007).
"""

from __future__ import annotations

import json
from functools import lru_cache
from importlib.resources import files

_COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "U": "A"}


@lru_cache(maxsize=1)
def _tables() -> tuple[dict[str, float], dict[str, float]]:
    raw = json.loads(files("popoff").joinpath("cfd_tables.json").read_text())
    return raw["MISMATCH_SCORES"], raw["PAM_SCORES"]


def revcomp(seq: str) -> str:
    return "".join(_COMP.get(base, "N") for base in reversed(seq.upper()))


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def _key(guide_base: str, off_base: str, position: int) -> str:
    """r<guide>:d<complement of the written protospacer base>, position.

    The written protospacer is the non-target strand, same sequence as the
    guide when they match. The DNA base in the hybrid is its complement, which
    is how the published calculator builds the key.
    """
    guide_base = guide_base.upper().replace("T", "U")
    off_base = off_base.upper().replace("T", "U")
    return f"r{guide_base}:d{_COMP[off_base]},{position}"


def cfd_score(guide: str, site23: str) -> float:
    """Score a 20 nt guide against a 23 nt site (protospacer + PAM)."""
    guide = guide.upper()
    site23 = site23.upper()
    if len(guide) != 20 or len(site23) != 23:
        raise ValueError("CFD expects a 20 nt guide and a 23 nt site")
    if any(base not in "ACGT" for base in guide + site23):
        return 0.0
    mismatch, pam = _tables()
    score = 1.0
    for i, (g, off) in enumerate(zip(guide, site23[:20]), start=1):
        if g == off:
            continue
        score *= mismatch[_key(g, off, i)]
    score *= pam[site23[-2:]]
    return score


def pair_weights(guide: str, ref23: str, alt23: str | None) -> list[dict]:
    """Per-position mismatch identity for the reference and alternate windows."""
    mismatch, _pam = _tables()
    guide = guide.upper()
    rows = []
    for i, g in enumerate(guide, start=1):
        ref_b = ref23[i - 1].upper()
        alt_b = alt23[i - 1].upper() if alt23 else None
        if g == ref_b and (alt_b is None or g == alt_b):
            continue
        row = {
            "position": i,
            "guide": g,
            "ref_base": ref_b,
            "alt_base": alt_b,
            "ref_weight": 1.0 if g == ref_b else mismatch[_key(g, ref_b, i)],
            "alt_weight": None,
        }
        if alt_b is not None:
            row["alt_weight"] = 1.0 if g == alt_b else mismatch[_key(g, alt_b, i)]
        if ref_b != g or (alt_b is not None and alt_b != g) or ref_b != alt_b:
            rows.append(row)
    return rows


def one_snp_envelope(guide: str, site23: str) -> tuple[float, float]:
    """Min and max CFD reachable by a single substitution in the 23 nt window."""
    site = list(site23.upper())
    scores = [cfd_score(guide, "".join(site))]
    for i in range(23):
        original = site[i]
        for base in "ACGT":
            if base == original:
                continue
            site[i] = base
            scores.append(cfd_score(guide, "".join(site)))
        site[i] = original
    return min(scores), max(scores)
