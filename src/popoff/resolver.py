"""Compare the engine with the second checker.

The result is a next action. It is not a safety verdict. "Accept" means the
two scorers do not justify another tool call. It does not mean a guide is
safe to use.
"""

from __future__ import annotations

REL_QUIET = 0.1
REL_SPLIT = 0.3
REL_HIGH = 0.5


def _rel(left: float | None, right: float | None) -> float:
    if left is None or right is None:
        return 0.0
    return abs(left - right) / max(abs(left), 1e-6)


def analyze(
    *,
    cfd_ref: float | None,
    cfd_alt: float | None,
    checker_ref: float | None,
    checker_alt: float | None,
    af: float | None,
    coordinates_agree: bool,
    windows_match: bool = True,
) -> dict:
    """One site. Inputs are scores the tools already returned."""
    if not coordinates_agree or not windows_match:
        return {
            "status": "DISAGREEMENT",
            "action": "INVESTIGATE",
            "reason": (
                "The engine and the second checker do not describe the same window. "
                "Check window start, strand orientation, and which allele was applied "
                "before any PAM expansion."
            ),
            "af": af,
        }
    cfd_rel = _rel(cfd_ref, cfd_alt)
    checker_rel = _rel(checker_ref, checker_alt)
    split = (cfd_rel > REL_SPLIT and checker_rel < REL_QUIET) or (
        checker_rel > REL_SPLIT and cfd_rel < REL_QUIET
    )
    if split:
        return {
            "status": "DISAGREEMENT",
            "action": "ESCALATE_GPU_ENUMERATION",
            "reason": (
                "One scorer's haplotype delta is large and the other's is small. "
                "An alternate PAM is a hypothesis. Search NAG and NGA for this guide only."
            ),
            "af": af,
            "cfd_relative_change": cfd_rel,
            "checker_relative_change": checker_rel,
        }
    if cfd_rel > REL_HIGH and checker_rel > REL_HIGH:
        return {
            "status": "HIGH_DELTA",
            "action": "ESCALATE_EVO2",
            "reason": (
                "Both scorers move on this haplotype. Evo 2 can score the window as a "
                "sequence-context prior. That score is not a cutting probability."
            ),
            "af": af,
            "cfd_relative_change": cfd_rel,
            "checker_relative_change": checker_rel,
        }
    if cfd_rel < REL_QUIET and checker_rel < REL_QUIET:
        return {
            "status": "CONSISTENT",
            "action": "ACCEPT",
            "reason": "Both haplotype deltas are under 10% relative. This is not a safety verdict.",
            "af": af,
            "cfd_relative_change": cfd_rel,
            "checker_relative_change": checker_rel,
        }
    return {
        "status": "AMBIGUOUS",
        "action": "ACCEPT_WITH_NOTE",
        "reason": "The scorers move by different moderate amounts. No extra tool is required.",
        "af": af,
        "cfd_relative_change": cfd_rel,
        "checker_relative_change": checker_rel,
    }


_RANK = {
    "INVESTIGATE": 0,
    "ESCALATE_GPU_ENUMERATION": 1,
    "ESCALATE_EVO2": 2,
    "ACCEPT_WITH_NOTE": 3,
    "ACCEPT": 4,
}


def resolve_sites(site_results: list[dict], *, coordinates_agree: bool) -> dict:
    """Wave-level action. The strongest disagreement wins."""
    if not coordinates_agree:
        return {
            "status": "DISAGREEMENT",
            "action": "INVESTIGATE",
            "reason": "The two enumerators returned different coordinates.",
            "n_sites": len(site_results),
        }
    if not site_results:
        return {
            "status": "CONSISTENT",
            "action": "ACCEPT",
            "reason": "No scored sites to compare.",
            "n_sites": 0,
        }
    chosen = min(site_results, key=lambda item: _RANK[item["action"]])
    return {
        "status": chosen["status"],
        "action": chosen["action"],
        "reason": chosen["reason"],
        "n_sites": len(site_results),
        "n_disagreement": sum(1 for item in site_results if item["status"] == "DISAGREEMENT"),
        "n_high_delta": sum(1 for item in site_results if item["status"] == "HIGH_DELTA"),
    }
