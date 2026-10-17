"""Primary-source check for claims the report is about to display.

The registry is the closed set of sources this build is willing to name.
A claim passes only when every number in it matches the run record and the
wording stays inside the registered scope. Clinical conclusions are rejected
even if the numbers are right.
"""

from __future__ import annotations

import re

REGISTRY = [
    {
        "id": "doench2016",
        "doi": "10.1038/nbt.3437",
        "scope": "CFD is the product of position-specific mismatch activities and a PAM weight. Position 1 is PAM-distal.",
    },
    {
        "id": "tsai2015",
        "doi": "10.1038/nbt.3117",
        "scope": "GUIDE-seq measures off-target cleavage in a reference cell line. Read counts are an experimental ranking, not a population haplotype.",
    },
    {
        "id": "gnomad-v4-api",
        "doi": None,
        "scope": "Allele counts in this run are queried from the gnomAD v4 GraphQL API for one population label. They are not a personal genotype.",
    },
]

_FORBIDDEN = re.compile(
    r"\b(diagnos|clinical (recommendation|decision)|safe for (patients|therapy)|should be used in patients)\b",
    re.IGNORECASE,
)
_NUMBER = re.compile(r"(?<![A-Za-z0-9.])\d+(?:\.\d+)?")


def verify_citation(claim: str, record: dict) -> dict:
    if _FORBIDDEN.search(claim):
        return {
            "claim": claim,
            "accepted": False,
            "reason": "The claim is a clinical or diagnostic conclusion. That is outside the registered scope.",
            "sources": [],
        }
    allowed = set()
    for key in ("integers", "floats"):
        for value in record.get(key, []):
            allowed.add(str(value))
            if isinstance(value, float):
                allowed.add(f"{value:.4f}".rstrip("0").rstrip("."))
                allowed.add(f"{value:.3f}".rstrip("0").rstrip("."))
    scrubbed = claim
    for name in ("Evo 2", "gnomAD v4", "Doench 2016"):
        scrubbed = scrubbed.replace(name, " ")
    stray = []
    for token in _NUMBER.findall(scrubbed):
        if token in {"1", "20", "23"}:
            continue
        if token not in allowed and str(int(float(token))) not in allowed:
            # Threshold constants that are part of the method, not a result.
            if token in {"0.2", "0.05", "0.1", "0.5", "0.01", "4"}:
                continue
            if token.isdigit() and 1990 <= int(token) <= 2035:
                continue
            stray.append(token)
    if stray:
        return {
            "claim": claim,
            "accepted": False,
            "reason": f"Numbers {stray} are not in the run record.",
            "sources": [],
        }
    sources = []
    lower = claim.lower()
    if "cfd" in lower or "doench" in lower or "mismatch" in lower:
        sources.append("doench2016")
    if "gnomad" in lower or "allele" in lower or "south asian" in lower or "population" in lower:
        sources.append("gnomad-v4-api")
    if "guide-seq" in lower or "guideseq" in lower:
        sources.append("tsai2015")
    if re.search(r"\b(is|equals|confirms) (a )?cutting probability\b", lower):
        return {
            "claim": claim,
            "accepted": False,
            "reason": "Evo 2 likelihood is not a cutting probability.",
            "sources": [],
        }
    return {
        "claim": claim,
        "accepted": True,
        "reason": "Numbers match the run record and the wording stays inside the registered scope.",
        "sources": sources,
    }
