"""Research-use policy. Personal genomes and clinical framings are refused."""

from __future__ import annotations

import re

# The gene name HBB and the words "sickle cell" describe the locus under study.
# They are not, by themselves, a request to diagnose a person.
_REFUSAL = [
    (r"\bmy (genome|dna|vcf|sample)\b", "personal genome"),
    (r"\bpersonal genome\b", "personal genome"),
    (r"\bpatient('s|s)?\b", "clinical framing"),
    (r"\bthis (individual|person|subject)\b", "personal genome"),
    (r"\bdiagnos(e|is|tic)\b", "clinical framing"),
    (r"\bshould (we|i) (treat|edit|inject|infuse)\b", "clinical framing"),
    (r"\bclinical (decision|use|recommendation)\b", "clinical framing"),
    (r"\bsafe to (edit|inject|infuse|treat)\b", "clinical framing"),
    (r"\b\.vcf\b", "personal genome"),
    (r"\bupload(ed)? (my )?genome\b", "personal genome"),
]


def review(text: str) -> dict | None:
    """Return a refusal record, or None if the request stays in scope."""
    haystack = text or ""
    hits = []
    for pattern, kind in _REFUSAL:
        if re.search(pattern, haystack, flags=re.IGNORECASE):
            hits.append(kind)
    if not hits:
        return None
    kinds = sorted(set(hits))
    return {
        "refused": True,
        "kinds": kinds,
        "reason": (
            "Refused. PopOff audits how off-target scores move under a named "
            "population's common variants. It does not accept a personal genome "
            "and it does not make a clinical or diagnostic claim."
        ),
    }
