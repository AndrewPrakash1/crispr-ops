"""Independent site check.

Two jobs, kept separate on purpose:

1. A second enumerator (`scan`) must return the same coordinates as the
   engine (`pam-skip`). A miss here is an enumeration bug.
2. A second window builder scores those coordinates. In legacy mode it has
   two real differences from the engine, both of which show up in published
   tool comparisons:
   - minus-strand windows are read off the + strand without a reverse
     complement (a strand-orientation bug)
   - every overlapping common SNV is applied at once, as if the alleles were
     known to sit on one haplotype (a phasing assumption gnomAD cannot support)

`corrected=True` rebuilds windows the same way the engine does. The agent
turns that flag on only after it has named which difference it found.
This module is not the CRISPRme+ binary. Set CRISPRME_BIN to record that an
external binary was requested; until that binary is wired, the trace says
the in-repo checker ran.
"""

from __future__ import annotations

import os

from popoff.cfd import cfd_score, revcomp
from popoff.genome import Genome
from popoff.native import enumerate_sites


def _apply(plus: str, window_pos: int, variants: list[dict], how: str) -> str:
    chars = list(plus)
    chosen = variants
    if how == "highest" and variants:
        chosen = [max(variants, key=lambda v: (v["af"], -v["pos"]))]
    for variant in chosen:
        index = variant["pos"] - window_pos
        if index < 0 or index >= len(chars):
            continue
        if chars[index] != variant["ref"]:
            continue
        if len(variant["ref"]) == 1 and len(variant["alt"]) == 1:
            chars[index] = variant["alt"]
    return "".join(chars)


def crispr_me_check(
    fasta: str,
    guide: str,
    pam_set: list[str],
    max_mm: int,
    sites: list[dict],
    genome: Genome,
    *,
    corrected: bool = False,
) -> dict:
    external = os.environ.get("CRISPRME_BIN")
    second = enumerate_sites(fasta, guide, pam_set, max_mm, algo="scan")
    engine_keys = {(s["chrom"], s["pos"], s["strand"], s["pam_class"]) for s in sites}
    second_keys = {(s["chrom"], s["pos"], s["strand"], s["pam_class"]) for s in second}
    rows = []
    for site in sites:
        plus = genome.slice(site["chrom"], site["pos"], 23)
        variants = site.get("common_variants") or []
        if corrected:
            plus_alt = _apply(plus, site["pos"], variants, "highest")
            ref_oriented = plus if site["strand"] == "+" else revcomp(plus)
            alt_oriented = plus_alt if site["strand"] == "+" else revcomp(plus_alt)
            phasing = "highest"
        else:
            plus_alt = _apply(plus, site["pos"], variants, "all")
            if site["strand"] == "-":
                ref_oriented = plus
                alt_oriented = plus_alt
            else:
                ref_oriented = plus
                alt_oriented = plus_alt
            phasing = "all-in-cis"
        ref_score = cfd_score(guide, ref_oriented) if all(b in "ACGT" for b in ref_oriented) else None
        alt_score = cfd_score(guide, alt_oriented) if all(b in "ACGT" for b in alt_oriented) else None
        rows.append({
            "chrom": site["chrom"],
            "pos": site["pos"],
            "strand": site["strand"],
            "pam_class": site["pam_class"],
            "ref_window": ref_oriented,
            "alt_window": alt_oriented,
            "cfd_ref": ref_score,
            "cfd_alt": alt_score,
            "phasing": phasing,
            "legacy_strand": (not corrected) and site["strand"] == "-",
        })
    return {
        "implementation": "external-binary-requested" if external else "in-repo-second-enumerator",
        "external_binary": external,
        "corrected": corrected,
        "engine_sites": len(engine_keys),
        "second_sites": len(second_keys),
        "only_engine": sorted(engine_keys - second_keys),
        "only_second": sorted(second_keys - engine_keys),
        "coordinates_agree": engine_keys == second_keys,
        "rows": rows,
    }
