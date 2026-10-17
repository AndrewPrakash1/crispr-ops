"""Default HBB panel.

HBB-HTY3 is a published exon-1 spacer (the reference allele; the sickle site is
the penultimate base). The other two are NGG sites read off GRCh38 at this
locus. They are candidates for the audit, not therapeutic recommendations.
"""

from __future__ import annotations

GUIDES = {
    "HBB": [
        {
            "id": "HBB-HTY3",
            "spacer": "GTAACGGCAGACTTCTCCTC",
            "note": "Published HBB exon-1 spacer against the reference allele. GRCh38 hit is chr11:5226984 with PAM AGG.",
        },
        {
            "id": "HBB-EX1-A",
            "spacer": "CTTGCCCCACAGGGCAGTAA",
            "note": "Plus-strand NGG site in HBB exon 1 on GRCh38 (PAM CGG). Locus-derived candidate.",
        },
        {
            "id": "HBB-PROM",
            "spacer": "AGGGCTGGGCATAAAAGTCA",
            "note": "Minus-strand NGG site at the HBB promoter on GRCh38 (window chr11:5227090, PAM GGG). Locus-derived candidate.",
        },
    ]
}


def panel_for(gene: str, guides: list[dict] | None) -> list[dict]:
    if guides:
        cleaned = []
        for guide in guides:
            spacer = guide["spacer"].upper().replace("U", "T")
            if len(spacer) == 23:
                spacer = spacer[:20]
            if len(spacer) != 20 or any(base not in "ACGT" for base in spacer):
                raise ValueError(f"guide {guide.get('id')} must be a 20 nt ACGT spacer")
            cleaned.append({"id": guide.get("id") or spacer, "spacer": spacer, "note": guide.get("note", "")})
        return cleaned
    gene = gene.upper()
    if gene not in GUIDES:
        raise ValueError(f"no default panel for {gene}; pass guide spacers")
    return [dict(item) for item in GUIDES[gene]]
