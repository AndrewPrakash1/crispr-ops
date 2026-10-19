#!/usr/bin/env python3
"""Design candidate SpCas9 guide spacers inside GRCh38 gene regions.

Output is a JSON list ready to paste into data/panel.json. The script avoids
random sequences by deriving every spacer from an NGG (or CCN) PAM site in the
provided FASTA.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from popoff.genome import load_fasta  # noqa: E402


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGT", "TGCA")
    return seq.translate(table)[::-1]


def gc_content(seq: str) -> float:
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq)


def has_homopolymer(seq: str, max_run: int = 4) -> bool:
    for base in "ACGT":
        if base * (max_run + 1) in seq:
            return True
    return False


def low_complexity(seq: str) -> bool:
    """Reject obvious low-complexity spacers."""
    if has_homopolymer(seq, 4):
        return True
    # Reject dinucleotide repeats that dominate the spacer.
    for i in range(2, 6):
        pattern = seq[:i]
        if len(pattern) == i and (pattern * (20 // i))[:20] == seq:
            return True
    # Reject longer simple repeats (e.g. GTTGTTGTTGTT...)
    for k in [3, 4]:
        for start in range(k):
            sub = seq[start::k]
            if has_homopolymer(sub, 4):
                return True
    return False


def design_guides(fasta_path: str, chrom: str, start: int, end: int, n: int = 10) -> list[dict]:
    genome = load_fasta(fasta_path)
    seq = genome.sequences[chrom]
    candidates = []
    seen = set()
    for i in range(max(0, start - 23), min(end + 23, len(seq) - 22)):
        # Plus strand NGG PAM (positions 21 and 22 of the 23-mer window)
        if seq[i + 21] == "G" and seq[i + 22] == "G":
            spacer = seq[i : i + 20]
            if "N" in spacer:
                continue
            if spacer not in seen:
                seen.add(spacer)
                candidates.append({
                    "spacer": spacer,
                    "chrom": chrom,
                    "pos": i + 1,
                    "strand": "+",
                    "pam": seq[i + 20 : i + 23],
                    "gc": gc_content(spacer),
                })
        # Minus strand NGG PAM -> plus strand CCN
        if seq[i] == "C" and seq[i + 1] == "C":
            protospacer = seq[i + 3 : i + 23]
            if "N" in protospacer:
                continue
            spacer = revcomp(protospacer)
            if spacer not in seen:
                seen.add(spacer)
                candidates.append({
                    "spacer": spacer,
                    "chrom": chrom,
                    "pos": i + 1,
                    "strand": "-",
                    "pam": revcomp(seq[i : i + 3]),
                    "gc": gc_content(spacer),
                })

    # Quality filter
    filtered = []
    for c in candidates:
        if c["gc"] < 0.25 or c["gc"] > 0.85:
            continue
        if low_complexity(c["spacer"]):
            continue
        # Avoid spacers ending in stretches of G/C that can mis-fold.
        if c["spacer"].endswith("GGGG") or c["spacer"].endswith("CCCC"):
            continue
        filtered.append(c)

    # Diversify by position across the locus.
    filtered.sort(key=lambda c: c["pos"])
    if not filtered:
        return []
    if len(filtered) <= n:
        return filtered
    # Pick n guides at evenly spaced indices.
    step = len(filtered) // n
    chosen = []
    for i in range(n):
        idx = min(i * step + step // 2, len(filtered) - 1)
        if i > 0 and filtered[idx]["spacer"] == filtered[chosen[-1]]["spacer"]:
            continue
        chosen.append(idx)
    return [filtered[i] for i in chosen]


def main() -> int:
    regions = {
        "CCR5": ("chr3", 46370142, 46376206),
        "PCSK9": ("chr1", 55039548, 55064853),
    }
    results = {}
    for gene, (chrom, start, end) in regions.items():
        fasta = ROOT / "data" / "genome" / f"{chrom}.fa"
        guides = design_guides(str(fasta), chrom, start, end, n=10)
        results[gene] = guides

    print(json.dumps(results, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
