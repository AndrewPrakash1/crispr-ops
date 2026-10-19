#!/usr/bin/env python3
"""Verify that every spacer in data/panel.json maps back to its chromosome."""

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


def parse_note(note: str) -> tuple[str, int, str] | None:
    # Designed-guide format: "GRCh38 chr1:55043961 minus strand, PAM TGG"
    m = re.search(r"(chr\d+):(\d+)\s+(minus|plus)\s+strand", note)
    if m:
        strand = "+" if m.group(3) == "plus" else "-"
        return m.group(1), int(m.group(2)), strand
    # Legacy HBB format: "... chr11:5226984 with PAM AGG"
    m = re.search(r"(chr\d+):(\d+)\s+with PAM (\w+)", note)
    if m:
        pam = m.group(3)
        strand = "+" if pam in {"AGG", "CGG", "TGG", "GGG", "GGA", "GGC", "GGT"} else "-"
        return m.group(1), int(m.group(2)), strand
    return None


def main() -> int:
    panel = json.loads((ROOT / "data" / "panel.json").read_text())
    loaded = {}
    errors = []
    total = 0
    for entry in panel["panel"]:
        gene = entry["gene"]
        for guide in entry["guides"]:
            total += 1
            parsed = parse_note(guide["note"])
            if parsed is None:
                print(f"SKIP {gene}/{guide['id']}: note not parseable")
                continue
            chrom, pos, strand = parsed
            if chrom not in loaded:
                loaded[chrom] = load_fasta(str(ROOT / "data" / "genome" / f"{chrom}.fa"))
            genome = loaded[chrom]
            spacer = guide["spacer"].upper()
            window = genome.slice(chrom, pos, 23)
            if strand == "+":
                expected = spacer
                protospacer = window[:20]
                pam = window[20:23]
            else:
                expected = revcomp(spacer)
                protospacer = window[3:23]
                pam = revcomp(window[:3])
            ok = protospacer == expected and len(pam) == 3 and pam[1:] == "GG"
            status = "OK" if ok else "FAIL"
            print(f"{status} {gene}/{guide['id']}: {spacer} {chrom}:{pos}{strand} {pam}")
            if not ok:
                errors.append(guide)
    print(f"\n{total - len(errors)}/{total} verified")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
