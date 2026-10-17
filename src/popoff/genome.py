"""FASTA reference, 1-based genomic coordinates."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field


@dataclass
class Genome:
    path: str
    sequences: dict[str, str]
    sha256: str

    def slice(self, chrom: str, pos: int, length: int) -> str:
        """Return `length` bases starting at 1-based `pos`."""
        seq = self.sequences[chrom]
        start = pos - 1
        if start < 0 or start + length > len(seq):
            raise ValueError(f"{chrom}:{pos} length {length} is outside the contig")
        return seq[start : start + length]

    def base(self, chrom: str, pos: int) -> str:
        return self.sequences[chrom][pos - 1]


def load_fasta(path: str) -> Genome:
    sequences: dict[str, str] = {}
    name = ""
    chunks: list[str] = []
    digest = hashlib.sha256()

    def flush() -> None:
        nonlocal chunks, name
        if name:
            sequences[name] = "".join(chunks)
        chunks = []

    with open(path, "rb") as handle:
        for line in handle:
            digest.update(line)
            if line.startswith(b">"):
                flush()
                header = line[1:].decode().strip()
                name = header.split()[0]
            else:
                chunks.append(line.decode().strip().upper())
        flush()
    if not sequences:
        raise ValueError(f"no sequences in {path}")
    return Genome(path=path, sequences=sequences, sha256=digest.hexdigest())


@dataclass
class Loaded:
    """Process-wide cache so the web app does not reparse chr11 per request."""

    by_path: dict[str, Genome] = field(default_factory=dict)

    def get(self, path: str) -> Genome:
        if path not in self.by_path:
            self.by_path[path] = load_fasta(path)
        return self.by_path[path]


GENOMES = Loaded()
