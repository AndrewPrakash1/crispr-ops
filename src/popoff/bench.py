"""CPU enumeration benchmark. The CUDA binary is timed when it has been built."""

from __future__ import annotations

import time
from pathlib import Path

from popoff.native import CPU_BIN, GPU_BIN, build_cpu, enumerate_sites, gpu_available

ROOT = Path(__file__).resolve().parents[2]


def run_bench(fasta: str, guide: str, pam: list[str] | None = None, max_mm: int = 4, threads: int = 4) -> dict:
    pam = pam or ["NGG"]
    build_cpu()
    t0 = time.perf_counter()
    cpu_hits = enumerate_sites(fasta, guide, pam, max_mm, algo="pam-skip", threads=threads)
    cpu_s = time.perf_counter() - t0
    t1 = time.perf_counter()
    scan_hits = enumerate_sites(fasta, guide, pam, max_mm, algo="scan", threads=threads)
    scan_s = time.perf_counter() - t1
    cpu_keys = [(h["chrom"], h["pos"], h["strand"], h["mismatches"], h["site23"], h["pam_class"]) for h in cpu_hits]
    scan_keys = [(h["chrom"], h["pos"], h["strand"], h["mismatches"], h["site23"], h["pam_class"]) for h in scan_hits]
    gpu = {"available": gpu_available(), "binary": str(GPU_BIN), "equal": None, "seconds": None, "hits": None}
    if gpu_available():
        t2 = time.perf_counter()
        gpu_hits = enumerate_sites(fasta, guide, pam, max_mm, algo="pam-skip", threads=1, binary=GPU_BIN)
        gpu["seconds"] = time.perf_counter() - t2
        gpu_keys = [(h["chrom"], h["pos"], h["strand"], h["mismatches"], h["site23"], h["pam_class"]) for h in gpu_hits]
        gpu["equal"] = gpu_keys == cpu_keys
        gpu["hits"] = len(gpu_hits)
    bases = _bases(fasta)
    return {
        "fasta": fasta,
        "guide": guide,
        "pam": pam,
        "max_mm": max_mm,
        "threads": threads,
        "bases": bases,
        "cpu_binary": str(CPU_BIN),
        "cpu_seconds": cpu_s,
        "scan_seconds": scan_s,
        "hits": len(cpu_hits),
        "algorithms_equal": cpu_keys == scan_keys,
        "guides_per_sec": 1 / cpu_s if cpu_s else None,
        "bases_per_sec": bases / cpu_s if cpu_s else None,
        "gpu": gpu,
    }


def _bases(fasta: str) -> int:
    n = 0
    with open(fasta) as handle:
        for line in handle:
            if not line.startswith(">"):
                n += len(line.strip())
    return n
