"""Call the C++ enumerator. The CUDA binary is used when it has been built."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CPU_SRC = ROOT / "cpp" / "enumerate.cpp"
CPU_BIN = ROOT / "cpp" / "enumerate_cpu"
GPU_BIN = ROOT / "cuda" / "enumerate_gpu"


def build_cpu() -> Path:
    CPU_BIN.parent.mkdir(parents=True, exist_ok=True)
    if CPU_BIN.exists() and CPU_BIN.stat().st_mtime >= CPU_SRC.stat().st_mtime:
        return CPU_BIN
    subprocess.run(
        ["c++", "-O3", "-std=c++17", "-pthread", str(CPU_SRC), "-o", str(CPU_BIN)],
        check=True,
    )
    return CPU_BIN


def gpu_available() -> bool:
    return GPU_BIN.exists() and os.access(GPU_BIN, os.X_OK)


def enumerate_sites(
    fasta: str,
    guide: str,
    pam_set: list[str],
    max_mm: int,
    *,
    algo: str = "pam-skip",
    threads: int | None = None,
    binary: Path | None = None,
) -> list[dict]:
    exe = binary or build_cpu()
    cmd = [
        str(exe),
        "--fasta", fasta,
        "--guide", guide.upper(),
        "--pam", ",".join(pam_set),
        "--max-mm", str(max_mm),
        "--algo", algo,
        "--threads", str(threads or os.cpu_count() or 1),
    ]
    proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    hits = []
    for line in proc.stdout.splitlines()[1:]:
        chrom, pos, strand, mm, pam_class, site23 = line.split("\t")
        hits.append({
            "chrom": chrom,
            "pos": int(pos),
            "strand": strand,
            "mismatches": int(mm),
            "pam_class": pam_class,
            "site23": site23,
            "protospacer": site23[:20],
            "pam": site23[20:],
        })
    return hits
