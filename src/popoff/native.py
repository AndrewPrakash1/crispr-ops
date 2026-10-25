"""Call the C++ enumerator. The CUDA binary is used when it has been built."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CPU_SRC = ROOT / "cpp" / "enumerate.cpp"
CPU_BIN = ROOT / "cpp" / "enumerate_cpu"
GPU_BIN = ROOT / "cuda" / "enumerate_gpu"


def build_cpu() -> Path:
    binary = CPU_BIN.with_suffix(".exe") if os.name == "nt" else CPU_BIN
    CPU_BIN.parent.mkdir(parents=True, exist_ok=True)
    if binary.exists() and binary.stat().st_mtime >= CPU_SRC.stat().st_mtime:
        return binary

    compiler = shutil.which("c++")
    if compiler:
        subprocess.run(
            [compiler, "-O3", "-std=c++17", "-pthread", str(CPU_SRC), "-o", str(binary)],
            check=True,
        )
    elif os.name == "nt":
        _build_cpu_msvc(binary)
    else:
        subprocess.run(
            ["c++", "-O3", "-std=c++17", "-pthread", str(CPU_SRC), "-o", str(binary)],
            check=True,
        )
    return binary


def _build_cpu_msvc(binary: Path) -> None:
    program_files_x86 = os.environ.get("ProgramFiles(x86)", r"C:\Program Files (x86)")
    vswhere = Path(program_files_x86) / "Microsoft Visual Studio" / "Installer" / "vswhere.exe"
    installations = subprocess.run(
        [
            str(vswhere),
            "-latest",
            "-products",
            "*",
            "-requires",
            "Microsoft.VisualStudio.Component.VC.Tools.x86.x64",
            "-property",
            "installationPath",
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip().splitlines()
    if not installations:
        raise FileNotFoundError("No Visual Studio C++ Build Tools installation was found.")
    vcvars = Path(installations[-1]) / "VC" / "Auxiliary" / "Build" / "vcvars64.bat"
    command = (
        f'& $env:ComSpec /c \'call "{vcvars}" >nul && cl /nologo /O2 /std:c++17 /EHsc '
        f'/Fe:"{binary}" "{CPU_SRC}"\''
    )
    subprocess.run(
        [shutil.which("powershell.exe") or "powershell.exe", "-NoProfile", "-Command", command],
        check=True,
        cwd=CPU_SRC.parent,
    )


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
