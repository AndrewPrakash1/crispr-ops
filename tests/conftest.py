from pathlib import Path

import pytest


@pytest.fixture
def root() -> Path:
    return Path(__file__).resolve().parents[1]


@pytest.fixture
def fasta_path(root: Path) -> Path:
    return root / "data" / "genome" / "chr11.fa"
