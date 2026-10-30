from pathlib import Path

import pytest


@pytest.fixture
def root() -> Path:
    return Path(__file__).resolve().parents[1]


@pytest.fixture(autouse=True)
def isolate_run_store(tmp_path: Path, monkeypatch):
    import popoff.store as store_module

    monkeypatch.setenv("POPOFF_DB", str(tmp_path / "test-runs.sqlite"))
    monkeypatch.setattr(store_module, "_store", None)


@pytest.fixture
def fasta_path(tmp_path: Path) -> Path:
    path = tmp_path / "test-reference.fa"
    path.write_text(">chr11\n" + ("N" * 6_000_000) + "\n")
    return path
