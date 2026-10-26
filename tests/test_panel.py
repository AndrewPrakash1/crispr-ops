import json
from pathlib import Path

import pytest

from popoff.panel import grid_rows, load_panel, run_panel, site_rows
from popoff.thresholds import GENE_LOCI


@pytest.fixture
def tiny_panel(tmp_path: Path) -> Path:
    panel_path = tmp_path / "panel.json"
    panel_path.write_text(
        json.dumps(
            {
                "name": "tiny",
                "populations": ["sas"],
                "panel": [
                    {
                        "gene": "HBB",
                        "guides": [
                            {"id": "HBB-HTY3", "spacer": "GTAACGGCAGACTTCTCCTC"}
                        ],
                    }
                ],
            }
        )
    )
    return panel_path


def test_load_panel(tiny_panel):
    panel = load_panel(tiny_panel)
    assert panel["panel"][0]["gene"] == "HBB"
    assert panel["panel"][0]["guides"][0]["spacer"] == "GTAACGGCAGACTTCTCCTC"


def test_demo_panel_loci_cover_all_guides():
    panel_path = Path(__file__).resolve().parents[1] / "data" / "panel.json"
    panel = load_panel(panel_path)
    assert sum(len(entry["guides"]) for entry in panel["panel"]) == 20
    assert len(panel["populations"]) == 3
    assert {
        entry["gene"]: GENE_LOCI[entry["gene"]]["chrom"]
        for entry in panel["panel"]
    } == {"HBB": "chr11", "CCR5": "chr3", "PCSK9": "chr1"}
    assert GENE_LOCI["CCR5"]["start"] <= 46371354 <= GENE_LOCI["CCR5"]["end"]
    assert GENE_LOCI["CCR5"]["start"] <= 46375211 <= GENE_LOCI["CCR5"]["end"]
    assert GENE_LOCI["PCSK9"]["start"] <= 55040729 <= GENE_LOCI["PCSK9"]["end"]
    assert GENE_LOCI["PCSK9"]["start"] <= 55061172 <= GENE_LOCI["PCSK9"]["end"]


def test_run_panel_creates_records_and_reuses_matching_run(tiny_panel, fasta_path, tmp_path, monkeypatch):
    panel = load_panel(tiny_panel)
    db_path = tmp_path / "panel.sqlite"
    monkeypatch.setenv("POPOFF_DB", str(db_path))
    results = run_panel(panel, fasta=str(fasta_path), db_path=db_path)
    assert len(results) == 1
    record = results[0]["record"]
    assert record["gene"] == "HBB"
    assert "guides" in record
    assert len(record["guides"]) == 1
    resumed = run_panel(panel, fasta=str(fasta_path), db_path=db_path)
    assert resumed[0]["run_id"] == results[0]["run_id"]


def test_grid_and_site_rows(tiny_panel, fasta_path):
    panel = load_panel(tiny_panel)
    results = run_panel(panel, fasta=str(fasta_path))
    grid = grid_rows(results)
    assert len(grid) == 1
    assert grid[0]["guide_id"] == "HBB-HTY3"
    sites = site_rows(results)
    assert all(isinstance(s["cfd_ref"], float) for s in sites)
