import json
from pathlib import Path

import pytest

from popoff.panel import grid_rows, load_panel, run_panel, site_rows


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


def test_run_panel_creates_records(tiny_panel, fasta_path):
    panel = load_panel(tiny_panel)
    results = run_panel(panel, fasta=str(fasta_path))
    assert len(results) == 1
    record = results[0]["record"]
    assert record["gene"] == "HBB"
    assert "guides" in record
    assert len(record["guides"]) == 1


def test_grid_and_site_rows(tiny_panel, fasta_path):
    panel = load_panel(tiny_panel)
    results = run_panel(panel, fasta=str(fasta_path))
    grid = grid_rows(results)
    assert len(grid) == 1
    assert grid[0]["guide_id"] == "HBB-HTY3"
    sites = site_rows(results)
    assert all(isinstance(s["cfd_ref"], float) for s in sites)
