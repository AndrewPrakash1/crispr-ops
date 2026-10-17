"""Run a multi-guide, multi-population audit panel and export an assumption grid.

The panel is a JSON file with a list of genes, each carrying its own guide
spacers.  Every gene x population combination is audited independently, the
final records are stored in the SQLite run store, and a CSV assumption grid is
written so the table can be inspected without the API server.
"""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

from popoff.agent import audit_events
from popoff.store import get_store


def load_panel(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def run_panel(panel: dict, fasta: str | None = None, db_path: str | Path | None = None) -> list[dict]:
    """Run every gene x population combination in *panel*.

    Returns a list of result dicts, each containing the gene, population,
    run_id, and the assembled audit record.
    """
    if db_path:
        os.environ["POPOFF_DB"] = str(db_path)

    populations = panel.get("populations", ["sas", "afr", "nfe"])
    default_fasta = panel.get("fasta")
    fasta = fasta or default_fasta or ""

    results: list[dict] = []
    for entry in panel.get("panel", []):
        gene = entry["gene"]
        guides = entry.get("guides", [])
        if not guides:
            continue
        for population in populations:
            request = {
                "gene": gene,
                "population": population,
                "note": entry.get("note", ""),
                "fasta": fasta,
                "guides": [dict(g) for g in guides],
            }
            events = list(audit_events(request))
            run_id = _last_run_id()
            record = {}
            if events and events[-1].get("stage") == "result":
                record = events[-1].get("payload") or {}
            results.append({
                "gene": gene,
                "population": population,
                "run_id": run_id,
                "record": record,
            })
    return results


def _last_run_id() -> str | None:
    recent = get_store().recent(1)
    if recent:
        return recent[0]["run_id"]
    return None


def grid_rows(results: list[dict]) -> list[dict]:
    """Flatten panel results into assumption-grid rows."""
    rows: list[dict] = []
    for result in results:
        record = result.get("record") or {}
        for guide in record.get("guides", []):
            rows.append({
                "gene": result["gene"],
                "population": result["population"],
                "guide_id": guide["id"],
                "spacer": guide["spacer"],
                "reference_verdict": guide["reference_verdict"],
                "population_verdict": guide["population_verdict"],
                "flipped": guide["flipped"],
                "escalated": guide["escalated"],
                "n_sites": guide["n_sites"],
            })
    return rows


def site_rows(results: list[dict]) -> list[dict]:
    """Flatten per-site rows from all panel results."""
    rows: list[dict] = []
    for result in results:
        record = result.get("record") or {}
        for site in record.get("sites", []):
            if site.get("on_target"):
                continue
            evo = site.get("evo2") or {}
            variant = site.get("applied_variant") or {}
            rows.append({
                "gene": result["gene"],
                "population": result["population"],
                "guide_id": site.get("guide_id"),
                "chrom": site.get("chrom"),
                "pos": site.get("pos"),
                "strand": site.get("strand"),
                "pam_class": site.get("pam_class"),
                "mm_ref": site.get("mm_ref"),
                "mm_alt": site.get("mm_alt"),
                "cfd_ref": site.get("cfd_ref"),
                "cfd_alt": site.get("cfd_alt"),
                "af": variant.get("af"),
                "variant": f"{variant.get('ref', '')}>{variant.get('alt', '')}" if variant else "",
                "evo2_direction": evo.get("direction"),
            })
    return rows


def write_csv(rows: list[dict], path: str | Path) -> None:
    if not rows:
        return
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(f=handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def export_grid(results: list[dict], grid_path: str | Path, site_path: str | Path | None = None) -> None:
    """Write the assumption grid and, optionally, the per-site detail CSVs."""
    write_csv(grid_rows(results), grid_path)
    if site_path:
        write_csv(site_rows(results), site_path)
