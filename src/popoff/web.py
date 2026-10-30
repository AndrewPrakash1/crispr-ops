"""Local audit server. The trace is a server-sent stream of the agent loop."""

from __future__ import annotations

import json
from pathlib import Path

from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import FileResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles

from popoff.agent import audit_events
from popoff.bench import run_bench
from popoff.guides import GUIDES
from popoff.panel import grid_rows, site_rows
from popoff.store import get_store
from popoff.validation import evaluate

ROOT = Path(__file__).resolve().parents[2]
STATIC = ROOT / "web" / "static"
FASTA = ROOT / "data" / "genome" / "chr11.fa"
PANEL_SNAPSHOT = ROOT / "data" / "grid.json"

app = FastAPI(title="PopOff")
app.mount("/static", StaticFiles(directory=STATIC), name="static")


@app.get("/")
def index():
    return FileResponse(STATIC / "index.html")


@app.get("/api/guides")
def guides():
    return {"guides": GUIDES, "fasta": str(FASTA), "fasta_present": FASTA.exists()}


@app.get("/api/runs")
def recent_runs(limit: int = Query(default=10, ge=1, le=50)):
    return {"runs": get_store().recent(limit)}


@app.get("/api/runs/{run_id}")
def replay_run(run_id: str):
    store = get_store()
    run = store.get_run(run_id)
    record = store.get_record(run_id)
    if run is None or record is None:
        raise HTTPException(status_code=404, detail="Saved run not found.")
    return {"run": run, "events": store.replay(run_id), "record": record}


@app.get("/api/panel")
def panel_snapshot():
    if not PANEL_SNAPSHOT.exists():
        raise HTTPException(status_code=503, detail="No frozen panel snapshot is available.")
    results = json.loads(PANEL_SNAPSHOT.read_text(encoding="utf-8"))
    grid = grid_rows(results)
    sites = site_rows(results)
    return {
        "kind": "frozen_snapshot",
        "source": "data/grid.json",
        "n_runs": len(results),
        "grid": grid,
        "sites": sites,
    }


@app.post("/api/audit")
def audit(body: dict):
    request = {
        "gene": body.get("gene") or "HBB",
        "population": body.get("population") or "sas",
        "note": body.get("note") or "",
        "fasta": body.get("fasta") or str(FASTA),
        "guides": body.get("guides") or None,
    }

    def generate():
        try:
            for event in audit_events(request):
                yield f"data: {json.dumps(event, default=str)}\n\n"
        except Exception as exc:
            yield f"data: {json.dumps({'stage': 'revise', 'title': 'Run stopped', 'detail': str(exc), 'tool': None, 'payload': {}})}\n\n"

    return StreamingResponse(generate(), media_type="text/event-stream")


@app.post("/api/bench")
def bench():
    if not FASTA.exists():
        return {"error": f"missing {FASTA}"}
    report = run_bench(str(FASTA), "GTAACGGCAGACTTCTCCTC", ["NGG"], 4, threads=4)
    report["gpu_note"] = (
        "CUDA binary timed and compared to the CPU hit list."
        if report["gpu"]["available"]
        else "No CUDA binary on this machine. cuda/enumerate.cu is the kernel; build it with nvcc -O3 -o enumerate_gpu enumerate.cu."
    )
    return report


@app.get("/api/validation")
def validation():
    path = ROOT / "data" / "guideseq.tsv"
    if not path.exists():
        return {
            "available": False,
            "reason": (
                "No GUIDE-seq table is shipped, because inventing read counts would be a fake result. "
                "Place a TSV at data/guideseq.tsv with columns guide, ot_seq, reads "
                "and optional alt_seq. Rank correlation is computed from that file, and a null is reported as a null."
            ),
        }
    return {"available": True, **evaluate(path)}
