from fastapi.testclient import TestClient

import popoff.store as store_module
from popoff.store import RunStore
from popoff.web import app


def test_saved_run_can_be_listed_and_replayed(tmp_path, monkeypatch):
    db_path = tmp_path / "runs.sqlite"
    monkeypatch.setenv("POPOFF_DB", str(db_path))
    monkeypatch.setattr(store_module, "_store", None)
    store = RunStore()
    run_id = store.start({"gene": "HBB", "population": "sas"})
    record = {"lead": "Stored result", "guides": [], "limitations": []}
    store.add(run_id, {
        "stage": "plan",
        "title": "Saved thought",
        "detail": "Offline replay test.",
        "tool": None,
        "payload": {},
    })
    store.add(run_id, {
        "stage": "result",
        "title": "Audit complete",
        "detail": "Stored result",
        "tool": None,
        "payload": record,
    })
    store.finish(run_id, record)

    client = TestClient(app)
    listed = client.get("/api/runs")
    assert listed.status_code == 200
    assert listed.json()["runs"][0]["run_id"] == run_id

    replay = client.get(f"/api/runs/{run_id}")
    assert replay.status_code == 200
    assert replay.json()["record"] == record
    assert [event["title"] for event in replay.json()["events"]] == ["Saved thought", "Audit complete"]


def test_replay_returns_404_for_unknown_run(tmp_path, monkeypatch):
    monkeypatch.setenv("POPOFF_DB", str(tmp_path / "runs.sqlite"))
    monkeypatch.setattr(store_module, "_store", None)
    response = TestClient(app).get("/api/runs/missing")
    assert response.status_code == 404


def test_frozen_panel_snapshot_covers_full_panel():
    response = TestClient(app).get("/api/panel")
    assert response.status_code == 200
    payload = response.json()
    assert payload["kind"] == "frozen_snapshot"
    assert len(payload["grid"]) == 60
    assert len(payload["sites"]) == 9639
