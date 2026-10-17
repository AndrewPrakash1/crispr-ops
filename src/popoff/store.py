"""SQLite run record.

Every event from an audit is stored with the run id, timestamp, and payload.
The final assembled record is also stored so the demo can replay without
spending API budget.
"""

from __future__ import annotations

import json
import os
import sqlite3
import uuid
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DB = ROOT / "data" / "popoff.sqlite"


def _db() -> Path:
    return Path(os.environ.get("POPOFF_DB", DEFAULT_DB))


def _connect() -> sqlite3.Connection:
    path = _db()
    path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(path)
    conn.row_factory = sqlite3.Row
    return conn


def init_db() -> None:
    with _connect() as conn:
        conn.executescript(
            """
            CREATE TABLE IF NOT EXISTS runs (
                run_id TEXT PRIMARY KEY,
                created_at TEXT NOT NULL,
                gene TEXT,
                population TEXT,
                provider TEXT,
                request_id TEXT,
                finished INTEGER DEFAULT 0
            );
            CREATE TABLE IF NOT EXISTS events (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                run_id TEXT NOT NULL,
                stage TEXT NOT NULL,
                title TEXT NOT NULL,
                detail TEXT,
                tool TEXT,
                payload TEXT,
                created_at TEXT NOT NULL,
                FOREIGN KEY (run_id) REFERENCES runs(run_id)
            );
            CREATE INDEX IF NOT EXISTS idx_events_run ON events(run_id);
            CREATE TABLE IF NOT EXISTS records (
                run_id TEXT PRIMARY KEY,
                record TEXT NOT NULL,
                created_at TEXT NOT NULL,
                FOREIGN KEY (run_id) REFERENCES runs(run_id)
            );
            """
        )


class RunStore:
    def __init__(self):
        init_db()

    def start(self, request: dict) -> str:
        run_id = str(uuid.uuid4())
        from popoff.token_factory import eligible

        env = eligible()
        with _connect() as conn:
            conn.execute(
                "INSERT INTO runs (run_id, created_at, gene, population, provider) VALUES (?, ?, ?, ?, ?)",
                (
                    run_id,
                    datetime.now(timezone.utc).isoformat(),
                    request.get("gene"),
                    request.get("population"),
                    env.get("provider"),
                ),
            )
        return run_id

    def add(self, run_id: str, event: dict) -> None:
        with _connect() as conn:
            conn.execute(
                "INSERT INTO events (run_id, stage, title, detail, tool, payload, created_at) VALUES (?, ?, ?, ?, ?, ?, ?)",
                (
                    run_id,
                    event.get("stage"),
                    event.get("title"),
                    event.get("detail"),
                    event.get("tool"),
                    json.dumps(event.get("payload") or {}, default=str),
                    datetime.now(timezone.utc).isoformat(),
                ),
            )

    def finish(self, run_id: str, record: dict, *, request_id: str | None = None) -> None:
        with _connect() as conn:
            conn.execute(
                "UPDATE runs SET finished = 1, request_id = ? WHERE run_id = ?",
                (request_id, run_id),
            )
            conn.execute(
                "INSERT OR REPLACE INTO records (run_id, record, created_at) VALUES (?, ?, ?)",
                (run_id, json.dumps(record, default=str), datetime.now(timezone.utc).isoformat()),
            )

    def replay(self, run_id: str) -> list[dict]:
        with _connect() as conn:
            rows = conn.execute(
                "SELECT stage, title, detail, tool, payload FROM events WHERE run_id = ? ORDER BY id",
                (run_id,),
            ).fetchall()
        events = []
        for row in rows:
            payload = json.loads(row["payload"])
            events.append({
                "stage": row["stage"],
                "title": row["title"],
                "detail": row["detail"],
                "tool": row["tool"],
                "payload": payload,
            })
        return events

    def recent(self, limit: int = 10) -> list[dict]:
        with _connect() as conn:
            rows = conn.execute(
                "SELECT run_id, created_at, gene, population, provider, finished FROM runs ORDER BY created_at DESC LIMIT ?",
                (limit,),
            ).fetchall()
        return [dict(row) for row in rows]

    def get_record(self, run_id: str) -> dict | None:
        with _connect() as conn:
            row = conn.execute("SELECT record FROM records WHERE run_id = ?", (run_id,)).fetchone()
        if row:
            return json.loads(row["record"])
        return None


_store: RunStore | None = None


def get_store() -> RunStore:
    global _store
    if _store is None:
        _store = RunStore()
    return _store


def save_events(request: dict, events: list[dict]) -> str:
    store = get_store()
    run_id = store.start(request)
    final = {}
    request_id = None
    for event in events:
        store.add(run_id, event)
        if event.get("stage") == "result":
            final = event.get("payload") or {}
            request_id = (event.get("payload") or {}).get("orchestrator", {}).get("request_id")
    store.finish(run_id, final, request_id=request_id)
    return run_id
