#!/usr/bin/env python3
"""Run the same audit with and without the orchestrator model and produce a diff.

When NEBIUS_API_KEY+NEBIUS_BASE_URL are set, the "model" branch calls Nemotron.
The "policy" branch forces the registered policy by unsetting the model keys in a
subprocess. The resulting traces are compared on:
- Evo 2 spend (number of sites sent to Evo 2)
- first hypothesis chosen during a checker disagreement
- overall site/hypothesis order
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PYTHON = ROOT / ".venv" / "bin" / "python"
if sys.platform == "win32":
    PYTHON = PYTHON.with_suffix(".exe")


def _run(*, policy: bool = False) -> list[dict]:
    cmd = [
        str(PYTHON),
        "-c",
        "import sys; sys.path.insert(0, 'src'); "
        "from popoff.agent import audit_events; "
        "import json; "
        "events = list(audit_events({'gene':'HBB','population':'sas','fasta':'data/genome/chr11.fa'})); "
        "print(json.dumps(events, default=str))",
    ]
    run_env = {**os.environ}
    if policy:
        run_env["POPOFF_NO_MODEL"] = "1"
        run_env.pop("NEBIUS_API_KEY", None)
        run_env.pop("NEBIUS_BASE_URL", None)
        run_env.pop("NVIDIA_API_KEY", None)
    else:
        run_env.pop("POPOFF_NO_MODEL", None)
    result = subprocess.run(cmd, cwd=ROOT, env=run_env, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(result.stderr)
    return json.loads(result.stdout)


def _summarize(events: list[dict]) -> dict:
    evo_spent = 0
    evo_deferred = 0
    first_hypothesis = None
    for event in events:
        if event.get("tool") == "evo2_score":
            evo_spent += 1
        if event.get("tool") == "defer_evo":
            evo_deferred += 1
        payload = event.get("payload") or {}
        if first_hypothesis is None and event.get("tool") == "crispr_me_check" and event.get("stage") == "revise":
            hyp = (payload.get("hypothesis") or {}).get("tool")
            if hyp:
                first_hypothesis = hyp
    return {
        "evo_spent": evo_spent,
        "evo_deferred": evo_deferred,
        "first_hypothesis": first_hypothesis,
        "n_events": len(events),
    }


def main() -> int:
    print("Running policy branch ...")
    policy_events = _run(policy=True)
    policy_summary = _summarize(policy_events)

    has_keys = bool(os.environ.get("NEBIUS_API_KEY") and os.environ.get("NEBIUS_BASE_URL"))
    if has_keys:
        print("Running model branch ...")
        model_events = _run(policy=False)
        model_summary = _summarize(model_events)
    else:
        print("No NEBIUS_API_KEY+NEBIUS_BASE_URL set; model branch skipped.")
        model_events = []
        model_summary = {"skipped": True}

    diff = {
        "model_keys_present": has_keys,
        "policy": policy_summary,
        "model": model_summary,
        "differ": (
            policy_summary.get("evo_spent") != model_summary.get("evo_spent") or
            policy_summary.get("first_hypothesis") != model_summary.get("first_hypothesis")
        ) if not model_summary.get("skipped") else False,
    }

    out = ROOT / "data" / "ab_diff.json"
    out.write_text(json.dumps(diff, indent=2, default=str))
    print("Wrote", out)
    print(json.dumps(diff, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
