#!/usr/bin/env python3
"""Prove eligibility: list Nemotron models and make one chat + one tool call.

Run with NEBIUS_API_KEY and NEBIUS_BASE_URL exported. The script writes a JSON
artifact containing the model id, request id, latency, and token usage. If the
credentials are missing, it exits with a clear message instead of inventing data.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from popoff import orchestrator, token_factory  # noqa: E402


def main() -> int:
    env = token_factory.eligible()
    if not env.get("has_nebius") and not env.get("has_legacy_key"):
        print("Set NEBIUS_API_KEY+NEBIUS_BASE_URL (Nebius track) or NVIDIA_API_KEY+ORCHESTRATOR_URL.")
        return 1
    if env.get("has_nebius"):
        print("Nebius track path (Token Factory / Studio).")
    else:
        print(
            "Dev orchestrator path via",
            env.get("provider"),
            "— official hackathon eligibility still requires Nebius Token Factory.",
        )

    print("Listing models ...")
    listing = token_factory.list_models()
    if not listing.get("available"):
        print("Model listing failed:", listing.get("reason"))
        return 1

    model_id = listing.get("nemotron_model")
    if not model_id:
        print("No Nemotron model found in /v1/models response.")
        print("Models returned:", listing.get("model_ids"))
        return 1

    print("Making one chat completion with", model_id, "...")
    chat_result = token_factory.chat(
        messages=[
            {"role": "system", "content": "You are a helpful assistant. Reply with a single word."},
            {"role": "user", "content": "Hello"},
        ],
        temperature=0.0,
        max_tokens=20,
    )
    if "error" in chat_result:
        print("Chat failed:", chat_result.get("error"))
        return 1

    print("Making one tool/orchestrator call ...")
    tool_result = orchestrator.choose(
        legal=["strand", "window", "phasing"],
        policy="strand",
        observation="Cross-check disagrees. Minus-strand windows may be reverse-complemented.",
    )

    artifact = {
        "nebius_track_eligible": env.get("nebius_track_eligible"),
        "listing": listing,
        "chat": {
            "model": chat_result.get("model"),
            "request_id": chat_result.get("request_id"),
            "latency_ms": chat_result.get("latency_ms"),
            "usage": chat_result.get("usage"),
            "provider": chat_result.get("provider"),
            "endpoint": chat_result.get("endpoint"),
        },
        "tool": {
            "tool": tool_result.get("tool"),
            "source": tool_result.get("source"),
            "model": tool_result.get("model"),
            "provider": tool_result.get("provider"),
            "request_id": tool_result.get("request_id"),
            "latency_ms": tool_result.get("latency_ms"),
            "usage": tool_result.get("usage"),
            "accepted": tool_result.get("accepted"),
        },
    }

    out = ROOT / "data" / "eligibility_trace.json"
    out.write_text(json.dumps(artifact, indent=2, default=str))
    print("Wrote", out)
    print(json.dumps(artifact, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
