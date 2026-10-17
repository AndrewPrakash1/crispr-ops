"""Choose the next legal tool.

Nemotron is the orchestrator when a Token Factory or NVIDIA API key is set. It
may return only a tool from the legal set handed to it. Any other tool is
rejected and the registered policy choice is kept. A thought that adds a
number the tools did not return is dropped, so the trace cannot grow an allele
frequency the run did not measure.

When the key is unset, the model is not called.
"""

from __future__ import annotations

import json
import os
import re

from popoff.policy import review
from popoff import token_factory

MODEL = os.environ.get("ORCHESTRATOR_MODEL") or "nvidia/nvidia-nemotron-nano-9b-v2"
_NUMBER = re.compile(r"\d+(?:\.\d+)?")
_ALLOWED = {"0.2", "0.05", "0.1", "0.3", "0.5", "4", "1", "20", "23"}


def choose(*, legal: list[str], policy: str, observation: str, caller=None) -> dict:
    if policy not in legal:
        raise ValueError(f"policy tool {policy} is not in {legal}")
    env = token_factory.eligible()
    force_policy = os.environ.get("POPOFF_NO_MODEL")
    if force_policy or not (env.get("has_nebius_key") or env.get("has_legacy_key")):
        reason = "POPOFF_NO_MODEL is set." if force_policy else "No eligible API key is set."
        return {
            "tool": policy,
            "thought": (
                f"{reason} The orchestrator model was not called. "
                f"The registered policy chose {policy}."
            ),
            "source": "policy",
            "model": None,
            "accepted": True,
            "provider": None,
        }
    caller = caller or _nim_chat
    try:
        raw = caller(_messages(legal, observation))
    except Exception:
        return _policy(policy, "The orchestrator call failed. The registered policy chose the next tool.")
    parsed = _parse(raw)
    meta = last_chat_meta()
    if parsed is None:
        return _policy(policy, "The orchestrator reply was not a tool choice. The registered policy chose the next tool.")
    if meta:
        parsed["request_id"] = meta.get("request_id")
        parsed["latency_ms"] = meta.get("latency_ms")
        parsed["usage"] = meta.get("usage")
        parsed["provider"] = meta.get("provider")
    tool = parsed.get("tool")
    thought = str(parsed.get("thought") or "").strip()
    if tool not in legal:
        return _policy(
            policy,
            f"The model named {tool}, which is not a legal next tool. The registered policy kept {policy}.",
        )
    if review(thought):
        thought = f"The model chose {tool}. Its explanation was refused by policy and was not kept."
    else:
        extra = [n for n in _NUMBER.findall(thought) if n not in _ALLOWED and n not in _NUMBER.findall(observation)]
        if extra:
            thought = (
                f"The model chose {tool}. Its explanation contained a number that is not in the "
                "tool observation, so the explanation was dropped."
            )
    return {
        "tool": tool,
        "thought": thought or f"The model chose {tool}.",
        "source": "nim",
        "model": MODEL,
        "accepted": True,
        "provider": env.get("provider"),
        "request_id": parsed.get("request_id"),
        "latency_ms": parsed.get("latency_ms"),
        "usage": parsed.get("usage"),
    }


def banner() -> dict:
    env = token_factory.eligible()
    if env.get("has_nebius_key") or env.get("has_legacy_key"):
        thought = (
            f"Orchestrator is {MODEL} via {env.get('provider', 'unknown')}. "
            "It chooses among legal tools. A choice outside that set is rejected."
        )
        source = "nim"
    else:
        thought = (
            "No eligible API key is set. The orchestrator model was not called. "
            "Branch decisions use the registered policy until a key is set. "
            "The trace will not fill in allele frequencies, CFD scores, or a GUIDE-seq correlation."
        )
        source = "policy"
    return {"tool": None, "thought": thought, "source": source, "model": MODEL if source == "nim" else None, "accepted": True, "provider": env.get("provider")}


def _policy(policy: str, thought: str) -> dict:
    return {"tool": policy, "thought": thought, "source": "policy", "model": MODEL, "accepted": False, "provider": None}


def _messages(legal: list[str], observation: str) -> list[dict]:
    return [
        {
            "role": "system",
            "content": (
                "You are the orchestrator of a CRISPR off-target audit. "
                "Reply with one JSON object and nothing else: "
                '{"tool": "<one of the legal tools>", "thought": "<one sentence>"}. '
                "Use only numbers that appear in the observation. "
                "Do not call a guide safe or unsafe. Do not diagnose. "
                "Do not invent an allele frequency, a CFD score, or an experimental correlation. "
                f"Legal tools: {', '.join(legal)}."
            ),
        },
        {"role": "user", "content": observation},
    ]


def _parse(raw: str) -> dict | None:
    text = re.sub(r"\s*<think>.*?</think>\s*", "", raw or "", flags=re.DOTALL)
    start = text.find("{")
    end = text.rfind("}")
    if start < 0 or end <= start:
        return None
    try:
        parsed = json.loads(text[start : end + 1])
    except json.JSONDecodeError:
        return None
    if not isinstance(parsed, dict) or "tool" not in parsed:
        return None
    return parsed


def _nim_chat(messages: list[dict]) -> str:
    result = token_factory.chat(messages, temperature=0.0, max_tokens=400)
    if "error" in result:
        raise RuntimeError(result["error"])
    # Attach metadata to the module-level result cache so callers can log it.
    _nim_chat._last_meta = {  # type: ignore[attr-defined]
        "model": result.get("model"),
        "request_id": result.get("request_id"),
        "latency_ms": result.get("latency_ms"),
        "usage": result.get("usage"),
        "provider": result.get("provider"),
    }
    return result["text"]


_nim_chat._last_meta = {}  # type: ignore[attr-defined]


def last_chat_meta() -> dict:
    """Return metadata from the most recent Nemotron chat call."""
    return getattr(_nim_chat, "_last_meta", {})
