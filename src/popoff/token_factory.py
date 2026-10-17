"""Nebius Token Factory client for eligible Nemotron calls.

This module prefers Nebius Token Factory when ``NEBIUS_API_KEY`` and
``NEBIUS_BASE_URL`` are set.  It falls back to the legacy NVIDIA integrate
endpoint only when the Nebius variables are absent, so existing demos keep
working but the eligible path is preferred.

The model id is discovered from ``GET /v1/models`` and the first matching
Nemotron model is used.  Every completion records request id, latency, and
token usage so the trace can prove eligibility.
"""

from __future__ import annotations

import json
import os
import ssl
import time
import urllib.request
from typing import Any

import certifi

_SSL = ssl.create_default_context(cafile=certifi.where())


def _nebius_key() -> str:
    return os.environ.get("NEBIUS_API_KEY", "")


def _nebius_base() -> str:
    return os.environ.get("NEBIUS_BASE_URL", "").rstrip("/")


def _legacy_key() -> str:
    return os.environ.get("NVIDIA_API_KEY", "")


def _has_nebius() -> bool:
    return bool(_nebius_key() and _nebius_base())


def _authorization() -> str:
    return f"Bearer {_nebius_key() or _legacy_key()}"


def _base_url() -> str:
    if _has_nebius():
        return _nebius_base()
    return os.environ.get("ORCHESTRATOR_URL", "https://integrate.api.nvidia.com")


def _pick_model(models: list[dict[str, Any]]) -> str | None:
    for model in models:
        name = model.get("id") or model.get("model") or model.get("name", "")
        if name and "nemotron" in name.lower():
            return name
    return None


def list_models() -> dict[str, Any]:
    """Return eligible Nemotron model ids from the configured endpoint."""
    if not _has_nebius():
        return {"available": False, "reason": "NEBIUS_API_KEY and NEBIUS_BASE_URL are not both set"}
    url = f"{_base_url()}/v1/models"
    req = urllib.request.Request(
        url,
        headers={"Authorization": _authorization(), "Accept": "application/json"},
        method="GET",
    )
    try:
        with urllib.request.urlopen(req, timeout=30, context=_SSL) as resp:
            payload = json.loads(resp.read().decode())
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode(errors="replace")[:400]
        return {"available": False, "reason": f"HTTP {exc.code}: {detail}", "endpoint": url}
    except Exception as exc:
        return {"available": False, "reason": str(exc), "endpoint": url}
    models = payload.get("data") or payload.get("models") or []
    picked = _pick_model(models)
    return {
        "available": True,
        "endpoint": url,
        "model_ids": [m.get("id") or m.get("model") or m.get("name") for m in models],
        "nemotron_model": picked,
    }


def _model_id() -> str:
    env_model = os.environ.get("ORCHESTRATOR_MODEL")
    if env_model:
        return env_model
    listing = list_models()
    if listing.get("available") and listing.get("nemotron_model"):
        return listing["nemotron_model"]
    return "nvidia/nvidia-nemotron-nano-9b-v2"


def chat(messages: list[dict[str, str]], *, temperature: float = 0.0, max_tokens: int = 400) -> dict[str, Any]:
    """Call the configured chat/completions endpoint and return a structured result.

    The returned dict always contains:
      - ``text``: the model's message content
      - ``model``: the model id used
      - ``endpoint``: the base URL
      - ``provider``: "nebius" or "nvidia"
      - ``request_id``: upstream request id
      - ``latency_ms``: round-trip latency
      - ``usage``: token usage dict (prompt, completion, total)
      - ``error``: only present if the call failed
    """
    base = _base_url()
    if not (_nebius_key() or _legacy_key()):
        return {"error": "NEBIUS_API_KEY and NVIDIA_API_KEY are both unset"}
    url = f"{base}/v1/chat/completions"
    model = _model_id()
    body = {
        "model": model,
        "messages": messages,
        "temperature": temperature,
        "max_tokens": max_tokens,
    }
    req = urllib.request.Request(
        url,
        data=json.dumps(body).encode(),
        headers={
            "Authorization": _authorization(),
            "Content-Type": "application/json",
            "Accept": "application/json",
        },
        method="POST",
    )
    start = time.perf_counter()
    try:
        with urllib.request.urlopen(req, timeout=60, context=_SSL) as resp:
            payload = json.loads(resp.read().decode())
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode(errors="replace")[:400]
        return {
            "error": f"HTTP {exc.code}: {detail}",
            "model": model,
            "endpoint": base,
            "provider": "nebius" if _has_nebius() else "nvidia",
            "latency_ms": int((time.perf_counter() - start) * 1000),
        }
    except Exception as exc:
        return {
            "error": str(exc),
            "model": model,
            "endpoint": base,
            "provider": "nebius" if _has_nebius() else "nvidia",
            "latency_ms": int((time.perf_counter() - start) * 1000),
        }
    latency = int((time.perf_counter() - start) * 1000)
    choice = (payload.get("choices") or [{}])[0]
    text = (choice.get("message") or {}).get("content", "")
    usage = payload.get("usage") or {}
    return {
        "text": text,
        "model": payload.get("model") or model,
        "endpoint": base,
        "provider": "nebius" if _has_nebius() else "nvidia",
        "request_id": payload.get("id") or payload.get("request_id"),
        "latency_ms": latency,
        "usage": {
            "prompt_tokens": usage.get("prompt_tokens") or usage.get("prompt_tokens", 0),
            "completion_tokens": usage.get("completion_tokens") or usage.get("completion_tokens", 0),
            "total_tokens": usage.get("total_tokens") or usage.get("total_tokens", 0),
        },
    }


def eligible() -> dict[str, Any]:
    """Summary of eligibility state for the trace."""
    return {
        "provider": "nebius" if _has_nebius() else "nvidia" if _legacy_key() else None,
        "base_url": _base_url() if _has_nebius() else None,
        "has_nebius_key": bool(_nebius_key()),
        "has_legacy_key": bool(_legacy_key()),
        "model_from_env": bool(os.environ.get("ORCHESTRATOR_MODEL")),
    }
