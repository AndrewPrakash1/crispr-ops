"""Evo 2 sequence-context prior.

The score is the mean per-base log-probability of a window under Evo 2, via
the NVIDIA NIM forward endpoint. The audit uses the difference between the
alternate haplotype and the reference haplotype as a second opinion. It is
not a cutting probability and it is not treated as ground truth.

Token indices are the ASCII codes of A, C, G and T, which is the Evo 2
vocabulary mapping published with the NIM.
"""

from __future__ import annotations

import json
import math
import os
import ssl
import time
import urllib.error
import urllib.request

import certifi

_SSL = ssl.create_default_context(cafile=certifi.where())

NIM_URL = os.environ.get(
    "EVO2_URL",
    "https://health.api.nvidia.com/v1/biology/arc/evo2-7b/forward",
)
_ASCII = {"A": 65, "C": 67, "G": 71, "T": 84}


def _authorization() -> str:
    nebius = os.environ.get("NEBIUS_API_KEY")
    legacy = os.environ.get("NVIDIA_API_KEY")
    key = nebius or legacy
    return f"Bearer {key}" if key else ""


def _log_softmax_at(logits_row: list[float], index: int) -> float:
    m = max(logits_row)
    log_z = m + math.log(sum(math.exp(v - m) for v in logits_row))
    return logits_row[index] - log_z


def mean_logprob_from_logits(sequence: str, logits: list[list[float]]) -> float:
    """Causal likelihood: the logit at position i predicts the base at i+1."""
    sequence = sequence.upper()
    if len(logits) < len(sequence):
        raise ValueError("logit rows are shorter than the sequence")
    total = 0.0
    n = 0
    for i, base in enumerate(sequence[1:], start=1):
        row = logits[i - 1]
        # NIM may return a 512-way vector. Guard short rows used in tests.
        if _ASCII[base] >= len(row):
            raise ValueError("logit vector is missing the DNA token index")
        total += _log_softmax_at(row, _ASCII[base])
        n += 1
    return total / max(n, 1)


def evo2_score(window: str) -> dict:
    window = "".join(base for base in window.upper() if base in "ACGT")
    if len(window) < 8:
        return {"available": False, "reason": "window shorter than 8 bp", "window_length": len(window)}
    key = os.environ.get("NEBIUS_API_KEY") or os.environ.get("NVIDIA_API_KEY")
    if not key:
        return {
            "available": False,
            "reason": "NEBIUS_API_KEY and NVIDIA_API_KEY are both unset, so the Evo 2 NIM was not called.",
            "endpoint": NIM_URL,
            "window_length": len(window),
        }
    body = json.dumps({"sequence": window, "output_layers": ["output_layer"]}).encode()
    req = urllib.request.Request(
        NIM_URL,
        data=body,
        headers={
            "content-type": "application/json",
            "authorization": _authorization(),
            "accept": "application/json",
        },
    )
    start = time.perf_counter()
    try:
        with urllib.request.urlopen(req, timeout=120, context=_SSL) as resp:
            payload = json.loads(resp.read().decode())
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode(errors="replace")[:400]
        return {"available": False, "reason": f"NIM HTTP {exc.code}: {detail}", "endpoint": NIM_URL, "latency_ms": int((time.perf_counter() - start) * 1000)}
    except Exception as exc:
        return {"available": False, "reason": str(exc), "endpoint": NIM_URL, "latency_ms": int((time.perf_counter() - start) * 1000)}
    latency_ms = int((time.perf_counter() - start) * 1000)
    logits = _extract_logits(payload)
    if logits is None:
        return {"available": False, "reason": "NIM response had no output_layer logits", "endpoint": NIM_URL, "latency_ms": latency_ms}
    return {
        "available": True,
        "mean_logprob": mean_logprob_from_logits(window, logits),
        "window_length": len(window),
        "endpoint": NIM_URL,
        "model": "arc/evo2-7b",
        "request_id": payload.get("id") or payload.get("request_id"),
        "latency_ms": latency_ms,
        "usage": payload.get("usage", {}),
    }


def _extract_logits(payload: dict) -> list[list[float]] | None:
    data = payload.get("data") or payload.get("outputs") or payload
    if isinstance(data, dict):
        layer = data.get("output_layer") or data.get("layers", {}).get("output_layer")
        if isinstance(layer, list):
            return layer
    if isinstance(payload.get("output_layer"), list):
        return payload["output_layer"]
    return None
