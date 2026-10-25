"""gnomAD v4 region queries. Allele counts are cached on disk."""

from __future__ import annotations

import json
import ssl
import threading
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import certifi

_SSL = ssl.create_default_context(cafile=certifi.where())
_REQUEST_INTERVAL_SECONDS = 6.0
_last_request_at = 0.0
_request_lock = threading.Lock()

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "data" / "cache" / "gnomad"
API = "https://gnomad.broadinstitute.org/api"
_BATCH_SIZE = 50

QUERY = """
query($chrom: String!, $start: Int!, $stop: Int!) {
  region(chrom: $chrom, start: $start, stop: $stop, reference_genome: GRCh38) {
    variants(dataset: gnomad_r4) {
      variant_id
      pos
      ref
      alt
      genome { ac an populations { id ac an } }
      exome { ac an populations { id ac an } }
    }
  }
}
"""

POPULATIONS = {
    "sas": "South Asian",
    "afr": "African / African American",
    "amr": "Admixed American",
    "eas": "East Asian",
    "nfe": "Non-Finnish European",
    "fin": "Finnish",
    "asj": "Ashkenazi Jewish",
    "mid": "Middle Eastern",
    "ami": "Amish",
    "remaining": "Remaining",
}


@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    af: float
    ac: int
    an: int
    source: str
    variant_id: str

    def as_dict(self) -> dict:
        return {
            "chrom": self.chrom,
            "pos": self.pos,
            "ref": self.ref,
            "alt": self.alt,
            "af": self.af,
            "ac": self.ac,
            "an": self.an,
            "source": self.source,
            "variant_id": self.variant_id,
        }


def _pop_af(block: dict | None, population: str) -> tuple[int, int] | None:
    if not block:
        return None
    for pop in block.get("populations") or []:
        if pop.get("id") == population and pop.get("an"):
            return int(pop["ac"]), int(pop["an"])
    return None


def _fetch_region(chrom: str, start: int, stop: int) -> dict:
    CACHE.mkdir(parents=True, exist_ok=True)
    key = _cache_path(chrom, start, stop)
    if key.exists():
        return json.loads(key.read_text())
    body = json.dumps({
        "query": QUERY,
        "variables": {"chrom": chrom, "start": start, "stop": stop},
    }).encode()
    req = urllib.request.Request(
        API, data=body, headers={"content-type": "application/json", "accept": "application/json"}
    )
    payload = _request_payload(req)
    key.write_text(json.dumps(payload))
    return payload


def prefetch_af(sites: list[dict]) -> int:
    """Fetch uncached 23-base site intervals in batched GraphQL requests."""
    CACHE.mkdir(parents=True, exist_ok=True)
    intervals = {}
    for site in sites:
        chrom = site["chrom"]
        api_chrom = chrom[3:] if chrom.startswith("chr") else chrom
        if not api_chrom.isdecimal():
            raise ValueError(f"unexpected chromosome name {chrom!r}")
        start = int(site["pos"])
        stop = start + 22
        intervals[(api_chrom, start, stop)] = _cache_path(api_chrom, start, stop)
    pending = [(interval, path) for interval, path in intervals.items() if not path.exists()]
    for offset in range(0, len(pending), _BATCH_SIZE):
        _prefetch_batch(pending[offset : offset + _BATCH_SIZE])
    return len(pending)


def _prefetch_batch(batch: list[tuple[tuple[str, int, int], Path]]) -> None:
    fields = []
    for index, ((chrom, start, stop), _path) in enumerate(batch):
        fields.append(
            f'  r{index}: region(chrom: "{chrom}", start: {start}, stop: {stop}, '
            "reference_genome: GRCh38) { variants(dataset: gnomad_r4) { "
            "variant_id pos ref alt genome { ac an populations { id ac an } } "
            "exome { ac an populations { id ac an } } } }"
        )
    query = "query {\n" + "\n".join(fields) + "\n}"
    req = urllib.request.Request(
        API,
        data=json.dumps({"query": query}).encode(),
        headers={"content-type": "application/json", "accept": "application/json"},
    )
    try:
        payload = _request_payload(req)
    except urllib.error.HTTPError as exc:
        if exc.code not in (400, 413) or len(batch) == 1:
            raise
        midpoint = len(batch) // 2
        _prefetch_batch(batch[:midpoint])
        _prefetch_batch(batch[midpoint:])
        return
    data = payload.get("data") or {}
    for index, (_interval, path) in enumerate(batch):
        path.write_text(json.dumps({"data": {"region": data.get(f"r{index}")}}))


def _cache_path(chrom: str, start: int, stop: int) -> Path:
    return CACHE / f"{chrom}_{start}_{stop}.json"


def _request_payload(req: urllib.request.Request) -> dict:
    payload = None
    for attempt in range(4):
        _wait_for_request_slot()
        try:
            with urllib.request.urlopen(req, timeout=60, context=_SSL) as resp:
                payload = json.loads(resp.read().decode())
            break
        except urllib.error.HTTPError as exc:
            if exc.code not in (429, 500, 502, 503) or attempt == 3:
                raise
            retry_after = exc.headers.get("Retry-After")
            delay = float(retry_after) if retry_after and retry_after.isdecimal() else 0.0
            time.sleep(max(delay, 1.5 * (attempt + 1)))
        except urllib.error.URLError:
            if attempt == 3:
                raise
            time.sleep(1.5 * (attempt + 1))
    if payload is None:
        raise RuntimeError("gnomAD request failed")
    if payload.get("errors"):
        raise RuntimeError(payload["errors"])
    return payload


def _wait_for_request_slot() -> None:
    global _last_request_at
    with _request_lock:
        elapsed = time.monotonic() - _last_request_at
        if elapsed < _REQUEST_INTERVAL_SECONDS:
            time.sleep(_REQUEST_INTERVAL_SECONDS - elapsed)
        _last_request_at = time.monotonic()


def fetch_af(chrom: str, start: int, end: int, population: str) -> list[Variant]:
    """Variants overlapping 1-based inclusive [start, end] for one population.

    Genome calls are preferred. Exome calls are used when the variant has no
    genome genotype counts. gnomAD chromosome names have no 'chr' prefix.
    """
    population = population.lower()
    if population not in POPULATIONS:
        raise ValueError(f"unknown population {population}")
    api_chrom = chrom[3:] if chrom.startswith("chr") else chrom
    payload = _fetch_region(api_chrom, start, end)
    region = (payload.get("data") or {}).get("region") or {}
    found: list[Variant] = []
    for row in region.get("variants") or []:
        pos = int(row["pos"])
        if pos < start or pos > end:
            continue
        ref, alt = row["ref"], row["alt"]
        if len(ref) != 1 or len(alt) != 1:
            continue
        counts = _pop_af(row.get("genome"), population)
        source = "genome"
        if counts is None:
            counts = _pop_af(row.get("exome"), population)
            source = "exome"
        if counts is None:
            continue
        ac, an = counts
        found.append(Variant(
            chrom=chrom if chrom.startswith("chr") else f"chr{chrom}",
            pos=pos,
            ref=ref,
            alt=alt,
            af=ac / an,
            ac=ac,
            an=an,
            source=source,
            variant_id=row.get("variant_id") or f"{api_chrom}-{pos}-{ref}-{alt}",
        ))
    found.sort(key=lambda v: (-v.af, v.pos))
    return found


def network_ok() -> bool:
    try:
        urllib.request.urlopen(API, timeout=5, context=_SSL)
        return True
    except urllib.error.HTTPError:
        return True
    except Exception:
        return False
