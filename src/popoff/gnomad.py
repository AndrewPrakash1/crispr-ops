"""gnomAD v4 region queries. Allele counts are cached on disk."""

from __future__ import annotations

import json
import ssl
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import certifi

_SSL = ssl.create_default_context(cafile=certifi.where())

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "data" / "cache" / "gnomad"
API = "https://gnomad.broadinstitute.org/api"

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
    key = CACHE / f"{chrom}_{start}_{stop}.json"
    if key.exists():
        return json.loads(key.read_text())
    # gnomAD GraphQL is rate-limited; a small pause keeps long panels from 429-ing.
    time.sleep(0.3)
    body = json.dumps({
        "query": QUERY,
        "variables": {"chrom": chrom, "start": start, "stop": stop},
    }).encode()
    req = urllib.request.Request(
        API, data=body, headers={"content-type": "application/json", "accept": "application/json"}
    )
    payload = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=60, context=_SSL) as resp:
                payload = json.loads(resp.read().decode())
            break
        except urllib.error.HTTPError as exc:
            if exc.code not in (429, 500, 502, 503) or attempt == 3:
                raise
            time.sleep(1.5 * (attempt + 1))
        except urllib.error.URLError:
            if attempt == 3:
                raise
            time.sleep(1.5 * (attempt + 1))
    if payload is None:
        raise RuntimeError("gnomAD request failed")
    if payload.get("errors"):
        raise RuntimeError(payload["errors"])
    key.write_text(json.dumps(payload))
    return payload


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
