import json
import urllib.error
from io import BytesIO

from popoff import gnomad


class _Response:
    def __init__(self, payload):
        self.payload = json.dumps(payload).encode()

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return None

    def read(self):
        return self.payload


def test_prefetch_batches_site_queries_and_caches_population_data(tmp_path, monkeypatch):
    monkeypatch.setattr(gnomad, "CACHE", tmp_path)
    monkeypatch.setattr(gnomad, "_BATCH_SIZE", 50)
    monkeypatch.setattr(gnomad, "_wait_for_request_slot", lambda: None)
    calls = []

    def urlopen(request, **_kwargs):
        query = json.loads(request.data)["query"]
        alias_count = query.count(": region(")
        calls.append(alias_count)
        if alias_count > 2:
            raise urllib.error.HTTPError(request.full_url, 400, "query too large", {}, BytesIO(b""))
        aliases = [f"r{index}" for index in range(alias_count)]
        data = {alias: {"variants": []} for alias in aliases}
        if len(calls) == 2:
            data["r0"]["variants"] = [
                {
                    "variant_id": "1-101-A-G",
                    "pos": 101,
                    "ref": "A",
                    "alt": "G",
                    "genome": {"ac": 8, "an": 100, "populations": [{"id": "sas", "ac": 4, "an": 20}]},
                    "exome": None,
                }
            ]
        return _Response({"data": data})

    monkeypatch.setattr(gnomad.urllib.request, "urlopen", urlopen)
    sites = [
        {"chrom": "chr1", "pos": 100},
        {"chrom": "chr1", "pos": 200},
        {"chrom": "chr3", "pos": 300},
    ]

    assert gnomad.prefetch_af(sites) == 3
    assert calls == [3, 1, 2]
    assert gnomad.fetch_af("chr1", 100, 122, "sas")[0].af == 0.2
    assert gnomad.fetch_af("chr1", 200, 222, "sas") == []
    assert gnomad.fetch_af("chr3", 300, 322, "sas") == []
    assert gnomad.prefetch_af(sites) == 0
    assert len(calls) == 3
