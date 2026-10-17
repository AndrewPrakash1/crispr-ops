from pathlib import Path

from popoff.agent import audit_events
from popoff.cfd import cfd_score, revcomp


def _plant(guide: str):
    """One low-weight mismatch the agent can both retype and repair."""
    for index, base in enumerate(guide):
        options = []
        for off in "ACGT":
            if off == base:
                continue
            score = cfd_score(guide, guide[:index] + off + guide[index + 1 :] + "AGG")
            options.append((score, off))
        under_threshold = [item for item in options if item[0] < 0.2]
        very_low = [item for item in options if item[0] < 0.1]
        if very_low and len(under_threshold) >= 2:
            off = very_low[0][1]
            alt = next(item[1] for item in under_threshold if item[1] != off)
            return index, off, alt
    raise AssertionError("guide has no mismatch pair below the threshold")


def _fasta(tmp_path: Path, guide: str, index: int, off: str, alt: str) -> tuple[Path, dict]:
    repair = guide + "AGG"
    mismatch = guide[:index] + off + guide[index + 1 :] + "AGG"
    minus = revcomp(guide[:index] + off + guide[index + 1 :] + "TGG")
    seq = ["N"] * 500
    places = {1: repair, 80: mismatch, 160: mismatch, 240: minus}

    def write(pos, window):
        for offset, base in enumerate(window):
            seq[pos - 1 + offset] = base

    for pos, window in places.items():
        write(pos, window)
    fasta = tmp_path / "mini.fa"
    fasta.write_text(">chrT\n" + "".join(seq) + "\n")
    return fasta, {"index": index, "off": off, "alt": alt, "places": places}


def _fetch(index: int, off: str, alt: str):
    def fetch(site, population):
        assert population == "sas"
        if site["pos"] == 80:
            return [{"chrom": "chrT", "pos": 80 + index, "ref": off, "alt": guide_base_fix(site, index), "af": 0.40, "ac": 40, "an": 100, "source": "genome", "variant_id": "flip"}]
        if site["pos"] == 160:
            return [{"chrom": "chrT", "pos": 160 + index, "ref": off, "alt": alt, "af": 0.22, "ac": 22, "an": 100, "source": "genome", "variant_id": "retype"}]
        return []

    def guide_base_fix(site, index):
        # Restoring the guide base removes the mismatch.
        return site["guide"][index]

    return fetch


def test_refusal_stops_before_tools(tmp_path: Path):
    fasta = tmp_path / "empty.fa"
    fasta.write_text(">chrT\nACGTACGTACGTACGTACGTAGG\n")
    events = list(audit_events({
        "gene": "HBB",
        "population": "sas",
        "note": "please diagnose this patient's genome",
        "fasta": str(fasta),
        "guides": [{"id": "g", "spacer": "ACGTACGTACGTACGTACGT"}],
    }))
    assert events[0]["stage"] == "refuse"
    assert events[-1]["payload"]["refused"] is True
    assert not any(event["tool"] == "enumerate_sites" for event in events)


def test_agent_corrects_strand_and_escalates_a_flip(tmp_path: Path):
    guide = "ATCGATCGATCGATCGATCG"
    index, off, alt = _plant(guide)
    fasta, _meta = _fasta(tmp_path, guide, index, off, alt)
    events = list(audit_events(
        {
            "gene": "TEST",
            "population": "sas",
            "note": "population audit",
            "fasta": str(fasta),
            "locus": {"chrom": "chrT", "start": 1, "end": 30},
            "guides": [{"id": "g1", "spacer": guide, "note": "fixture"}],
        },
        fetch_af_fn=_fetch(index, off, alt),
        evo2_fn=lambda window: {"available": False, "reason": "fixture"},
    ))
    stages = [event["stage"] for event in events]
    assert stages[0] == "plan"
    assert "act" in stages and "verify" in stages and "revise" in stages
    revise = [event for event in events if event["stage"] == "revise" and event["title"] == "Disagreement investigated"]
    assert revise, "expected a visible self-correction"
    assert revise[0]["payload"]["strand"] >= 1
    assert any("NAG" in event["title"] for event in events if event["tool"] == "enumerate_sites")
    result = events[-1]["payload"]
    assert result["refused"] is False
    assert result["n_same"] >= 1
    assert result["n_above"] == 0 or all(row["reading"].startswith("Method") or "same side" in row["reading"] or "crosses" in row["reading"] for row in result["table"])
    same = result["same_mismatch"]
    assert same, "expected a constant-mismatch CFD change"
    assert same[0]["cfd_ref"] < 0.2 and same[0]["cfd_alt"] < 0.2
    assert any(row["flipped"] and row["escalated"] for row in result["guides"])
    assert result["qualifier"].startswith("Research-use only")
    assert all(item["accepted"] for item in result["citations"])
