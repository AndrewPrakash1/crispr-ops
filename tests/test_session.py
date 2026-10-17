from pathlib import Path

from popoff.nat_tools import AuditSession


def _fasta(tmp_path: Path) -> str:
    path = tmp_path / "t.fa"
    path.write_text(">chr11\n" + ("N" * 80) + "\n")
    return str(path)


def test_session_refuses_a_clinical_question(tmp_path: Path):
    session = AuditSession(_fasta(tmp_path), "HBB", "sas", note="please diagnose this patient's genome")
    assert session.refused
    assert session.enumerate("HBB-HTY3")["refused"] is True


def test_gpu_enumerate_is_blocked_until_a_flip_or_unresolved_disagreement(tmp_path: Path):
    session = AuditSession(_fasta(tmp_path), "HBB", "sas", note="population audit")
    guide = session.guides[0]
    guide.phase = "expand"
    result = session.gpu_enumerate(guide.spec["id"])
    assert "error" in result
    assert "not legal" in result["error"]
