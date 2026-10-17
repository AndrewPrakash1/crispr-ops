from pathlib import Path

from popoff.native import enumerate_sites


def test_pam_skip_matches_full_scan(tmp_path: Path):
    guide = "ATCGATCGATCGATCGATCG"
    # Plus-strand hit with one mismatch, and a minus-strand perfect hit.
    plus = "ATCGATCGATCGATAGATCGAGG"
    minus_oriented = guide + "TGG"
    comp = str.maketrans("ACGT", "TGCA")
    minus = minus_oriented.translate(comp)[::-1]
    seq = "N" * 30 + plus + "N" * 30 + minus + "N" * 10
    fasta = tmp_path / "mini.fa"
    fasta.write_text(">chrT\n" + seq + "\n")
    skip = enumerate_sites(str(fasta), guide, ["NGG"], 4, algo="pam-skip", threads=1)
    scan = enumerate_sites(str(fasta), guide, ["NGG"], 4, algo="scan", threads=1)
    keys = lambda rows: [(r["pos"], r["strand"], r["mismatches"], r["site23"]) for r in rows]
    assert keys(skip) == keys(scan)
    assert {row["strand"] for row in skip} == {"+", "-"}
    assert any(row["mismatches"] == 0 and row["strand"] == "-" for row in skip)
    assert any(row["mismatches"] == 1 and row["strand"] == "+" for row in skip)
