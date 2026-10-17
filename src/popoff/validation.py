"""Ranking check against an experimental off-target table, when one is supplied.

The shipped default does not invent GUIDE-seq read counts. Put a TSV at
data/guideseq.tsv with columns guide, ot_seq, reads (and optionally af, alt_seq)
to score it. A null result is reported as a null.
"""

from __future__ import annotations

from pathlib import Path

from popoff.cfd import cfd_score, hamming


def _rank(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        average = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranks[order[k]] = average
        i = j + 1
    return ranks


def spearman(left: list[float], right: list[float]) -> float | None:
    if len(left) < 3 or len(left) != len(right):
        return None
    a, b = _rank(left), _rank(right)
    mean_a = sum(a) / len(a)
    mean_b = sum(b) / len(b)
    num = sum((x - mean_a) * (y - mean_b) for x, y in zip(a, b))
    den_a = sum((x - mean_a) ** 2 for x in a) ** 0.5
    den_b = sum((y - mean_b) ** 2 for y in b) ** 0.5
    if den_a == 0 or den_b == 0:
        return None
    return num / (den_a * den_b)


def load_table(path: Path) -> list[dict]:
    lines = [line for line in path.read_text().splitlines() if line.strip() and not line.startswith("#")]
    header = [cell.strip() for cell in lines[0].split("\t")]
    rows = []
    for line in lines[1:]:
        cells = line.split("\t")
        rows.append({key: cells[i] if i < len(cells) else "" for i, key in enumerate(header)})
    return rows


def evaluate(path: Path) -> dict:
    rows = load_table(path)
    guides = []
    reads = []
    cfd = []
    mm_score = []
    pop = []
    pop_reads = []
    for row in rows:
        guide = row["guide"].upper()[:20]
        ot = row["ot_seq"].upper()
        if len(ot) < 23:
            continue
        guides.append(guide)
        count = float(row["reads"])
        reads.append(count)
        cfd.append(cfd_score(guide, ot[:23]))
        mm_score.append(-float(hamming(guide, ot[:20])))
        if row.get("alt_seq"):
            pop.append(cfd_score(guide, row["alt_seq"].upper()[:23]))
            pop_reads.append(count)
    cfd_rho = spearman(cfd, reads)
    mm_rho = spearman(mm_score, reads)
    pop_rho = spearman(pop, pop_reads) if len(pop) >= 3 else None
    if cfd_rho is None or mm_rho is None:
        comparison = "Too few rows to compare rankings."
    elif abs(cfd_rho - mm_rho) < 0.02:
        comparison = "Null: CFD and mismatch count rank the experimental reads about equally."
    elif cfd_rho > mm_rho:
        comparison = "CFD ranks the experimental reads better than mismatch count."
    else:
        comparison = "Mismatch count ranks the experimental reads better than CFD. Reported as a null for the CFD claim."
    population = "No population-haplotype sequences in the table."
    if pop_rho is not None and cfd_rho is not None:
        if abs(pop_rho - cfd_rho) < 0.02:
            population = "Null: population-adjusted CFD does not improve rank agreement with these reads."
        elif pop_rho > cfd_rho:
            population = "Population-adjusted CFD agrees with these reads better than reference CFD."
        else:
            population = "Population-adjusted CFD agrees with these reads worse than reference CFD."
    return {
        "n": len(reads),
        "n_population": len(pop),
        "spearman_cfd": cfd_rho,
        "spearman_mismatch": mm_rho,
        "spearman_population": pop_rho,
        "comparison": comparison,
        "population": population,
        "source": str(path),
    }
