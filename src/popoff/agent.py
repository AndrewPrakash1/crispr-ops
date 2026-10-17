"""Plan, act, verify, revise.

The next tool is chosen from what the previous tool returned. That is the
difference between this loop and a script that always runs the same stages:

- triage sends Evo 2 only to sites whose CFD delta is large or whose score
  sits near the registered threshold
- a CFD disagreement is investigated (window, strand, phasing) and the check
  is run again
- a guide whose verdict flips is re-enumerated at NAG and NGA, and only that
  guide
- a personal genome or a clinical question is refused before any tool runs
"""

from __future__ import annotations

import hashlib
import os
import platform
import sys
from collections.abc import Callable, Iterator

from popoff import thresholds as T
from popoff.cfd import cfd_score, hamming, one_snp_envelope, pair_weights, revcomp
from popoff.citations import verify_citation
from popoff.crisprme import crispr_me_check
from popoff.evo2 import evo2_score
from popoff.genome import GENOMES, Genome
from popoff.gnomad import POPULATIONS, fetch_af, network_ok
from popoff.guides import panel_for
from popoff.native import GPU_BIN, enumerate_sites, gpu_available
from popoff.orchestrator import banner, choose
from popoff.policy import review
from popoff.resolver import analyze, resolve_sites

QUALIFIER = (
    "Research-use only. This is an audit of how off-target scores behave under "
    "population variation. It is not a clinical screen, not a diagnosis, and it "
    "does not accept personal genomes. A CFD change below the registered "
    "threshold of 0.2 is a statement about the scoring method, not a statement "
    "that a site is safe or unsafe. Evo 2 likelihood is a sequence-context prior, "
    "not a cutting probability."
)


def _event(stage: str, title: str, detail: str = "", tool: str | None = None, payload: dict | None = None) -> dict:
    return {"stage": stage, "title": title, "detail": detail, "tool": tool, "payload": payload or {}}


def _orient(plus: str, strand: str) -> str:
    return plus if strand == "+" else revcomp(plus)


def _apply(plus: str, start: int, variants: list[dict]) -> str:
    chars = list(plus)
    for variant in variants:
        index = variant["pos"] - start
        if 0 <= index < len(chars) and chars[index] == variant["ref"] and len(variant["alt"]) == 1:
            chars[index] = variant["alt"]
    return "".join(chars)


def _near(score: float | None) -> bool:
    return score is not None and abs(score - T.ACTIONABLE_CFD) <= 0.05


class _Guide:
    def __init__(self, spec: dict):
        self.spec = spec
        self.waves = [("NGG",)]
        self.wave_i = 0
        self.phase = "enumerate"
        self.sites: list[dict] = []
        self.cursor = 0
        self.escalated = False
        self.force_expand = False
        self.resolution: dict | None = None
        self.resolved_waves: set[str] = set()
        self.triaged_waves: set[str] = set()
        self.evo_queue: list[dict] = []

    @property
    def wave(self) -> tuple[str, ...]:
        return self.waves[self.wave_i]

    @property
    def wave_name(self) -> str:
        return ",".join(self.wave)

    def wave_sites(self) -> list[dict]:
        return [site for site in self.sites if site["wave"] == self.wave_name]


def audit_events(
    request: dict,
    *,
    fetch_af_fn: Callable | None = None,
    evo2_fn: Callable | None = None,
    store_run: bool = True,
) -> Iterator[dict]:
    fetch_af_fn = fetch_af_fn or _default_fetch
    evo2_fn = evo2_fn or evo2_score

    def _gen():
        nonlocal fetch_af_fn, evo2_fn
        yield from _audit_events_inner(request, fetch_af_fn=fetch_af_fn, evo2_fn=evo2_fn)

    events = []
    for event in _gen():
        events.append(event)
        yield event
    if store_run:
        from popoff.store import save_events
        save_events(request, events)


def _audit_events_inner(
    request: dict,
    *,
    fetch_af_fn: Callable | None = None,
    evo2_fn: Callable | None = None,
) -> Iterator[dict]:
    fetch_af_fn = fetch_af_fn or _default_fetch
    evo2_fn = evo2_fn or evo2_score
    text = " ".join([
        request.get("gene", ""),
        request.get("population", ""),
        request.get("note", ""),
        " ".join(g.get("spacer", "") for g in request.get("guides") or []),
    ])
    refusal = review(text)
    if refusal:
        yield _event("refuse", "Policy refusal", refusal["reason"], "policy", refusal)
        yield _event("result", "Refused", refusal["reason"], payload={"refused": True, "refusal": refusal, "qualifier": QUALIFIER})
        return

    population = request.get("population", "sas").lower()
    if population not in POPULATIONS:
        yield _event("refuse", "Unknown population", f"{population} is not a gnomAD label.", payload={"refused": True})
        return
    try:
        guides = [_Guide(spec) for spec in panel_for(request["gene"], request.get("guides"))]
    except ValueError as exc:
        yield _event("refuse", "Guides rejected", str(exc), payload={"refused": True, "qualifier": QUALIFIER})
        return

    genome = GENOMES.get(request["fasta"])
    locus = request.get("locus") or T.GENE_LOCI.get(request["gene"].upper())
    env = _environment(genome)
    yield _event(
        "plan",
        "Initial plan",
        f"Audit {request['gene']} in {POPULATIONS[population]}. "
        f"Enumerate NGG with at most {T.MAX_MISMATCHES} mismatches, score CFD, "
        f"then decide. Threshold τ = {T.ACTIONABLE_CFD} was registered before this run. "
        f"Common variant means allele frequency ≥ {T.COMMON_AF}.",
        payload={"environment": env, "guides": [g.spec for g in guides]},
    )
    opening = banner()
    yield _event("plan", "Orchestrator", opening["thought"], payload={"orchestrator": opening})

    cited = False
    claims: list[str] = []
    citation_cursor = 0
    result: dict = {}

    while True:
        action = _next(guides, cited, claims, citation_cursor, result)
        if action is None:
            break
        action, branch = _branch(action)
        if branch is not None:
            yield _event("plan", "Orchestrator", branch["thought"], action.get("tool"), {"orchestrator": branch})
        if action.get("deferred"):
            continue
        why = action.get("why", "")
        if branch is not None and branch.get("tool") == action.get("tool"):
            why = f"{branch['thought']} {why}".strip()
        yield _event("plan", action["plan"], why, action["tool"], action.get("plan_payload"))
        if action["tool"] == "finish":
            result = _assemble(request, population, guides, genome, env, claims)
            claims = result.pop("_claims")
            citation_cursor = 0
            cited = False
            yield _event("verify", "Run record is ready for citation checks", f"{len(claims)} claims.", payload={"n_claims": len(claims)})
            continue
        if action["tool"] == "verify_citation":
            verdict = verify_citation(claims[citation_cursor], result["record"])
            yield _event("act", "verify_citation", claims[citation_cursor], "verify_citation", {"claim": claims[citation_cursor]})
            yield _event(
                "verify" if verdict["accepted"] else "revise",
                "Citation accepted" if verdict["accepted"] else "Citation rejected",
                verdict["reason"],
                "verify_citation",
                verdict,
            )
            result.setdefault("citations", []).append(verdict)
            citation_cursor += 1
            if citation_cursor >= len(claims):
                cited = True
                result["qualifier"] = QUALIFIER
                yield _event("result", "Audit complete", result["lead"], payload=result)
            continue

        observations = _act(action, request, population, genome, locus, fetch_af_fn, evo2_fn)
        yield _event("act", action["title"], action.get("detail", ""), action["tool"], observations["public"])
        verification = _verify(action, observations, genome)
        yield _event("verify", verification["title"], verification["detail"], action["tool"], verification["payload"])
        revision = _revise(action, observations, verification, genome)
        if revision:
            yield _event("revise", revision["title"], revision["detail"], action["tool"], revision["payload"])


def _next(guides: list[_Guide], cited: bool, claims: list[str], citation_cursor: int, result: dict) -> dict | None:
    for guide in guides:
        if guide.phase == "done":
            continue
        spec = guide.spec
        if guide.phase == "enumerate":
            return {
                "tool": "enumerate_sites",
                "guide": guide,
                "plan": f"Enumerate {spec['id']} at PAM {guide.wave_name}",
                "why": "This guide has no sites yet for the current PAM set.",
                "title": f"enumerate_sites({spec['id']}, {guide.wave_name})",
                "detail": f"max_mm={T.MAX_MISMATCHES}",
            }
        pending_fetch = [s for s in guide.wave_sites() if not s.get("fetched")]
        if pending_fetch:
            site = pending_fetch[0]
            return {
                "tool": "fetch_af",
                "guide": guide,
                "site": site,
                "plan": f"fetch_af {site['chrom']}:{site['pos']}{site['strand']}",
                "why": "Allele frequency is still missing, so the haplotype is unknown.",
                "title": "fetch_af",
                "detail": f"{site['chrom']}:{site['pos']}{site['strand']} mm={site['mismatches']}",
            }
        pending_score = [s for s in guide.wave_sites() if s.get("needs_score") and not s.get("scored")]
        if pending_score:
            site = pending_score[0]
            return {
                "tool": "score_cfd",
                "guide": guide,
                "site": site,
                "plan": f"score_cfd at {site['chrom']}:{site['pos']}{site['strand']}",
                "why": "A common SNV overlaps this window. Reference and alternate haplotypes get separate CFD scores.",
                "title": "score_cfd",
                "detail": site.get("applied_variant", {}).get("variant_id", ""),
            }
        if guide.phase == "check":
            return {
                "tool": "crispr_me_check",
                "guide": guide,
                "corrected": False,
                "plan": f"Cross-check {spec['id']} {guide.wave_name}",
                "why": "A second enumerator and a second window builder have not looked at this wave.",
                "title": "crispr_me_check",
                "detail": "legacy window builder",
            }
        if guide.phase == "recheck":
            return {
                "tool": "crispr_me_check",
                "guide": guide,
                "corrected": True,
                "plan": f"Re-run the cross-check for {spec['id']}",
                "why": "The investigation named a fix. The check runs again with that fix.",
                "title": "crispr_me_check",
                "detail": "strand-corrected, highest-AF haplotype",
            }
        if guide.phase == "evo":
            if guide.wave_name not in guide.triaged_waves:
                _triage(guide)
            if guide.evo_queue:
                site = guide.evo_queue[0]
                return {
                    "tool": "evo2_score",
                    "guide": guide,
                    "site": site,
                    "plan": f"Evo 2 second opinion at {site['chrom']}:{site['pos']}",
                    "why": site.get("triage_reason", "Selected by triage."),
                    "title": "evo2_score",
                    "detail": "ref window and alt window",
                    "legal": ["evo2_score", "defer_evo"],
                }
            guide.phase = "expand"
        if guide.phase == "expand":
            if (_flipped(guide) or guide.force_expand) and not guide.escalated and guide.wave == ("NGG",):
                guide.escalated = True
                guide.waves.append(T.EXPANDED_PAM)
                guide.wave_i += 1
                guide.phase = "enumerate"
                return {
                    "tool": "enumerate_sites",
                    "guide": guide,
                    "plan": f"Escalate {spec['id']} to NAG and NGA",
                    "why": "The guide's verdict flipped, or the corrected checker still disagrees. "
                           "The PAM search widens for this guide only.",
                    "title": f"enumerate_sites({spec['id']}, NAG,NGA)",
                    "detail": "escalation",
                    "legal": ["enumerate_sites", "hold"],
                }
            guide.phase = "done"
    if not result:
        return {"tool": "finish", "plan": "Close the loop and check citations", "why": "Every guide has a verdict.", "title": "assemble"}
    if not cited and citation_cursor < len(claims):
        return {
            "tool": "verify_citation",
            "plan": "Check the next claim against the run record",
            "why": "A number that is not in the record, or a clinical conclusion, is not displayed.",
            "title": "verify_citation",
        }
    return None


def _branch(action: dict) -> tuple[dict, dict | None]:
    """Ask the orchestrator when more than one tool is legal."""
    legal = action.get("legal") or []
    if len(legal) < 2:
        return action, None
    observation = f"{action.get('plan', '')}. {action.get('why', '')} Legal tools: {', '.join(legal)}."
    decision = choose(legal=legal, policy=action["tool"], observation=observation)
    if decision["tool"] == "defer_evo":
        guide: _Guide = action["guide"]
        guide.evo_queue.clear()
        for site in guide.wave_sites():
            if not site.get("evo2"):
                site["evo2"] = {"available": False, "reason": "Orchestrator deferred Evo 2 for this wave."}
        guide.phase = "expand"
        return {**action, "deferred": True}, decision
    if decision["tool"] == "hold":
        decision = {
            **decision,
            "accepted": False,
            "tool": action["tool"],
            "thought": decision["thought"] + " The registered rule still expands NAG and NGA for this guide.",
        }
    return action, decision


def _triage(guide: _Guide) -> None:
    guide.triaged_waves.add(guide.wave_name)
    ranked = []
    for site in guide.wave_sites():
        if site.get("on_target"):
            continue
        delta = None
        if site.get("cfd_alt") is not None:
            delta = abs(site["cfd_alt"] - site["cfd_ref"])
        near = _near(site.get("cfd_ref")) or _near(site.get("cfd_alt"))
        if (delta is not None and delta >= 0.05) or near:
            ranked.append((0 if near else 1, -(delta or 0), site))
            site["triage_reason"] = (
                f"|ΔCFD|={delta:.4f}" if delta is not None else "near threshold"
            ) + (" and within 0.05 of τ" if near else "")
        else:
            site["evo2"] = {"available": False, "reason": "Deferred. |ΔCFD| < 0.05 and the score is not within 0.05 of τ."}
        site["_triaged"] = True
    ranked.sort(key=lambda item: (item[0], item[1]))
    guide.evo_queue = [site for _, _, site in ranked[: T.EVO2_BUDGET]]
    skipped = max(0, len(ranked) - T.EVO2_BUDGET)
    for _, _, site in ranked[T.EVO2_BUDGET :]:
        site["evo2"] = {"available": False, "reason": f"Over the Evo 2 budget of {T.EVO2_BUDGET}."}
    guide.triage_note = f"{len(guide.evo_queue)} sites selected, {skipped} over budget, {sum(1 for s in guide.sites if s.get('evo2', {}).get('reason', '').startswith('Deferred'))} deferred."


def _act(action, request, population, genome: Genome, locus, fetch_af_fn, evo2_fn) -> dict:
    guide: _Guide = action["guide"]
    if action["tool"] == "enumerate_sites":
        device = "cpu"
        binary = None
        gpu_reason = ""
        if guide.wave != ("NGG",):
            if gpu_available():
                device = "gpu"
                binary = GPU_BIN
            else:
                gpu_reason = (
                    "cuda/enumerate_gpu is not built on this machine. "
                    "The expanded PAM search used the CPU enumerator and the same hit contract."
                )
        hits = enumerate_sites(
            request["fasta"], guide.spec["spacer"], list(guide.wave), T.MAX_MISMATCHES,
            algo="pam-skip", binary=binary,
        )
        for hit in hits:
            hit.update({
                "guide_id": guide.spec["id"],
                "guide": guide.spec["spacer"],
                "wave": guide.wave_name,
                "fetched": False,
                "on_target": bool(
                    hit["mismatches"] == 0 and locus and hit["chrom"] == locus["chrom"]
                    and locus["start"] <= hit["pos"] <= locus["end"]
                ),
            })
            lo, hi = one_snp_envelope(guide.spec["spacer"], hit["site23"])
            hit["envelope"] = [lo, hi]
        guide.sites.extend(hits)
        guide.phase = "fetch" if hits else "check"
        return {
            "public": {
                "n": len(hits),
                "pam": guide.wave_name,
                "on_target": sum(1 for h in hits if h["on_target"]),
                "device": device,
                "gpu_reason": gpu_reason,
            },
            "hits": hits,
        }

    if action["tool"] == "fetch_af":
        site = action["site"]
        variants = fetch_af_fn(site, population)
        conflicts = []
        applicable = []
        for variant in variants:
            if genome.base(site["chrom"], variant["pos"]) != variant["ref"]:
                conflicts.append(variant)
                continue
            applicable.append(variant)
        common = [v for v in applicable if v["af"] >= T.COMMON_AF]
        site["variants"] = applicable
        site["conflicts"] = conflicts
        site["common_variants"] = common
        site["fetched"] = True
        plus = genome.slice(site["chrom"], site["pos"], 23)
        site["ref_window"] = _orient(plus, site["strand"])
        site["cfd_ref"] = cfd_score(guide.spec["spacer"], site["ref_window"])
        site["mm_ref"] = hamming(guide.spec["spacer"], site["ref_window"][:20])
        if common:
            top = max(common, key=lambda v: (v["af"], -v["pos"]))
            site["applied_variant"] = top
            site["needs_score"] = True
        else:
            site["scored"] = True
            site["cfd_alt"] = None
        if all(s.get("fetched") for s in guide.wave_sites()) and not any(s.get("needs_score") and not s.get("scored") for s in guide.wave_sites()):
            guide.phase = "check"
        return {"public": {"n_variants": len(applicable), "n_common": len(common), "conflicts": len(conflicts), "pos": site["pos"]}, "site": site}

    if action["tool"] == "score_cfd":
        site = action["site"]
        plus = genome.slice(site["chrom"], site["pos"], 23)
        alt_plus = _apply(plus, site["pos"], [site["applied_variant"]])
        alt = _orient(alt_plus, site["strand"])
        scored = {
            "cfd_ref": site["cfd_ref"],
            "cfd_alt": cfd_score(guide.spec["spacer"], alt),
            "mm_ref": site["mm_ref"],
            "mm_alt": hamming(guide.spec["spacer"], alt[:20]),
            "pairs": pair_weights(guide.spec["spacer"], site["ref_window"], alt),
        }
        site["alt_window"] = alt
        site.update(scored)
        site["scored"] = True
        if not any(s.get("needs_score") and not s.get("scored") for s in guide.wave_sites()):
            guide.phase = "check"
        return {"public": scored | {"pos": site["pos"], "strand": site["strand"]}, "site": site}

    if action["tool"] == "crispr_me_check":
        report = crispr_me_check(
            request["fasta"], guide.spec["spacer"], list(guide.wave), T.MAX_MISMATCHES,
            guide.wave_sites(), genome, corrected=action["corrected"],
        )
        key = guide.wave_name + ("#fixed" if action["corrected"] else "")
        guide.resolved_waves.add(key)
        if action["corrected"]:
            guide.phase = "evo"
        return {"public": {k: report[k] for k in ("coordinates_agree", "engine_sites", "second_sites", "corrected", "implementation")}, "report": report, "corrected": action["corrected"]}

    if action["tool"] == "evo2_score":
        site = guide.evo_queue.pop(0) if guide.evo_queue else action["site"]
        ref_ctx, alt_ctx = _contexts(genome, site)
        ref = evo2_fn(ref_ctx)
        alt = evo2_fn(alt_ctx) if site.get("alt_window") else {"available": False, "reason": "no alternate haplotype"}
        opinion = {"ref": ref, "alt": alt}
        if ref.get("available") and alt.get("available") and site.get("cfd_alt") is not None:
            evo_delta = alt["mean_logprob"] - ref["mean_logprob"]
            cfd_delta = site["cfd_alt"] - site["cfd_ref"]
            opinion["evo_delta"] = evo_delta
            opinion["cfd_delta"] = cfd_delta
            opinion["direction"] = "same direction" if evo_delta * cfd_delta > 0 else "disagreement"
            opinion["reading"] = "Second opinion only. This is not a cutting probability."
        site["evo2"] = opinion
        if not guide.evo_queue:
            guide.phase = "expand"
        return {"public": {"pos": site["pos"], "opinion": {k: opinion[k] for k in opinion if k in ("direction", "reading", "evo_delta", "cfd_delta")} | {"ref_available": ref.get("available"), "alt_available": alt.get("available"), "reason": ref.get("reason")}}, "site": site}

    raise RuntimeError(action["tool"])


def _contexts(genome: Genome, site: dict) -> tuple[str, str]:
    length = 256
    start = max(1, site["pos"] - 116)
    chrom_len = len(genome.sequences[site["chrom"]])
    if start + length - 1 > chrom_len:
        start = max(1, chrom_len - length + 1)
    plus = genome.slice(site["chrom"], start, min(length, chrom_len - start + 1))
    alt = plus
    if site.get("applied_variant"):
        alt = _apply(plus, start, [site["applied_variant"]])
    return plus, alt


def _verify(action, observations, genome: Genome) -> dict:
    if action["tool"] == "enumerate_sites":
        n = observations["public"]["n"]
        device = observations["public"].get("device", "cpu")
        detail = f"Hit list is sorted and restricted to the requested PAM set. Device: {device}."
        if observations["public"].get("gpu_reason"):
            detail += " " + observations["public"]["gpu_reason"]
        return {"title": f"{n} sites", "detail": detail, "payload": observations["public"]}
    if action["tool"] == "fetch_af":
        site = observations["site"]
        mismatch = site["mm_ref"] != site["mismatches"]
        detail = "Reference base matches the FASTA." if not site["conflicts"] else f"{len(site['conflicts'])} alleles did not match the FASTA and were dropped."
        if mismatch:
            detail += " Mismatch count does not match the enumerator."
        return {"title": "Allele frequencies checked", "detail": detail, "payload": {"mm_ok": not mismatch, "conflicts": site["conflicts"]}}
    if action["tool"] == "score_cfd":
        site = observations["site"]
        same = site["mm_ref"] == site["mm_alt"]
        changed = abs(site["cfd_ref"] - site["cfd_alt"]) > 1e-12
        kind = "mismatch type changed, count did not" if same and changed else "mismatch count changed" if not same else "no CFD change"
        return {"title": "CFD pair scored", "detail": f"ref {site['cfd_ref']:.4f} → alt {site['cfd_alt']:.4f} ({kind}).", "payload": {"kind": kind, "cfd_ref": site["cfd_ref"], "cfd_alt": site["cfd_alt"]}}
    if action["tool"] == "crispr_me_check":
        report = observations["report"]
        n_bad = 0
        if report["coordinates_agree"]:
            by_key = {(s["chrom"], s["pos"], s["strand"], s["pam_class"]): s for s in action["guide"].wave_sites()}
            for row in report["rows"]:
                site = by_key[(row["chrom"], row["pos"], row["strand"], row["pam_class"])]
                if row["cfd_ref"] is None or abs(row["cfd_ref"] - site["cfd_ref"]) > 1e-9:
                    n_bad += 1
                elif site.get("cfd_alt") is not None and row["cfd_alt"] is not None and abs(row["cfd_alt"] - site["cfd_alt"]) > 1e-9:
                    n_bad += 1
        title = "Cross-check agrees" if report["coordinates_agree"] and n_bad == 0 else "Cross-check disagrees"
        resolution = _resolution(action, observations)
        action["guide"].resolution = resolution
        return {
            "title": title,
            "detail": (
                f"coordinate agreement={report['coordinates_agree']}, CFD disagreements={n_bad}. "
                f"Resolver: {resolution['action']}. {resolution['reason']}"
            ),
            "payload": {"n_bad": n_bad, "coordinates_agree": report["coordinates_agree"], "resolution": resolution},
        }
    if action["tool"] == "evo2_score":
        opinion = observations["site"]["evo2"]
        if opinion.get("direction"):
            return {"title": f"Evo 2 and CFD: {opinion['direction']}", "detail": opinion["reading"], "payload": opinion}
        reason = opinion["ref"].get("reason") or opinion["alt"].get("reason")
        return {"title": "Evo 2 unavailable", "detail": reason or "No score returned.", "payload": opinion}
    return {"title": "Checked", "detail": "", "payload": {}}


def _revise(action, observations, verification, genome: Genome) -> dict | None:
    if action["tool"] == "fetch_af" and observations["site"]["conflicts"]:
        return {"title": "Dropped alleles that contradict the reference", "detail": "A gnomAD ref base that is not in the FASTA is an assembly conflict, not a haplotype.", "payload": {"dropped": len(observations["site"]["conflicts"])}}
    if action["tool"] != "crispr_me_check" or action.get("corrected"):
        if action["tool"] == "crispr_me_check" and action.get("corrected"):
            resolution = verification["payload"].get("resolution") or {}
            if verification["payload"]["n_bad"] == 0 and verification["payload"]["coordinates_agree"]:
                return {
                    "title": "Re-run agrees with the engine",
                    "detail": "The corrected window builder now matches the engine CFD on this wave.",
                    "payload": {"agreed": True, "resolution": resolution},
                }
            if resolution.get("action") == "ESCALATE_GPU_ENUMERATION":
                action["guide"].force_expand = True
            return {
                "title": "Re-run still disagrees",
                "detail": "Engine scores are kept. "
                + ("NAG and NGA will be searched for this guide." if action["guide"].force_expand else "No PAM expansion from this disagreement."),
                "payload": {"agreed": False, "resolution": resolution},
            }
        return None
    if verification["payload"]["n_bad"] == 0 and verification["payload"]["coordinates_agree"]:
        action["guide"].phase = "evo"
        return None
    report = observations["report"]
    notes = []
    by_key = {(s["chrom"], s["pos"], s["strand"], s["pam_class"]): s for s in action["guide"].wave_sites()}
    strand_n = 0
    phase_n = 0
    window_n = 0
    for row in report["rows"]:
        site = by_key.get((row["chrom"], row["pos"], row["strand"], row["pam_class"]))
        if site is None or row["ref_window"] == site.get("ref_window"):
            if site and site.get("alt_window") and row["alt_window"] != site["alt_window"]:
                phase_n += _phasing_note(site, row, False, notes)
            continue
        if revcomp(row["ref_window"]) == site.get("ref_window"):
            strand_n += 1
            _phasing_note(site, row, True, notes)
            if site.get("alt_window") and revcomp(row["alt_window"]) != site["alt_window"]:
                phase_n += 1
            continue
        if _window_shift(genome, site):
            window_n += 1
    action["guide"].phase = "recheck"
    n_minus = sum(1 for site in action["guide"].wave_sites() if site["strand"] == "-")
    hypothesis = choose(
        legal=["strand", "window", "phasing"],
        policy="strand" if n_minus else "window",
        observation=(
            f"coordinate agreement={report['coordinates_agree']}. "
            f"minus-strand windows={n_minus}. "
            f"Resolver action: {action['guide'].resolution['action'] if action['guide'].resolution else 'unset'}."
        ),
    )
    explained = "strand" if strand_n else "window" if window_n else "phasing" if phase_n else "none"
    detail = (
        f"Hypothesis: {hypothesis['tool']} ({hypothesis['source']}). {hypothesis['thought']} "
        f"Tested window shift ({window_n} explained), strand ({strand_n} minus-strand windows "
        f"were not reverse-complemented), and phasing ({phase_n} compound haplotypes). "
        f"The tests point at {explained}. "
        "The checker was wrong on strand wherever a reverse-complement restores the engine window. "
        "Where it stacked every SNV, phase is not known; the engine keeps the single highest-AF allele."
    )
    return {
        "title": "Disagreement investigated",
        "detail": detail,
        "payload": {
            "strand": strand_n,
            "phasing": phase_n,
            "window": window_n,
            "notes": notes[:8],
            "hypothesis": hypothesis,
            "explained": explained,
        },
    }


def _resolution(action, observations) -> dict:
    report = observations["report"]
    by_key = {(s["chrom"], s["pos"], s["strand"], s["pam_class"]): s for s in action["guide"].wave_sites()}
    results = []
    for row in report["rows"]:
        site = by_key.get((row["chrom"], row["pos"], row["strand"], row["pam_class"]))
        if site is None:
            continue
        windows_match = row.get("ref_window") == site.get("ref_window")
        alt = site.get("alt_window")
        if alt and row.get("alt_window") and row["alt_window"] != alt and revcomp(row["alt_window"]) != alt:
            windows_match = False
        variant = site.get("applied_variant") or {}
        results.append(analyze(
            cfd_ref=site.get("cfd_ref"),
            cfd_alt=site.get("cfd_alt"),
            checker_ref=row.get("cfd_ref"),
            checker_alt=row.get("cfd_alt"),
            af=variant.get("af"),
            coordinates_agree=bool(report["coordinates_agree"]),
            windows_match=windows_match,
        ))
    return resolve_sites(results, coordinates_agree=bool(report["coordinates_agree"]))


def _phasing_note(site, row, strand_fixed: bool, notes: list) -> int:
    if not site.get("common_variants"):
        return 0
    aligned = revcomp(row["alt_window"]) if strand_fixed else row["alt_window"]
    if aligned == site.get("alt_window"):
        return 0
    plus = None
    # The compound haplotype is what the legacy checker builds. Record its CFD
    # only as a sensitivity, after orienting it the way the engine does.
    notes.append({"pos": site["pos"], "checker_phasing": row["phasing"], "engine_phasing": "highest-AF"})
    return 1 if plus is None else 1


def _window_shift(genome: Genome, site: dict) -> bool:
    for delta in range(-6, 7):
        if delta == 0:
            continue
        try:
            window = genome.slice(site["chrom"], site["pos"] + delta, 23)
        except ValueError:
            continue
        if window == site["ref_window"] or revcomp(window) == site["ref_window"]:
            return True
    return False


def _population_cfd(site: dict) -> float:
    if site.get("cfd_alt") is None:
        return site["cfd_ref"]
    return site["cfd_alt"]


def _verdict(sites: list[dict], key: str) -> str:
    flagged = False
    for site in sites:
        if site.get("on_target"):
            continue
        score = site["cfd_ref"] if key == "ref" else _population_cfd(site)
        if score >= T.ACTIONABLE_CFD:
            flagged = True
    return "FLAGGED" if flagged else "CLEAR"


def _flipped(guide: _Guide) -> bool:
    ngg = [s for s in guide.sites if s["wave"] == "NGG"]
    if not ngg:
        return False
    return _verdict(ngg, "ref") != _verdict(ngg, "pop")


def _assemble(request, population, guides: list[_Guide], genome: Genome, env: dict, claims: list[str]) -> dict:
    rows = []
    off = []
    for guide in guides:
        for site in guide.sites:
            if site.get("on_target"):
                continue
            off.append(site)
            if site.get("cfd_alt") is None:
                continue
            same = site["mm_ref"] == site["mm_alt"]
            changed = abs(site["cfd_ref"] - site["cfd_alt"]) > 1e-12
            if not changed:
                changed_at = [i for i in range(23) if site["ref_window"][i] != site["alt_window"][i]]
                if changed_at == [20]:
                    reading = (
                        "The substitution is the N of the PAM. CFD multiplies mismatch weights "
                        "and the last two PAM bases, so this base does not move the score."
                    )
                else:
                    continue
                above_ref = site["cfd_ref"] >= T.ACTIONABLE_CFD
                above_alt = site["cfd_alt"] >= T.ACTIONABLE_CFD
                variant = site.get("applied_variant") or {}
                rows.append({
                    "guide_id": site["guide_id"],
                    "chrom": site["chrom"],
                    "pos": site["pos"],
                    "strand": site["strand"],
                    "pam_class": site["pam_class"],
                    "mm_ref": site["mm_ref"],
                    "mm_alt": site["mm_alt"],
                    "cfd_ref": site["cfd_ref"],
                    "cfd_alt": site["cfd_alt"],
                    "delta": 0.0,
                    "af": variant.get("af"),
                    "variant": f"{variant.get('ref', '')}>{variant.get('alt', '')}" if variant else "",
                    "variant_pos": variant.get("pos"),
                    "same_mismatch_count": True,
                    "above_ref": above_ref,
                    "above_alt": above_alt,
                    "reading": reading,
                    "ref_window": site["ref_window"],
                    "alt_window": site.get("alt_window"),
                    "pairs": site.get("pairs") or [],
                    "evo2": site.get("evo2"),
                })
                continue
            above_ref = site["cfd_ref"] >= T.ACTIONABLE_CFD
            above_alt = site["cfd_alt"] >= T.ACTIONABLE_CFD
            if same and not above_ref and not above_alt:
                reading = "Method behavior: the mismatch count is unchanged and both scores stay below 0.2."
            elif above_ref != above_alt:
                reading = "The haplotype crosses the registered 0.2 threshold."
            else:
                reading = "CFD moved and the mismatch count changed. Both sides are on the same side of 0.2."
            variant = site.get("applied_variant") or {}
            rows.append({
                "guide_id": site["guide_id"],
                "chrom": site["chrom"],
                "pos": site["pos"],
                "strand": site["strand"],
                "pam_class": site["pam_class"],
                "mm_ref": site["mm_ref"],
                "mm_alt": site["mm_alt"],
                "cfd_ref": site["cfd_ref"],
                "cfd_alt": site["cfd_alt"],
                "delta": site["cfd_alt"] - site["cfd_ref"],
                "af": variant.get("af"),
                "variant": f"{variant.get('ref', '')}>{variant.get('alt', '')}" if variant else "",
                "variant_pos": variant.get("pos"),
                "same_mismatch_count": same,
                "above_ref": above_ref,
                "above_alt": above_alt,
                "reading": reading,
                "ref_window": site["ref_window"],
                "alt_window": site.get("alt_window"),
                "pairs": site.get("pairs") or [],
                "evo2": site.get("evo2"),
            })
    same_rows = [r for r in rows if r["same_mismatch_count"] and abs(r["delta"]) > 1e-12]
    above = [r for r in same_rows if r["above_ref"] or r["above_alt"]]
    movers = [r for r in rows if abs(r["delta"]) > 1e-12]
    minority = [r for r in movers if (r.get("af") or 0) > 0.5]
    if minority:
        top = max(minority, key=lambda r: r["af"])
        lead = (
            f"The reference base is the minority allele at {top['chrom']}:{top['variant_pos']} "
            f"in {POPULATIONS[population]} (AF {top['af']:.2f})."
        )
    elif movers:
        top = max(movers, key=lambda r: r["af"] or 0)
        lead = (
            f"The highest overlapping {POPULATIONS[population]} frequency in this run is AF {top['af']:.3f} "
            f"at {top['chrom']}:{top['variant_pos']}. It is not a majority allele."
        )
    elif rows:
        top = max(rows, key=lambda r: r["af"] or 0)
        lead = (
            f"A {POPULATIONS[population]} variant at AF {top['af']:.2f} overlaps "
            f"{top['chrom']}:{top['variant_pos']}, and the CFD score does not move. "
            f"{top['reading']}"
        )
    else:
        lead = (
            f"No common {POPULATIONS[population]} SNV (AF ≥ {T.COMMON_AF}) overlaps an off-target window "
            f"in the loaded reference for this panel."
        )
    n_off = len(off)
    n_same = len(same_rows)
    n_above = len(above)
    n_fetched_common = sum(1 for s in off if s.get("common_variants"))
    guide_rows = []
    for guide in guides:
        ngg = [s for s in guide.sites if s["wave"] == "NGG" and not s.get("on_target")]
        guide_rows.append({
            "id": guide.spec["id"],
            "spacer": guide.spec["spacer"],
            "note": guide.spec["note"],
            "reference_verdict": _verdict(ngg, "ref"),
            "population_verdict": _verdict(ngg, "pop"),
            "flipped": _verdict(ngg, "ref") != _verdict(ngg, "pop"),
            "escalated": guide.escalated,
            "n_sites": len([s for s in guide.sites if not s.get("on_target")]),
        })
    record = {
        "integers": [n_off, n_same, n_above, n_fetched_common, len(guides), T.MAX_MISMATCHES],
        "floats": [T.ACTIONABLE_CFD, T.COMMON_AF] + [round(r["cfd_ref"], 4) for r in rows] + [round(r["cfd_alt"], 4) for r in rows] + [round(r["af"], 4) for r in rows if r.get("af") is not None],
    }
    # Integers that appear in claims must be in the record. Counts are enough;
    # individual CFD values are floats.
    claims = [
        f"{n_same} off-target sites keep the same mismatch count and change CFD. {n_above} of them reach 0.2 on either haplotype.",
        f"{n_off} off-target sites were enumerated. {n_fetched_common} "
        f"{'carries' if n_fetched_common == 1 else 'carry'} a common variant at allele frequency 0.05 or higher. "
        f"Frequencies are gnomAD v4 {POPULATIONS[population]} calls.",
        "CFD is the Doench 2016 product of mismatch weights and a PAM weight. Position 1 is PAM-distal.",
        "Evo 2, when a score is returned, is a sequence prior. A disagreement with the CFD delta is not a cutting probability.",
    ]
    limitations = [
        f"The loaded reference is {genome.path}. Sites outside that FASTA are not in this run.",
        "Substitutions only. Bulges and indels are out of scope.",
        "The primary haplotype is the single highest-AF SNV. gnomAD does not phase sites inside the window.",
        "The cross-check is an in-repo second enumerator unless CRISPRME_BIN points at an external binary.",
        "NAG and NGA are searched only after a guide's NGG verdict flips.",
        "τ = 0.2 is a registered method cutoff, not a measured cutting rate.",
        "The orchestrator model runs only when NEBIUS_API_KEY+NEBIUS_BASE_URL or NVIDIA_API_KEY is set. Without an eligible key the model is not called.",
        "No GUIDE-seq correlation is reported unless a real data/guideseq.tsv is loaded.",
        "CRISPRme+ is not called unless CRISPRME_BIN is set. The second checker is in-repo.",
    ]
    return {
        "refused": False,
        "lead": lead,
        "gene": request["gene"],
        "population": population,
        "population_name": POPULATIONS[population],
        "threshold": T.ACTIONABLE_CFD,
        "common_af": T.COMMON_AF,
        "guides": guide_rows,
        "table": rows,
        "same_mismatch": same_rows,
        "n_off": n_off,
        "n_same": n_same,
        "n_above": n_above,
        "environment": env,
        "limitations": limitations,
        "qualifier": QUALIFIER,
        "record": record,
        "sites": [_public_site(s) for s in off],
        "_claims": claims,
    }


def _public_site(site: dict) -> dict:
    keep = ("guide_id", "chrom", "pos", "strand", "mismatches", "pam_class", "site23", "on_target",
            "cfd_ref", "cfd_alt", "mm_ref", "mm_alt", "ref_window", "alt_window", "envelope", "wave",
            "evo2", "pairs", "applied_variant", "common_variants")
    return {key: site.get(key) for key in keep}


def _default_fetch(site: dict, population: str) -> list[dict]:
    found = fetch_af(site["chrom"], site["pos"], site["pos"] + 22, population)
    return [v.as_dict() for v in found]


def _environment(genome: Genome) -> dict:
    from popoff.token_factory import eligible
    env = eligible()
    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "threads": os.cpu_count(),
        "genome": genome.path,
        "genome_sha256": genome.sha256,
        "genome_sha256_12": genome.sha256[:12],
        "network_gnomad": network_ok(),
        "nvidia_api_key": bool(os.environ.get("NVIDIA_API_KEY")),
        "nebius_api_key": env.get("has_nebius_key"),
        "eligible_provider": env.get("provider"),
        "crisprme_bin": os.environ.get("CRISPRME_BIN"),
    }


def file_sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()
