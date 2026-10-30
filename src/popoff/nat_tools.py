"""Tools the NeMo Agent Toolkit ReAct workflow can call.

The model chooses the tool. Each function is deterministic and returns an
error when the call is not legal yet, so a bad choice comes back as an
observation instead of a silent skip. One process runs one audit.
"""

from __future__ import annotations

import json
from pathlib import Path

try:
    from nat.builder.function_info import FunctionInfo
    from nat.cli.register_workflow import register_function
    from nat.data_models.function import FunctionBaseConfig
except ModuleNotFoundError as exc:  # pragma: no cover - optional dependency in test environments
    if exc.name != "nat":
        raise

    class FunctionInfo:
        @staticmethod
        def from_fn(fn, description=None, **kwargs):
            return {"fn": fn, "description": description, **kwargs}

    def register_function(config_type=None):
        def decorator(fn):
            return fn
        return decorator

    class FunctionBaseConfig:
        def __init_subclass__(cls, **kwargs):
            super().__init_subclass__()
            for key, value in kwargs.items():
                setattr(cls, key, value)

        def __init__(self, **kwargs):
            for key, value in kwargs.items():
                setattr(self, key, value)

    class _NoopField:
        def __init__(self, default=None, **kwargs):
            self.default = default

    def Field(default=None, **kwargs):
        return _NoopField(default=default, **kwargs)
else:
    from pydantic import Field

from popoff import thresholds as T
from popoff.agent import (
    QUALIFIER,
    _Guide,
    _act as run_tool,
    _assemble,
    _default_fetch,
    _environment,
    _flipped,
    _revise,
    _triage,
    _verify,
)
from popoff.evo2 import evo2_score
from popoff.citations import verify_citation
from popoff.genome import GENOMES
from popoff.gnomad import POPULATIONS
from popoff.guides import panel_for
from popoff.policy import review

ROOT = Path(__file__).resolve().parents[2]
SESSION: "AuditSession | None" = None


class _FastaConfig(FunctionBaseConfig):
    genome_fasta: str = Field(default="data/genome/chr11.fa")


class BeginAuditConfig(_FastaConfig, name="popoff_begin_audit"):
    pass


class EnumerateConfig(FunctionBaseConfig, name="popoff_enumerate_sites"):
    pass


class FetchScoreConfig(FunctionBaseConfig, name="popoff_fetch_and_score"):
    pass


class CheckConfig(FunctionBaseConfig, name="popoff_crispr_me_check"):
    pass


class ResolveConfig(FunctionBaseConfig, name="popoff_resolve_disagreement"):
    pass


class InvestigateConfig(FunctionBaseConfig, name="popoff_investigate"):
    pass


class EvoConfig(FunctionBaseConfig, name="popoff_evo2_score"):
    pass


class GpuConfig(FunctionBaseConfig, name="popoff_gpu_enumerate"):
    pass


class CiteConfig(FunctionBaseConfig, name="popoff_verify_citation"):
    pass


class ReportConfig(FunctionBaseConfig, name="popoff_audit_report"):
    pass


class AuditSession:
    def __init__(self, fasta: str, gene: str, population: str, note: str = ""):
        self.fasta = fasta
        self.gene = gene
        self.population = population
        self.note = note
        text = f"{gene} {population} {note}"
        self.refusal = review(text)
        self.guides: list[_Guide] = []
        self.genome = None
        self.locus = None
        self.request: dict = {}
        self.report_body: dict | None = None
        if self.refusal or population not in POPULATIONS:
            return
        self.guides = [_Guide(spec) for spec in panel_for(gene, None)]
        self.genome = GENOMES.get(fasta)
        self.locus = T.GENE_LOCI.get(gene.upper())
        self.request = {"gene": gene, "population": population, "note": note, "fasta": fasta}

    @property
    def refused(self) -> bool:
        return self.refusal is not None or self.population not in POPULATIONS

    def guide(self, guide_id: str) -> _Guide:
        for guide in self.guides:
            if guide.spec["id"] == guide_id:
                return guide
        raise KeyError(guide_id)

    def summary(self) -> dict:
        if self.refusal:
            return {"refused": True, "reason": self.refusal["reason"], "qualifier": QUALIFIER}
        if self.population not in POPULATIONS:
            return {"refused": True, "reason": f"{self.population} is not a gnomAD population label."}
        return {
            "refused": False,
            "gene": self.gene,
            "population": self.population,
            "guides": [guide.spec["id"] for guide in self.guides],
            "qualifier": QUALIFIER,
        }

    def enumerate(self, guide_id: str) -> dict:
        guide = self._ready(guide_id)
        if isinstance(guide, dict):
            return guide
        if guide.phase != "enumerate":
            return {"error": f"{guide_id} is in phase {guide.phase}, not enumerate."}
        obs = self._act({"tool": "enumerate_sites", "guide": guide, "title": "enumerate_sites"})
        return obs["public"]

    def fetch_and_score(self, guide_id: str) -> dict:
        guide = self._ready(guide_id)
        if isinstance(guide, dict):
            return guide
        if not guide.wave_sites():
            return {"error": "Enumerate this guide before fetching allele frequencies."}
        steps = 0
        while steps < 500:
            steps += 1
            pending = [site for site in guide.wave_sites() if not site.get("fetched")]
            if pending:
                self._act({"tool": "fetch_af", "guide": guide, "site": pending[0]})
                continue
            scoring = [site for site in guide.wave_sites() if site.get("needs_score") and not site.get("scored")]
            if scoring:
                self._act({"tool": "score_cfd", "guide": guide, "site": scoring[0]})
                continue
            break
        common = sum(1 for site in guide.wave_sites() if site.get("common_variants"))
        return {"guide_id": guide_id, "sites": len(guide.wave_sites()), "common_variant_sites": common}

    def cross_check(self, guide_id: str) -> dict:
        guide = self._ready(guide_id)
        if isinstance(guide, dict):
            return guide
        if any(not site.get("fetched") for site in guide.wave_sites()):
            return {"error": "Fetch allele frequencies before the cross-check."}
        obs = self._act({"tool": "crispr_me_check", "guide": guide, "corrected": False})
        verification = _verify({"tool": "crispr_me_check", "guide": guide, "corrected": False}, obs, self.genome)
        revision = _revise({"tool": "crispr_me_check", "guide": guide, "corrected": False}, obs, verification, self.genome)
        return {
            "coordinates_agree": obs["public"]["coordinates_agree"],
            "implementation": obs["public"]["implementation"],
            "resolution": verification["payload"].get("resolution"),
            "revision": None if revision is None else revision["title"],
        }

    def resolve(self, guide_id: str) -> dict:
        guide = self._ready(guide_id)
        if isinstance(guide, dict):
            return guide
        if not guide.resolution:
            return {"error": "Cross-check this guide before resolving a disagreement."}
        return guide.resolution

    def investigate(self, guide_id: str) -> dict:
        guide = self._ready(guide_id)
        if isinstance(guide, dict):
            return guide
        if guide.phase != "recheck":
            return {"error": "No open disagreement. The cross-check has to name one before investigate."}
        obs = self._act({"tool": "crispr_me_check", "guide": guide, "corrected": True})
        verification = _verify({"tool": "crispr_me_check", "guide": guide, "corrected": True}, obs, self.genome)
        revision = _revise({"tool": "crispr_me_check", "guide": guide, "corrected": True}, obs, verification, self.genome)
        return {"agreed": revision["payload"]["agreed"], "detail": revision["detail"]}

    def evo(self, guide_id: str) -> dict:
        guide = self._ready(guide_id)
        if isinstance(guide, dict):
            return guide
        if guide.phase not in ("evo", "expand", "check"):
            return {"error": f"Evo 2 is not legal in phase {guide.phase}."}
        if guide.wave_name not in guide.triaged_waves:
            _triage(guide)
            guide.phase = "evo"
        scored = 0
        while guide.evo_queue and scored < T.EVO2_BUDGET:
            self._act({"tool": "evo2_score", "guide": guide, "site": guide.evo_queue[0]})
            scored += 1
        guide.phase = "expand"
        return {"scored": scored, "note": guide.triage_note}

    def gpu_enumerate(self, guide_id: str) -> dict:
        guide = self._ready(guide_id)
        if isinstance(guide, dict):
            return guide
        if guide.escalated:
            return {"error": "This guide was already expanded to NAG and NGA."}
        if not (_flipped(guide) or guide.force_expand):
            return {
                "error": (
                    "PAM expansion is not legal. The NGG verdict did not flip, and the "
                    "checker disagreement was not left unresolved."
                )
            }
        guide.escalated = True
        guide.waves.append(T.EXPANDED_PAM)
        guide.wave_i += 1
        guide.phase = "enumerate"
        obs = self._act({"tool": "enumerate_sites", "guide": guide})
        return obs["public"]

    def cite(self, claim: str) -> dict:
        if self.report_body is None:
            return {"error": "Call audit_report before verify_citation."}
        return verify_citation(claim, self.report_body["record"])

    def report(self) -> dict:
        blocked = self._blocked()
        if blocked:
            return blocked
        body = _assemble(self.request, self.population, self.guides, self.genome, _environment(self.genome), [])
        body.pop("_claims", None)
        self.report_body = body
        return {
            "lead": body["lead"],
            "qualifier": body["qualifier"],
            "guides": body["guides"],
            "limitations": body["limitations"],
        }

    def _ready(self, guide_id: str):
        blocked = self._blocked()
        if blocked:
            return blocked
        try:
            return self.guide(guide_id)
        except KeyError:
            return {"error": f"Unknown guide {guide_id}.", "guides": [g.spec["id"] for g in self.guides]}

    def _blocked(self) -> dict | None:
        if self.refusal:
            return {"refused": True, "reason": self.refusal["reason"]}
        if self.population not in POPULATIONS:
            return {"refused": True, "reason": f"{self.population} is not a gnomAD population label."}
        return None

    def _act(self, action: dict) -> dict:
        return run_tool(
            action, self.request, self.population, self.genome, self.locus, _default_fetch, evo2_score
        )


def _open(config: BeginAuditConfig, query: str) -> AuditSession:
    global SESSION
    parts = query.split()
    gene = parts[0] if parts else "HBB"
    population = parts[1] if len(parts) > 1 else "sas"
    note = " ".join(parts[2:])
    fasta = config.genome_fasta
    if not Path(fasta).is_absolute():
        fasta = str(ROOT / fasta)
    SESSION = AuditSession(fasta, gene, population, note)
    return SESSION


def _dumps(payload: dict) -> str:
    return json.dumps(payload, default=str)


def _session_or_error() -> AuditSession | str:
    if SESSION is None:
        return _dumps({"error": "Call begin_audit first."})
    return SESSION


@register_function(config_type=BeginAuditConfig)
async def popoff_begin_audit(config: BeginAuditConfig, builder):
    del builder

    async def _run(query: str) -> str:
        return _dumps(_open(config, query).summary())

    yield FunctionInfo.from_fn(
        _run,
        description=(
            "Start an audit. Input is 'GENE POPULATION' plus an optional note. "
            "Refuses a personal genome or a clinical question and returns the reason. "
            "Call this before any other tool."
        ),
    )


@register_function(config_type=EnumerateConfig)
async def popoff_enumerate_sites(_config: EnumerateConfig, builder):
    del builder

    async def _run(guide_id: str) -> str:
        session = _session_or_error()
        if isinstance(session, str):
            return session
        return _dumps(session.enumerate(guide_id))

    yield FunctionInfo.from_fn(
        _run,
        description="Enumerate NGG off-target sites for one guide id at mismatch <= 4. CPU enumerator.",
    )


@register_function(config_type=FetchScoreConfig)
async def popoff_fetch_and_score(_config: FetchScoreConfig, builder):
    del builder

    async def _run(guide_id: str) -> str:
        session = _session_or_error()
        if isinstance(session, str):
            return session
        return _dumps(session.fetch_and_score(guide_id))

    yield FunctionInfo.from_fn(
        _run,
        description="Fetch gnomAD allele frequencies and score CFD for every site of one guide.",
    )


@register_function(config_type=CheckConfig)
async def popoff_crispr_me_check(_config: CheckConfig, builder):
    del builder

    async def _run(guide_id: str) -> str:
        session = _session_or_error()
        if isinstance(session, str):
            return session
        return _dumps(session.cross_check(guide_id))

    yield FunctionInfo.from_fn(
        _run,
        description=(
            "Cross-check one guide with the second enumerator. "
            "This is the in-repo checker unless CRISPRME_BIN is set. It is not CRISPRme+ by default."
        ),
    )


@register_function(config_type=ResolveConfig)
async def popoff_resolve_disagreement(_config: ResolveConfig, builder):
    del builder

    async def _run(guide_id: str) -> str:
        session = _session_or_error()
        if isinstance(session, str):
            return session
        return _dumps(session.resolve(guide_id))

    yield FunctionInfo.from_fn(
        _run,
        description=(
            "Read the disagreement resolver for one guide. "
            "INVESTIGATE means call investigate. ESCALATE_EVO2 means call evo2_score. "
            "ESCALATE_GPU_ENUMERATION means call gpu_enumerate. ACCEPT is not a safety verdict."
        ),
    )


@register_function(config_type=InvestigateConfig)
async def popoff_investigate(_config: InvestigateConfig, builder):
    del builder

    async def _run(guide_id: str) -> str:
        session = _session_or_error()
        if isinstance(session, str):
            return session
        return _dumps(session.investigate(guide_id))

    yield FunctionInfo.from_fn(
        _run,
        description="Re-run the cross-check after testing strand, window start, and variant phasing.",
    )


@register_function(config_type=EvoConfig)
async def popoff_evo2_score(_config: EvoConfig, builder):
    del builder

    async def _run(guide_id: str) -> str:
        session = _session_or_error()
        if isinstance(session, str):
            return session
        return _dumps(session.evo(guide_id))

    yield FunctionInfo.from_fn(
        _run,
        description=(
            "Score triaged windows with Evo 2. The result is a sequence-context prior, "
            "not a cutting probability. Sites far from the CFD threshold are not sent."
        ),
    )


@register_function(config_type=GpuConfig)
async def popoff_gpu_enumerate(_config: GpuConfig, builder):
    del builder

    async def _run(guide_id: str) -> str:
        session = _session_or_error()
        if isinstance(session, str):
            return session
        return _dumps(session.gpu_enumerate(guide_id))

    yield FunctionInfo.from_fn(
        _run,
        description=(
            "Expand one guide to NAG and NGA. Legal only after a verdict flip or an unresolved "
            "checker disagreement. Uses the CUDA binary when it exists, otherwise the CPU enumerator."
        ),
    )


@register_function(config_type=ReportConfig)
async def popoff_audit_report(_config: ReportConfig, builder):
    del builder

    async def _run(unused: str) -> str:
        del unused
        session = _session_or_error()
        if isinstance(session, str):
            return session
        return _dumps(session.report())

    yield FunctionInfo.from_fn(
        _run,
        description="Return the lead, qualifier, verdicts, and limitations. Do not add numbers that are not in this JSON.",
    )


@register_function(config_type=CiteConfig)
async def popoff_verify_citation(_config: CiteConfig, builder):
    del builder

    async def _run(claim: str) -> str:
        session = _session_or_error()
        if isinstance(session, str):
            return session
        return _dumps(session.cite(claim))

    yield FunctionInfo.from_fn(
        _run,
        description="Accept a claim only when its numbers are in the audit report and it is not a clinical conclusion.",
    )
