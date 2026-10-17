"""NeMo Agent Toolkit registration.

`nat run --config_file workflows/crispr_ops.yml --input "HBB sas"` is the
ReAct workflow: Nemotron chooses among the registered tools. It needs
NVIDIA_API_KEY. `workflows/popoff_audit.yml` runs the same audit with the
registered policy choosing every branch, which is what the demo uses when
that key is unset.
"""

from __future__ import annotations

import json
from pathlib import Path

from nat.builder.function_info import FunctionInfo
from nat.cli.register_workflow import register_function
from nat.data_models.function import FunctionBaseConfig
from pydantic import Field

from popoff.agent import audit_events

# Registers the ReAct tools as an import side effect.
import popoff.nat_tools  # noqa: F401, E402

ROOT = Path(__file__).resolve().parents[2]


class PopoffAuditConfig(FunctionBaseConfig, name="popoff_audit"):
    genome_fasta: str = Field(description="Reference FASTA to search")
    population: str = Field(default="sas", description="gnomAD population label")
    gene: str = Field(default="HBB")


@register_function(config_type=PopoffAuditConfig)
async def popoff_audit(config: PopoffAuditConfig, builder):
    del builder

    async def _run(query: str) -> str:
        gene = config.gene
        population = config.population
        note = query or ""
        parts = query.split()
        if parts:
            gene = parts[0]
        if len(parts) > 1:
            population = parts[1]
        events = list(audit_events({
            "gene": gene,
            "population": population,
            "note": note,
            "fasta": str(ROOT / config.genome_fasta) if not Path(config.genome_fasta).is_absolute() else config.genome_fasta,
        }))
        payload = events[-1].get("payload") or {}
        payload["trace_stages"] = [event["stage"] for event in events]
        return json.dumps(payload, default=str)

    yield FunctionInfo.from_fn(
        _run,
        description="Population-aware CRISPR off-target audit. Input is 'GENE POPULATION'.",
    )
