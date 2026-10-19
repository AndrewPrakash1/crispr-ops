"""Command line entry points."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from popoff.agent import audit_events
from popoff.bench import run_bench
from popoff.panel import export_grid, grid_rows, load_panel, run_panel, site_rows, write_csv

ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    parser = argparse.ArgumentParser(prog="popoff")
    sub = parser.add_subparsers(dest="cmd", required=True)

    audit = sub.add_parser("audit")
    audit.add_argument("--gene", default="HBB")
    audit.add_argument("--population", default="sas")
    audit.add_argument("--fasta", default=str(ROOT / "data" / "genome" / "chr11.fa"))
    audit.add_argument("--note", default="")
    audit.add_argument("--guides", default="", help="Comma-separated 20 nt spacers. Default is the HBB panel.")

    bench = sub.add_parser("bench")
    bench.add_argument("--fasta", default=str(ROOT / "data" / "genome" / "chr11.fa"))
    bench.add_argument("--guide", default="GTAACGGCAGACTTCTCCTC")
    bench.add_argument("--threads", type=int, default=4)

    serve = sub.add_parser("serve")
    serve.add_argument("--host", default="127.0.0.1")
    serve.add_argument("--port", type=int, default=8000)

    panel = sub.add_parser("panel")
    panel.add_argument("--file", default=str(ROOT / "data" / "panel.json"))
    panel.add_argument("--fasta", default=str(ROOT / "data" / "genome" / "chr1_chr3_chr11.fa"))
    panel.add_argument("--db", default=str(ROOT / "data" / "popoff.sqlite"))
    panel.add_argument("--grid", default=str(ROOT / "data" / "grid.csv"))
    panel.add_argument("--sites", default=str(ROOT / "data" / "grid_sites.csv"))
    panel.add_argument("--json", dest="json_out", default="")

    export = sub.add_parser("export")
    export.add_argument("--run-id", action="append", default=[])
    export.add_argument("--grid", required=True)
    export.add_argument("--sites", default="")

    args = parser.parse_args()
    if args.cmd == "audit":
        guides = None
        if args.guides:
            guides = [{"id": f"g{i+1}", "spacer": s.strip()} for i, s in enumerate(args.guides.split(",")) if s.strip()]
        events = list(audit_events({
            "gene": args.gene,
            "population": args.population,
            "note": args.note,
            "fasta": args.fasta,
            "guides": guides,
        }))
        print(json.dumps(events[-1]["payload"], indent=2, default=str))
    elif args.cmd == "bench":
        print(json.dumps(run_bench(args.fasta, args.guide, threads=args.threads), indent=2))
    elif args.cmd == "serve":
        import uvicorn
        uvicorn.run("popoff.web:app", host=args.host, port=args.port, reload=False)
    elif args.cmd == "panel":
        panel_def = load_panel(args.file)
        results = run_panel(panel_def, fasta=args.fasta, db_path=args.db)
        export_grid(results, args.grid, args.sites)
        if args.json_out:
            with open(args.json_out, "w", encoding="utf-8") as handle:
                json.dump(results, handle, indent=2, default=str)
        print(json.dumps({
            "n_runs": len(results),
            "grid": args.grid,
            "sites": args.sites,
            "json": args.json_out,
        }, indent=2))
    elif args.cmd == "export":
        from popoff.store import get_store

        store = get_store()
        records = []
        for run_id in args.run_id:
            record = store.get_record(run_id)
            if record:
                records.append({"run_id": run_id, "record": record})
        if not records:
            raise SystemExit("No records found for the supplied run ids.")
        write_csv(grid_rows(records), args.grid)
        if args.sites:
            write_csv(site_rows(records), args.sites)
        print(json.dumps({"grid": args.grid, "sites": args.sites, "n_runs": len(records)}, indent=2))


if __name__ == "__main__":
    main()
