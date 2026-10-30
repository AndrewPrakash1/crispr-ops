#!/usr/bin/env bash
# Load orchestrator credentials without committing secrets.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

if [[ -f "$ROOT/.env.local" ]]; then
  set -a
  # shellcheck disable=SC1091
  source "$ROOT/.env.local"
  set +a
elif [[ -f "$HOME/.commandcode/auth.json" ]]; then
  export NVIDIA_API_KEY="$(
    "$ROOT/.venv/bin/python" -c 'import json, pathlib; print(json.loads(pathlib.Path.home().joinpath(".commandcode/auth.json").read_text())["apiKey"])'
  )"
  export ORCHESTRATOR_URL="${ORCHESTRATOR_URL:-https://api.commandcode.ai/provider}"
  export ORCHESTRATOR_MODEL="${ORCHESTRATOR_MODEL:-nvidia/nemotron-3-ultra-550b-a55b}"
fi

if [[ -z "${NVIDIA_API_KEY:-}" && -z "${NEBIUS_API_KEY:-}" ]]; then
  echo "No API key. Create $ROOT/.env.local or sign in with the Command Code CLI." >&2
  exit 1
fi

exec "$@"
