# Development Workflow

- Uses `pytest` as the test runner. Confidence: 0.95
- Expects tests to cover refusal, illegal-choice rejection, fallback/cached mode, and scientific edge cases. Confidence: 0.85
- Prefers reproducible artifacts saved to the repo (JSON traces, CSV grids, SQLite run records). Confidence: 0.9
- Wants one-command re-runs for grids and eligibility proofs (e.g. `popoff panel`, `scripts/prove_eligibility.py`). Confidence: 0.85
- Comfortable with environment-variable driven configuration for secrets and endpoints. Confidence: 0.8
