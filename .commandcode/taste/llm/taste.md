# LLM Eligibility & Orchestration

- Requires eligible Nemotron calls to go through Nebius Token Factory (`NEBIUS_API_KEY` + `NEBIUS_BASE_URL` from the console), not `NVIDIA_API_KEY` / `integrate.api.nvidia.com` / OpenRouter. Confidence: 0.95
- Wants eligibility proven and stored: exact model id, request id, latency, token usage in a JSON artifact. Confidence: 0.9
- Wants a model-vs-policy A/B diff captured, comparing Evo 2 spend, hypothesis order, and site order. Confidence: 0.9
- Expects graceful fallback to policy/cached mode with a clear banner when API keys or the model are unavailable. Confidence: 0.9
- Wants illegal model choices caught, rejected, and logged. Confidence: 0.9
