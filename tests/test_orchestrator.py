import popoff.orchestrator as orchestrator
from popoff.orchestrator import choose


def test_missing_key_does_not_call_the_model(monkeypatch):
    monkeypatch.delenv("NVIDIA_API_KEY", raising=False)

    def caller(_messages):
        raise AssertionError("model was called")

    decision = choose(legal=["evo2_score", "defer_evo"], policy="evo2_score", observation="delta 0.01", caller=caller)
    assert decision["source"] == "policy"
    assert decision["tool"] == "evo2_score"
    assert decision["accepted"] is True
    assert "not called" in decision["thought"]


def test_illegal_tool_and_invented_number_are_dropped(monkeypatch):
    monkeypatch.setenv("NVIDIA_API_KEY", "test-key")

    def caller(_messages):
        return '{"tool": "gpu_enumerate", "thought": "AF 0.72 makes this dangerous for the patient"}'

    decision = choose(
        legal=["evo2_score", "defer_evo"],
        policy="evo2_score",
        observation="cfd delta 0.01",
        caller=caller,
    )
    assert decision["tool"] == "evo2_score"
    assert decision["accepted"] is False
    assert decision["source"] == "policy"
    assert "0.72" not in decision["thought"]
    assert "patient" not in decision["thought"]


def test_legal_choice_keeps_a_grounded_thought(monkeypatch):
    monkeypatch.setenv("NVIDIA_API_KEY", "test-key")

    def caller(_messages):
        return 'note {"tool": "defer_evo", "thought": "delta 0.01 is under the bar"}'

    decision = choose(
        legal=["evo2_score", "defer_evo"],
        policy="evo2_score",
        observation="delta 0.01",
        caller=caller,
    )
    assert decision["tool"] == "defer_evo"
    assert decision["source"] == "nim"
    assert decision["accepted"] is True
    assert decision["model"] == orchestrator.MODEL
