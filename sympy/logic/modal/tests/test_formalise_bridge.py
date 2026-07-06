from __future__ import annotations
from sympy.logic.modal.formalise import LLMPromptBuilder

def test_llm_prompt_builder_stub():
    builder = LLMPromptBuilder(api_key="")
    prompt = builder.formalise_prompt("Necessarily p implies q")
    assert "Necessarily p implies q" in prompt
    assert "sympy.logic.modal" in prompt

    res = builder.call_gemini(prompt)
    assert res == "# LLM API not configured or requests not installed."
