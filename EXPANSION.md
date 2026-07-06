# sympy_modal Expansion Plan

This document tracks the implementation of new features to expand the capabilities of `sympy_modal`.

## 1. Enhancing Expressiveness (Broadening the Logic)
- [x] **Multi-Agent Epistemic Logic:**
  - [x] Implement `AgentBox` (e.g., $\Box_A$) representing knowledge of specific agents.
  - [x] Implement `CommonKnowledge` operator.
  - [x] Add tests for multi-agent epistemic properties.
- [x] **Rich Temporal Logic (LTL/CTL):**
  - [x] Implement `Next` ($X$) operator.
  - [x] Implement `Until` ($U$) operator.
  - [x] Add tests for LTL/CTL operators.
- [x] **Opt-in Classical Logic:**
  - [x] Add `allow_classical=True` flag to `ProofContext`.
  - [x] Automatically inject classical axioms (e.g., Law of Excluded Middle) when the flag is true.
  - [x] Silence warnings for classical axioms when the flag is true.
  - [x] Add tests verifying classical logic proofs work when enabled.

## 2. Upgrading Proof Automation & Tactics
- [x] **Interactive Tactics Engine:**
  - [x] Implement basic scriptable tactics (e.g., `apply`, `rewrite`, `induction`).
  - [x] Add tests for the tactics engine.
- [x] **SMT Solver Integration (SymPy Native):**
  - [x] Integrate `sympy.logic.inference` for solving intermediate non-modal propositional goals.
  - [x] Add tests verifying improved proof search speed and capability.

## 3. Practical Verification Features
- [x] **Semantic Model Checking:**
  - [x] Implement a "Semantic Evaluator" taking finite Kripke models.
  - [x] Evaluate formulas against the specific models.
  - [x] Add tests for semantic model checking on defined worlds and relations.
- [x] **LLM-Assisted Formalisation Bridge:**
  - [x] Implement a utility class for calling the Gemini 2.5 Flash API via OpenRouter.
  - [x] Add prompt generation for converting natural language to `sympy_modal` syntax.
  - [x] Include comments about potential future expansion to other models.
  - [x] Add tests using mocked LLM responses.
- [x] **Export to External Proof Assistants:**
  - [x] Implement `export_lean()` robust string generator via AST on `ProofTerm`.
  - [x] Implement `export_coq()` robust string generator via AST on `ProofTerm`.
  - [x] Ensure the strings are structured to be interpreted robustly.
  - [x] Add tests for exporting proof terms to Lean and Coq formats.

## 4. Documentation
- [x] **Update `sympy_modal.md`:**
  - [x] Document all the newly added features.
  - [x] Replace silent `assert` statements in tutorial examples with meaningful `print` statements.
  - [x] Expand example explanations to include structured two-paragraph introductions and step-by-step code walkthroughs.
