# Coding Agent Prompt: Second-Order Modal Extension for SymPy

## Overview

You are tasked with designing and implementing `sympy-modal`: a second-order modal predicate calculus extension to SymPy. This is a research-grade system combining a computer algebra interface with a proof-theoretic kernel. Read this entire prompt before writing any code. The architecture is precise and the order of implementation matters.

---

## Background and Motivation

SymPy's existing logic module (`sympy.logic`) can represent propositional and first-order formulae syntactically but provides no proof calculus, no modal operators, no Kripke semantics, and no type-theoretic foundation. It is a computer algebra system that manipulates logical syntax; it is not a proof assistant.

This extension adds the components necessary to make SymPy capable of:

1. Expressing and verifying second-order modal propositions
2. Maintaining Kripke frame semantics at runtime
3. Producing proof terms whose well-typedness constitutes a validity certificate
4. Supporting provability logic (GL) sufficient to express and verify Löb's theorem

The theoretical foundation is the Curry–Howard–Lambek correspondence: propositions are types, proofs are terms, proof checking is type checking. The system must enforce this correspondence through a trusted kernel, not merely represent it syntactically.

---

## Architecture

The system has five layers. Implement them in order; each layer depends on the one below it.

### Layer 1: Modal Type System (`sympy_modal/types.py`)

Replace SymPy's untyped predicate logic with a stratified type universe supporting second-order quantification.

**Requirements:**

- Implement a `TypeUniverse` class with cumulative levels (Type₀, Type₁, ...) sufficient to prevent Russell-style paradoxes while permitting controlled self-reference
- Implement `ModalPredicate(name, domain, arity, modality)` — a typed predicate whose validity is relative to a Kripke frame
- Implement `PredicateVariable(name, type)` — a second-order variable ranging over predicates of a given type, for use in universal and existential second-order quantification
- Implement `GuardedFixedPoint(operator, name)` — a controlled fixed-point constructor permitting the diagonal constructions required by Löb's theorem, guarded to ensure well-foundedness
- The type universe must be predicative at base level; impredicative quantification is permitted only at explicitly elevated universe levels

**Test:** The following must be expressible as a well-typed term:

```python
P = PredicateVariable('P', type=FunctionType(Universe(0), BoolType()))
lob_type = ForAllPredicates(P,
    Box(Box(P(x)) >> P(x)) >> Box(P(x))
)
assert lob_type.is_well_typed()
```

---

### Layer 2: Kripke Frame Semantics (`sympy_modal/frames.py`)

Implement Kripke frames as runtime objects that govern inference, not merely data structures.

**Requirements:**

- Implement `KripkeFrame(worlds, accessibility, axioms)` where:
  - `worlds` is a set of world identifiers
  - `accessibility` is a dict mapping each world to the set of worlds it can access
  - `axioms` is a list drawn from `{Axiom.K, Axiom.T, Axiom.B, Axiom.Four, Axiom.Five, Axiom.Lob}`
- Implement the following named constructors as class methods:
  - `KripkeFrame.K()` — minimal normal modal logic
  - `KripkeFrame.T()` — reflexive (alethic)
  - `KripkeFrame.S4()` — reflexive and transitive
  - `KripkeFrame.S5()` — reflexive, transitive, symmetric (equivalence relation)
  - `KripkeFrame.GL()` — Gödel-Löb provability logic; transitive and converse well-founded; required for Löb's theorem
  - `KripkeFrame.D()` — deontic; serial accessibility
  - `KripkeFrame.K45()` — epistemic; transitive and euclidean
- Implement `frame.validates(formula)` — checks whether a formula is valid in all worlds of this frame under its accessibility relation
- Implement `frame.is_valid_inference(premises, conclusion)` — checks whether the inference is sound in this frame
- The frame must be a live runtime object: inferences submitted to the proof context must be checked against the active frame, and unsound inferences must raise `FrameViolationError` with a precise statement of which frame condition was violated

**Test:**

```python
s4  = KripkeFrame.S4()
gl  = KripkeFrame.GL()

# □P → P is valid in S4 (reflexivity) but not in GL
assert s4.validates(Box(P) >> P)
assert not gl.validates(Box(P) >> P)

# Löb's axiom □(□P → P) → □P is valid in GL but not S4
assert gl.validates(Box(Box(P) >> P) >> Box(P))
assert not s4.validates(Box(Box(P) >> P) >> Box(P))
```

---

### Layer 3: Trusted Proof Kernel (`sympy_modal/kernel.py`)

This is the most critical component. It must be small, auditable, and complete. Every inference in the system passes through it. Nothing bypasses it.

**Requirements:**

- The kernel implements natural deduction rules for second-order modal logic:
  - Standard propositional rules: →I, →E (modus ponens), ∧I, ∧E, ∨I, ∨E, ⊥E
  - First-order rules: ∀I, ∀E, ∃I, ∃E with eigenvariable conditions enforced
  - Second-order rules: ∀²I, ∀²E ranging over predicate variables; ∃²I, ∃²E
  - Modal rules: Necessitation (if ⊢ P then ⊢ □P — only for theorems, never for hypotheses; kernel enforces this distinction strictly), K axiom distribution, and frame-specific axioms loaded from the active `KripkeFrame`
- Implement `ProofTerm(formula, derivation, source)` — a term whose well-typedness is the certificate of validity
- Implement `TrustedKernel(frame)` with the following methods:
  - `kernel.check_axiom(formula)` → `ProofTerm` or raises `NotAnAxiomError`
  - `kernel.verify_rule(rule, premises)` → `ProofTerm` or raises `InvalidInferenceError`
  - `kernel.necessitate(proof_term)` → `ProofTerm` for □P given proof of P; raises `NecessitationError` if proof_term depends on undischarged hypotheses
  - `kernel.check_term(term, type)` → `bool`; the type checker; this is the trusted base
- The kernel source must not exceed 600 lines. If it grows beyond this, the design is wrong. Complexity belongs in the proof search layer above, not in the kernel.
- The kernel must have 100% test coverage. It is the trusted base; untested paths are unsound paths.

**Test:**

```python
kernel = TrustedKernel(frame=KripkeFrame.GL())

# Modus ponens
ab   = kernel.check_axiom(A >> B)
a    = kernel.check_axiom(A)
b    = kernel.verify_rule(ModusPonens, ab, a)
assert b.formula == B

# Necessitation refused on hypothesis
hyp  = ProofTerm(A, derivation=None, source='hypothesis')
with pytest.raises(NecessitationError):
    kernel.necessitate(hyp)

# Necessitation accepted on theorem
thm  = kernel.check_axiom(A >> A)  # tautology
box  = kernel.necessitate(thm)
assert box.formula == Box(A >> A)
```

---

### Layer 4: Proof Context and Search (`sympy_modal/context.py`)

The stateful proof environment that accumulates lemmas, manages hypotheses, and orchestrates proof search.

**Requirements:**

- Implement `ProofContext(frame, axioms)` with:
  - `ctx.assume(formula)` → `ProofTerm` with source `'hypothesis'`; adds to open hypothesis set
  - `ctx.discharge(hypothesis, proof_term)` → `ProofTerm` for the implication; removes hypothesis from open set
  - `ctx.apply(rule, *premises)` → `ProofTerm`; delegates to kernel
  - `ctx.necessitate(proof_term)` → `ProofTerm`; delegates to kernel with hygiene check
  - `ctx.lemma(name, proof_term)` → registers a proved theorem for reuse
  - `ctx.prove(formula, strategy=None)` → `ProofTerm | ProofFailure`; attempts proof search
  - `ctx.save()` / `ctx.restore()` → checkpoint and rollback for backtracking proof search
- Implement `ProofFailure(formula, obstacle, missing_axioms)` — a precise account of why proof search failed, including which additional axioms would suffice
- Implement at minimum the following proof search strategies:
  - `Strategy.Backward` — goal-directed backward chaining
  - `Strategy.ForwardChain` — forward chaining from hypotheses
  - `Strategy.ModalInduction` — specialised for modal fixed-point arguments

**Test:**

```python
ctx = ProofContext(frame=KripkeFrame.GL())

# Löb's theorem should be provable in GL
P   = PredicateVariable('P', type=FunctionType(Universe(0), BoolType()))
lob = ForAllPredicates(P,
        ctx.Box(ctx.Box(P) >> P) >> ctx.Box(P)
      )

result = ctx.prove(lob)
assert isinstance(result, ProofTerm)
assert result.formula == lob
assert ctx.kernel.check_term(result.term, lob)  # kernel verifies certificate

# Same theorem should fail in S4
ctx_s4 = ProofContext(frame=KripkeFrame.S4())
failure = ctx_s4.prove(lob)
assert isinstance(failure, ProofFailure)
assert Axiom.Lob in failure.missing_axioms
```

---

### Layer 5: Formalisation Interface (`sympy_modal/formalise.py`)

The natural language front end. This is explicitly the hardest layer and the one most likely to be incomplete in a first implementation. Design it for extensibility.

**Requirements:**

- Implement `FormalisationInterface` with:
  - `fi.resolve_modality(text)` → `ModalSignature` identifying which operators are present and which frame they require. Must distinguish at minimum: alethic (S5), deontic (D), epistemic (K45), provability (GL), temporal (S4). Ambiguous cases must return `AmbiguousModalityError` with the candidate readings listed explicitly.
  - `fi.resolve_quantifier_scope(text)` → `ScopeResolution` distinguishing de re from de dicto readings. For "necessarily someone wins": must return both `□∃x Wins(x)` (de dicto) and `∃x □Wins(x)` (de re) as distinct candidates with a flag indicating they are not equivalent in any normal modal logic.
  - `fi.resolve_order(text)` → `QuantifierOrder` identifying whether quantification is first-order (over individuals) or second-order (over predicates). Ambiguous cases must be flagged.
  - `fi.infer_frame(modal_signature)` → `KripkeFrame`; selects the most conservative frame consistent with the detected modality
  - `fi.formalise(text, frame=None)` → `ModalFormula | FormalisationError`; composes the above into a complete formalisation or returns a structured error with the precise point of failure

- This layer must fail loudly and precisely rather than silently produce incorrect formalisations. A wrong formalisation that passes into the proof kernel is worse than a formalisation error.
- Stub implementations returning `NotImplemented` with a descriptive message are acceptable for complex cases in a first version. Mark all stubs with `# TODO: requires LLM integration` so a downstream system can identify where LLM assistance is needed.

---

## Implementation Constraints

**Language and dependencies:**

- Python 3.11+
- SymPy 1.13+ as the base; extend, do not replace
- No dependencies on existing proof assistants (Lean, Rocq) in the core; these may be added as optional backends in a later phase
- `pytest` for testing; `mypy` for type checking; all public interfaces must be fully typed

**Code organisation:**

```
sympy_modal/
    __init__.py
    types.py        # Layer 1
    frames.py       # Layer 2
    kernel.py       # Layer 3 — trusted base; keep small
    context.py      # Layer 4
    formalise.py    # Layer 5
    operators.py    # Box, Diamond, ForAllPredicates, 
                    # ExistsPredicates, FixedPoint
    errors.py       # All custom exceptions
    tests/
        test_types.py
        test_frames.py
        test_kernel.py    # Must achieve 100% coverage
        test_context.py
        test_formalise.py
        test_integration.py
```

**Quality requirements:**

- The kernel (Layer 3) must achieve 100% test coverage; all other layers 80% minimum
- All public methods must have docstrings stating: the logical rule or semantic condition being implemented, preconditions, postconditions, and which layer of the Curry–Howard correspondence the method inhabits
- `mypy --strict` must pass on all files
- Every `FrameViolationError`, `NecessitationError`, `InvalidInferenceError`, and `FormalisationError` must include a human-readable explanation of the precise logical condition that was violated, suitable for display to a user who understands modal logic but may not know the implementation

---

## Validation: The Löb Integration Test

The following integration test is the primary acceptance criterion. It must pass before the implementation is considered complete:

```python
from sympy_modal import ProofContext, KripkeFrame, PredicateVariable
from sympy_modal import ForAllPredicates, Box, Universe, FunctionType, BoolType
from sympy_modal.errors import ProofFailure
from sympy import symbols

x = symbols('x')

# Löb's theorem in GL
ctx = ProofContext(frame=KripkeFrame.GL())
P   = PredicateVariable('P', type=FunctionType(Universe(0), BoolType()))

lob = ForAllPredicates(P,
        Box(Box(P(x)) >> P(x)) >> Box(P(x))
      )

proof = ctx.prove(lob)

# Certificate exists
assert proof.is_valid()

# Certificate is independently verifiable by kernel alone
assert ctx.kernel.check_term(proof.term, proof.formula)

# Certificate fails in S4 — wrong frame
ctx_s4  = ProofContext(frame=KripkeFrame.S4())
failure = ctx_s4.prove(lob)
assert isinstance(failure, ProofFailure)
assert failure.missing_axioms  # precise statement of what is missing

# The proof term from GL is rejected by the S4 kernel
assert not ctx_s4.kernel.check_term(proof.term, proof.formula)
```

---

## What This System Is Not

Do not implement:

- A full dependent type theory (Martin-Löf ITT or Calculus of Constructions in full generality) — implement only the fragment required for second-order modal logic
- A general-purpose proof assistant — this is a SymPy extension, not a replacement for Rocq or Lean
- An LLM integration in the core — the formalisation interface (Layer 5) should be designed so an LLM can populate it, but the LLM call itself is out of scope for this implementation
- Completeness for classical second-order logic — the target is the intuitionistic fragment, which is computationally tractable and sufficient for the Curry–Howard correspondence

---

## Suggested Implementation Order

1. `errors.py` — define all exceptions first; they are referenced everywhere
2. `types.py` — the type universe; nothing works without this
3. `operators.py` — Box, Diamond, ForAllPredicates, FixedPoint as SymPy-compatible symbolic objects
4. `frames.py` — Kripke frames; depends on operators
5. `kernel.py` — trusted kernel; depends on types and frames; write tests in parallel
6. `context.py` — proof context and search; depends on kernel
7. `formalise.py` — formalisation interface; depends on all layers; implement stubs first
8. `tests/test_integration.py` — Löb integration test; the acceptance criterion

Do not proceed to the next layer until the current layer's tests pass. The kernel in particular must be complete and fully tested before the context layer is begun, since the context layer's correctness depends entirely on the kernel's soundness.
