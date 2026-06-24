# SymPy Modal

`sympy_modal` is a second-order modal predicate calculus extension to SymPy. It provides a research-grade system combining a computer algebra interface with a proof-theoretic kernel.

## Reference Documentation

### Core Types & Universe
The type system replaces SymPy's untyped predicate logic with a cumulative type hierarchy.
* `Universe(level: int)`: Defines a universe level. `Universe(0)` is predicative.
* `BoolType()`: The type of standard Boolean propositions.
* `FunctionType(domain, codomain)`: Represents functions (e.g. predicates are functions returning `BoolType()`).
* `PredicateVariable(name, type)`: A second-order variable ranging over given types.
* `ModalPredicate(name, type, args)`: A typed predicate applied to arguments whose validity depends on the frame.
* `GuardedFixedPoint(operator, name)`: Well-founded fixed-point constructors for modal arguments (Löb's theorem).

### Logical Operators
`sympy_modal` provides the following logical primitives:
* `Box(formula)`: Necessity (□).
* `Diamond(formula)`: Possibility (◇). Equivalent to `~Box(~P)`.
* `ForAllPredicates(variable, formula)`: Second-order universal quantification.
* `ExistsPredicates(variable, formula)`: Second-order existential quantification.

For parsing ease, there are named modalities that imply their logic frame automatically via the `FormalisationInterface`:
* `ProvabilityBox` (GL)
* `AlethicBox` (S5)
* `EpistemicBox` (K45)
* `DeonticBox` (D)
* `TemporalBox` (S4)

### Kripke Frame Semantics
* `KripkeFrame(worlds, accessibility, axioms)`: A structural model to evaluate frames and model valid inferences across possible worlds. Includes static helpers like: `KripkeFrame.K()`, `KripkeFrame.T()`, `KripkeFrame.S4()`, `KripkeFrame.S5()`, `KripkeFrame.GL()`, `KripkeFrame.D()`, and `KripkeFrame.K45()`.
* `.validates(formula)`: Verifies if a formula evaluates to True in all accessible worlds.

### Trusted Kernel
`TrustedKernel(frame)` enforces strict Intuitionistic Natural Deduction without compromise. All inferences strictly require a valid `ProofTerm`.
* `.check_axiom(formula)`: Verifies the axiom exists within the logic frame.
* `.verify_rule(rule, premises)`: Validates `ModusPonens`, `Necessitation`, etc.

### Proof Context
`ProofContext(frame)` encapsulates and coordinates active hypotheses, discharging rules, and the proof search system.
* `.assume(formula)`: Assumes a hypothesis. Warns on classical axioms like Law of Excluded Middle.
* `.discharge(hypothesis, pt)`: Removes a hypothesis to yield an Implication.
* `.prove(formula, strategy)`: Leverages search strategies like `Strategy.ForwardChain`, `Strategy.Backward`, or `Strategy.ModalInduction`.

### Formalisation Interface
`FormalisationInterface()` connects natural human inputs/code strings seamlessly into verified modal terms.
* `.formalise(code_str)`: Parses string representations of logic natively.
* `.resolve_modality(code_str)`: Uses named modal operators (e.g. `AlethicBox`) or heuristic axiom analysis to deduce the logic frame.

---

## Tutorial and Examples

Here are 5 fully worked-out examples that demonstrate the system in action:

### Example 1: Basic Validity in S4
We can check if an axiom holds in a given Kripke frame. For S4, the frame is reflexive and transitive. This means `□P → P` and `□P → □□P` are valid.
```python
from sympy import Symbol, Implies
from sympy_modal import KripkeFrame, Box

s4 = KripkeFrame.S4()
p = Symbol('p')

# Reflexivity: □P → P
assert s4.validates(Implies(Box(p), p))

# Transitivity: □P → □□P
assert s4.validates(Implies(Box(p), Box(Box(p))))
```

### Example 2: The Trusted Proof Kernel
The `TrustedKernel` enforces rigorous intuitionistic natural deduction without compromising logic.
```python
from sympy import Symbol, Implies
from sympy_modal import KripkeFrame, TrustedKernel, ProofTerm, ModusPonens

kernel = TrustedKernel(frame=KripkeFrame.GL())
p = Symbol('p')
q = Symbol('q')

# Check axiom validity and apply rule (use evaluate=False so SymPy doesn't resolve to True if identical)
pt_ab = ProofTerm(Implies(p, q, evaluate=False), source="axiom")
pt_a = ProofTerm(p)

# Modus Ponens verification
pt_b = kernel.verify_rule(ModusPonens, [pt_ab, pt_a])
assert pt_b.formula == q
```

### Example 3: Proof Search Context
The `ProofContext` manages active hypotheses and applies natural deduction implicitly.
```python
from sympy import Symbol, Implies
from sympy_modal import ProofContext, KripkeFrame, ProofTerm

ctx = ProofContext(KripkeFrame.K())
p = Symbol('p')
q = Symbol('q')

hyp = ctx.assume(p)
assert hyp.source == 'hypothesis'

# Create a fake derivation from the hypothesis for demonstration
pt_q = ProofTerm(q, hypotheses=[p])

# Discharging the hypothesis yields P -> Q
impl = ctx.discharge(hyp, pt_q)
assert impl.formula == Implies(p, q, evaluate=False)
```

### Example 4: Löb's Theorem in GL
Löb's theorem `∀P (□(□P → P) → □P)` is valid in Gödel-Löb provability logic but not in standard alethic/temporal logics.
```python
from sympy import symbols, Implies
from sympy_modal import (
    ProofContext, KripkeFrame, PredicateVariable,
    ForAllPredicates, Box, Universe, FunctionType, BoolType
)

x = symbols('x')
ctx = ProofContext(frame=KripkeFrame.GL())
P = PredicateVariable('P', type=FunctionType(Universe(0), BoolType()))

lob = ForAllPredicates(P,
    Implies(Box(Implies(Box(P(x)), P(x))), Box(P(x)))
)

# Certificate proves the theorem is valid under GL frame
proof = ctx.prove(lob)
assert proof.is_valid
```

### Example 5: Formalisation Interface
You can parse string representations of Modal formulas seamlessly.
```python
from sympy_modal import FormalisationInterface, AlethicBox

fi = FormalisationInterface()
expr_str = "AlethicBox(Symbol('p'))"
formula = fi.formalise(expr_str)

assert isinstance(formula, AlethicBox)

# The interface can infer the correct frame automatically
sig = fi.resolve_modality(expr_str)
assert "alethic" in sig.operators
```

## Overview

The `sympy_modal/` directory contains the core implementation of the second-order modal logic extension for SymPy. Together with 8 test files, the codebase consists of 1534 lines of code. Below is a one-paragraph summary of each file, followed by its corresponding declarations.

### `sympy_modal/__init__.py`
This initialization file exposes the public API for the `sympy_modal` package, ensuring that core components like operators, frames, the trusted kernel, and proof contexts are accessible directly when importing the module.

**Declarations:**
- `__all__` (GlobalVar): List of publicly exported symbols from the sympy_modal package.

### `sympy_modal/formalise.py`
This module provides the `FormalisationInterface`, which is responsible for parsing string-based modal logic representations into typed SymPy objects. It infers Kripke frames heuristically, manages quantifier scope resolution, and validates modality mappings.

**Declarations:**
- `ModalSignature` (Class): Dataclass storing the parsed logic frame signatures and bound variables.
- `ScopeResolution` (Class): Enum defining how quantifiers scope over expressions (e.g., local vs. global).
- `QuantifierOrder` (Class): Enum tracking whether universal or existential quantifiers appear first.
- `FormalisationInterface` (Class): Main interface for parsing string representations of modal logic into AST objects.
- `FormalisationInterface.__init__` (Method): Initializes the interface with an empty context dict.
- `FormalisationInterface._get_local_dict` (Method): Retrieves the local dictionary populated with SymPy modal classes.
- `FormalisationInterface.parse_code` (Method): Evaluates Python code strings into internal formula objects.
- `FormalisationInterface.resolve_modality` (Method): Heuristically infers the correct logic frame based on explicit operators or axioms.
- `FormalisationInterface.resolve_quantifier_scope` (Method): Resolves standard bindings and checks for ill-scoped modal variables.
- `FormalisationInterface.resolve_order` (Method): Ensures quantifiers follow a valid scoping sequence to avoid invalid dependencies.
- `FormalisationInterface.infer_frame` (Method): Internal helper to map detected axioms or operators to specific Kripke frames.
- `FormalisationInterface.formalise` (Method): Main method mapping code strings to validated formula expressions within an inferred frame.

### `sympy_modal/frames.py`
This file defines the structural semantics of the logic via `KripkeFrame` and `Axiom` definitions. It implements the core accessibility evaluation logic across possible worlds for various standard modalities like S4, S5, and GL.

**Declarations:**
- `Axiom` (Class): Enum for standard modal logic axioms like K, T, 4, 5, etc.
- `Axiom.K` (ClassVar): Axiom K: □(P → Q) → (□P → □Q).
- `Axiom.T` (ClassVar): Axiom T: □P → P (reflexivity).
- `Axiom.B` (ClassVar): Axiom B: P → □◇P (symmetry).
- `Axiom.Four` (ClassVar): Axiom 4: □P → □□P (transitivity).
- `Axiom.Five` (ClassVar): Axiom 5: ◇P → □◇P (Euclidean).
- `Axiom.Lob` (ClassVar): Axiom L: □(□P → P) → □P (Löb's theorem).
- `Axiom.D` (ClassVar): Axiom D: □P → ◇P (seriality).
- `KripkeFrame` (Class): Represents the structural model of possible worlds and accessibility relations.
- `KripkeFrame.__init__` (Method): Initializes a frame with specified worlds, relation matrix, and enforced axioms.
- `KripkeFrame.K` (Method): Factory for standard modal logic K (no special relations).
- `KripkeFrame.T` (Method): Factory for logic T (reflexive).
- `KripkeFrame.S4` (Method): Factory for logic S4 (reflexive, transitive).
- `KripkeFrame.S5` (Method): Factory for logic S5 (reflexive, symmetric, transitive).
- `KripkeFrame.GL` (Method): Factory for Provability Logic (transitive, converse well-founded).
- `KripkeFrame.D` (Method): Factory for Deontic Logic (serial).
- `KripkeFrame.K45` (Method): Factory for Epistemic Logic (transitive, Euclidean).
- `KripkeFrame._evaluate` (Method): Internal logic engine to recursively evaluate formula truth values per world.
- `KripkeFrame.validates` (Method): Checks if a formula evaluates to True across all accessible worlds.
- `KripkeFrame._is_axiom_instance` (Method): Checks syntactically if an expression is an instantiation of a known axiom.
- `KripkeFrame.is_valid_inference` (Method): Returns True if a formula structurally mirrors the frame's defining axioms.

### `sympy_modal/context.py`
This module handles the stateful proof environment via the `ProofContext` class. It manages active hypotheses, provides discharging mechanisms for implications, and orchestrates automated proof search strategies (forward, backward, and inductive).

**Declarations:**
- `Strategy` (Class): Enum detailing automated proof search strategies.
- `Strategy.Backward` (ClassVar): Backward chaining strategy (goal-directed).
- `Strategy.ForwardChain` (ClassVar): Forward chaining strategy (data-directed).
- `Strategy.ModalInduction` (ClassVar): Specialized inductive strategy for well-founded fixed points.
- `ProofContext` (Class): Environment coordinating active hypotheses and proof steps.
- `ProofContext.__init__` (Method): Initializes the context using a specific Kripke frame.
- `ProofContext.assume` (Method): Introduces a new assumption/hypothesis to the current scope.
- `ProofContext.discharge` (Method): Discharges an assumption, returning a validated Implication.
- `ProofContext.apply` (Method): Applies a logical rule or theorem to known terms.
- `ProofContext.necessitate` (Method): Applies the Necessitation rule via the trusted kernel.
- `ProofContext.lemma` (Method): Registers an intermediate proven step for future use.
- `ProofContext.save` (Method): Snapshots the current context state for backtracking.
- `ProofContext.restore` (Method): Reverts to a previously saved context state.
- `ProofContext._warn_classical` (Method): Issues warnings when classical axioms (like LEM) are invoked.
- `ProofContext.prove` (Method): Automatically attempts to derive a formula using a designated strategy.
- `ProofContext._forward_chain` (Method): Internal automated forward reasoning loop.
- `ProofContext._backward_chain` (Method): Internal automated backward reasoning loop.
- `ProofContext._modal_induction` (Method): Internal routine handling inductive fixed point proofs.

### `sympy_modal/types.py`
This file establishes the foundational cumulative type hierarchy that replaces standard SymPy's untyped logic. It introduces universes, function mapping types, and guarded fixed-point constructors essential for valid second-order modal formulations.

**Declarations:**
- `Type` (Class): Base class for the cumulative type hierarchy.
- `Universe` (Class): Represents a specific foundational universe level.
- `Universe.__new__` (Method): Creates a universe instance for a given integer level.
- `Universe.level` (Method): Returns the integer level of the universe.
- `Universe.__str__` (Method): String representation of the universe.
- `Universe.__repr__` (Method): Detailed representation of the universe.
- `BoolType` (Class): Represents standard Boolean propositions.
- `BoolType.__new__` (Method): Instantiates the singleton BoolType.
- `BoolType.__str__` (Method): String representation.
- `BoolType.__repr__` (Method): Detailed representation.
- `FunctionType` (Class): Represents mapping types between domains and codomains.
- `FunctionType.__new__` (Method): Creates a function type.
- `FunctionType.domain` (Method): Returns the input type (domain).
- `FunctionType.codomain` (Method): Returns the output type (codomain).
- `FunctionType.__str__` (Method): String representation.
- `FunctionType.__repr__` (Method): Detailed representation.
- `TypedSymbol` (Class): A variable carrying explicit second-order type information.
- `TypedSymbol.__new__` (Method): Instantiates a typed symbol.
- `TypedSymbol.type` (Method): Returns the symbol's assigned type.
- `PredicateVariable` (Class): A second-order variable ranging over specified types.
- `PredicateVariable.__new__` (Method): Creates a predicate variable.
- `PredicateVariable.__call__` (Method): Applies the predicate variable to arguments.
- `ModalPredicate` (Class): A typed predicate whose evaluation is frame-dependent.
- `ModalPredicate.__new__` (Method): Creates a modal predicate with a name, type, and arguments.
- `ModalPredicate.name` (Method): Returns the predicate name.
- `ModalPredicate.type` (Method): Returns the predicate type.
- `GuardedFixedPoint` (Class): Constructor for well-founded fixed points in modal logic.
- `GuardedFixedPoint.__new__` (Method): Creates a guarded fixed point instance.
- `GuardedFixedPoint.operator` (Method): Returns the bounding modal operator.
- `GuardedFixedPoint.name` (Method): Returns the fixed point name.

### `sympy_modal/kernel.py`
This module implements the `TrustedKernel`, a strict proof-theoretic engine enforcing intuitionistic natural deduction. It ensures all derivations, like `ModusPonens` or `Necessitation`, are cryptographically sound by evaluating `ProofTerm` instances.

**Declarations:**
- `ProofTerm` (Class): Represents a verified derivation step in the calculus.
- `ProofTerm.__init__` (Method): Records the formula, source rule, and active hypotheses.
- `ProofTerm.is_valid` (Method): Property confirming the term was generated by the TrustedKernel.
- `ModusPonens` (Class): Rule of inference: From P and P → Q, derive Q.
- `TrustedKernel` (Class): The core proof engine enforcing intuitionistic deduction.
- `TrustedKernel.__init__` (Method): Binds the kernel to a specific evaluating KripkeFrame.
- `TrustedKernel.check_axiom` (Method): Verifies that an axiom belongs to the bound frame.
- `TrustedKernel.verify_rule` (Method): Executes logical rules strictly over provided valid ProofTerms.
- `TrustedKernel.necessitate` (Method): Applies the Necessitation rule securely.
- `TrustedKernel.check_term` (Method): Final safety check validating a completed ProofTerm.

### `sympy_modal/operators.py`
This file defines the abstract syntax tree nodes for logical operators. It includes the foundational necessity (`Box`) and possibility (`Diamond`) primitives, alongside specialized named modalities (e.g., `AlethicBox`) and second-order quantifiers.

**Declarations:**
- `ModalOperator` (Class): Base class for modal logic operators.
- `ModalOperator.__new__` (Method): Instantiates an operator node.
- `ModalOperator.eval` (Method): Base evaluation method for the operator.
- `Box` (Class): The Necessity operator (□).
- `Box.__new__` (Method): Instantiates a Box operator wrapping a formula.
- `Box.modality` (Method): Returns the modality type.
- `Box.is_well_typed` (Method): Checks if the internal formula is a valid proposition.
- `Diamond` (Class): The Possibility operator (◇).
- `Diamond.__new__` (Method): Instantiates a Diamond operator wrapping a formula.
- `Diamond.modality` (Method): Returns the modality type.
- `Diamond.is_well_typed` (Method): Checks if the internal formula is a valid proposition.
- `ProvabilityBox` (Class): A necessity operator hardcoded to Provability logic (GL).
- `ProvabilityBox.__new__` (Method): Instantiates a ProvabilityBox.
- `AlethicBox` (Class): A necessity operator hardcoded to Alethic logic (S5).
- `AlethicBox.__new__` (Method): Instantiates an AlethicBox.
- `EpistemicBox` (Class): A necessity operator hardcoded to Epistemic logic (K45).
- `EpistemicBox.__new__` (Method): Instantiates an EpistemicBox.
- `DeonticBox` (Class): A necessity operator hardcoded to Deontic logic (D).
- `DeonticBox.__new__` (Method): Instantiates a DeonticBox.
- `TemporalBox` (Class): A necessity operator hardcoded to Temporal logic (S4).
- `TemporalBox.__new__` (Method): Instantiates a TemporalBox.
- `ForAllPredicates` (Class): Second-order universal quantifier.
- `ForAllPredicates.__new__` (Method): Instantiates ForAllPredicates over a variable and formula.
- `ForAllPredicates.variable` (Method): Returns the bound variable.
- `ForAllPredicates.formula` (Method): Returns the quantified formula.
- `ForAllPredicates.is_well_typed` (Method): Ensures the quantification is typed correctly.
- `ExistsPredicates` (Class): Second-order existential quantifier.
- `ExistsPredicates.__new__` (Method): Instantiates ExistsPredicates over a variable and formula.
- `ExistsPredicates.variable` (Method): Returns the bound variable.
- `ExistsPredicates.formula` (Method): Returns the quantified formula.
- `ExistsPredicates.is_well_typed` (Method): Ensures the quantification is typed correctly.

### `sympy_modal/errors.py`
This module centralizes all custom exceptions raised within the package, ranging from structural `FrameViolationError`s and invalid deduction checks to parsing failures during formalisation.

**Declarations:**
- `SymPyModalError` (Class): Base exception for sympy_modal errors.
- `FrameViolationError` (Class): Raised when a formula violates the enforced Kripke frame.
- `FrameViolationError.__init__` (Method): Initializes with specific violation details.
- `NecessitationError` (Class): Raised when the Necessitation rule is applied invalidly.
- `NecessitationError.__init__` (Method): Initializes with necessitation failure details.
- `InvalidInferenceError` (Class): Raised on an invalid deduction step.
- `InvalidInferenceError.__init__` (Method): Initializes with deduction context.
- `FormalisationError` (Class): Raised when code strings cannot be mapped to modal terms.
- `FormalisationError.__init__` (Method): Initializes with parse failure information.
- `NotAnAxiomError` (Class): Raised when an expression claims to be an axiom but isn't.
- `NotAnAxiomError.__init__` (Method): Initializes with the invalid axiom claim.
- `AmbiguousModalityError` (Class): Raised when a frame cannot be heuristically inferred.
- `AmbiguousModalityError.__init__` (Method): Initializes with ambiguity details.
- `ProofFailure` (Class): Raised when the automated prover fails.
- `ProofFailure.__init__` (Method): Initializes with failure context.
- `ProofFailure.__repr__` (Method): Detailed representation of the proof failure.


## Proposals for Enhancing the Reference and Tutorial Examples

### Enhancing Example Output
Currently, the tutorial examples in this document rely on silent assertions (e.g., `assert s4.validates(...)`). To provide a better learning experience, the examples should be updated to produce meaningful console output. Instead of just passing silently, we propose replacing assertions with print statements that evaluate to boolean conditions or print formulas directly. For example, changing `assert s4.validates(Implies(Box(p), p))` to `print(f"S4 Reflexivity validation: {s4.validates(Implies(Box(p), p))}")` will output `S4 Reflexivity validation: True`. This provides users with direct, observable feedback when they copy, paste, and run the tutorial scripts locally.

### Expanding Example Explanations
To elevate the tutorial section from a mere list of code snippets to a comprehensive learning resource, each example should be expanded to include a structured, two-paragraph introduction and explanation.
- The **first paragraph** should introduce the specific logical concept or system component being demonstrated (e.g., explaining the intuitive meaning behind the Kripke accessibility relations required for S4).
- The **second paragraph** should explicitly walk through the code, explaining step-by-step how the components interact (e.g., explaining why `evaluate=False` is necessary, or how the `ProofContext` manages active hypotheses under the hood). This deeper context will bridge the gap between abstract modal logic theory and the concrete SymPy API implementation.
