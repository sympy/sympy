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

In modal logic, different systems are characterized by the properties of their accessibility relations between possible worlds. The S4 system is defined by Kripke frames that are both reflexive and transitive. Reflexivity ensures that whatever is necessary in a given world must be true in that same world (□P → P), while transitivity guarantees that if something is necessarily true, it is necessarily necessarily true (□P → □□P). This system forms the foundation for many temporal and epistemic logical structures.

The code below demonstrates how to programmatically construct an S4 Kripke frame and verify these core axioms using SymPy. We first instantiate the `s4` frame and a basic boolean proposition `p`. We then use the `.validates()` method to evaluate implications across the frame. By wrapping these checks in `print` statements, we can directly observe that the frame confirms both the reflexivity and transitivity axioms as valid inferences.

```python
from sympy import Symbol, Implies
from sympy_modal import KripkeFrame, Box

s4 = KripkeFrame.S4()
p = Symbol('p')

# Reflexivity: □P → P
print(f"S4 Reflexivity validation: {s4.validates(Implies(Box(p), p))}")
# Output: S4 Reflexivity validation: True

# Transitivity: □P → □□P
print(f"S4 Transitivity validation: {s4.validates(Implies(Box(p), Box(Box(p))))}")
# Output: S4 Transitivity validation: True
```

### Example 2: The Trusted Proof Kernel

At the core of the proof-theoretic system is the `TrustedKernel`, which enforces strict intuitionistic natural deduction. Rather than trusting human input or heuristic checks, the kernel demands cryptographically sound derivations for every logical step. It operates by verifying instances of `ProofTerm`, ensuring that logical operations such as Modus Ponens or Necessitation are only applied to valid premises within the bounds of a specified Kripke frame.

In this example, we configure the `TrustedKernel` within the context of Provability Logic (GL). We manually construct two `ProofTerm` instances: one representing an implication axiom (`P → Q`) and another representing the premise (`P`). It is crucial to use `evaluate=False` when constructing the implication so SymPy maintains the structural integrity of the formula instead of eagerly resolving it. The kernel then strictly verifies the application of the `ModusPonens` rule against these terms, returning a new derived `ProofTerm` which we print to verify its formula.

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
print(f"Derived formula via Modus Ponens: {pt_b.formula}")
# Output: Derived formula via Modus Ponens: q
```

### Example 3: Proof Search Context

The `ProofContext` serves as an intelligent environment that manages the stateful progression of formal proofs. It orchestrates active hypotheses, keeps track of dependencies, and facilitates automated proof search strategies. A fundamental operation in natural deduction managed by the context is assumption and discharging: temporarily assuming a hypothesis to derive a consequence, and then formally removing the assumption to yield a proven implication.

This script demonstrates how to utilize `ProofContext` initialized with the basic modal logic frame K. We introduce `p` as an active hypothesis using `ctx.assume(p)`. After simulating a derivation that relies on this hypothesis resulting in `q`, we use the `discharge` method. This action formally removes the hypothesis from the active context and correctly binds it into an implication (`P → Q`). Printing the results clearly illustrates how the context transitions from managing an assumption to finalizing a logical statement.

```python
from sympy import Symbol, Implies
from sympy_modal import ProofContext, KripkeFrame, ProofTerm

ctx = ProofContext(KripkeFrame.K())
p = Symbol('p')
q = Symbol('q')

hyp = ctx.assume(p)
print(f"Hypothesis source: {hyp.source}")
# Output: Hypothesis source: hypothesis

# Create a fake derivation from the hypothesis for demonstration
pt_q = ProofTerm(q, hypotheses=[p])

# Discharging the hypothesis yields P -> Q
impl = ctx.discharge(hyp, pt_q)
print(f"Discharged implication: {impl.formula}")
# Output: Discharged implication: Implies(p, q)
```

### Example 4: Löb's Theorem in GL

Löb's theorem is a profound result in provability logic, expressing that if a formal system can prove that its own provability of a statement implies the statement itself, then it can just prove the statement unconditionally. Formally, it is written as `∀P (□(□P → P) → □P)`. While this holds true under the specific conditions of Gödel-Löb (GL) provability logic—characterized by transitive and converse well-founded accessibility relations—it fails in standard alethic or temporal frameworks.

In the provided example, we leverage `sympy_modal`'s typed second-order capabilities to formalize Löb's theorem. We define a `PredicateVariable` `P` that strictly maps elements of a base universe to boolean propositions. Constructing the formula requires precisely nesting the `Box` operators and implications inside a universal quantifier. When we pass this constructed formula to `ctx.prove(lob)` configured for a GL frame, the automated prover attempts to derive a valid proof certificate. The print output confirms whether the automated strategy successfully authenticated the theorem's validity.

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
print(f"Löb's theorem proved in GL: {proof.is_valid}")
# Output: Löb's theorem proved in GL: True
```

### Example 5: Formalisation Interface

Constructing deeply nested logical ASTs manually can be tedious and prone to syntax errors. The `FormalisationInterface` provides a streamlined bridge, allowing users to express complex modal formulas as natural string inputs. It processes code strings seamlessly into properly typed SymPy objects. Additionally, the interface possesses heuristic capabilities to resolve quantifier scopes and deduce appropriate underlying Kripke frames just by analyzing the named modal operators present in the text.

In this concluding example, we instantiate the `FormalisationInterface` and supply it with a string defining a proposition wrapped in an `AlethicBox`. The `formalise` method parses this string and returns the corresponding SymPy AST object. To showcase its heuristic power, we then call `resolve_modality`, which analyzes the string to deduce the logical signature. By printing the type of the formalized object and the contents of the detected signature, we can verify that the system correctly inferred the intention of utilizing an alethic (S5) structural frame.

```python
from sympy_modal import FormalisationInterface, AlethicBox

fi = FormalisationInterface()
expr_str = "AlethicBox(Symbol('p'))"
formula = fi.formalise(expr_str)

print(f"Parsed formula type: {type(formula).__name__}")
# Output: Parsed formula type: AlethicBox

# The interface can infer the correct frame automatically
sig = fi.resolve_modality(expr_str)
print(f"Inferred modality signatures: {sig.operators}")
# Output: Inferred modality signatures: ['alethic']
```

## Overview

The `sympy_modal/` directory contains the core implementation of the second-order modal logic extension for SymPy. Together with 8 test files, the codebase consists of 1534 lines of code. Below is a detailed breakdown.

### `sympy_modal/__init__.py`
(74 lines of code) This initialization file exposes the public API for the `sympy_modal` package, ensuring that core components like operators, frames, the trusted kernel, and proof contexts are accessible directly when importing the module.

**Declarations:**
- `__all__` (GlobalVar): List of publicly exported symbols from the sympy_modal package. (18 lines of code)

### `sympy_modal/formalise.py`
(171 lines of code) This module provides the `FormalisationInterface`, which is responsible for parsing string-based modal logic representations into typed SymPy objects. It infers Kripke frames heuristically, manages quantifier scope resolution, and validates modality mappings.

**Declarations:**
- `ModalSignature` (Class): Dataclass storing the parsed logic frame signatures and bound variables. (3 lines of code)
- `ScopeResolution` (Class): Enum defining how quantifiers scope over expressions (e.g., local vs. global). (3 lines of code)
- `QuantifierOrder` (Class): Enum tracking whether universal or existential quantifiers appear first. (3 lines of code)
- `FormalisationInterface` (Class): Main interface for parsing string representations of modal logic into AST objects. (136 lines of code)
- `FormalisationInterface.__init__` (Method): Initializes the interface with an empty context dict. (9 lines of code)
- `FormalisationInterface._get_local_dict` (Method): Retrieves the local dictionary populated with SymPy modal classes. (19 lines of code)
- `FormalisationInterface.parse_code` (Method): Evaluates Python code strings into internal formula objects. (8 lines of code)
- `FormalisationInterface.resolve_modality` (Method): Heuristically infers the correct logic frame based on explicit operators or axioms. (41 lines of code)
- `FormalisationInterface.resolve_quantifier_scope` (Method): Resolves standard bindings and checks for ill-scoped modal variables. (11 lines of code)
- `FormalisationInterface.resolve_order` (Method): Ensures quantifiers follow a valid scoping sequence to avoid invalid dependencies. (18 lines of code)
- `FormalisationInterface.infer_frame` (Method): Internal helper to map detected axioms or operators to specific Kripke frames. (5 lines of code)
- `FormalisationInterface.formalise` (Method): Main method mapping code strings to validated formula expressions within an inferred frame. (17 lines of code)

### `sympy_modal/frames.py`
(198 lines of code) This file defines the structural semantics of the logic via `KripkeFrame` and `Axiom` definitions. It implements the core accessibility evaluation logic across possible worlds for various standard modalities like S4, S5, and GL.

**Declarations:**
- `Axiom` (Class): Enum for standard modal logic axioms like K, T, 4, 5, etc. (8 lines of code)
- `Axiom.K` (ClassVar): Axiom K: □(P → Q) → (□P → □Q). (1 lines of code)
- `Axiom.T` (ClassVar): Axiom T: □P → P (reflexivity). (1 lines of code)
- `Axiom.B` (ClassVar): Axiom B: P → □◇P (symmetry). (1 lines of code)
- `Axiom.Four` (ClassVar): Axiom 4: □P → □□P (transitivity). (1 lines of code)
- `Axiom.Five` (ClassVar): Axiom 5: ◇P → □◇P (Euclidean). (1 lines of code)
- `Axiom.Lob` (ClassVar): Axiom L: □(□P → P) → □P (Löb's theorem). (1 lines of code)
- `Axiom.D` (ClassVar): Axiom D: □P → ◇P (seriality). (1 lines of code)
- `KripkeFrame` (Class): Represents the structural model of possible worlds and accessibility relations. (178 lines of code)
- `KripkeFrame.__init__` (Method): Initializes a frame with specified worlds, relation matrix, and enforced axioms. (4 lines of code)
- `KripkeFrame.K` (Method): Factory for standard modal logic K (no special relations). (3 lines of code)
- `KripkeFrame.T` (Method): Factory for logic T (reflexive). (3 lines of code)
- `KripkeFrame.S4` (Method): Factory for logic S4 (reflexive, transitive). (3 lines of code)
- `KripkeFrame.S5` (Method): Factory for logic S5 (reflexive, symmetric, transitive). (3 lines of code)
- `KripkeFrame.GL` (Method): Factory for Provability Logic (transitive, converse well-founded). (6 lines of code)
- `KripkeFrame.D` (Method): Factory for Deontic Logic (serial). (3 lines of code)
- `KripkeFrame.K45` (Method): Factory for Epistemic Logic (transitive, Euclidean). (3 lines of code)
- `KripkeFrame._evaluate` (Method): Internal logic engine to recursively evaluate formula truth values per world. (43 lines of code)
- `KripkeFrame.validates` (Method): Checks if a formula evaluates to True across all accessible worlds. (47 lines of code)
- `KripkeFrame._is_axiom_instance` (Method): Checks syntactically if an expression is an instantiation of a known axiom. (26 lines of code)
- `KripkeFrame.is_valid_inference` (Method): Returns True if a formula structurally mirrors the frame's defining axioms. (12 lines of code)

### `sympy_modal/context.py`
(257 lines of code) This module handles the stateful proof environment via the `ProofContext` class. It manages active hypotheses, provides discharging mechanisms for implications, and orchestrates automated proof search strategies (forward, backward, and inductive).

**Declarations:**
- `Strategy` (Class): Enum detailing automated proof search strategies. (4 lines of code)
- `Strategy.Backward` (ClassVar): Backward chaining strategy (goal-directed). (1 lines of code)
- `Strategy.ForwardChain` (ClassVar): Forward chaining strategy (data-directed). (1 lines of code)
- `Strategy.ModalInduction` (ClassVar): Specialized inductive strategy for well-founded fixed points. (1 lines of code)
- `ProofContext` (Class): Environment coordinating active hypotheses and proof steps. (236 lines of code)
- `ProofContext.__init__` (Method): Initializes the context using a specific Kripke frame. (13 lines of code)
- `ProofContext.assume` (Method): Introduces a new assumption/hypothesis to the current scope. (8 lines of code)
- `ProofContext.discharge` (Method): Discharges an assumption, returning a validated Implication. (14 lines of code)
- `ProofContext.apply` (Method): Applies a logical rule or theorem to known terms. (5 lines of code)
- `ProofContext.necessitate` (Method): Applies the Necessitation rule via the trusted kernel. (5 lines of code)
- `ProofContext.lemma` (Method): Registers an intermediate proven step for future use. (5 lines of code)
- `ProofContext.save` (Method): Snapshots the current context state for backtracking. (5 lines of code)
- `ProofContext.restore` (Method): Reverts to a previously saved context state. (6 lines of code)
- `ProofContext._warn_classical` (Method): Issues warnings when classical axioms (like LEM) are invoked. (23 lines of code)
- `ProofContext.prove` (Method): Automatically attempts to derive a formula using a designated strategy. (42 lines of code)
- `ProofContext._forward_chain` (Method): Internal automated forward reasoning loop. (50 lines of code)
- `ProofContext._backward_chain` (Method): Internal automated backward reasoning loop. (34 lines of code)
- `ProofContext._modal_induction` (Method): Internal routine handling inductive fixed point proofs. (10 lines of code)

### `sympy_modal/types.py`
(134 lines of code) This file establishes the foundational cumulative type hierarchy that replaces standard SymPy's untyped logic. It introduces universes, function mapping types, and guarded fixed-point constructors essential for valid second-order modal formulations.

**Declarations:**
- `Type` (Class): Base class for the cumulative type hierarchy. (3 lines of code)
- `Universe` (Class): Represents a specific foundational universe level. (23 lines of code)
- `Universe.__new__` (Method): Creates a universe instance for a given integer level. (7 lines of code)
- `Universe.level` (Method): Returns the integer level of the universe. (2 lines of code)
- `Universe.__str__` (Method): String representation of the universe. (2 lines of code)
- `Universe.__repr__` (Method): Detailed representation of the universe. (2 lines of code)
- `BoolType` (Class): Represents standard Boolean propositions. (10 lines of code)
- `BoolType.__new__` (Method): Instantiates the singleton BoolType. (2 lines of code)
- `BoolType.__str__` (Method): String representation. (2 lines of code)
- `BoolType.__repr__` (Method): Detailed representation. (2 lines of code)
- `FunctionType` (Class): Represents mapping types between domains and codomains. (20 lines of code)
- `FunctionType.__new__` (Method): Creates a function type. (4 lines of code)
- `FunctionType.domain` (Method): Returns the input type (domain). (2 lines of code)
- `FunctionType.codomain` (Method): Returns the output type (codomain). (2 lines of code)
- `FunctionType.__str__` (Method): String representation. (2 lines of code)
- `FunctionType.__repr__` (Method): Detailed representation. (2 lines of code)
- `TypedSymbol` (Class): A variable carrying explicit second-order type information. (10 lines of code)
- `TypedSymbol.__new__` (Method): Instantiates a typed symbol. (4 lines of code)
- `TypedSymbol.type` (Method): Returns the symbol's assigned type. (2 lines of code)
- `PredicateVariable` (Class): A second-order variable ranging over specified types. (12 lines of code)
- `PredicateVariable.__new__` (Method): Creates a predicate variable. (4 lines of code)
- `PredicateVariable.__call__` (Method): Applies the predicate variable to arguments. (3 lines of code)
- `ModalPredicate` (Class): A typed predicate whose evaluation is frame-dependent. (18 lines of code)
- `ModalPredicate.__new__` (Method): Creates a modal predicate with a name, type, and arguments. (5 lines of code)
- `ModalPredicate.name` (Method): Returns the predicate name. (2 lines of code)
- `ModalPredicate.type` (Method): Returns the predicate type. (2 lines of code)
- `GuardedFixedPoint` (Class): Constructor for well-founded fixed points in modal logic. (17 lines of code)
- `GuardedFixedPoint.__new__` (Method): Creates a guarded fixed point instance. (4 lines of code)
- `GuardedFixedPoint.operator` (Method): Returns the bounding modal operator. (2 lines of code)
- `GuardedFixedPoint.name` (Method): Returns the fixed point name. (2 lines of code)

### `sympy_modal/kernel.py`
(97 lines of code) This module implements the `TrustedKernel`, a strict proof-theoretic engine enforcing intuitionistic natural deduction. It ensures all derivations, like `ModusPonens` or `Necessitation`, are cryptographically sound by evaluating `ProofTerm` instances.

**Declarations:**
- `ProofTerm` (Class): Represents a verified derivation step in the calculus. (15 lines of code)
- `ProofTerm.__init__` (Method): Records the formula, source rule, and active hypotheses. (6 lines of code)
- `ProofTerm.is_valid` (Method): Property confirming the term was generated by the TrustedKernel. (3 lines of code)
- `ModusPonens` (Class): Rule of inference: From P and P → Q, derive Q. (3 lines of code)
- `TrustedKernel` (Class): The core proof engine enforcing intuitionistic deduction. (65 lines of code)
- `TrustedKernel.__init__` (Method): Binds the kernel to a specific evaluating KripkeFrame. (2 lines of code)
- `TrustedKernel.check_axiom` (Method): Verifies that an axiom belongs to the bound frame. (10 lines of code)
- `TrustedKernel.verify_rule` (Method): Executes logical rules strictly over provided valid ProofTerms. (26 lines of code)
- `TrustedKernel.necessitate` (Method): Applies the Necessitation rule securely. (12 lines of code)
- `TrustedKernel.check_term` (Method): Final safety check validating a completed ProofTerm. (7 lines of code)

### `sympy_modal/operators.py`
(130 lines of code) This file defines the abstract syntax tree nodes for logical operators. It includes the foundational necessity (`Box`) and possibility (`Diamond`) primitives, alongside specialized named modalities (e.g., `AlethicBox`) and second-order quantifiers.

**Declarations:**
- `ModalOperator` (Class): Base class for modal logic operators. (9 lines of code)
- `ModalOperator.__new__` (Method): Instantiates an operator node. (2 lines of code)
- `ModalOperator.eval` (Method): Base evaluation method for the operator. (2 lines of code)
- `Box` (Class): The Necessity operator (□). (16 lines of code)
- `Box.__new__` (Method): Instantiates a Box operator wrapping a formula. (4 lines of code)
- `Box.modality` (Method): Returns the modality type. (2 lines of code)
- `Box.is_well_typed` (Method): Checks if the internal formula is a valid proposition. (3 lines of code)
- `Diamond` (Class): The Possibility operator (◇). (15 lines of code)
- `Diamond.__new__` (Method): Instantiates a Diamond operator wrapping a formula. (4 lines of code)
- `Diamond.modality` (Method): Returns the modality type. (2 lines of code)
- `Diamond.is_well_typed` (Method): Checks if the internal formula is a valid proposition. (2 lines of code)
- `ProvabilityBox` (Class): A necessity operator hardcoded to Provability logic (GL). (4 lines of code)
- `ProvabilityBox.__new__` (Method): Instantiates a ProvabilityBox. (2 lines of code)
- `AlethicBox` (Class): A necessity operator hardcoded to Alethic logic (S5). (4 lines of code)
- `AlethicBox.__new__` (Method): Instantiates an AlethicBox. (2 lines of code)
- `EpistemicBox` (Class): A necessity operator hardcoded to Epistemic logic (K45). (4 lines of code)
- `EpistemicBox.__new__` (Method): Instantiates an EpistemicBox. (2 lines of code)
- `DeonticBox` (Class): A necessity operator hardcoded to Deontic logic (D). (4 lines of code)
- `DeonticBox.__new__` (Method): Instantiates a DeonticBox. (2 lines of code)
- `TemporalBox` (Class): A necessity operator hardcoded to Temporal logic (S4). (4 lines of code)
- `TemporalBox.__new__` (Method): Instantiates a TemporalBox. (2 lines of code)
- `ForAllPredicates` (Class): Second-order universal quantifier. (22 lines of code)
- `ForAllPredicates.__new__` (Method): Instantiates ForAllPredicates over a variable and formula. (6 lines of code)
- `ForAllPredicates.variable` (Method): Returns the bound variable. (2 lines of code)
- `ForAllPredicates.formula` (Method): Returns the quantified formula. (2 lines of code)
- `ForAllPredicates.is_well_typed` (Method): Ensures the quantification is typed correctly. (3 lines of code)
- `ExistsPredicates` (Class): Second-order existential quantifier. (21 lines of code)
- `ExistsPredicates.__new__` (Method): Instantiates ExistsPredicates over a variable and formula. (6 lines of code)
- `ExistsPredicates.variable` (Method): Returns the bound variable. (2 lines of code)
- `ExistsPredicates.formula` (Method): Returns the quantified formula. (2 lines of code)
- `ExistsPredicates.is_well_typed` (Method): Ensures the quantification is typed correctly. (2 lines of code)

### `sympy_modal/errors.py`
(74 lines of code) This module centralizes all custom exceptions raised within the package, ranging from structural `FrameViolationError`s and invalid deduction checks to parsing failures during formalisation.

**Declarations:**
- `SymPyModalError` (Class): Base exception for sympy_modal errors. (3 lines of code)
- `FrameViolationError` (Class): Raised when a formula violates the enforced Kripke frame. (6 lines of code)
- `FrameViolationError.__init__` (Method): Initializes with specific violation details. (2 lines of code)
- `NecessitationError` (Class): Raised when the Necessitation rule is applied invalidly. (7 lines of code)
- `NecessitationError.__init__` (Method): Initializes with necessitation failure details. (2 lines of code)
- `InvalidInferenceError` (Class): Raised on an invalid deduction step. (6 lines of code)
- `InvalidInferenceError.__init__` (Method): Initializes with deduction context. (2 lines of code)
- `FormalisationError` (Class): Raised when code strings cannot be mapped to modal terms. (6 lines of code)
- `FormalisationError.__init__` (Method): Initializes with parse failure information. (2 lines of code)
- `NotAnAxiomError` (Class): Raised when an expression claims to be an axiom but isn't. (7 lines of code)
- `NotAnAxiomError.__init__` (Method): Initializes with the invalid axiom claim. (2 lines of code)
- `AmbiguousModalityError` (Class): Raised when a frame cannot be heuristically inferred. (7 lines of code)
- `AmbiguousModalityError.__init__` (Method): Initializes with ambiguity details. (3 lines of code)
- `ProofFailure` (Class): Raised when the automated prover fails. (12 lines of code)
- `ProofFailure.__init__` (Method): Initializes with failure context. (4 lines of code)
- `ProofFailure.__repr__` (Method): Detailed representation of the proof failure. (2 lines of code)


## Proposals for Enhancing the Reference and Tutorial Examples

### Enhancing Example Output
Currently, the tutorial examples in this document rely on silent assertions (e.g., `assert s4.validates(...)`). To provide a better learning experience, the examples should be updated to produce meaningful console output. Instead of just passing silently, we propose replacing assertions with print statements that evaluate to boolean conditions or print formulas directly. For example, changing `assert s4.validates(Implies(Box(p), p))` to `print(f"S4 Reflexivity validation: {s4.validates(Implies(Box(p), p))}")` will output `S4 Reflexivity validation: True`. This provides users with direct, observable feedback when they copy, paste, and run the tutorial scripts locally.

### Expanding Example Explanations
To elevate the tutorial section from a mere list of code snippets to a comprehensive learning resource, each example should be expanded to include a structured, two-paragraph introduction and explanation.
- The **first paragraph** should introduce the specific logical concept or system component being demonstrated (e.g., explaining the intuitive meaning behind the Kripke accessibility relations required for S4).
- The **second paragraph** should explicitly walk through the code, explaining step-by-step how the components interact (e.g., explaining why `evaluate=False` is necessary, or how the `ProofContext` manages active hypotheses under the hood). This deeper context will bridge the gap between abstract modal logic theory and the concrete SymPy API implementation.
