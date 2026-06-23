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
