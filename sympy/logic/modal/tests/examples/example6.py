from sympy import Symbol, Or, Not, Implies
from sympy.logic.modal import ProofContext, KripkeFrame

print("Example 6: Classical Logic Injection and SMT Integration")

# Standard context rejects LEM without the flag
ctx_intuitionistic = ProofContext(KripkeFrame.K())
p = Symbol('p')
lem = Or(p, Not(p))
proof_fail = ctx_intuitionistic.prove(lem)
print(f"Intuitionistic proof valid: {isinstance(proof_fail, type(ctx_intuitionistic.assume(p)))}") # Will be false (ProofFailure)

# With classical logic enabled
ctx_classical = ProofContext(KripkeFrame.K(), allow_classical=True)
proof_success = ctx_classical.prove(lem)
print(f"Classical proof valid (LEM): {proof_success.is_valid}")
print(f"Classical derivation source: {proof_success.derivation}")

# SMT solver solving a complex propositional tautology instantly
q = Symbol('q')
complex_tautology = Or(Implies(p, q), Implies(q, p))
smt_proof = ctx_classical.prove(complex_tautology)
print(f"SMT proof valid: {smt_proof.is_valid}")
print(f"SMT derivation source: {smt_proof.derivation[0]}")
