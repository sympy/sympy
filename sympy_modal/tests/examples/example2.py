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
