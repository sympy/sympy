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
