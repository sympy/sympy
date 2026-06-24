from sympy import Symbol, Implies
from sympy_modal import KripkeFrame, Box

s4 = KripkeFrame.S4()
p = Symbol('p')

# Reflexivity: □P → P
print(f"S4 Reflexivity validation: {s4.validates(Implies(Box(p), p))}")

# Transitivity: □P → □□P
print(f"S4 Transitivity validation: {s4.validates(Implies(Box(p), Box(Box(p))))}")
