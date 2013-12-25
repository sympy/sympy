#!/usr/bin/env python

# equations contain abs-sign

from sympy import *
from sympy.utilities.solution import *

x = Symbol('x')

eqs = [
    abs(x) - 5,
    abs(3 * x + 4) - 7,
    abs(2 - 5 * x) + 3,
    abs(2 * x - 5) - abs(4 - x) + 18,
    abs(abs(x) - 3) - 15,
    abs(x**2 - 1*x) - 2
]

for eq in eqs:
    print '===================================================='
    print '=== Equation: ' + latex(eq)

    reset_solution()
    res = solve(eq, x)
    R = last_solution()
    for r in R:
        print r
    print '=== Answer:'
    for r in res:
        print latex(r)