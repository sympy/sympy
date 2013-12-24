#!/usr/bin/env python

# Exponential equations

from sympy import *
from sympy.utilities.solution import *

x = Symbol('x')

eqs = [
    2**(pi*x + E) - 4,
    2**x - 8,
    5**(x + 2) - 125,
    2**(2 * x) - 8**(x + 1),
    3**(2 * x + 4) - 11 * 9**x - 210,
    4**x - 3 * 2**x + 2,
    2**(5 * x - 1) * 3**(3 * x - 1) * 5**(2 * x - 1) - 720**x
]

for eq in eqs:
    print '===================================================='
    print '=== Equation: ' + str(eq)

    reset_solution()
    res = solve(eq, x)
    R = last_solution()
    for r in R:
        print r
    print '=== Answer: ' + repr(res)
