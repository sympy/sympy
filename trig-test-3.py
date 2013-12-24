#!/usr/bin/env python

# Trigonometric equations reducible to algebraic equations.

from sympy import *
from sympy.utilities.solution import *

x = Symbol('x')

eqs = [
    5 * cos(x) ** 2 - 5 * cos(x) + 1,
    8 * cos(x) ** 2 + 6 * sin(x) - 3,
    3 * tan(x) ** 3 + tan(x),
    sin(3 * x) * cos(4 * x)
]

for eq in eqs:
    print '===================================================='
    print '=== Equation: ' + str(eq)

    reset_solution()
    res = solve(eq, x)
    R = last_solution()
    for r in R: 
        print r
    print '=== Answer: ' + str(res)
