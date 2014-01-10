#!/usr/bin/env python

# Trigonometric equations reducible to algebraic equations.

from sympy import *
from sympy.utilities.solution import *

x = Symbol('x')

eqs = [
    asin(x) - 1,
    asin(2 * x + 1) - 0,
    asin(x) + 4,
    acos(x - 3) - pi / 2,
    acos(x) + 1,
    atan(x) - pi,
    acot(2 *x - 4) - pi / 3,
    5 * cos(x) ** 2 - 5 * cos(x) + 1,
    8 * cos(x) ** 2 + 6 * sin(x) - 3,
    3 * tan(x) ** 3 + tan(x),
    sin(3 * x) * cos(4 * x)
]

for eq in eqs:
    print '===================================================='
    print '=== Equation: ' + latex(eq) + ' = 0'

    reset_solution()
    res = solve(eq, x)
    R = last_solution()
    for r in R: 
        print r
    print '=== Answer:'
    for r in res:
        print latex(r)