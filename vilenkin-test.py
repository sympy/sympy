#!/usr/bin/env python

# Equations from Vilenkin book

from sympy import *
from sympy.utilities.solution import *

x = Symbol('x')

eqs = [
    x**sqrt(x) - x**(x/2),
    4 - log(x, 10) - 3 * sqrt(log(x, 10)),
    log(x, Rational(1, 2)) + log(x, 3) - 1,
    x**log(x, 10) - x**100,
    sin(3*x) * cos(2 * x) * tan(7 * x),
    cos(x**2) + cos(5*x**2),
    sqrt(3) * sin(x) + cos(x) - sqrt(2),
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