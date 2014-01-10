#!/usr/bin/env python

# Logarithmic equations

from sympy import *
from sympy.utilities.solution import *

x = Symbol('x')

eqs = [
    log(x, 2) - 10,
    log(50 * x - 1, 7) - 5,
    log(x, Rational(1, 3)) - 2,
    log(2 * x - 1, Rational(1, 3)) - 2,
    log(8, x - 1) - 1,
    ln(E**2 + 2*x - 3) - 2,
    log(x, 3) - log(9, 3),
    log(x**2 - 3, 3) - log(2*x, 3),
    2 * log((x - 1)**2, 7) + log((2*x + 9) / (7 * x + 9), sqrt(7)),
    log(x + 1)**2 + 10 - 11 * log(x + 1),
    log(x**2 + 9*x, 10) + log((x + 9) / x, 10),
    log(6 * sin(x) + 4, 3) * log(6 * sin(x) + 4, 5) - log(6 * sin(x) + 4, 3) - log(6 * sin(x) + 4, 5),
    log(x**2 + 5*x - 6, 2) - log(4 * x, 2),
    log((x**3 - 5 * x**2) / (x - 5), 5) - 2,
    log(2*x)**2 + 3 * log(2 * x) + 2
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