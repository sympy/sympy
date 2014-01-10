#!/usr/bin/env python

# Equations of the form trig(x) = m

from sympy import *
from sympy.utilities.solution import *

x = Symbol('x')

eqs = [
    sin(x),
    sin(x) - 1,
    sin(x) + 1,
    3 * sin(x) - 1,
    sin(x) + 5 * sin(x),
    sin(x) - Rational(1, 2),
    sin(x) - sqrt(2) / 2,
    sin(x) - Rational(1, 3),
    cos(x),
    cos(x) + 1,
    cos(x) - 1,
    tan(x) + -sqrt(3),
    cot(x) - 1,
    sin(2*x) - Rational(1, 2),
    sin(2*x / 3),
    sin(2*x / 5 - 1),
    sin(3*x + pi / 4) + 1,
    cos(x / 3) + sqrt(2) / 2,
    tan(pi/4 - x/2) + 1,
    cot(pi/6 - x) - sqrt(3) / 3
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