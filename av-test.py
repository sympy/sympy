#!/usr/bin/env python

# admissible values test

from sympy import *
from sympy.solvers.solvers import sub_trig_solution, _k
from sympy.utilities.solution import *

x = Symbol('x')

eqs = [
    sin(3*x) * cot(4*x),
    sin(6 * x) / sin(4 * x),
    sin(2*x) * sin(4*x) * sin(6*x) / sin(x)
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