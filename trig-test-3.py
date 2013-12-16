#!/usr/bin/env python
# vim: set smartindent expandtab shiftwidth=4 softtabstop=4:

# Equations of the form trig(x) = m

from sympy import *
from sympy.utilities.solution import *
from sympy.simplify.fu import *

x = Symbol('x')
y = Symbol('y')
print srepr(Add(x, x, evaluate=False))

z = sin(x)**2 + cos(x)
print TR6(z)

eqs = [
    5 * cos(x) ** 2 - 5 * cos(x) + 1,
    8 * cos(x) ** 2 + 6 * sin(x) - 3,
    3 * tan(x) ** 3 + tan(x),
    sin(3*x) * cos(4 * x)
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
