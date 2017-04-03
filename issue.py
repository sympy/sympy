from __future__ import division, print_function

from sympy import *
from sympy.solvers.solveset import *

x = Symbol('x')

expr = sqrt(x - 2) + 2

ar = solveset_real(expr, x)

a = solveset(expr, x)


print(ar)

print(a)

# q = Poly(y**4 - 3*y**3 + (-3*sqrt(2) + 1)*y**2 + 2*sqrt(2)*y + 2, y, domain='QQ<sqrt(2)>')
