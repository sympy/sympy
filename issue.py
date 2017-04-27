from __future__ import division, print_function

from sympy import *

x = Symbol('x', real = True)

expr = x**2 + x*sin(x) + cos(x)
l = limit(expr, x, -oo)

expr1 = x**2 + x*sin(x)
expr2 = x * (x + sin(x))

l1 = limit(expr1, x, -oo)
l2 = limit(expr2, x, -oo)


l1 = limit(expr1, x, -oo)

expr3 = sin(1/x)/x + x**(-2)
expr4 = (sin(1/x) + 1/x)/x

l3 = limit(expr3, x, 0)
l4 = limit(expr4, x, 0)
