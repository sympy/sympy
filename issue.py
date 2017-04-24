from __future__ import division, print_function

from sympy import *

x = Symbol('x', real = True)

expr = x**2 + x*sin(x) + cos(x)
l = limit(expr, x, -oo)

