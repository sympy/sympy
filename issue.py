from __future__ import division, print_function

from sympy import *

x = Symbol('x', real = True)

expr = (3**x + 2* x**10) / (x**10 + E**x)
l = limit(expr, x, -oo)
