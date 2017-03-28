#!/usr/bin/env python

from __future__ import division, print_function

from sympy import *
from sympy.solvers.solveset import invert_real

x = Symbol('x')
expr = log(x)/x 

dom = Interval(0, oo)

print(invert_real(expr, 0, x, dom))

