#!/usr/bin/env python

from __future__ import division, print_function

from sympy import *

x = Symbol('x')


expr = log(x)/x <= 0

solveset(expr,x,S.Reals)

print(solveset(expr,x,S.Reals))
