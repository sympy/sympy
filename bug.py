#!/usr/bin/env python

from __future__ import division, print_function

from sympy import *
from sympy.physics.vector import vlatex

x = symbols('x')
J = symbols('J')

f = Function('f')
g = Function('g')
h = Function('h')


# print(J*f(x).diff(x))

# print(J*f(x).diff(x).subs(f(x), g(x)-h(x)))

tmp = J*f(x).diff(x).subs(f(x), g(x)-h(x))

# print(tmp)



# print(latex(tmp))

# print(vlatex(tmp))

latex(tmp)
vlatex(tmp)
