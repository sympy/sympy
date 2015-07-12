from __future__ import division, print_function

from sympy.core import sympify, Symbol

x = Symbol('x')


def timeit_sympify_1():
    sympify(1)


def timeit_sympify_x():
    sympify(x)
