from __future__ import division, print_function

from sympy import Symbol, sin
from sympy.integrals.trigonometry import trigintegrate

x = Symbol('x')


def timeit_trigintegrate_sin3x():
    trigintegrate(sin(x)**3, x)


def timeit_trigintegrate_x2():
    trigintegrate(x**2, x)  # -> None
