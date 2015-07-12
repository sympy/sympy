from __future__ import division, print_function

from sympy import symbols
from sympy.functions.special.spherical_harmonics import Ynm

x, y = symbols('x,y')


def timeit_Ynm_xy():
    Ynm(1, 1, x, y)
