from __future__ import print_function, division

from sympy.core import symbols, S, C

x, y = symbols('x,y')


def timeit_Symbol_meth_lookup():
    x.diff  # no call, just method lookup


def timeit_S_lookup():
    S.Exp1


def timeit_C_lookup():
    C.Add


def timeit_Symbol_eq_xy():
    x == y
