"""Numerical Methods for Holonomic Functions"""

from __future__ import print_function, division

from sympy.core.sympify import sympify
from mpmath import mp
from sympy.holonomic.holonomic import DMFsubs


def euler(func, points):
    ann = func.annihilator
    a = ann.order
    R = ann.parent.base
    K = R.get_field()

    dmf = []
    for j in ann.listofpoly:
        dmf.append(K.new(j.rep))

    red = [-dmf[i] / dmf[a] for i in range(a)]

    y0 = func.y0
    if len(y0) < a:
        raise TypeError("Not Enough Initial Conditions")
    x0 = func.x0
    sol = _eular(red, x0, points[0], y0, a)

    for i, j in enumerate(points[1:]):
        sol = _eular(red, points[i], j, sol, a)

    return sol


def _eular(red, x0, x1, y0, a):
    A = sympify(x0)._to_mpmath(mp.prec)
    B = sympify(x1)._to_mpmath(mp.prec)
    y_0 = [sympify(i)._to_mpmath(mp.prec) for i in y0]
    h = B - A
    f_0 = y_0[1:]
    f_0_n = 0

    for i in range(a):
        f_0_n += sympify(DMFsubs(red[i], A))._to_mpmath(mp.prec) * y_0[i]
    f_0.append(f_0_n)

    sol = []
    for i in range(a):
        sol.append(y_0[i] + h * f_0[i])

    return sol
