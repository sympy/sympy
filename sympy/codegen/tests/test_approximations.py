# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from sympy import symbols, exp
from sympy.codegen.rewriting import optimize
from sympy.codegen.approximations import SumApprox


def test_SumApprox_trivial():
    x = symbols('x')
    expr1 = 1 + x
    sum_approx = SumApprox(bounds={x: (-1e-20, 1e-20)}, reltol=1e-16)
    apx1 = optimize(expr1, [sum_approx])
    assert apx1 - 1 == 0


def test_SumApprox_monotone_terms():
    x, y, z = symbols('x y z')
    expr1 = exp(z)*(x**2 + y**2 + 1)
    bnds1 = {x: (0, 1e-3), y: (100, 1000)}
    sum_approx_m2 = SumApprox(bounds=bnds1, reltol=1e-2)
    sum_approx_m5 = SumApprox(bounds=bnds1, reltol=1e-5)
    sum_approx_m11 = SumApprox(bounds=bnds1, reltol=1e-11)
    assert (optimize(expr1, [sum_approx_m2])/exp(z) - (y**2)).simplify() == 0
    assert (optimize(expr1, [sum_approx_m5])/exp(z) - (y**2 + 1)).simplify() == 0
    assert (optimize(expr1, [sum_approx_m11])/exp(z) - (y**2 + 1 + x**2)).simplify() == 0
