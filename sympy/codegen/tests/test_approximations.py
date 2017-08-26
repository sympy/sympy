# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from sympy import symbols
from sympy.codegen.rewriting import optimize
from sympy.codegen.approximations import SumApprox


def test_sum_approx():
    x, y = symbols('x y')
    expr1 = 1 + x
    sum_approx = SumApprox(bounds={x: (-1e-20, 1e-20)}, reltol=1e-16)
    apx1 = optimize(expr1, [sum_approx])
    assert apx1 - 1 == 0
