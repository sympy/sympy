# -*- coding: utf-8 -*-

from __future__ import division

from sympy import Add, Mul, Pow
from sympy.physics.unitsystems.dimensions import Dimension


def dim_simplify(expr):
    """
    Simplify expression by recursively evaluating the dimension arguments.

    This function proceeds to a very rough dimensional analysis. It tries to
    simplify expression with dimensions, and it deletes all what multiplies a
    dimension without being a dimension. This is necessary to avoid strange
    behavior when Add(L, L) be transformed into Mul(2, L).
    """

    args = []
    for arg in expr.args:
        if isinstance(arg, (Mul, Pow, Add)):
            arg = dim_simplify(arg)
        args.append(arg)

    if isinstance(expr, Pow):
        return args[0]**args[1]
    elif isinstance(expr, Add):
        return reduce(lambda x, y: x+y, args)
    elif isinstance(expr, Mul):
        args = [arg for arg in args if isinstance(arg, Dimension)]
        return reduce(lambda x, y: x*y, args)
    else:
        return expr
