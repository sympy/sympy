# -*- coding: utf-8 -*-

"""
Several methods to simplify expressions involving unit objects.
"""

from __future__ import division

from sympy import Add, Mul, Pow, Function
from sympy.core.compatibility import reduce
from sympy.physics.unitsystems import Dimension, Unit, Quantity


def dim_simplify(expr):
    """
    Simplify expression by recursively evaluating the dimension arguments.

    This function proceeds to a very rough dimensional analysis. It tries to
    simplify expression with dimensions, and it deletes all what multiplies a
    dimension without being a dimension. This is necessary to avoid strange
    behavior when Add(L, L) be transformed into Mul(2, L).
    """

    if isinstance(expr, Dimension):
        return expr

    if isinstance(expr, Pow):
        return dim_simplify(expr.base)**dim_simplify(expr.exp)
    elif isinstance(expr, Function):
        return dim_simplify(expr.args[0])
    elif isinstance(expr, Add):
        if (all(isinstance(arg, Dimension) for arg in args) or
            all(arg.is_dimensionless for arg in args if isinstance(arg, Dimension))):
            return reduce(lambda x, y: x.add(y), args)
        else:
            raise ValueError("Dimensions cannot be added: %s" % expr)
    elif isinstance(expr, Mul):
        return Dimension(Mul(*[dim_simplify(i).name for i in expr.args if isinstance(i, Dimension)]))

    raise ValueError("Cannot be simplifed: %s", expr)
