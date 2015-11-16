# -*- coding: utf-8 -*-

"""
Several methods to simplify expressions involving unit objects.
"""

from __future__ import division

from sympy import Add, Mul, Pow
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

    args = []
    for arg in expr.args:
        if isinstance(arg, (Mul, Pow, Add)):
            arg = dim_simplify(arg)
        args.append(arg)

    if all([arg.is_number or (isinstance(arg, Dimension) and arg.is_dimensionless) for arg in args]):
        return Dimension({})

    if isinstance(expr, Pow):
        if isinstance(args[0], Dimension):
            return args[0].pow(args[1])
        else:
            raise ValueError("Basis of Pow is not a Dimension: %s" % args[0])
    elif isinstance(expr, Add):
        if (all(isinstance(arg, Dimension) for arg in args) or
            all(arg.is_dimensionless() for arg in args if isinstance(arg, Dimension))):
            return reduce(lambda x, y: x.add(y), args)
        else:
            raise ValueError("Dimensions cannot be added: %s" % expr)
    elif isinstance(expr, Mul):
        args = [arg for arg in args if isinstance(arg, Dimension)]
        return reduce(lambda x, y: x.mul(y), args)

    raise ValueError("Cannot be simplifed: %s", expr)

def qsimplify(expr):
    """
    Simplify expression by recursively evaluating the quantity arguments.

    If units are encountered, as it can be when using Constant, they are
    converted to quantity.
    """

    def redmul(x, y):
        """
        Function used to combine args in multiplications.

        This is necessary because the previous computation was not commutative,
        and Mul(3, u) was not simplified; but Mul(u, 3) was.
        """
        if isinstance(x, Quantity):
            return x.mul(y)
        elif isinstance(y, Quantity):
            return y.mul(x)
        else:
            return x*y

    args = []
    for arg in expr.args:
        arg = arg.evalf()
        if isinstance(arg, (Mul, Pow, Add)):
            arg = qsimplify(arg)
        args.append(arg)

    q_args, o_args = [], []

    for arg in args:
        if isinstance(arg, Quantity):
            q_args.append(arg)
        elif isinstance(arg, Unit):
            # replace unit by a quantity to make the simplification
            q_args.append(arg.as_quantity)
        else:
            o_args.append(arg)

    if isinstance(expr, Pow):
        return args[0].pow(args[1])
    elif isinstance(expr, Add):
        if q_args != []:
            quantities = reduce(lambda x, y: x.add(y), q_args)
        else:
            quantities = []
        return reduce(lambda x, y: x+y, o_args, quantities)
    elif isinstance(expr, Mul):
        if q_args != []:
            quantities = reduce(lambda x, y: x.mul(y), q_args)
        else:
            quantities = []
        return reduce(redmul, o_args, quantities)
    else:
        return expr
