# -*- coding: utf-8 -*-

"""
Dimensional analysis and unit systems.

This module defines dimension/unit systems and physical quantities. It is
based on a group-theoretical construction where dimensions are represented as
vectors (coefficient being the exponents), and units are defined as a dimension
to which we added a scale.

Quantities are built from a factor and a unit, and are the basic objects that
one will use when doing computations.

All objects except systems and prefixes can be used in sympy expressions.

Details about the implementation can be found in the documentation, and we
will not repeat all the explanations we gave there concerning our approach.
Ideas about future developments can be found on the `Github wiki
<https://github.com/sympy/sympy/wiki/Unit-systems>`_, and you should consult
this page if you are willing to help.
"""

from __future__ import division

from sympy import Add, Mul, Pow
from sympy.core.compatibility import reduce
from sympy.physics.unitsystems.dimensions import Dimension
from sympy.physics.unitsystems.units import Unit, Constant
from sympy.physics.unitsystems.quantities import Quantity


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


def qsimplify(expr):
    """
    Simplify expression by recursively evaluating the quantity arguments.

    If units are encountered, as it can be when using Constant, they are
    converted to quantity.
    """

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
        return args[0]**args[1]
    elif isinstance(expr, Add):
        if q_args != []:
            quantities = reduce(lambda x, y: x+y, q_args)
        else:
            quantities = []
        return reduce(lambda x, y: x+y, o_args, quantities)
    elif isinstance(expr, Mul):
        if q_args != []:
            quantities = reduce(lambda x, y: x*y, q_args)
        else:
            quantities = []
        return reduce(lambda x, y: x*y, o_args, quantities)
    else:
        return expr
