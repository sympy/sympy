# -*- coding: utf-8 -*-

"""
Several methods to simplify expressions involving unit objects.
"""

from __future__ import division

from sympy.physics.units.quantities import Quantity
from sympy import Add, Mul, Pow, Function, Rational
from sympy.core.compatibility import reduce
from sympy.physics.units.dimensions import Dimension


def dim_simplify(expr):
    """
    NOTE: this function could be deprecated in the future.

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
        if (all(isinstance(arg, Dimension) for arg in expr.args) or
            all(arg.is_dimensionless for arg in expr.args if isinstance(arg, Dimension))):
            return reduce(lambda x, y: x.add(y), expr.args)
        else:
            raise ValueError("Dimensions cannot be added: %s" % expr)
    elif isinstance(expr, Mul):
        return Dimension(Mul(*[dim_simplify(i).name for i in expr.args if isinstance(i, Dimension)]))

    raise ValueError("Cannot be simplifed: %s", expr)


def convert_to(expr, quantity):
    """
    Convert `expr` to the same expression with all of its units and quantities
    represented as factors of `quantity`, whenever the dimension is compatible.

    Examples
    ========

    >>> from sympy.physics.units import speed_of_light, meter, gram, \
        second, day, mile, newton, kilogram, inch, centimeter, atomic_mass_constant
    >>> from sympy.physics.units.definitions import kilometer
    >>> from sympy.physics.units import convert_to
    >>> convert_to(mile, kilometer)
    25146*kilometer/15625
    >>> convert_to(mile, kilometer).n()
    1.609344*kilometer
    >>> convert_to(speed_of_light, meter/second)
    299792458*meter/second
    >>> convert_to(day, second)
    86400*second
    >>> 3*newton
    3*newton
    >>> convert_to(3*newton, kilogram*meter/second**2)
    3*kilogram*meter/second**2
    >>> convert_to(atomic_mass_constant, gram)
    1.66053904e-24*gram
    """

    def get_total_scale_factor(expr):
        if isinstance(expr, Mul):
            return reduce(lambda x, y: x*y, [get_total_scale_factor(i) for i in expr.args])
        elif isinstance(expr, Pow):
            return get_total_scale_factor(expr.base)**expr.exp
        elif isinstance(expr, Quantity):
            return expr.scale_factor
        return 1

    def get_units(expr):
        if isinstance(expr, Mul):
            return reduce(lambda x, y: x*y, [get_units(i) for i in expr.args])
        elif isinstance(expr, Pow):
            return get_units(expr.base)**expr.exp
        elif isinstance(expr, Quantity):
            return expr
        return 1

    if isinstance(quantity, Quantity):
        backup_quantity = None
    else:
        backup_quantity = quantity
        quantity = Quantity("_temp", Dimension(Quantity.get_dimensional_expr(quantity)), get_total_scale_factor(quantity))

    def _convert_to(expr, quantity):
        if isinstance(expr, Add):
            return Add(*[_convert_to(i, quantity) for i in expr.args])
        elif isinstance(expr, Mul):
            new_args = [_convert_to(i, quantity) for i in expr.args]
            edim = Dimension(Quantity.get_dimensional_expr(expr))
            if edim == quantity.dimension:
                scale_factor_old = get_total_scale_factor(expr)
                return expr / get_units(expr) * scale_factor_old / quantity.scale_factor * quantity
            return Mul(*new_args)
        elif isinstance(expr, Pow):
            base = _convert_to(expr.base, quantity)
            edim = Dimension(Quantity.get_dimensional_expr(base))**expr.exp
            if edim == quantity.dimension:
                scale_factor_old = get_total_scale_factor(expr)
                return expr / get_units(expr) * scale_factor_old / quantity.scale_factor * quantity
            return base**expr.exp
        elif isinstance(expr, Quantity):
            edim = Dimension(Quantity.get_dimensional_expr(expr))
            edep1 = edim.get_dimensional_dependencies()
            edep2 = quantity.dimension.get_dimensional_dependencies()
            if edim == quantity.dimension:
                return expr.scale_factor / quantity.scale_factor * quantity
            if set(edep1.keys()) == set(edep2.keys()):
                fracs = [Rational(v1, v2) for v1, v2 in zip(edep1.values(), edep2.values())]
                powers = list(set(fracs))
                if len(powers) == 1:
                    return expr.scale_factor / quantity.scale_factor**powers[0] * quantity**powers[0]
            else:
                return expr
        return expr

    res = _convert_to(expr, quantity)
    if backup_quantity:
        res = res.subs(quantity, backup_quantity)
    return res
