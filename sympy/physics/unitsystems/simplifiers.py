# -*- coding: utf-8 -*-

"""
Several methods to simplify expressions involving unit objects.
"""

from __future__ import division

from sympy import Add, Mul, Pow, Function
from sympy.core.compatibility import reduce
from sympy.physics.unitsystems.quantities import Quantity
from sympy.physics.unitsystems.dimensions import Dimension


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

    >>> from sympy.physics.unitsystems import speed_of_light, meter, hour, minute, second, day, mile, newton, kilogram
    >>> from sympy.physics.unitsystems.definitions import kilometer
    >>> from sympy.physics.unitsystems import convert_to

    >>> from sympy import solve, symbols, pi, Eq
    >>> from sympy.physics.unitsystems import Unit, Quantity, length, mass
    >>> from sympy.physics.unitsystems.systems import mks
    >>> from sympy.physics.unitsystems import day, gravitational_constant as G
    >>> T = symbols("T")
    >>> a = Quantity("venus_a", length, 108208000e3)
    >>> M = Quantity("solar_mass", mass, 1.9891e30)
    >>> Eq(T**2 / a**3, 4*pi**2 / G / M)
    Eq(T**2/venus_a**3, 4*pi**2/(gravitational_constant*solar_mass))
    >>> q = solve(T**2 / a**3 - 4*pi**2 / G / M, T)[1]
    >>> q
    6.28318530717959*sqrt(venus_a**3/(gravitational_constant*solar_mass))
    >>> convert_to(q, day)
    1.59123003109442e-10*sqrt(1993480633430344427209462)*day
    >>> convert_to(q, day).n()
    224.666840070582*day

    # TODO: move these to the tests:

    >>> convert_to(mile, kilometer)
    1.609344*kilometer
    >>> convert_to(meter/second, speed_of_light)
    speed_of_light/299792458
    >>> convert_to(299792458*meter/second, speed_of_light)
    speed_of_light
    >>> convert_to(2*299792458*meter/second, speed_of_light)
    2*speed_of_light
    >>> convert_to(speed_of_light, meter/second)
    299792458*meter/second
    >>> convert_to(2*speed_of_light, meter/second)
    599584916*meter/second
    >>> convert_to(day, second)
    86400*second
    >>> convert_to(2*hour, minute)
    120*minute
    >>> convert_to(mile, meter)
    1609.344*meter
    >>> convert_to(mile/hour, kilometer/hour)
    1.609344*kilometer/hour
    >>> convert_to(3*newton, meter/second)
    3*newton
    >>> convert_to(3*newton, kilogram*meter/second**2)
    3*meter*kilogram/second**2
    """
    def get_dimensional_expr(expr):
        if isinstance(expr, Mul):
            return Mul(*[get_dimensional_expr(i) for i in expr.args])
        elif isinstance(expr, Pow):
            return get_dimensional_expr(expr.base)**expr.exp
        elif isinstance(expr, Add):
            #return get_dimensional_expr()
            raise NotImplementedError
        elif isinstance(expr, Quantity):
            return expr.dimension.name
        return 1

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
        quantity = Quantity("_temp", Dimension(get_dimensional_expr(quantity)), get_total_scale_factor(quantity))

    def _convert_to(expr, quantity):
        if isinstance(expr, Add):
            return Add(*[_convert_to(i, quantity) for i in expr.args])
        elif isinstance(expr, Mul):
            new_args = [_convert_to(i, quantity) for i in expr.args]
            edim = Dimension(get_dimensional_expr(expr))
            if edim == quantity.dimension:
                scale_factor_old = get_total_scale_factor(expr)
                return expr / get_units(expr) * scale_factor_old / quantity.scale_factor * quantity
            return Mul(*new_args)
        elif isinstance(expr, Pow):
            base = _convert_to(expr.base, quantity)
            edim = Dimension(get_dimensional_expr(base))**expr.exp
            if edim == quantity.dimension:
                scale_factor_old = get_total_scale_factor(expr)
                return expr / get_units(expr) * scale_factor_old / quantity.scale_factor * quantity
            return base**expr.exp
        elif isinstance(expr, Quantity):
            edim = Dimension(get_dimensional_expr(expr))
            if edim == quantity.dimension:
                return expr.scale_factor / quantity.scale_factor * quantity
            else:
                return expr
        return expr

    res = _convert_to(expr, quantity)
    if backup_quantity:
        res = res.subs(quantity, backup_quantity)
    return res
