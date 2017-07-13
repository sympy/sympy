# -*- coding: utf-8 -*-

"""
Several methods to simplify expressions involving unit objects.
"""

from __future__ import division

import collections

from sympy.physics.units.quantities import Quantity
from sympy import Add, Mul, Pow, Function, Rational, Tuple, sympify
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


def _get_conversion_matrix_for_expr(expr, target_units):
    from sympy import Matrix

    expr_dim = Dimension(Quantity.get_dimensional_expr(expr))
    dim_dependencies = expr_dim.get_dimensional_dependencies(mark_dimensionless=True)
    target_dims = [Dimension(Quantity.get_dimensional_expr(x)) for x in target_units]
    canon_dim_units = {i for x in target_dims for i in x.get_dimensional_dependencies(mark_dimensionless=True)}
    canon_expr_units = {i for i in dim_dependencies}

    if not canon_expr_units.issubset(canon_dim_units):
        return None

    canon_dim_units = sorted(canon_dim_units)

    camat = Matrix([[i.get_dimensional_dependencies(mark_dimensionless=True).get(j, 0)  for i in target_dims] for j in canon_dim_units])
    exprmat = Matrix([dim_dependencies.get(k, 0) for k in canon_dim_units])

    res_exponents = camat.solve_least_squares(exprmat, method=None)
    return res_exponents


def convert_to(expr, target_units):
    """
    Convert ``expr`` to the same expression with all of its units and quantities
    represented as factors of ``target_units``, whenever the dimension is compatible.

    ``target_units`` may be a single unit/quantity, or a collection of
    units/quantities.

    Examples
    ========

    >>> from sympy.physics.units import speed_of_light, meter, gram, second, day
    >>> from sympy.physics.units import mile, newton, kilogram, atomic_mass_constant
    >>> from sympy.physics.units import kilometer, centimeter
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

    Conversion to multiple units:

    >>> convert_to(speed_of_light, [meter, second])
    299792458*meter/second
    >>> convert_to(3*newton, [centimeter, gram, second])
    300000*centimeter*gram/second**2

    Conversion to Planck units:

    >>> from sympy.physics.units import gravitational_constant, hbar
    >>> convert_to(atomic_mass_constant, [gravitational_constant, speed_of_light, hbar]).n()
    7.62950196312651e-20*gravitational_constant**(-0.5)*hbar**0.5*speed_of_light**0.5

    """
    if not isinstance(target_units, (collections.Iterable, Tuple)):
        target_units = [target_units]

    if isinstance(expr, Add):
        return Add.fromiter(convert_to(i, target_units) for i in expr.args)

    expr = sympify(expr)

    if not isinstance(expr, Quantity) and expr.has(Quantity):
        expr = expr.replace(lambda x: isinstance(x, Quantity), lambda x: x.convert_to(target_units))

    def get_total_scale_factor(expr):
        if isinstance(expr, Mul):
            return reduce(lambda x, y: x * y, [get_total_scale_factor(i) for i in expr.args])
        elif isinstance(expr, Pow):
            return get_total_scale_factor(expr.base) ** expr.exp
        elif isinstance(expr, Quantity):
            return expr.scale_factor
        return expr

    depmat = _get_conversion_matrix_for_expr(expr, target_units)
    if depmat is None:
        return expr

    expr_scale_factor = get_total_scale_factor(expr)
    return expr_scale_factor * Mul.fromiter((1/get_total_scale_factor(u) * u) ** p for u, p in zip(target_units, depmat))
