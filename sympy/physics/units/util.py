# -*- coding: utf-8 -*-

"""
Several methods to simplify expressions involving unit objects.
"""

from __future__ import division

from sympy.utilities.exceptions import SymPyDeprecationWarning

from sympy import Add, Function, Mul, Pow, Rational, Tuple, sympify
from sympy.core.compatibility import reduce, Iterable
from sympy.physics.units.dimensions import Dimension, dimsys_default
from sympy.physics.units.quantities import Quantity
from sympy.physics.units.prefixes import Prefix
from sympy.utilities.iterables import sift


def dim_simplify(expr):
    """
    NOTE: this function could be deprecated in the future.

    Simplify expression by recursively evaluating the dimension arguments.

    This function proceeds to a very rough dimensional analysis. It tries to
    simplify expression with dimensions, and it deletes all what multiplies a
    dimension without being a dimension. This is necessary to avoid strange
    behavior when Add(L, L) be transformed into Mul(2, L).
    """
    SymPyDeprecationWarning(
        deprecated_since_version="1.2",
        feature="dimensional simplification function",
        issue=13336,
        useinstead="don't use",
    ).warn()
    _, expr = Quantity._collect_factor_and_dimension(expr)
    return expr


def _get_conversion_matrix_for_expr(expr, target_units):
    from sympy import Matrix

    expr_dim = Dimension(Quantity.get_dimensional_expr(expr))
    dim_dependencies = dimsys_default.get_dimensional_dependencies(expr_dim, mark_dimensionless=True)
    target_dims = [Dimension(Quantity.get_dimensional_expr(x)) for x in target_units]
    canon_dim_units = {i for x in target_dims for i in dimsys_default.get_dimensional_dependencies(x, mark_dimensionless=True)}
    canon_expr_units = {i for i in dim_dependencies}

    if not canon_expr_units.issubset(canon_dim_units):
        return None

    canon_dim_units = sorted(canon_dim_units)

    camat = Matrix([[dimsys_default.get_dimensional_dependencies(i, mark_dimensionless=True).get(j, 0)  for i in target_dims] for j in canon_dim_units])
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
    if not isinstance(target_units, (Iterable, Tuple)):
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


def quantity_simplify(expr):
    if expr.is_Atom:
        return expr
    if not expr.is_Mul:
        return expr.func(*map(quantity_simplify, expr.args))

    if expr.has(Prefix):
        coeff, args = expr.as_coeff_mul(Prefix)
        args = list(args)
        for arg in args:
            if isinstance(arg, Pow):
                coeff = coeff * (arg.base.scale_factor ** arg.exp)
            else:
                coeff = coeff * arg.scale_factor
        expr = coeff

    coeff, args = expr.as_coeff_mul(Quantity)
    args_pow = [arg.as_base_exp() for arg in args]
    quantity_pow, other_pow = sift(args_pow, lambda x: isinstance(x[0], Quantity), binary=True)
    quantity_pow_by_dim = sift(quantity_pow, lambda x: x[0].dimension)
    # Just pick the first quantity:
    ref_quantities = [i[0][0] for i in quantity_pow_by_dim.values()]
    new_quantities = [
        Mul.fromiter(
            (quantity*i.scale_factor/quantity.scale_factor)**p for i, p in v)
            if len(v) > 1 else v[0][0]**v[0][1]
        for quantity, (k, v) in zip(ref_quantities, quantity_pow_by_dim.items())]
    return coeff*Mul.fromiter(other_pow)*Mul.fromiter(new_quantities)
