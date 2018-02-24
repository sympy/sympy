# -*- coding: utf-8 -*-

"""
Physical quantities.
"""

from __future__ import division

from sympy import (Abs, Add, AtomicExpr, Basic, Derivative, Function, Mul,
    Pow, S, Symbol, sympify)
from sympy.core.compatibility import string_types
from sympy.physics.units import Dimension, dimensions
from sympy.physics.units.dimensions import dimsys_default, DimensionSystem
from sympy.physics.units.prefixes import Prefix
from sympy.utilities.exceptions import SymPyDeprecationWarning


class Quantity(AtomicExpr):
    """
    Physical quantity: can be a unit of measure, a constant or a generic quantity.
    """

    is_commutative = True
    is_real = True
    is_number = False
    is_nonzero = True
    _diff_wrt = True

    def __new__(cls, name, abbrev=None, dimension=None, scale_factor=None, **assumptions):

        if not isinstance(name, Symbol):
            name = Symbol(name)

        if dimension is not None:
            SymPyDeprecationWarning(
                deprecated_since_version="1.3",
                issue=14319,
                feature="Quantity arguments",
                useinstead="SI_quantity_dimension_map",
            ).warn()

        if scale_factor is not None:
            SymPyDeprecationWarning(
                deprecated_since_version="1.3",
                issue=14319,
                feature="Quantity arguments",
                useinstead="SI_quantity_scale_factors",
            ).warn()

        #if not isinstance(dim_sys, DimensionSystem):
            #raise TypeError("%s is not a DimensionSystem" % dim_sys)

        if abbrev is None:
            abbrev = name
        elif isinstance(abbrev, string_types):
            abbrev = Symbol(abbrev)

        obj = AtomicExpr.__new__(cls, name, abbrev)
        obj._name = name
        obj._abbrev = abbrev

        if scale_factor is not None:
            # TODO: remove after deprecation:
            from sympy.physics.units.definitions import SI_quantity_scale_factors
            scale_factor = process_scale_factor(scale_factor)
            SI_quantity_scale_factors[obj] = scale_factor
        return obj

    @property
    def name(self):
        return self._name

    @property
    def dimension(self):
        # TODO: add support for units other than SI:
        from sympy.physics.units.definitions import SI_quantity_dimension_map
        return SI_quantity_dimension_map[self]

    @property
    def abbrev(self):
        """
        Symbol representing the unit name.

        Prepend the abbreviation with the prefix symbol if it is defines.
        """
        return self._abbrev

    @property
    def scale_factor(self):
        """
        Overall magnitude of the quantity as compared to the canonical units.
        """
        from sympy.physics.units.definitions import SI_quantity_scale_factors
        return SI_quantity_scale_factors.get(self, S.One)

    def _eval_is_positive(self):
        return self.scale_factor.is_positive

    def _eval_is_constant(self):
        return self.scale_factor.is_constant()

    def _eval_Abs(self):
        # FIXME prefer usage of self.__class__ or type(self) instead
        q = self.func(self.name, self.abbrev)
        # TODO: add support for unit systems other that SI:
        from sympy.physics.units.definitions import SI_quantity_scale_factors
        SI_quantity_scale_factors[q] = process_scale_factor(Abs(self.scale_factor))

    def _eval_subs(self, old, new):
        if isinstance(new, Quantity) and self != old:
            return self

    @staticmethod
    def get_dimensional_expr(expr):
        if isinstance(expr, Mul):
            return Mul(*[Quantity.get_dimensional_expr(i) for i in expr.args])
        elif isinstance(expr, Pow):
            return Quantity.get_dimensional_expr(expr.base) ** expr.exp
        elif isinstance(expr, Add):
            return Quantity.get_dimensional_expr(expr.args[0])
        elif isinstance(expr, Derivative):
            dim = Quantity.get_dimensional_expr(expr.expr)
            for independent, count in expr.variable_count:
                dim /= Quantity.get_dimensional_expr(independent)**count
            return dim
        elif isinstance(expr, Function):
            args = [Quantity.get_dimensional_expr(arg) for arg in expr.args]
            if all(i == 1 for i in args):
                return S.One
            return expr.func(*args)
        elif isinstance(expr, Quantity):
            return expr.dimension.name
        return S.One

    @staticmethod
    def _collect_factor_and_dimension(expr):
        """Return tuple with factor expression and dimension expression."""
        if isinstance(expr, Quantity):
            return expr.scale_factor, expr.dimension
        elif isinstance(expr, Mul):
            factor = 1
            dimension = Dimension(1)
            for arg in expr.args:
                arg_factor, arg_dim = Quantity._collect_factor_and_dimension(arg)
                factor *= arg_factor
                dimension *= arg_dim
            return factor, dimension
        elif isinstance(expr, Pow):
            factor, dim = Quantity._collect_factor_and_dimension(expr.base)
            exp_factor, exp_dim = Quantity._collect_factor_and_dimension(expr.exp)
            if exp_dim.is_dimensionless:
               exp_dim = 1
            return factor ** exp_factor, dim ** (exp_factor * exp_dim)
        elif isinstance(expr, Add):
            factor, dim = Quantity._collect_factor_and_dimension(expr.args[0])
            for addend in expr.args[1:]:
                addend_factor, addend_dim = \
                    Quantity._collect_factor_and_dimension(addend)
                if dim != addend_dim:
                    raise ValueError(
                        'Dimension of "{0}" is {1}, '
                        'but it should be {2}'.format(
                            addend, addend_dim.name, dim.name))
                factor += addend_factor
            return factor, dim
        elif isinstance(expr, Derivative):
            factor, dim = Quantity._collect_factor_and_dimension(expr.args[0])
            for independent, count in expr.variable_count:
                ifactor, idim = Quantity._collect_factor_and_dimension(independent)
                factor /= ifactor**count
                dim /= idim**count
            return factor, dim
        elif isinstance(expr, Function):
            fds = [Quantity._collect_factor_and_dimension(
                arg) for arg in expr.args]
            return (expr.func(*(f[0] for f in fds)),
                    expr.func(*(d[1] for d in fds)))
        elif isinstance(expr, Dimension):
            return 1, expr
        else:
            return expr, Dimension(1)

    def convert_to(self, other):
        """
        Convert the quantity to another quantity of same dimensions.

        Examples
        ========

        >>> from sympy.physics.units import speed_of_light, meter, second
        >>> speed_of_light
        speed_of_light
        >>> speed_of_light.convert_to(meter/second)
        299792458*meter/second

        >>> from sympy.physics.units import liter
        >>> liter.convert_to(meter**3)
        meter**3/1000
        """
        from .util import convert_to
        return convert_to(self, other)

    @property
    def free_symbols(self):
        """Return free symbols from quantity."""
        return self.scale_factor.free_symbols


def _Quantity_constructor_postprocessor_Add(expr):
    # Construction postprocessor for the addition,
    # checks for dimension mismatches of the addends, thus preventing
    # expressions like `meter + second` to be created.

    deset = {
        tuple(sorted(dimsys_default.get_dimensional_dependencies(
            Dimension(Quantity.get_dimensional_expr(i) if not i.is_number else 1
        )).items()))
        for i in expr.args
        if i.free_symbols == set()  # do not raise if there are symbols
                    # (free symbols could contain the units corrections)
    }
    # If `deset` has more than one element, then some dimensions do not
    # match in the sum:
    if len(deset) > 1:
        raise ValueError("summation of quantities of incompatible dimensions")
    return expr


Basic._constructor_postprocessor_mapping[Quantity] = {
    "Add" : [_Quantity_constructor_postprocessor_Add],
}


def process_dimension(dimension):
    from sympy.physics.units.dimensions import dimsys_default, DimensionSystem
    # TODO: add support for more dimension systems:
    dim_sys = dimsys_default

    if not isinstance(dimension, dimensions.Dimension):
        if dimension == 1:
            dimension = Dimension(1)
        else:
            raise ValueError("expected dimension or 1")
    else:
        for dim_sym in dimension.name.atoms(Dimension):
            if dim_sym not in [i.name for i in dim_sys._dimensional_dependencies]:
                raise ValueError("Dimension %s is not registered in the "
                                 "dimensional dependency tree." % dim_sym)
    return dimension


def process_scale_factor(scale_factor):
    scale_factor = sympify(scale_factor)
    # replace all prefixes by their ratio to canonical units:
    scale_factor = scale_factor.replace(lambda x: isinstance(x, Prefix), lambda x: x.scale_factor)
    # replace all quantities by their ratio to canonical units:
    scale_factor = scale_factor.replace(lambda x: isinstance(x, Quantity), lambda x: x.scale_factor)
    return scale_factor
