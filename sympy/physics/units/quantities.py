# -*- coding: utf-8 -*-

"""
Physical quantities.
"""

from __future__ import division

from sympy.core.compatibility import string_types
from sympy import Abs, sympify, Mul, Pow, S, Symbol, Add, AtomicExpr, Basic, Function
from sympy.physics.units import Dimension
from sympy.physics.units import dimensions
from sympy.physics.units.prefixes import Prefix


class Quantity(AtomicExpr):
    """
    Physical quantity.
    """

    is_commutative = True
    is_real = True
    is_number = False
    is_nonzero = True
    _diff_wrt = True

    def __new__(cls, name, dimension, scale_factor=S.One, abbrev=None, **assumptions):

        if not isinstance(name, Symbol):
            name = Symbol(name)

        if not isinstance(dimension, dimensions.Dimension):
            if dimension == 1:
                dimension = Dimension(1)
            else:
                raise ValueError("expected dimension or 1")
        scale_factor = sympify(scale_factor)

        dimex = Quantity.get_dimensional_expr(scale_factor)
        if dimex != 1:
            if dimension != Dimension(dimex):
                raise ValueError("quantity value and dimension mismatch")

        # replace all prefixes by their ratio to canonical units:
        scale_factor = scale_factor.replace(lambda x: isinstance(x, Prefix), lambda x: x.scale_factor)
        # replace all quantities by their ratio to canonical units:
        scale_factor = scale_factor.replace(lambda x: isinstance(x, Quantity), lambda x: x.scale_factor)

        if abbrev is None:
            abbrev = name
        elif isinstance(abbrev, string_types):
            abbrev = Symbol(abbrev)

        obj = AtomicExpr.__new__(cls, name, dimension, scale_factor, abbrev)
        obj._name = name
        obj._dimension = dimension
        obj._scale_factor = scale_factor
        obj._abbrev = abbrev
        return obj

    @property
    def name(self):
        return self._name

    @property
    def dimension(self):
        return self._dimension

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
        return self._scale_factor

    def _eval_is_positive(self):
       return self.scale_factor.is_positive

    def _eval_is_constant(self):
        return self.scale_factor.is_constant()

    def _eval_Abs(self):
        # FIXME prefer usage of self.__class__ or type(self) instead
        return self.func(self.name, self.dimension, Abs(self.scale_factor),
                         self.abbrev)

    @staticmethod
    def get_dimensional_expr(expr):
        if isinstance(expr, Mul):
            return Mul(*[Quantity.get_dimensional_expr(i) for i in expr.args])
        elif isinstance(expr, Pow):
            return Quantity.get_dimensional_expr(expr.base) ** expr.exp
        elif isinstance(expr, Add):
            return Quantity.get_dimensional_expr(expr.args[0])
        elif isinstance(expr, Function):
            fds = [Quantity.get_dimensional_expr(arg) for arg in expr.args]
            return expr.func(*fds)
        elif isinstance(expr, Quantity):
            return expr.dimension.name
        return 1

    @staticmethod
    def _collect_factor_and_dimension(expr):

        if isinstance(expr, Quantity):
            return expr.scale_factor, expr.dimension
        elif isinstance(expr, Mul):
            factor = 1
            dimension = 1
            for arg in expr.args:
                arg_factor, arg_dim = Quantity._collect_factor_and_dimension(arg)
                factor *= arg_factor
                dimension *= arg_dim
            return factor, dimension
        elif isinstance(expr, Pow):
            factor, dim = Quantity._collect_factor_and_dimension(expr.base)
            return factor ** expr.exp, dim ** expr.exp
        elif isinstance(expr, Add):
            raise NotImplementedError
        else:
            return 1, 1

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
        return set([])


def _Quantity_constructor_postprocessor_Add(expr):
    # Construction postprocessor for the addition,
    # checks for dimension mismatches of the addends, thus preventing
    # expressions like `meter + second` to be created.

    deset = {
        tuple(sorted(Dimension(
            Quantity.get_dimensional_expr(i) if not i.is_number else 1
        ).get_dimensional_dependencies().items()))
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
