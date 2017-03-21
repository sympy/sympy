# -*- coding: utf-8 -*-

"""
Physical quantities.
"""

from __future__ import division
import numbers

from sympy.core.compatibility import string_types
from sympy import sympify, Expr, Number, Mul, Pow, S, Symbol, Add
from sympy.physics.unitsystems import dimensions


class Quantity(Expr):
    """
    Physical quantity.
    """

    is_commutative = True
    is_number = False

    def __new__(cls, name, dimension, factor=S.One, abbrev=None, prefix=None, **assumptions):

        if not isinstance(name, Symbol):
            name = Symbol(name)

        if not isinstance(dimension, dimensions.Dimension):
            dimension = getattr(dimensions, str(dimension))
        factor = sympify(factor)

        if abbrev is None:
            abbrev = name
        elif isinstance(abbrev, string_types):
            abbrev = Symbol(abbrev)

        if prefix is None:
            obj = Expr.__new__(cls, name, dimension, factor, abbrev)
        else:
            prefix = sympify(prefix)
            obj = Expr.__new__(cls, name, dimension, factor, abbrev, prefix)
        obj._name = name
        obj._dimension = dimension
        obj._factor_without_prefix = factor
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
    def factor(self):
        """
        Overall magnitude of the quantity as compared to the canonical units.
        """
        return self._factor_without_prefix

    def __str__(self):
        return "%s" % (self.name)

    def __repr__(self):
        return self.__str__()

    @staticmethod
    def _collect_factor_and_dimension(expr):

        if isinstance(expr, Quantity):
            return expr.factor, expr.dimension
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

        >>> from sympy.physics.unitsystems import speed_of_light, meter, second
        >>> speed_of_light
        speed_of_light
        >>> speed_of_light.convert_to(meter/second)
        299792458*meter/second

        >>> from sympy.physics.unitsystems import liter
        >>> liter.convert_to(meter**3)
        meter**3/1000
        """
        factor, dimension = self._collect_factor_and_dimension(other)

        if self.dimension != dimension:
            return self

        return self.factor/factor*other
