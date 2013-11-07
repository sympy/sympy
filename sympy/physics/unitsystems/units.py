# -*- coding: utf-8 -*-

"""
Unit system for physical quantities.
"""

from __future__ import division
import numbers

from sympy import sympify, AtomicExpr, Number, Pow, Mul
from .dimensions import Dimension, DimensionSystem


class Unit(AtomicExpr):
    """
    Class for the units.

    A unit is defined by two things:

    - a dimension;
    - a factor.

    The factor represents the position of the unit with respect to the
    canonical unit of this dimension. For example if we choose the gram to be
    the canonical dimension for the mass, then by definition its factor is 1;
    on the other hand the factor defined here for kilogram is 1000, even when
    it is a base unit. The explanation is that here we do not have defined any
    system that we could use as a reference: here the canonical unit is the
    only scale, and thus the only available origin.

    Additionnaly one can add a prefix and an abbreviation. The only utility of
    the former is to provide a shorthand for some units, but it is never used
    among computations; it appears only when defining and priting units. The
    same reamrl applies to the abbreviation.

    All operations (pow, mul, etc.) are defined as the corresponding ones
    acting on the factor (a number) and the dimension.
    """

    is_commutative = True
    is_number = False
    # make sqrt(m**2) --> m
    is_positive = True

    def __new__(cls, dim, abbrev="", factor=1, prefix=None, **assumptions):
        """
        Create a new unit instance.

        ``dim`` can be a Dimension or Unit object. The latter allows to
        construct derived units and constants. Note that the argument prefix
        is ignored if ``dim`` is a Unit instance.
        """

        factor = sympify(factor)

        obj = AtomicExpr.__new__(cls, **assumptions)

        if isinstance(dim, Dimension):
            obj._abbrev = abbrev
            obj._factor = factor
            obj.dim = dim
            obj.prefix = prefix
        elif isinstance(dim, Unit):
            obj._abbrev = abbrev or dim.abbrev
            obj._factor = factor * dim.factor
            obj.dim = dim.dim
            obj.prefix = None
        else:
            raise TypeError("'dim' object should be Unit or Dimension "
                            "instance; found %s" % type(dim))

        return obj

    @property
    def abbrev(self):
        """
        Symbol representing the unit name.

        Prepend the abbreviation with the prefix symbol if it is defines.
        """

        if self._abbrev == "":
            return ""
        if self.prefix is not None:
            return self.prefix.abbrev + self._abbrev
        else:
            return self._abbrev

    @property
    def abbrev_dim(self):
        """
        Abbreviation which use only intrinsinc properties of the unit.
        """

        return '(%g %s)' % (self.factor, self.dim)

    def __str__(self):
        if self.abbrev != "":
            return self.abbrev
        else:
            return self.abbrev_dim

    def __repr__(self):
        return self.abbrev_dim

    @property
    def factor(self):
        """
        Overall magnitude of the unit.
        """

        if self.prefix is not None:
            return self.prefix.factor * self._factor
        else:
            return self._factor

    def is_compatible(self, other):
        """
        Test if argument is a unit and has the same dimension as self.

        This function is used to verify that some operations can be done.
        """

        if isinstance(other, Unit):
            if self.dim == other.dim:
                return True

        return False

    def __eq__(self, other):
        return (isinstance(other, Unit) and self.factor == other.factor
                and self.dim == other.dim)

    def __add__(self, other):
        if not isinstance(other, Unit):
            raise TypeError("Only unit can be added; '%s' is not valid"
                            % type(other))
        else:
            if self.is_compatible(other):
                return Unit(self.dim, factor=self.factor + other.factor)
            else:
                raise ValueError("Only dimension which are equal can be "
                                 "added; '%s' and '%s' are different"
                                 % (self, other))

    def __sub__(self, other):

        if not isinstance(other, Unit):
            raise TypeError("Only unit can be added; '%s' is not valid"
                            % type(other))
        else:
            if self.is_compatible(other):
                return Unit(self.dim, factor=self.factor - other.factor)
            else:
                raise ValueError("Only dimension which are equal can be "
                                 "subtracted; '%s' and '%s' are different"
                                 % (self, other))

    def __pow__(self, other):

        other = sympify(other)
        #TODO: check consistency when having rational, float...
        if isinstance(other, (numbers.Real, Number)):
            if other == 0:
                return sympify(1)
            elif other == 1:
                return self
            else:
                factor = (self.factor**other).evalf()
                dim = self.dim**other
                if dim == 1:
                    return factor
                else:
                    return Unit(dim, factor=factor)
        else:
            return Pow(self, other)

    def __mul__(self, other):
        other = sympify(other)

        if other == 1:
            return self
        elif isinstance(other, Unit):
            factor = self.factor * other.factor
            dim = self.dim * other.dim
            if dim == 1:
                return factor
            else:
                return Unit(dim, factor=factor)
        #TODO: what to do when other is a number? return a unit or a quantity?
        #      or nothing special?
        #elif isinstance(other, Number):
        #    return Quantity(other, self)
        #elif other.is_number:
        #    factor = self.factor * other
        #    return Unit(self.dim, factor=factor)
        else:
            return Mul(self, other)

    def __rmul__(self, other):
        return self * other

    def __div__(self, other):
        other = sympify(other)

        if other == 1:
            return self
        elif isinstance(other, Unit):
            factor = self.factor / other.factor
            dim = self.dim / other.dim
            if dim == 1:
                return factor
            else:
                return Unit(dim, factor=factor)
        #TODO same remark as in __mul__
        #elif isinstance(other, Number):
        #    return Quantity(1/other, unit)
        #elif other.is_number:
            #factor = self.factor / other
            #return Unit(self.dimension, factor=factor, system=system)
        else:
            return Mul(self, Pow(other, -1))

    __truediv__ = __div__

    def __rdiv__(self, other):

        return self**-1 * other

    __rtruediv__ = __rdiv__


class UnitSystem(object):
    """
    UnitSystem represents a coherent set of units.

    A unit system is basically a dimension system with notions of scales. Many
    of the methods are defined in the same way.

    It is much better if all base units have a symbol.
    """

    def __init__(self, base, units=(), name="", descr=""):
        self.name = name
        self.descr = descr


        self._base_units = base
        self._units = tuple(list(units) + [u for u in base if u not in units])

        # construct the associated dimension system
        self._system = DimensionSystem([u.dim for u in base],
                                       [u.dim for u in units])

        if self.is_consistent is False:
            raise ValueError("The system with basis '%s' is not consistent"
                             % str(self._base_units))

    @property
    def is_consistent(self):
        return self._system.is_consistent
