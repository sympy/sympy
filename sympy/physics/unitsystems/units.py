# -*- coding: utf-8 -*-

"""
Unit system for physical quantities.
"""

from __future__ import division

from sympy import sympify, AtomicExpr
from .dimensions import Dimension


class Unit(AtomicExpr):
    """
    Class for the units.
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

        This factor represents the position of the unit with respect to the
        canonical unit of this dimension. For example if we choose the gram
        to be the canonical dimension for the mass, then by definition its
        factor is 1; on the other hand the factor defined here for kilogram is
        1000, even when it is a base unit. The explanation is that here we do
        not have defined any system that we could use as a reference: here the
        canonical unit is the only scale, and thus the only available origin.
        """

        if self.prefix is not None:
            return self.prefix.factor * self._factor
        else:
            return self._factor
