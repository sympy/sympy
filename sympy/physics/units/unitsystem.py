# -*- coding: utf-8 -*-

"""
Unit system for physical quantities; include definition of constants.
"""

from __future__ import division

from sympy.core.decorators import deprecated
from sympy.physics.units.quantities import Quantity

from sympy import S
from .dimensions import DimensionSystem


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

        # construct the associated dimension system
        self._system = DimensionSystem([u.dimension for u in base],
                                       [u.dimension for u in units])

        if self.is_consistent is False:
            raise ValueError("The system with basis '%s' is not consistent"
                             % str(self._base_units))

        self._units = tuple(set(base) | set(units))

        # create a dict linkin
        # this is possible since we have already verified that the base units
        # form a coherent system
        base_dict = dict((u.dimension, u) for u in base)
        # order the base units in the same order than the dimensions in the
        # associated system, in order to ensure that we get always the same
        self._base_units = tuple(base_dict[d] for d in self._system._base_dims)

    def __str__(self):
        """
        Return the name of the system.

        If it does not exist, then it makes a list of symbols (or names) of
        the base dimensions.
        """

        if self.name != "":
            return self.name
        else:
            return "UnitSystem((%s))" % ", ".join(str(d) for d in self._base_units)

    def __repr__(self):
        return '<UnitSystem: %s>' % repr(self._base_units)

    def __getitem__(self, key):
        """
        Shortcut to the get_unit method, using key access.
        """

        u = self.get_unit(key)

        #TODO: really want to raise an error?
        if u is None:
            raise KeyError(key)

        return u

    def extend(self, base, units=(), name="", description=""):
        """
        Extend the current system into a new one.

        Take the base and normal units of the current system to merge
        them to the base and normal units given in argument.
        If not provided, name and description are overriden by empty strings.
        """

        base = self._base_units + tuple(base)
        units = self._units + tuple(units)

        return UnitSystem(base, units, name, description)

    def print_unit_base(self, unit):
        """
        Give the string expression of a unit in term of the basis.

        Units are displayed by decreasing power.
        """

        res = S.One

        factor = unit.scale_factor
        vec = self._system.dim_vector(unit.dimension)

        for (u, p) in sorted(zip(self._base_units, vec), key=lambda x: x[1],
                             reverse=True):

            factor /= u.scale_factor ** p
            if p == 0:
                continue
            elif p == 1:
                res *= u
            else:
                res *= u**p

        return factor * res

    @property
    def dim(self):
        """
        Give the dimension of the system.

        That is return the number of units forming the basis.
        """

        return self._system.dim

    @property
    def is_consistent(self):
        """
        Check if the underlying dimension system is consistent.
        """
        return self._system.is_consistent
