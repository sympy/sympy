# -*- coding: utf-8 -*-

"""
Unit system for physical quantities; include definition of constants.
"""

from __future__ import division

from sympy.physics.unitsystems.prefixes import PREFIXES

from sympy.physics.unitsystems.quantities import Quantity

from sympy import sympify
from sympy.core.compatibility import string_types
from .dimensions import DimensionSystem


class Unit(Quantity):
    """
    Class for the units.

    A Unit is a quantity that support prefixes.
    """

    def __new__(cls, name, dimension, factor, abbrev="", prefix=None, **assumptions):
        """
        Create a new unit instance.

        ``dim`` can be a Dimension or Unit object. The latter allows to
        construct derived units and constants. Note that the argument prefix
        is ignored if ``dim`` is a Unit instance and already has a prefix.
        """
        factor = sympify(factor)

        if isinstance(prefix, string_types):
            prefix = PREFIXES[prefix]

        if prefix is not None:
            factor *= prefix
            prefix = sympify(prefix)

        obj = Quantity.__new__(cls, name, dimension, factor, abbrev, prefix, **assumptions)
        obj._prefix = prefix
        return obj

    @property
    def prefix(self):
        return self._prefix

    @property
    def abbrev(self):
        if self.prefix is not None:
            return "%s%s" % (self.prefix.abbrev, self.abbrev)
        else:
            return self._abbrev

    def is_compatible(self, other):
        """
        Test if argument is a unit and has the same dimension as self.

        This function is used to verify that some operations can be done.
        """

        if isinstance(other, Unit):
            if self.dim == other.dim:
                return True

        return False


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
            return "(%s)" % ", ".join(str(d) for d in self._base_units)

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

    def get_unit(self, unit):
        """
        Find a specific unit which is part of the system.

        unit can be a string or a dimension object. If no unit is found, then
        return None.
        """

        #TODO: if the argument is a list, return a list of all matching dims

        found_unit = None

        #TODO: use copy instead of direct assignment for found_dim?
        if isinstance(unit, str):
            for u in self._units:
                #TODO: verify not only abbrev
                if unit in (u.abbrev,):
                    found_unit = u
                    break
        elif isinstance(unit, Unit):
            try:
                i = self._units.index(unit)
                found_unit = self._units[i]
            except ValueError:
                pass

        return found_unit

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

        res = ""

        factor = unit.factor
        vec = self._system.dim_vector(unit.dim)

        for (u, p) in sorted(zip(self._base_units, vec), key=lambda x: x[1],
                             reverse=True):

            factor /= u.factor**p
            if p == 0:
                continue
            elif p == 1:
                res += "%s " % str(u)
            else:
                res += "%s^%d " % (str(u), p)

        return "%g %s" % (factor, res.strip())

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
