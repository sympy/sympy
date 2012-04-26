# -*- coding:utf-8 -*-

"""
Units for physical quantities.

Also define prefixes.
We speak of canonical basis concerning the original dimensions/units/etc. used
to define all objects before the existence of the unit system.
"""

from __future__ import division
from copy import copy

from sympy import Rational, Matrix, AtomicExpr

# Is it a good idea to combine prefixes with between them, instead of just
# keeping them to better print results when the user asks for it?
# So do we want not to care about prefixes until the end?

class Prefix(object):
    """
    This class represent prefixes, with their name, symbol and factor.
    """

    def __init__(self, name, abbrev, exponent):

        self.name = name
        self.abbrev = abbrev

        self.factor = Rational(10)**exponent

    def __repr__(self):

        return self.name

    __str__ = __repr__

    def __mul__(self, other):
        fact = self.factor * other.factor
        if fact == 1:
            return 1
        elif isinstance(other, Prefix):
            # simplify prefix
            for p in PREFIXES:
                if PREFIXES[p].factor == fact:
                    return PREFIXES[p]
        #elif isinstance(other, Unit):
            #return PrefixedUnit(self, other)
        return self.factor * other

    def __div__(self, other):
        fact = self.factor / other.factor
        if fact == 1:
            return 1
        elif isinstance(other, Prefix):
            for p in PREFIXES:
                if PREFIXES[p].factor == fact:
                    return PREFIXES[p]
        return self.factor / other

    __truediv__ = __div__

    def __rdiv__(self, other):
        if other == 1:
            for p in PREFIXES:
                if PREFIXES[p].factor == 1/self.factor:
                    return PREFIXES[p]
        return other / self.factor

    __rtruediv__ = __rdiv__

PREFIXES = {
    'Y': Prefix('yotta', 'Y', 24),
    'Z': Prefix('zetta', 'Z', 21),
    'E': Prefix('exa', 'E', 18),
    'P': Prefix('peta', 'P', 15),
    'T': Prefix('tera', 'T', 12),
    'G': Prefix('giga', 'G', 9),
    'M': Prefix('mega', 'M', 6),
    'k': Prefix('kilo', 'k', 3),
    'h': Prefix('hecto', 'h', 2),
    'da': Prefix('deca', 'da', 1),
    'd': Prefix('deci', 'd', -1),
    'c': Prefix('centi', 'c', -2),
    'm': Prefix('milli', 'm', -3),
    'µ': Prefix('micro', 'µ', -6),
    'n': Prefix('nano', 'n', -9),
    'p': Prefix('pico', 'p', -12),
    'f': Prefix('femto', 'f', -15),
    'a': Prefix('atto', 'a', -18),
    'z': Prefix('zepto', 'z', -21),
    'y': Prefix('yocto', 'y', -24)
}

class Unit(AtomicExpr):

    #used to check if the unit is part of an unit system, and if it's one of
    #the base unit in it; should not be modified by hand
    _system = None

    def __init__(self, abbrev, dimension, value=1):
        self.abbrev = abbrev
        self._value = value
        self.dimension = dimension

    @property
    def value(self):
        return self._value

    @property
    def is_base_unit(self):
        """
        Check if the unit is part of the base unit of the system which is
        currently used.
        """
        if self._system is not None:
            if self in self._system._base_units:
                return True

        return False

    @property
    def has_system(self):
        if self._system is not None:
            return True

        return False


class UnitSystem():
    """
    Representation of an unit system.

    _base_units is a tuple of the units which form the basis. _units is also
    a tuple.
    """

    def __init__(self, base, units=(), name='', description=''):
        self.name = name
        self.description = description

        self._base_units = base
        self._units = list(units) + [u for u in base if u not in units]

        self._list_dim()
        self._compute_dim_matrix()

        if self._system_is_well_defined() is False:
            # TODO: check precisely which test failed
            raise AttributeError('The base units can not well-defined an unit '
                                 'system.')

    def __str__(self):
        return self.name

    def __repr__(self):
        return '<UnitSystem: %s>' % self.name

    def _list_dim(self):
        """
        Do a list of all the base dimensions, and sort them. It is then
        interpreted as a vector basis.
        """

        dims = [unit.dimension for unit in self._base_units]
        gen = reduce(lambda x,y: x*y, dims)
        self._list_dim = sorted(gen.keys())

    def _compute_dim_matrix(self):
        dim_tab = []
        for unit in self._base_units:
            row = []
            dim_tab.append(self.can_dim_vector(unit.dimension))
        self._dim_matrix = Matrix(dim_tab)

    def can_dim_vector(self, dimension):
        """
        Canonical vector representation of a dimension.
        """
        vec = []
        for dim in self._list_dim:
            vec.append(dimension.get(dim, 0))
        return vec

    def _system_is_well_defined(self):
        # check redundancy between base units, i.e. if we can invert the
        # matrices:
        if self._dim_matrix.det() == 0:
            return False

        return True

    def get_unit(self, unit):
        """
        Find a specific unit which is part of the system. Modify its attribute
        to say in which system it is.

        unit can be a string or an unit object.
        """
        found_unit = None

        if isinstance(unit, str):
            for u in self._units:
                if unit in (u.abbrev,):
                    found_unit = copy(u)

        if isinstance(unit, Unit):
            if unit in self._units:
                found_unit = copy(unit)

        if found_unit is not None:
            found_unit._system = self

        return found_unit

