# -*- coding:utf-8 -*-

"""
Units for physical quantities.

Also define prefixes.

The philosophy of this module is the following. We can consider an unit system
as a vector space on the integer, where the unit are the vector and
multiplication is the composition law. For example, if you choose the meter m
and the second s as the basic units, then you can express any unit as
u = m^a s^b where a and b are two integers, and it can be represented by the
column vector (a, b). Then the multiplication of two units u1 and u2 correspond
to the addition of the column vector. Exponentiation of units correspond to
scalar multiplication.

We speak of canonical basis concerning the original dimensions/units/etc. used
to define all objects before the existence of the unit system. Then we can
use linear algebra (transformation matrices...) to go from one system to
another one.
"""

from __future__ import division
from copy import copy

from sympy import Rational, Matrix, AtomicExpr, Mul

# Is it a good idea to combine prefixes with between them, instead of just
# keeping them to better print results when the user asks for it?
# So do we want not to care about prefixes until the end?

_UNIT_SYSTEM = None

def set_system(system):
    """
    Define the default unit system used in computations.

    Results are expressed in units of this system. If it occurs that quantities
    are expressed in an unit not defined in this system, it will be converted.

    In arithmetic operations, the unit system used is determined in this order:
    1. general system, 2. left unit system and 3. right unit system.
    """

    global _UNIT_SYSTEM
    _UNIT_SYSTEM = system


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

    def __init__(self, abbrev, dimension, factor=1):
        self.abbrev = abbrev
        self._factor = factor
        self.dimension = dimension

    @property
    def factor(self):
        return self._factor

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

    def __eq__(self, other):
        return (isinstance(other, Unit) and self.factor == other.factor
                and self.dimension == other.dimension)

    def __mul__(self, other):
        system = _UNIT_SYSTEM or self._system

        if isinstance(other, Unit):
            system = system or other.system

            factor = self.factor * other.factor
            dim = self.dimension * other.dimension
            unit = Unit(abbrev, dim, factor)

            if system is None:
                abbrev = '%s %s' % (self.abbrev, other.abbrev)
                return unit
            else:
                u = system.get_unit(unit)
                if u != []:
                    return u
                else:
                    system.can_dim_vector()

        else:
            return Mul(self, other)


class UnitSystem():
    """
    Representation of an unit system.

    _base_units is a tuple of the units which form the basis. _units is also
    a tuple.

    The base matrix is a multiple of the identity.
    """

    def __init__(self, base, units=(), name='', description=''):
        self.name = name
        self.description = description

        self._base_units = base
        self._units = list(units) + [u for u in base if u not in units]

        self._list_dim()
        self._compute_can_transf_matrix()

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
        self._list_dim = sorted(gen.names)

    def _compute_can_transf_matrix(self):
        """
        Compute the canonical transformation matrix from the canonical to the
        base unit basis.

        It's simply the matrix where columns are the vector of base units in
        canonical basis.
        """

        self._transf_matrix = reduce(lambda x, y: x.row_join(y),
                                     [self.can_dim_vector(unit.dimension)
                                            for unit in self._base_units])

    def can_dim_vector(self, dimension):
        """
        Canonical vector representation of a dimension.
        """
        vec = []
        for dim in self._list_dim:
            vec.append(dimension.get(dim, 0))
        return Matrix(vec)

    def dim_vector(self, dimension):
        """
        Vector representation in terms of the base dimensions.
        """

        return self._transf_matrix * self.can_dim_vector(dimension)

    def _system_is_well_defined(self):
        # check redundancy between base units, i.e. if we can invert the
        # matrices
        #if self._transf_matrix.det() == 0:
        if self._transf_matrix.is_square is False:
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
