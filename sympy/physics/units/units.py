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

from sympy import (sympify, Number, Integer, Rational, Matrix, AtomicExpr,
                   Mul, Pow)

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
    if system is None or isinstance(system, UnitSystem):
        _UNIT_SYSTEM = system
    else:
        raise TypeError('The argument should be None or UnitSystem instance.')

def get_system():
    return _UNIT_SYSTEM

class Prefix(object):
    """
    This class represent prefixes, with their name, symbol and factor.
    """

    def __init__(self, name, abbrev, exponent):

        self.name = name
        self.abbrev = abbrev

        self.factor = sympify(10)**exponent

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

    is_positive = True # make sqrt(m**2) --> m
    is_commutative = True
    is_number = False

    def __new__(cls, dimension, abbrev='', factor=1, prefix=None, system=None,
                **assumptions):

        obj = AtomicExpr.__new__(cls, **assumptions)
        obj._abbrev = abbrev
        obj._factor = sympify(factor)
        obj.dimension = dimension
        obj.prefix = prefix

        #used to check if the unit is part of an unit system, and if it's one of
        #the base unit in it; should not be modified by hand
        obj._system = system

        return obj

    @property
    def factor(self):
        if self.prefix is not None:
            return self.prefix.factor * self._factor
        else:
            return self._factor

    @property
    def abbrev(self):
        if self._abbrev == '':
            return ''
        if self.prefix is not None:
            return self.prefix.abbrev + self._abbrev
        else:
            return self._abbrev

    @property
    def abbrev_base(self):
        '''
        Compute the abbreviation in terms of base units.
        '''
        if self.has_system is False or _UNIT_SYSTEM is None:
            return ''

        string = ''
        syst = self._system

        ratio = self.factor
        # sort by order of decreasing power
        l = zip(syst._base_units, syst.dim_vector(self.dimension))
        for u, d in sorted(l, key=lambda x: x[1], reverse=True):
            if d == 0:
                continue
            elif d == 1:
                ratio /= u.factor
                string += '%s ' % u
            elif d != 0 and d != 1:
                string += '%s**%s ' % (u, d)
                ratio /= u.factor**d

        if ratio != 1:
            string = '%s * %s' % (ratio, string)

        return string.strip()
        # cache result?

    def __str__(self):
        if self.abbrev != '':
            return self.abbrev
        elif self.abbrev_base != '':
            return self.abbrev_base
        else:
            return '%s %s' % (self.factor, self.dimension)

    __repr__ = __str__

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

    @staticmethod
    def _comptute_unit(unit, system):
        '''
        Given an unit and a system, try to find the unit in the system.
        '''
        if system is None:
            return unit
        else:
            u = system.get_unit(unit)
            if u is not None:
                return u
            else:
                return unit

    def __pow__(self, other):
        system = _UNIT_SYSTEM or self._system

        other = sympify(other)
        #TODO: check consistency when having rational, float...
        if isinstance(other, Number):
            if other == 0:
                return sympify(1)
            else:
                factor = self.factor**other
                dim = self.dimension**other
                if dim == 1:
                    return factor
                else:
                    unit = Unit(dim, factor=factor, system=system)
                    return self._comptute_unit(unit, system)
        else:
            return Pow(self, other)

    def __mul__(self, other):
        system = _UNIT_SYSTEM or self._system

        if isinstance(other, Unit):
            system = system or other._system

            factor = self.factor * other.factor
            dim = self.dimension * other.dimension
            if dim == 1:
                return factor
            else:
                unit = Unit(dim, factor=factor, system=system)
                return self._comptute_unit(unit, system)
        else:
            return Mul(self, other)

    def __div__(self, other):
        system = _UNIT_SYSTEM or self._system

        if isinstance(other, Unit):
            system = system or other._system

            factor = self.factor / other.factor
            dim = self.dimension / other.dimension
            if dim == 1:
                return factor
            else:
                unit = Unit(dim, factor=factor, system=system)
                return self._comptute_unit(unit, system)
        else:
            return Mul(self, Pow(other, -1))

    __truediv__ = __div__

    #TODO: should m+m be interpreted as 1*m + 1*m?
    #def __add__(self, other):
    #    if self != other:
    #        raise TypeError('Only identical units can be added.')
    #    else:
    #        return self

    #def __sub__(self, other):
    #    if self != other:
    #        raise TypeError('Only identical units can be subtracted.')
    #    else:
    #        return self


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
        self._list_dim = sorted(gen.keys())

    @property
    def base_factor(self):
        factors = [unit.factor for unit in self._base_units]
        return reduce(lambda x,y: x*y, factors)

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

    def __getitem__(self, key):
        '''
        Shortcut to the get_unit method, using key access.
        '''
        u = self.get_unit(key)
        if u is None:
            raise KeyError(key)
        return u

    def dim_vector(self, dimension):
        """
        Vector representation in terms of the base dimensions.
        """

        return self._transf_matrix * self.can_dim_vector(dimension)

    def _system_is_well_defined(self):
        #TODO: check redundancy between base units, i.e. if we can invert the
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
                    #TODO: use copy instead of direct assignment?
                    found_unit = u
        elif isinstance(unit, Unit):
            try:
                i = self._units.index(unit)
                #TODO: use copy instead of direct assignment?
                found_unit = self._units[i]
            except ValueError:
                pass
        if found_unit is not None:
            found_unit._system = self

        return found_unit

class Quantity(AtomicExpr):

    is_commutative = True

    def __new__(cls, factor=1, unit=None, **assumptions):

        factor = sympify(factor)

        if (unit is None or unit == 1) and isinstance(factor, Number):
            return factor

        obj = AtomicExpr.__new__(cls, **assumptions)

        #TODO: if factor is of the form "1 m", use the current defined system
        #      to get the unit
        if isinstance(factor, Number):
            if isinstance(unit, Unit):
                obj.unit = unit
                obj.factor = factor
            else:
                raise TypeError('"unit" should be an Unit instance.')
        else:
            raise NotImplementedError

        return obj

    def __str__(self):
        return '%s %s' % (self.factor, self.unit)

    __repr__ = __str__

    def __neg__(self):
        return Quantity(-self.factor, self.unit)

    def __add__(self, other):

        if isinstance(other, Quantity):
            if self.unit == other.unit:
                factor = self.factor + other.factor
                return Quantity(factor, self.unit)
            else:
                raise ValueError('Units should be the same to add quantities.')
        elif isinstance(other, Unit):
            if self.unit == other:
                factor = self.factor + 1
                return Quantity(factor, self.unit)
            else:
                raise ValueError('Units should be the same to add quantities.')
        else:
            raise TypeError('Only quantities and units can be added.')

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):

        if isinstance(other, Quantity):
            if self.unit == other.unit:
                return Quantity(self.factor - other.factor, self.unit)
            else:
                raise ValueError('Units should be the same to subtract '
                                 'quantities.')
        elif isinstance(other, Unit):
            if self.unit == other:
                return Quantity(self.factor - 1, self.unit)
            else:
                raise ValueError('Units should be the same to subtract '
                                 'quantities.')
        else:
            raise TypeError('Only quantities and units can be subtracted.')

    def __rsub__(self, other):
       return -self + other

    def __mul__(self, other):

        other = sympify(other)
        if isinstance(other, Quantity):
            return Quantity(self.factor*other.factor, self.unit*other.unit)
        elif isinstance(other, Unit):
            return Quantity(self.factor, self.unit*other)
        elif isinstance(other, Number):
            return Quantity(self.factor*other, self.unit)
        else:
            return Mul(self, other)

    def __rmul__(self, other):

        return self*other

    def __div__(self, other):

        other = sympify(other)
        if isinstance(other, Quantity):
            return Quantity(self.factor/other.factor, self.unit/other.unit)
        elif isinstance(other, Unit):
            return Quantity(self.factor, self.unit/other)
        elif isinstance(other, Number):
            return Quantity(self.factor/other, self.unit)
        else:
            return Mul(self, other)

    __truediv__ = __div__

    def __rdiv__(self, other):

        other = sympify(other)
        if isinstance(other, Quantity):
            return Quantity(other.factor/self.factor, other.unit/self.unit)
        elif isinstance(other, Unit):
            return Quantity(1/self.factor, other/self.unit)
        elif isinstance(other, Number):
            return Quantity(other/self.factor, self.unit**-1)
        else:
            return Mul(self, other)

    __rtruediv__ = __rdiv__

    def __pow__(self, other):

        other = sympify(other)
        if isinstance(other, Number):
            return Quantity(self.factor**other, self.unit**other)
        else:
            return Pow(self, other)

    def __eq__(self, other):
        #TODO: interpret an Unit as a quantity with factor 1
        return (isinstance(other, Quantity) and self.factor == other.factor
                and self.unit == other.unit)
