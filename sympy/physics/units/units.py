# -*- coding:utf-8 -*-

"""
Units for physical quantities.

Also define prefixes.
We speak of canonical basis concerning the original dimensions/units/etc. used
to define all objects before the existence of the unit system.
"""

from __future__ import division

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

    def __init__(self, abbrev, dimension):
        self.abbrev = abbrev
        self.dimension = dimension

class UnitSystem():

    def __init__(self, base, units=None):
        self._base_units = base
        self._units = units

        self._list_dim()
        self._compute_dim_matrix()

        if self._system_is_well_defined() is False:
            # TODO: check precisely which test failed
            raise AttributeError('The base units can not well-defined an unit '
                                 'system.')

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
