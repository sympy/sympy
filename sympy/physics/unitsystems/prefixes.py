# -*- coding: utf-8 -*-

"""
Module defining unit prefixe class and some constants.

Constant dict for SI and binary prefixes are defined as PREFIXES and
BIN_PREFIXES.
"""

from sympy import sympify, Expr


class Prefix(Expr):
    """
    This class represent prefixes, with their name, symbol and factor.

    Prefixes are used to create derived units from a given unit. They should
    always be encapsulated into units.

    The factor is constructed from a base (default is 10) to some power, and
    it gives the total multiple or fraction. For example the kilometer km
    is constructed from the meter (factor 1) and the kilo (10 to the power 3,
    i.e. 1000). The base can be changed to allow e.g. binary prefixes.

    A prefix multiplied by something will always return the product of this
    other object times the factor, except if the other object:

    - is a prefix and they can be combined into a new prefix;
    - defines multiplication with prefixes (which is the case for the Unit
      class).
    """
    _op_priority = 13.0

    def __new__(cls, name, abbrev, exponent, base=sympify(10)):

        name = sympify(name)
        abbrev = sympify(abbrev)
        exponent = sympify(exponent)
        base = sympify(base)

        obj = Expr.__new__(cls, name, abbrev, exponent, base)
        obj._name = name
        obj._abbrev = abbrev
        obj._factor = base**exponent
        return obj

    @property
    def name(self):
        return self._name

    @property
    def abbrev(self):
        return self._abbrev

    def __str__(self):
        if self.base == 10:
            return "Prefix(%s, %s, %s)" % (self.name, self.abbrev, self.exponent)
        else:
            return "Prefix(%s, %s, %s, %s)" % (self.name, self.abbrev, self.exponent, self.base)

    __repr__ = __str__

    def __mul__(self, other):
        fact = self.factor * other.factor

        if fact == 1:
            return 1
        elif isinstance(other, Prefix):
            # simplify prefix
            for p in PREFIXES:
                if PREFIXES[p].factor == fact:
                    return PREFIXES[p]
            return fact

        return self.factor * other

    def __div__(self, other):
        fact = self.factor / other.factor

        if fact == 1:
            return 1
        elif isinstance(other, Prefix):
            for p in PREFIXES:
                if PREFIXES[p].factor == fact:
                    return PREFIXES[p]
            return fact

        return self.factor / other

    __truediv__ = __div__

    def __rdiv__(self, other):
        if other == 1:
            for p in PREFIXES:
                if PREFIXES[p].factor == 1 / self.factor:
                    return PREFIXES[p]
        return other / self.factor

    __rtruediv__ = __rdiv__


def prefix_unit(unit, prefixes):
    """
    Return a list of all units formed by unit and the given prefixes.

    You can use the predefined PREFIXES or BIN_PREFIXES, but you can also
    pass as argument a subdict of them if you don't want all prefixed units.

        >>> from sympy.physics.unitsystems.prefixes import (PREFIXES,
        ...                                                 prefix_unit)
        >>> from sympy.physics.unitsystems.systems import mks
        >>> m = mks["m"]
        >>> pref = {"m": PREFIXES["m"], "c": PREFIXES["c"], "d": PREFIXES["d"]}
        >>> prefix_unit(m, pref)  #doctest: +SKIP
        [cm, dm, mm]
    """

    from sympy.physics.unitsystems.units import Unit

    prefixed_units = []

    for prefix_abbr, prefix in prefixes.items():
        prefixed_units.append(Unit("%s%s" % (prefix.name, unit.name), unit.dimension, unit.factor, unit.abbrev,
                                   prefix=prefix))

    return prefixed_units


# http://physics.nist.gov/cuu/Units/prefixes.html
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

# http://physics.nist.gov/cuu/Units/binary.html
BIN_PREFIXES = {
    'Ki': Prefix('kibi', 'Y', 10, 2),
    'Mi': Prefix('mebi', 'Y', 20, 2),
    'Gi': Prefix('gibi', 'Y', 30, 2),
    'Ti': Prefix('tebi', 'Y', 40, 2),
    'Pi': Prefix('pebi', 'Y', 50, 2),
    'Ei': Prefix('exbi', 'Y', 60, 2)
}
