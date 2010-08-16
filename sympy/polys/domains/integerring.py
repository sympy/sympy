"""Implementation of :class:`IntegerRing` class. """

from sympy.polys.domains.ring import Ring
from sympy.polys.domains.simpledomain import SimpleDomain
from sympy.polys.domains.characteristiczero import CharacteristicZero

import math

class IntegerRing(Ring, CharacteristicZero, SimpleDomain):
    """General class for integer rings. """

    is_ZZ = True
    rep   = 'ZZ'

    is_Numerical = True

    has_assoc_Ring         = True
    has_assoc_Field        = True

    def get_field(self):
        """Returns a field associated with `self`. """
        from sympy.polys.domains import QQ
        return QQ

    def from_AlgebraicField(K1, a, K0):
        """Convert a `ANP` object to `dtype`. """
        if a.is_ground:
            return K1.convert(a.LC(), K0.dom)

    def log(self, a, b):
        """Returns b-base logarithm of `a`. """
        return self.dtype(math.log(a, b))

