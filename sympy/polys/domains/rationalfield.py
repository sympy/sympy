"""Implementation of :class:`RationalField` class. """

from sympy.polys.domains.field import Field
from sympy.polys.domains.simpledomain import SimpleDomain
from sympy.polys.domains.characteristiczero import CharacteristicZero

class RationalField(Field, CharacteristicZero, SimpleDomain):
    """General class for rational fields. """

    is_QQ = True
    rep   = 'QQ'

    is_Numerical = True

    has_assoc_Ring         = True
    has_assoc_Field        = True

    def get_ring(self):
        """Returns a ring associated with `self`. """
        from sympy.polys.domains import ZZ
        return ZZ

    def algebraic_field(self, *extension):
        """Returns an algebraic field, i.e. `QQ(alpha, ...)`. """
        from sympy.polys.domains import AlgebraicField
        return AlgebraicField(self, *extension)

    def from_AlgebraicField(K1, a, K0):
        """Convert a `ANP` object to `dtype`. """
        if a.is_ground:
            return K1.convert(a.LC(), K0.dom)

