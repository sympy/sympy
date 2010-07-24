"""Implementation of :class:`Ring` class. """

from sympy.polys.domains.domain import Domain
from sympy.polys.polyerrors import ExactQuotientFailed, NotInvertible, NotReversible

class Ring(Domain):
    """Represents a ring domain. """

    has_Ring = True

    def get_ring(self):
        """Returns a ring associated with `self`. """
        return self

    def exquo(self, a, b):
        """Exact quotient of `a` and `b`, implies `__floordiv__`.  """
        return a // b

    def quo(self, a, b):
        """Quotient of `a` and `b`, implies `__floordiv__`. """
        if a % b:
            raise ExactQuotientFailed(a, b, self)
        else:
            return a // b

    def rem(self, a, b):
        """Remainder of `a` and `b`, implies `__mod__`.  """
        return a % b

    def div(self, a, b):
        """Division of `a` and `b`, implies `__divmod__`. """
        return divmod(a, b)

    def invert(self, a, b):
        """Returns inversion of `a mod b`. """
        s, t, h = self.gcdex(a, b)

        if self.is_one(h):
            return s % b
        else:
            raise NotInvertible("zero divisor")

    def revert(self, a):
        """Returns `a**(-1)` if possible. """
        if self.is_one(a):
            return a
        else:
            raise NotReversible('only unity is reversible in a ring')

    def numer(self, a):
        """Returns numerator of `a`. """
        return a

    def denom(self, a):
        """Returns denominator of `a`. """
        return self.one

