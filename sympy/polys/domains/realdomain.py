"""Implementation of :class:`RealDomain` class. """

from sympy.polys.domains.simpledomain import SimpleDomain
from sympy.polys.domains.characteristiczero import CharacteristicZero

from sympy.polys.polyerrors import DomainError

import math

class RealDomain(CharacteristicZero, SimpleDomain): # XXX: should be a field
    """Abstract domain for real numbers. """

    rep   = 'RR'

    is_Exact     = False
    is_Numerical = True

    def as_integer_ratio(self, a, **args):
        """Convert real number to a (numer, denom) pair. """
        v, n = math.frexp(a) # XXX: hack, will work only for floats

        for i in xrange(300):
            if v != math.floor(v):
                v, n = 2*v, n-1
            else:
                break

        numer, denom = int(v), 1

        m = 1 << abs(n)

        if n > 0:
            numer *= m
        else:
            denom = m

        n, d = self.limit_denom(numer, denom, **args)

        if a and not n:
            return numer, denom
        else:
            return n, d

    def limit_denom(self, n, d, **args):
        """Find closest rational to `n/d` (up to `max_denom`). """
        max_denom = args.get('max_denom', 1000000)

        if d <= max_denom:
            return n, d

        from sympy.polys.domains import QQ
        self = QQ(n, d)

        p0, q0, p1, q1 = 0, 1, 1, 0

        while True:
            a  = n//d
            q2 = q0 + a*q1

            if q2 > max_denom:
                break

            p0, q0, p1, q1, n, d = \
                p1, q1, p0 + a*p1, q2, d, n - a*d

        k = (max_denom - q0)//q1

        P1, Q1 = p0 + k*p1, q0 + k*q1

        bound1 = QQ(P1, Q1)
        bound2 = QQ(p1, q1)

        if abs(bound2 - self) <= abs(bound1 - self):
            return p1, q1
        else:
            return P1, Q1

    def get_ring(self):
        """Returns a ring associated with `self`. """
        raise DomainError('there is no ring associated with %s' % self)

    def get_field(self):
        """Returns a field associated with `self`. """
        raise DomainError('there is no field associated with %s' % self)

    def get_exact(self):
        """Returns an exact domain associated with `self`. """
        from sympy.polys.domains import QQ
        return QQ

    def exquo(self, a, b):
        """Exact quotient of `a` and `b`, implies `__div__`.  """
        return a / b

    def quo(self, a, b):
        """Quotient of `a` and `b`, implies `__div__`. """
        return a / b

    def rem(self, a, b):
        """Remainder of `a` and `b`, implies nothing.  """
        return self.zero

    def div(self, a, b):
        """Division of `a` and `b`, implies `__div__`. """
        return a / b, self.zero

