"""Implementation of :class:`Field` class. """

from sympy.polys.domains.ring import Ring
from sympy.polys.polyerrors import NotReversible, DomainError

class Field(Ring):
    """Represents a field domain. """

    has_Field = True

    def get_ring(self):
        """Returns a ring associated with ``self``. """
        raise DomainError('there is no ring associated with %s' % self)

    def get_field(self):
        """Returns a field associated with ``self``. """
        return self

    def exquo(self, a, b):
        """Exact quotient of ``a`` and ``b``, implies ``__div__``.  """
        return a / b

    def quo(self, a, b):
        """Quotient of ``a`` and ``b``, implies ``__div__``. """
        return a / b

    def rem(self, a, b):
        """Remainder of ``a`` and ``b``, implies nothing.  """
        return self.zero

    def div(self, a, b):
        """Division of ``a`` and ``b``, implies ``__div__``. """
        return a / b, self.zero

    def gcd(self, a, b):
        """Returns GCD of ``a`` and ``b`` that is consistent with the core
        implementation. In addition, this allows the primitive of an
        expression to be cleared of Rationals.

        >>> from sympy.polys.domains import QQ
        >>> from sympy import gcd, Rational, primitive
        >>> from sympy.abc import x
        >>> QQ.gcd(QQ(2, 3), QQ(4, 9))
        2/9
        >>> gcd(Rational(2, 3), Rational(4, 9))
        2/9
        >>> primitive(2*x/3 + Rational(4, 9))
        (2/9, 3*x + 2)

        """
        try:
            ring = self.get_ring()
        except DomainError:
            return self.one

        p = ring.gcd(self.numer(a), self.numer(b))
        q = ring.lcm(self.denom(a), self.denom(b))

        return self.convert(p, ring)/q

    def lcm(self, a, b):
        """Returns LCM of ``a`` and ``b`` that is consistent with the core
        implementation.

        >>> from sympy.polys.domains import QQ
        >>> from sympy import lcm, Rational, primitive
        >>> from sympy.abc import x
        >>> QQ.lcm(QQ(2, 3), QQ(4, 9))
        4/3
        >>> lcm(Rational(2, 3), Rational(4, 9))
        4/3
        """

        try:
            ring = self.get_ring()
        except DomainError:
            return a*b

        p = ring.lcm(self.numer(a), self.numer(b))
        q = ring.gcd(self.denom(a), self.denom(b))

        return self.convert(p, ring)/q

    def revert(self, a):
        """Returns ``a**(-1)`` if possible. """
        if a:
            return 1/a
        else:
            raise NotReversible('zero is not reversible')
