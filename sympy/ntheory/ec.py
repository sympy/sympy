from sympy.abc import x, y
from sympy.core.numbers import oo
from sympy.core.relational import Eq
from sympy.polys.domains import FiniteField, QQ, RationalField
from sympy.solvers.solvers import solve

from .factor_ import divisors
from .residue_ntheory import sqrt_mod


class EllipticCurve():
    """
    Create the following Elliptic Curve over domain.

    `y^{2} + a_{1} x y + a_{3} y = x^{3} + a_{2} x^{2} + a_{4} x + a_{6}`

    The default domain is ``QQ``. If no coefficient ``a1``, ``a2``, ``a3``,
    it create curve as following form.

    `y^{2} = x^{3} + a_{4} x + a_{6}`

    Examples
    ========

    References
    ==========

    [1] J. Silverman "A Friendly Introduction to Number Theory" Third Edition
    [2] http://mathworld.wolfram.com/EllipticDiscriminant.html
    [3] G. Hardy, E. Wright "An Introduction to the Theory of Numbers" Sixth Edition

    """

    O = (0, 1, 0)

    def __init__(self, a4, a6, a1=0, a2=0, a3=0, domain=QQ):
        self._domain = domain
        # Calculate discriminant
        self._b2 = a1**2 + 4 * a2
        self._b4 = 2 * a4 + a1 * a3
        self._b6 = a3**2 + 4 * a6
        self._b8 = a1**2 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3**2 - a4**2
        self._discrim = self._domain(-self._b2**2 * self._b8 - 8 * self._b4**3 - 27 * self._b6**2 + 9 * self._b2 * self._b4 * self._b6)
        self._a1 = self._domain(a1)
        self._a2 = self._domain(a2)
        self._a3 = self._domain(a3)
        self._a4 = self._domain(a4)
        self._a6 = self._domain(a6)
        self._eq = Eq(y**2 + self._a1*x*y + self._a3*y, x**3 + self._a2*x**2 + self._a4*x + self._a6)
        if isinstance(self._domain, FiniteField):
            self._rank = 0
        elif isinstance(self._domain, RationalField):
            self._rank = None

    def __contains__(self, point):
        if self.characteristic == 0 and len(point) == 3 and point[2] == 0:
            return True
        return self._eq.subs({x: point[0], y: point[1]})

    def __repr__(self):
        return 'E({}): {}'.format(self._domain, self._eq)

    def add(self, p1, p2, to_sympy=True):
        """
        Return new point R = p1 + p2 in curve.

        Examples
        ========

        >>> from sympy.ntheory.ec import EllipticCurve
        >>> e1 = EllipticCurve(-17, 16)
        >>> e1.add((0, -4), (1, 0))
        (15, -56, 1)

        """
        if len(p1) == 3 and p1[2] == 0:
            return p2
        if len(p2) == 3 and p2[2] == 0:
            return p1
        x1, y1 = p1[:2]
        x2, y2 = p2[:2]
        if x1 != x2:
            slope = (y1 - y2) / (x1 - x2)
        else:
            if (y1 + y2) == 0:
                return self.O
            slope = (3 * x1**2 + self._a4) / (2 * y1)
        x3 = slope**2 - x1 - x2
        y3 = -y1 - slope * (x3 - x1)
        if to_sympy:
            try:
                return self._domain.to_sympy(x3), self._domain.to_sympy(y3), 1
            except TypeError:
                pass
        return x3, y3, 1

    def minimal(self):
        """
        Return minimal Weierstrass equation Elliptic Curve.

        Examples
        ========

        >>> from sympy.ntheory.ec import EllipticCurve
        >>> e1 = EllipticCurve(-10, -20, 0, -1, 1)
        >>> e1.minimal()
        E(QQ): y**2 == x**3 - 13392*x - 1080432

        """
        char = self.characteristic
        if char == 2:
            return self
        if char == 3:
            return EllipticCurve(self._b4/2, self._b6/4, a2=self._b2/4, domain=self._domain)
        c4 = self._b2**2 - 24*self._b4
        c6 = -self._b2**3 + 36*self._b2*self._b4 - 216*self._b6
        return EllipticCurve(-27*c4, -54*c6, domain=self._domain)

    def mul(self, p, n):
        """
        Return new point R = nP in curve.

        Compute nP using Double-and-Add method.

        Examples
        ========

        >>> from sympy.ntheory.ec import EllipticCurve
        >>> e3 = EllipticCurve(-1, 9)
        >>> e3.mul((1, -3), 3)
        (664/169, 17811/2197, 1)

        """
        if n < 1:
            return p
        r = self.O
        while n:
            if n & 1:
                r = self.add(r, p, to_sympy=False)
            n >>= 1
            p = self.add(p, p, to_sympy=False)
        try:
            return self._domain.to_sympy(r[0]), self._domain.to_sympy(r[1]), r[2]
        except TypeError:
            return r

    def point_order(self, p):
        """
        Return order of point p.

        Order of point is integer n that nP = O.
        """
        if p not in self:
            raise ValueError('Invalid point.')
        if p == self.O:
            return 1
        if p[1] == 0:
            return 2
        n = self.add(p, p)
        if n[1] == -p[1]:
            return 3
        i = 2
        while int(n[0]) == n[0]:
            n = self.add(p, n)
            i += 1
            if n == self.O:
                return i
        return oo

    def points(self):
        """
        Return points of curve over Finite Field.

        Examples
        ========

        >>> from sympy.polys.domains import FF
        >>> from sympy.ntheory.ec import EllipticCurve
        >>> e2 = EllipticCurve(1, 0, domain=FF(2))
        >>> list(e2.points())
        [(0, 0), (1, 0)]

        """
        char = self.characteristic
        if char > 1:
            for i in range(char):
                y = sqrt_mod(i**3 + self._a2*i**2 + self._a4*i + self._a6, char)
                if y is not None:
                    yield i, y
                    if y != 0:
                        yield i, char - y
        else:
            raise NotImplementedError("Still not implemented")

    def torsion_points(self):
        """
        Return torsion points of curve over Rational number.

        According to Nagell-Lutz theorem, torsion point p(x, y)
        x and y are integers, either y = 0 or y**2 is divisor
        of discriminent. According to Mazur's theorem, there are
        at most 15 points in torsion collection.

        Examples
        ========

        >>> from sympy.ntheory.ec import EllipticCurve
        >>> e2 = EllipticCurve(-43, 166)
        >>> sorted(e2.torsion_points())
        [(-5, -16), (-5, 16), (0, 1, 0), (3, -8), (3, 8), (11, -32), (11, 32)]

        """
        if self.characteristic > 0:
            raise ValueError("No torsion point for Finite Field.")
        l = [self.O]
        for x in solve(self._eq.subs(y, 0)):
            if x.is_rational and self.point_order((x, 0,)) != oo:
                l.append((x, 0,))
        for i in divisors(self.discriminent, generator=True):
            j = int(i**.5)
            if j**2 == i:
                for x in solve(self._eq.subs(y, j)):
                    if x.is_rational and self.point_order((x, j,)) != oo:
                        l.extend([(x, j,), (x, -j,)])
        return l

    @property
    def characteristic(self):
        """
        Return domain characteristic.

        """
        return self._domain.characteristic()

    @property
    def discriminent(self):
        """
        Return curve discriminent.

        """
        return self._discrim

    @property
    def is_singular(self):
        """
        Return True if curve discriminent is equal to zero.
        """
        return self.discriminent == 0

    @property
    def j_invariant(self):
        """
        Return curve j-invariant.

        Examples
        ========

        >>> from sympy.ntheory.ec import EllipticCurve
        >>> e1 = EllipticCurve(-10, -20, 0, -1, 1)
        >>> e2 = e1.minimal()
        >>> e1.j_invariant == e2.j_invariant
        True

        """
        c4 = self._b2**2 - 24*self._b4
        return c4**3 / self.discriminent

    @property
    def order(self):
        """
        Number of points in Finite field.

        Examples
        ========

        >>> from sympy.polys.domains import FF
        >>> from sympy.ntheory.ec import EllipticCurve
        >>> e2 = EllipticCurve(1, 0, domain=FF(19))
        >>> e2.order
        19

        """
        if self.characteristic == 0:
            raise NotImplementedError("Still not implemented")
        return len(list(self.points()))

    @property
    def rank(self):
        """
        Number of independent points of infinite order.

        For Finite field, it must be 0.
        """
        if self._rank is not None:
            return self._rank
        raise NotImplementedError("Still not implemented")
