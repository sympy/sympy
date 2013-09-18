from sympy import QQ
from sympy.abc import x, y
from sympy.core.relational import Eq
from .factor_ import divisors
from .residue_ntheory import sqrt_mod
from sympy.polys.domains import FiniteField, RationalField
from sympy.solvers.solvers import solve


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

    def __init__(self, a4, a6, a1=0, a2=0, a3=0, domain=QQ):
        self._domain = domain
        # Calculate discriminant
        b2 = a1**2 + 4 * a2
        b4 = 2 * a4 + a1 * a3
        b6 = a3**2 + 4 * a6
        b8 = a1**2 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3**2 - a4**2
        self._discrim = int(self._domain(-b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6))
        self._eq = Eq(y**2 + a1*x*y + a3*y, x**3 + a2*x**2 + a4*x + a6)
        self._a1 = self._domain(a1)
        self._a2 = self._domain(a2)
        self._a3 = self._domain(a3)
        self._a4 = self._domain(a4)
        self._a6 = self._domain(a6)
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
                return 0, 1, 0
            slope = (3 * x1**2 + self._a4) / (2 * y1)
        x3 = slope**2 - x1 - x2
        y3 = -y1 - slope * (x3 - x1)
        if to_sympy:
            return self._domain.to_sympy(x3), self._domain.to_sympy(y3), 1
        return x3, y3, 1

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
        r = (0, 1, 0)
        while n:
            if n & 1:
                r = self.add(r, p, to_sympy=False)
            n >>= 1
            p = self.add(p, p, to_sympy=False)
        return self._domain.to_sympy(r[0]), self._domain.to_sympy(r[1]), r[2]

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

    def torsion_list(self):
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
        >>> sorted(e2.torsion_list())
        [(-5, -16), (-5, 16), (3, -8), (3, 8), (11, -32), (11, 32)]

        """
        if self.characteristic > 0:
            raise ValueError("No torsion point for Finite Field.")
        l = []
        for x in solve(self._eq.subs(y, 0)):
            if x.is_rational:
                l.append((x, 0,))
        for i in divisors(self.discriminent, generator=True):
            j = int(i**.5)
            if j**2 == i:
                for x in solve(self._eq.subs(y, j)):
                    if x.is_rational:
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
