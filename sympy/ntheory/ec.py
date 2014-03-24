from sympy.abc import x, y
from sympy.core.compatibility import as_int, is_sequence
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

    def __call__(self, x, y, z=1):
        return EllipticCurvePoint(x, y, z, self)

    def __contains__(self, point):
        if is_sequence(point):
            if len(point) == 2:
                z1 = 1
            else:
                z1 = point[2]
            x1, y1 = point[:2]
        elif isinstance(point, EllipticCurvePoint):
            x1, y1, z1 = point.x, point.y, point.z
        else:
            raise ValueError('Invalid point.')
        if self.characteristic == 0 and z1 == 0:
            return True
        return self._eq.subs({x: x1, y: y1})

    def __repr__(self):
        return 'E({}): {}'.format(self._domain, self._eq)

    def minimal(self):
        """
        Return minimal Weierstrass equation.

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
                    yield self(i, y)
                    if y != 0:
                        yield self(i, char - y)
        else:
            raise NotImplementedError("Still not implemented")

    def torsion_points(self):
        """
        Return torsion points of curve over Rational number.

        Return point objects those are finite order.
        According to Nagell-Lutz theorem, torsion point p(x, y)
        x and y are integers, either y = 0 or y**2 is divisor
        of discriminent. According to Mazur's theorem, there are
        at most 15 points in torsion collection.

        Examples
        ========

        >>> from sympy.ntheory.ec import EllipticCurve
        >>> e2 = EllipticCurve(-43, 166)
        >>> sorted(e2.torsion_points())
        [(-5, -16), (-5, 16), O, (3, -8), (3, 8), (11, -32), (11, 32)]

        """
        if self.characteristic > 0:
            raise ValueError("No torsion point for Finite Field.")
        l = [EllipticCurvePoint.point_at_infinity(self)]
        for x in solve(self._eq.subs(y, 0)):
            if x.is_rational:
                l.append(self(x, 0))
        for i in divisors(self.discriminant, generator=True):
            j = int(i**.5)
            if j**2 == i:
                for x in solve(self._eq.subs(y, j)):
                    p = self(x, j)
                    if x.is_rational and p.order() != oo:
                        l.extend([p, -p])
        return l

    @property
    def characteristic(self):
        """
        Return domain characteristic.

        Examples
        ========

        >>> from sympy.ntheory.ec import EllipticCurve
        >>> e2 = EllipticCurve(-43, 166)
        >>> e2.characteristic
        0

        """
        return self._domain.characteristic()

    @property
    def discriminant(self):
        """
        Return curve discriminant.

        Examples
        ========

        >>> from sympy.ntheory.ec import EllipticCurve
        >>> e2 = EllipticCurve(0, 17)
        >>> e2.discriminant
        -124848

        """
        return int(self._discrim)

    @property
    def is_singular(self):
        """
        Return True if curve discriminant is equal to zero.
        """
        return self.discriminant == 0

    @property
    def j_invariant(self):
        """
        Return curve j-invariant.

        Examples
        ========

        >>> from sympy.ntheory.ec import EllipticCurve
        >>> e1 = EllipticCurve(-2, 0, 0, 1, 1)
        >>> e1.j_invariant
        1404928/389

        """
        c4 = self._b2**2 - 24*self._b4
        return self._domain.to_sympy(c4**3 / self._discrim)

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

EC = EllipticCurve

class EllipticCurvePoint():
    """
    Point of Elliptic Curve

    Examples
    ========

    >>> from sympy.ntheory.ec import EllipticCurve
    >>> e1 = EllipticCurve(-17, 16)
    >>> p1 = e1(0, -4, 1)
    >>> p2 = e1(1, 0)
    >>> p1 + p2
    (15, -56)
    >>> e3 = EllipticCurve(-1, 9)
    >>> e3(1, -3) * 3
    (664/169, 17811/2197)
    >>> e2 = EC(-2, 0, 0, 1, 1)
    >>> p = e2(-1,1)
    >>> q = e2(0, -1)
    >>> p+q
    (4, 8)
    >>> p-q
    (1, 0)
    >>> 3*p-5*q
    (328/361, -2800/6859)
    """

    @staticmethod
    def point_at_infinity(curve):
        return EllipticCurvePoint(0, 1, 0, curve)

    def __init__(self, x, y, z, curve):
        self.x = x
        self.y = y
        self.z = z
        self._curve = curve

    def __add__(self, p):
        if self.z == 0:
            return p
        if p.z == 0:
            return self
        x1, y1 = self.x, self.y
        x2, y2 = p.x, p.y
        a1 = self._curve._a1
        a2 = self._curve._a2
        a3 = self._curve._a3
        a4 = self._curve._a4
        a6 = self._curve._a6
        if x1 != x2:
            slope = (y1 - y2) / (x1 - x2)
            yint = (y1 * x2 - y2 * x1) / (x2 - x1)
        else:
            if (y1 + y2) == 0:
                return self.point_at_infinity(self._curve)
            slope = (3 * x1**2 + 2*a2*x1 + a4 - a1*y1) / (a1 * x1 + a3 + 2 * y1)
            yint = (-x1**3 + a4*x1 + 2*a6 - a3*y1) / (a1*x1 + a3 + 2*y1)
        x3 = slope**2 + a1*slope - a2 - x1 - x2
        y3 = -(slope + a1) * x3 - yint - a3
        return EllipticCurvePoint(x3, y3, 1, self._curve)

    def __lt__(self, other):
        return (self.x, self.y, self.z) < (other.x, other.y, other.z)

    def __mul__(self, n):
        n = as_int(n)
        r = self.point_at_infinity(self._curve)
        if n == 0:
            return r
        if n < 0:
            return -self * -n
        p = self
        while n:
            if n & 1:
                r = r + p
            n >>= 1
            p = p + p
        return r

    def __rmul__(self, n):
        return self * n

    def __neg__(self):
        return EllipticCurvePoint(self.x, -self.y - self._curve._a1*self.x - self._curve._a3, self.z, self._curve)

    def __repr__(self):
        if self.z == 0:
            return 'O'
        dom = self._curve._domain
        try:
            return '({}, {})'.format(dom.to_sympy(self.x), dom.to_sympy(self.y))
        except TypeError:
            pass
        return '({}, {})'.format(self.x, self.y)

    def __sub__(self, other):
        return self + -other

    def order(self):
        """
        Return point order n where nP = 0.

        """
        if self.z == 0:
            return 1
        if self.y == 0:  # P = -P
            return 2
        p = self * 2
        if p.y == -self.y:  # 2P = -P
            return 3
        i = 2
        while int(p.x) == p.x:
            p = self + p
            i += 1
            if p.z == 0:
                return i
        return oo
