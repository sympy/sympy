from sympy import QQ
from sympy.abc import x, y
from sympy.core.relational import Eq
from .residue_ntheory import sqrt_mod
from sympy.polys.domains import FiniteField, RationalField

class EllipticCurve():

    def __init__(self, a4, a6, a1=0, a2=0, a3=0, domain=QQ):
        self._coeff = [a4, a6, a1, a2, a3]
        self._domain = domain
        self._discrim = int(self._domain(-16)*(4*a4**3+27*a6**2))
        self._eq = Eq(y**2 + a1*x*y + a3*y, x**3 + a2*x**2 + a4*x + a6)
        if isinstance(self._domain, FiniteField):
            self._char = self._domain.mod
            self._rank = 0
        elif isinstance(self._domain, RationalField):
            self._char = 0
            self._rank = None
        self._point = None
        self._points = None

    def __add__(self, other):
        p1 = self._point
        if p1[2] == 0:
            return other
        p2 = other
        if p2[2] == 0:
            return p1
        x1 = self._domain(p1[0])
        y1 = self._domain(p1[1])
        x2, y2 = p2[:2]
        if x1 != x2:
            slope = (y1 - y2) / (x1 - x2)
        else:
            if (y1 + y2) == 0:
                return 0, 1, 0
            slope = (3 * x1**2 + self._coeff[0]) / (2 * y1)
        x3 = slope**2 - x1 - x2
        y3 = -y1 - slope * (x3 - x1)
        return x3, y3, 1

    def __call__(self, x, y, z=1):
        self._point = (x, y, z)
        return self

    def __contains__(self, point):
        if self.characteristic == 0 and len(point) == 3 and point[2] == 0:
            return True
        return self._eq.subs({x:point[0], y:point[1]})

    def __mul__(self, other):
        if other < 1:
            return self._point
        r = (0, 1, 0)
        p = self._point
        while other:
            if other & 1:
                self._point = r
                r = self.__add__(p)
            other >>= 1
            self._point = p
            p = self.__add__(p)
        return r

    def __repr__(self):
        return 'E({}): y**2 = x**3 + {}x + {}'.format(self._domain, *self._coeff)

    def points(self):
        if self._points is not None:
            return self._points
        char = self.characteristic
        if char > 1:
            self._points = []
            for i in range(1, char + 1):
                y = sqrt_mod(i**3 + self._coeff[0]*i + self._coeff[1], char)
                if y is not None:
                    self._points.extend([(i, y), (i, char - y)])
        raise NotImplementedError("Still not implemented")

    @property
    def characteristic(self):
        return self._char

    @property
    def discriminent(self):
        return self._discrim

    @property
    def order(self):
        if self.characteristic == 0:
            raise NotImplementedError("Still not implemented")
        return len(self.points())

    @property
    def rank(self):
        if self._rank is not None:
            return self._rank
        raise NotImplementedError("Still not implemented")
