from functools import reduce
import random
from sympy import QQ
from sympy.abc import x, y
from sympy.core.numbers import igcdex, ilcm
from sympy.core.relational import Eq
from sympy.core.sympify import sympify
from sympy.polys.domains import FiniteField, RationalField

class EllipticCurve():

    def __init__(self, a4, a6, a1=0, a2=0, a3=0, domain=QQ):
        self._coeff = [a4, a6, a1, a2, a3]
        self._domain = domain
        self._discrim = int(self._domain(-16)*(4*a4**3+27*a6**2))
        self._eq = Eq(y**2 + a1*x*y + a3*y, x**3 + a2*x**2 + a4*x + a6)
        if isinstance(self._domain, FiniteField):
            self._char = self._domain.mod
        elif isinstance(self._domain, RationalField):
            self._char = 0

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
        if len(point) == 3 and point[2] == 0:
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
        return self._eq.__repr__()

    @property
    def characteristic(self):
        return self._char

    @property
    def discriminent(self):
        return self._discrim
