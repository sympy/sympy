from utils import memoizer_immutable_args
from basic import Basic, sympify
from number import Rational

def makefraction(p,q):
    return tuple.__new__(Fraction, (p, q))

def makefraction_from_man_exp(man, exp):
    if exp > 0:
        return makefraction(man * 2 ** exp, 1)
    obj = Fraction(man, 2** -exp)
    if obj.is_Fraction: return obj
    return obj.as_Fraction

class Fraction(Rational, tuple):
    """
    Represents a ratio p/q of two integers.

    This implementation relies completely on the implementation
    of Integer.
    """

    @memoizer_immutable_args('Fraction.__new__')
    def __new__(cls, p, q):
        if q<0:
            p, q = -p, -q
        r = Basic.Integer.gcd(abs(p), q)
        if r>1:
            p //= r
            q //= r
        if q==1:
            return Basic.Integer(p)
        return tuple.__new__(cls, (p, q))

    make = staticmethod(makefraction)
    make_from_man_exp = staticmethod(makefraction_from_man_exp)

    @property
    def p(self): return self[0]
    
    @property
    def q(self): return self[1]

    def __pos__(self):
        return self

    def __neg__(self):
        return self.make(-self.p, self.q)

    def __add__(self, other):
        other = sympify(other)
        if other.is_Integer:
            other = other.as_Fraction
        if other.is_Fraction:
            return Fraction(self.p * other.q + self.q * other.p,
                            self.q * other.q)
        return NotImplemented

    def __sub__(self, other):
        other = sympify(other)
        if other.is_Integer:
            other = other.as_Fraction
        if other.is_Fraction:
            return Fraction(self.p * other.q - self.q * other.p,
                            self.q * other.q)
        return NotImplemented

    def __mul__(self, other):
        other = sympify(other)
        if other.is_Integer:
            other = other.as_Fraction
        if other.is_Fraction:
            return Fraction(self.p * other.p, self.q * other.q)
        return NotImplemented

    def __div__(self, other):
        other = sympify(other)
        if other.is_Integer:
            other = other.as_Fraction
        if other.is_Fraction:
            return Fraction(self.p * other.q, self.q * other.p)
        return NotImplemented

    def __pow__(self, other):
        other = sympify(other)
        if other.is_Integer:
            if other.is_negative:
                p = -other.p
                return Fraction.make(self.q ** p, self.p ** p)
            p = other.p
            return Fraction.make(self.p ** p, self.q ** p)
        if other.is_Fraction:
            return Basic.Pow(self, other)
        return NotImplemented

    def __radd__(self, other):
        if isinstance(other, Basic):
            if other.is_Integer:
                other = other.as_Fraction
            if other.is_Fraction:
                return Fraction(self.p * other.q + self.q * other.p,
                                self.q * other.q)
            return Basic.Add(other, self)
        return sympify(other) + self

    def __rsub__(self, other):
        if isinstance(other, Basic):
            if other.is_Integer:
                other = other.as_Fraction
            if other.is_Fraction:
                return Fraction(-self.p * other.q + self.q * other.p,
                                self.q * other.q)
            return Basic.Add(other, -self)
        return sympify(other) - self

    def __rmul__(self, other):
        if isinstance(other, Basic):
            if other.is_Integer:
                other = other.as_Fraction
            if other.is_Fraction:
                return Fraction(self.p * other.p, self.q * other.q)
            return Basic.Mul(other, self)
        return sympify(other) * self

    def __rdiv__(self, other):
        if isinstance(other, Basic):
            if other.is_Integer:
                other = other.as_Fraction
            if other.is_Fraction:
                return Fraction(self.q * other.p, self.p * other.q)
            return Basic.Mul(other, 1/self)
        return sympify(other) / self

    def __rpow__(self, other):
        if isinstance(other, Basic):
            return Basic.Pow(other, self)
        return sympify(other) ** self
