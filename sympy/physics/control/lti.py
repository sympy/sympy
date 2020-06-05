from sympy import Basic, Mul, degree, Symbol, sympify, expand
from sympy.core.numbers import Integer

__all__ = ['TransferFunction',]


class TransferFunction(Basic):

    def __new__(cls, num, den, var):

        if not isinstance(var, Symbol):
            raise TypeError("Last argument must be a complex Symbol.")
        if den == 0:
            raise ValueError("TransferFunction can't have a zero denominator.")

        num, den = sympify(num), sympify(den)
        obj = Basic.__new__(cls, num, den, var)
        obj._num = num
        obj._den = den
        obj._var = var

        return obj

    @property
    def num(self):
        return self._num

    @property
    def den(self):
        return self._den

    @property
    def var(self):
        return self._var

    def __add__(self, other):
        if not self.var == other.var:
            raise ValueError("Last argument for both Transfer Functions must be same.")
        p = self.num * other.den + other.num * self.den
        q = self.den * other.den
        return TransferFunction(expand(p), expand(q), self.var)

    def __sub__(self, other):
        if not self.var == other.var:
            raise ValueError("Last argument for both Transfer Functions must be same.")
        p = self.num * other.den - other.num * self.den
        q = self.den * other.den
        return TransferFunction(expand(p), expand(q), self.var)

    def __mul__(self, other):
        if not self.var == other.var:
            raise ValueError("Last argument for both Transfer Functions must be same.")
        p = self.num * other.num
        q = self.den * other.den
        return TransferFunction(expand(p), expand(q), self.var)

    def __div__(self, other):
        if not self.var == other.var:
            raise ValueError("Last argument for both Transfer Functions must be same.")
        p = self.num * other.den
        q = self.den * other.num
        return TransferFunction(expand(p), expand(q), self.var)

    __truediv__ = __div__

    def __pow__(self, p):
        p = sympify(p)
        if not isinstance(p, Integer):
            raise ValueError("Exponent must be an Integer.")
        if p == 0:
            return TransferFunction(1, 1, self.var)
        if p < 0:
            p = abs(p)
            num_, den_ = expand(self.den**p), expand(self.num**p)
            return TransferFunction(num_, den_, self.var)
        if p > 0:
            num_, den_ = expand(self.num**p), expand(self.den**p)
            return TransferFunction(num_, den_, self.var)

    def __neg__(self):
        return TransferFunction(-self.num, self.den, self.var)

    @property
    def is_proper(self):
        return degree(self.num) <= degree(self.den)

    @property
    def is_strictly_proper(self):
        return degree(self.num) < degree(self.den)

    @property
    def is_biproper(self):
        return degree(self.num) == degree(self.den)
