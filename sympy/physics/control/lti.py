from sympy import Basic, Mul, degree

__all__ = ['TransferFunction',]


class TransferFunction(Mul):
    # this will subclass from `Basic`.
    def __new__(cls, num, den):
        obj = Mul.__new__(cls, num, 1/den)
        obj.num = num
        obj.den = den
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
        pass

    def __sub__(self, other):
        pass

    def __mul__(self, other):
        pass

    def __div__(self, other):
        pass

    def __pow__(self, p):
        pass

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
