
from basic import Basic

class ArithMeths:

    def __pos__(self):
        return self

    def __neg__(self):
        return Basic.Mul(-1, self)

    def __add__(self, other):
        return Basic.Add(self, other)

    def __sub__(self, other):
        return self + (-Basic.sympify(other))

    def __mul__(self, other):
        return Basic.Mul(self, other)

    def __div__(self, other):
        return self * (Basic.sympify(other) ** (-1))

    def _eval_power(self, exponent):
        return

    def __pow__(self, other):
        return Basic.Pow(self, other)

class ImmutableMeths:

    def __setitem__(self, k, v):
        raise TypeError('%s instance is immutable' % (self.__class__.__name__))

    def __delitem__(self, k):
        raise TypeError('%s instance is immutable' % (self.__class__.__name__))

    def popitem(self):
        raise TypeError('%s instance is immutable' % (self.__class__.__name__))
