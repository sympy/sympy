
from basic import Basic

class RelationalMeths:

    def __eq__(self, other):
        if self is other: return True
        return Basic.Equality(self, other)

    def __ne__(self, other):
        if self is other: return False
        return Basic.Unequality(self, other)

    def __lt__(self, other):
        return Basic.StrictInequality(self, other)

    def __gt__(self, other):
        return Basic.StrictInequality(other, self)

    def __le__(self, other):
        return Basic.Inequality(self, other)

    def __ge__(self, other):
        return Basic.Inequality(other, self)


class ArithMeths:

    def __pos__(self):
        return self

    def __neg__(self):
        return Basic.Mul(-1, self)

    def __add__(self, other):
        return Basic.Add(self, other)

    __radd__ = __add__

    def __sub__(self, other):
        return self + (-Basic.sympify(other))

    def __mul__(self, other):
        return Basic.Mul(self, other)

    __rmul__ = __mul__

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
