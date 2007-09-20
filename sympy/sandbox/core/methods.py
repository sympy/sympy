
from basic import Basic, sympify

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
        return self * (-1)

    def __add__(self, other):
        return Basic.Add(self, other)

    def __radd__(self, other):
        return sympify(other) + self

    def __sub__(self, other):
        return self + (-sympify(other))

    def __rsub__(self, other):
        return sympify(other) - self

    def __mul__(self, other):
        return Basic.Mul(self, other)

    def __rmul__(self, other):
        return sympify(other) * self

    def __div__(self, other):
        return self * (sympify(other) ** (-1))

    def __rdiv__(self, other):
        return sympify(other) / self

    def __pow__(self, other):
        return Basic.Pow(self, other)

    def __rpow__(self, other):
        return sympify(other) ** self

    def _eval_power(self, exponent):
        return



class ImmutableMeths:

    def __setitem__(self, k, v):
        raise TypeError('%s instance is immutable' % (self.__class__.__name__))

    def __delitem__(self, k):
        raise TypeError('%s instance is immutable' % (self.__class__.__name__))

    def popitem(self):
        raise TypeError('%s instance is immutable' % (self.__class__.__name__))

class NumberMeths(ArithMeths, RelationalMeths):

    def __eq__(self, other):
        raise NotImplementedError('%s must implement __eq__ method' % (self.__class__.__name__))

    def __add__(self, other):
        raise NotImplementedError('%s must implement __add__ method' % (self.__class__.__name__))

    def __mul__(self, other):
        raise NotImplementedError('%s must implement __mul__ method' % (self.__class__.__name__))

    def __div__(self, other):
        raise NotImplementedError('%s must implement __div__ method' % (self.__class__.__name__))

    def __pow__(self, other):
        raise NotImplementedError('%s must implement __pow__ method' % (self.__class__.__name__))
