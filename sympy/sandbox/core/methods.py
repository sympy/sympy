
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
        return Basic.Mul(Basic.Integer(-1), self)

    def __add__(self, other):
        return Basic.Add(self, other)

    def __sub__(self, other):
        return Basic.Add(self, (-sympify(other)))

    def __mul__(self, other):
        return Basic.Mul(self, other)

    def __div__(self, other):
        return Basic.Mul(self, (sympify(other) ** (-1)))

    def __pow__(self, other):
        return Basic.Pow(self, other)

    def _eval_power(self, exponent):
        return

    def __radd__(self, other):
        if isinstance(other, Basic):
            return Basic.Add(other, self)
        return sympify(other) + self

    def __rsub__(self, other):
        if isinstance(other, Basic):
            return Basic.Add(other, -self)            
        return sympify(other) - self

    def __rmul__(self, other):
        if isinstance(other, Basic):
            return Basic.Mul(other, self)            
        return sympify(other) * self

    def __rdiv__(self, other):
        if isinstance(other, Basic):
            return Basic.Mul(other, self ** (-1))
        return sympify(other) / self

    def __rpow__(self, other):
        if isinstance(other, Basic):
            return Basic.Pow(other, self)            
        return sympify(other) ** self


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
