from utils import memoizer_immutable_args
from basic import Basic, sympify
from number import Real

class Float(Real, float):

    @memoizer_immutable_args('Float.__new__')
    def __new__(cls, f):
        if isinstance(f, Basic):
            return f.evalf()
        return float.__new__(cls, f)

    def torepr(self):
        return '%s(%r)' % (self.__class__.__name__, float(self))

    def evalf(self):
        return self

    def compare(self, other):
        if self is other: return 0
        c = cmp(self.__class__, other.__class__)
        if c: return c
        # float does not have __cmp__ method.
        s,o = float(self), float(other)
        if s==o: return 0
        if s<o: return -1
        return 1

    def __eq__(self, other):
        other = sympify(other)
        if self is other: return True
        if other.is_Number:
            return self.compare(other.evalf())==0
        return super(Real, self).__eq__(other)

    # float has __int__, __float__

    def __add__(self, other):
        other = sympify(other)
        if other.is_Number:
            if other.is_Float:
                return Float(float.__add__(self, float(other)))
            elif other.is_Integer or other.is_Fraction:
                return Float(float.__add__(self, float(other)))
            else:
                raise NotImplementedError('%s.__add__(%s)' \
                                          % (self.__class__.__name__,
                                             other.__class__.__name__))
        return Basic.Add(self, other)

    def __mul__(self, other):
        other = sympify(other)
        if other.is_Number:
            if other.is_Float:
                return Float(float.__mul__(self, other))
            elif other.is_Integer or other.is_Fraction:
                return Float(float.__mul__(self, float(other)))
            else:
                raise NotImplementedError('%s.__mul__(%s)' \
                                          % (self.__class__.__name__,
                                             other.__class__.__name__))
        return Basic.Mul(self, other)

    def __div__(self, other):
        other = sympify(other)
        if other.is_Number:
            if other.is_Float:
                return Float(float.__div__(self, other))
            elif other.is_Integer or other.is_Fraction:
                return Float(float.__div__(self, float(other)))
            else:
                raise NotImplementedError('%s.__div__(%s)' \
                                          % (self.__class__.__name__,
                                             other.__class__.__name__))
        return super(Real, self).__div__(other)

    def __pow__(self, other):
        other = sympify(other)
        if other.is_Number:
            if other.is_Float:
                return Float(float.__pow__(self, other))
            elif other.is_Integer or other.is_Fraction:
                return Float(float.__pow__(self, float(other)))
            else:
                raise NotImplementedError('%s.__pow__(%s)' \
                                          % (self.__class__.__name__,
                                             other.__class__.__name__))
        return Basic.Pow(self, other)
