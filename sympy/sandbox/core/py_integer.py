from utils import memoizer_immutable_args
from basic import Basic, sympify
from number import Rational


class Integer(Rational, int):

    is_integer = True
    is_even = None
    is_odd = None

    @memoizer_immutable_args('Integer.__new__')
    def __new__(cls, p):
        obj = int.__new__(cls, p)
        if p==0:
            obj.is_zero = True            
        elif p==1:
            obj.is_one = True
        elif p==2:
            obj.is_two = True
        return obj

    @property
    def p(self): return int(self)

    @property
    def q(self): return 1

    def torepr(self):
        return '%s(%r)' % (self.__class__.__name__, self.p)

    # relational methods

    def compare(self, other):
        if self is other: return 0
        c = cmp(self.__class__, other.__class__)
        if c: return c
        #
        return int.__cmp__(self, int(other))

    # converter methods

    # int has __int__, __float__

    def evalf(self):
        return Basic.Float(int(self))

    # arithmetic methods

    def __add__(self, other):
        other = sympify(other)
        if other.is_Number:
            if other.is_Integer:
                return Integer(int.__add__(self, other))
            elif other.is_Fraction:
                return other + self
            elif other.is_Float:
                return self.evalf() + other
            else:
                raise NotImplementedError('%s.__add__(%s)' \
                                          % (self.__class__.__name__,
                                             other.__class__.__name__))
        return Basic.Add(self, other)

    def __mul__(self, other):
        other = sympify(other)
        if other.is_Number:
            if other.is_Integer:
                return Integer(int.__mul__(self, other))
            elif other.is_Fraction:
                return other * self
            elif other.is_Float:
                return self.evalf() * other
            else:
                raise NotImplementedError('%s.__mul__(%s)' \
                                          % (self.__class__.__name__,
                                             other.__class__.__name__))
        return Basic.Mul(self, other)

    def __div__(self, other):
        other = sympify(other)
        if other.is_Number:
            if other.is_Integer or other.is_Fraction:
                return Basic.Fraction(self, other)
            elif other.is_Float:
                return self.evalf() / other
            else:
                raise NotImplementedError('%s.__div__(%s)' \
                                          % (self.__class__.__name__,
                                             other.__class__.__name__))
        return Basic.Mul(self, other)

    def __pow__(self, other):
        other = sympify(other)
        if other.is_Number:
            if other.is_Integer:
                if other.is_nonnegative:
                    return Integer(int.__pow__(self, other))
                else:
                    return Basic.Fraction(1, int.__pow__(self, -other))
            elif other.is_Fraction:
                pass
            elif other.is_Float:
                return self.evalf() ** other
            else:
                raise NotImplementedError('%s.__pow__(%s)' \
                                          % (self.__class__.__name__,
                                             other.__class__.__name__))
        return Basic.Pow(self, other)

    # mathematical properties

    @property
    def is_even(self):
        return int.__mod__(self,2)==0

    @property
    def is_odd(self):
        return int.__mod__(self,2)==1

    @property
    def is_positive(self):
        return int(self)>0

    @property
    def is_negative(self):
        return int(self)<0
