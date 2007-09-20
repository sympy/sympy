from utils import memoizer_immutable_args
from basic import Basic, sympify
from number import Rational

class Fraction(Rational, tuple):

    @memoizer_immutable_args('Fraction.__new__')
    def __new__(cls, p, q):
        p, q = sympify(p), sympify(q)
        if q==1:
            return p
        if p.is_Integer and q.is_Integer:
            return tuple.__new__(cls, (p,q))
        if q.is_Fraction:
            iq = Fraction(q.q, q.p)
            return p * iq
        raise TypeError(`p,q`)
        return p/q

    @property
    def p(self): return self[0]
    
    @property
    def q(self): return self[1]

    def torepr(self):
        return '%s(%r, %r)' % (self.__class__.__name__, self.p, self.q)

    def evalf(self):
        return self.p.evalf() / self.q.evalf()

    def compare(self, other):
        if self is other: return 0
        c = cmp(self.__class__, other.__class__)
        if c: return c
        #
        return cmp(self.p*other.q, self.q*other.p)

    @property
    def is_positive(self):
        return self.p > 0

    @property
    def is_negative(self):
        return self.p < 0

    def __int__(self):
        raise NotImplementedError

    def __float__(self):
        return float(self.p) / float(self.q)

    def __add__(self, other):
        other = Basic.sympify(other)
        if other.is_Number:
            if other.is_Integer:
                return Fraction(self.p + self.q * other, self.q)
            elif other.is_Fraction:
                return Fraction(self.p * other.q + self.q * other.p,
                                self.q * other.q)
            elif other.is_Float:
                return self.evalf() + other
            else:
                raise NotImplementedError('%s.__add__(%s)' \
                                          % (self.__class__.__name__,
                                             other.__class__.__name__))
        return Basic.Add(self, other)

    def __mul__(self, other):
        other = Basic.sympify(other)
        if other.is_Number:
            if other.is_Integer:
                return Fraction(self.p * other, self.q)
            elif other.is_Fraction:
                return Fraction(self.p * other.p, self.q * other.q)
            elif other.is_Float:
                return self.evalf() * other
            else:                  
                raise NotImplementedError('%s.__mul__(%s)' \
                                          % (self.__class__.__name__,
                                             other.__class__.__name__))
        return Basic.Mul(self, other)

    def __div__(self, other):
        other = Basic.sympify(other)
        if other.is_Number:
            if other.is_Integer:
                return Fraction(self.p, self.q * other)
            elif other.is_Fraction:
                return Fraction(self.p * other.q, self.q * other.p)
            elif other.is_Float:
                return self.evalf() / other
            else:                  
                raise NotImplementedError('%s.__div__(%s)' \
                                          % (self.__class__.__name__,
                                             other.__class__.__name__))
        return super(Rational, self).__div__(other)

    def __pow__(self, other):
        other = Basic.sympify(other)
        if other.is_Number:
            if other.is_Integer:
                return Fraction(self.p ** other, self.q ** other)
            elif other.is_Fraction:
                pass
            elif other.is_Float:
                return self.evalf() ** other
            else:             
                raise NotImplementedError('%s.__pow__(%s)' \
                                          % (self.__class__.__name__,
                                             other.__class__.__name__))
        return Basic.Pow(self, other)
