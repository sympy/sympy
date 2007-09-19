from utils import memoizer_immutable_args
from basic import Basic, Atom
from methods import ArithMeths, RelationalMeths

Basic.is_zero = None
Basic.is_one = None
Basic.is_half = None
Basic.is_two = None
Basic.is_even = None
Basic.is_odd = None
Basic.is_negative = None
Basic.is_positive = None
Basic.is_nonnegative = None
Basic.is_nonpositive = None
Basic.is_real = None
Basic.is_integer = None
Basic.is_rational = None
Basic.is_finite = None
Basic.is_bounded = None
Basic.is_commutative = None
Basic.is_prime = None

class Number(ArithMeths, RelationalMeths, Atom):
    """A Number is an atomic object with a definite numerical value.
    Examples include rational numbers (-25, 2/3, ...) via the Rational
    class, floating-point numbers via the Real class, the imaginary
    unit I, and some special constants like pi."""

    is_negative = None
    is_positive = None
    is_real = True
    is_finite = True
    is_bounded = True
    is_commutative = True
    
    def __new__(cls, x, **options):
        if isinstance(x, Basic): return x
        if isinstance(x, (int, long)): return Integer(x, **options)
        return Basic.__new__(cls, x, **options)

    @property
    def is_nonpositive(self):
        return self.is_negative or self.is_zero

    @property
    def is_nonnegative(self):
        return self.is_positive or self.is_zero


class Real(Number):

    pass

class Rational(Number):

    is_rational = True

    def __new__(cls, p, q=1):
        return Fraction(p, q)

class Fraction(Rational, tuple):

    @memoizer_immutable_args('Fraction.__new__')
    def __new__(cls, p, q):
        if q==1:
            return Integer(p)
        assert isinstance(p, int),`p`
        assert isinstance(q, int),`q`
        return tuple.__new__(cls, (p,q))

    @property
    def p(self): return self[0]
    
    @property
    def q(self): return self[1]

    def torepr(self):
        return '%s(%r, %r)' % (self.__class__.__name__, self.p, self.q)

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

    def compare(self, other):
        if self is other: return 0
        c = cmp(self.__class__, other.__class__)
        if c: return c
        #
        return int.__cmp__(self, int(other))

    def __neg__(self):
        return Integer(-int(self))

    def __add__(self, other):
        if isinstance(other, int):
            return Integer(int(self) + int(other))
        return Basic.Add(self, other)

    def __mul__(self, other):
        if isinstance(other, int):
            return Integer(int(self) * int(other))
        return Basic.Mul(self, other)

    def __pow__(self, other):
        if isinstance(other, int) and other>0:
            return Integer(int(self) ** int(other))
        return Basic.Pow(self, other)

    __iadd__ = __add__

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
