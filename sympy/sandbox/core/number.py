from utils import memoizer_immutable_args
from basic import Basic, Atom, sympify
from methods import NumberMeths

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

class Number(NumberMeths, Atom):
    """A Number is an atomic object with a definite numerical value.
    Examples include rational numbers (-25, 2/3, ...) via the Rational
    class, floating-point numbers via the Real class, the imaginary
    unit I, and some special constants like pi."""

    is_negative = None
    is_positive = None
    is_finite = True
    is_bounded = True
    is_commutative = True
    
    def __new__(cls, x, **options):
        if isinstance(x, Basic): return x
        if isinstance(x, (int, long)): return Integer(x, **options)
        return Basic.__new__(cls, x, **options)

    @property
    def as_Interval(self):
        if self.is_Interval: return self
        return Interval.make(self, self)

    @property
    def as_Float(self):
        if self.is_Float: return self
        if self.is_Interval: return self.mid.as_Float
        if self.is_Rational: return Float.make_from_fraction(self.p, self.q)
        raise NotImplementedError(`self`)

    @property
    def as_Fraction(self):
        if self.is_Fraction: return self
        if self.is_Integer: return Fraction.make(self.p, self.q)
        if self.is_Float: return Fraction.make_from_man_exp(self.man, self.exp)
        if self.is_Interval: return self.mid.as_Fraction
        raise NotImplementedError(`self`)

    @property
    def as_Integer(self):
        if self.is_Integer: return self
        raise NotImplementedError(`self`)

    def as_native(self):
        """
        Return internal representation of number implementation.
        """
        raise NotImplementedError(`self`)

    def compare(self, other):
        if self is other: return 0
        c = cmp(self.__class__, other.__class__)
        if c: return c
        return cmp(self.as_native(), other.as_native())

class Real(Number):

    """
    The implementation of a floating point number must
    derive from Real class and have the methods defined
    in the following template class:

    class Float(Real, <floating point implementation class>):

        def __new__(cls, f):
            if isinstance(f, Basic):
                return f.evalf()
            # f can be python int/long, float, str and
            # any other object that the implementation
            # can handle.
            obj = <floating point implementation class>.__new__(cls,..)
            return obj

        def as_native(self):
            return <floating point implementation instance>

        def __float__(self):
            return <python float instance>

        def __int__(self):
            return <python int instance>
    """
    is_real = True

    def __new__(cls, f):
        return Float(f)

    def torepr(self):
        return '%s(%r)' % (self.__class__.__name__, self.as_native())

    @property
    def is_nonpositive(self):
        return self.is_negative or self.is_zero

    @property
    def is_nonnegative(self):
        return self.is_positive or self.is_zero

    def __eq__(self, other):
        other = sympify(other)
        if self is other: return True
        if other.is_Number:
            return self.compare(other.evalf())==0
        return super(Number, self).__eq__(other)


class Rational(Real):

    """ Base class for Integer and Fraction.
    
    Rational subclasses must define attributes p and q that
    hold integer instances of numerator and denominator.

    The type of p and q must support arithmetics and comparsions
    with python integer types and they should return instances
    of integer implementation classes.

    Integer must define staticmethod gcd(a,b) and probably
    also redefine __int__,__long__,__float__ methods.
    """

    is_rational = True

    def __new__(cls, p, q=1):
        return Fraction(p, q)

    def torepr(self):
        if self.q==1:
            return '%s(%r)' % (self.__class__.__name__, self.p)
        return '%s(%r, %r)' % (self.__class__.__name__, self.p, self.q)

    @property
    def is_positive(self):
        return self.p > 0

    @property
    def is_negative(self):
        return self.p < 0

    def compare(self, other):
        if self is other: return 0
        c = cmp(self.__class__, other.__class__)
        if c: return c
        return cmp(self.p*other.q, self.q*other.p)

    def evalf(self):
        return Basic.Float(self.p) / Basic.Float(self.q)

    def __int__(self):
        return int(self.p // self.q)

    def __long__(self):
        return long(self.p // self.q)

    def __float__(self):
        return float(self.evalf())

    def __eq__(self, other):
        other = sympify(other)
        if self is other: return True
        if other.is_Number:
            if other.is_Float:
                return self.evalf()==other
            return self.compare(other)==0
        return super(Real, self).__eq__(other)

from py_integer import Integer
from py_fraction import Fraction
#from py_float import Float
from numerics_float import Float
#from decimal_float import Float
from interval import Interval
