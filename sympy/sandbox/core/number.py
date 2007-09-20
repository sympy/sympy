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

    def __new__(cls, f):
        return Float(f)

class Rational(Real):

    is_rational = True

    def __new__(cls, p, q=1):
        return Fraction(p, q)

    def __eq__(self, other):
        other = sympify(other)
        if self is other: return True
        if other.is_Number:
            if other.is_Rational:
                return self.compare(other)==0
            if other.is_Float:
                return self.evalf()==other
        return super(Real, self).__eq__(other)



from py_integer import Integer
from py_fraction import Fraction
from py_float import Float
