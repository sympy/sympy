from utils import memoizer_immutable_args
from basic import Basic, sympify
from number import Rational

class Fraction(Rational, tuple):
    """
    Represents a ratio p/q of two integers.

    This implementation relies completely on the implementation
    of Integer.
    """

    @memoizer_immutable_args('Fraction.__new__')
    def __new__(cls, p, q):
        if q<0:
            p, q = -p, -q
        r = Basic.Integer.gcd(abs(p), q)
        if r>1:
            p //= r
            q //= r
        if q==1:
            return Basic.Integer(p)
        return tuple.__new__(cls, (p, q))

    @property
    def p(self): return self[0]
    
    @property
    def q(self): return self[1]

