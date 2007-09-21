from utils import memoizer_immutable_args
from basic import Basic, sympify
from number import Rational

pyint = long
makeinteger = lambda p: pyint.__new__(Integer, p)

class Integer(Rational, pyint):

    is_integer = True

    @memoizer_immutable_args('Integer.__new__')
    def __new__(cls, p):
        obj = pyint.__new__(cls, p)
        if p==0:
            obj.is_zero = True            
        elif p==1:
            obj.is_one = True
        elif p==2:
            obj.is_two = True
        return obj

    make = staticmethod(makeinteger)

    @property
    def p(self): return pyint(self)

    @property
    def q(self): return 1

    # relational methods

    def compare(self, other):
        if self is other: return 0
        c = cmp(self.__class__, other.__class__)
        if c: return c
        return pyint.__cmp__(self, pyint(other))

    # converter methods

    __int__ = pyint.__int__
    __long__ = pyint.__long__
    __float__ = pyint.__float__

    def evalf(self):
        return Basic.Float(pyint(self))

    # mathematical properties

    @property
    def is_even(self):
        return pyint.__mod__(self,2)==0

    @property
    def is_odd(self):
        return pyint.__mod__(self,2)==1

    @property
    def is_positive(self):
        return pyint.__cmp__(self, 0L)==1

    @property
    def is_negative(self):
        return pyint.__cmp__(self, 0L)==-1

    # algorithms

    @staticmethod
    def gcd(a, b):
        """
        Returns the Greatest Common Divisor, implementing Euclid's algorithm.
        """
        while a:
            a, b = b%a, a
        return b

    @staticmethod
    def factor_trial_division(n):
        """
        Factor any integer into a product of primes, 0, 1, and -1.
        Returns a dictionary {<prime: exponent>}.
        """
        if not n:
            return {0:1}
        factors = {}
        if n < 0:
            factors[-1] = 1
            n = -n
        if n==1:
            factors[1] = 1
            return factors
        d = 2
        while n % d == 0:
            try:
                factors[d] += 1
            except KeyError:
                factors[d] = 1
            n //= d
        d = 3
        while n > 1 and d*d <= n:
            if n % d:
                d += 2
            else:
                try:
                    factors[d] += 1
                except KeyError:
                    factors[d] = 1
                n //= d
        if n>1:
            try:
                factors[n] += 1
            except KeyError:
                factors[n] = 1
        return factors

    def __pos__(self):
        return self

    def __neg__(self):
        return Integer(pyint.__neg__(self))

    def __add__(self, other):
        other = sympify(other)
        if other.is_Integer:
            return Integer(pyint.__add__(self, other))
        return NotImplemented

    def __sub__(self, other):
        other = sympify(other)
        if other.is_Integer:
            return Integer(pyint.__sub__(self, other))
        return NotImplemented

    def __mul__(self, other):
        other = sympify(other)
        if other.is_Integer:
            return Integer(pyint.__mul__(self, other))
        return NotImplemented

    def __div__(self, other):
        other = sympify(other)
        if other.is_Integer:
            return Basic.Fraction(self.p, other.p)
        return NotImplemented

    def __pow__(self, other):
        other = sympify(other)
        if other.is_Integer:
            if other.is_negative:
                return Basic.Fraction(1, pyint.__pow__(self, -other.p))
            return Integer(pyint.__pow__(self, other))
        return NotImplemented

    # __r*__ methods are defined in methods.py
