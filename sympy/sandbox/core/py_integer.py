from utils import memoizer_immutable_args
from basic import Basic, sympify
from number import Rational


class Integer(Rational, long):

    is_integer = True

    @memoizer_immutable_args('Integer.__new__')
    def __new__(cls, p):
        obj = long.__new__(cls, p)
        if p==0:
            obj.is_zero = True            
        elif p==1:
            obj.is_one = True
        elif p==2:
            obj.is_two = True
        return obj

    @property
    def p(self): return long(self)

    @property
    def q(self): return 1

    # relational methods

    def compare(self, other):
        if self is other: return 0
        c = cmp(self.__class__, other.__class__)
        if c: return c
        return long.__cmp__(self, long(other))

    # converter methods

    __int__ = long.__int__
    __long__ = long.__long__
    __float__ = long.__float__

    def evalf(self):
        return Basic.Float(long(self))

    # mathematical properties

    @property
    def is_even(self):
        return long.__mod__(self,2)==0

    @property
    def is_odd(self):
        return long.__mod__(self,2)==1

    @property
    def is_positive(self):
        return long(self)>0

    @property
    def is_negative(self):
        return long(self)<0

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
