"""Tools for manipulating of large commutative expressions. """

from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.core.numbers import Rational
from sympy.core.singleton import S
from sympy.core.coreerrors import NonCommutativeExpression

def decompose_power(expr):
    """
    Decompose power into symbolic base and integer exponent.

    Example
    =======

    >>> from sympy.core.exprtools import decompose_power
    >>> from sympy.abc import x, y

    >>> decompose_power(x)
    (x, 1)
    >>> decompose_power(x**2)
    (x, 2)
    >>> decompose_power(x**(2*y))
    (x**y, 2)
    >>> decompose_power(x**(2*y/3))
    (x**(y/3), 2)

    """
    base, exp = expr.as_base_exp()

    if exp.is_Number:
        if exp.is_Rational:
            if not exp.is_Integer:
                base = Pow(base, Rational(1, exp.q))

            exp = exp.p
        else:
            base, exp = expr, 1
    else:
        exp, tail = exp.as_coeff_mul()

        if exp.is_Rational:
            if not exp.is_Integer:
                tail += (Rational(1, exp.q),)

            base, exp = Pow(base, Mul(*tail)), exp.p
        else:
            base, exp = expr, 1

    return base, exp

class Factors(object):
    """Efficient representation of ``f_1*f_2*...*f_n``. """

    __slots__ = ['factors', 'gens']

    def __init__(self, factors=None):
        if factors is None:
            factors = {}

        self.factors = factors
        self.gens = frozenset(factors.keys())

    def __repr__(self):
        return "Factors(%s)" % self.factors

    def as_expr(self):
        return Mul(*[ factor**exp for factor, exp in self.factors.iteritems() ])

    def normal(self, other):
        self_factors = dict(self.factors)
        other_factors = dict(other.factors)

        for factor, self_exp in self.factors.iteritems():
            try:
                other_exp = other.factors[factor]
            except KeyError:
                continue

            exp = self_exp - other_exp

            if not exp:
                del self_factors[factor]
                del other_factors[factor]
            else:
                if exp > 0:
                    self_factors[factor] = exp
                    del other_factors[factor]
                else:
                    del self_factors[factor]
                    other_factors[factor] = -exp

        return Factors(self_factors), Factors(other_factors)

    def mul(self, other):
        factors = dict(self.factors)

        for factor, exp in other.factors.iteritems():
            if factor in factors:
                exp = factors[factor] + exp

                if not exp:
                    del factors[factor]
                    continue

            factors[factor] = exp

        return Factors(factors)

    def div(self, other):
        quo, rem = dict(self.factors), {}

        for factor, exp in other.factors.iteritems():
            if factor in quo:
                exp = quo[factor] - exp

                if exp <= 0:
                    del quo[factor]

                if exp >= 0:
                    if exp:
                        quo[factor] = exp

                    continue

                exp = -exp

            rem[factor] = exp

        return Factors(quo), Factors(rem)

    def quo(self, other):
        return self.div(other)[0]

    def rem(self, other):
        return self.div(other)[1]

    def pow(self, other):
        if type(other) is int and other >= 0:
            factors = {}

            if other:
                for factor, exp in self.factors.iteritems():
                    factors[factor] = exp*other

            return Factors(factors)
        else:
            raise ValueError("expected non-negative integer, got %s" % other)

    def gcd(self, other):
        factors = {}

        for factor, exp in self.factors.iteritems():
            if factor in other.factors:
                exp = min(exp, other.factors[factor])
                factors[factor] = exp

        return Factors(factors)

    def lcm(self, other):
        factors = dict(self.factors)

        for factor, exp in other.factors.iteritems():
            if factor in factors:
                exp = max(exp, factors[factor])

            factors[factor] = exp

        return Factors(factors)

    def __mul__(self, other):
        if isinstance(other, Factors):
            return self.mul(other)
        else:
            return NotImplemented

    def __divmod__(self, other):
        if isinstance(other, Factors):
            return self.div(other)
        else:
            return NotImplemented

    def __div__(self, other):
        if isinstance(other, Factors):
            return self.quo(other)
        else:
            return NotImplemented

    __truediv__ = __div__

    def __mod__(self, other):
        if isinstance(other, Factors):
            return self.rem(other)
        else:
            return NotImplemented

    def __pow__(self, other):
        if type(other) is int:
            return self.pow(other)
        else:
            return NotImplemented

    def __eq__(self, other):
        return self.factors == other.factors

    def __ne__(self, other):
        return not self.__eq__(other)

class Term(object):
    """Efficient representation of ``coeff*(numer/denom)``. """

    __slots__ = ['coeff', 'numer', 'denom']

    def __init__(self, term, numer=None, denom=None):
        if numer is None and denom is None:
            if not term.is_commutative:
                raise NonCommutativeExpression('commutative expression expected')

            coeff, factors = term.as_coeff_mul()
            numer, denom = {}, {}

            for factor in factors:
                base, exp = decompose_power(factor)

                if base.is_Add:
                    cont, base = base.primitive()
                    coeff *= cont

                if exp > 0:
                    numer[base] = exp
                else:
                    denom[base] = -exp

            numer = Factors(numer)
            denom = Factors(denom)
        else:
            coeff = term

            if numer is None:
                numer = Factors()

            if denom is None:
                denom = Factors()

        self.coeff = coeff
        self.numer = numer
        self.denom = denom

    def __repr__(self):
        return "Term(%s, %s, %s)" % (self.coeff, self.numer, self.denom)

    def as_expr(self):
        return self.coeff*(self.numer.as_expr()/self.denom.as_expr())

    def mul(self, other):
        coeff = self.coeff*other.coeff
        numer = self.numer.mul(other.numer)
        denom = self.denom.mul(other.denom)

        numer, denom = numer.normal(denom)

        return Term(coeff, numer, denom)

    def inv(self):
        return Term(1/self.coeff, self.denom, self.numer)

    def quo(self, other):
        return self.mul(other.inv())

    def pow(self, other):
        if other < 0:
            return self.inv().pow(-other)
        else:
            return Term(self.coeff **  other,
                        self.numer.pow(other),
                        self.denom.pow(other))

    def gcd(self, other):
        return Term(self.coeff.gcd(other.coeff),
                    self.numer.gcd(other.numer),
                    self.denom.gcd(other.denom))

    def lcm(self, other):
        return Term(self.coeff.lcm(other.coeff),
                    self.numer.lcm(other.numer),
                    self.denom.lcm(other.denom))

    def __mul__(self, other):
        if isinstance(other, Term):
            return self.mul(other)
        else:
            return NotImplemented

    def __div__(self, other):
        if isinstance(other, Term):
            return self.quo(other)
        else:
            return NotImplemented

    __truediv__ = __div__

    def __pow__(self, other):
        if type(other) is int:
            return self.pow(other)
        else:
            return NotImplemented

    def __eq__(self, other):
        return (self.coeff == other.coeff and
                self.numer == other.numer and
                self.denom == other.denom)

    def __ne__(self, other):
        return not self.__eq__(other)
