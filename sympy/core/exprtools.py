"""Tools for manipulating of large commutative expressions. """

from sympy.core.add import Add
from sympy.core.compatibility import iterable
from sympy.core.mul import Mul, _keep_coeff
from sympy.core.power import Pow
from sympy.core.basic import Basic
from sympy.core.expr import Expr
from sympy.core.sympify import sympify
from sympy.core.numbers import Rational, Integer
from sympy.core.singleton import S
from sympy.core.symbol import Dummy
from sympy.core.coreerrors import NonCommutativeExpression
from sympy.core.containers import Tuple

def decompose_power(expr):
    """
    Decompose power into symbolic base and integer exponent.

    Examples
    ========

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
            tail = _keep_coeff(Rational(1, exp.q), Mul(*tail))

            base, exp = Pow(base, tail), exp.p
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

    def __hash__(self):
        return hash((tuple(self.factors), self.gens))

    def __repr__(self):
        return "Factors(%s)" % self.factors

    def as_expr(self):
        args = []
        for factor, exp in self.factors.iteritems():
            if exp != 1:
                b, e = factor.as_base_exp()
                e = _keep_coeff(Integer(exp), e)
                args.append(b**e)
            else:
                args.append(factor)
        return Mul(*args)

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
                    coeff *= cont**exp

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

    def __hash__(self):
        return hash((self.coeff, self.numer, self.denom))

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

def _gcd_terms(terms, isprimitive=False):
    """Helper function for :func:`gcd_terms`. If `isprimitive` is True then the
    call to primitive for an Add will be skipped. This is useful when the
    content has already been extrated."""
    if isinstance(terms, Basic) and not isinstance(terms, Tuple):
        terms = Add.make_args(terms)

    if len(terms) <= 1:
        if not terms:
            return S.Zero, S.Zero, S.One
        else:
            return terms[0], S.One, S.One

    terms = map(Term, terms)
    cont = terms[0]

    for term in terms[1:]:
        cont = cont.gcd(term)

    for i, term in enumerate(terms):
        terms[i] = term.quo(cont)

    denom = terms[0].denom

    for term in terms[1:]:
        denom = denom.lcm(term.denom)

    numers = []

    for term in terms:
        numer = term.numer.mul(denom.quo(term.denom))
        numers.append(term.coeff*numer.as_expr())

    cont = cont.as_expr()
    numer = Add(*numers)
    denom = denom.as_expr()
    if not isprimitive and numer.is_Add:
        _cont, numer = numer.primitive()
        cont *= _cont

    return cont, numer, denom

def gcd_terms(terms, isprimitive=False):
    """
    Compute the GCD of ``terms`` and put them together. If ``isprimitive`` is
    True the _gcd_terms will not run the primitive method on the terms.

    Examples
    ========

    >>> from sympy.core import gcd_terms
    >>> from sympy.abc import x, y

    >>> gcd_terms((x + 1)**2*y + (x + 1)*y**2)
    y*(x + 1)*(x + y + 1)

    """
    terms = sympify(terms)
    isexpr = isinstance(terms, Expr)
    if not isexpr or terms.is_Add:
        cont, numer, denom = _gcd_terms(terms, isprimitive)
        coeff, factors = cont.as_coeff_Mul()
        return _keep_coeff(coeff, factors*numer/denom)

    if terms.is_Atom:
        return terms

    if terms.is_Mul:
        c, args = terms.as_coeff_mul()
        return _keep_coeff(c, Mul(*[gcd_terms(i, isprimitive) for i in args]))

    def handle(a):
        if iterable(a):
            if isinstance(a, Basic):
                return a.func(*[gcd_terms(i, isprimitive) for i in a.args])
            return type(a)([gcd_terms(i, isprimitive) for i in a])
        return gcd_terms(a, isprimitive)
    return terms.func(*[handle(i) for i in terms.args])


def factor_terms(expr, radical=False):
    """Remove common factors from terms in all arguments without
    changing the underlying structure of the expr. No expansion or
    simplification (and no processing of non-commutatives) is performed.

    If radical=True then a radical common to all terms will be factored
    out of any Add sub-expressions of the expr.

    Examples
    ========

    >>> from sympy import factor_terms, Symbol
    >>> from sympy.abc import x, y
    >>> factor_terms(x + x*(2 + 4*y)**3)
    x*(8*(2*y + 1)**3 + 1)
    >>> A = Symbol('A', commutative=False)
    >>> factor_terms(x*A + x*A + x*y*A)
    x*(y*A + 2*A)

    """

    expr = sympify(expr)
    is_iterable = iterable(expr)

    if not isinstance(expr, Basic) or expr.is_Atom:
        if is_iterable:
            return type(expr)([factor_terms(i, radical=radical) for i in expr])
        return expr

    if expr.is_Function or is_iterable or not hasattr(expr, 'args_cnc'):
        args = expr.args
        newargs = tuple([factor_terms(i, radical=radical) for i in args])
        if newargs == args:
            return expr
        return expr.func(*newargs)

    cont, p = expr.as_content_primitive(radical=radical)
    list_args, nc = zip(*[ai.args_cnc() for ai in Add.make_args(p)])
    list_args = list(list_args)
    nc = [((Dummy(), Mul._from_args(i)) if i else None) for i in nc]
    ncreps = dict([i for i in nc if i is not None])
    for i, a in enumerate(list_args):
        if nc[i] is not None:
            a.append(nc[i][0])
        a = Mul._from_args(a) # gcd_terms will fix up ordering
        list_args[i] = gcd_terms(a, isprimitive=True)
        # cancel terms that may not have cancelled
    p = Add._from_args(list_args) # gcd_terms will fix up ordering
    p = gcd_terms(p, isprimitive=True).xreplace(ncreps)
    return _keep_coeff(cont, p)
