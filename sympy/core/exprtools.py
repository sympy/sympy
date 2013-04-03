"""Tools for manipulating of large commutative expressions. """

from sympy.core.add import Add
from sympy.core.compatibility import iterable, is_sequence, SYMPY_INTS
from sympy.core.mul import Mul, _keep_coeff
from sympy.core.power import Pow
from sympy.core.basic import Basic, preorder_traversal
from sympy.core.expr import Expr
from sympy.core.sympify import sympify
from sympy.core.numbers import Rational, Integer, Number
from sympy.core.singleton import S
from sympy.core.symbol import Dummy
from sympy.core.coreerrors import NonCommutativeExpression
from sympy.core.containers import Tuple, Dict
from sympy.utilities import default_sort_key
from sympy.utilities.iterables import (common_prefix, common_suffix,
        variations, ordered)

from collections import defaultdict


def decompose_power(expr):
    """
    Decompose power into symbolic base and integer exponent.

    This is strictly only valid if the exponent from which
    the integer is extracted is itself an integer or the
    base is positive. These conditions are assumed and not
    checked here.

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
        exp, tail = exp.as_coeff_Mul(rational=True)

        if exp is S.NegativeOne:
            base, exp = Pow(base, tail), -1
        elif exp is not S.One:
            tail = _keep_coeff(Rational(1, exp.q), tail)
            base, exp = Pow(base, tail), exp.p
        else:
            base, exp = expr, 1

    return base, exp


class Factors(object):
    """Efficient representation of ``f_1*f_2*...*f_n``. """

    __slots__ = ['factors', 'gens']

    def __init__(self, factors=None):  # Factors
        """Initialize Factors from dict or expr.

        Examples
        ========

        >>> from sympy.core.exprtools import Factors
        >>> from sympy.abc import x
        >>> from sympy import I
        >>> e = 2*x**3
        >>> Factors(e)
        Factors({2: 1, x: 3})
        >>> Factors(e.as_powers_dict())
        Factors({2: 1, x: 3})
        >>> f = _
        >>> f.factors  # underlying dictionary
        {2: 1, x: 3}
        >>> f.gens  # base of each factor
        frozenset([2, x])
        >>> Factors(0)
        Factors({0: 1})
        >>> Factors(I)
        Factors({I: 1})

        Notes
        =====

        Although a dictionary can be passed, no checking will be done to
        ensure that -1 is split from Rationals. This will then affect
        gcd and normal calculations:

        >>> Factors(-2*x).factors
        {-1: 1, 2: 1, x: 1}
        >>> Factors((-2*x).as_powers_dict()).factors
        {-2: 1, x: 1}

        """
        if factors is None or factors == 1:
            factors = {}
        elif factors is S.Zero or factors == 0:
            factors = {S.Zero: S.One}
        elif isinstance(factors, Number) and factors < 0 and not factors is S.NegativeOne:
            factors = {-factors: S.One, S.NegativeOne: S.One}
        elif isinstance(factors, Expr):
            c, nc = factors.args_cnc()
            factors = dict(Mul(*c, **dict(evaluate=False)).as_powers_dict())
            if nc:
                factors[Mul(*nc, **dict(evaluate=False))] = S.One

        # sanitizing
        if factors.get(S.NegativeOne, S.One) is not S.One:
            e = factors.pop(S.NegativeOne)
            factors[S.ImaginaryUnit] = 1
            if e == S(3)/2:
                factors[S.NegativeOne] = S.One
            elif e is not S.Half:
                raise ValueError('unanticipated value for exponent of -1')

        self.factors = factors
        try:
            self.gens = frozenset(factors.keys())
        except AttributeError:
            raise TypeError('expecting Expr or dictionary')

    def __hash__(self):  # Factors
        keys = tuple(ordered(self.factors.keys()))
        values = [self.factors[k] for k in keys]
        return hash((keys, values))

    def __repr__(self):  # Factors
        return "Factors({%s})" % ', '.join(
            ['%s: %s' % (k, v) for k, v in ordered(self.factors.items())])

    @property
    def is_zero(self):
        """
        >>> from sympy.core.exprtools import Factors
        >>> Factors(1).div(Factors(1))[1].is_zero
        True
        """
        f = self.factors
        return len(f) == 1 and S.Zero in f

    @property
    def is_one(self):
        """
        >>> from sympy.core.exprtools import Factors
        >>> Factors(1).is_one
        True
        """
        return not self.factors

    def as_expr(self):  # Factors
        """Return the underlying expression.

        Examples
        ========

        >>> from sympy.core.exprtools import Factors
        >>> from sympy.abc import x, y, z
        >>> Factors((x*y**2).as_powers_dict()).as_expr()
        x*y**2

        """

        args = []
        for factor, exp in self.factors.iteritems():
            if exp != 1:
                b, e = factor.as_base_exp()
                if isinstance(exp, int):
                    e = _keep_coeff(Integer(exp), e)
                elif isinstance(exp, Rational):
                    e = _keep_coeff(exp, e)
                else:
                    e *= exp
                args.append(b**e)
            else:
                args.append(factor)
        return Mul(*args)

    def normal(self, other):  # Factors
        """Return unique Factors of self and other.

        Examples
        ========

        >>> from sympy.core.exprtools import Factors
        >>> from sympy.abc import x, y, z
        >>> a = Factors(2*x*y**2)
        >>> b = Factors(-2*x*y/z)
        >>> a.normal(b)
        (Factors({y: 1}), Factors({-1: 1, z: -1}))
        >>> a.normal(a)
        (Factors({}), Factors({}))

        """
        if not isinstance(other, Factors):
            other = Factors(other)
        if other.is_zero:
            return Factors(), Factors(S.Zero)

        self_factors = dict(self.factors)
        other_factors = dict(other.factors)
        if self_factors == other_factors:
            return Factors(), Factors()

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

    def mul(self, other):  # Factors
        """Return Factors of ``self * other``.

        Examples
        ========

        >>> from sympy.core.exprtools import Factors
        >>> from sympy.abc import x, y, z
        >>> a = Factors((x*y**2).as_powers_dict())
        >>> b = Factors((x*y/z).as_powers_dict())
        >>> a.mul(b)
        Factors({x: 2, y: 3, z: -1})
        >>> a*b
        Factors({x: 2, y: 3, z: -1})
        """
        if not isinstance(other, Factors):
            other = Factors(other)
        if any(f.is_zero for f in (self, other)):
            return Factors(S.Zero)
        factors = dict(self.factors)

        for factor, exp in other.factors.iteritems():
            if factor in factors:
                exp = factors[factor] + exp

                if not exp:
                    del factors[factor]
                    continue

                if factor is S.NegativeOne:
                    if not exp % 2:
                        del factors[factor]
                        continue
                    exp = S.One

            factors[factor] = exp

        return Factors(factors)

    def div(self, other):  # Factors
        """Return the result of ``divmod(self, other)`` as two Factors,
        ``quo`` and ``rem`` such that ``self = other*(quo + rem)``.

        Notes
        =====

        If there is no remainder, the factors will not be empty, it will
        be {0: 1}; the empty factors corresponds to {1: 1}. This is similar
        to the result of string.find returning 0 for a hit and -1 for a
        non-hit. (Normally {} and 0 mean "null" but with Factors and find
        they have different semantics.)

        Examples
        ========

        >>> from sympy.core.exprtools import Factors
        >>> from sympy.abc import x, y, z
        >>> a = Factors((x*y**2).as_powers_dict())
        >>> b = Factors((x*y/z).as_powers_dict())
        >>> a.div(b)
        (Factors({y: 1}), Factors({z: -1}))
        >>> q, r = a.div(a)
        >>> a*q, a*r
        (Factors({x: 1, y: 2}), Factors({0: 1}))

        The ``/`` operator only gives the quotient:

        >>> a/b
        Factors({y: 1})
        """
        if not isinstance(other, Factors):
            other = Factors(other)
            if other.is_zero:
                raise ZeroDivisionError
            if self.is_zero:
                return Factors(S.Zero)
        quo, rem = dict(self.factors), {}

        for factor, exp in other.factors.iteritems():
            if factor in quo:
                if factor is S.NegativeOne:
                    del quo[factor]
                    continue

                exp = quo[factor] - exp

                if exp <= 0:
                    del quo[factor]

                if exp >= 0:
                    if exp:
                        quo[factor] = exp

                    continue

                exp = -exp

            rem[factor] = exp

        rem = rem or S.Zero

        return Factors(quo), Factors(rem)

    def quo(self, other):  # Factors
        """Return quotient Factors of ``self / other``.

        Examples
        ========

        >>> from sympy.core.exprtools import Factors
        >>> from sympy.abc import x, y, z
        >>> a = Factors((x*y**2).as_powers_dict())
        >>> b = Factors((x*y/z).as_powers_dict())
        >>> a.quo(b)  # same as a/b
        Factors({y: 1})
        """
        return self.div(other)[0]

    def rem(self, other):  # Factors
        """Return remainder Factors of ``self / other``.

        Examples
        ========

        >>> from sympy.core.exprtools import Factors
        >>> from sympy.abc import x, y, z
        >>> a = Factors((x*y**2).as_powers_dict())
        >>> b = Factors((x*y/z).as_powers_dict())
        >>> a.rem(b)
        Factors({z: -1})
        >>> a.rem(a)
        Factors({0: 1})
        """
        return self.div(other)[1]

    def pow(self, other):  # Factors
        """Return self raised to a non-negative integer power.

        Examples
        ========

        >>> from sympy.core.exprtools import Factors
        >>> from sympy.abc import x, y, z
        >>> a = Factors((x*y**2).as_powers_dict())
        >>> a**2
        Factors({x: 2, y: 4})

        """
        if isinstance(other, Factors):
            other = other.as_expr()
            if other.is_Integer:
                other = int(other)
        if isinstance(other, SYMPY_INTS) and other >= 0:
            factors = {}

            if other:
                for factor, exp in self.factors.iteritems():
                    factors[factor] = exp*other
                    if factor is S.NegativeOne:
                        if factors[factor] % 2 == 1:
                            factors[factor] = S.One
                        else:
                            del factors[factor]

            return Factors(factors)
        else:
            raise ValueError("expected non-negative integer, got %s" % other)

    def gcd(self, other):  # Factors
        """Return Factors of ``gcd(self, other)``. The keys are
        the intersection of factors with the minimum exponent for
        each factor.

        Examples
        ========

        >>> from sympy.core.exprtools import Factors
        >>> from sympy.abc import x, y, z
        >>> a = Factors((x*y**2).as_powers_dict())
        >>> b = Factors((x*y/z).as_powers_dict())
        >>> a.gcd(b)
        Factors({x: 1, y: 1})
        """
        if not isinstance(other, Factors):
            other = Factors(other)
            if other.is_zero:
                return Factors(self.factors)

        factors = {}

        for factor, exp in self.factors.iteritems():
            if factor in other.factors:
                exp = min(exp, other.factors[factor])
                factors[factor] = exp

        return Factors(factors)

    def lcm(self, other):  # Factors
        """Return Factors of ``lcm(self, other)`` which are
        the union of factors with the maximum exponent for
        each factor.

        Examples
        ========

        >>> from sympy.core.exprtools import Factors
        >>> from sympy.abc import x, y, z
        >>> a = Factors((x*y**2).as_powers_dict())
        >>> b = Factors((x*y/z).as_powers_dict())
        >>> a.lcm(b)
        Factors({x: 1, y: 2, z: -1})
        """
        if not isinstance(other, Factors):
            other = Factors(other)
            if any(f.is_zero for f in (self, other)):
                return Factors(S.Zero)

        factors = dict(self.factors)

        for factor, exp in other.factors.iteritems():
            if factor in factors:
                exp = max(exp, factors[factor])

            factors[factor] = exp

        return Factors(factors)

    def __mul__(self, other):  # Factors
        return self.mul(other)

    def __divmod__(self, other):  # Factors
        return self.div(other)

    def __div__(self, other):  # Factors
        return self.quo(other)

    __truediv__ = __div__

    def __mod__(self, other):  # Factors
        return self.rem(other)

    def __pow__(self, other):  # Factors
        return self.pow(other)

    def __eq__(self, other):  # Factors
        if not isinstance(other, Factors):
            other = Factors(other)
        return self.factors == other.factors

    def __ne__(self, other):  # Factors
        return not self.__eq__(other)


class Term(object):
    """Efficient representation of ``coeff*(numer/denom)``. """

    __slots__ = ['coeff', 'numer', 'denom']

    def __init__(self, term, numer=None, denom=None):  # Term
        if numer is None and denom is None:
            if not term.is_commutative:
                raise NonCommutativeExpression(
                    'commutative expression expected')

            coeff, factors = term.as_coeff_mul()
            numer, denom = defaultdict(int), defaultdict(int)

            for factor in factors:
                base, exp = decompose_power(factor)

                if base.is_Add:
                    cont, base = base.primitive()
                    coeff *= cont**exp

                if exp > 0:
                    numer[base] += exp
                else:
                    denom[base] += -exp

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

    def __hash__(self):  # Term
        return hash((self.coeff, self.numer, self.denom))

    def __repr__(self):  # Term
        return "Term(%s, %s, %s)" % (self.coeff, self.numer, self.denom)

    def as_expr(self):  # Term
        return self.coeff*(self.numer.as_expr()/self.denom.as_expr())

    def mul(self, other):  # Term
        coeff = self.coeff*other.coeff
        numer = self.numer.mul(other.numer)
        denom = self.denom.mul(other.denom)

        numer, denom = numer.normal(denom)

        return Term(coeff, numer, denom)

    def inv(self):  # Term
        return Term(1/self.coeff, self.denom, self.numer)

    def quo(self, other):  # Term
        return self.mul(other.inv())

    def pow(self, other):  # Term
        if other < 0:
            return self.inv().pow(-other)
        else:
            return Term(self.coeff ** other,
                        self.numer.pow(other),
                        self.denom.pow(other))

    def gcd(self, other):  # Term
        return Term(self.coeff.gcd(other.coeff),
                    self.numer.gcd(other.numer),
                    self.denom.gcd(other.denom))

    def lcm(self, other):  # Term
        return Term(self.coeff.lcm(other.coeff),
                    self.numer.lcm(other.numer),
                    self.denom.lcm(other.denom))

    def __mul__(self, other):  # Term
        if isinstance(other, Term):
            return self.mul(other)
        else:
            return NotImplemented

    def __div__(self, other):  # Term
        if isinstance(other, Term):
            return self.quo(other)
        else:
            return NotImplemented

    __truediv__ = __div__

    def __pow__(self, other):  # Term
        if isinstance(other, SYMPY_INTS):
            return self.pow(other)
        else:
            return NotImplemented

    def __eq__(self, other):  # Term
        return (self.coeff == other.coeff and
                self.numer == other.numer and
                self.denom == other.denom)

    def __ne__(self, other):  # Term
        return not self.__eq__(other)


def _gcd_terms(terms, isprimitive=False, fraction=True):
    """Helper function for :func:`gcd_terms`.

    If ``isprimitive`` is True then the call to primitive
    for an Add will be skipped. This is useful when the
    content has already been extrated.

    If ``fraction`` is True then the expression will appear over a common
    denominator, the lcm of all term denominators.
    """

    if isinstance(terms, Basic) and not isinstance(terms, Tuple):
        terms = Add.make_args(terms)

    terms = map(Term, [t for t in terms if t])

    # there is some simplification that may happen if we leave this
    # here rather than duplicate it before the mapping of Term onto
    # the terms
    if len(terms) == 0:
        return S.Zero, S.Zero, S.One

    if len(terms) == 1:
        cont = terms[0].coeff
        numer = terms[0].numer.as_expr()
        denom = terms[0].denom.as_expr()

    else:
        cont = terms[0]
        for term in terms[1:]:
            cont = cont.gcd(term)

        for i, term in enumerate(terms):
            terms[i] = term.quo(cont)

        if fraction:
            denom = terms[0].denom

            for term in terms[1:]:
                denom = denom.lcm(term.denom)

            numers = []
            for term in terms:
                numer = term.numer.mul(denom.quo(term.denom))
                numers.append(term.coeff*numer.as_expr())
        else:
            numers = [t.as_expr() for t in terms]
            denom = Term(S(1)).numer

        cont = cont.as_expr()
        numer = Add(*numers)
        denom = denom.as_expr()

    if not isprimitive and numer.is_Add:
        _cont, numer = numer.primitive()
        cont *= _cont

    return cont, numer, denom


def gcd_terms(terms, isprimitive=False, clear=True, fraction=True):
    """Compute the GCD of ``terms`` and put them together.

    ``terms`` can be an expression or a non-Basic sequence of expressions
    which will be handled as though they are terms from a sum.

    If ``isprimitive`` is True the _gcd_terms will not run the primitive
    method on the terms.

    ``clear`` controls the removal of integers from the denominator of an Add
    expression. When True (default), all numerical denominator will be cleared;
    when False the denominators will be cleared only if all terms had numerical
    denominators other than 1.

    ``fraction``, when True (default), will put the expression over a common
    denominator.

    Examples
    ========

    >>> from sympy.core import gcd_terms
    >>> from sympy.abc import x, y

    >>> gcd_terms((x + 1)**2*y + (x + 1)*y**2)
    y*(x + 1)*(x + y + 1)
    >>> gcd_terms(x/2 + 1)
    (x + 2)/2
    >>> gcd_terms(x/2 + 1, clear=False)
    x/2 + 1
    >>> gcd_terms(x/2 + y/2, clear=False)
    (x + y)/2
    >>> gcd_terms(x/2 + 1/x)
    (x**2 + 2)/(2*x)
    >>> gcd_terms(x/2 + 1/x, fraction=False)
    (x + 2/x)/2
    >>> gcd_terms(x/2 + 1/x, fraction=False, clear=False)
    x/2 + 1/x

    >>> gcd_terms(x/2/y + 1/x/y)
    (x**2 + 2)/(2*x*y)
    >>> gcd_terms(x/2/y + 1/x/y, fraction=False, clear=False)
    (x + 2/x)/(2*y)

    See Also
    ========
    factor_terms, sympy.polys.polytools.terms_gcd

    """
    def mask(terms):
        """replace nc portions of each term with a unique Dummy symbols
        and return the replacements to restore them"""
        args = [(a, []) if a.is_commutative else a.args_cnc() for a in terms]
        reps = []
        for i, (c, nc) in enumerate(args):
            if nc:
                nc = Mul._from_args(nc)
                d = Dummy()
                reps.append((d, nc))
                c.append(d)
                args[i] = Mul._from_args(c)
            else:
                args[i] = c
        return args, dict(reps)

    isadd = isinstance(terms, Add)
    addlike = isadd or not isinstance(terms, Basic) and \
        is_sequence(terms, include=set) and \
        not isinstance(terms, Dict)

    if addlike:
        if isadd:  # i.e. an Add
            terms = list(terms.args)
        else:
            terms = sympify(terms)
        terms, reps = mask(terms)
        cont, numer, denom = _gcd_terms(terms, isprimitive, fraction)
        numer = numer.xreplace(reps)
        coeff, factors = cont.as_coeff_Mul()
        return _keep_coeff(coeff, factors*numer/denom, clear=clear)

    if not isinstance(terms, Basic):
        return terms

    if terms.is_Atom:
        return terms

    if terms.is_Mul:
        c, args = terms.as_coeff_mul()
        return _keep_coeff(c, Mul(*[gcd_terms(i, isprimitive, clear, fraction)
            for i in args]), clear=clear)

    def handle(a):
        # don't treat internal args like terms of an Add
        if not isinstance(a, Expr):
            if isinstance(a, Basic):
                return a.func(*[handle(i) for i in a.args])
            return type(a)([handle(i) for i in a])
        return gcd_terms(a, isprimitive, clear, fraction)

    if isinstance(terms, Dict):
        return Dict(*[(k, handle(v)) for k, v in terms.args])
    return terms.func(*[handle(i) for i in terms.args])


def factor_terms(expr, radical=False, clear=False, fraction=False, sign=True):
    """Remove common factors from terms in all arguments without
    changing the underlying structure of the expr. No expansion or
    simplification (and no processing of non-commutatives) is performed.

    If radical=True then a radical common to all terms will be factored
    out of any Add sub-expressions of the expr.

    If clear=False (default) then coefficients will not be separated
    from a single Add if they can be distributed to leave one or more
    terms with integer coefficients.

    If fraction=True (default is False) then a common denominator will be
    constructed for the expression.

    If sign=True (default) then even if the only factor in common is a -1,
    it will be factored out of the expression.

    Examples
    ========

    >>> from sympy import factor_terms, Symbol, Mul, primitive
    >>> from sympy.abc import x, y
    >>> factor_terms(x + x*(2 + 4*y)**3)
    x*(8*(2*y + 1)**3 + 1)
    >>> A = Symbol('A', commutative=False)
    >>> factor_terms(x*A + x*A + x*y*A)
    x*(y*A + 2*A)

    When ``clear`` is False, a rational will only be factored out of an
    Add expression if all terms of the Add have coefficients that are
    fractions:

    >>> factor_terms(x/2 + 1, clear=False)
    x/2 + 1
    >>> factor_terms(x/2 + 1, clear=True)
    (x + 2)/2

    This only applies when there is a single Add that the coefficient
    multiplies:

    >>> factor_terms(x*y/2 + y, clear=True)
    y*(x + 2)/2
    >>> factor_terms(x*y/2 + y, clear=False) == _
    True

    If a -1 is all that can be factored out, to *not* factor it out, the
    flag ``sign`` must be False:

    >>> factor_terms(-x - y)
    -(x + y)
    >>> factor_terms(-x - y, sign=False)
    -x - y
    >>> factor_terms(-2*x - 2*y, sign=False)
    -2*(x + y)

    See Also
    ========
    gcd_terms, sympy.polys.polytools.terms_gcd

    """

    expr = sympify(expr)
    is_iterable = iterable(expr)

    if not isinstance(expr, Basic) or expr.is_Atom:
        if is_iterable:
            return type(expr)([factor_terms(i,
                radical=radical,
                clear=clear,
                fraction=fraction) for i in expr])
        return expr

    if expr.is_Pow or expr.is_Function or is_iterable or not hasattr(expr, 'args_cnc'):
        args = expr.args
        newargs = tuple([factor_terms(i,
            radical=radical,
            clear=clear,
            fraction=fraction) for i in args])
        if newargs == args:
            return expr
        return expr.func(*newargs)

    cont, p = expr.as_content_primitive(radical=radical)
    if p.is_Add:
        list_args = [gcd_terms(a,
        isprimitive=True,
        clear=clear,
        fraction=fraction) for a in Add.make_args(p)]
        # get a common negative (if there) which gcd_terms does not remove
        if all(a.as_coeff_Mul()[0] < 0 for a in list_args):
            cont = -cont
            list_args = [-a for a in list_args]
        p = Add._from_args(list_args)  # gcd_terms will fix up ordering
    elif p.args:
        p = p.func(
            *[factor_terms(a, radical, clear, fraction) for a in p.args])
    p = gcd_terms(p,
        isprimitive=True,
        clear=clear,
        fraction=fraction)
    return _keep_coeff(cont, p, clear=clear, sign=sign)


def _mask_nc(eq):
    """Return ``eq`` with non-commutative objects replaced with dummy
    symbols. A dictionary that can be used to restore the original
    values is returned: if it is None, the expression is
    noncommutative and cannot be made commutative. The third value
    returned is a list of any non-commutative symbols that appear
    in the returned equation.

    Notes
    =====
    All non-commutative objects other than Symbols are replaced with
    a non-commutative Symbol. Identical objects will be identified
    by identical symbols.

    If there is only 1 non-commutative object in an expression it will
    be replaced with a commutative symbol. Otherwise, the non-commutative
    entities are retained and the calling routine should handle
    replacements in this case since some care must be taken to keep
    track of the ordering of symbols when they occur within Muls.

    Examples
    ========
    >>> from sympy.physics.secondquant import Commutator, NO, F, Fd
    >>> from sympy import Dummy, symbols, Sum, Mul, Basic
    >>> from sympy.abc import x, y
    >>> from sympy.core.exprtools import _mask_nc
    >>> A, B, C = symbols('A,B,C', commutative=False)
    >>> Dummy._count = 0 # reset for doctest purposes

    One nc-symbol:

    >>> _mask_nc(A**2 - x**2)
    (_0**2 - x**2, {_0: A}, [])

    Multiple nc-symbols:

    >>> _mask_nc(A**2 - B**2)
    (A**2 - B**2, None, [A, B])

    An nc-object with nc-symbols but no others outside of it:

    >>> _mask_nc(1 + x*Commutator(A, B))
    (_1*x + 1, {_1: Commutator(A, B)}, [])
    >>> _mask_nc(NO(Fd(x)*F(y)))
    (_2, {_2: NO(CreateFermion(x)*AnnihilateFermion(y))}, [])

    Multiple nc-objects:

    >>> eq = x*Commutator(A, B) + x*Commutator(A, C)*Commutator(A, B)
    >>> _mask_nc(eq)
    (x*_3 + x*_4*_3, {_3: Commutator(A, B), _4: Commutator(A, C)}, [_3, _4])

    Multiple nc-objects and nc-symbols:

    >>> eq = A*Commutator(A, B) + B*Commutator(A, C)
    >>> _mask_nc(eq)
    (A*_5 + B*_6, {_5: Commutator(A, B), _6: Commutator(A, C)}, [_5, _6, A, B])

    If there is an object that:

        - doesn't contain nc-symbols
        - but has arguments which derive from Basic, not Expr
        - and doesn't define an _eval_is_commutative routine

    then it will give False (or None?) for the is_commutative test. Such
    objects are also removed by this routine:

    >>> from sympy import Basic, Mul
    >>> eq = (1 + Mul(Basic(), Basic(), evaluate=False))
    >>> eq.is_commutative
    False
    >>> _mask_nc(eq)
    (_7**2 + 1, {_7: Basic()}, [])

    """
    expr = eq
    if expr.is_commutative:
        return eq, {}, []

    # identify nc-objects; symbols and other
    rep = []
    nc_obj = set()
    nc_syms = set()
    pot = preorder_traversal(expr, keys=default_sort_key)
    for i, a in enumerate(pot):
        if any(a == r[0] for r in rep):
            pot.skip()
        elif not a.is_commutative:
            if a.is_Symbol:
                nc_syms.add(a)
            elif not (a.is_Add or a.is_Mul or a.is_Pow):
                if all(s.is_commutative for s in a.free_symbols):
                    rep.append((a, Dummy()))
                else:
                    nc_obj.add(a)
                pot.skip()

    # If there is only one nc symbol or object, it can be factored regularly
    # but polys is going to complain, so replace it with a Dummy.
    if len(nc_obj) == 1 and not nc_syms:
        rep.append((nc_obj.pop(), Dummy()))
    elif len(nc_syms) == 1 and not nc_obj:
        rep.append((nc_syms.pop(), Dummy()))

    # Any remaining nc-objects will be replaced with an nc-Dummy and
    # identified as an nc-Symbol to watch out for
    nc_obj = sorted(nc_obj, key=default_sort_key)
    for n in nc_obj:
        nc = Dummy(commutative=False)
        rep.append((n, nc))
        nc_syms.add(nc)
    expr = expr.subs(rep)

    nc_syms = list(nc_syms)
    nc_syms.sort(key=default_sort_key)
    return expr, dict([(v, k) for k, v in rep]) or None, nc_syms


def factor_nc(expr):
    """Return the factored form of ``expr`` while handling non-commutative
    expressions.

    **examples**
    >>> from sympy.core.exprtools import factor_nc
    >>> from sympy import Symbol
    >>> from sympy.abc import x
    >>> A = Symbol('A', commutative=False)
    >>> B = Symbol('B', commutative=False)
    >>> factor_nc((x**2 + 2*A*x + A**2).expand())
    (x + A)**2
    >>> factor_nc(((x + A)*(x + B)).expand())
    (x + A)*(x + B)
    """
    from sympy.simplify.simplify import powsimp
    from sympy.polys import gcd, factor

    def _pemexpand(expr):
        "Expand with the minimal set of hints necessary to check the result."
        return expr.expand(deep=True, mul=True, power_exp=True,
            power_base=False, basic=False, multinomial=True, log=False)

    expr = sympify(expr)
    if not isinstance(expr, Expr) or not expr.args:
        return expr
    if not expr.is_Add:
        return expr.func(*[factor_nc(a) for a in expr.args])

    expr, rep, nc_symbols = _mask_nc(expr)
    if rep:
        return factor(expr).subs(rep)
    else:
        args = [a.args_cnc() for a in Add.make_args(expr)]
        c = g = l = r = S.One
        hit = False
        # find any commutative gcd term
        for i, a in enumerate(args):
            if i == 0:
                c = Mul._from_args(a[0])
            elif a[0]:
                c = gcd(c, Mul._from_args(a[0]))
            else:
                c = S.One
        if c is not S.One:
            hit = True
            c, g = c.as_coeff_Mul()
            if g is not S.One:
                for i, (cc, _) in enumerate(args):
                    cc = list(Mul.make_args(Mul._from_args(list(cc))/g))
                    args[i][0] = cc
            else:
                for i, (cc, _) in enumerate(args):
                    cc[0] = cc[0]/c
                    args[i][0] = cc
        # find any noncommutative common prefix
        for i, a in enumerate(args):
            if i == 0:
                n = a[1][:]
            else:
                n = common_prefix(n, a[1])
            if not n:
                # is there a power that can be extracted?
                if not args[0][1]:
                    break
                b, e = args[0][1][0].as_base_exp()
                ok = False
                if e.is_Integer:
                    for t in args:
                        if not t[1]:
                            break
                        bt, et = t[1][0].as_base_exp()
                        if et.is_Integer and bt == b:
                            e = min(e, et)
                        else:
                            break
                    else:
                        ok = hit = True
                        l = b**e
                        il = b**-e
                        for i, a in enumerate(args):
                            args[i][1][0] = il*args[i][1][0]
                        break
                if not ok:
                    break
        else:
            hit = True
            lenn = len(n)
            l = Mul(*n)
            for i, a in enumerate(args):
                args[i][1] = args[i][1][lenn:]
        # find any noncommutative common suffix
        for i, a in enumerate(args):
            if i == 0:
                n = a[1][:]
            else:
                n = common_suffix(n, a[1])
            if not n:
                # is there a power that can be extracted?
                if not args[0][1]:
                    break
                b, e = args[0][1][-1].as_base_exp()
                ok = False
                if e.is_Integer:
                    for t in args:
                        if not t[1]:
                            break
                        bt, et = t[1][-1].as_base_exp()
                        if et.is_Integer and bt == b:
                            e = min(e, et)
                        else:
                            break
                    else:
                        ok = hit = True
                        r = b**e
                        il = b**-e
                        for i, a in enumerate(args):
                            args[i][1][-1] = args[i][1][-1]*il
                        break
                if not ok:
                    break
        else:
            hit = True
            lenn = len(n)
            r = Mul(*n)
            for i, a in enumerate(args):
                args[i][1] = a[1][:len(a[1]) - lenn]
        if hit:
            mid = Add(*[Mul(*cc)*Mul(*nc) for cc, nc in args])
        else:
            mid = expr

        # sort the symbols so the Dummys would appear in the same
        # order as the original symbols, otherwise you may introduce
        # a factor of -1, e.g. A**2 - B**2) -- {A:y, B:x} --> y**2 - x**2
        # and the former factors into two terms, (A - B)*(A + B) while the
        # latter factors into 3 terms, (-1)*(x - y)*(x + y)
        rep1 = [(n, Dummy()) for n in sorted(nc_symbols, key=default_sort_key)]
        unrep1 = [(v, k) for k, v in rep1]
        unrep1.reverse()
        new_mid, r2, _ = _mask_nc(mid.subs(rep1))
        new_mid = powsimp(factor(new_mid))

        new_mid = new_mid.subs(r2).subs(unrep1)

        if new_mid.is_Pow:
            return _keep_coeff(c, g*l*new_mid*r)

        if new_mid.is_Mul:
            # XXX TODO there should be a way to inspect what order the terms
            # must be in and just select the plausible ordering without
            # checking permutations
            cfac = []
            ncfac = []
            for f in new_mid.args:
                if f.is_commutative:
                    cfac.append(f)
                else:
                    b, e = f.as_base_exp()
                    if e.is_Integer:
                        ncfac.extend([b]*e)
                    else:
                        ncfac.append(f)
            pre_mid = g*Mul(*cfac)*l
            target = _pemexpand(expr/c)
            for s in variations(ncfac, len(ncfac)):
                ok = pre_mid*Mul(*s)*r
                if _pemexpand(ok) == target:
                    return _keep_coeff(c, ok)

        # mid was an Add that didn't factor successfully
        return _keep_coeff(c, g*l*mid*r)
