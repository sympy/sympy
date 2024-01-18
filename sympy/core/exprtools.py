"""Tools for manipulating of large commutative expressions. """

from .add import Add
from .mul import Mul, _keep_coeff
from .power import Pow
from .basic import Basic
from .expr import Expr
from .function import expand_power_exp
from .sympify import sympify
from .numbers import Rational, Integer, Number, I, equal_valued
from .singleton import S
from .sorting import default_sort_key, ordered
from .symbol import Dummy
from .traversal import preorder_traversal
from .coreerrors import NonCommutativeExpression
from .containers import Tuple, Dict
from sympy.external.gmpy import SYMPY_INTS
from sympy.utilities.iterables import (common_prefix, common_suffix,
        variations, iterable, is_sequence)

from collections import defaultdict
from typing import Tuple as tTuple


_eps = Dummy(positive=True)


def _isnumber(i):
    return isinstance(i, (SYMPY_INTS, float)) or i.is_Number


def _monotonic_sign(self):
    """Return the value closest to 0 that ``self`` may have if all symbols
    are signed and the result is uniformly the same sign for all values of symbols.
    If a symbol is only signed but not known to be an
    integer or the result is 0 then a symbol representative of the sign of self
    will be returned. Otherwise, None is returned if a) the sign could be positive
    or negative or b) self is not in one of the following forms:

    - L(x, y, ...) + A: a function linear in all symbols x, y, ... with an
      additive constant; if A is zero then the function can be a monomial whose
      sign is monotonic over the range of the variables, e.g. (x + 1)**3 if x is
      nonnegative.
    - A/L(x, y, ...) + B: the inverse of a function linear in all symbols x, y, ...
      that does not have a sign change from positive to negative for any set
      of values for the variables.
    - M(x, y, ...) + A: a monomial M whose factors are all signed and a constant, A.
    - A/M(x, y, ...) + B: the inverse of a monomial and constants A and B.
    - P(x): a univariate polynomial

    Examples
    ========

    >>> from sympy.core.exprtools import _monotonic_sign as F
    >>> from sympy import Dummy
    >>> nn = Dummy(integer=True, nonnegative=True)
    >>> p = Dummy(integer=True, positive=True)
    >>> p2 = Dummy(integer=True, positive=True)
    >>> F(nn + 1)
    1
    >>> F(p - 1)
    _nneg
    >>> F(nn*p + 1)
    1
    >>> F(p2*p + 1)
    2
    >>> F(nn - 1)  # could be negative, zero or positive
    """
    if not self.is_extended_real:
        return

    if (-self).is_Symbol:
        rv = _monotonic_sign(-self)
        return rv if rv is None else -rv

    if not self.is_Add and self.as_numer_denom()[1].is_number:
        s = self
        if s.is_prime:
            if s.is_odd:
                return Integer(3)
            else:
                return Integer(2)
        elif s.is_composite:
            if s.is_odd:
                return Integer(9)
            else:
                return Integer(4)
        elif s.is_positive:
            if s.is_even:
                if s.is_prime is False:
                    return Integer(4)
                else:
                    return Integer(2)
            elif s.is_integer:
                return S.One
            else:
                return _eps
        elif s.is_extended_negative:
            if s.is_even:
                return Integer(-2)
            elif s.is_integer:
                return S.NegativeOne
            else:
                return -_eps
        if s.is_zero or s.is_extended_nonpositive or s.is_extended_nonnegative:
            return S.Zero
        return None

    # univariate polynomial
    free = self.free_symbols
    if len(free) == 1:
        if self.is_polynomial():
            from sympy.polys.polytools import real_roots
            from sympy.polys.polyroots import roots
            from sympy.polys.polyerrors import PolynomialError
            x = free.pop()
            x0 = _monotonic_sign(x)
            if x0 in (_eps, -_eps):
                x0 = S.Zero
            if x0 is not None:
                d = self.diff(x)
                if d.is_number:
                    currentroots = []
                else:
                    try:
                        currentroots = real_roots(d)
                    except (PolynomialError, NotImplementedError):
                        currentroots = [r for r in roots(d, x) if r.is_extended_real]
                y = self.subs(x, x0)
                if x.is_nonnegative and all(
                        (r - x0).is_nonpositive for r in currentroots):
                    if y.is_nonnegative and d.is_positive:
                        if y:
                            return y if y.is_positive else Dummy('pos', positive=True)
                        else:
                            return Dummy('nneg', nonnegative=True)
                    if y.is_nonpositive and d.is_negative:
                        if y:
                            return y if y.is_negative else Dummy('neg', negative=True)
                        else:
                            return Dummy('npos', nonpositive=True)
                elif x.is_nonpositive and all(
                        (r - x0).is_nonnegative for r in currentroots):
                    if y.is_nonnegative and d.is_negative:
                        if y:
                            return Dummy('pos', positive=True)
                        else:
                            return Dummy('nneg', nonnegative=True)
                    if y.is_nonpositive and d.is_positive:
                        if y:
                            return Dummy('neg', negative=True)
                        else:
                            return Dummy('npos', nonpositive=True)
        else:
            n, d = self.as_numer_denom()
            den = None
            if n.is_number:
                den = _monotonic_sign(d)
            elif not d.is_number:
                if _monotonic_sign(n) is not None:
                    den = _monotonic_sign(d)
            if den is not None and (den.is_positive or den.is_negative):
                v = n*den
                if v.is_positive:
                    return Dummy('pos', positive=True)
                elif v.is_nonnegative:
                    return Dummy('nneg', nonnegative=True)
                elif v.is_negative:
                    return Dummy('neg', negative=True)
                elif v.is_nonpositive:
                    return Dummy('npos', nonpositive=True)
        return None

    # multivariate
    c, a = self.as_coeff_Add()
    v = None
    if not a.is_polynomial():
        # F/A or A/F where A is a number and F is a signed, rational monomial
        n, d = a.as_numer_denom()
        if not (n.is_number or d.is_number):
            return
        if (
                a.is_Mul or a.is_Pow) and \
                a.is_rational and \
                all(p.exp.is_Integer for p in a.atoms(Pow) if p.is_Pow) and \
                (a.is_positive or a.is_negative):
            v = S.One
            for ai in Mul.make_args(a):
                if ai.is_number:
                    v *= ai
                    continue
                reps = {}
                for x in ai.free_symbols:
                    reps[x] = _monotonic_sign(x)
                    if reps[x] is None:
                        return
                v *= ai.subs(reps)
    elif c:
        # signed linear expression
        if not any(p for p in a.atoms(Pow) if not p.is_number) and (a.is_nonpositive or a.is_nonnegative):
            free = list(a.free_symbols)
            p = {}
            for i in free:
                v = _monotonic_sign(i)
                if v is None:
                    return
                p[i] = v or (_eps if i.is_nonnegative else -_eps)
            v = a.xreplace(p)
    if v is not None:
        rv = v + c
        if v.is_nonnegative and rv.is_positive:
            return rv.subs(_eps, 0)
        if v.is_nonpositive and rv.is_negative:
            return rv.subs(_eps, 0)


def decompose_power(expr: Expr) -> tTuple[Expr, int]:
    """
    Decompose power into symbolic base and integer exponent.

    Examples
    ========

    >>> from sympy.core.exprtools import decompose_power
    >>> from sympy.abc import x, y
    >>> from sympy import exp

    >>> decompose_power(x)
    (x, 1)
    >>> decompose_power(x**2)
    (x, 2)
    >>> decompose_power(exp(2*y/3))
    (exp(y/3), 2)

    """
    base, exp = expr.as_base_exp()

    if exp.is_Number:
        if exp.is_Rational:
            if not exp.is_Integer:
                base = Pow(base, Rational(1, exp.q))  # type: ignore
            e = exp.p  # type: ignore
        else:
            base, e = expr, 1
    else:
        exp, tail = exp.as_coeff_Mul(rational=True)

        if exp is S.NegativeOne:
            base, e = Pow(base, tail), -1
        elif exp is not S.One:
            # todo: after dropping python 3.7 support, use overload and Literal
            #  in as_coeff_Mul to make exp Rational, and remove these 2 ignores
            tail = _keep_coeff(Rational(1, exp.q), tail)  # type: ignore
            base, e = Pow(base, tail), exp.p  # type: ignore
        else:
            base, e = expr, 1

    return base, e


def decompose_power_rat(expr: Expr) -> tTuple[Expr, Rational]:
    """
    Decompose power into symbolic base and rational exponent;
    if the exponent is not a Rational, then separate only the
    integer coefficient.

    Examples
    ========

    >>> from sympy.core.exprtools import decompose_power_rat
    >>> from sympy.abc import x
    >>> from sympy import sqrt, exp

    >>> decompose_power_rat(sqrt(x))
    (x, 1/2)
    >>> decompose_power_rat(exp(-3*x/2))
    (exp(x/2), -3)

    """
    _ = base, exp = expr.as_base_exp()
    return _ if exp.is_Rational else decompose_power(expr)


class Factors:
    """Efficient representation of ``f_1*f_2*...*f_n``."""

    __slots__ = ('factors', 'gens')

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
        frozenset({2, x})
        >>> Factors(0)
        Factors({0: 1})
        >>> Factors(I)
        Factors({I: 1})

        Notes
        =====

        Although a dictionary can be passed, only minimal checking is
        performed: powers of -1 and I are made canonical.

        """
        if isinstance(factors, (SYMPY_INTS, float)):
            factors = S(factors)
        if isinstance(factors, Factors):
            factors = factors.factors.copy()
        elif factors in (None, S.One):
            factors = {}
        elif factors is S.Zero or factors == 0:
            factors = {S.Zero: S.One}
        elif isinstance(factors, Number):
            n = factors
            factors = {}
            if n < 0:
                factors[S.NegativeOne] = S.One
                n = -n
            if n is not S.One:
                if n.is_Float or n.is_Integer or n is S.Infinity:
                    factors[n] = S.One
                elif n.is_Rational:
                    # since we're processing Numbers, the denominator is
                    # stored with a negative exponent; all other factors
                    # are left .
                    if n.p != 1:
                        factors[Integer(n.p)] = S.One
                    factors[Integer(n.q)] = S.NegativeOne
                else:
                    raise ValueError('Expected Float|Rational|Integer, not %s' % n)
        elif isinstance(factors, Basic) and not factors.args:
            factors = {factors: S.One}
        elif isinstance(factors, Expr):
            c, nc = factors.args_cnc()
            i = c.count(I)
            for _ in range(i):
                c.remove(I)
            factors = dict(Mul._from_args(c).as_powers_dict())
            # Handle all rational Coefficients
            for f in list(factors.keys()):
                if isinstance(f, Rational) and not isinstance(f, Integer):
                    p, q = Integer(f.p), Integer(f.q)
                    factors[p] = (factors[p] if p in factors else S.Zero) + factors[f]
                    factors[q] = (factors[q] if q in factors else S.Zero) - factors[f]
                    factors.pop(f)
            if i:
                factors[I] = factors.get(I, S.Zero) + i
            if nc:
                factors[Mul(*nc, evaluate=False)] = S.One
        else:
            factors = factors.copy()  # /!\ should be dict-like

            # tidy up -/+1 and I exponents if Rational

            handle = [k for k in factors if k is I or k in (-1, 1)]
            if handle:
                i1 = S.One
                for k in handle:
                    if not _isnumber(factors[k]):
                        continue
                    i1 *= k**factors.pop(k)
                if i1 is not S.One:
                    for a in i1.args if i1.is_Mul else [i1]:  # at worst, -1.0*I*(-1)**e
                        if a is S.NegativeOne:
                            factors[a] = S.One
                        elif a is I:
                            factors[I] = S.One
                        elif a.is_Pow:
                            factors[a.base] = factors.get(a.base, S.Zero) + a.exp
                        elif equal_valued(a, 1):
                            factors[a] = S.One
                        elif equal_valued(a, -1):
                            factors[-a] = S.One
                            factors[S.NegativeOne] = S.One
                        else:
                            raise ValueError('unexpected factor in i1: %s' % a)

        self.factors = factors
        keys = getattr(factors, 'keys', None)
        if keys is None:
            raise TypeError('expecting Expr or dictionary')
        self.gens = frozenset(keys())

    def __hash__(self):  # Factors
        keys = tuple(ordered(self.factors.keys()))
        values = [self.factors[k] for k in keys]
        return hash((keys, values))

    def __repr__(self):  # Factors
        return "Factors({%s})" % ', '.join(
            ['%s: %s' % (k, v) for k, v in ordered(self.factors.items())])

    @property
    def is_zero(self):  # Factors
        """
        >>> from sympy.core.exprtools import Factors
        >>> Factors(0).is_zero
        True
        """
        f = self.factors
        return len(f) == 1 and S.Zero in f

    @property
    def is_one(self):  # Factors
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
        >>> from sympy.abc import x, y
        >>> Factors((x*y**2).as_powers_dict()).as_expr()
        x*y**2

        """

        args = []
        for factor, exp in self.factors.items():
            if exp != 1:
                if isinstance(exp, Integer):
                    b, e = factor.as_base_exp()
                    e = _keep_coeff(exp, e)
                    args.append(b**e)
                else:
                    args.append(factor**exp)
            else:
                args.append(factor)
        return Mul(*args)

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

        for factor, exp in other.factors.items():
            if factor in factors:
                exp = factors[factor] + exp

                if not exp:
                    del factors[factor]
                    continue

            factors[factor] = exp

        return Factors(factors)

    def normal(self, other):
        """Return ``self`` and ``other`` with ``gcd`` removed from each.
        The only differences between this and method ``div`` is that this
        is 1) optimized for the case when there are few factors in common and
        2) this does not raise an error if ``other`` is zero.

        See Also
        ========
        div

        """
        if not isinstance(other, Factors):
            other = Factors(other)
            if other.is_zero:
                return (Factors(), Factors(S.Zero))
            if self.is_zero:
                return (Factors(S.Zero), Factors())

        self_factors = dict(self.factors)
        other_factors = dict(other.factors)

        for factor, self_exp in self.factors.items():
            try:
                other_exp = other.factors[factor]
            except KeyError:
                continue

            exp = self_exp - other_exp

            if not exp:
                del self_factors[factor]
                del other_factors[factor]
            elif _isnumber(exp):
                if exp > 0:
                    self_factors[factor] = exp
                    del other_factors[factor]
                else:
                    del self_factors[factor]
                    other_factors[factor] = -exp
            else:
                r = self_exp.extract_additively(other_exp)
                if r is not None:
                    if r:
                        self_factors[factor] = r
                        del other_factors[factor]
                    else:  # should be handled already
                        del self_factors[factor]
                        del other_factors[factor]
                else:
                    sc, sa = self_exp.as_coeff_Add()
                    if sc:
                        oc, oa = other_exp.as_coeff_Add()
                        diff = sc - oc
                        if diff > 0:
                            self_factors[factor] -= oc
                            other_exp = oa
                        elif diff < 0:
                            self_factors[factor] -= sc
                            other_factors[factor] -= sc
                            other_exp = oa - diff
                        else:
                            self_factors[factor] = sa
                            other_exp = oa
                    if other_exp:
                        other_factors[factor] = other_exp
                    else:
                        del other_factors[factor]

        return Factors(self_factors), Factors(other_factors)

    def div(self, other):  # Factors
        """Return ``self`` and ``other`` with ``gcd`` removed from each.
        This is optimized for the case when there are many factors in common.

        Examples
        ========

        >>> from sympy.core.exprtools import Factors
        >>> from sympy.abc import x, y, z
        >>> from sympy import S

        >>> a = Factors((x*y**2).as_powers_dict())
        >>> a.div(a)
        (Factors({}), Factors({}))
        >>> a.div(x*z)
        (Factors({y: 2}), Factors({z: 1}))

        The ``/`` operator only gives ``quo``:

        >>> a/x
        Factors({y: 2})

        Factors treats its factors as though they are all in the numerator, so
        if you violate this assumption the results will be correct but will
        not strictly correspond to the numerator and denominator of the ratio:

        >>> a.div(x/z)
        (Factors({y: 2}), Factors({z: -1}))

        Factors is also naive about bases: it does not attempt any denesting
        of Rational-base terms, for example the following does not become
        2**(2*x)/2.

        >>> Factors(2**(2*x + 2)).div(S(8))
        (Factors({2: 2*x + 2}), Factors({8: 1}))

        factor_terms can clean up such Rational-bases powers:

        >>> from sympy import factor_terms
        >>> n, d = Factors(2**(2*x + 2)).div(S(8))
        >>> n.as_expr()/d.as_expr()
        2**(2*x + 2)/8
        >>> factor_terms(_)
        2**(2*x)/2

        """
        quo, rem = dict(self.factors), {}

        if not isinstance(other, Factors):
            other = Factors(other)
            if other.is_zero:
                raise ZeroDivisionError
            if self.is_zero:
                return (Factors(S.Zero), Factors())

        for factor, exp in other.factors.items():
            if factor in quo:
                d = quo[factor] - exp
                if _isnumber(d):
                    if d <= 0:
                        del quo[factor]

                    if d >= 0:
                        if d:
                            quo[factor] = d

                        continue

                    exp = -d

                else:
                    r = quo[factor].extract_additively(exp)
                    if r is not None:
                        if r:
                            quo[factor] = r
                        else:  # should be handled already
                            del quo[factor]
                    else:
                        other_exp = exp
                        sc, sa = quo[factor].as_coeff_Add()
                        if sc:
                            oc, oa = other_exp.as_coeff_Add()
                            diff = sc - oc
                            if diff > 0:
                                quo[factor] -= oc
                                other_exp = oa
                            elif diff < 0:
                                quo[factor] -= sc
                                other_exp = oa - diff
                            else:
                                quo[factor] = sa
                                other_exp = oa
                        if other_exp:
                            rem[factor] = other_exp
                        else:
                            assert factor not in rem
                    continue

            rem[factor] = exp

        return Factors(quo), Factors(rem)

    def quo(self, other):  # Factors
        """Return numerator Factor of ``self / other``.

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
        """Return denominator Factors of ``self / other``.

        Examples
        ========

        >>> from sympy.core.exprtools import Factors
        >>> from sympy.abc import x, y, z
        >>> a = Factors((x*y**2).as_powers_dict())
        >>> b = Factors((x*y/z).as_powers_dict())
        >>> a.rem(b)
        Factors({z: -1})
        >>> a.rem(a)
        Factors({})
        """
        return self.div(other)[1]

    def pow(self, other):  # Factors
        """Return self raised to a non-negative integer power.

        Examples
        ========

        >>> from sympy.core.exprtools import Factors
        >>> from sympy.abc import x, y
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
                for factor, exp in self.factors.items():
                    factors[factor] = exp*other

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

        for factor, exp in self.factors.items():
            factor, exp = sympify(factor), sympify(exp)
            if factor in other.factors:
                lt = (exp - other.factors[factor]).is_negative
                if lt == True:
                    factors[factor] = exp
                elif lt == False:
                    factors[factor] = other.factors[factor]

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

        for factor, exp in other.factors.items():
            if factor in factors:
                exp = max(exp, factors[factor])

            factors[factor] = exp

        return Factors(factors)

    def __mul__(self, other):  # Factors
        return self.mul(other)

    def __divmod__(self, other):  # Factors
        return self.div(other)

    def __truediv__(self, other):  # Factors
        return self.quo(other)

    def __mod__(self, other):  # Factors
        return self.rem(other)

    def __pow__(self, other):  # Factors
        return self.pow(other)

    def __eq__(self, other):  # Factors
        if not isinstance(other, Factors):
            other = Factors(other)
        return self.factors == other.factors

    def __ne__(self, other):  # Factors
        return not self == other


class Term:
    """Efficient representation of ``coeff*(numer/denom)``. """

    __slots__ = ('coeff', 'numer', 'denom')

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

    def __truediv__(self, other):  # Term
        if isinstance(other, Term):
            return self.quo(other)
        else:
            return NotImplemented

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
        return not self == other


def _gcd_terms(terms, isprimitive=False, fraction=True):
    """Helper function for :func:`gcd_terms`.

    Parameters
    ==========

    isprimitive : boolean, optional
        If ``isprimitive`` is True then the call to primitive
        for an Add will be skipped. This is useful when the
        content has already been extracted.

    fraction : boolean, optional
        If ``fraction`` is True then the expression will appear over a common
        denominator, the lcm of all term denominators.
    """

    if isinstance(terms, Basic) and not isinstance(terms, Tuple):
        terms = Add.make_args(terms)

    terms = list(map(Term, [t for t in terms if t]))

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
            denom = Term(S.One).numer

        cont = cont.as_expr()
        numer = Add(*numers)
        denom = denom.as_expr()

    if not isprimitive and numer.is_Add:
        _cont, numer = numer.primitive()
        cont *= _cont

    return cont, numer, denom


def gcd_terms(terms, isprimitive=False, clear=True, fraction=True):
    """Compute the GCD of ``terms`` and put them together.

    Parameters
    ==========

    terms : Expr
        Can be an expression or a non-Basic sequence of expressions
        which will be handled as though they are terms from a sum.

    isprimitive : bool, optional
        If ``isprimitive`` is True the _gcd_terms will not run the primitive
        method on the terms.

    clear : bool, optional
        It controls the removal of integers from the denominator of an Add
        expression. When True (default), all numerical denominator will be cleared;
        when False the denominators will be cleared only if all terms had numerical
        denominators other than 1.

    fraction : bool, optional
        When True (default), will put the expression over a common
        denominator.

    Examples
    ========

    >>> from sympy import gcd_terms
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
    >>> gcd_terms(x/2/y + 1/x/y, clear=False)
    (x**2/2 + 1)/(x*y)
    >>> gcd_terms(x/2/y + 1/x/y, clear=False, fraction=False)
    (x/2 + 1/x)/y

    The ``clear`` flag was ignored in this case because the returned
    expression was a rational expression, not a simple sum.

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
                nc = Mul(*nc)
                d = Dummy()
                reps.append((d, nc))
                c.append(d)
                args[i] = Mul(*c)
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
        if not clear:
            c, _coeff = coeff.as_coeff_Mul()
            if not c.is_Integer and not clear and numer.is_Add:
                n, d = c.as_numer_denom()
                _numer = numer/d
                if any(a.as_coeff_Mul()[0].is_Integer
                        for a in _numer.args):
                    numer = _numer
                    coeff = n*_coeff
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
                if not a.args:
                    return a
                return a.func(*[handle(i) for i in a.args])
            return type(a)([handle(i) for i in a])
        return gcd_terms(a, isprimitive, clear, fraction)

    if isinstance(terms, Dict):
        return Dict(*[(k, handle(v)) for k, v in terms.args])
    return terms.func(*[handle(i) for i in terms.args])


def _factor_sum_int(expr, **kwargs):
    """Return Sum or Integral object with factors that are not
    in the wrt variables removed. In cases where there are additive
    terms in the function of the object that are independent, the
    object will be separated into two objects.

    Examples
    ========

    >>> from sympy import Sum, factor_terms
    >>> from sympy.abc import x, y
    >>> factor_terms(Sum(x + y, (x, 1, 3)))
    y*Sum(1, (x, 1, 3)) + Sum(x, (x, 1, 3))
    >>> factor_terms(Sum(x*y, (x, 1, 3)))
    y*Sum(x, (x, 1, 3))

    Notes
    =====

    If a function in the summand or integrand is replaced
    with a symbol, then this simplification should not be
    done or else an incorrect result will be obtained when
    the symbol is replaced with an expression that depends
    on the variables of summation/integration:

    >>> eq = Sum(y, (x, 1, 3))
    >>> factor_terms(eq).subs(y, x).doit()
    3*x
    >>> eq.subs(y, x).doit()
    6
    """
    result = expr.function
    if result == 0:
        return S.Zero
    limits = expr.limits

    # get the wrt variables
    wrt = {i.args[0] for i in limits}

    # factor out any common terms that are independent of wrt
    f = factor_terms(result, **kwargs)
    i, d = f.as_independent(*wrt)
    if isinstance(f, Add):
        return i * expr.func(1, *limits) + expr.func(d, *limits)
    else:
        return i * expr.func(d, *limits)


def factor_terms(expr, radical=False, clear=False, fraction=False, sign=True):
    """Remove common factors from terms in all arguments without
    changing the underlying structure of the expr. No expansion or
    simplification (and no processing of non-commutatives) is performed.

    Parameters
    ==========

    radical: bool, optional
        If radical=True then a radical common to all terms will be factored
        out of any Add sub-expressions of the expr.

    clear : bool, optional
        If clear=False (default) then coefficients will not be separated
        from a single Add if they can be distributed to leave one or more
        terms with integer coefficients.

    fraction : bool, optional
        If fraction=True (default is False) then a common denominator will be
        constructed for the expression.

    sign : bool, optional
        If sign=True (default) then even if the only factor in common is a -1,
        it will be factored out of the expression.

    Examples
    ========

    >>> from sympy import factor_terms, Symbol
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
    def do(expr):
        from sympy.concrete.summations import Sum
        from sympy.integrals.integrals import Integral
        is_iterable = iterable(expr)

        if not isinstance(expr, Basic) or expr.is_Atom:
            if is_iterable:
                return type(expr)([do(i) for i in expr])
            return expr

        if expr.is_Pow or expr.is_Function or \
                is_iterable or not hasattr(expr, 'args_cnc'):
            args = expr.args
            newargs = tuple([do(i) for i in args])
            if newargs == args:
                return expr
            return expr.func(*newargs)

        if isinstance(expr, (Sum, Integral)):
            return _factor_sum_int(expr,
                radical=radical, clear=clear,
                fraction=fraction, sign=sign)

        cont, p = expr.as_content_primitive(radical=radical, clear=clear)
        if p.is_Add:
            list_args = [do(a) for a in Add.make_args(p)]
            # get a common negative (if there) which gcd_terms does not remove
            if not any(a.as_coeff_Mul()[0].extract_multiplicatively(-1) is None
                       for a in list_args):
                cont = -cont
                list_args = [-a for a in list_args]
            # watch out for exp(-(x+2)) which gcd_terms will change to exp(-x-2)
            special = {}
            for i, a in enumerate(list_args):
                b, e = a.as_base_exp()
                if e.is_Mul and e != Mul(*e.args):
                    list_args[i] = Dummy()
                    special[list_args[i]] = a
            # rebuild p not worrying about the order which gcd_terms will fix
            p = Add._from_args(list_args)
            p = gcd_terms(p,
                isprimitive=True,
                clear=clear,
                fraction=fraction).xreplace(special)
        elif p.args:
            p = p.func(
                *[do(a) for a in p.args])
        rv = _keep_coeff(cont, p, clear=clear, sign=sign)
        return rv
    expr = sympify(expr)
    return do(expr)


def factor_nc(expr):
    """Return the factored form of ``expr`` while handling non-commutative
    expressions.

    Examples
    ========

    >>> from sympy import factor_nc, Symbol
    >>> from sympy.abc import x
    >>> A = Symbol('A', commutative=False)
    >>> B = Symbol('B', commutative=False)
    >>> factor_nc((x**2 + 2*A*x + A**2).expand())
    (x + A)**2
    >>> factor_nc(((x + A)*(x + B)).expand())
    (x + A)*(x + B)
    """
    expr = sympify(expr)
    if not isinstance(expr, Expr) or not expr.args:
        return expr
    if not expr.is_Add:
        return expr.func(*[factor_nc(a) for a in expr.args])
    expr = expr.func(*[expand_power_exp(i) for i in expr.args])

    return _NCPoly(expr).factor_one()


def _nc_generators(expr):
    """Return noncommutative polynomial generators of an expression.

    Commutative symbols are excluded, and functions of generators other than
    positive integer powers are treated as generators, as with the polynomials module
    """
    if expr.is_Add or expr.is_Mul:
        generators = set()
        for arg in expr.args:
            generators |= _nc_generators(arg)
        return generators
    if expr.is_commutative:
        return set()
    if expr.is_Atom:
        return set([expr])
    if expr.is_Pow:
        term, pw = expr.args
        if pw.is_positive and pw.is_Integer:
            return _nc_generators(term)
    return set([expr])

def _nc_degree(expr, generators):
    """Return the degree of an expression in the given generators"""
    if expr.is_Add:
        return max([_nc_degree(arg, generators) for arg in expr.args])
    if expr.is_Mul:
        return sum([_nc_degree(arg, generators) for arg in expr.args])
    if expr.is_commutative:
        return 0
    if expr in generators:
        return 1
    if expr.is_Pow:
        term, pw = expr.args
        if term not in generators:
            raise ValueError(f"_nc_degree: Expected {term} to be in generators")
        if pw.is_positive and pw.is_Integer:
            return pw
    raise ValueError(f"_nc_degree: Expected {expr} to be in generators")

class _NCMonomial():
    """Class that represents a monomial with noncommutative generators"""
    def __init__(self, coeff, term, generators, degree=None):
        self.coeff = coeff
        self.term = term
        self.generators = generators
        self.degree = degree
        if self.degree is None:
            self.degree = _nc_degree(self.term, self.generators)

    def __add__(self, other):
        if not isinstance(other, _NCMonomial):
            return NotImplemented
        if other.term != self.term:
            raise ValueError(f"Cannot add {self} and {other} as monomials")
        if other.generators != self.generators:
            raise ValueError(f"Cannot add {self} and {other} with different generators")
        return _NCMonomial(self.coeff + other.coeff, self.term, self.generators, self.degree)

    def __mul__(self, other):
        if isinstance(other, _NCMonomial):
            if other.generators != self.generators:
                raise ValueError(f"Cannot multiply {self} and {other} with different generators")
            return _NCMonomial(self.coeff * other.coeff, self.term * other.term, self.generators, self.degree + other.degree)
        if isinstance(other, _NCPoly):
            return NotImplemented
        return _NCMonomial((self.coeff * other).cancel(), self.term, self.generators, self.degree)

    def __str__(self):
        return str(self.coeff * self.term)

    def __repr__(self):
        return str(self)
            
    def as_expr(self):
        """Convert to Expr"""
        return self.coeff * self.term

    def split(self, ldegree):
        """Split into two monomials of ldegree and self.degree - ldegree

        The left monomial will have the original coefficient.
        The right monomial will always have a coefficient of 1.
        """
        if ldegree < 0:
            raise ValueError(f"ldegree must be nonnegative")
        if ldegree == 0:
            return (
                _NCMonomial(self.coeff, Integer(1), self.generators, 0),
                _NCMonomial(Integer(1), self.term, self.generators, self.degree),
            )
        if ldegree == self.degree:
            return (
                self,
                _NCMonomial(Integer(1), Integer(1), self.generators, 0)
            )
            
        if ldegree > self.degree:
            raise ValueError(f"Monomial {self} of degree {self.degree} can't be split with left term of degree {ldegree}")
        rdegree = self.degree - ldegree
        
        if self.term.is_Pow:
            subterms = [self.term]
        elif self.term.is_Mul:
            subterms = self.term.args

        acc = 0
        for idx, subterm in enumerate(subterms):
            if subterm.is_Pow:
                trm = subterm.args[0]
                pw = subterm.args[1]
            else:
                trm = subterm
                pw = 1
            acc += pw
            if acc < ldegree:
                continue
            gap = acc - ldegree
            lterm = Mul(*subterms[:idx], trm**(pw - gap))
            rterm = Mul(trm**gap, *subterms[idx+1:])
            return (
                _NCMonomial(self.coeff, lterm, self.generators, degree=ldegree),
                _NCMonomial(Integer(1), rterm, self.generators, degree=rdegree)
            )


    def ldivide(self, other):
        """Return the monomial m such that self == other * m

        If no such monomial exists, return None
        """
        if other.degree > self.degree:
            return None
        if other.degree == self.degree:
            if self.term == other.term:
                return _NCMonomial(self.coeff / other.coeff, Integer(1), self.generators, degree=0)
            return None
        # Could detect inequality earlier if I don't use split
        left, right = self.split(other.degree)
        if left.term == other.term:
            return _NCMonomial(self.coeff / other.coeff, right.term, self.generators, self.degree - other.degree)
        return None


    def rdivide(self, other):
        """Return the monomial m such that self == m * other

        If no such monomial exists, return None
        """
        if other.degree > self.degree:
            return None
        if other.degree == self.degree:
            if self.term == other.term:
                return _NCMonomial(self.coeff / other.coeff, Integer(1), self.generators, degree=0)
            return None
        # Could detect inequality earlier if I don't use split
        left, right = self.split(self.degree - other.degree)
        if right.term == other.term:
            return _NCMonomial(self.coeff / other.coeff, left.term, self.generators, self.degree - other.degree)
        return None



class _NCPoly():
    """Class that represents a polynomial with noncommutative generators

    It is used only for factoring expressions with noncommutative terms and does not have the
    same interface as the Poly class (which supports only commutative terms).

    A few functions (add, subtract) update an instance in-place, unlike most sympy classes."""
    def __init__(self, expr, generators=None):
        if expr is None:
            expr = {}
        if isinstance(expr, dict):
            self.rep = expr
            self.generators = generators
        else:
            expr = expr.expand()
            self.generators = generators
            if generators is None:
                self.generators = _nc_generators(expr)
            coeffs_dict = expr.as_coefficients_dict(*(self.generators or [None]))

            # Representation: Dict[degree, Dict[term, _NCMonomial]]
            self.rep = {}
            for term, coeff in coeffs_dict.items():
                degree = _nc_degree(term, self.generators)
                if degree not in self.rep:
                    self.rep[degree] = {}
                self.rep[degree][term] = _NCMonomial(coeff, term, self.generators, degree)
        self._eval_degree()

    def __add__(self, other):
        out = self.copy()
        if isinstance(other, _NCPoly):
            for monom in other.monomials():
                out.add(monom)
        elif isinstance(other, _NCMonomial):
            out.add(other)
        else:
            raise ValueError(f"Can't add {other} to _NCPoly")
        return out

    def __sub__(self, other):
        out = self.copy()
        if isinstance(other, _NCPoly):
            for monom in other.monomials():
                out.subtract(monom)
        elif isinstance(other, _NCMonomial):
            out.subtract(other)
        else:
            raise ValueError(f"Can't add {other} to _NCPoly")
        return out

    def __mul__(self, other):
        out = _NCPoly(None, self.generators)
        if isinstance(other, _NCPoly):
            if self.generators != other.generators:
                raise ValueError("Incompatible generators")
            for lterm in self.monomials():
                for rterm in other.monomials():
                    out.add(lterm * rterm)
        else:
            for term in self.monomials():
                out.add(term * other)
        return out

    def __rmul__(self, other):
        out = _NCPoly(None, self.generators)
        for term in self.monomials():
            out.add(other * term)
        return out

    def __len__(self):
        """Number of monomials"""
        if len(self.rep) == 0:
            return 0
        return sum([len(terms) for terms in self.rep.values()])

    def _eval_degree(self):
        """Update the degree field"""
        if len(self.rep) == 0:
            self.degree = 0
        else:
            self.degree = max(self.rep.keys())

    def copy(self, degree=None):
        """Return a copy, optionally with only the monomials of the specified degree"""
        if degree is None:
            new_rep = {k: v.copy() for k, v in self.rep.items()}
        else:
            new_rep = {}
            if degree in self.rep:
                new_rep[degree] = self.rep[degree].copy()
        return _NCPoly(new_rep, self.generators)

    def as_expr(self):
        """Convert to Expr"""
        return Add(*[monom.as_expr() for monom in self.monomials()])
        
    def monomials(self, degree=None):
        """Return iterator of all monomials (with specified degree)"""
        if degree is None:
            for monomials in self.rep.values():
                for monomial in monomials.values():
                    yield monomial
        else:
            for monomial in self.rep.get(degree, {}).values():
                yield monomial

    def extract_coeff_factor(self, degree=None):
        from sympy.simplify.radsimp import fraction
        from sympy.polys.polytools import gcd, lcm
        numerators = []
        denominators = []
        for monom in self.monomials(degree):
            num, denom = fraction(monom.coeff.cancel())
            numerators.append(num)
            denominators.append(denom)
        return gcd(numerators) / lcm(denominators)
 
    def add(self, monom):
        """Add a monomial to this polynomial, updating in-place"""
        degree = monom.degree
        term = monom.term
        if degree > self.degree:
            self.degree = degree

        if degree not in self.rep:
            self.rep[degree] = {}

        if term in self.rep[degree]:
            new_monom = self.rep[degree][term] + monom

            if new_monom.coeff == 0:
                del self.rep[degree][term]
                if len(self.rep[degree]) == 0:
                    del self.rep[degree]
                    if degree == self.degree:
                        self._eval_degree()
            else:
                self.rep[degree][term] = new_monom
        else:
            self.rep[degree][term] = monom

    def subtract(self, monom):
        """Subtract a monomial from this polynomial, updating in-place"""
        return self.add(monom * -1)


    def _factor_homogeneous(self, ldegree):
        """Factor the maximum-degree terms into two polynomials
        of degree ldegree and self.degree - ldegree, if possible.

        Return (None, None) if no such factorization exists.

        Implementation based on Algorithm 1 from https://arxiv.org/abs/1002.3180
        """
        rdegree = self.degree - ldegree
        terms = list(self.monomials(self.degree))
        Ghat, Hhat = terms[0].split(ldegree)
        G = _NCPoly(None, self.generators)
        H = _NCPoly(None, self.generators)
        G.add(Ghat)
        H.add(Hhat)
        for monomial in terms[1:]:
            R = monomial.ldivide(Ghat)
            if R is None:
                L = monomial.rdivide(Hhat)
                if L is not None:
                    G.add(L)
            else:
                H.add(R)
        if len(G) * len(H) != len(terms):
            return None, None
        GH = G * H
        if GH.as_expr() == Add(*[term.as_expr() for term in terms]):
            return G, H
        return None, None


    def _factor_once(self, ldegree):
        """Factor into pairs of polynomials of degree ldegree and self.degree - ldegree, if possible.
        Coefficients will not be factored or otherwise simplified.

        Return a generator that yields tuples of polynomials.

        Implementation based on Algorithm 2 from https://arxiv.org/abs/1002.3180
        """
        from sympy.solvers.polysys import solve_poly_system
        if self.degree == 0:
            return
        if len(self) < 2:
            return

        if ldegree == 0 or ldegree >= self.degree:
            return
            
        G, H = self._factor_homogeneous(ldegree)
        if G is None:
            return

        n = self.degree
        h = ldegree
        k = n - ldegree
        # Future optimization: select the monomials with the least overlap potential
        Ghat = next(G.monomials())
        Hhat = next(H.monomials())
        
        dummies = []
        Fhat = self.copy()

        for j in range(1, n + 1):
            for i in range(1, j):
                Fhat -= G.copy(h - i) * H.copy(k - j + i)

            if j <= min(h, k):
                Ghat_L, Ghat_R = Ghat.split(h - j)
                Hhat_L, Hhat_R = Hhat.split(j)
                if Ghat_R.term == Hhat_L.term:
                    # Lines 7-8 of the published algorithm suggest that the term should be
                    # Ghat_L.term * Hhat_R.term. This doesn't make sense, however - the term is
                    # supposed to have degree n - j, not n - 2j.
                    term = Ghat.term * Hhat_R.term
                    c = 0
                    if term in Fhat.rep.get(n - j, {}):
                        c = Fhat.rep[n - j][term].coeff
                        Fhat.subtract(_NCMonomial(c, term, self.generators, n - j))
                    dummy = Dummy()
                    dummies.append(dummy)
                    # Lines 10-11 of the published algorithm suggest that we should add a monomial
                    # with term Ghat_L * Hhat_R to both G and H. If this were the case, then in
                    # step 3 of the example in 3.3 from the paper, we would add -alpha * y to G instead of -alpha.
                    # This seems to be a second minor mistake in the published algorithm.
                    #
                    # There is also a subtlety on handling coefficients not mentioned in the algorithm.
                    # We have to divide the new H term by the coefficient of Ghat
                    # (We would have to do the same to the new G term, except for the invariant
                    # that Hhat never has a coefficient)
                    G.add(_NCMonomial(dummy, Ghat_L.term, self.generators, h - j))
                    H.add(_NCMonomial((c - dummy)/Ghat.coeff, Hhat_R.term, self.generators, k - j))

            for monom in Fhat.monomials(n - j):
                R = monom.ldivide(Ghat)
                if R is not None:
                    H.add(R)
                L = monom.rdivide(Hhat)
                if L is not None:
                    G.add(L)

        diff = self - G * H
        if len(dummies) == 0:
            if len(diff) == 0:
                yield G, H
            return
        eqs = [term.coeff for term in diff.monomials()]
        solns = solve_poly_system([eqs], *dummies)
        if len(solns) == 0:
            return
        Gexpr = G.as_expr()
        Hexpr = H.as_expr()
        for soln in solns:
            for sub in soln:
                if not sub.is_rational:
                    import pdb; pdb.set_trace()
            substitutes = list(zip(dummies, soln))
            yield (
                _NCPoly(Gexpr.subs(substitutes), self.generators),
                _NCPoly(Hexpr.subs(substitutes), self.generators),
            )
                    
    def _factor_recursive(self):
        """Return a generator that yields all unique, valid, irreducible factorizations.
        Does not fully simplify coefficients.
        """
        self_expr = self.as_expr()
        out = set()
        for degree in range(1, self.degree):
            for G, H in self._factor_once(degree):
                Gcoeff, G = G.factor_coeff()
                Hcoeff, H = H.factor_coeff()
                coeff = Gcoeff * Hcoeff
                Gs = G._factor_recursive()
                Hs = H._factor_recursive()
                H_list = []
                G_first = next(Gs)
                coeff_G = coeff * G_first
                for H_expr in Hs:
                    H_list.append(H_expr)
                    expr = coeff_G * H_expr
                    if expr not in out:
                        yield expr
                        out.add(expr)
                for G_expr in Gs:
                    coeff_G = coeff * G_expr
                    for H_expr in H_list:
                        expr = coeff_G * H_expr
                        if expr not in out:
                            yield expr
                            out.add(expr)
        if len(out) == 0:
            yield self_expr

    def factor_one(self):
        """Return a valid factorization. Also handles factorization of commutative coefficients.

        Factorizations of noncommutative polynomials are not unique in general. This method is
        deterministic, and it is guaranteed to return a valid, irreducible factorization. There are
        no further guarantees as to *which* valid, irreducible factorization is returned.
        """
        coeff_factor, F = self.factor_coeff()
        out = next(F._factor_recursive())
        if coeff_factor != 1:
            out = coeff_factor * out
        return out

    def factor_all(self):
        """Return a generator that yields all unique, valid, irreducible factorizations.
        Also handles factorizations of commutative coefficients.
        """
        coeff_factor, F = self.factor_coeff()
        factorizations = F._factor_recursive()
        if coeff_factor != 1:
            factorizations = (coeff_factor * f for f in factorizations)
        return factorizations

    def factor_coeff(self):
        coeff_factor = self.extract_coeff_factor()
        if coeff_factor == 1:
            return (Integer(1), self)
        return (coeff_factor.factor(), self * (1/coeff_factor))


    def __str__(self):
        return str(self.as_expr())

    def __repr__(self):
        return str(self)
