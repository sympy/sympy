from __future__ import print_function, division

from math import log as _log

from .sympify import _sympify
from .cache import cacheit
from .singleton import S
from .expr import Expr
from .evalf import PrecisionExhausted
from .function import (_coeff_isneg, expand_complex, expand_multinomial,
    expand_mul)
from .logic import fuzzy_bool
from .compatibility import as_int, range
from .evaluate import global_evaluate

from mpmath.libmp import sqrtrem as mpmath_sqrtrem
from sympy.utilities.iterables import sift


def integer_nthroot(y, n):
    """
    Return a tuple containing x = floor(y**(1/n))
    and a boolean indicating whether the result is exact (that is,
    whether x**n == y).

    >>> from sympy import integer_nthroot
    >>> integer_nthroot(16,2)
    (4, True)
    >>> integer_nthroot(26,2)
    (5, False)

    """
    y, n = int(y), int(n)
    if y < 0:
        raise ValueError("y must be nonnegative")
    if n < 1:
        raise ValueError("n must be positive")
    if y in (0, 1):
        return y, True
    if n == 1:
        return y, True
    if n == 2:
        x, rem = mpmath_sqrtrem(y)
        return int(x), not rem
    if n > y:
        return 1, False
    # Get initial estimate for Newton's method. Care must be taken to
    # avoid overflow
    try:
        guess = int(y**(1./n) + 0.5)
    except OverflowError:
        exp = _log(y, 2)/n
        if exp > 53:
            shift = int(exp - 53)
            guess = int(2.0**(exp - shift) + 1) << shift
        else:
            guess = int(2.0**exp)
    if guess > 2**50:
        # Newton iteration
        xprev, x = -1, guess
        while 1:
            t = x**(n - 1)
            xprev, x = x, ((n - 1)*x + y//t)//n
            if abs(x - xprev) < 2:
                break
    else:
        x = guess
    # Compensate
    t = x**n
    while t < y:
        x += 1
        t = x**n
    while t > y:
        x -= 1
        t = x**n
    return x, t == y


class Pow(Expr):
    """
    Defines the expression x**y as "x raised to a power y"

    Singleton definitions involving (0, 1, -1, oo, -oo):

    +--------------+---------+-----------------------------------------------+
    | expr         | value   | reason                                        |
    +==============+=========+===============================================+
    | z**0         | 1       | Although arguments over 0**0 exist, see [2].  |
    +--------------+---------+-----------------------------------------------+
    | z**1         | z       |                                               |
    +--------------+---------+-----------------------------------------------+
    | (-oo)**(-1)  | 0       |                                               |
    +--------------+---------+-----------------------------------------------+
    | (-1)**-1     | -1      |                                               |
    +--------------+---------+-----------------------------------------------+
    | S.Zero**-1   | zoo     | This is not strictly true, as 0**-1 may be    |
    |              |         | undefined, but is convenient in some contexts |
    |              |         | where the base is assumed to be positive.     |
    +--------------+---------+-----------------------------------------------+
    | 1**-1        | 1       |                                               |
    +--------------+---------+-----------------------------------------------+
    | oo**-1       | 0       |                                               |
    +--------------+---------+-----------------------------------------------+
    | 0**oo        | 0       | Because for all complex numbers z near        |
    |              |         | 0, z**oo -> 0.                                |
    +--------------+---------+-----------------------------------------------+
    | 0**-oo       | zoo     | This is not strictly true, as 0**oo may be    |
    |              |         | oscillating between positive and negative     |
    |              |         | values or rotating in the complex plane.      |
    |              |         | It is convenient, however, when the base      |
    |              |         | is positive.                                  |
    +--------------+---------+-----------------------------------------------+
    | 1**oo        | nan     | Because there are various cases where         |
    | 1**-oo       |         | lim(x(t),t)=1, lim(y(t),t)=oo (or -oo),       |
    | 1**zoo       |         | but lim( x(t)**y(t), t) != 1.  See [3].       |
    +--------------+---------+-----------------------------------------------+
    | (-1)**oo     | nan     | Because of oscillations in the limit.         |
    | (-1)**(-oo)  |         |                                               |
    +--------------+---------+-----------------------------------------------+
    | oo**oo       | oo      |                                               |
    +--------------+---------+-----------------------------------------------+
    | oo**-oo      | 0       |                                               |
    +--------------+---------+-----------------------------------------------+
    | (-oo)**oo    | nan     |                                               |
    | (-oo)**-oo   |         |                                               |
    +--------------+---------+-----------------------------------------------+

    Because symbolic computations are more flexible that floating point
    calculations and we prefer to never return an incorrect answer,
    we choose not to conform to all IEEE 754 conventions.  This helps
    us avoid extra test-case code in the calculation of limits.

    See Also
    ========

    sympy.core.numbers.Infinity
    sympy.core.numbers.NegativeInfinity
    sympy.core.numbers.NaN

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Exponentiation
    .. [2] http://en.wikipedia.org/wiki/Exponentiation#Zero_to_the_power_of_zero
    .. [3] http://en.wikipedia.org/wiki/Indeterminate_forms

    """
    is_Pow = True

    __slots__ = ['is_commutative']

    @cacheit
    def __new__(cls, b, e, evaluate=None):
        if evaluate is None:
            evaluate = global_evaluate[0]
        from sympy.functions.elementary.exponential import exp_polar

        b = _sympify(b)
        e = _sympify(e)
        if evaluate:
            if e is S.Zero:
                return S.One
            elif e is S.One:
                return b
            elif e.is_integer and _coeff_isneg(b):
                if e.is_even:
                    b = -b
                elif e.is_odd:
                    return -Pow(-b, e)
            if S.NaN in (b, e):  # XXX S.NaN**x -> S.NaN under assumption that x != 0
                return S.NaN
            elif b is S.One:
                if abs(e).is_infinite:
                    return S.NaN
                return S.One
            else:
                # recognize base as E
                if not e.is_Atom and b is not S.Exp1 and b.func is not exp_polar:
                    from sympy import numer, denom, log, sign, im, factor_terms
                    c, ex = factor_terms(e, sign=False).as_coeff_Mul()
                    den = denom(ex)
                    if den.func is log and den.args[0] == b:
                        return S.Exp1**(c*numer(ex))
                    elif den.is_Add:
                        s = sign(im(b))
                        if s.is_Number and s and den == \
                                log(-factor_terms(b, sign=False)) + s*S.ImaginaryUnit*S.Pi:
                            return S.Exp1**(c*numer(ex))

                obj = b._eval_power(e)
                if obj is not None:
                    return obj
        obj = Expr.__new__(cls, b, e)
        obj.is_commutative = (b.is_commutative and e.is_commutative)
        return obj

    @property
    def base(self):
        return self._args[0]

    @property
    def exp(self):
        return self._args[1]

    @classmethod
    def class_key(cls):
        return 3, 2, cls.__name__

    def _eval_power(self, other):
        from sympy import Abs, arg, exp, floor, im, log, re, sign
        b, e = self.as_base_exp()
        if b is S.NaN:
            return (b**e)**other  # let __new__ handle it

        s = None
        if other.is_integer:
            s = 1
        elif b.is_polar:  # e.g. exp_polar, besselj, var('p', polar=True)...
            s = 1
        elif e.is_real is not None:
            # helper functions ===========================
            def _half(e):
                """Return True if the exponent has a literal 2 as the
                denominator, else None."""
                if getattr(e, 'q', None) == 2:
                    return True
                n, d = e.as_numer_denom()
                if n.is_integer and d == 2:
                    return True
            def _n2(e):
                """Return ``e`` evaluated to a Number with 2 significant
                digits, else None."""
                try:
                    rv = e.evalf(2, strict=True)
                    if rv.is_Number:
                        return rv
                except PrecisionExhausted:
                    pass
            # ===================================================
            if e.is_real:
                # we need _half(other) with constant floor or
                # floor(S.Half - e*arg(b)/2/pi) == 0

                # handle -1 as special case
                if (e == -1) == True:
                    # floor arg. is 1/2 + arg(b)/2/pi
                    if _half(other):
                        if b.is_negative is True:
                            return S.NegativeOne**other*Pow(-b, e*other)
                        if b.is_real is False:
                            return Pow(b.conjugate()/Abs(b)**2, other)
                elif e.is_even:
                    if b.is_real:
                        b = abs(b)
                    if b.is_imaginary:
                        b = abs(im(b))*S.ImaginaryUnit

                if (abs(e) < 1) == True or (e == 1) == True:
                    s = 1  # floor = 0
                elif b.is_nonnegative:
                    s = 1  # floor = 0
                elif re(b).is_nonnegative and (abs(e) < 2) == True:
                    s = 1  # floor = 0
                elif im(b).is_nonzero and (abs(e) == 2) == True:
                    s = 1  # floor = 0
                elif _half(other):
                    s = exp(2*S.Pi*S.ImaginaryUnit*other*floor(
                        S.Half - e*arg(b)/(2*S.Pi)))
                    if s.is_real and _n2(sign(s) - s) == 0:
                        s = sign(s)
                    else:
                        s = None
            else:
                # e.is_real is False requires:
                #     _half(other) with constant floor or
                #     floor(S.Half - im(e*log(b))/2/pi) == 0
                try:
                    s = exp(2*S.ImaginaryUnit*S.Pi*other*
                        floor(S.Half - im(e*log(b))/2/S.Pi))
                    # be careful to test that s is -1 or 1 b/c sign(I) == I:
                    # so check that s is real
                    if s.is_real and _n2(sign(s) - s) == 0:
                        s = sign(s)
                    else:
                        s = None
                except PrecisionExhausted:
                    s = None

        if s is not None:
            return s*Pow(b, e*other)

    def _eval_is_even(self):
        if self.exp.is_integer and self.exp.is_positive:
            return self.base.is_even

    def _eval_is_positive(self):
        from sympy import log
        if self.base == self.exp:
            if self.base.is_nonnegative:
                return True
        elif self.base.is_positive:
            if self.exp.is_real:
                return True
        elif self.base.is_negative:
            if self.exp.is_even:
                return True
            if self.exp.is_odd:
                return False
        elif self.base.is_nonpositive:
            if self.exp.is_odd:
                return False
        elif self.base.is_imaginary:
            if self.exp.is_integer:
                m = self.exp % 4
                if m.is_zero:
                    return True
                if m.is_integer and m.is_zero is False:
                    return False
            if self.exp.is_imaginary:
                return log(self.base).is_imaginary

    def _eval_is_negative(self):
        if self.base.is_negative:
            if self.exp.is_odd:
                return True
            if self.exp.is_even:
                return False
        elif self.base.is_positive:
            if self.exp.is_real:
                return False
        elif self.base.is_nonnegative:
            if self.exp.is_nonnegative:
                return False
        elif self.base.is_nonpositive:
            if self.exp.is_even:
                return False
        elif self.base.is_real:
            if self.exp.is_even:
                return False

    def _eval_is_zero(self):
        if self.base.is_zero:
            if self.exp.is_positive:
                return True
            elif self.exp.is_nonpositive:
                return False
        elif self.base.is_nonzero:
            if self.exp.is_finite:
                return False
            elif self.exp.is_infinite:
                if (1 - abs(self.base)).is_positive:
                    return self.exp.is_positive
                elif (1 - abs(self.base)).is_negative:
                    return self.exp.is_negative

    def _eval_is_integer(self):
        b, e = self.args
        if b.is_rational:
            if b.is_integer is False and e.is_positive:
                return False  # rat**nonneg
        if b.is_integer and e.is_integer:
            if b is S.NegativeOne:
                return True
            if e.is_nonnegative or e.is_positive:
                return True
        if b.is_integer and e.is_negative and (e.is_finite or e.is_integer):
            if (b - 1).is_nonzero and (b + 1).is_nonzero:
                return False
        if b.is_Number and e.is_Number:
            check = self.func(*self.args)
            return check.is_Integer

    def _eval_is_real(self):
        from sympy import arg, exp, log, Mul
        real_b = self.base.is_real
        if real_b is None:
            if self.base.func == exp and self.base.args[0].is_imaginary:
                return self.exp.is_imaginary
            return
        real_e = self.exp.is_real
        if real_e is None:
            return
        if real_b and real_e:
            if self.base.is_positive:
                return True
            elif self.base.is_nonnegative:
                if self.exp.is_nonnegative:
                    return True
            else:
                if self.exp.is_integer:
                    return True
                elif self.base.is_negative:
                    if self.exp.is_Rational:
                        return False
        if real_e and self.exp.is_negative:
            return Pow(self.base, -self.exp).is_real
        im_b = self.base.is_imaginary
        im_e = self.exp.is_imaginary
        if im_b:
            if self.exp.is_integer:
                if self.exp.is_even:
                    return True
                elif self.exp.is_odd:
                    return False
            elif im_e and log(self.base).is_imaginary:
                return True
            elif self.exp.is_Add:
                c, a = self.exp.as_coeff_Add()
                if c and c.is_Integer:
                    return Mul(
                        self.base**c, self.base**a, evaluate=False).is_real
            elif self.base in (-S.ImaginaryUnit, S.ImaginaryUnit):
                if (self.exp/2).is_integer is False:
                    return False
        if real_b and im_e:
            if self.base is S.NegativeOne:
                return True
            c = self.exp.coeff(S.ImaginaryUnit)
            if c:
                ok = (c*log(self.base)/S.Pi).is_Integer
                if ok is not None:
                    return ok

        if real_b is False:  # we already know it's not imag
            i = arg(self.base)*self.exp/S.Pi
            return i.is_integer

    def _eval_is_complex(self):
        if all(a.is_complex for a in self.args):
            return True

    def _eval_is_imaginary(self):
        from sympy import arg, log
        if self.base.is_imaginary:
            if self.exp.is_integer:
                odd = self.exp.is_odd
                if odd is not None:
                    return odd
                return

        if self.exp.is_imaginary:
            imlog = log(self.base).is_imaginary
            if imlog is not None:
                return False  # I**i -> real; (2*I)**i -> complex ==> not imaginary

        if self.base.is_real and self.exp.is_real:
            if self.base.is_positive:
                return False
            else:
                rat = self.exp.is_rational
                if not rat:
                    return rat
                if self.exp.is_integer:
                    return False
                else:
                    half = (2*self.exp).is_integer
                    if half:
                        return self.base.is_negative
                    return half

        if self.base.is_real is False:  # we already know it's not imag
            i = arg(self.base)*self.exp/S.Pi
            return (2*i).is_odd

    def _eval_is_odd(self):
        if self.exp.is_integer:
            if self.exp.is_positive:
                return self.base.is_odd
            elif self.exp.is_nonnegative and self.base.is_odd:
                return True
            elif self.base is S.NegativeOne:
                return True

    def _eval_is_finite(self):
        if self.exp.is_negative:
            if self.base.is_zero:
                return False
            if self.base.is_infinite:
                return True
        c1 = self.base.is_finite
        if c1 is None:
            return
        c2 = self.exp.is_finite
        if c2 is None:
            return
        if c1 and c2:
            if self.exp.is_nonnegative or self.base.is_nonzero:
                return True

    def _eval_is_prime(self):
        if self.exp == S.One:
            return self.base.is_prime
        if self.is_number:
            return self.doit().is_prime

        if self.is_integer and self.is_positive:
            """
            a Power will be non-prime only if both base and exponent
            are greater than 1
            """
            if (self.base-1).is_positive or (self.exp-1).is_positive:
                return False

    def _eval_is_polar(self):
        return self.base.is_polar

    def _eval_subs(self, old, new):
        from sympy import exp, log, Symbol
        def _check(ct1, ct2, old):
            """Return bool, pow where, if bool is True, then the exponent of
            Pow `old` will combine with `pow` so the substitution is valid,
            otherwise bool will be False,

            cti are the coefficient and terms of an exponent of self or old
            In this _eval_subs routine a change like (b**(2*x)).subs(b**x, y)
            will give y**2 since (b**x)**2 == b**(2*x); if that equality does
            not hold then the substitution should not occur so `bool` will be
            False.
            """
            coeff1, terms1 = ct1
            coeff2, terms2 = ct2
            if terms1 == terms2:
                pow = coeff1/coeff2
                try:
                    pow = as_int(pow)
                    combines = True
                except ValueError:
                    combines = Pow._eval_power(
                        Pow(*old.as_base_exp(), evaluate=False),
                        pow) is not None
                return combines, pow
            return False, None

        if old == self.base:
            return new**self.exp._subs(old, new)

        if old.func is self.func and self.base == old.base:
            if self.exp.is_Add is False:
                ct1 = self.exp.as_independent(Symbol, as_Add=False)
                ct2 = old.exp.as_independent(Symbol, as_Add=False)
                ok, pow = _check(ct1, ct2, old)
                if ok:
                    # issue 5180: (x**(6*y)).subs(x**(3*y),z)->z**2
                    return self.func(new, pow)
            else:  # b**(6*x+a).subs(b**(3*x), y) -> y**2 * b**a
                # exp(exp(x) + exp(x**2)).subs(exp(exp(x)), w) -> w * exp(exp(x**2))
                oarg = old.exp
                new_l = []
                o_al = []
                ct2 = oarg.as_coeff_mul()
                for a in self.exp.args:
                    newa = a._subs(old, new)
                    ct1 = newa.as_coeff_mul()
                    ok, pow = _check(ct1, ct2, old)
                    if ok:
                        new_l.append(new**pow)
                        continue
                    o_al.append(newa)
                if new_l:
                    new_l.append(Pow(self.base, Add(*o_al), evaluate=False))
                    return Mul(*new_l)

        if old.func is exp and self.exp.is_real and self.base.is_positive:
            ct1 = old.args[0].as_independent(Symbol, as_Add=False)
            ct2 = (self.exp*log(self.base)).as_independent(
                Symbol, as_Add=False)
            ok, pow = _check(ct1, ct2, old)
            if ok:
                return self.func(new, pow)  # (2**x).subs(exp(x*log(2)), z) -> z

    def as_base_exp(self):
        """Return base and exp of self.

        If base is 1/Integer, then return Integer, -exp. If this extra
        processing is not needed, the base and exp properties will
        give the raw arguments

        Examples
        ========

        >>> from sympy import Pow, S
        >>> p = Pow(S.Half, 2, evaluate=False)
        >>> p.as_base_exp()
        (2, -2)
        >>> p.args
        (1/2, 2)

        """

        b, e = self.args
        if b.is_Rational and b.p == 1 and b.q != 1:
            return Integer(b.q), -e
        return b, e

    def _eval_adjoint(self):
        from sympy.functions.elementary.complexes import adjoint
        i, p = self.exp.is_integer, self.base.is_positive
        if i:
            return adjoint(self.base)**self.exp
        if p:
            return self.base**adjoint(self.exp)
        if i is False and p is False:
            expanded = expand_complex(self)
            if expanded != self:
                return adjoint(expanded)

    def _eval_conjugate(self):
        from sympy.functions.elementary.complexes import conjugate as c
        i, p = self.exp.is_integer, self.base.is_positive
        if i:
            return c(self.base)**self.exp
        if p:
            return self.base**c(self.exp)
        if i is False and p is False:
            expanded = expand_complex(self)
            if expanded != self:
                return c(expanded)

    def _eval_transpose(self):
        from sympy.functions.elementary.complexes import transpose
        i, p = self.exp.is_integer, self.base.is_complex
        if p:
            return self.base**self.exp
        if i:
            return transpose(self.base)**self.exp
        if i is False and p is False:
            expanded = expand_complex(self)
            if expanded != self:
                return transpose(expanded)

    def _eval_expand_power_exp(self, **hints):
        """a**(n+m) -> a**n*a**m"""
        b = self.base
        e = self.exp
        if e.is_Add and e.is_commutative:
            expr = []
            for x in e.args:
                expr.append(self.func(self.base, x))
            return Mul(*expr)
        return self.func(b, e)

    def _eval_expand_power_base(self, **hints):
        """(a*b)**n -> a**n * b**n"""
        force = hints.get('force', False)

        b = self.base
        e = self.exp
        if not b.is_Mul:
            return self

        cargs, nc = b.args_cnc(split_1=False)

        # expand each term - this is top-level-only
        # expansion but we have to watch out for things
        # that don't have an _eval_expand method
        if nc:
            nc = [i._eval_expand_power_base(**hints)
                if hasattr(i, '_eval_expand_power_base') else i
                for i in nc]

            if e.is_Integer:
                if e.is_positive:
                    rv = Mul(*nc*e)
                else:
                    rv = 1/Mul(*nc*-e)
                if cargs:
                    rv *= Mul(*cargs)**e
                return rv

            if not cargs:
                return self.func(Mul(*nc), e, evaluate=False)

            nc = [Mul(*nc)]

        # sift the commutative bases
        def pred(x):
            if x is S.ImaginaryUnit:
                return S.ImaginaryUnit
            polar = x.is_polar
            if polar:
                return True
            if polar is None:
                return fuzzy_bool(x.is_nonnegative)
        sifted = sift(cargs, pred)
        nonneg = sifted[True]
        other = sifted[None]
        neg = sifted[False]
        imag = sifted[S.ImaginaryUnit]
        if imag:
            I = S.ImaginaryUnit
            i = len(imag) % 4
            if i == 0:
                pass
            elif i == 1:
                other.append(I)
            elif i == 2:
                if neg:
                    nonn = -neg.pop()
                    if nonn is not S.One:
                        nonneg.append(nonn)
                else:
                    neg.append(S.NegativeOne)
            else:
                if neg:
                    nonn = -neg.pop()
                    if nonn is not S.One:
                        nonneg.append(nonn)
                else:
                    neg.append(S.NegativeOne)
                other.append(I)
            del imag

        # bring out the bases that can be separated from the base

        if force or e.is_integer:
            # treat all commutatives the same and put nc in other
            cargs = nonneg + neg + other
            other = nc
        else:
            # this is just like what is happening automatically, except
            # that now we are doing it for an arbitrary exponent for which
            # no automatic expansion is done

            assert not e.is_Integer

            # handle negatives by making them all positive and putting
            # the residual -1 in other
            if len(neg) > 1:
                o = S.One
                if not other and neg[0].is_Number:
                    o *= neg.pop(0)
                if len(neg) % 2:
                    o = -o
                for n in neg:
                    nonneg.append(-n)
                if o is not S.One:
                    other.append(o)
            elif neg and other:
                if neg[0].is_Number and neg[0] is not S.NegativeOne:
                    other.append(S.NegativeOne)
                    nonneg.append(-neg[0])
                else:
                    other.extend(neg)
            else:
                other.extend(neg)
            del neg

            cargs = nonneg
            other += nc

        rv = S.One
        if cargs:
            rv *= Mul(*[self.func(b, e, evaluate=False) for b in cargs])
        if other:
            rv *= self.func(Mul(*other), e, evaluate=False)
        return rv

    def _eval_expand_multinomial(self, **hints):
        """(a+b+..) ** n -> a**n + n*a**(n-1)*b + .., n is nonzero integer"""

        base, exp = self.args
        result = self

        if exp.is_Rational and exp.p > 0 and base.is_Add:
            if not exp.is_Integer:
                n = Integer(exp.p // exp.q)

                if not n:
                    return result
                else:
                    radical, result = self.func(base, exp - n), []

                    expanded_base_n = self.func(base, n)
                    if expanded_base_n.is_Pow:
                        expanded_base_n = \
                            expanded_base_n._eval_expand_multinomial()
                    for term in Add.make_args(expanded_base_n):
                        result.append(term*radical)

                    return Add(*result)

            n = int(exp)

            if base.is_commutative:
                order_terms, other_terms = [], []

                for b in base.args:
                    if b.is_Order:
                        order_terms.append(b)
                    else:
                        other_terms.append(b)

                if order_terms:
                    # (f(x) + O(x^n))^m -> f(x)^m + m*f(x)^{m-1} *O(x^n)
                    f = Add(*other_terms)
                    o = Add(*order_terms)

                    if n == 2:
                        return expand_multinomial(f**n, deep=False) + n*f*o
                    else:
                        g = expand_multinomial(f**(n - 1), deep=False)
                        return expand_mul(f*g, deep=False) + n*g*o

                if base.is_number:
                    # Efficiently expand expressions of the form (a + b*I)**n
                    # where 'a' and 'b' are real numbers and 'n' is integer.
                    a, b = base.as_real_imag()

                    if a.is_Rational and b.is_Rational:
                        if not a.is_Integer:
                            if not b.is_Integer:
                                k = self.func(a.q * b.q, n)
                                a, b = a.p*b.q, a.q*b.p
                            else:
                                k = self.func(a.q, n)
                                a, b = a.p, a.q*b
                        elif not b.is_Integer:
                            k = self.func(b.q, n)
                            a, b = a*b.q, b.p
                        else:
                            k = 1

                        a, b, c, d = int(a), int(b), 1, 0

                        while n:
                            if n & 1:
                                c, d = a*c - b*d, b*c + a*d
                                n -= 1
                            a, b = a*a - b*b, 2*a*b
                            n //= 2

                        I = S.ImaginaryUnit

                        if k == 1:
                            return c + I*d
                        else:
                            return Integer(c)/k + I*d/k

                p = other_terms
                # (x+y)**3 -> x**3 + 3*x**2*y + 3*x*y**2 + y**3
                # in this particular example:
                # p = [x,y]; n = 3
                # so now it's easy to get the correct result -- we get the
                # coefficients first:
                from sympy import multinomial_coefficients
                from sympy.polys.polyutils import basic_from_dict
                expansion_dict = multinomial_coefficients(len(p), n)
                # in our example: {(3, 0): 1, (1, 2): 3, (0, 3): 1, (2, 1): 3}
                # and now construct the expression.
                return basic_from_dict(expansion_dict, *p)
            else:
                if n == 2:
                    return Add(*[f*g for f in base.args for g in base.args])
                else:
                    multi = (base**(n - 1))._eval_expand_multinomial()
                    if multi.is_Add:
                        return Add(*[f*g for f in base.args
                            for g in multi.args])
                    else:
                        # XXX can this ever happen if base was an Add?
                        return Add(*[f*multi for f in base.args])
        elif (exp.is_Rational and exp.p < 0 and base.is_Add and
                abs(exp.p) > exp.q):
            return 1 / self.func(base, -exp)._eval_expand_multinomial()
        elif exp.is_Add and base.is_Number:
            #  a + b      a  b
            # n      --> n  n  , where n, a, b are Numbers

            coeff, tail = S.One, S.Zero
            for term in exp.args:
                if term.is_Number:
                    coeff *= self.func(base, term)
                else:
                    tail += term

            return coeff * self.func(base, tail)
        else:
            return result

    def as_real_imag(self, deep=True, **hints):
        from sympy import atan2, cos, im, re, sin
        from sympy.polys.polytools import poly

        if self.exp.is_Integer:
            exp = self.exp
            re, im = self.base.as_real_imag(deep=deep)
            if not im:
                return self, S.Zero
            a, b = symbols('a b', cls=Dummy)
            if exp >= 0:
                if re.is_Number and im.is_Number:
                    # We can be more efficient in this case
                    expr = expand_multinomial(self.base**exp)
                    return expr.as_real_imag()

                expr = poly(
                    (a + b)**exp)  # a = re, b = im; expr = (a + b*I)**exp
            else:
                mag = re**2 + im**2
                re, im = re/mag, -im/mag
                if re.is_Number and im.is_Number:
                    # We can be more efficient in this case
                    expr = expand_multinomial((re + im*S.ImaginaryUnit)**-exp)
                    return expr.as_real_imag()

                expr = poly((a + b)**-exp)

            # Terms with even b powers will be real
            r = [i for i in expr.terms() if not i[0][1] % 2]
            re_part = Add(*[cc*a**aa*b**bb for (aa, bb), cc in r])
            # Terms with odd b powers will be imaginary
            r = [i for i in expr.terms() if i[0][1] % 4 == 1]
            im_part1 = Add(*[cc*a**aa*b**bb for (aa, bb), cc in r])
            r = [i for i in expr.terms() if i[0][1] % 4 == 3]
            im_part3 = Add(*[cc*a**aa*b**bb for (aa, bb), cc in r])

            return (re_part.subs({a: re, b: S.ImaginaryUnit*im}),
            im_part1.subs({a: re, b: im}) + im_part3.subs({a: re, b: -im}))

        elif self.exp.is_Rational:
            re, im = self.base.as_real_imag(deep=deep)

            if im.is_zero and self.exp is S.Half:
                if re.is_nonnegative:
                    return self, S.Zero
                if re.is_nonpositive:
                    return S.Zero, (-self.base)**self.exp

            # XXX: This is not totally correct since for x**(p/q) with
            #      x being imaginary there are actually q roots, but
            #      only a single one is returned from here.
            r = self.func(self.func(re, 2) + self.func(im, 2), S.Half)
            t = atan2(im, re)

            rp, tp = self.func(r, self.exp), t*self.exp

            return (rp*cos(tp), rp*sin(tp))
        else:

            if deep:
                hints['complex'] = False

                expanded = self.expand(deep, **hints)
                if hints.get('ignore') == expanded:
                    return None
                else:
                    return (re(expanded), im(expanded))
            else:
                return (re(self), im(self))

    def _eval_derivative(self, s):
        from sympy import log
        dbase = self.base.diff(s)
        dexp = self.exp.diff(s)
        return self * (dexp * log(self.base) + dbase * self.exp/self.base)

    def _eval_evalf(self, prec):
        base, exp = self.as_base_exp()
        base = base._evalf(prec)
        if not exp.is_Integer:
            exp = exp._evalf(prec)
        if exp.is_negative and base.is_number and base.is_real is False:
            base = base.conjugate() / (base * base.conjugate())._evalf(prec)
            exp = -exp
            return self.func(base, exp).expand()
        return self.func(base, exp)

    def _eval_is_polynomial(self, syms):
        if self.exp.has(*syms):
            return False

        if self.base.has(*syms):
            return bool(self.base._eval_is_polynomial(syms) and
                self.exp.is_Integer and (self.exp >= 0))
        else:
            return True

    def _eval_is_rational(self):
        p = self.func(*self.as_base_exp())  # in case it's unevaluated
        if not p.is_Pow:
            return p.is_rational
        b, e = p.as_base_exp()
        if e.is_Rational and b.is_Rational:
            # we didn't check that e is not an Integer
            # because Rational**Integer autosimplifies
            return False
        if e.is_integer:
            if b.is_rational:
                if b.is_nonzero or e.is_nonnegative:
                    return True
                if b == e:  # always rational, even for 0**0
                    return True
            elif b.is_irrational:
                return e.is_zero

    def _eval_is_algebraic(self):
        if self.base.is_zero or (self.base - 1).is_zero:
            return True
        elif self.exp.is_rational:
            return self.base.is_algebraic
        elif self.base.is_algebraic and self.exp.is_algebraic:
            if ((self.base.is_nonzero and (self.base - 1).is_nonzero)
                or self.base.is_integer is False
                or self.base.is_irrational):
                return self.exp.is_rational

    def _eval_is_rational_function(self, syms):
        if self.exp.has(*syms):
            return False

        if self.base.has(*syms):
            return self.base._eval_is_rational_function(syms) and \
                self.exp.is_Integer
        else:
            return True

    def _eval_is_algebraic_expr(self, syms):
        if self.exp.has(*syms):
            return False

        if self.base.has(*syms):
            return self.base._eval_is_algebraic_expr(syms) and \
                self.exp.is_Rational
        else:
            return True

    def as_numer_denom(self):
        if not self.is_commutative:
            return self, S.One
        base, exp = self.as_base_exp()
        n, d = base.as_numer_denom()
        # this should be the same as ExpBase.as_numer_denom wrt
        # exponent handling
        neg_exp = exp.is_negative
        if not neg_exp and not (-exp).is_negative:
            neg_exp = _coeff_isneg(exp)
        int_exp = exp.is_integer
        # the denominator cannot be separated from the numerator if
        # its sign is unknown unless the exponent is an integer, e.g.
        # sqrt(a/b) != sqrt(a)/sqrt(b) when a=1 and b=-1. But if the
        # denominator is negative the numerator and denominator can
        # be negated and the denominator (now positive) separated.
        if not (d.is_real or int_exp):
            n = base
            d = S.One
        dnonpos = d.is_nonpositive
        if dnonpos:
            n, d = -n, -d
        elif dnonpos is None and not int_exp:
            n = base
            d = S.One
        if neg_exp:
            n, d = d, n
            exp = -exp
        return self.func(n, exp), self.func(d, exp)

    def matches(self, expr, repl_dict={}, old=False):
        expr = _sympify(expr)

        # special case, pattern = 1 and expr.exp can match to 0
        if expr is S.One:
            d = repl_dict.copy()
            d = self.exp.matches(S.Zero, d)
            if d is not None:
                return d

        # make sure the expression to be matched is an Expr
        if not isinstance(expr, Expr):
            return None

        b, e = expr.as_base_exp()

        # special case number
        sb, se = self.as_base_exp()
        if sb.is_Symbol and se.is_Integer and expr:
            if e.is_rational:
                return sb.matches(b**(e/se), repl_dict)
            return sb.matches(expr**(1/se), repl_dict)

        d = repl_dict.copy()
        d = self.base.matches(b, d)
        if d is None:
            return None

        d = self.exp.xreplace(d).matches(e, d)
        if d is None:
            return Expr.matches(self, expr, repl_dict)
        return d

    def _eval_nseries(self, x, n, logx):
        # NOTE! This function is an important part of the gruntz algorithm
        #       for computing limits. It has to return a generalized power
        #       series with coefficients in C(log, log(x)). In more detail:
        # It has to return an expression
        #     c_0*x**e_0 + c_1*x**e_1 + ... (finitely many terms)
        # where e_i are numbers (not necessarily integers) and c_i are
        # expressions involving only numbers, the log function, and log(x).
        from sympy import ceiling, collect, exp, log, O, Order, powsimp
        b, e = self.args
        if e.is_Integer:
            if e > 0:
                # positive integer powers are easy to expand, e.g.:
                # sin(x)**4 = (x-x**3/3+...)**4 = ...
                return expand_multinomial(self.func(b._eval_nseries(x, n=n,
                    logx=logx), e), deep=False)
            elif e is S.NegativeOne:
                # this is also easy to expand using the formula:
                # 1/(1 + x) = 1 - x + x**2 - x**3 ...
                # so we need to rewrite base to the form "1+x"

                nuse = n
                cf = 1

                try:
                    ord = b.as_leading_term(x)
                    cf = Order(ord, x).getn()
                    if cf and cf.is_Number:
                        nuse = n + 2*ceiling(cf)
                    else:
                        cf = 1
                except NotImplementedError:
                    pass

                b_orig, prefactor = b, O(1, x)
                while prefactor.is_Order:
                    nuse += 1
                    b = b_orig._eval_nseries(x, n=nuse, logx=logx)
                    prefactor = b.as_leading_term(x)

                # express "rest" as: rest = 1 + k*x**l + ... + O(x**n)
                rest = expand_mul((b - prefactor)/prefactor)

                if rest.is_Order:
                    return 1/prefactor + rest/prefactor + O(x**n, x)

                k, l = rest.leadterm(x)
                if l.is_Rational and l > 0:
                    pass
                elif l.is_number and l > 0:
                    l = l.evalf()
                elif l == 0:
                    k = k.simplify()
                    if k == 0:
                        # if prefactor == w**4 + x**2*w**4 + 2*x*w**4, we need to
                        # factor the w**4 out using collect:
                        return 1/collect(prefactor, x)
                    else:
                        raise NotImplementedError()
                else:
                    raise NotImplementedError()

                if cf < 0:
                    cf = S.One/abs(cf)

                try:
                    dn = Order(1/prefactor, x).getn()
                    if dn and dn < 0:
                        pass
                    else:
                        dn = 0
                except NotImplementedError:
                    dn = 0

                terms = [1/prefactor]
                for m in range(1, ceiling((n - dn)/l*cf)):
                    new_term = terms[-1]*(-rest)
                    if new_term.is_Pow:
                        new_term = new_term._eval_expand_multinomial(
                            deep=False)
                    else:
                        new_term = expand_mul(new_term, deep=False)
                    terms.append(new_term)
                terms.append(O(x**n, x))
                return powsimp(Add(*terms), deep=True, combine='exp')
            else:
                # negative powers are rewritten to the cases above, for
                # example:
                # sin(x)**(-4) = 1/( sin(x)**4) = ...
                # and expand the denominator:
                nuse, denominator = n, O(1, x)
                while denominator.is_Order:
                    denominator = (b**(-e))._eval_nseries(x, n=nuse, logx=logx)
                    nuse += 1
                if 1/denominator == self:
                    return self
                # now we have a type 1/f(x), that we know how to expand
                return (1/denominator)._eval_nseries(x, n=n, logx=logx)

        if e.has(Symbol):
            return exp(e*log(b))._eval_nseries(x, n=n, logx=logx)

        # see if the base is as simple as possible
        bx = b
        while bx.is_Pow and bx.exp.is_Rational:
            bx = bx.base
        if bx == x:
            return self

        # work for b(x)**e where e is not an Integer and does not contain x
        # and hopefully has no other symbols

        def e2int(e):
            """return the integer value (if possible) of e and a
            flag indicating whether it is bounded or not."""
            n = e.limit(x, 0)
            infinite = n.is_infinite
            if not infinite:
                # XXX was int or floor intended? int used to behave like floor
                # so int(-Rational(1, 2)) returned -1 rather than int's 0
                try:
                    n = int(n)
                except TypeError:
                    #well, the n is something more complicated (like 1+log(2))
                    try:
                        n = int(n.evalf()) + 1  # XXX why is 1 being added?
                    except TypeError:
                        pass  # hope that base allows this to be resolved
                n = _sympify(n)
            return n, infinite

        order = O(x**n, x)
        ei, infinite = e2int(e)
        b0 = b.limit(x, 0)
        if infinite and (b0 is S.One or b0.has(Symbol)):
            # XXX what order
            if b0 is S.One:
                resid = (b - 1)
                if resid.is_positive:
                    return S.Infinity
                elif resid.is_negative:
                    return S.Zero
                raise ValueError('cannot determine sign of %s' % resid)

            return b0**ei

        if (b0 is S.Zero or b0.is_infinite):
            if infinite is not False:
                return b0**e  # XXX what order

            if not ei.is_number:  # if not, how will we proceed?
                raise ValueError(
                    'expecting numerical exponent but got %s' % ei)

            nuse = n - ei

            if e.is_real and e.is_positive:
                lt = b.as_leading_term(x)

                # Try to correct nuse (= m) guess from:
                # (lt + rest + O(x**m))**e =
                # lt**e*(1 + rest/lt + O(x**m)/lt)**e =
                # lt**e + ... + O(x**m)*lt**(e - 1) = ... + O(x**n)
                try:
                    cf = Order(lt, x).getn()
                    nuse = ceiling(n - cf*(e - 1))
                except NotImplementedError:
                    pass

            bs = b._eval_nseries(x, n=nuse, logx=logx)
            terms = bs.removeO()
            if terms.is_Add:
                bs = terms
                lt = terms.as_leading_term(x)

                # bs -> lt + rest -> lt*(1 + (bs/lt - 1))
                return ((self.func(lt, e) * self.func((bs/lt).expand(), e).nseries(
                    x, n=nuse, logx=logx)).expand() + order)

            if bs.is_Add:
                from sympy import O
                # So, bs + O() == terms
                c = Dummy('c')
                res = []
                for arg in bs.args:
                    if arg.is_Order:
                        arg = c*arg.expr
                    res.append(arg)
                bs = Add(*res)
                rv = (bs**e).series(x).subs(c, O(1, x))
                rv += order
                return rv

            rv = bs**e
            if terms != bs:
                rv += order
            return rv

        # either b0 is bounded but neither 1 nor 0 or e is infinite
        # b -> b0 + (b-b0) -> b0 * (1 + (b/b0-1))
        o2 = order*(b0**-e)
        z = (b/b0 - 1)
        o = O(z, x)
        if o is S.Zero or o2 is S.Zero:
            infinite = True
        else:
            if o.expr.is_number:
                e2 = log(o2.expr*x)/log(x)
            else:
                e2 = log(o2.expr)/log(o.expr)
            n, infinite = e2int(e2)
        if infinite:
            # requested accuracy gives infinite series,
            # order is probably non-polynomial e.g. O(exp(-1/x), x).
            r = 1 + z
        else:
            l = []
            g = None
            for i in range(n + 2):
                g = self._taylor_term(i, z, g)
                g = g.nseries(x, n=n, logx=logx)
                l.append(g)
            r = Add(*l)
        return expand_mul(r*b0**e) + order

    def _eval_as_leading_term(self, x):
        from sympy import exp, log
        if not self.exp.has(x):
            return self.func(self.base.as_leading_term(x), self.exp)
        return exp(self.exp * log(self.base)).as_leading_term(x)

    @cacheit
    def _taylor_term(self, n, x, *previous_terms): # of (1+x)**e
        from sympy import binomial
        return binomial(self.exp, n) * self.func(x, n)

    def _sage_(self):
        return self.args[0]._sage_()**self.args[1]._sage_()

    def as_content_primitive(self, radical=False):
        """Return the tuple (R, self/R) where R is the positive Rational
        extracted from self.

        Examples
        ========

        >>> from sympy import sqrt
        >>> sqrt(4 + 4*sqrt(2)).as_content_primitive()
        (2, sqrt(1 + sqrt(2)))
        >>> sqrt(3 + 3*sqrt(2)).as_content_primitive()
        (1, sqrt(3)*sqrt(1 + sqrt(2)))

        >>> from sympy import expand_power_base, powsimp, Mul
        >>> from sympy.abc import x, y

        >>> ((2*x + 2)**2).as_content_primitive()
        (4, (x + 1)**2)
        >>> (4**((1 + y)/2)).as_content_primitive()
        (2, 4**(y/2))
        >>> (3**((1 + y)/2)).as_content_primitive()
        (1, 3**((y + 1)/2))
        >>> (3**((5 + y)/2)).as_content_primitive()
        (9, 3**((y + 1)/2))
        >>> eq = 3**(2 + 2*x)
        >>> powsimp(eq) == eq
        True
        >>> eq.as_content_primitive()
        (9, 3**(2*x))
        >>> powsimp(Mul(*_))
        3**(2*x + 2)

        >>> eq = (2 + 2*x)**y
        >>> s = expand_power_base(eq); s.is_Mul, s
        (False, (2*x + 2)**y)
        >>> eq.as_content_primitive()
        (1, (2*(x + 1))**y)
        >>> s = expand_power_base(_[1]); s.is_Mul, s
        (True, 2**y*(x + 1)**y)

        See docstring of Expr.as_content_primitive for more examples.
        """

        b, e = self.as_base_exp()
        b = _keep_coeff(*b.as_content_primitive(radical=radical))
        ce, pe = e.as_content_primitive(radical=radical)
        if b.is_Rational:
            #e
            #= ce*pe
            #= ce*(h + t)
            #= ce*h + ce*t
            #=> self
            #= b**(ce*h)*b**(ce*t)
            #= b**(cehp/cehq)*b**(ce*t)
            #= b**(iceh+r/cehq)*b**(ce*t)
            #= b**(iceh)*b**(r/cehq)*b**(ce*t)
            #= b**(iceh)*b**(ce*t + r/cehq)
            h, t = pe.as_coeff_Add()
            if h.is_Rational:
                ceh = ce*h
                c = self.func(b, ceh)
                r = S.Zero
                if not c.is_Rational:
                    iceh, r = divmod(ceh.p, ceh.q)
                    c = self.func(b, iceh)
                return c, self.func(b, _keep_coeff(ce, t + r/ce/ceh.q))
        e = _keep_coeff(ce, pe)
        # b**e = (h*t)**e = h**e*t**e = c*m*t**e
        if e.is_Rational and b.is_Mul:
            h, t = b.as_content_primitive(radical=radical)  # h is positive
            c, m = self.func(h, e).as_coeff_Mul()  # so c is positive
            m, me = m.as_base_exp()
            if m is S.One or me == e:  # probably always true
                # return the following, not return c, m*Pow(t, e)
                # which would change Pow into Mul; we let sympy
                # decide what to do by using the unevaluated Mul, e.g
                # should it stay as sqrt(2 + 2*sqrt(5)) or become
                # sqrt(2)*sqrt(1 + sqrt(5))
                return c, self.func(_keep_coeff(m, t), e)
        return S.One, self.func(b, e)

    def is_constant(self, *wrt, **flags):
        expr = self
        if flags.get('simplify', True):
            expr = expr.simplify()
        b, e = expr.as_base_exp()
        bz = b.equals(0)
        if bz:  # recalculate with assumptions in case it's unevaluated
            new = b**e
            if new != expr:
                return new.is_constant()
        econ = e.is_constant(*wrt)
        bcon = b.is_constant(*wrt)
        if bcon:
            if econ:
                return True
            bz = b.equals(0)
            if bz is False:
                return False
        elif bcon is None:
            return None

        return e.equals(0)

from .add import Add
from .numbers import Integer
from .mul import Mul, _keep_coeff
from .symbol import Symbol, Dummy, symbols
