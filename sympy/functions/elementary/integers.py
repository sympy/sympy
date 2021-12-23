from typing import Tuple as tTuple

from sympy.core.basic import Basic
from sympy.core.expr import Expr

from sympy.core import Add, Mul, S
from sympy.core.evalf import get_integer_part, PrecisionExhausted
from sympy.core.exprtools import gcd_terms
from sympy.core.function import Function
from sympy.core.kind import NumberKind
from sympy.core.logic import fuzzy_or, fuzzy_and, fuzzy_not
from sympy.core.numbers import Integer
from sympy.core.relational import Gt, Lt, Ge, Le, Relational, is_eq
from sympy.core.symbol import Symbol
from sympy.core.sympify import _sympify
from sympy.multipledispatch import dispatch

###############################################################################
######################### FLOOR and CEILING FUNCTIONS #########################
###############################################################################


class RoundFunction(Function):
    """The base class for rounding functions."""

    args: tTuple[Expr]

    @classmethod
    def eval(cls, arg):
        from sympy.functions.elementary.complexes import im
        v = cls._eval_number(arg)
        if v is not None:
            return v

        if arg.is_integer or arg.is_finite is False:
            return arg
        if arg.is_imaginary or (S.ImaginaryUnit*arg).is_real:
            i = im(arg)
            if not i.has(S.ImaginaryUnit):
                return cls(i)*S.ImaginaryUnit
            return cls(arg, evaluate=False)

        # Integral, numerical, symbolic part
        ipart = npart = spart = S.Zero

        # Extract integral (or complex integral) terms
        terms = Add.make_args(arg)

        for t in terms:
            if t.is_integer or (t.is_imaginary and im(t).is_integer):
                ipart += t
            elif t.has(Symbol):
                spart += t
            else:
                npart += t

        if not (npart or spart):
            return ipart

        # Evaluate npart numerically if independent of spart
        if npart and (
            not spart or
            npart.is_real and (spart.is_imaginary or (S.ImaginaryUnit*spart).is_real) or
                npart.is_imaginary and spart.is_real):
            try:
                r, i = get_integer_part(
                    npart, cls._dir, {}, return_ints=True)
                ipart += Integer(r) + Integer(i)*S.ImaginaryUnit
                npart = S.Zero
            except (PrecisionExhausted, NotImplementedError):
                pass

        spart += npart
        if not spart:
            return ipart
        elif spart.is_imaginary or (S.ImaginaryUnit*spart).is_real:
            return ipart + cls(im(spart), evaluate=False)*S.ImaginaryUnit
        elif isinstance(spart, (floor, ceiling)):
            return ipart + spart
        else:
            return ipart + cls(spart, evaluate=False)

    def _eval_is_finite(self):
        return self.args[0].is_finite

    def _eval_is_real(self):
        return self.args[0].is_real

    def _eval_is_integer(self):
        return self.args[0].is_real


class floor(RoundFunction):
    """
    Floor is a univariate function which returns the largest integer
    value not greater than its argument. This implementation
    generalizes floor to complex numbers by taking the floor of the
    real and imaginary parts separately.

    Examples
    ========

    >>> from sympy import floor, E, I, S, Float, Rational
    >>> floor(17)
    17
    >>> floor(Rational(23, 10))
    2
    >>> floor(2*E)
    5
    >>> floor(-Float(0.567))
    -1
    >>> floor(-I/2)
    -I
    >>> floor(S(5)/2 + 5*I/2)
    2 + 2*I

    See Also
    ========

    sympy.functions.elementary.integers.ceiling

    References
    ==========

    .. [1] "Concrete mathematics" by Graham, pp. 87
    .. [2] http://mathworld.wolfram.com/FloorFunction.html

    """
    _dir = -1

    @classmethod
    def _eval_number(cls, arg):
        if arg.is_Number:
            return arg.floor()
        elif any(isinstance(i, j)
                for i in (arg, -arg) for j in (floor, ceiling)):
            return arg
        if arg.is_NumberSymbol:
            return arg.approximation_interval(Integer)[0]

    def _eval_as_leading_term(self, x, logx=None, cdir=0):
        arg = self.args[0]
        arg0 = arg.subs(x, 0)
        r = self.subs(x, 0)
        if arg0.is_finite:
            if arg0 == r:
                if cdir == 0:
                    ndirl = arg.dir(x, cdir=-1)
                    ndir = arg.dir(x, cdir=1)
                    if ndir != ndirl:
                        raise ValueError("Two sided limit of %s around 0"
                                    "does not exist" % self)
                else:
                    ndir = arg.dir(x, cdir=cdir)
                return r - 1 if ndir.is_negative else r
            else:
                return r
        return arg.as_leading_term(x, logx=logx, cdir=cdir)

    def _eval_nseries(self, x, n, logx, cdir=0):
        arg = self.args[0]
        arg0 = arg.subs(x, 0)
        if arg0.is_infinite:
            from sympy.calculus.accumulationbounds import AccumBounds
            from sympy.series.order import Order
            s = arg._eval_nseries(x, n, logx, cdir)
            o = Order(1, (x, 0)) if n <= 0 else AccumBounds(-1, 0)
            return s + o
        r = self.subs(x, 0)
        if arg0 == r:
            ndir = arg.dir(x, cdir=cdir if cdir != 0 else 1)
            return r - 1 if ndir.is_negative else r
        else:
            return r

    def _eval_is_negative(self):
        return self.args[0].is_negative

    def _eval_is_nonnegative(self):
        return self.args[0].is_nonnegative

    def _eval_rewrite_as_ceiling(self, arg, **kwargs):
        return -ceiling(-arg)

    def _eval_rewrite_as_frac(self, arg, **kwargs):
        return arg - frac(arg)

    def __le__(self, other):
        other = S(other)
        if self.args[0].is_real:
            if other.is_integer:
                return self.args[0] < other + 1
            if other.is_number and other.is_real:
                return self.args[0] < ceiling(other)
        if self.args[0] == other and other.is_real:
            return S.true
        if other is S.Infinity and self.is_finite:
            return S.true

        return Le(self, other, evaluate=False)

    def __ge__(self, other):
        other = S(other)
        if self.args[0].is_real:
            if other.is_integer:
                return self.args[0] >= other
            if other.is_number and other.is_real:
                return self.args[0] >= ceiling(other)
        if self.args[0] == other and other.is_real:
            return S.false
        if other is S.NegativeInfinity and self.is_finite:
            return S.true

        return Ge(self, other, evaluate=False)

    def __gt__(self, other):
        other = S(other)
        if self.args[0].is_real:
            if other.is_integer:
                return self.args[0] >= other + 1
            if other.is_number and other.is_real:
                return self.args[0] >= ceiling(other)
        if self.args[0] == other and other.is_real:
            return S.false
        if other is S.NegativeInfinity and self.is_finite:
            return S.true

        return Gt(self, other, evaluate=False)

    def __lt__(self, other):
        other = S(other)
        if self.args[0].is_real:
            if other.is_integer:
                return self.args[0] < other
            if other.is_number and other.is_real:
                return self.args[0] < ceiling(other)
        if self.args[0] == other and other.is_real:
            return S.false
        if other is S.Infinity and self.is_finite:
            return S.true

        return Lt(self, other, evaluate=False)


@dispatch(floor, Expr)
def _eval_is_eq(lhs, rhs): # noqa:F811
   return is_eq(lhs.rewrite(ceiling), rhs) or \
        is_eq(lhs.rewrite(frac),rhs)


class ceiling(RoundFunction):
    """
    Ceiling is a univariate function which returns the smallest integer
    value not less than its argument. This implementation
    generalizes ceiling to complex numbers by taking the ceiling of the
    real and imaginary parts separately.

    Examples
    ========

    >>> from sympy import ceiling, E, I, S, Float, Rational
    >>> ceiling(17)
    17
    >>> ceiling(Rational(23, 10))
    3
    >>> ceiling(2*E)
    6
    >>> ceiling(-Float(0.567))
    0
    >>> ceiling(I/2)
    I
    >>> ceiling(S(5)/2 + 5*I/2)
    3 + 3*I

    See Also
    ========

    sympy.functions.elementary.integers.floor

    References
    ==========

    .. [1] "Concrete mathematics" by Graham, pp. 87
    .. [2] http://mathworld.wolfram.com/CeilingFunction.html

    """
    _dir = 1

    @classmethod
    def _eval_number(cls, arg):
        if arg.is_Number:
            return arg.ceiling()
        elif any(isinstance(i, j)
                for i in (arg, -arg) for j in (floor, ceiling)):
            return arg
        if arg.is_NumberSymbol:
            return arg.approximation_interval(Integer)[1]

    def _eval_as_leading_term(self, x, logx=None, cdir=0):
        arg = self.args[0]
        arg0 = arg.subs(x, 0)
        r = self.subs(x, 0)
        if arg0.is_finite:
            if arg0 == r:
                if cdir == 0:
                    ndirl = arg.dir(x, cdir=-1)
                    ndir = arg.dir(x, cdir=1)
                    if ndir != ndirl:
                        raise ValueError("Two sided limit of %s around 0"
                                    "does not exist" % self)
                else:
                    ndir = arg.dir(x, cdir=cdir)
                return r if ndir.is_negative else r + 1
            else:
                return r
        return arg.as_leading_term(x, logx=logx, cdir=cdir)

    def _eval_nseries(self, x, n, logx, cdir=0):
        arg = self.args[0]
        arg0 = arg.subs(x, 0)
        if arg0.is_infinite:
            from sympy.calculus.accumulationbounds import AccumBounds
            from sympy.series.order import Order
            s = arg._eval_nseries(x, n, logx, cdir)
            o = Order(1, (x, 0)) if n <= 0 else AccumBounds(0, 1)
            return s + o
        r = self.subs(x, 0)
        if arg0 == r:
            ndir = arg.dir(x, cdir=cdir if cdir != 0 else 1)
            return r if ndir.is_negative else r + 1
        else:
            return r

    def _eval_rewrite_as_floor(self, arg, **kwargs):
        return -floor(-arg)

    def _eval_rewrite_as_frac(self, arg, **kwargs):
        return arg + frac(-arg)

    def _eval_is_positive(self):
        return self.args[0].is_positive

    def _eval_is_nonpositive(self):
        return self.args[0].is_nonpositive

    def __lt__(self, other):
        other = S(other)
        if self.args[0].is_real:
            if other.is_integer:
                return self.args[0] <= other - 1
            if other.is_number and other.is_real:
                return self.args[0] <= floor(other)
        if self.args[0] == other and other.is_real:
            return S.false
        if other is S.Infinity and self.is_finite:
            return S.true

        return Lt(self, other, evaluate=False)

    def __gt__(self, other):
        other = S(other)
        if self.args[0].is_real:
            if other.is_integer:
                return self.args[0] > other
            if other.is_number and other.is_real:
                return self.args[0] > floor(other)
        if self.args[0] == other and other.is_real:
            return S.false
        if other is S.NegativeInfinity and self.is_finite:
            return S.true

        return Gt(self, other, evaluate=False)

    def __ge__(self, other):
        other = S(other)
        if self.args[0].is_real:
            if other.is_integer:
                return self.args[0] > other - 1
            if other.is_number and other.is_real:
                return self.args[0] > floor(other)
        if self.args[0] == other and other.is_real:
            return S.true
        if other is S.NegativeInfinity and self.is_finite:
            return S.true

        return Ge(self, other, evaluate=False)

    def __le__(self, other):
        other = S(other)
        if self.args[0].is_real:
            if other.is_integer:
                return self.args[0] <= other
            if other.is_number and other.is_real:
                return self.args[0] <= floor(other)
        if self.args[0] == other and other.is_real:
            return S.false
        if other is S.Infinity and self.is_finite:
            return S.true

        return Le(self, other, evaluate=False)


@dispatch(ceiling, Basic)  # type:ignore
def _eval_is_eq(lhs, rhs): # noqa:F811
    return is_eq(lhs.rewrite(floor), rhs) or is_eq(lhs.rewrite(frac),rhs)


class frac(Function):
    r"""Represents the fractional part of x

    For real numbers it is defined [1]_ as

    .. math::
        x - \left\lfloor{x}\right\rfloor

    Examples
    ========

    >>> from sympy import Symbol, frac, Rational, floor, I
    >>> frac(Rational(4, 3))
    1/3
    >>> frac(-Rational(4, 3))
    2/3

    returns zero for integer arguments

    >>> n = Symbol('n', integer=True)
    >>> frac(n)
    0

    rewrite as floor

    >>> x = Symbol('x')
    >>> frac(x).rewrite(floor)
    x - floor(x)

    for complex arguments

    >>> r = Symbol('r', real=True)
    >>> t = Symbol('t', real=True)
    >>> frac(t + I*r)
    I*frac(r) + frac(t)

    See Also
    ========

    sympy.functions.elementary.integers.floor
    sympy.functions.elementary.integers.ceiling

    References
    ===========

    .. [1] https://en.wikipedia.org/wiki/Fractional_part
    .. [2] http://mathworld.wolfram.com/FractionalPart.html

    """
    @classmethod
    def eval(cls, arg):
        from sympy.calculus.accumulationbounds import AccumBounds
        from sympy.functions.elementary.complexes import im

        def _eval(arg):
            if arg in (S.Infinity, S.NegativeInfinity):
                return AccumBounds(0, 1)
            if arg.is_integer:
                return S.Zero
            if arg.is_number:
                if arg is S.NaN:
                    return S.NaN
                elif arg is S.ComplexInfinity:
                    return S.NaN
                else:
                    return arg - floor(arg)
            return cls(arg, evaluate=False)

        terms = Add.make_args(arg)
        real, imag = S.Zero, S.Zero
        for t in terms:
            # Two checks are needed for complex arguments
            # see issue-7649 for details
            if t.is_imaginary or (S.ImaginaryUnit*t).is_real:
                i = im(t)
                if not i.has(S.ImaginaryUnit):
                    imag += i
                else:
                    real += t
            else:
                real += t

        real = _eval(real)
        imag = _eval(imag)
        return real + S.ImaginaryUnit*imag

    def _eval_rewrite_as_floor(self, arg, **kwargs):
        return arg - floor(arg)

    def _eval_rewrite_as_ceiling(self, arg, **kwargs):
        return arg + ceiling(-arg)

    def _eval_is_finite(self):
        return True

    def _eval_is_real(self):
        return self.args[0].is_extended_real

    def _eval_is_imaginary(self):
        return self.args[0].is_imaginary

    def _eval_is_integer(self):
        return self.args[0].is_integer

    def _eval_is_zero(self):
        return fuzzy_or([self.args[0].is_zero, self.args[0].is_integer])

    def _eval_is_negative(self):
        return False

    def __ge__(self, other):
        if self.is_extended_real:
            other = _sympify(other)
            # Check if other <= 0
            if other.is_extended_nonpositive:
                return S.true
            # Check if other >= 1
            res = self._value_one_or_more(other)
            if res is not None:
                return not(res)
        return Ge(self, other, evaluate=False)

    def __gt__(self, other):
        if self.is_extended_real:
            other = _sympify(other)
            # Check if other < 0
            res = self._value_one_or_more(other)
            if res is not None:
                return not(res)
            # Check if other >= 1
            if other.is_extended_negative:
                return S.true
        return Gt(self, other, evaluate=False)

    def __le__(self, other):
        if self.is_extended_real:
            other = _sympify(other)
            # Check if other < 0
            if other.is_extended_negative:
                return S.false
            # Check if other >= 1
            res = self._value_one_or_more(other)
            if res is not None:
                return res
        return Le(self, other, evaluate=False)

    def __lt__(self, other):
        if self.is_extended_real:
            other = _sympify(other)
            # Check if other <= 0
            if other.is_extended_nonpositive:
                return S.false
            # Check if other >= 1
            res = self._value_one_or_more(other)
            if res is not None:
                return res
        return Lt(self, other, evaluate=False)

    def _value_one_or_more(self, other):
        if other.is_extended_real:
            if other.is_number:
                res = other >= 1
                if res and not isinstance(res, Relational):
                    return S.true
            if other.is_integer and other.is_positive:
                return S.true


@dispatch(frac, Basic)  # type:ignore
def _eval_is_eq(lhs, rhs): # noqa:F811
    if (lhs.rewrite(floor) == rhs) or \
        (lhs.rewrite(ceiling) == rhs):
        return True
    # Check if other < 0
    if rhs.is_extended_negative:
        return False
    # Check if other >= 1
    res = lhs._value_one_or_more(rhs)
    if res is not None:
        return False


class Mod(Function):
    """Represents a modulo operation on symbolic expressions.

    Parameters
    ==========

    p : Expr
        Dividend.

    q : Expr
        Divisor.

    Notes
    =====

    The convention used is the same as Python's: the remainder always has the
    same sign as the divisor.

    Examples
    ========

    >>> from sympy.abc import x, y
    >>> x**2 % y
    Mod(x**2, y)
    >>> _.subs({x: 5, y: 6})
    1

    """

    kind = NumberKind

    @classmethod
    def eval(cls, p, q):
        def doit(p, q):
            """Try to return p % q if both are numbers or +/-p is known
            to be less than or equal q.
            """

            if q.is_zero:
                raise ZeroDivisionError("Modulo by zero")
            if p is S.NaN or q is S.NaN or p.is_finite is False or q.is_finite is False:
                return S.NaN
            if p is S.Zero or p in (q, -q) or (p.is_integer and q == 1):
                return S.Zero

            if q.is_Number:
                if p.is_Number:
                    return p%q
                if q == 2:
                    if p.is_even:
                        return S.Zero
                    elif p.is_odd:
                        return S.One

            if hasattr(p, '_eval_Mod'):
                rv = getattr(p, '_eval_Mod')(q)
                if rv is not None:
                    return rv

            # by ratio
            r = p/q
            if r.is_integer:
                return S.Zero
            try:
                d = int(r)
            except TypeError:
                pass
            else:
                if isinstance(d, int):
                    rv = p - d*q
                    if (rv*q < 0) == True:
                        rv += q
                    return rv

            # by difference
            # -2|q| < p < 2|q|
            d = abs(p)
            for _ in range(2):
                d -= abs(q)
                if d.is_negative:
                    if q.is_positive:
                        if p.is_positive:
                            return d + q
                        elif p.is_negative:
                            return -d
                    elif q.is_negative:
                        if p.is_positive:
                            return d
                        elif p.is_negative:
                            return -d + q
                    break

        rv = doit(p, q)
        if rv is not None:
            return rv

        # denest
        if isinstance(p, cls):
            qinner = p.args[1]
            if qinner % q == 0:
                return cls(p.args[0], q)
            elif (qinner*(q - qinner)).is_nonnegative:
                # |qinner| < |q| and have same sign
                return p
        elif isinstance(-p, cls):
            qinner = (-p).args[1]
            if qinner % q == 0:
                return cls(-(-p).args[0], q)
            elif (qinner*(q + qinner)).is_nonpositive:
                # |qinner| < |q| and have different sign
                return p
        elif isinstance(p, Add):
            # separating into modulus and non modulus
            both_l = non_mod_l, mod_l = [], []
            for arg in p.args:
                both_l[isinstance(arg, cls)].append(arg)
            # if q same for all
            if mod_l and all(inner.args[1] == q for inner in mod_l):
                net = Add(*non_mod_l) + Add(*[i.args[0] for i in mod_l])
                return cls(net, q)

        elif isinstance(p, Mul):
            # separating into modulus and non modulus
            both_l = non_mod_l, mod_l = [], []
            for arg in p.args:
                both_l[isinstance(arg, cls)].append(arg)

            if mod_l and all(inner.args[1] == q for inner in mod_l):
                # finding distributive term
                non_mod_l = [cls(x, q) for x in non_mod_l]
                mod = []
                non_mod = []
                for j in non_mod_l:
                    if isinstance(j, cls):
                        mod.append(j.args[0])
                    else:
                        non_mod.append(j)
                prod_mod = Mul(*mod)
                prod_non_mod = Mul(*non_mod)
                prod_mod1 = Mul(*[i.args[0] for i in mod_l])
                net = prod_mod1*prod_mod
                return prod_non_mod*cls(net, q)

            if q.is_Integer and q is not S.One:
                _ = []
                for i in non_mod_l:
                    if i.is_Integer and (i % q is not S.Zero):
                        _.append(i%q)
                    else:
                        _.append(i)
                non_mod_l = _

            p = Mul(*(non_mod_l + mod_l))

        # XXX other possibilities?

        from sympy.polys.polyerrors import PolynomialError
        from sympy.polys.polytools import gcd

        # extract gcd; any further simplification should be done by the user
        try:
            G = gcd(p, q)
            if G != 1:
                p, q = [gcd_terms(i/G, clear=False, fraction=False)
                        for i in (p, q)]
        except PolynomialError:  # issue 21373
            G = S.One
        pwas, qwas = p, q

        # simplify terms
        # (x + y + 2) % x -> Mod(y + 2, x)
        if p.is_Add:
            args = []
            for i in p.args:
                a = cls(i, q)
                if a.count(cls) > i.count(cls):
                    args.append(i)
                else:
                    args.append(a)
            if args != list(p.args):
                p = Add(*args)

        else:
            # handle coefficients if they are not Rational
            # since those are not handled by factor_terms
            # e.g. Mod(.6*x, .3*y) -> 0.3*Mod(2*x, y)
            cp, p = p.as_coeff_Mul()
            cq, q = q.as_coeff_Mul()
            ok = False
            if not cp.is_Rational or not cq.is_Rational:
                r = cp % cq
                if r == 0:
                    G *= cq
                    p *= int(cp/cq)
                    ok = True
            if not ok:
                p = cp*p
                q = cq*q

        # simple -1 extraction
        if p.could_extract_minus_sign() and q.could_extract_minus_sign():
            G, p, q = [-i for i in (G, p, q)]

        # check again to see if p and q can now be handled as numbers
        rv = doit(p, q)
        if rv is not None:
            return rv*G

        # put 1.0 from G on inside
        if G.is_Float and G == 1:
            p *= G
            return cls(p, q, evaluate=False)
        elif G.is_Mul and G.args[0].is_Float and G.args[0] == 1:
            p = G.args[0]*p
            G = Mul._from_args(G.args[1:])
        return G*cls(p, q, evaluate=(p, q) != (pwas, qwas))

    def _eval_is_integer(self):
        p, q = self.args
        if fuzzy_and([p.is_integer, q.is_integer, fuzzy_not(q.is_zero)]):
            return True

    def _eval_is_nonnegative(self):
        if self.args[1].is_positive:
            return True

    def _eval_is_nonpositive(self):
        if self.args[1].is_negative:
            return True

    def _eval_rewrite_as_floor(self, a, b, **kwargs):
        from sympy.functions.elementary.integers import floor
        return a - b*floor(a/b)
