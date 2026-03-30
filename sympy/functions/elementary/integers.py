from __future__ import annotations

from sympy.core.basic import Basic
from sympy.core.expr import Expr

from sympy.core import Add, S
from sympy.core.evalf import get_integer_part, PrecisionExhausted
from sympy.core.function import DefinedFunction
from sympy.core.logic import fuzzy_or, fuzzy_and
from sympy.core.numbers import Integer, int_valued
from sympy.core.relational import Relational, is_eq, is_le, is_lt, is_ge, is_gt
from sympy.functions.elementary.complexes import im, re
from sympy.multipledispatch import dispatch

###############################################################################
######################### FLOOR and CEILING FUNCTIONS #########################
###############################################################################



class RoundFunction(DefinedFunction):
    """Abstract base class for rounding functions."""

    args: tuple[Expr]

    @classmethod
    def eval(cls, arg):
        if (v := cls._eval_number(arg)) is not None:
            return v
        if (v := cls._eval_const_number(arg)) is not None:
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
        intof = lambda x: int(x) if int_valued(x) else (
            x if x.is_integer else None)
        for t in Add.make_args(arg):
            if t.is_imaginary and (i := intof(im(t))) is not None:
                ipart += i*S.ImaginaryUnit
            elif (i := intof(t)) is not None:
                ipart += i
            elif t.is_number:
                npart += t
            else:
                spart += t

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

    @classmethod
    def _eval_number(cls, arg):
        raise NotImplementedError()

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
    .. [2] https://mathworld.wolfram.com/FloorFunction.html

    """
    _dir = -1

    @classmethod
    def _eval_number(cls, arg):
        if arg.is_Number:
            return arg.floor()
        if any(isinstance(i, j)
                for i in (arg, -arg) for j in (floor, ceiling)):
            return arg
        if arg.is_NumberSymbol:
            return arg.approximation_interval(Integer)[0]

    @classmethod
    def _eval_const_number(cls, arg):
        if arg.is_real:
            if arg.is_zero:
                return S.Zero
            if arg.is_positive:
                num, den = arg.as_numer_denom()
                s = den.is_negative
                if s is None:
                    return None
                if s:
                    num, den = -num, -den
                # 0 <= num/den < 1 -> 0
                if is_lt(num, den):
                    return S.Zero
                # 1 <= num/den < 2 -> 1
                if fuzzy_and([is_le(den, num), is_lt(num, 2*den)]):
                    return S.One
            if arg.is_negative:
                num, den = arg.as_numer_denom()
                s = den.is_negative
                if s is None:
                    return None
                if s:
                    num, den = -num, -den
                # -1 <= num/den < 0 -> -1
                if is_le(-den, num):
                    return S.NegativeOne
                # -2 <= num/den < -1 -> -2
                if fuzzy_and([is_le(-2*den, num), is_lt(num, -den)]):
                    return Integer(-2)

    def _eval_as_leading_term(self, x, logx, cdir):
        from sympy.calculus.accumulationbounds import AccumBounds
        arg = self.args[0]
        arg0 = arg.subs(x, 0)
        r = self.subs(x, 0)
        if arg0 is S.NaN or isinstance(arg0, AccumBounds):
            arg0 = arg.limit(x, 0, dir='-' if re(cdir).is_negative else '+')
            r = floor(arg0)
        if arg0.is_finite:
            if arg0 == r:
                ndir = arg.dir(x, cdir=cdir if cdir != 0 else 1)
                if ndir.is_negative:
                    return r - 1
                elif ndir.is_positive:
                    return r
                else:
                    raise NotImplementedError("Not sure of sign of %s" % ndir)
            else:
                return r
        return arg.as_leading_term(x, logx=logx, cdir=cdir)

    def _eval_nseries(self, x, n, logx, cdir=0):
        arg = self.args[0]
        arg0 = arg.subs(x, 0)
        r = self.subs(x, 0)
        if arg0 is S.NaN:
            arg0 = arg.limit(x, 0, dir='-' if re(cdir).is_negative else '+')
            r = floor(arg0)
        if arg0.is_infinite:
            from sympy.calculus.accumulationbounds import AccumBounds
            from sympy.series.order import Order
            s = arg._eval_nseries(x, n, logx, cdir)
            o = Order(1, (x, 0)) if n <= 0 else AccumBounds(-1, 0)
            return s + o
        if arg0 == r:
            ndir = arg.dir(x, cdir=cdir if cdir != 0 else 1)
            if ndir.is_negative:
                return r - 1
            elif ndir.is_positive:
                return r
            else:
                raise NotImplementedError("Not sure of sign of %s" % ndir)
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

@dispatch(floor, Expr)
def _eval_is_eq(lhs, rhs): # noqa:F811
    return is_eq(lhs.rewrite(ceiling), rhs) or \
        is_eq(lhs.rewrite(frac),rhs)

@dispatch(floor, floor)
def _eval_is_ge(lhs,rhs):  # noqa:F811
    x = lhs.args[0]
    y = rhs.args[0]
    if x.is_real and y.is_real:
        if x == y:
            return True
        if x.is_integer:
            if y.is_integer:
                return is_ge(x, y)
            if y.is_number:
                return is_ge(x, floor(y))
        if x.is_number:
            if y.is_integer:
                return is_ge(floor(x), y)
            if y.is_number:
                return is_ge(floor(x), floor(y))
    return None

@dispatch(floor, Expr)
def _eval_is_ge(lhs, rhs): # noqa:F811
    if lhs.args[0].is_real:
        if rhs.is_integer:
            return is_ge(lhs.args[0], rhs)
        if rhs.is_number and rhs.is_real:
            return is_ge(lhs.args[0], ceiling(rhs))
    if lhs.args[0] == rhs and rhs.is_real:
        if rhs.is_integer:
            return True
        if rhs.is_noninteger:
            return False
    if rhs is S.NegativeInfinity and lhs.is_finite:
        return True
    return None

@dispatch(Expr, floor)
def _eval_is_ge(lhs, rhs): # noqa:F811
    if rhs.args[0].is_real and lhs.is_real:
        if lhs.is_integer:
            return is_lt(rhs.args[0], lhs + 1)
        if lhs.is_number:
            return is_lt(rhs.args[0], ceiling(lhs))
        if rhs.args[0] == lhs:
            return True
    if lhs is S.Infinity and rhs.is_finite:
        return True
    return None

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
    .. [2] https://mathworld.wolfram.com/CeilingFunction.html

    """
    _dir = 1

    @classmethod
    def _eval_number(cls, arg):
        if arg.is_Number:
            return arg.ceiling()
        if any(isinstance(i, j)
                for i in (arg, -arg) for j in (floor, ceiling)):
            return arg
        if arg.is_NumberSymbol:
            return arg.approximation_interval(Integer)[1]

    @classmethod
    def _eval_const_number(cls, arg):
        if arg.is_real:
            if arg.is_zero:
                return S.Zero
            if arg.is_positive:
                num, den = arg.as_numer_denom()
                s = den.is_negative
                if s is None:
                    return None
                if s:
                    num, den = -num, -den
                # 0 < num/den <= 1 -> 1
                if is_le(num, den):
                    return S.One
                # 1 < num/den <= 2 -> 2
                if fuzzy_and([is_lt(den, num), is_le(num, 2*den)]):
                    return Integer(2)
            if arg.is_negative:
                num, den = arg.as_numer_denom()
                s = den.is_negative
                if s is None:
                    return None
                if s:
                    num, den = -num, -den
                # -1 < num/den <= 0 -> 0
                if is_lt(-den, num):
                    return S.Zero
                # -2 < num/den <= -1 -> -1
                if fuzzy_and([is_lt(-2*den, num), is_le(num, -den)]):
                    return S.NegativeOne

    def _eval_as_leading_term(self, x, logx, cdir):
        from sympy.calculus.accumulationbounds import AccumBounds
        arg = self.args[0]
        arg0 = arg.subs(x, 0)
        r = self.subs(x, 0)
        if arg0 is S.NaN or isinstance(arg0, AccumBounds):
            arg0 = arg.limit(x, 0, dir='-' if re(cdir).is_negative else '+')
            r = ceiling(arg0)
        if arg0.is_finite:
            if arg0 == r:
                ndir = arg.dir(x, cdir=cdir if cdir != 0 else 1)
                if ndir.is_negative:
                    return r
                elif ndir.is_positive:
                    return r + 1
                else:
                    raise NotImplementedError("Not sure of sign of %s" % ndir)
            else:
                return r
        return arg.as_leading_term(x, logx=logx, cdir=cdir)

    def _eval_nseries(self, x, n, logx, cdir=0):
        arg = self.args[0]
        arg0 = arg.subs(x, 0)
        r = self.subs(x, 0)
        if arg0 is S.NaN:
            arg0 = arg.limit(x, 0, dir='-' if re(cdir).is_negative else '+')
            r = ceiling(arg0)
        if arg0.is_infinite:
            from sympy.calculus.accumulationbounds import AccumBounds
            from sympy.series.order import Order
            s = arg._eval_nseries(x, n, logx, cdir)
            o = Order(1, (x, 0)) if n <= 0 else AccumBounds(0, 1)
            return s + o
        if arg0 == r:
            ndir = arg.dir(x, cdir=cdir if cdir != 0 else 1)
            if ndir.is_negative:
                return r
            elif ndir.is_positive:
                return r + 1
            else:
                raise NotImplementedError("Not sure of sign of %s" % ndir)
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

@dispatch(ceiling, Basic)  # type:ignore
def _eval_is_eq(lhs, rhs): # noqa:F811
    return is_eq(lhs.rewrite(floor), rhs) or is_eq(lhs.rewrite(frac),rhs)

@dispatch(ceiling, ceiling)
def _eval_is_ge(lhs,rhs): # noqa:F811
    x = lhs.args[0]
    y = rhs.args[0]
    if x.is_real and y.is_real:
        if x == y:
            return True
        if x.is_integer:
            if y.is_integer:
                return is_ge(x, y)
            if y.is_number:
                return is_ge(x, ceiling(y))
        if x.is_number:
            if y.is_integer:
                return is_ge(ceiling(x), y)
            if y.is_number:
                return is_ge(ceiling(x), ceiling(y))
    return None

@dispatch(ceiling, floor)
def _eval_is_ge(lhs, rhs): # noqa:F811
    if lhs.args[0].is_real and rhs.args[0].is_real:
        return is_ge(lhs.args[0], rhs.args[0])
    return None

@dispatch(floor, ceiling)
def _eval_is_ge(lhs, rhs): # noqa:F811
    if lhs.args[0].is_real and rhs.args[0].is_real:
        if lhs.args[0].is_integer:
            return is_ge(lhs.args[0], rhs)
        if rhs.args[0].is_integer:
            return is_ge(lhs, rhs.args[0])
        if lhs.args[0].is_number:
            return is_ge(floor(lhs.args[0]), rhs)
        if rhs.args[0].is_number:
            return is_ge(lhs, ceiling(rhs.args[0]))
    return is_ge(floor(lhs.args[0]) - floor(rhs.args[0]), frac(lhs.args[0]) - frac(rhs.args[0]))

@dispatch(ceiling, Expr)
def _eval_is_ge(lhs, rhs): # noqa:F811
    if lhs.args[0].is_real and rhs.is_real:
        if rhs.is_integer:
            return is_gt(lhs.args[0], rhs - 1)
        if rhs.is_number:
            return is_gt(lhs.args[0], floor(rhs))
        if lhs.args[0] == rhs:
            return True
    if rhs is S.NegativeInfinity and lhs.is_finite:
        return True
    return None

@dispatch(Expr, ceiling)
def _eval_is_ge(lhs, rhs): # noqa:F811
    if rhs.args[0].is_real and lhs.is_real:
        if lhs.is_integer:
            return is_le(rhs.args[0], lhs)
        if lhs.is_number:
            return is_le(rhs.args[0], floor(lhs))
        if rhs.args[0] == lhs  and lhs.is_integer is False:
            return False
    if lhs is S.NegativeInfinity and rhs.is_finite:
        return False
    return None

class frac(DefinedFunction):
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
    .. [2] https://mathworld.wolfram.com/FractionalPart.html

    """
    @classmethod
    def eval(cls, arg):
        from sympy.calculus.accumulationbounds import AccumBounds

        def _eval(arg):
            if arg in (S.Infinity, S.NegativeInfinity):
                return AccumBounds(0, 1)
            if arg.is_integer:
                return S.Zero
            if arg.is_number:
                if arg is S.NaN or arg is S.ComplexInfinity:
                    return S.NaN
                return arg - floor(arg)
            return cls(arg, evaluate=False)

        real, imag = S.Zero, S.Zero
        for t in Add.make_args(arg):
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

    def _value_one_or_more(self, other):
        if other.is_extended_real:
            if other.is_number:
                res = other >= 1
                if res and not isinstance(res, Relational):
                    return S.true
            if other.is_integer and other.is_positive:
                return S.true

    def _eval_as_leading_term(self, x, logx, cdir):
        from sympy.calculus.accumulationbounds import AccumBounds
        arg = self.args[0]
        arg0 = arg.subs(x, 0)
        r = self.subs(x, 0)

        if arg0.is_finite:
            if r.is_zero:
                ndir = arg.dir(x, cdir=cdir)
                if ndir.is_negative:
                    return S.One
                return (arg - arg0).as_leading_term(x, logx=logx, cdir=cdir)
            else:
                return r
        elif arg0 in (S.ComplexInfinity, S.Infinity, S.NegativeInfinity):
            return AccumBounds(0, 1)
        return arg.as_leading_term(x, logx=logx, cdir=cdir)

    def _eval_nseries(self, x, n, logx, cdir=0):
        from sympy.series.order import Order
        arg = self.args[0]
        arg0 = arg.subs(x, 0)
        r = self.subs(x, 0)

        if arg0.is_infinite:
            from sympy.calculus.accumulationbounds import AccumBounds
            o = Order(1, (x, 0)) if n <= 0 else AccumBounds(0, 1) + Order(x**n, (x, 0))
            return o
        else:
            res = (arg - arg0)._eval_nseries(x, n, logx=logx, cdir=cdir)
            if r.is_zero:
                ndir = arg.dir(x, cdir=cdir)
                res += S.One if ndir.is_negative else S.Zero
            else:
                res += r
            return res

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

@dispatch(frac, frac)
def _eval_is_ge(lhs, rhs): # noqa:F811
    if lhs.args[0].is_real and rhs.args[0].is_real:
        if lhs.args[0].is_integer:
            return is_ge(0, rhs)
        if rhs.args[0].is_integer:
            return is_ge(lhs, 0)
    return None

@dispatch(frac, floor)
def _eval_is_ge(lhs, rhs): # noqa:F811
    if lhs.args[0].is_real and rhs.args[0].is_real:
        return is_lt(rhs.args[0],1)
    return None

@dispatch(floor, frac)
def _eval_is_ge(lhs, rhs): # noqa:F811
    if lhs.args[0].is_real and rhs.args[0].is_real:
        return is_ge(lhs.args[0], 1)
    return None

@dispatch(frac, ceiling)
def _eval_is_ge(lhs, rhs): # noqa:F811
    if lhs.args[0].is_real and rhs.args[0].is_real:
        return is_le(rhs.args[0],0)
    return None

@dispatch(ceiling, frac)
def _eval_is_ge(lhs, rhs): # noqa:F811
    if lhs.args[0].is_real and rhs.args[0].is_real:
        return is_ge(lhs.args[0],0)
    return None

@dispatch(frac, Expr)
def _eval_is_ge(lhs, rhs): # noqa:F811
    # Check if other <= 0
    if rhs.is_extended_nonpositive:
        return True
    # Check if other >= 1
    res = lhs._value_one_or_more(rhs)
    if res is not None:
        return not(res)
    return None

@dispatch(Expr, frac)
def _eval_is_ge(lhs, rhs): # noqa:F811
    # Check if other < 0
    if lhs.is_extended_negative:
        return False
    # Check if other >= 1
    res = rhs._value_one_or_more(lhs)
    if res is not None:
        return res
    return None
