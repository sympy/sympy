from __future__ import annotations
from itertools import product

from sympy.core.add import Add
from sympy.core.cache import cacheit
from sympy.core.function import (DefinedFunction, ArgumentIndexError, expand_log,
    expand_mul, FunctionClass, PoleError, expand_multinomial, expand_complex)
from sympy.core.logic import fuzzy_and, fuzzy_not, fuzzy_or
from sympy.core.mul import Mul
from sympy.core.numbers import Integer, Rational, pi, I
from sympy.core.parameters import global_parameters
from sympy.core.power import Pow
from sympy.core.singleton import S
from sympy.core.symbol import Wild, Dummy
from sympy.core.sympify import sympify
from sympy.functions.combinatorial.factorials import factorial
from sympy.functions.elementary.complexes import arg, unpolarify, im, re, Abs
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.ntheory import multiplicity, perfect_power
from sympy.ntheory.factor_ import factorint
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from sympy.core.expr import Expr

# NOTE IMPORTANT
# The series expansion code in this file is an important part of the gruntz
# algorithm for determining limits. _eval_nseries has to return a generalized
# power series with coefficients in C(log(x), log).
# In more detail, the result of _eval_nseries(self, x, n) must be
#   c_0*x**e_0 + ... (finitely many terms)
# where e_i are numbers (not necessarily integers) and c_i involve only
# numbers, the function log, and log(x). [This also means it must not contain
# log(x(1+p)), this *has* to be expanded to log(x)+log(1+p) if x.is_positive and
# p.is_positive.]


class ExpBase(DefinedFunction):

    unbranched = True
    _singularities = (S.ComplexInfinity,)

    @property
    def kind(self):
        return self.exp.kind

    def inverse(self, argindex=1):
        """
        Returns the inverse function of ``exp(x)``.
        """
        return log

    def as_numer_denom(self):
        """
        Returns this with a positive exponent as a 2-tuple (a fraction).

        Examples
        ========

        >>> from sympy import exp
        >>> from sympy.abc import x
        >>> exp(-x).as_numer_denom()
        (1, exp(x))
        >>> exp(x).as_numer_denom()
        (exp(x), 1)
        """
        # this should be the same as Pow.as_numer_denom wrt
        # exponent handling
        if not self.is_commutative:
            return self, S.One
        exp = self.exp
        neg_exp = exp.is_negative
        if not neg_exp and not (-exp).is_negative:
            neg_exp = exp.could_extract_minus_sign()
        if neg_exp:
            return S.One, self.func(-exp)
        return self, S.One

    @property
    def exp(self):
        """
        Returns the exponent of the function.
        """
        return self.args[0]

    def as_base_exp(self):
        """
        Returns the 2-tuple (base, exponent).
        """
        return self.func(1), Mul(*self.args)

    def _eval_adjoint(self):
        return self.func(self.exp.adjoint())

    def _eval_conjugate(self):
        return self.func(self.exp.conjugate())

    def _eval_transpose(self):
        return self.func(self.exp.transpose())

    def _eval_is_finite(self):
        arg = self.exp
        if arg.is_infinite:
            if arg.is_extended_negative:
                return True
            if arg.is_extended_positive:
                return False
        if arg.is_finite:
            return True

    def _eval_is_rational(self):
        s = self.func(*self.args)
        if s.func == self.func:
            z = s.exp.is_zero
            if z:
                return True
            elif s.exp.is_rational and fuzzy_not(z):
                return False
        else:
            return s.is_rational

    def _eval_is_zero(self):
        return self.exp is S.NegativeInfinity

    def _eval_power(self, other):
        """exp(arg)**e -> exp(arg*e) if assumptions allow it.
        """
        b, e = self.as_base_exp()
        return Pow._eval_power(Pow(b, e, evaluate=False), other)

    def _eval_expand_power_exp(self, **hints):
        from sympy.concrete.products import Product
        from sympy.concrete.summations import Sum
        arg = self.args[0]
        if arg.is_Add and arg.is_commutative:
            return Mul.fromiter(self.func(x) for x in arg.args)
        elif isinstance(arg, Sum) and arg.is_commutative:
            return Product(self.func(arg.function), *arg.limits)
        return self.func(arg)


class exp_polar(ExpBase):
    r"""
    Represent a *polar number* (see g-function Sphinx documentation).

    Explanation
    ===========

    ``exp_polar`` represents the function
    `Exp: \mathbb{C} \rightarrow \mathcal{S}`, sending the complex number
    `z = a + bi` to the polar number `r = exp(a), \theta = b`. It is one of
    the main functions to construct polar numbers.

    Examples
    ========

    >>> from sympy import exp_polar, pi, I, exp

    The main difference is that polar numbers do not "wrap around" at `2 \pi`:

    >>> exp(2*pi*I)
    1
    >>> exp_polar(2*pi*I)
    exp_polar(2*I*pi)

    apart from that they behave mostly like classical complex numbers:

    >>> exp_polar(2)*exp_polar(3)
    exp_polar(5)

    See Also
    ========

    sympy.simplify.powsimp.powsimp
    polar_lift
    periodic_argument
    principal_branch
    """

    is_polar = True
    is_comparable = False  # cannot be evalf'd

    def _eval_Abs(self):   # Abs is never a polar number
        return exp(re(self.args[0]))

    def _eval_evalf(self, prec):
        """ Careful! any evalf of polar numbers is flaky """
        i = im(self.args[0])
        try:
            bad = (i <= -pi or i > pi)
        except TypeError:
            bad = True
        if bad:
            return self  # cannot evalf for this argument
        res = exp(self.args[0])._eval_evalf(prec)
        if i > 0 and im(res) < 0:
            # i ~ pi, but exp(I*i) evaluated to argument slightly bigger than pi
            return re(res)
        return res

    def _eval_power(self, other):
        return self.func(self.args[0]*other)

    def _eval_is_extended_real(self):
        if self.args[0].is_extended_real:
            return True

    def as_base_exp(self):
        # XXX exp_polar(0) is special!
        if self.args[0] == 0:
            return self, S.One
        return ExpBase.as_base_exp(self)


class ExpMeta(FunctionClass):
    def __instancecheck__(cls, instance):
        if exp in instance.__class__.__mro__:
            return True
        return isinstance(instance, Pow) and instance.base is S.Exp1


class exp(ExpBase, metaclass=ExpMeta):
    """
    The exponential function, :math:`e^x`.

    Examples
    ========

    >>> from sympy import exp, I, pi
    >>> from sympy.abc import x
    >>> exp(x)
    exp(x)
    >>> exp(x).diff(x)
    exp(x)
    >>> exp(I*pi)
    -1

    Parameters
    ==========

    arg : Expr

    See Also
    ========

    sympy.functions.elementary.exponential.log
    """

    def fdiff(self, argindex=1):
        """
        Returns the first derivative of this function.
        """
        if argindex == 1:
            return self
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_refine(self, assumptions):
        from sympy.assumptions import ask, Q
        arg = self.args[0]
        if arg.is_Mul:
            Ioo = I*S.Infinity
            if arg in [Ioo, -Ioo]:
                return S.NaN

            coeff = arg.as_coefficient(pi*I)
            if coeff:
                if ask(Q.integer(2*coeff)):
                    if ask(Q.even(coeff)):
                        return S.One
                    elif ask(Q.odd(coeff)):
                        return S.NegativeOne
                    elif ask(Q.even(coeff + S.Half)):
                        return -I
                    elif ask(Q.odd(coeff + S.Half)):
                        return I

    @classmethod
    def eval(cls, arg):
        from sympy.calculus import AccumBounds
        from sympy.matrices.matrixbase import MatrixBase
        from sympy.sets.setexpr import SetExpr
        from sympy.simplify.simplify import logcombine
        if isinstance(arg, MatrixBase):
            return arg.exp()
        elif global_parameters.exp_is_pow:
            return Pow(S.Exp1, arg)
        elif arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg.is_zero:
                return S.One
            elif arg is S.One:
                return S.Exp1
            elif arg is S.Infinity:
                return S.Infinity
            elif arg is S.NegativeInfinity:
                return S.Zero
        elif arg is S.ComplexInfinity:
            return S.NaN
        elif isinstance(arg, log):
            return arg.args[0]
        elif isinstance(arg, AccumBounds):
            return AccumBounds(exp(arg.min), exp(arg.max))
        elif isinstance(arg, SetExpr):
            return arg._eval_func(cls)
        elif arg.is_Mul:
            coeff = arg.as_coefficient(pi*I)
            if coeff:
                if (2*coeff).is_integer:
                    if coeff.is_even:
                        return S.One
                    elif coeff.is_odd:
                        return S.NegativeOne
                    elif (coeff + S.Half).is_even:
                        return -I
                    elif (coeff + S.Half).is_odd:
                        return I
                elif coeff.is_Rational:
                    ncoeff = coeff % 2 # restrict to [0, 2pi)
                    if ncoeff > 1: # restrict to (-pi, pi]
                        ncoeff -= 2
                    if ncoeff != coeff:
                        return cls(ncoeff*pi*I)

            # Warning: code in risch.py will be very sensitive to changes
            # in this (see DifferentialExtension).

            # look for a single log factor

            coeff, terms = arg.as_coeff_Mul()

            # but it can't be multiplied by oo
            if coeff in [S.NegativeInfinity, S.Infinity]:
                if terms.is_number:
                    if coeff is S.NegativeInfinity:
                        terms = -terms
                    if re(terms).is_zero and terms is not S.Zero:
                        return S.NaN
                    if re(terms).is_positive and im(terms) is not S.Zero:
                        return S.ComplexInfinity
                    if re(terms).is_negative:
                        return S.Zero
                return None

            coeffs, log_term = [coeff], None
            for term in Mul.make_args(terms):
                term_ = logcombine(term)
                if isinstance(term_, log):
                    if log_term is None:
                        log_term = term_.args[0]
                    else:
                        return None
                elif term.is_comparable:
                    coeffs.append(term)
                else:
                    return None

            return log_term**Mul(*coeffs) if log_term else None

        elif arg.is_Add:
            out = []
            add = []
            argchanged = False
            for a in arg.args:
                if a is S.One:
                    add.append(a)
                    continue
                newa = cls(a)
                if isinstance(newa, cls):
                    if newa.args[0] != a:
                        add.append(newa.args[0])
                        argchanged = True
                    else:
                        add.append(a)
                else:
                    out.append(newa)
            if out or argchanged:
                return Mul(*out)*cls(Add(*add), evaluate=False)

        if arg.is_zero:
            return S.One

    @property
    def base(self):
        """
        Returns the base of the exponential function.
        """
        return S.Exp1

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        """
        Calculates the next term in the Taylor series expansion.
        """
        if n < 0:
            return S.Zero
        if n == 0:
            return S.One
        x = sympify(x)
        if previous_terms:
            p = previous_terms[-1]
            if p is not None:
                return p * x / n
        return x**n/factorial(n)

    def as_real_imag(self, deep=True, **hints):
        """
        Returns this function as a 2-tuple representing a complex number.

        Examples
        ========

        >>> from sympy import exp, I
        >>> from sympy.abc import x
        >>> exp(x).as_real_imag()
        (exp(re(x))*cos(im(x)), exp(re(x))*sin(im(x)))
        >>> exp(1).as_real_imag()
        (E, 0)
        >>> exp(I).as_real_imag()
        (cos(1), sin(1))
        >>> exp(1+I).as_real_imag()
        (E*cos(1), E*sin(1))

        See Also
        ========

        sympy.functions.elementary.complexes.re
        sympy.functions.elementary.complexes.im
        """
        from sympy.functions.elementary.trigonometric import cos, sin
        re, im = self.args[0].as_real_imag()
        if deep:
            re = re.expand(deep, **hints)
            im = im.expand(deep, **hints)
        cos, sin = cos(im), sin(im)
        return (exp(re)*cos, exp(re)*sin)

    def _eval_subs(self, old, new):
        # keep processing of power-like args centralized in Pow
        if old.is_Pow:  # handle (exp(3*log(x))).subs(x**2, z) -> z**(3/2)
            old = exp(old.exp*log(old.base))
        elif old is S.Exp1 and new.is_Function:
            old = exp
        if isinstance(old, exp) or old is S.Exp1:
            f = lambda a: Pow(*a.as_base_exp(), evaluate=False) if (
                a.is_Pow or isinstance(a, exp)) else a
            return Pow._eval_subs(f(self), f(old), new)

        if old is exp and not new.is_Function:
            return new**self.exp._subs(old, new)
        return super()._eval_subs(old, new)

    def _eval_is_extended_real(self):
        if self.args[0].is_extended_real:
            return True
        elif self.args[0].is_imaginary:
            arg2 = -S(2) * I * self.args[0] / pi
            return arg2.is_even

    def _eval_is_complex(self):
        def complex_extended_negative(arg):
            yield arg.is_complex
            yield arg.is_extended_negative
        return fuzzy_or(complex_extended_negative(self.args[0]))

    def _eval_is_algebraic(self):
        if (self.exp / pi / I).is_rational:
            return True
        if fuzzy_not(self.exp.is_zero):
            if self.exp.is_algebraic:
                return False
            elif (self.exp / pi).is_rational:
                return False

    def _eval_is_extended_positive(self):
        if self.exp.is_extended_real:
            return self.args[0] is not S.NegativeInfinity
        elif self.exp.is_imaginary:
            arg2 = -I * self.args[0] / pi
            return arg2.is_even

    def _eval_nseries(self, x, n, logx, cdir=0):
        # NOTE Please see the comment at the beginning of this file, labelled
        #      IMPORTANT.
        from sympy.functions.elementary.complexes import sign
        from sympy.functions.elementary.integers import ceiling
        from sympy.series.limits import limit
        from sympy.series.order import Order
        from sympy.simplify.powsimp import powsimp
        arg = self.exp
        arg_series = arg._eval_nseries(x, n=n, logx=logx)
        if arg_series.is_Order:
            return 1 + arg_series
        arg0 = limit(arg_series.removeO(), x, 0)
        if arg0 is S.NegativeInfinity:
            return Order(x**n, x)
        if arg0 is S.Infinity:
            return self
        if arg0.is_infinite:
            raise PoleError("Cannot expand %s around 0" % (self))
        # checking for indecisiveness/ sign terms in arg0
        if any(isinstance(arg, sign) for arg in arg0.args):
            return self
        t = Dummy("t")
        nterms = n
        try:
            cf = Order(arg.as_leading_term(x, logx=logx), x).getn()
        except (NotImplementedError, PoleError):
            cf = 0
        if cf and cf > 0:
            nterms = ceiling(n/cf)
        exp_series = exp(t)._taylor(t, nterms)
        r = exp(arg0)*exp_series.subs(t, arg_series - arg0)
        rep = {logx: log(x)} if logx is not None else {}
        if r.subs(rep) == self:
            return r
        if cf and cf > 1:
            r += Order((arg_series - arg0)**n, x)/x**((cf-1)*n)
        else:
            r += Order((arg_series - arg0)**n, x)
        r = r.expand()
        r = powsimp(r, deep=True, combine='exp')
        # powsimp may introduce unexpanded (-1)**Rational; see PR #17201
        simplerat = lambda x: x.is_Rational and x.q in [3, 4, 6]
        w = Wild('w', properties=[simplerat])
        r = r.replace(S.NegativeOne**w, expand_complex(S.NegativeOne**w))
        return r

    def _taylor(self, x, n):
        l = []
        g = None
        for i in range(n):
            g = self.taylor_term(i, self.args[0], g)
            g = g.nseries(x, n=n)
            l.append(g.removeO())
        return Add(*l)

    def _eval_as_leading_term(self, x, logx, cdir):
        from sympy.calculus.util import AccumBounds
        arg = self.args[0].cancel().as_leading_term(x, logx=logx)
        arg0 = arg.subs(x, 0)
        if arg is S.NaN:
            return S.NaN
        if isinstance(arg0, AccumBounds):
            # This check addresses a corner case involving AccumBounds.
            # if isinstance(arg, AccumBounds) is True, then arg0 can either be 0,
            # AccumBounds(-oo, 0) or AccumBounds(-oo, oo).
            # Check out function: test_issue_18473() in test_exponential.py and
            # test_limits.py for more information.
            if re(cdir) < S.Zero:
                return exp(-arg0)
            return exp(arg0)
        if arg0 is S.NaN:
            arg0 = arg.limit(x, 0)
        if arg0.is_infinite is False:
            return exp(arg0)
        raise PoleError("Cannot expand %s around 0" % (self))

    def _eval_rewrite_as_sin(self, arg, **kwargs):
        from sympy.functions.elementary.trigonometric import sin
        return sin(I*arg + pi/2) - I*sin(I*arg)

    def _eval_rewrite_as_cos(self, arg, **kwargs):
        from sympy.functions.elementary.trigonometric import cos
        return cos(I*arg) + I*cos(I*arg + pi/2)

    def _eval_rewrite_as_tanh(self, arg, **kwargs):
        from sympy.functions.elementary.hyperbolic import tanh
        return (1 + tanh(arg/2))/(1 - tanh(arg/2))

    def _eval_rewrite_as_sqrt(self, arg, **kwargs):
        from sympy.functions.elementary.trigonometric import sin, cos
        if arg.is_Mul:
            coeff = arg.coeff(pi*I)
            if coeff and coeff.is_number:
                cosine, sine = cos(pi*coeff), sin(pi*coeff)
                if not isinstance(cosine, cos) and not isinstance (sine, sin):
                    return cosine + I*sine

    def _eval_rewrite_as_Pow(self, arg, **kwargs):
        if arg.is_Mul:
            logs = [a for a in arg.args if isinstance(a, log) and len(a.args) == 1]
            if logs:
                return Pow(logs[0].args[0], arg.coeff(logs[0]))

    def _eval_rewrite_as_EML(self, arg, **kwargs):
        return EML(arg, S.One, evaluate=False)


def match_real_imag(expr):
    r"""
    Try to match expr with $a + Ib$ for real $a$ and $b$.

    ``match_real_imag`` returns a tuple containing the real and imaginary
    parts of expr or ``(None, None)`` if direct matching is not possible. Contrary
    to :func:`~.re`, :func:`~.im``, and ``as_real_imag()``, this helper will not force things
    by returning expressions themselves containing ``re()`` or ``im()`` and it
    does not expand its argument either.

    """
    r_, i_ = expr.as_independent(I, as_Add=True)
    if i_ == 0 and r_.is_real:
        return (r_, i_)
    i_ = i_.as_coefficient(I)
    if i_ and i_.is_real and r_.is_real:
        return (r_, i_)
    else:
        return (None, None) # simpler to check for than None


class log(DefinedFunction):
    r"""
    The natural logarithm function `\ln(x)` or `\log(x)`.

    Explanation
    ===========

    Logarithms are taken with the natural base, `e`. To get
    a logarithm of a different base ``b``, use ``log(x, b)``,
    which is essentially short-hand for ``log(x)/log(b)``.

    ``log`` represents the principal branch of the natural
    logarithm. As such it has a branch cut along the negative
    real axis and returns values having a complex argument in
    `(-\pi, \pi]`.

    Examples
    ========

    >>> from sympy import log, sqrt, S, I
    >>> log(8, 2)
    3
    >>> log(S(8)/3, 2)
    -log(3)/log(2) + 3
    >>> log(-1 + I*sqrt(3))
    log(2) + 2*I*pi/3

    See Also
    ========

    sympy.functions.elementary.exponential.exp

    """

    args: tuple[Expr]

    _singularities = (S.Zero, S.ComplexInfinity)

    def fdiff(self, argindex=1):
        """
        Returns the first derivative of the function.
        """
        if argindex == 1:
            return 1/self.args[0]
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        r"""
        Returns `e^x`, the inverse function of `\log(x)`.
        """
        return exp

    @classmethod
    def eval(cls, arg, base=None):
        from sympy.calculus import AccumBounds
        from sympy.sets.setexpr import SetExpr

        arg = sympify(arg)

        if base is not None:
            base = sympify(base)
            if base == 1:
                if arg == 1:
                    return S.NaN
                else:
                    return S.ComplexInfinity
            try:
                # handle extraction of powers of the base now
                # or else expand_log in Mul would have to handle this
                n = multiplicity(base, arg)
                if n:
                    return n + log(arg / base**n) / log(base)
                else:
                    return log(arg)/log(base)
            except ValueError:
                pass
            if base is not S.Exp1:
                return cls(arg)/cls(base)
            else:
                return cls(arg)

        if arg.is_Number:
            if arg.is_zero:
                return S.ComplexInfinity
            elif arg is S.One:
                return S.Zero
            elif arg is S.Infinity or arg is S.NegativeInfinity:
                return S.Infinity
            elif arg is S.NaN:
                return S.NaN
            elif arg.is_Rational and arg.p == 1:
                return -cls(arg.q)

        if arg.is_Pow and arg.base is S.Exp1 and arg.exp.is_extended_real:
            return arg.exp
        if isinstance(arg, exp) and arg.exp.is_extended_real:
            return arg.exp
        elif isinstance(arg, exp) and arg.exp.is_number:
            r_, i_ = match_real_imag(arg.exp)
            if i_ and i_.is_comparable:
                i_ %= 2*pi
                if i_ > pi:
                    i_ -= 2*pi
                return r_ + expand_mul(i_ * I, deep=False)
        elif isinstance(arg, exp_polar):
            return unpolarify(arg.exp)
        elif isinstance(arg, AccumBounds):
            if arg.min.is_positive:
                return AccumBounds(log(arg.min), log(arg.max))
            elif arg.min.is_zero:
                return AccumBounds(S.NegativeInfinity, log(arg.max))
            else:
                return S.NaN
        elif isinstance(arg, SetExpr):
            return arg._eval_func(cls)

        if arg.is_number:
            if arg.is_negative:
                return pi * I + cls(-arg)
            elif arg is S.ComplexInfinity:
                return S.ComplexInfinity
            elif arg is S.Exp1:
                return S.One

        if arg.is_zero:
            return S.ComplexInfinity

        # don't autoexpand Pow or Mul (see the issue 3351):
        if not arg.is_Add:
            coeff = arg.as_coefficient(I)

            if coeff is not None:
                if coeff is S.Infinity or coeff is S.NegativeInfinity:
                    return S.Infinity
                elif coeff.is_Rational:
                    if coeff.is_nonnegative:
                        return pi * I * S.Half + cls(coeff)
                    else:
                        return -pi * I * S.Half + cls(-coeff)

        if arg.is_number and arg.is_algebraic:
            # Match arg = coeff*(r_ + i_*I) with coeff>0, r_ and i_ real.
            coeff, arg_ = arg.as_independent(I, as_Add=False)
            if coeff.is_negative:
                coeff *= -1
                arg_ *= -1
            arg_ = expand_mul(arg_, deep=False)
            r_, i_ = arg_.as_independent(I, as_Add=True)
            i_ = i_.as_coefficient(I)
            if coeff.is_real and i_ and i_.is_real and r_.is_real:
                if r_.is_zero:
                    if i_.is_positive:
                        return pi * I * S.Half + cls(coeff * i_)
                    elif i_.is_negative:
                        return -pi * I * S.Half + cls(coeff * -i_)
                else:
                    from sympy.simplify import ratsimp
                    # Check for arguments involving rational multiples of pi
                    t = (i_/r_).cancel()
                    t1 = (-t).cancel()
                    atan_table = _log_atan_table()
                    if t in atan_table:
                        modulus = ratsimp(coeff * Abs(arg_))
                        if r_.is_positive:
                            return cls(modulus) + I * atan_table[t]
                        else:
                            return cls(modulus) + I * (atan_table[t] - pi)
                    elif t1 in atan_table:
                        modulus = ratsimp(coeff * Abs(arg_))
                        if r_.is_positive:
                            return cls(modulus) + I * (-atan_table[t1])
                        else:
                            return cls(modulus) + I * (pi - atan_table[t1])

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):  # of log(1+x)
        r"""
        Returns the next term in the Taylor series expansion of `\log(1+x)`.
        """
        from sympy.simplify.powsimp import powsimp
        if n < 0:
            return S.Zero
        x = sympify(x)
        if n == 0:
            return x
        if previous_terms:
            p = previous_terms[-1]
            if p is not None:
                return powsimp((-n) * p * x / (n + 1), deep=True, combine='exp')
        return (1 - 2*(n % 2)) * x**(n + 1)/(n + 1)

    def _eval_expand_log(self, deep=True, **hints):
        from sympy.concrete import Sum, Product
        force = hints.get('force', False)
        factor = hints.get('factor', False)
        if (len(self.args) == 2):
            return expand_log(self.func(*self.args), deep=deep, force=force)
        arg = self.args[0]
        if arg.is_Integer:
            # remove perfect powers
            p = perfect_power(arg)
            logarg = None
            coeff = 1
            if p is not False:
                arg, coeff = p
                logarg = self.func(arg)
            # expand as product of its prime factors if factor=True
            if factor:
                p = factorint(arg)
                if arg not in p.keys():
                    logarg = sum(n*log(val) for val, n in p.items())
            if logarg is not None:
                return coeff*logarg
        elif arg.is_Rational:
            return log(arg.p) - log(arg.q)
        elif arg.is_Mul:
            expr = []
            nonpos = []
            for x in arg.args:
                if force or x.is_positive or x.is_polar:
                    a = self.func(x)
                    if isinstance(a, log):
                        expr.append(self.func(x)._eval_expand_log(**hints))
                    else:
                        expr.append(a)
                elif x.is_negative:
                    a = self.func(-x)
                    expr.append(a)
                    nonpos.append(S.NegativeOne)
                else:
                    nonpos.append(x)
            return Add(*expr) + log(Mul(*nonpos))
        elif arg.is_Pow or isinstance(arg, exp):
            if force or (arg.exp.is_extended_real and (arg.base.is_positive or ((arg.exp+1)
                .is_positive and (arg.exp-1).is_nonpositive))) or arg.base.is_polar:
                b = arg.base
                e = arg.exp
                a = self.func(b)
                if isinstance(a, log):
                    return unpolarify(e) * a._eval_expand_log(**hints)
                else:
                    return unpolarify(e) * a
        elif isinstance(arg, Product):
            if force or arg.function.is_positive:
                return Sum(log(arg.function), *arg.limits)

        return self.func(arg)

    def _eval_simplify(self, **kwargs):
        from sympy.simplify.simplify import expand_log, simplify, inversecombine
        if len(self.args) == 2:  # it's unevaluated
            return simplify(self.func(*self.args), **kwargs)

        expr = self.func(simplify(self.args[0], **kwargs))
        if kwargs['inverse']:
            expr = inversecombine(expr)
        expr = expand_log(expr, deep=True)
        return min([expr, self], key=kwargs['measure'])

    def as_real_imag(self, deep=True, **hints):
        """
        Returns this function as a complex coordinate.

        Examples
        ========

        >>> from sympy import I, log
        >>> from sympy.abc import x
        >>> log(x).as_real_imag()
        (log(Abs(x)), arg(x))
        >>> log(I).as_real_imag()
        (0, pi/2)
        >>> log(1 + I).as_real_imag()
        (log(sqrt(2)), pi/4)
        >>> log(I*x).as_real_imag()
        (log(Abs(x)), arg(I*x))

        """
        sarg = self.args[0]
        if deep:
            sarg = self.args[0].expand(deep, **hints)
        sarg_abs = Abs(sarg)
        if sarg_abs == sarg:
            return self, S.Zero
        sarg_arg = arg(sarg)
        if hints.get('log', False):  # Expand the log
            hints['complex'] = False
            return (log(sarg_abs).expand(deep, **hints), sarg_arg)
        else:
            return log(sarg_abs), sarg_arg

    def _eval_is_rational(self):
        s = self.func(*self.args)
        if s.func == self.func:
            if (self.args[0] - 1).is_zero:
                return True
            if s.args[0].is_rational and fuzzy_not((self.args[0] - 1).is_zero):
                return False
        else:
            return s.is_rational

    def _eval_is_algebraic(self):
        s = self.func(*self.args)
        if s.func == self.func:
            if (self.args[0] - 1).is_zero:
                return True
            elif fuzzy_not((self.args[0] - 1).is_zero):
                if self.args[0].is_algebraic:
                    return False
        else:
            return s.is_algebraic

    def _eval_is_extended_real(self):
        return self.args[0].is_extended_positive

    def _eval_is_complex(self):
        z = self.args[0]
        return fuzzy_and([z.is_complex, fuzzy_not(z.is_zero)])

    def _eval_is_finite(self):
        arg = self.args[0]
        if arg.is_zero:
            return False
        return arg.is_finite

    def _eval_is_extended_positive(self):
        return (self.args[0] - 1).is_extended_positive

    def _eval_is_zero(self):
        return (self.args[0] - 1).is_zero

    def _eval_is_extended_nonnegative(self):
        return (self.args[0] - 1).is_extended_nonnegative

    def _eval_nseries(self, x, n, logx, cdir=0):
        # NOTE Please see the comment at the beginning of this file, labelled
        #      IMPORTANT.
        from sympy.series.order import Order
        from sympy.simplify.simplify import logcombine
        from sympy.core.symbol import Dummy

        if self.args[0] == x:
            return log(x) if logx is None else logx
        arg = self.args[0]
        t = Dummy('t', positive=True)
        if cdir == 0:
            cdir = 1
        z = arg.subs(x, cdir*t)

        k, l = Wild("k"), Wild("l")
        r = z.match(k*t**l)
        if r is not None:
            k, l = r[k], r[l]
            if l != 0 and not l.has(t) and not k.has(t):
                r = l*log(x) if logx is None else l*logx
                r += log(k) - l*log(cdir) # XXX true regardless of assumptions?
                return r

        def coeff_exp(term, x):
            coeff, exp = S.One, S.Zero
            for factor in Mul.make_args(term):
                if factor.has(x):
                    base, exp = factor.as_base_exp()
                    if base != x:
                        try:
                            return term.leadterm(x)
                        except ValueError:
                            return term, S.Zero
                else:
                    coeff *= factor
            return coeff, exp

        # TODO new and probably slow
        try:
            a, b = z.leadterm(t, logx=logx, cdir=1)
        except (ValueError, NotImplementedError, PoleError):
            s = z._eval_nseries(t, n=n, logx=logx, cdir=1)
            while s.is_Order:
                n += 1
                s = z._eval_nseries(t, n=n, logx=logx, cdir=1)
            try:
                a, b = s.removeO().leadterm(t, cdir=1)
            except ValueError:
                a, b = s.removeO().as_leading_term(t, cdir=1), S.Zero

        p = (z/(a*t**b) - 1).cancel()._eval_nseries(t, n=n, logx=logx, cdir=1)
        if p.has(exp):
            p = logcombine(p)
        if isinstance(p, Order):
            n = p.getn()
        _, d = coeff_exp(p, t)
        logx = log(x) if logx is None else logx

        if not d.is_positive:
            res = log(a) - b*log(cdir) + b*logx
            _res = res
            logflags = {"deep": True, "log": True, "mul": False, "power_exp": False,
                "power_base": False, "multinomial": False, "basic": False, "force": True,
                "factor": False}
            expr = self.expand(**logflags)
            if (not a.could_extract_minus_sign() and
                logx.could_extract_minus_sign()):
                _res = _res.subs(-logx, -log(x)).expand(**logflags)
            else:
                _res = _res.subs(logx, log(x)).expand(**logflags)
            if _res == expr:
                return res
            return res + Order(x**n, x)

        def mul(d1, d2):
            res = {}
            for e1, e2 in product(d1, d2):
                ex = e1 + e2
                if ex < n:
                    res[ex] = res.get(ex, S.Zero) + d1[e1]*d2[e2]
            return res

        pterms = {}

        for term in Add.make_args(p.removeO()):
            co1, e1 = coeff_exp(term, t)
            pterms[e1] = pterms.get(e1, S.Zero) + co1

        k = S.One
        terms = {}
        pk = pterms

        while k*d < n:
            coeff = -S.NegativeOne**k/k
            for ex in pk:
                terms[ex] = terms.get(ex, S.Zero) + coeff*pk[ex]
            pk = mul(pk, pterms)
            k += S.One

        res = log(a) - b*log(cdir) + b*logx
        for ex in terms:
            res += terms[ex].cancel()*t**(ex)

        if a.is_negative and im(z) != 0:
            from sympy.functions.special.delta_functions import Heaviside
            for i, term in enumerate(z.lseries(t)):
                if not term.is_real or i == 5:
                    break
            if i < 5:
                coeff, _ = term.as_coeff_exponent(t)
                res += -2*I*pi*Heaviside(-im(coeff), 0)

        res = res.subs(t, x/cdir)
        return res + Order(x**n, x)

    def _eval_as_leading_term(self, x, logx, cdir):
        # NOTE
        # Refer https://github.com/sympy/sympy/pull/23592 for more information
        # on each of the following steps involved in this method.
        arg0 = self.args[0].together()

        # STEP 1
        t = Dummy('t', positive=True)
        if cdir == 0:
            cdir = 1
        z = arg0.subs(x, cdir*t)

        # STEP 2
        try:
            c, e = z.leadterm(t, logx=logx, cdir=1)
        except ValueError:
            arg = arg0.as_leading_term(x, logx=logx, cdir=cdir)
            return log(arg)
        if c.has(t):
            c = c.subs(t, x/cdir)
            if e != 0:
                raise PoleError("Cannot expand %s around 0" % (self))
            return log(c)

        # STEP 3
        if c == S.One and e == S.Zero:
            return (arg0 - S.One).as_leading_term(x, logx=logx)

        # STEP 4
        res = log(c) - e*log(cdir)
        logx = log(x) if logx is None else logx
        res += e*logx

        # STEP 5
        if c.is_negative and im(z) != 0:
            from sympy.functions.special.delta_functions import Heaviside
            for i, term in enumerate(z.lseries(t)):
                if not term.is_real or i == 5:
                    break
            if i < 5:
                coeff, _ = term.as_coeff_exponent(t)
                res += -2*I*pi*Heaviside(-im(coeff), 0)
        return res

    def _eval_derivative_n_times(self, s, n):
        if self.args[0] == s and n.is_integer and n.is_positive:
            return S.NegativeOne**(n-1) * factorial(n - 1) / s**n
        return super()._eval_derivative_n_times(s, n)

    def _eval_rewrite_as_EML(self, arg, **kwargs):
        # Pure-EML grammar (Odrzywolek 2026): the result is a binary tree whose
        # only constant leaf is 1 and whose only internal node is EML, namely
        #
        #     log(x) = EML(1, EML(EML(1, x), 1)).
        #
        # Each EML below is constructed with ``evaluate=False`` to preserve
        # the tree structure; otherwise the inner ``EML(_, 1)`` would auto-
        # collapse to ``exp(_)`` and the recursive form would be lost.
        inner = EML(S.One, arg, evaluate=False)
        middle = EML(inner, S.One, evaluate=False)
        return EML(S.One, middle, evaluate=False)


class EML(DefinedFunction):
    r"""
    The two-argument function ``EML(x, y) = exp(x) - log(y)``.

    Explanation
    ===========

    ``EML`` combines the exponential and (natural) logarithm into a single
    binary primitive. It is defined by

    .. math::

        \operatorname{EML}(x, y) = e^{x} - \log(y),

    where :math:`\log` denotes the principal branch of the natural logarithm.

    Following Odrzywolek (arXiv:2603.21852), every elementary expression
    can be written as a binary tree whose internal nodes are ``EML`` and whose
    only constant leaf is ``1`` (together with the free symbols of the
    expression). The recursive grammar of these trees is

    .. math::

        S \;\to\; 1 \;\mid\; \operatorname{EML}(S, S).

    The conversion is obtained with ``expr.rewrite(EML)``. The key building
    blocks given by the paper are

    .. math::

        e^{x} = \operatorname{EML}(x, 1), \qquad
        \log(x) = \operatorname{EML}(1, \operatorname{EML}(\operatorname{EML}(1, x), 1)).

    Examples
    ========

    >>> from sympy import EML, exp, log, symbols, simplify
    >>> x, y = symbols('x y')

    ``EML`` does not auto-evaluate into ``exp``/``log`` form; the conversion
    is explicit, through ``rewrite``, ``expand_func`` or :func:`from_eml`:

    >>> EML(x, 1)
    EML(x, 1)
    >>> EML(x, 1).rewrite(exp)
    exp(x)
    >>> EML(0, 1).rewrite(exp)
    1
    >>> EML(1, 1).rewrite(exp)
    E

    Round-tripping through ``exp`` and ``log``:

    >>> EML(x, y).rewrite(exp)
    exp(x) - log(y)
    >>> EML(x, y).rewrite(log)
    exp(x) - log(y)

    Differentiation in either argument:

    >>> EML(x, y).diff(x)
    exp(x)
    >>> EML(x, y).diff(y)
    -1/y

    Rewriting elementary expressions in pure-EML recursive form:

    >>> exp(x).rewrite(EML)
    EML(x, 1)
    >>> log(x).rewrite(EML)
    EML(1, EML(EML(1, x), 1))

    The recursive form round-trips back to the original (for positive ``x``):

    >>> from sympy import Symbol
    >>> xp = Symbol('x', positive=True)
    >>> simplify(log(xp).rewrite(EML).rewrite(exp) - log(xp))
    0

    See Also
    ========

    sympy.functions.elementary.exponential.exp
    sympy.functions.elementary.exponential.log
    """

    @property
    def exp_arg(self):
        """Argument that enters the function under ``exp``."""
        return self.args[0]

    @property
    def log_arg(self):
        """Argument that enters the function under ``log``."""
        return self.args[1]

    @classmethod
    def eval(cls, x, y):
        # ``EML`` deliberately does not auto-evaluate into ``exp``/``log``
        # form.  That conversion is explicit, via ``rewrite``, ``expand_func``
        # or :func:`from_eml`, so that the pure-EML tree structure is
        # preserved (e.g. ``EML(x, 1)`` stays ``EML(x, 1)`` rather than
        # collapsing to ``exp(x)``).  Only genuine ``NaN`` results are folded
        # here -- ``exp`` is ``NaN`` at ``zoo``/``NaN`` and ``log`` is ``NaN``
        # at ``NaN`` -- so that invalid expressions do not survive as live
        # ``EML`` nodes.
        if x is S.NaN or x is S.ComplexInfinity or y is S.NaN:
            return S.NaN

    def fdiff(self, argindex=1):
        """
        Returns the first derivative of this function with respect to the
        argument selected by *argindex*.
        """
        x, y = self.args
        if argindex == 1:
            return exp(x)
        elif argindex == 2:
            return -1/y
        raise ArgumentIndexError(self, argindex)

    def _eval_derivative_n_times(self, s, n):
        """
        Closed-form nth derivative when *s* is exactly one of the arguments
        and the other does not depend on it.
        """
        x, y = self.args
        if n.is_integer and n.is_positive:
            if x == s and not y.has(s):
                return exp(x)
            if y == s and not x.has(s):
                return S.NegativeOne**n * factorial(n - 1) / y**n
        return super()._eval_derivative_n_times(s, n)

    def _eval_rewrite_as_exp(self, x, y, **kwargs):
        return exp(x) - log(y)

    def _eval_rewrite_as_log(self, x, y, **kwargs):
        return exp(x) - log(y)

    def _eval_rewrite_as_sinh(self, x, y, **kwargs):
        from sympy.functions.elementary.hyperbolic import sinh, cosh
        return sinh(x) + cosh(x) - log(y)

    def _eval_rewrite_as_cosh(self, x, y, **kwargs):
        from sympy.functions.elementary.hyperbolic import sinh, cosh
        return cosh(x) + sinh(x) - log(y)

    def _eval_rewrite_as_Pow(self, x, y, **kwargs):
        return S.Exp1**x - log(y)

    def _eval_rewrite_as_tractable(self, x, y, limitvar=None, **kwargs):
        return exp(x) - log(y)

    def _eval_expand_func(self, **hints):
        x, y = self.args
        return exp(x) - log(y)

    def _eval_conjugate(self):
        x, y = self.args
        return self.func(x.conjugate(), y.conjugate())

    def as_real_imag(self, deep=True, **hints):
        """
        Return the real and imaginary parts of ``EML(x, y)``.

        The decomposition follows from ``exp(x) - log(y)`` together with the
        standard real/imag splits of ``exp`` and ``log``.

        Examples
        ========

        >>> from sympy import EML, symbols
        >>> x, y = symbols('x y', real=True)
        >>> EML(x, y).as_real_imag()
        (exp(x) - log(Abs(y)), -arg(y))
        """
        x, y = self.args
        return (exp(x) - log(y)).as_real_imag(deep=deep, **hints)

    def _eval_as_leading_term(self, x, logx=None, cdir=0):
        xa, ya = self.args
        return (exp(xa) - log(ya)).as_leading_term(x, logx=logx, cdir=cdir)

    def _eval_nseries(self, x, n, logx, cdir=0):
        xa, ya = self.args
        return (exp(xa) - log(ya))._eval_nseries(x, n=n, logx=logx, cdir=cdir)

    def _eval_is_extended_real(self):
        x, y = self.args
        return fuzzy_and([x.is_extended_real, y.is_extended_positive])

    def _eval_is_extended_positive(self):
        x, y = self.args
        if x.is_extended_real and (y - 1).is_extended_nonpositive and y.is_extended_positive:
            # exp(x) > 0 always; -log(y) >= 0 when 0 < y <= 1.
            return True
        # exp(x) > log(y) iff y < exp(exp(x)); hard to test in general.
        return None

    def _eval_is_extended_negative(self):
        x, y = self.args
        if x is S.NegativeInfinity and (y - 1).is_extended_positive:
            # -log(y) < 0 when y > 1.
            return True
        return None

    def _eval_is_zero(self):
        # EML(x, y) = 0 iff log(y) = exp(x), i.e. y = exp(exp(x)).
        # Handle the easy case of two concrete numbers via evalf-free equality.
        x, y = self.args
        if x.is_extended_real and y.is_extended_real:
            # Detect the only real-real solution: y == exp(exp(x)).
            if y == exp(exp(x)):
                return True
            if (y - exp(exp(x))).is_extended_nonzero:
                return False
        return None

    def _eval_is_finite(self):
        x, y = self.args
        return fuzzy_and([x.is_finite, y.is_finite, fuzzy_not(y.is_zero)])

    def _eval_is_complex(self):
        x, y = self.args
        return fuzzy_and([x.is_complex, y.is_complex, fuzzy_not(y.is_zero)])

    def _eval_is_rational(self):
        # exp(x) - log(y) is rational only in degenerate cases
        # (e.g. EML(0, 1) = 1). Otherwise transcendental.
        x, y = self.args
        if x.is_zero and (y - 1).is_zero:
            return True
        if (x.is_rational and fuzzy_not(x.is_zero)) or \
           (y.is_rational and fuzzy_not((y - 1).is_zero)):
            return False
        return None

    def _eval_is_algebraic(self):
        # By Lindemann-Weierstrass, exp(x) is transcendental for nonzero
        # algebraic x; log(y) is transcendental for algebraic y != 1.
        x, y = self.args
        if x.is_zero and (y - 1).is_zero:
            return True
        if (x.is_algebraic and fuzzy_not(x.is_zero)) or \
           (y.is_algebraic and fuzzy_not((y - 1).is_zero)):
            return False
        return None

    def _eval_evalf(self, prec):
        x, y = self.args
        return (exp(x) - log(y))._eval_evalf(prec)


def to_eml(expr, deep=True, numbers=False):
    r"""Rewrite the elementary functions in *expr* using the ``EML`` primitive.

    ``EML(x, y) = exp(x) - log(y)`` is the binary primitive of Odrzywolek
    (arXiv:2603.21852).  ``to_eml`` is the public entry point for the
    ``expr.rewrite(EML)`` machinery: it walks *expr* and rewrites every
    elementary function it understands -- ``exp``, ``log``, powers, and the
    trigonometric, hyperbolic and inverse functions -- into ``EML`` form.

    The surrounding arithmetic structure (``Add``, ``Mul``, ``Pow`` with the
    base/exponent themselves converted) is preserved: ``EML`` generates
    ``exp``, ``log`` and subtraction, but not the addition or multiplication of
    independent terms, so an expression collapses to ``EML`` *nodes* embedded
    in its ordinary arithmetic, not to a single ``EML`` tree.

    Numeric values are kept as leaves by default.  The only constant the
    grammar singles out is ``1``; on top of it ``EML(1, 1) = e`` is the one
    numeric value with a finite ``EML`` encoding, which *numbers=True*
    substitutes in.  Other integers and rationals (``2``, ``3``, ``1/2``, ...)
    have no finite ``EML`` encoding -- building them would require addition,
    which ``EML`` composition cannot express -- so they are left untouched.

    Parameters
    ==========

    expr : Expr or iterable of Expr
        The expression (or a list/tuple of expressions) to convert.
    deep : bool, optional
        If ``True`` (the default) the rewrite descends into all
        subexpressions; if ``False`` only the top level is rewritten.
    numbers : bool, optional
        If ``True``, replace Euler's number ``e`` by its ``EML`` encoding
        ``EML(1, 1)``.  Defaults to ``False``.

    Examples
    ========

    >>> from sympy import to_eml, exp, log, E, symbols
    >>> x, y = symbols('x y', positive=True)

    Elementary functions, including powers, are rewritten:

    >>> to_eml(exp(x))
    EML(x, 1)
    >>> to_eml(log(x))
    EML(1, EML(EML(1, x), 1))
    >>> to_eml(x**y)
    EML(y*EML(1, EML(EML(1, x), 1)), 1)

    Arithmetic structure is kept while each function becomes ``EML``:

    >>> to_eml(exp(x) + log(y))
    EML(1, EML(EML(1, y), 1)) + EML(x, 1)

    A list (or tuple) of expressions is converted element-wise:

    >>> to_eml([exp(x), log(x)])
    [EML(x, 1), EML(1, EML(EML(1, x), 1))]

    With ``numbers=True`` the constant ``e`` is encoded as ``EML(1, 1)``:

    >>> to_eml(E, numbers=True)
    EML(1, 1)

    See Also
    ========

    EML
    sympy.functions.elementary.exponential.exp
    sympy.functions.elementary.exponential.log
    """
    from sympy.core.basic import Basic
    from sympy.utilities.iterables import iterable
    if not isinstance(expr, Basic) and iterable(expr):
        return type(expr)(to_eml(e, deep=deep, numbers=numbers) for e in expr)
    result = sympify(expr).rewrite(EML, deep=deep)
    if numbers:
        result = result.xreplace({S.Exp1: EML(S.One, S.One, evaluate=False)})
    return result


def from_eml(expr):
    r"""Expand every ``EML`` node in *expr* back to ``exp(x) - log(y)`` form.

    ``from_eml`` is the inverse of :func:`to_eml`: it replaces each
    ``EML(x, y)`` by ``exp(x) - log(y)`` (recursively), recovering an ordinary
    expression in terms of ``exp`` and ``log``.  It is a thin, named wrapper
    around ``expr.rewrite(exp)`` that also accepts lists and tuples.

    Round-tripping ``to_eml`` then ``from_eml`` returns a mathematically
    equivalent expression (it is literally an identity on the EML nodes; only
    the usual ``exp``/``log`` simplifications such as ``log(1) = 0`` are
    applied).

    Parameters
    ==========

    expr : Expr or iterable of Expr
        The expression (or a list/tuple of expressions) to expand.

    Examples
    ========

    >>> from sympy import EML, to_eml, from_eml, exp, log, simplify, symbols
    >>> x, y = symbols('x y', positive=True)

    >>> from_eml(EML(x, y))
    exp(x) - log(y)
    >>> from_eml(to_eml(exp(x)))
    exp(x)

    No ``EML`` nodes remain; the result equals the original up to the branch
    simplifications SymPy applies (so a ``log`` round-trip matches only after
    ``simplify`` for positive arguments):

    >>> simplify(from_eml(to_eml(log(x))) - log(x))
    0

    A list (or tuple) is expanded element-wise:

    >>> from_eml([EML(x, 1, evaluate=False), EML(1, y, evaluate=False)])
    [exp(x), E - log(y)]

    See Also
    ========

    to_eml
    EML
    """
    from sympy.core.basic import Basic
    from sympy.utilities.iterables import iterable
    if not isinstance(expr, Basic) and iterable(expr):
        return type(expr)(from_eml(e) for e in expr)
    expr = sympify(expr)
    # Convert *only* ``EML`` nodes back to ``exp(x) - log(y)`` (bottom-up),
    # leaving every other function untouched -- unlike ``rewrite(exp)``, which
    # would also expand e.g. ``sin`` into its exponential form.
    expr = expr.replace(
        lambda e: isinstance(e, EML),
        lambda e: exp(e.args[0]) - log(e.args[1]))
    # The pure-EML encoding of ``log`` introduces nested ``log(exp(...))``
    # terms; ``expand_log`` folds them so the recursive log/exp grammar
    # collapses back to its closed form (e.g. ``EML(1, EML(EML(1, x), 1))``
    # becomes ``log(x)``).
    return expand_log(expr, force=True)


class LambertW(DefinedFunction):
    r"""
    The Lambert W function $W(z)$ is defined as the inverse
    function of $w \exp(w)$ [1]_.

    Explanation
    ===========

    In other words, the value of $W(z)$ is such that $z = W(z) \exp(W(z))$
    for any complex number $z$.  The Lambert W function is a multivalued
    function with infinitely many branches $W_k(z)$, indexed by
    $k \in \mathbb{Z}$.  Each branch gives a different solution $w$
    of the equation $z = w \exp(w)$.

    The Lambert W function has two partially real branches: the
    principal branch ($k = 0$) is real for real $z > -1/e$, and the
    $k = -1$ branch is real for $-1/e < z < 0$. All branches except
    $k = 0$ have a logarithmic singularity at $z = 0$.

    Examples
    ========

    >>> from sympy import LambertW
    >>> LambertW(1.2)
    0.635564016364870
    >>> LambertW(1.2, -1).n()
    -1.34747534407696 - 4.41624341514535*I
    >>> LambertW(-1).is_real
    False

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Lambert_W_function
    """
    _singularities = (-Pow(S.Exp1, -1, evaluate=False), S.ComplexInfinity)

    @classmethod
    def eval(cls, x, k=None):
        if k == S.Zero:
            return cls(x)
        elif k is None:
            k = S.Zero

        if k.is_zero:
            if x.is_zero:
                return S.Zero
            if x is S.Exp1:
                return S.One
            if x == -1/S.Exp1:
                return S.NegativeOne
            if x == -log(2)/2:
                return -log(2)
            if x == 2*log(2):
                return log(2)
            if x == -pi/2:
                return I*pi/2
            if x == exp(1 + S.Exp1):
                return S.Exp1
            if x is S.Infinity:
                return S.Infinity

        if fuzzy_not(k.is_zero):
            if x.is_zero:
                return S.NegativeInfinity
        if k is S.NegativeOne:
            if x == -pi/2:
                return -I*pi/2
            elif x == -1/S.Exp1:
                return S.NegativeOne
            elif x == -2*exp(-2):
                return -Integer(2)

    def fdiff(self, argindex=1):
        """
        Return the first derivative of this function.
        """
        x = self.args[0]

        if len(self.args) == 1:
            if argindex == 1:
                return LambertW(x)/(x*(1 + LambertW(x)))
        else:
            k = self.args[1]
            if argindex == 1:
                return LambertW(x, k)/(x*(1 + LambertW(x, k)))

        raise ArgumentIndexError(self, argindex)

    def _eval_is_extended_real(self):
        x = self.args[0]
        if len(self.args) == 1:
            k = S.Zero
        else:
            k = self.args[1]
        if k.is_zero:
            if (x + 1/S.Exp1).is_positive:
                return True
            elif (x + 1/S.Exp1).is_nonpositive:
                return False
        elif (k + 1).is_zero:
            if x.is_negative and (x + 1/S.Exp1).is_positive:
                return True
            elif x.is_nonpositive or (x + 1/S.Exp1).is_nonnegative:
                return False
        elif fuzzy_not(k.is_zero) and fuzzy_not((k + 1).is_zero):
            if x.is_extended_real:
                return False

    def _eval_is_finite(self):
        return self.args[0].is_finite

    def _eval_is_algebraic(self):
        s = self.func(*self.args)
        if s.func == self.func:
            if fuzzy_not(self.args[0].is_zero) and self.args[0].is_algebraic:
                return False
        else:
            return s.is_algebraic

    def _eval_as_leading_term(self, x, logx, cdir):
        if len(self.args) == 1:
            arg = self.args[0]
            arg0 = arg.subs(x, 0).cancel()
            if not arg0.is_zero:
                return self.func(arg0)
            return arg.as_leading_term(x)

    def _eval_nseries(self, x, n, logx, cdir=0):
        if len(self.args) == 1:
            from sympy.functions.elementary.integers import ceiling
            from sympy.series.order import Order
            arg = self.args[0].nseries(x, n=n, logx=logx)
            lt = arg.as_leading_term(x, logx=logx)
            lte = 1
            if lt.is_Pow:
                lte = lt.exp
            if ceiling(n/lte) >= 1:
                s = Add(*[(-S.One)**(k - 1)*Integer(k)**(k - 2)/
                          factorial(k - 1)*arg**k for k in range(1, ceiling(n/lte))])
                s = expand_multinomial(s)
            else:
                s = S.Zero

            return s + Order(x**n, x)
        return super()._eval_nseries(x, n, logx)

    def _eval_is_zero(self):
        x = self.args[0]
        if len(self.args) == 1:
            return x.is_zero
        else:
            return fuzzy_and([x.is_zero, self.args[1].is_zero])


@cacheit
def _log_atan_table():
    return {
        # first quadrant only
        sqrt(3): pi / 3,
        1: pi / 4,
        sqrt(5 - 2 * sqrt(5)): pi / 5,
        sqrt(2) * sqrt(5 - sqrt(5)) / (1 + sqrt(5)): pi / 5,
        sqrt(5 + 2 * sqrt(5)): pi * Rational(2, 5),
        sqrt(2) * sqrt(sqrt(5) + 5) / (-1 + sqrt(5)): pi * Rational(2, 5),
        sqrt(3) / 3: pi / 6,
        sqrt(2) - 1: pi / 8,
        sqrt(2 - sqrt(2)) / sqrt(sqrt(2) + 2): pi / 8,
        sqrt(2) + 1: pi * Rational(3, 8),
        sqrt(sqrt(2) + 2) / sqrt(2 - sqrt(2)): pi * Rational(3, 8),
        sqrt(1 - 2 * sqrt(5) / 5): pi / 10,
        (-sqrt(2) + sqrt(10)) / (2 * sqrt(sqrt(5) + 5)): pi / 10,
        sqrt(1 + 2 * sqrt(5) / 5): pi * Rational(3, 10),
        (sqrt(2) + sqrt(10)) / (2 * sqrt(5 - sqrt(5))): pi * Rational(3, 10),
        2 - sqrt(3): pi / 12,
        (-1 + sqrt(3)) / (1 + sqrt(3)): pi / 12,
        2 + sqrt(3): pi * Rational(5, 12),
        (1 + sqrt(3)) / (-1 + sqrt(3)): pi * Rational(5, 12)
    }
