from __future__ import print_function, division

from sympy.core import C, sympify
from sympy.core.add import Add
from sympy.core.function import Lambda, Function, ArgumentIndexError
from sympy.core.cache import cacheit
from sympy.core.singleton import S
from sympy.core.symbol import Wild, Dummy
from sympy.core.mul import Mul

from sympy.functions.elementary.miscellaneous import sqrt
from sympy.ntheory import multiplicity, perfect_power
from sympy.core.compatibility import xrange

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


class ExpBase(Function):

    unbranched = True

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

        >>> from sympy.functions import exp
        >>> from sympy.abc import x
        >>> exp(-x).as_numer_denom()
        (1, exp(x))
        >>> exp(x).as_numer_denom()
        (exp(x), 1)
        """
        # this should be the same as Pow.as_numer_denom wrt
        # exponent handling
        exp = self.exp
        neg_exp = exp.is_negative
        if not neg_exp and not (-exp).is_negative:
            neg_exp = _coeff_isneg(exp)
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

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def _eval_is_bounded(self):
        arg = self.args[0]
        if arg.is_unbounded:
            if arg.is_negative:
                return True
            if arg.is_positive:
                return False
        if arg.is_bounded:
            return True

    def _eval_is_rational(self):
        s = self.func(*self.args)
        if s.func == self.func:
            if s.args[0].is_rational:
                return False
        else:
            return s.is_rational

    def _eval_is_zero(self):
        return (self.args[0] is S.NegativeInfinity)

    def _eval_power(b, e):
        """exp(arg)**e -> exp(arg*e) if assumptions allow it.
        """
        f = b.func
        be = b.exp
        rv = f(be*e)
        if e.is_integer:
            return rv
        if be.is_real:
            return rv
        # "is True" needed below; exp.is_polar returns <property object ...>
        if f.is_polar is True:
            return rv
        if e.is_polar:
            return rv
        if be.is_polar:
            return rv
        besmall = abs(be) <= S.Pi
        if besmall is True:
            return rv
        elif besmall is False and e.is_Rational and e.q == 2:
            return -rv

    def _eval_expand_power_exp(self, **hints):
        arg = self.args[0]
        if arg.is_Add and arg.is_commutative:
            expr = 1
            for x in arg.args:
                expr *= self.func(x)
            return expr
        return self.func(arg)


class exp_polar(ExpBase):
    r"""
    Represent a 'polar number' (see g-function Sphinx documentation).

    ``exp_polar`` represents the function
    `Exp: \mathbb{C} \rightarrow \mathcal{S}`, sending the complex number
    `z = a + bi` to the polar number `r = exp(a), \theta = b`. It is one of
    the main functions to construct polar numbers.

    >>> from sympy import exp_polar, pi, I, exp

    The main difference is that polar numbers don't "wrap around" at `2 \pi`:

    >>> exp(2*pi*I)
    1
    >>> exp_polar(2*pi*I)
    exp_polar(2*I*pi)

    apart from that they behave mostly like classical complex numbers:

    >>> exp_polar(2)*exp_polar(3)
    exp_polar(5)

    See also
    ========

    sympy.simplify.simplify.powsimp
    sympy.functions.elementary.complexes.polar_lift
    sympy.functions.elementary.complexes.periodic_argument
    sympy.functions.elementary.complexes.principal_branch
    """

    is_polar = True
    is_comparable = False  # cannot be evalf'd

    def _eval_Abs(self):
        from sympy import expand_mul
        return sqrt( expand_mul(self * self.conjugate()) )

    def _eval_evalf(self, prec):
        """ Careful! any evalf of polar numbers is flaky """
        from sympy import im, pi, re
        i = im(self.args[0])
        if i <= -pi or i > pi:
            return self  # cannot evalf for this argument
        res = exp(self.args[0])._eval_evalf(prec)
        if i > 0 and im(res) < 0:
            # i ~ pi, but exp(I*i) evaluated to argument slightly bigger than pi
            return re(res)
        return res

    def _eval_is_real(self):
        if self.args[0].is_real:
            return True

    def as_base_exp(self):
        # XXX exp_polar(0) is special!
        if self.args[0] == 0:
            return self, S(1)
        return ExpBase.as_base_exp(self)


class exp(ExpBase):
    """
    The exponential function, :math:`e^x`.

    See Also
    ========

    log
    """

    def fdiff(self, argindex=1):
        """
        Returns the first derivative of this function.
        """
        if argindex == 1:
            return self
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Zero:
                return S.One
            elif arg is S.One:
                return S.Exp1
            elif arg is S.Infinity:
                return S.Infinity
            elif arg is S.NegativeInfinity:
                return S.Zero
        elif arg.func is log:
            return arg.args[0]
        elif arg.is_Mul:
            Ioo = S.ImaginaryUnit*S.Infinity
            if arg in [Ioo, -Ioo]:
                return S.NaN

            coeff = arg.coeff(S.Pi*S.ImaginaryUnit)
            if coeff:
                if (2*coeff).is_integer:
                    if coeff.is_even:
                        return S.One
                    elif coeff.is_odd:
                        return S.NegativeOne
                    elif (coeff + S.Half).is_even:
                        return -S.ImaginaryUnit
                    elif (coeff + S.Half).is_odd:
                        return S.ImaginaryUnit

            # Warning: code in risch.py will be very sensitive to changes
            # in this (see DifferentialExtension).

            # look for a single log factor

            coeff, terms = arg.as_coeff_Mul()

            # but it can't be multiplied by oo
            if coeff in [S.NegativeInfinity, S.Infinity]:
                return None

            coeffs, log_term = [coeff], None
            for term in Mul.make_args(terms):
                if term.func is log:
                    if log_term is None:
                        log_term = term.args[0]
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
            for a in arg.args:
                if a is S.One:
                    add.append(a)
                    continue
                newa = cls(a)
                if newa.func is cls:
                    add.append(a)
                else:
                    out.append(newa)
            if out:
                return Mul(*out)*cls(Add(*add), evaluate=False)

        elif arg.is_Matrix:
            from sympy import Matrix
            return arg.exp()

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
        return x**n/C.factorial()(n)

    def as_real_imag(self, deep=True, **hints):
        """
        Returns this function as a 2-tuple representing a complex number.

        Examples
        ========

        >>> from sympy import I
        >>> from sympy.abc import x
        >>> from sympy.functions import exp
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
        re, im = self.args[0].as_real_imag()
        if deep:
            re = re.expand(deep, **hints)
            im = im.expand(deep, **hints)
        cos, sin = C.cos(im), C.sin(im)
        return (exp(re)*cos, exp(re)*sin)

    def _eval_subs(self, old, new):
        arg = self.args[0]
        o = old
        if old.is_Pow:  # handle (exp(3*log(x))).subs(x**2, z) -> z**(3/2)
            o = exp(o.exp*log(o.base))
        if o.func is exp:
            # exp(a*expr) .subs( exp(b*expr), y )  ->  y ** (a/b)
            a, expr_terms = self.args[0].as_independent(
                C.Symbol, as_Add=False)
            b, expr_terms_ = o.args[0].as_independent(
                C.Symbol, as_Add=False)

            if expr_terms == expr_terms_:
                return new**(a/b)

            if arg.is_Add:  # exp(2*x+a).subs(exp(3*x),y) -> y**(2/3) * exp(a)
                # exp(exp(x) + exp(x**2)).subs(exp(exp(x)), w) -> w * exp(exp(x**2))
                oarg = o.args[0]
                new_l = []
                o_al = []
                coeff2, terms2 = oarg.as_coeff_mul()
                for a in arg.args:
                    a = a._subs(o, new)
                    coeff1, terms1 = a.as_coeff_mul()
                    if terms1 == terms2:
                        new_l.append(new**(coeff1/coeff2))
                    else:
                        o_al.append(a._subs(o, new))
                if new_l:
                    new_l.append(self.func(Add(*o_al)))
                    r = Mul(*new_l)
                    return r
        if o is S.Exp1:
            # treat this however Pow is being treated
            u = C.Dummy('u', positive=True)
            return (u**self.args[0]).xreplace({u: new})

        return Function._eval_subs(self, o, new)

    def _eval_is_real(self):
        if self.args[0].is_real:
            return True
        elif self.args[0].is_imaginary:
            arg2 = -S(2) * S.ImaginaryUnit * self.args[0] / S.Pi
            return arg2.is_even

    def _eval_is_positive(self):
        if self.args[0].is_real:
            return not self.args[0] is S.NegativeInfinity
        elif self.args[0].is_imaginary:
            arg2 = -S.ImaginaryUnit * self.args[0] / S.Pi
            return arg2.is_even

    def _eval_nseries(self, x, n, logx):
        # NOTE Please see the comment at the beginning of this file, labelled
        #      IMPORTANT.
        from sympy import limit, oo, powsimp
        arg = self.args[0]
        arg_series = arg._eval_nseries(x, n=n, logx=logx)
        if arg_series.is_Order:
            return 1 + arg_series
        arg0 = limit(arg_series.removeO(), x, 0)
        if arg0 in [-oo, oo]:
            return self
        t = Dummy("t")
        exp_series = exp(t)._taylor(t, n)
        o = exp_series.getO()
        exp_series = exp_series.removeO()
        r = exp(arg0)*exp_series.subs(t, arg_series - arg0)
        r += C.Order(o.expr.subs(t, (arg_series - arg0)), x)
        r = r.expand()
        return powsimp(r, deep=True, combine='exp')

    def _taylor(self, x, n):
        l = []
        g = None
        for i in xrange(n):
            g = self.taylor_term(i, self.args[0], g)
            g = g.nseries(x, n=n)
            l.append(g)
        return Add(*l) + C.Order(x**n, x)

    def _eval_as_leading_term(self, x):
        arg = self.args[0]
        if arg.is_Add:
            return Mul(*[exp(f).as_leading_term(x) for f in arg.args])
        arg = self.args[0].as_leading_term(x)
        if C.Order(1, x).contains(arg):
            return S.One
        return exp(arg)

    def _eval_rewrite_as_sin(self, arg):
        I = S.ImaginaryUnit
        return C.sin(I*arg + S.Pi/2) - I*C.sin(I*arg)

    def _eval_rewrite_as_cos(self, arg):
        I = S.ImaginaryUnit
        return C.cos(I*arg) + I*C.cos(I*arg + S.Pi/2)

    def _sage_(self):
        import sage.all as sage
        return sage.exp(self.args[0]._sage_())


class log(Function):
    """
    The natural logarithm function `\ln(x)` or `\log(x)`.
    Logarithms are taken with the natural base, `e`. To get
    a logarithm of a different base ``b``, use ``log(x, b)``,
    which is essentially short-hand for ``log(x)/log(b)``.

    See Also
    ========

    exp
    """

    def fdiff(self, argindex=1):
        """
        Returns the first derivative of the function.
        """
        if argindex == 1:
            return 1/self.args[0]
            s = C.Dummy('x')
            return Lambda(s**(-1), s)
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        """
        Returns `e^x`, the inverse function of `\log(x)`.
        """
        return exp

    @classmethod
    def eval(cls, arg, base=None):
        from sympy import unpolarify
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
                    den = base**n
                    if den.is_Integer:
                        return n + log(arg // den) / log(base)
                    else:
                        return n + log(arg / den) / log(base)
                else:
                    return log(arg)/log(base)
            except ValueError:
                pass
            if base is not S.Exp1:
                return cls(arg)/cls(base)
            else:
                return cls(arg)

        if arg.is_Number:
            if arg is S.Zero:
                return S.ComplexInfinity
            elif arg is S.One:
                return S.Zero
            elif arg is S.Infinity:
                return S.Infinity
            elif arg is S.NegativeInfinity:
                return S.Infinity
            elif arg is S.NaN:
                return S.NaN
            elif arg.is_negative:
                return S.Pi * S.ImaginaryUnit + cls(-arg)
            elif arg.is_Rational:
                if arg.q != 1:
                    return cls(arg.p) - cls(arg.q)
        elif arg is S.ComplexInfinity:
            return S.ComplexInfinity
        elif arg is S.Exp1:
            return S.One
        elif arg.func is exp and arg.args[0].is_real:
            return arg.args[0]
        elif arg.func is exp_polar:
            return unpolarify(arg.exp)
        #don't autoexpand Pow or Mul (see the issue 252):
        elif not arg.is_Add:
            coeff = arg.as_coefficient(S.ImaginaryUnit)

            if coeff is not None:
                if coeff is S.Infinity:
                    return S.Infinity
                elif coeff is S.NegativeInfinity:
                    return S.Infinity
                elif coeff.is_Rational:
                    if coeff.is_nonnegative:
                        return S.Pi * S.ImaginaryUnit * S.Half + cls(coeff)
                    else:
                        return -S.Pi * S.ImaginaryUnit * S.Half + cls(-coeff)

    def as_base_exp(self):
        """
        Returns this function in the form (base, exponent).
        """
        return self, S.One

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):  # of log(1+x)
        """
        Returns the next term in the Taylor series expansion of `\log(1+x)`.
        """
        from sympy import powsimp
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
        from sympy import unpolarify
        from sympy.concrete import Sum, Product
        force = hints.get('force', False)
        arg = self.args[0]
        if arg.is_Integer:
            # remove perfect powers
            p = perfect_power(int(arg))
            if p is not False:
                return p[1]*self.func(p[0])
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
        elif arg.is_Pow:
            if force or (arg.exp.is_real and arg.base.is_positive) or \
                    arg.base.is_polar:
                b = arg.base
                e = arg.exp
                a = self.func(b)
                if isinstance(a, log):
                    return unpolarify(e) * a._eval_expand_log(**hints)
                else:
                    return unpolarify(e) * a
        elif isinstance(arg, Product):
            if arg.function.is_positive:
                return Sum(log(arg.function), *arg.limits)

        return self.func(arg)

    def _eval_simplify(self, ratio, measure):
        from sympy.simplify.simplify import expand_log, logcombine, simplify
        expr = self.func(simplify(self.args[0], ratio=ratio, measure=measure))
        expr = expand_log(expr, deep=True)
        return min([expr, self], key=measure)

    def as_real_imag(self, deep=True, **hints):
        """
        Returns this function as a complex coordinate.

        Examples
        ========

        >>> from sympy import I
        >>> from sympy.abc import x
        >>> from sympy.functions import log
        >>> log(x).as_real_imag()
        (log(Abs(x)), arg(x))
        >>> log(I).as_real_imag()
        (0, pi/2)
        >>> log(1+I).as_real_imag()
        (log(sqrt(2)), pi/4)
        >>> log(I*x).as_real_imag()
        (log(Abs(x)), arg(I*x))

        """
        if deep:
            abs = C.Abs(self.args[0].expand(deep, **hints))
            arg = C.arg(self.args[0].expand(deep, **hints))
        else:
            abs = C.Abs(self.args[0])
            arg = C.arg(self.args[0])
        if hints.get('log', False):  # Expand the log
            hints['complex'] = False
            return (log(abs).expand(deep, **hints), arg)
        else:
            return (log(abs), arg)

    def _eval_is_rational(self):
        s = self.func(*self.args)
        if s.func == self.func:
            if s.args[0].is_rational:
                return False
        else:
            return s.is_rational

    def _eval_is_real(self):
        return self.args[0].is_positive

    def _eval_is_bounded(self):
        arg = self.args[0]
        if arg.is_infinitesimal:
            return False
        return arg.is_bounded

    def _eval_is_positive(self):
        arg = self.args[0]
        if arg.is_positive:
            if arg.is_unbounded:
                return True
            if arg.is_infinitesimal:
                return False
            return (arg - 1).is_positive

    def _eval_is_zero(self):
        # XXX This is not quite useless. Try evaluating log(0.5).is_negative
        #     without it. There's probably a nicer way though.
        if self.args[0] is S.One:
            return True
        elif self.args[0].is_number:
            return self.args[0].expand() is S.One
        elif self.args[0].is_negative:
            return False

    def _eval_nseries(self, x, n, logx):
        # NOTE Please see the comment at the beginning of this file, labelled
        #      IMPORTANT.
        from sympy import cancel
        if not logx:
            logx = log(x)
        if self.args[0] == x:
            return logx
        arg = self.args[0]
        k, l = Wild("k"), Wild("l")
        r = arg.match(k*x**l)
        if r is not None:
            #k = r.get(r, S.One)
            #l = r.get(l, S.Zero)
            k, l = r[k], r[l]
            if l != 0 and not l.has(x) and not k.has(x):
                r = log(k) + l*logx  # XXX true regardless of assumptions?
                return r

        # TODO new and probably slow
        s = self.args[0].nseries(x, n=n, logx=logx)
        while s.is_Order:
            n += 1
            s = self.args[0].nseries(x, n=n, logx=logx)
        a, b = s.leadterm(x)
        p = cancel(s/(a*x**b) - 1)
        g = None
        l = []
        for i in xrange(n + 2):
            g = log.taylor_term(i, p, g)
            g = g.nseries(x, n=n, logx=logx)
            l.append(g)
        return log(a) + b*logx + Add(*l) + C.Order(p**n, x)

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)
        if arg is S.One:
            return (self.args[0] - 1).as_leading_term(x)
        return self.func(arg)

    def _sage_(self):
        import sage.all as sage
        return sage.log(self.args[0]._sage_())


class LambertW(Function):
    """Lambert W function, defined as the inverse function of
    x*exp(x). This function represents the principal branch
    of this inverse function, which like the natural logarithm
    is multivalued.

    For more information, see:
    http://en.wikipedia.org/wiki/Lambert_W_function
    """

    @classmethod
    def eval(cls, x):
        if x == S.Zero:
            return S.Zero
        if x == S.Exp1:
            return S.One
        if x == -1/S.Exp1:
            return S.NegativeOne
        if x == -log(2)/2:
            return -log(2)
        if x == S.Infinity:
            return S.Infinity

    def fdiff(self, argindex=1):
        """
        Return the first derivative of this function.
        """
        if argindex == 1:
            x = self.args[0]
            return LambertW(x)/(x*(1 + LambertW(x)))
        else:
            raise ArgumentIndexError(self, argindex)

from sympy.core.function import _coeff_isneg
