""" This module contains various functions that are special cases
    of incomplete gamma functions. It should probably be renamed. """

from sympy.core import Add, S, C, sympify, cacheit, pi, I
from sympy.core.function import Function, ArgumentIndexError
from sympy.functions.elementary.miscellaneous import sqrt, root
from sympy.functions.elementary.complexes import polar_lift
from sympy.functions.special.hyper import hyper, meijerg

# TODO series expansions
# TODO see the "Note:" in Ei

###############################################################################
################################ ERROR FUNCTION ###############################
###############################################################################

class erf(Function):
    """
    The Gauss error function.

    This function is defined as:

    :math:`\\mathrm{erf}(x)=\\frac{2}{\\sqrt{\\pi}} \\int_0^x e^{-t^2} \\, \\mathrm{d}x`

    Or, in ASCII::

                x
            /
           |
           |     2
           |   -t
        2* |  e    dt
           |
          /
          0
        -------------
              ____
            \/ pi

    Examples
    ========

    >>> from sympy import I, oo, erf
    >>> from sympy.abc import z

    Several special values are known:

    >>> erf(0)
    0
    >>> erf(oo)
    1
    >>> erf(-oo)
    -1
    >>> erf(I*oo)
    oo*I
    >>> erf(-I*oo)
    -oo*I

    In general one can pull out factors of -1 and I from the argument:

    >>> erf(-z)
    -erf(z)

    The error function obeys the mirror symmetry:

    >>> from sympy import conjugate
    >>> conjugate(erf(z))
    erf(conjugate(z))

    Differentiation with respect to z is supported:

    >>> from sympy import diff
    >>> diff(erf(z), z)
    2*exp(-z**2)/sqrt(pi)

    We can numerically evaluate the error function to arbitrary precision
    on the whole complex plane:

    >>> erf(4).evalf(30)
    0.999999984582742099719981147840

    >>> erf(-4*I).evalf(30)
    -1296959.73071763923152794095062*I

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Error_function
    .. [2] http://dlmf.nist.gov/7
    .. [3] http://mathworld.wolfram.com/Erf.html
    .. [4] http://functions.wolfram.com/GammaBetaErf/Erf
    """

    nargs = 1
    unbranched = True

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 2*C.exp(-self.args[0]**2)/sqrt(S.Pi)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.One
            elif arg is S.NegativeInfinity:
                return S.NegativeOne
            elif arg is S.Zero:
                return S.Zero

        t = arg.extract_multiplicatively(S.ImaginaryUnit)
        if t == S.Infinity or t == S.NegativeInfinity:
            return arg

        if arg.could_extract_minus_sign():
            return -cls(-arg)

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)
            k = C.floor((n - 1)/S(2))
            if len(previous_terms) > 2:
                return -previous_terms[-2] * x**2 * (n-2)/(n*k)
            else:
                return 2*(-1)**k * x**n/(n*C.factorial(k)*sqrt(S.Pi))

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def _eval_is_real(self):
        return self.args[0].is_real

    def _eval_rewrite_as_uppergamma(self, z):
        return sqrt(z**2)/z*(S.One - C.uppergamma(S.Half, z**2)/sqrt(S.Pi))

    def _eval_rewrite_as_tractable(self, z):
        return S.One - _erfs(z)*C.exp(-z**2)

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and C.Order(1, x).contains(arg):
            return 2*x/sqrt(pi)
        else:
            return self.func(arg)


###############################################################################
#################### EXPONENTIAL INTEGRALS ####################################
###############################################################################

class Ei(Function):
    r"""
    The classical exponential integral.

    For the use in SymPy, this function is defined as

    .. math:: \operatorname{Ei}(x) = \sum_{n=1}^\infty \frac{x^n}{n\, n!}
                                     + \log(x) + \gamma,

    where :math:`\gamma` is the Euler-Mascheroni constant.

    If :math:`x` is a polar number, this defines an analytic function on the
    riemann surface of the logarithm. Otherwise this defines an analytic
    function in the cut plane :math:`\mathbb{C} \setminus (-\infty, 0]`.

    **Background**

    The name 'exponential integral' comes from the following statement:

    .. math:: \operatorname{Ei}(x) = \int_{-\infty}^x \frac{e^t}{t} \mathrm{d}t

    If the integral is interpreted as a Cauchy principal value, this statement
    holds for :math:`x > 0` and :math:`\operatorname{Ei}(x)` as defined above.

    Note that we carefully avoided defining :math:`\operatorname{Ei}(x)` for
    negative real x. This is because above integral formula does not hold for
    any polar lift of such :math:`x`, indeed all branches of
    :math:`\operatorname{Ei}(x)` above the negative reals are imaginary.

    However, the following statement holds for all :math:`x \in \mathbb{R}^*`:

    .. math:: \int_{-\infty}^x \frac{e^t}{t} \mathrm{d}t =
              \frac{\operatorname{Ei}\left(|x|e^{i \arg(x)}\right) +
                    \operatorname{Ei}\left(|x|e^{- i \arg(x)}\right)}{2},

    where the integral is again understood to be a principal value if
    :math:`x > 0`, and :math:`|x|e^{i \arg(x)}`,
    :math:`|x|e^{- i \arg(x)}` denote two conjugate polar lifts of :math:`x`.

    See Also
    ========

    expint, sympy.functions.special.gamma_functions.uppergamma

    References
    ==========

    - Abramowitz & Stegun, section 5: http://www.math.sfu.ca/~cbm/aands/page_228.htm
    - http://en.wikipedia.org/wiki/Exponential_integral

    Examples
    ========

    >>> from sympy import Ei, polar_lift, exp_polar, I, pi
    >>> from sympy.abc import x

    The exponential integral in SymPy is strictly undefined for negative values
    of the argument. For convenience, exponential integrals with negative
    arguments are immediately converted into an expression that agrees with
    the classical integral definition:

    >>> Ei(-1)
    -I*pi + Ei(exp_polar(I*pi))

    This yields a real value:

    >>> Ei(-1).n(chop=True)
    -0.219383934395520

    On the other hand the analytic continuation is not real:

    >>> Ei(polar_lift(-1)).n(chop=True)
    -0.21938393439552 + 3.14159265358979*I

    The exponential integral has a logarithmic branch point at the origin:

    >>> Ei(x*exp_polar(2*I*pi))
    Ei(x) + 2*I*pi

    Differentiation is supported:

    >>> Ei(x).diff(x)
    exp(x)/x

    The exponential integral is related to many other special functions.
    For example:

    >>> from sympy import uppergamma, expint, Shi
    >>> Ei(x).rewrite(expint)
    -expint(1, x*exp_polar(I*pi)) - I*pi
    >>> Ei(x).rewrite(Shi)
    Chi(x) + Shi(x)

    """

    nargs = 1

    @classmethod
    def eval(cls, z):
        if not z.is_polar and z.is_negative:
            # Note: is this a good idea?
            return Ei(polar_lift(z)) - pi*I
        nz, n = z.extract_branch_factor()
        if n:
            return Ei(nz) + 2*I*pi*n

    def fdiff(self, argindex=1):
        from sympy import unpolarify
        arg = unpolarify(self.args[0])
        if argindex == 1:
            return C.exp(arg)/arg
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_evalf(self, prec):
        if (self.args[0]/polar_lift(-1)).is_positive:
            return Function._eval_evalf(self, prec) + (I*pi)._eval_evalf(prec)
        return Function._eval_evalf(self, prec)

    def _eval_rewrite_as_uppergamma(self, z):
        from sympy import uppergamma
        # XXX this does not currently work usefully because uppergamma
        #     immediately turns into expint
        return -uppergamma(0, polar_lift(-1)*z) - I*pi

    def _eval_rewrite_as_expint(self, z):
        return -expint(1, polar_lift(-1)*z) - I*pi

    def _eval_rewrite_as_Si(self, z):
        return Shi(z) + Chi(z)
    _eval_rewrite_as_Ci = _eval_rewrite_as_Si
    _eval_rewrite_as_Chi = _eval_rewrite_as_Si
    _eval_rewrite_as_Shi = _eval_rewrite_as_Si

class expint(Function):
    r"""
    Generalized exponential integral.

    This function is defined as

    .. math:: \operatorname{E}_\nu(z) = z^{\nu - 1} \Gamma(1 - \nu, z),

    where `\Gamma(1 - \nu, z)` is the upper incomplete gamma function
    (``uppergamma``).

    Hence for :math:`z` with positive real part we have

    .. math:: \operatorname{E}_\nu(z)
              =   \int_1^\infty \frac{e^{-zt}}{z^\nu} \mathrm{d}t,

    which explains the name.

    The representation as an incomplete gamma function provides an analytic
    continuation for :math:`\operatorname{E}_\nu(z)`. If :math:`\nu` is a
    non-positive integer the exponential integral is thus an unbranched
    function of :math:`z`, otherwise there is a branch point at the origin.
    Refer to the incomplete gamma function documentation for details of the
    branching behavior.

    See Also
    ========

    E1: The classical case, returns expint(1, z).
    Ei: Another related function called exponential integral.
    sympy.functions.special.gamma_functions.uppergamma

    References
    ==========

    - http://dlmf.nist.gov/8.19
    - http://functions.wolfram.com/GammaBetaErf/ExpIntegralE/
    - http://en.wikipedia.org/wiki/Exponential_integral

    Examples
    ========

    >>> from sympy import expint, S
    >>> from sympy.abc import nu, z

    Differentiation is supported. Differentiation with respect to z explains
    further the name: for integral orders, the exponential integral is an
    iterated integral of the exponential function.

    >>> expint(nu, z).diff(z)
    -expint(nu - 1, z)

    Differentiation with respect to nu has no classical expression:

    >>> expint(nu, z).diff(nu)
    -z**(nu - 1)*meijerg(((), (1, 1)), ((0, 0, -nu + 1), ()), z)

    At non-postive integer orders, the exponential integral reduces to the
    exponential function:

    >>> expint(0, z)
    exp(-z)/z
    >>> expint(-1, z)
    exp(-z)/z + exp(-z)/z**2

    At half-integers it reduces to error functions:

    >>> expint(S(1)/2, z)
    -sqrt(pi)*erf(sqrt(z))/sqrt(z) + sqrt(pi)/sqrt(z)

    At positive integer orders it can be rewritten in terms of exponentials
    and expint(1, z). Use expand_func() to do this:

    >>> from sympy import expand_func
    >>> expand_func(expint(5, z))
    z**4*expint(1, z)/24 + (-z**3 + z**2 - 2*z + 6)*exp(-z)/24

    The generalised exponential integral is essentially equivalent to the
    incomplete gamma function:

    >>> from sympy import uppergamma
    >>> expint(nu, z).rewrite(uppergamma)
    z**(nu - 1)*uppergamma(-nu + 1, z)

    As such it is branched at the origin:

    >>> from sympy import exp_polar, pi, I
    >>> expint(4, z*exp_polar(2*pi*I))
    I*pi*z**3/3 + expint(4, z)
    >>> expint(nu, z*exp_polar(2*pi*I))
    z**(nu - 1)*(exp(2*I*pi*nu) - 1)*gamma(-nu + 1) + expint(nu, z)

    """

    nargs = 2

    @classmethod
    def eval(cls, nu, z):
        from sympy import (unpolarify, expand_mul, uppergamma, exp, gamma,
                           factorial)
        nu2 = unpolarify(nu)
        if nu != nu2:
            return expint(nu2, z)
        if nu.is_Integer and nu <= 0 or (not nu.is_Integer and (2*nu).is_Integer):
            return unpolarify(expand_mul(z**(nu - 1)*uppergamma(1 - nu, z)))

        # Extract branching information. This can be deduced from what is
        # explained in lowergamma.eval().
        z, n = z.extract_branch_factor()
        if n == 0:
            return
        if nu.is_integer:
            if (nu > 0) is not True:
                return
            return expint(nu, z) \
               - 2*pi*I*n*(-1)**(nu-1)/factorial(nu-1)*unpolarify(z)**(nu - 1)
        else:
            return (exp(2*I*pi*nu*n) - 1)*z**(nu-1)*gamma(1 - nu) + expint(nu, z)

    def fdiff(self, argindex):
        from sympy import meijerg
        nu, z = self.args
        if argindex == 1:
            return -z**(nu - 1)*meijerg([], [1, 1], [0, 0, 1 - nu], [], z)
        elif argindex == 2:
            return -expint(nu - 1, z)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_rewrite_as_uppergamma(self, nu, z):
        from sympy import uppergamma
        return z**(nu - 1)*uppergamma(1 - nu, z)

    def _eval_rewrite_as_Ei(self, nu, z):
        from sympy import exp_polar, unpolarify, exp, factorial
        if nu == 1:
            return -Ei(z*exp_polar(-I*pi)) - I*pi
        elif nu.is_Integer and nu > 1:
            # DLMF, 8.19.7
            x = -unpolarify(z)
            return x**(nu-1)/factorial(nu - 1)*E1(z).rewrite(Ei) + \
                   exp(x)/factorial(nu - 1) * \
                   Add(*[factorial(nu - k - 2)*x**k for k in range(nu - 1)])
        else:
            return self

    def _eval_expand_func(self, **hints):
        return self.rewrite(Ei).rewrite(expint, **hints)

    def _eval_rewrite_as_Si(self, nu, z):
        if nu != 1:
            return self
        return Shi(z) - Chi(z)
    _eval_rewrite_as_Ci = _eval_rewrite_as_Si
    _eval_rewrite_as_Chi = _eval_rewrite_as_Si
    _eval_rewrite_as_Shi = _eval_rewrite_as_Si

def E1(z):
    """
    Classical case of the generalized exponential integral.

    This is equivalent to ``expint(1, z)``.

    """
    return expint(1, z)


###############################################################################
#################### TRIGONOMETRIC INTEGRALS ##################################
###############################################################################

class TrigonometricIntegral(Function):
    """ Base class for trigonometric integrals. """

    nargs = 1

    @classmethod
    def eval(cls, z):
        if z == 0:
            return cls._atzero
        nz = z.extract_multiplicatively(polar_lift(I))
        if nz is None and cls._trigfunc(0) == 0:
            nz = z.extract_multiplicatively(I)
        if nz is not None:
            return cls._Ifactor(nz, 1)
        nz = z.extract_multiplicatively(polar_lift(-I))
        if nz is not None:
            return cls._Ifactor(nz, -1)

        nz = z.extract_multiplicatively(polar_lift(-1))
        if nz is None and cls._trigfunc(0) == 0:
            nz = z.extract_multiplicatively(-1)
        if nz is not None:
            return cls._minusfactor(nz)

        nz, n = z.extract_branch_factor()
        if n == 0 and nz == z:
            return
        return 2*pi*I*n*cls._trigfunc(0) + cls(nz)

    def fdiff(self, argindex=1):
        from sympy import unpolarify
        arg = unpolarify(self.args[0])
        if argindex == 1:
            return self._trigfunc(arg)/arg

    def _eval_rewrite_as_Ei(self, z):
        return self._eval_rewrite_as_expint(z).rewrite(Ei)

    def _eval_rewrite_as_uppergamma(self, z):
        from sympy import uppergamma
        return self._eval_rewrite_as_expint(z).rewrite(uppergamma)

    def _eval_nseries(self, x, n, logx):
        # NOTE this is fairly inefficient
        from sympy import log, EulerGamma, Pow
        n += 1
        if self.args[0].subs(x, 0) != 0:
            return super(TrigonometricIntegral, self)._eval_nseries(x, n, logx)
        baseseries = self._trigfunc(x)._eval_nseries(x, n, logx)
        if self._trigfunc(0) != 0:
            baseseries -= 1
        baseseries = baseseries.replace(Pow, lambda t, n: t**n/n)
        if self._trigfunc(0) != 0:
            baseseries += EulerGamma + log(x)
        return baseseries.subs(x, self.args[0])._eval_nseries(x, n, logx)

class Si(TrigonometricIntegral):
    r"""
    Sine integral.

    This function is defined by

    .. math:: \operatorname{Si}(z) = \int_0^z \frac{\sin{t}}{t} \mathrm{d}t.

    It is an entire function.

    See Also
    ========

    Ci: Cosine integral.
    Shi: Sinh integral.
    Chi: Cosh integral.
    expint: The generalised exponential integral.

    References
    ==========

    - http://en.wikipedia.org/wiki/Trigonometric_integral

    Examples
    ========

    >>> from sympy import Si
    >>> from sympy.abc import z

    The sine integral is an antiderivative of sin(z)/z:

    >>> Si(z).diff(z)
    sin(z)/z

    It is unbranched:

    >>> from sympy import exp_polar, I, pi
    >>> Si(z*exp_polar(2*I*pi))
    Si(z)

    Sine integral behaves much like ordinary sine under multiplication by I:

    >>> Si(I*z)
    I*Shi(z)
    >>> Si(-z)
    -Si(z)

    It can also be expressed in terms of exponential integrals, but beware
    that the latter is branched:

    >>> from sympy import expint
    >>> Si(z).rewrite(expint)
    -I*(-expint(1, z*exp_polar(-I*pi/2))/2 + expint(1, z*exp_polar(I*pi/2))/2) + pi/2

    """

    _trigfunc = C.sin
    _atzero = S(0)

    @classmethod
    def _minusfactor(cls, z):
        return -Si(z)

    @classmethod
    def _Ifactor(cls, z, sign):
        return I*Shi(z)*sign

    def _eval_rewrite_as_expint(self, z):
        # XXX should we polarify z?
        return pi/2 + (E1(polar_lift(I)*z) - E1(polar_lift(-I)*z))/2/I

class Ci(TrigonometricIntegral):
    r"""
    Cosine integral.

    This function is defined for positive :math:`x` by

    .. math:: \operatorname{Ci}(x) = \gamma + \log{x}
                         + \int_0^x \frac{\cos{t} - 1}{t} \mathrm{d}t
           = -\int_x^\infty \frac{\cos{t}}{t} \mathrm{d}t,

    where :math:`\gamma` is the Euler-Mascheroni constant.

    We have

    .. math:: \operatorname{Ci}(z) =
        -\frac{\operatorname{E}_1\left(e^{i\pi/2} z\right)
               + \operatorname{E}_1\left(e^{-i \pi/2} z\right)}{2}

    which holds for all polar :math:`z` and thus provides an analytic
    continuation to the Riemann surface of the logarithm.

    The formula also holds as stated
    for :math:`z \in \mathbb{C}` with :math:`Re(z) > 0`.
    By lifting to the principal branch we obtain an analytic function on the
    cut complex plane.

    See Also
    ========

    Si: Sine integral.
    Shi: Sinh integral.
    Chi: Cosh integral.
    expint: The generalised exponential integral.

    References
    ==========

    - http://en.wikipedia.org/wiki/Trigonometric_integral

    Examples
    ========

    >>> from sympy import Ci
    >>> from sympy.abc import z

    The cosine integral is a primitive of cos(z)/z:

    >>> Ci(z).diff(z)
    cos(z)/z

    It has a logarithmic branch point at the origin:

    >>> from sympy import exp_polar, I, pi
    >>> Ci(z*exp_polar(2*I*pi))
    Ci(z) + 2*I*pi

    Cosine integral behaves somewhat like ordinary cos under multiplication by I:

    >>> from sympy import polar_lift
    >>> Ci(polar_lift(I)*z)
    Chi(z) + I*pi/2
    >>> Ci(polar_lift(-1)*z)
    Ci(z) + I*pi

    It can also be expressed in terms of exponential integrals:

    >>> from sympy import expint
    >>> Ci(z).rewrite(expint)
    -expint(1, z*exp_polar(-I*pi/2))/2 - expint(1, z*exp_polar(I*pi/2))/2

    """

    _trigfunc = C.cos
    _atzero = S.ComplexInfinity

    @classmethod
    def _minusfactor(cls, z):
        return Ci(z) + I*pi

    @classmethod
    def _Ifactor(cls, z, sign):
        return Chi(z) + I*pi/2*sign

    def _eval_rewrite_as_expint(self, z):
        return -(E1(polar_lift(I)*z) + E1(polar_lift(-I)*z))/2

class Shi(TrigonometricIntegral):
    r"""
    Sinh integral.

    This function is defined by

    .. math:: \operatorname{Shi}(z) = \int_0^z \frac{\sinh{t}}{t} \mathrm{d}t.

    It is an entire function.

    See Also
    ========

    Si: Sine integral.
    Ci: Cosine integral.
    Chi: Cosh integral.
    expint: The generalised exponential integral.

    References
    ==========

    - http://en.wikipedia.org/wiki/Trigonometric_integral

    Examples
    ========

    >>> from sympy import Shi
    >>> from sympy.abc import z

    The Sinh integral is a primitive of sinh(z)/z:

    >>> Shi(z).diff(z)
    sinh(z)/z

    It is unbranched:

    >>> from sympy import exp_polar, I, pi
    >>> Shi(z*exp_polar(2*I*pi))
    Shi(z)

    Sinh integral behaves much like ordinary sinh under multiplication by I:

    >>> Shi(I*z)
    I*Si(z)
    >>> Shi(-z)
    -Shi(z)

    It can also be expressed in terms of exponential integrals, but beware
    that the latter is branched:

    >>> from sympy import expint
    >>> Shi(z).rewrite(expint)
    expint(1, z)/2 - expint(1, z*exp_polar(I*pi))/2 - I*pi/2

    """

    _trigfunc = C.sinh
    _atzero = S(0)

    @classmethod
    def _minusfactor(cls, z):
        return -Shi(z)

    @classmethod
    def _Ifactor(cls, z, sign):
        return I*Si(z)*sign

    def _eval_rewrite_as_expint(self, z):
        from sympy import exp_polar
        # XXX should we polarify z?
        return (E1(z)-E1(exp_polar(I*pi)*z))/2 - I*pi/2

class Chi(TrigonometricIntegral):
    r"""
    Cosh integral.

    This function is defined for positive :math:`x` by

    .. math:: \operatorname{Chi}(x) = \gamma + \log{x}
                         + \int_0^x \frac{\cosh{t} - 1}{t} \mathrm{d}t,

    where :math:`\gamma` is the Euler-Mascheroni constant.

    We have

    .. math:: \operatorname{Chi}(z) = \operatorname{Ci}\left(e^{i \pi/2}z\right)
                         - i\frac{\pi}{2},

    which holds for all polar :math:`z` and thus provides an analytic
    continuation to the Riemann surface of the logarithm.
    By lifting to the principal branch we obtain an analytic function on the
    cut complex plane.

    See Also
    ========

    Si: Sine integral.
    Ci: Cosine integral.
    Shi: Sinh integral.
    expint: The generalised exponential integral.

    References
    ==========

    - http://en.wikipedia.org/wiki/Trigonometric_integral

    Examples
    ========

    >>> from sympy import Chi
    >>> from sympy.abc import z

    The cosh integral is a primitive of cosh(z)/z:

    >>> Chi(z).diff(z)
    cosh(z)/z

    It has a logarithmic branch point at the origin:

    >>> from sympy import exp_polar, I, pi
    >>> Chi(z*exp_polar(2*I*pi))
    Chi(z) + 2*I*pi

    Cosh integral behaves somewhat like ordinary cosh under multiplication by I:

    >>> from sympy import polar_lift
    >>> Chi(polar_lift(I)*z)
    Ci(z) + I*pi/2
    >>> Chi(polar_lift(-1)*z)
    Chi(z) + I*pi

    It can also be expressed in terms of exponential integrals:

    >>> from sympy import expint
    >>> Chi(z).rewrite(expint)
    -expint(1, z)/2 - expint(1, z*exp_polar(I*pi))/2 - I*pi/2

    """

    _trigfunc = C.cosh
    _atzero = S.ComplexInfinity

    @classmethod
    def _minusfactor(cls, z):
        return Chi(z) + I*pi

    @classmethod
    def _Ifactor(cls, z, sign):
        return Ci(z) + I*pi/2*sign

    def _eval_rewrite_as_expint(self, z):
        from sympy import exp_polar
        return -I*pi/2 - (E1(z) + E1(exp_polar(I*pi)*z))/2


###############################################################################
#################### FRESNEL INTEGRALS ########################################
###############################################################################

class FresnelIntegral(Function):
    """ Base class for the Fresnel integrals."""

    nargs = 1
    unbranched = True

    @classmethod
    def eval(cls, z):
        # Value at zero
        if z is S.Zero:
            return S(0)

        # Try to pull out factors of -1 and I
        prefact = S.One
        newarg = z
        changed = False

        nz = newarg.extract_multiplicatively(-1)
        if nz is not None:
            prefact = -prefact
            newarg = nz
            changed = True

        nz = newarg.extract_multiplicatively(I)
        if nz is not None:
            prefact = cls._sign*I*prefact
            newarg = nz
            changed = True

        if changed:
            return prefact*cls(newarg)

        # Values at positive infinities signs
        # if any were extracted automatically
        if z is S.Infinity:
            return S.Half
        elif z is I*S.Infinity:
            return cls._sign*I*S.Half

    def fdiff(self, argindex=1):
        if argindex == 1:
            return self._trigfunc(S.Half*pi*self.args[0]**2)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_is_real(self):
        return self.args[0].is_real

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def _as_real_imag(self, deep=True, **hints):
        if self.args[0].is_real:
            if deep:
                hints['complex'] = False
                return (self.expand(deep, **hints), S.Zero)
            else:
                return (self, S.Zero)
        if deep:
            re, im = self.args[0].expand(deep, **hints).as_real_imag()
        else:
            re, im = self.args[0].as_real_imag()
        return (re, im)

    def as_real_imag(self, deep=True, **hints):
        # Fresnel S
        # http://functions.wolfram.com/06.32.19.0003.01
        # http://functions.wolfram.com/06.32.19.0006.01
        # Fresnel C
        # http://functions.wolfram.com/06.33.19.0003.01
        # http://functions.wolfram.com/06.33.19.0006.01
        x, y = self._as_real_imag(deep=deep, **hints)
        sq = -y**2/x**2
        re = S.Half*(self.func(x+x*sqrt(sq))+self.func(x-x*sqrt(sq)))
        im = x/(2*y) * sqrt(sq) * (self.func(x-x*sqrt(sq)) - self.func(x+x*sqrt(sq)))
        return (re, im)


class fresnels(FresnelIntegral):
    r"""
    Fresnel integral S.

    This function is defined by

    .. math:: \operatorname{S}(z) = \int_0^z \sin{\frac{\pi}{2} t^2} \mathrm{d}t.

    It is an entire function.

    Examples
    ========

    >>> from sympy import I, oo, fresnels
    >>> from sympy.abc import z

    Several special values are known:

    >>> fresnels(0)
    0
    >>> fresnels(oo)
    1/2
    >>> fresnels(-oo)
    -1/2
    >>> fresnels(I*oo)
    -I/2
    >>> fresnels(-I*oo)
    I/2

    In general one can pull out factors of -1 and I from the argument:
    >>> fresnels(-z)
    -fresnels(z)

    >>> fresnels(I*z)
    -I*fresnels(z)

    The Fresnel S integral obeys the mirror symmetry:

    >>> from sympy import conjugate
    >>> conjugate(fresnels(z))
    fresnels(conjugate(z))

    Differentiation with respect to z is supported:

    >>> from sympy import diff
    >>> diff(fresnels(z), z)
    sin(pi*z**2/2)

    Defining the Fresnel functions via an integral

    >>> from sympy import integrate, pi, sin, gamma, expand_func
    >>> integrate(sin(pi*z**2/2), z)
    3*fresnels(z)*gamma(3/4)/(4*gamma(7/4))
    >>> expand_func(integrate(sin(pi*z**2/2), z))
    fresnels(z)

    We can numerically evaluate the Fresnel integral to arbitrary precision
    on the whole complex plane:

    >>> fresnels(2).evalf(30)
    0.343415678363698242195300815958

    >>> fresnels(-2*I).evalf(30)
    0.343415678363698242195300815958*I

    See Also
    ========

    fresnelc

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Fresnel_integral
    .. [2] http://dlmf.nist.gov/7
    .. [3] http://mathworld.wolfram.com/FresnelIntegrals.html
    .. [4] http://functions.wolfram.com/GammaBetaErf/FresnelS
    """

    _trigfunc = C.sin
    _sign = -S.One

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0:
            return S.Zero
        else:
            x = sympify(x)
            if len(previous_terms) > 1:
                p = previous_terms[-1]
                return (-pi**2*x**4*(4*n - 1)/(8*n*(2*n + 1)*(4*n + 3))) * p
            else:
                return x**3 * (-x**4)**n * (S(2)**(-2*n - 1)*pi**(2*n + 1)) / ((4*n + 3)*C.factorial(2*n + 1))

    def _eval_rewrite_as_erf(self, z):
        return (S.One + I)/4 * (erf((S.One + I)/2*sqrt(pi)*z) - I*erf((S.One - I)/2*sqrt(pi)*z))

    def _eval_rewrite_as_hyper(self, z):
        return pi*z**3/6 * hyper([S(3)/4], [S(3)/2, S(7)/4], -pi**2*z**4/16)

    def _eval_rewrite_as_meijerg(self, z):
        return (pi*z**(S(9)/4) / (sqrt(2)*(z**2)**(S(3)/4)*(-z)**(S(3)/4))
                * meijerg([],[1],[S(3)/4],[S(1)/4,0],-pi**2*z**4/16))

class fresnelc(FresnelIntegral):
    r"""
    Fresnel integral C.

    This function is defined by

    .. math:: \operatorname{C}(z) = \int_0^z \cos{\frac{\pi}{2} t^2} \mathrm{d}t.

    It is an entire function.

    Examples
    ========

    >>> from sympy import I, oo, fresnelc
    >>> from sympy.abc import z

    Several special values are known:

    >>> fresnelc(0)
    0
    >>> fresnelc(oo)
    1/2
    >>> fresnelc(-oo)
    -1/2
    >>> fresnelc(I*oo)
    I/2
    >>> fresnelc(-I*oo)
    -I/2

    In general one can pull out factors of -1 and I from the argument:
    >>> fresnelc(-z)
    -fresnelc(z)

    >>> fresnelc(I*z)
    I*fresnelc(z)

    The Fresnel C integral obeys the mirror symmetry:

    >>> from sympy import conjugate
    >>> conjugate(fresnelc(z))
    fresnelc(conjugate(z))

    Differentiation with respect to z is supported:

    >>> from sympy import diff
    >>> diff(fresnelc(z), z)
    cos(pi*z**2/2)

    Defining the Fresnel functions via an integral

    >>> from sympy import integrate, pi, cos, gamma, expand_func
    >>> integrate(cos(pi*z**2/2), z)
    fresnelc(z)*gamma(1/4)/(4*gamma(5/4))
    >>> expand_func(integrate(cos(pi*z**2/2), z))
    fresnelc(z)

    We can numerically evaluate the Fresnel integral to arbitrary precision
    on the whole complex plane:

    >>> fresnelc(2).evalf(30)
    0.488253406075340754500223503357

    >>> fresnelc(-2*I).evalf(30)
    -0.488253406075340754500223503357*I

    See Also
    ========

    fresnels

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Fresnel_integral
    .. [2] http://dlmf.nist.gov/7
    .. [3] http://mathworld.wolfram.com/FresnelIntegrals.html
    .. [4] http://functions.wolfram.com/GammaBetaErf/FresnelC
    """

    _trigfunc = C.cos
    _sign = S.One

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0:
            return S.Zero
        else:
            x = sympify(x)
            if len(previous_terms) > 1:
                p = previous_terms[-1]
                return (-pi**2*x**4*(4*n - 3)/(8*n*(2*n - 1)*(4*n + 1))) * p
            else:
                return x * (-x**4)**n * (S(2)**(-2*n)*pi**(2*n)) / ((4*n + 1)*C.factorial(2*n))

    def _eval_rewrite_as_erf(self, z):
        return (S.One - I)/4 * (erf((S.One + I)/2*sqrt(pi)*z) + I*erf((S.One - I)/2*sqrt(pi)*z))

    def _eval_rewrite_as_hyper(self, z):
        return z * hyper([S.One/4], [S.One/2, S(5)/4], -pi**2*z**4/16)

    def _eval_rewrite_as_meijerg(self, z):
        return (pi*z**(S(3)/4) / (sqrt(2)*root(z**2,4)*root(-z,4))
                * meijerg([],[1],[S(1)/4],[S(3)/4,0],-pi**2*z**4/16))

###############################################################################
#################### HELPER FUNCTIONS #########################################
###############################################################################

class _erfs(Function):
    """
    Helper function to make the :math:`erf(z)` function
    tractable for the Gruntz algorithm.
    """

    nargs = 1

    def _eval_aseries(self, n, args0, x, logx):
        if args0[0] != S.Infinity:
            return super(_erfs, self)._eval_aseries(n, args0, x, logx)

        z = self.args[0]
        l = [ 1/sqrt(S.Pi) * C.factorial(2*k)*(-S(4))**(-k)/C.factorial(k) * (1/z)**(2*k+1) for k in xrange(0,n) ]
        o = C.Order(1/z**(2*n+1), x)
        # It is very inefficient to first add the order and then do the nseries
        return (Add(*l))._eval_nseries(x, n, logx) + o

    def fdiff(self, argindex=1):
        if argindex == 1:
            z = self.args[0]
            return -2/sqrt(S.Pi)+2*z*_erfs(z)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_rewrite_as_intractable(self, z):
        return (S.One-erf(z))*C.exp(z**2)
