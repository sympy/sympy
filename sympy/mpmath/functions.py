"""
This module defines most special functions and mathematical constants
provided by mpmath. [Exception: elliptic functions are currently
in elliptic.py]

Most of the actual computational code is located in the lib* modules
(libelefun, libhyper, ...); this module simply wraps this code to
handle precision management in a user friendly way, provide type
conversions, etc.

In addition, this module defines a number of functions that would
be inconvenient to define in the lib* modules, due to requiring
high level operations (e.g. numerical quadrature) for the computation,
or the need to support multiple arguments of mixed types.

"""

import libmpf
import libelefun
import libmpc
import libmpi
import gammazeta
import libhyper

from settings import dps_to_prec

from mptypes import (\
    mpnumeric, mpmathify,
    mpf, make_mpf,
    mpc, make_mpc,
    mpi, make_mpi,
    constant,
    prec_rounding, mp,
    extraprec,
    zero, one, inf, ninf, nan, j, isnan, isinf, isint, eps,
    ComplexResult,
)

class _pi(constant):
    r"""
    `\pi`, roughly equal to 3.141592654, represents the area of the unit
    circle, the half-period of trigonometric functions, and many other
    things in mathematics.

    Mpmath can evaluate `\pi` to arbitrary precision::

        >>> from sympy.mpmath import *
        >>> mp.dps = 50
        >>> print pi
        3.1415926535897932384626433832795028841971693993751

    This shows digits 99991-100000 of `\pi`::

        >>> mp.dps = 100000
        >>> str(pi)[-10:]
        '5549362464'

    **Possible issues**

    :data:`pi` always rounds to the nearest floating-point
    number when used. This means that exact mathematical identities
    involving `\pi` will generally not be preserved in floating-point
    arithmetic. In particular, multiples of :data:`pi` (except for
    the trivial case ``0*pi``) are *not* the exact roots of
    :func:`sin`, but differ roughly by the current epsilon::

        >>> mp.dps = 15
        >>> sin(pi)
        mpf('1.2246467991473532e-16')

    One solution is to use the :func:`sinpi` function instead::

        >>> sinpi(1)
        mpf('0.0')

    See the documentation of trigonometric functions for additional
    details.
    """

class _degree(constant):
    """
    Represents one degree of angle, `1^{\circ} = \pi/180`, or
    about 0.01745329. This constant may be evaluated to arbitrary
    precision::

        >>> from sympy.mpmath import *
        >>> mp.dps = 50
        >>> print degree
        0.017453292519943295769236907684886127134428718885417

    The :data:`degree` object is convenient for conversion
    to radians::

        >>> print sin(30 * degree)
        0.5
        >>> print asin(0.5) / degree
        30.0
    """

class _e(constant):
    """
    The transcendental number `e` = 2.718281828... is the base of the
    natural logarithm (:func:`ln`) and of the exponential function
    (:func:`exp`).

    Mpmath can be evaluate `e` to arbitrary precision::

        >>> mp.dps = 50
        >>> print e
        2.7182818284590452353602874713526624977572470937

    This shows digits 99991-100000 of `e`::

        >>> mp.dps = 100000
        >>> str(e)[-10:]
        '2100427165'

    **Possible issues**

    :data:`e` always rounds to the nearest floating-point number
    when used, and mathematical identities involving `e` may not
    hold in floating-point arithmetic. For example, ``ln(e)``
    might not evaluate exactly to 1.

    In particular, don't use ``e**x`` to compute the exponential
    function. Use ``exp(x)`` instead; this is both faster and more
    accurate.
    """

class _ln2(constant): pass

class _ln10(constant): pass

class _phi(constant):
    r"""
    Represents the golden ratio `\phi = (1+\sqrt 5)/2`,
    approximately equal to 1.6180339887. To high precision,
    its value is::

        >>> from sympy.mpmath import *
        >>> mp.dps = 50
        >>> print phi
        1.6180339887498948482045868343656381177203091798058

    Formulas for the golden ratio include the following::

        >>> print (1+sqrt(5))/2
        1.6180339887498948482045868343656381177203091798058
        >>> print findroot(lambda x: x**2-x-1, 1)
        1.6180339887498948482045868343656381177203091798058
        >>> print limit(lambda n: fib(n+1)/fib(n), inf)
        1.6180339887498948482045868343656381177203091798058
    """

class _euler(constant):
    r"""
    Euler's constant or the Euler-Mascheroni constant `\gamma`
    = 0.57721566... is a number of central importance to
    number theory and special functions. It is defined as the limit

    .. math ::

        \gamma = \lim_{n\to\infty} H_n - \log n

    where `H_n = 1 + \frac{1}{2} + \ldots + \frac{1}{n}` is a harmonic
    number (see :func:`harmonic`).

    Evaluation of `\gamma` is supported at arbitrary precision::

        >>> from sympy.mpmath import *
        >>> mp.dps = 50
        >>> print euler
        0.57721566490153286060651209008240243104215933593992

    We can also compute `\gamma` directly from the definition,
    although this is less efficient::

        >>> print limit(lambda n: harmonic(n)-log(n), inf)
        0.57721566490153286060651209008240243104215933593992

    This shows digits 9991-10000 of `\gamma`::

        >>> mp.dps = 10000
        >>> str(euler)[-10:]
        '4679858165'

    Integrals, series, and representations for `\gamma` in terms of
    special functions include the following (there are many others)::

        >>> mp.dps = 25
        >>> print -quad(lambda x: exp(-x)*log(x), [0,inf])
        0.5772156649015328606065121
        >>> print quad(lambda x,y: (x-1)/(1-x*y)/log(x*y), [0,1], [0,1])
        0.5772156649015328606065121
        >>> print nsum(lambda k: 1/k-log(1+1/k), [1,inf])
        0.5772156649015328606065121
        >>> print nsum(lambda k: (-1)**k*zeta(k)/k, [2,inf])
        0.5772156649015328606065121
        >>> print -diff(gamma, 1)
        0.5772156649015328606065121
        >>> print limit(lambda x: 1/x-gamma(x), 0)
        0.5772156649015328606065121
        >>> print limit(lambda x: zeta(x)-1/(x-1), 1)
        0.5772156649015328606065121
        >>> print (log(2*pi*nprod(lambda n:
        ...     exp(-2+2/n)*(1+2/n)**n, [1,inf]))-3)/2
        0.5772156649015328606065121

    For generalizations of the identities `\gamma = -\Gamma'(1)`
    and `\gamma = \lim_{x\to1} \zeta(x)-1/(x-1)`, see
    :func:`psi` and :func:`stieltjes` respectively.
    """

class _catalan(constant):
    r"""
    Catalan's constant `K` = 0.91596559... is given by the infinite
    series

    .. math ::

        K = \sum_{k=0}^{\infty} \frac{(-1)^k}{(2k+1)^2}.

    Mpmath can evaluate it to arbitrary precision::

        >>> from sympy.mpmath import *
        >>> mp.dps = 50
        >>> print catalan
        0.91596559417721901505460351493238411077414937428167

    One can also compute `K` directly from the definition, although
    this is significantly less efficient::

        >>> print nsum(lambda k: (-1)**k/(2*k+1)**2, [0, inf])
        0.91596559417721901505460351493238411077414937428167

    This shows digits 9991-10000 of `K`::

        >>> mp.dps = 10000
        >>> str(catalan)[-10:]
        '9537871503'

    Catalan's constant has numerous integral representations::

        >>> mp.dps = 50
        >>> print quad(lambda x: -log(x)/(1+x**2), [0, 1])
        0.91596559417721901505460351493238411077414937428167
        >>> print quad(lambda x: atan(x)/x, [0, 1])
        0.91596559417721901505460351493238411077414937428167
        >>> print quad(lambda x: ellipk(x**2)/2, [0, 1])
        0.91596559417721901505460351493238411077414937428167
        >>> print quad(lambda x,y: 1/(1+(x*y)**2), [0, 1], [0, 1])
        0.91596559417721901505460351493238411077414937428167

    As well as series representations::

        >>> print pi*log(sqrt(3)+2)/8 + 3*nsum(lambda n:
        ...  (fac(n)/(2*n+1))**2/fac(2*n), [0, inf])/8
        0.91596559417721901505460351493238411077414937428167
        >>> print 1-nsum(lambda n: n*zeta(2*n+1)/16**n, [1,inf])
        0.91596559417721901505460351493238411077414937428167
    """

class _khinchin(constant):
    r"""
    Khinchin's constant `K` = 2.68542... is a number that
    appears in the theory of continued fractions. Mpmath can evaluate
    it to arbitrary precision::

        >>> from sympy.mpmath import *
        >>> mp.dps = 50
        >>> print khinchin
        2.6854520010653064453097148354817956938203822939945

    An integral representation is::

        >>> I = quad(lambda x: log((1-x**2)/sincpi(x))/x/(1+x), [0, 1])
        >>> print 2*exp(1/log(2)*I)
        2.6854520010653064453097148354817956938203822939945

    The computation of ``khinchin`` is based on an efficient
    implementation of the following series::

        >>> f = lambda n: (zeta(2*n)-1)/n*sum((-1)**(k+1)/mpf(k)
        ...     for k in range(1,2*n))
        >>> print exp(nsum(f, [1,inf])/log(2))
        2.6854520010653064453097148354817956938203822939945
    """

class _glaisher(constant):
    r"""
    Glaisher's constant `A`, also known as the Glaisher-Kinkelin
    constant, is a number approximately equal to 1.282427129 that
    sometimes appears in formulas related to gamma and zeta functions.
    It is also related to the Barnes G-function (see :func:`barnesg`).

    The constant is defined  as `A = \exp(1/12-\zeta'(-1))` where
    `\zeta'(s)` denotes the derivative of the Riemann zeta function
    (see :func:`zeta`).

    Mpmath can evaluate Glaisher's constant to arbitrary precision:

        >>> from sympy.mpmath import *
        >>> mp.dps = 50
        >>> print glaisher
        1.282427129100622636875342568869791727767688927325

    We can verify that the value computed by :data:`glaisher` is
    correct using mpmath's facilities for numerical
    differentiation and arbitrary evaluation of the zeta function:

        >>> print exp(mpf(1)/12 - diff(zeta, -1))
        1.282427129100622636875342568869791727767688927325

    Here is an example of an integral that can be evaluated in
    terms of Glaisher's constant:

        >>> mp.dps = 15
        >>> print quad(lambda x: log(gamma(x)), [1, 1.5])
        -0.0428537406502909
        >>> print -0.5 - 7*log(2)/24 + log(pi)/4 + 3*log(glaisher)/2
        -0.042853740650291

    Mpmath computes Glaisher's constant by applying Euler-Maclaurin
    summation to a slowly convergent series. The implementation is
    reasonably efficient up to about 10,000 digits. See the source
    code for additional details.

    References:
    http://mathworld.wolfram.com/Glaisher-KinkelinConstant.html
    """

class _apery(constant):
    r"""
    Represents Apery's constant, which is the irrational number
    approximately equal to 1.2020569 given by

    .. math ::

        \zeta(3) = \sum_{k=1}^\infty\frac{1}{k^3}.

    The calculation is based on an efficient hypergeometric
    series. To 50 decimal places, the value is given by::

        >>> from sympy.mpmath import *
        >>> mp.dps = 50
        >>> print apery
        1.2020569031595942853997381615114499907649862923405

    Other ways to evaluate Apery's constant using mpmath
    include::

        >>> print zeta(3)
        1.2020569031595942853997381615114499907649862923405
        >>> print -diff(trigamma, 1)/2
        1.2020569031595942853997381615114499907649862923405
        >>> print 8*nsum(lambda k: 1/(2*k+1)**3, [0,inf])/7
        1.2020569031595942853997381615114499907649862923405
        >>> f = lambda k: 2/k**3/(exp(2*pi*k)-1)
        >>> print 7*pi**3/180 - nsum(f, [1,inf])
        1.2020569031595942853997381615114499907649862923405

    This shows digits 9991-10000 of Apery's constant::

        >>> mp.dps = 10000
        >>> str(apery)[-10:]
        '3189504235'
    """

# Mathematical constants
pi = _pi(libelefun.mpf_pi, "pi")
degree = _degree(libelefun.mpf_degree, "degree")
e = _e(libelefun.mpf_e, "e")
ln2 = _ln2(libelefun.mpf_ln2, "ln(2)")
ln10 = _ln10(libelefun.mpf_ln10, "ln(10)")
phi = _phi(libelefun.mpf_phi, "Golden ratio (phi)")
euler = _euler(gammazeta.mpf_euler, "Euler's constant (gamma)")
catalan = _catalan(gammazeta.mpf_catalan, "Catalan's constant")
khinchin = _khinchin(gammazeta.mpf_khinchin, "Khinchin's constant")
glaisher = _glaisher(gammazeta.mpf_glaisher, "Glaisher's constant")
apery = _apery(gammazeta.mpf_apery, "Apery's constant")


def funcwrapper(f):
    def g(*args, **kwargs):
        orig = mp.prec
        try:
            args = [mpmathify(z) for z in args]
            mp.prec = orig + 10
            v = f(*args, **kwargs)
        finally:
            mp.prec = orig
        return +v
    g.__name__ = f.__name__
    g.__doc__ = f.__doc__
    return g

def mpfunc(name, real_f, complex_f, doc, interval_f=None):
    def f(x, **kwargs):
        if not isinstance(x, mpnumeric):
            x = mpmathify(x)
        prec, rounding = prec_rounding
        if kwargs:
            prec = kwargs.get('prec', prec)
            if 'dps' in kwargs:
                prec = dps_to_prec(kwargs['dps'])
            rounding = kwargs.get('rounding', rounding)
        if isinstance(x, mpf):
            try:
                return make_mpf(real_f(x._mpf_, prec, rounding))
            except ComplexResult:
                # Handle propagation to complex
                if mp.trap_complex:
                    raise
                return make_mpc(complex_f((x._mpf_, libmpf.fzero), prec, rounding))
        elif isinstance(x, mpc):
            return make_mpc(complex_f(x._mpc_, prec, rounding))
        elif isinstance(x, mpi):
            if interval_f:
                return make_mpi(interval_f(x._val, prec))
        raise NotImplementedError("%s of a %s" % (name, type(x)))

    f.__name__ = name
    f.__doc__ = "Computes the %s of x" % doc
    return f

def altfunc(f, name, desc):
    def g(x):
        orig = mp.prec
        try:
            mp.prec = orig + 10
            return one/f(x)
        finally:
            mp.prec = orig
    g.__name__ = name
    g.__doc__ = "Computes the %s of x, 1/%s(x)" % (desc, f.__name__)
    return g

def altinvfunc(f, name, desc):
    def g(x):
        orig = mp.prec
        try:
            mp.prec = orig + 10
            return f(one/x)
        finally:
            mp.prec = orig
    g.__name__ = name
    g.__doc__ = "Computes the inverse %s of x, %s(1/x)" % (desc, f.__name__)
    return g

sqrt = mpfunc('sqrt', libelefun.mpf_sqrt, libmpc.mpc_sqrt, "principal square root", libmpi.mpi_sqrt)
sqrt.__doc__ = r"""
``sqrt(x)`` computes the principal square root of `x`, `\sqrt x`.
For positive real numbers, the principal root is simply the
positive square root. For arbitrary complex numbers, the principal
square root is defined to satisfy `\sqrt x = \exp(\log(x)/2)`.
The function thus has a branch cut along the negative half real axis.

For all mpmath numbers ``x``, calling ``sqrt(x)`` is equivalent to
performing ``x**0.5``.

**Examples**

Basic examples and limits::

    >>> from sympy.mpmath import *
    >>> mp.dps = 15
    >>> print sqrt(10)
    3.16227766016838
    >>> print sqrt(100)
    10.0
    >>> print sqrt(-4)
    (0.0 + 2.0j)
    >>> print sqrt(1+1j)
    (1.09868411346781 + 0.455089860562227j)
    >>> print sqrt(inf)
    +inf

Square root evaluation is fast at huge precision::

    >>> mp.dps = 50000
    >>> a = sqrt(3)
    >>> str(a)[-10:]
    '9329332814'

:func:`sqrt` supports interval arguments::

    >>> mp.dps = 15
    >>> print sqrt(mpi(16, 100))
    [4.0, 10.0]
    >>> print sqrt(mpi(2))
    [1.4142135623730949234, 1.4142135623730951455]
    >>> print sqrt(mpi(2)) ** 2
    [1.9999999999999995559, 2.0000000000000004441]

"""

cbrt = mpfunc('cbrt', libelefun.mpf_cbrt, libmpc.mpc_cbrt, "principal cubic root")
cbrt.__doc__ = """
``cbrt(x)`` computes the cube root of `x`, `x^{1/3}`. This
function is faster and more accurate than raising to a floating-point
fraction::

    >>> from sympy.mpmath import *
    >>> mp.dps = 15
    >>> 125**(mpf(1)/3)
    mpf('4.9999999999999991')
    >>> cbrt(125)
    mpf('5.0')

Every nonzero complex number has three cube roots. This function
returns the cube root defined by `\exp(\log(x)/3)` where the
principal branch of the natural logarithm is used. Note that this
does not give a real cube root for negative real numbers::

    >>> print cbrt(-1)
    (0.5 + 0.866025403784439j)

"""

exp = mpfunc('exp', libelefun.mpf_exp, libmpc.mpc_exp, "exponential function", libmpi.mpi_exp)

exp.__doc__ = r"""
Computes the exponential function,

.. math ::

    \exp(x) = e^x = \sum_{k=0}^{\infty} \frac{x^k}{k!}.

For complex numbers, the exponential function also satisfies

.. math ::

    \exp(x+yi) = e^x (\cos y + i \sin y).

**Basic examples**

Some values of the exponential function::

    >>> from sympy.mpmath import *
    >>> mp.dps = 25
    >>> print exp(0)
    1.0
    >>> print exp(1)
    2.718281828459045235360287
    >>> print exp(-1)
    0.3678794411714423215955238
    >>> print exp(inf)
    +inf
    >>> print exp(-inf)
    0.0

Arguments can be arbitrarily large::

    >>> print exp(10000)
    8.806818225662921587261496e+4342
    >>> print exp(-10000)
    1.135483865314736098540939e-4343

Evaluation is supported for interval arguments::

    >>> print exp(mpi(-inf,0))
    [0.0, 1.0]
    >>> print exp(mpi(0,1))
    [1.0, 2.71828182845904523536028749558]

The exponential function can be evaluated efficiently to arbitrary
precision::

    >>> mp.dps = 10000
    >>> print exp(pi)  #doctest: +ELLIPSIS
    23.140692632779269005729...8984304016040616

**Functional properties**

Numerical verification of Euler's identity for the complex
exponential function::

    >>> mp.dps = 15
    >>> print exp(j*pi)+1
    (0.0 + 1.22464679914735e-16j)
    >>> print chop(exp(j*pi)+1)
    0.0

This recovers the coefficients (reciprocal factorials) in the
Maclaurin series expansion of exp::

    >>> nprint(taylor(exp, 0, 5))
    [1.0, 1.0, 0.5, 0.166667, 4.16667e-2, 8.33333e-3]

The exponential function is its own derivative and antiderivative::

    >>> print exp(pi)
    23.1406926327793
    >>> print diff(exp, pi)
    23.1406926327793
    >>> print quad(exp, [-inf, pi])
    23.1406926327793

The exponential function can be evaluated using various methods,
including direct summation of the series, limits, and solving
the defining differential equation::

    >>> print nsum(lambda k: pi**k/fac(k), [0,inf])
    23.1406926327793
    >>> print limit(lambda k: (1+pi/k)**k, inf)
    23.1406926327793
    >>> print odefun(lambda t, x: x, 0, 1)(pi)
    23.1406926327793


"""


ln = mpfunc('ln', libelefun.mpf_log, libmpc.mpc_log, "natural logarithm", libmpi.mpi_log)

ln.__doc__ = r"""Computes the natural logarithm of `x`, `\ln x`.
See :func:`log` for additional documentation."""

cos = mpfunc('cos', libelefun.mpf_cos, libmpc.mpc_cos, "cosine", libmpi.mpi_cos)
sin = mpfunc('sin', libelefun.mpf_sin, libmpc.mpc_sin, "sine", libmpi.mpi_sin)
tan = mpfunc('tan', libelefun.mpf_tan, libmpc.mpc_tan, "tangent", libmpi.mpi_tan)
cosh = mpfunc('cosh', libelefun.mpf_cosh, libmpc.mpc_cosh, "hyperbolic cosine")
sinh = mpfunc('sinh', libelefun.mpf_sinh, libmpc.mpc_sinh, "hyperbolic sine")
tanh = mpfunc('tanh', libelefun.mpf_tanh, libmpc.mpc_tanh, "hyperbolic tangent")

acos = mpfunc('acos', libelefun.mpf_acos, libmpc.mpc_acos, "inverse cosine")
asin = mpfunc('asin', libelefun.mpf_asin, libmpc.mpc_asin, "inverse sine")
atan = mpfunc('atan', libelefun.mpf_atan, libmpc.mpc_atan, "inverse tangent")
asinh = mpfunc('asinh', libelefun.mpf_asinh, libmpc.mpc_asinh, "inverse hyperbolic sine")
acosh = mpfunc('acosh', libelefun.mpf_acosh, libmpc.mpc_acosh, "inverse hyperbolic cosine")
atanh = mpfunc('atanh', libelefun.mpf_atanh, libmpc.mpc_atanh, "inverse hyperbolic tangent")

sec = altfunc(cos, 'sec', 'secant')
csc = altfunc(sin, 'csc', 'cosecant')
cot = altfunc(tan, 'cot', 'cotangent')
sech = altfunc(cosh, 'sech', 'hyperbolic secant')
csch = altfunc(sinh, 'csch', 'hyperbolic cosecant')
coth = altfunc(tanh, 'coth', 'hyperbolic cotangent')

asec = altinvfunc(acos, 'asec', 'secant')
acsc = altinvfunc(asin, 'acsc', 'cosecant')
acot = altinvfunc(atan, 'acot', 'cotangent')
asech = altinvfunc(acosh, 'asech', 'hyperbolic secant')
acsch = altinvfunc(asinh, 'acsch', 'hyperbolic cosecant')
acoth = altinvfunc(atanh, 'acoth', 'hyperbolic cotangent')

cospi = mpfunc('cospi', libelefun.mpf_cos_pi, libmpc.mpc_cos_pi, "")
sinpi = mpfunc('sinpi', libelefun.mpf_sin_pi, libmpc.mpc_sin_pi, "")

cosh.__doc__ = r"""
Computes the hyperbolic cosine of `x`,
`\cosh(x) = (e^x + e^{-x})/2`. Values and limits include::

    >>> from sympy.mpmath import *
    >>> mp.dps = 25
    >>> print cosh(0)
    1.0
    >>> print cosh(1)
    1.543080634815243778477906
    >>> print cosh(-inf), cosh(+inf)
    +inf +inf

The hyperbolic cosine is an even, convex function with
a global minimum at `x = 0`, having a Maclaurin series
that starts::

    >>> nprint(chop(taylor(cosh, 0, 5)))
    [1.0, 0.0, 0.5, 0.0, 4.16667e-2, 0.0]

Generalized to complex numbers, the hyperbolic cosine is
equivalent to a cosine with the argument rotated
in the imaginary direction, or `\cosh x = \cos ix`::

    >>> print cosh(2+3j)
    (-3.724545504915322565473971 + 0.5118225699873846088344638j)
    >>> print cos(3-2j)
    (-3.724545504915322565473971 + 0.5118225699873846088344638j)
"""

sinh.__doc__ = r"""
Computes the hyperbolic sine of `x`,
`\sinh(x) = (e^x - e^{-x})/2`. Values and limits include::

    >>> from sympy.mpmath import *
    >>> mp.dps = 25
    >>> print sinh(0)
    0.0
    >>> print sinh(1)
    1.175201193643801456882382
    >>> print sinh(-inf), sinh(+inf)
    -inf +inf

The hyperbolic sine is an odd function, with a Maclaurin
series that starts::

    >>> nprint(chop(taylor(sinh, 0, 5)))
    [0.0, 1.0, 0.0, 0.166667, 0.0, 8.33333e-3]

Generalized to complex numbers, the hyperbolic sine is
essentially a sine with a rotation `i` applied to
the argument; more precisely, `\sinh x = -i \sin ix`::

    >>> print sinh(2+3j)
    (-3.590564589985779952012565 + 0.5309210862485198052670401j)
    >>> print j*sin(3-2j)
    (-3.590564589985779952012565 + 0.5309210862485198052670401j)
"""

tanh.__doc__ = r"""
Computes the hyperbolic tangent of `x`,
`\tanh(x) = \sinh(x)/\cosh(x)`. Values and limits include::

    >>> from sympy.mpmath import *
    >>> mp.dps = 25
    >>> print tanh(0)
    0.0
    >>> print tanh(1)
    0.7615941559557648881194583
    >>> print tanh(-inf), tanh(inf)
    -1.0 1.0

The hyperbolic tangent is an odd, sigmoidal function, similar
to the inverse tangent and error function. Its Maclaurin
series is::

    >>> nprint(chop(taylor(tanh, 0, 5)))
    [0.0, 1.0, 0.0, -0.333333, 0.0, 0.133333]

Generalized to complex numbers, the hyperbolic tangent is
essentially a tangent with a rotation `i` applied to
the argument; more precisely, `\tanh x = -i \tan ix`::

    >>> print tanh(2+3j)
    (0.9653858790221331242784803 - 0.009884375038322493720314034j)
    >>> print j*tan(3-2j)
    (0.9653858790221331242784803 - 0.009884375038322493720314034j)
"""

cos.__doc__ = r"""
Computes the cosine of `x`, `\cos(x)`.

    >>> from sympy.mpmath import *
    >>> mp.dps = 25
    >>> print cos(pi/3)
    0.5
    >>> print cos(100000001)
    -0.9802850113244713353133243
    >>> print cos(2+3j)
    (-4.189625690968807230132555 - 9.109227893755336597979197j)
    >>> print cos(inf)
    nan
    >>> nprint(chop(taylor(cos, 0, 6)))
    [1.0, 0.0, -0.5, 0.0, 4.16667e-2, 0.0, -1.38889e-3]
    >>> print cos(mpi(0,1))
    [0.540302305868139717400936602301, 1.0]
    >>> print cos(mpi(0,2))
    [-0.41614683654714238699756823214, 1.0]
"""

sin.__doc__ = r"""
Computes the sine of `x`, `\sin(x)`.

    >>> from sympy.mpmath import *
    >>> mp.dps = 25
    >>> print sin(pi/3)
    0.8660254037844386467637232
    >>> print sin(100000001)
    0.1975887055794968911438743
    >>> print sin(2+3j)
    (9.1544991469114295734673 - 4.168906959966564350754813j)
    >>> print sin(inf)
    nan
    >>> nprint(chop(taylor(sin, 0, 6)))
    [0.0, 1.0, 0.0, -0.166667, 0.0, 8.33333e-3, 0.0]
    >>> print sin(mpi(0,1))
    [0.0, 0.841470984807896506652502331201]
    >>> print sin(mpi(0,2))
    [0.0, 1.0]
"""

tan.__doc__ = r"""
Computes the tangent of `x`, `\tan(x) = \frac{\sin(x)}{\cos(x)}`.
The tangent function is singular at `x = (n+1/2)\pi`, but
``tan(x)`` always returns a finite result since `(n+1/2)\pi`
cannot be represented exactly using floating-point arithmetic.

    >>> from sympy.mpmath import *
    >>> mp.dps = 25
    >>> print tan(pi/3)
    1.732050807568877293527446
    >>> print tan(100000001)
    -0.2015625081449864533091058
    >>> print tan(2+3j)
    (-0.003764025641504248292751221 + 1.003238627353609801446359j)
    >>> print tan(inf)
    nan
    >>> nprint(chop(taylor(tan, 0, 6)))
    [0.0, 1.0, 0.0, 0.333333, 0.0, 0.133333, 0.0]
    >>> print tan(mpi(0,1))
    [0.0, 1.55740772465490223050697482944]
    >>> print tan(mpi(0,2))  # Interval includes a singularity
    [-inf, +inf]
"""

sec.__doc__ = r"""
Computes the secant of `x`, `\mathrm{sec}(x) = \frac{1}{\cos(x)}`.
The secant function is singular at `x = (n+1/2)\pi`, but
``sec(x)`` always returns a finite result since `(n+1/2)\pi`
cannot be represented exactly using floating-point arithmetic.

    >>> from sympy.mpmath import *
    >>> mp.dps = 25
    >>> print sec(pi/3)
    2.0
    >>> print sec(10000001)
    -1.184723164360392819100265
    >>> print sec(2+3j)
    (-0.04167496441114427004834991 + 0.0906111371962375965296612j)
    >>> print sec(inf)
    nan
    >>> nprint(chop(taylor(sec, 0, 6)))
    [1.0, 0.0, 0.5, 0.0, 0.208333, 0.0, 8.47222e-2]
    >>> print sec(mpi(0,1))
    [1.0, 1.85081571768092561791175324143]
    >>> print sec(mpi(0,2))  # Interval includes a singularity
    [-inf, +inf]
"""

csc.__doc__ = r"""
Computes the cosecant of `x`, `\mathrm{csc}(x) = \frac{1}{\sin(x)}`.
This cosecant function is singular at `x = n \pi`, but with the
exception of the point `x = 0`, ``csc(x)`` returns a finite result
since `n \pi` cannot be represented exactly using floating-point
arithmetic.

    >>> from sympy.mpmath import *
    >>> mp.dps = 25
    >>> print csc(pi/3)
    1.154700538379251529018298
    >>> print csc(10000001)
    -1.864910497503629858938891
    >>> print csc(2+3j)
    (0.09047320975320743980579048 + 0.04120098628857412646300981j)
    >>> print csc(inf)
    nan
    >>> print csc(mpi(0,1))  # Interval includes a singularity
    [1.18839510577812121626159945235, +inf]
    >>> print csc(mpi(0,2))
    [1.0, +inf]
"""

cot.__doc__ = r"""
Computes the cotangent of `x`,
`\mathrm{cot}(x) = \frac{1}{\tan(x)} = \frac{\cos(x)}{\sin(x)}`.
This cotangent function is singular at `x = n \pi`, but with the
exception of the point `x = 0`, ``cot(x)`` returns a finite result
since `n \pi` cannot be represented exactly using floating-point
arithmetic.

    >>> from sympy.mpmath import *
    >>> mp.dps = 25
    >>> print cot(pi/3)
    0.5773502691896257645091488
    >>> print cot(10000001)
    1.574131876209625656003562
    >>> print cot(2+3j)
    (-0.003739710376336956660117409 - 0.9967577965693583104609688j)
    >>> print cot(inf)
    nan
    >>> print cot(mpi(0,1))  # Interval includes a singularity
    [0.642092615934330703006419986575, +inf]
    >>> print cot(mpi(1,2))
    [-inf, +inf]
"""

acos.__doc__ = r"""
Computes the inverse cosine or arccosine of `x`, `\cos^{-1}(x)`.
Since `-1 \le \cos(x) \le 1` for real `x`, the inverse
cosine is real-valued only for `-1 \le x \le 1`. On this interval,
:func:`acos` is defined to be a monotonically decreasing
function assuming values between `+\pi` and `0`.

Basic values are::

    >>> from sympy.mpmath import *
    >>> mp.dps = 25
    >>> print acos(-1)
    3.141592653589793238462643
    >>> print acos(0)
    1.570796326794896619231322
    >>> print acos(1)
    0.0
    >>> nprint(chop(taylor(acos, 0, 6)))
    [1.5708, -1.0, 0.0, -0.166667, 0.0, -7.5e-2, 0.0]

:func:`acos` is defined so as to be a proper inverse function of
`\cos(\theta)` for `0 \le \theta < \pi`.
We have `\cos(\cos^{-1}(x)) = x` for all `x`, but
`\cos^{-1}(\cos(x)) = x` only for `0 \le \Re[x] < \pi`::

    >>> for x in [1, 10, -1, 2+3j, 10+3j]:
    ...     print cos(acos(x)), acos(cos(x))
    ...
    1.0 1.0
    (10.0 + 0.0j) 2.566370614359172953850574
    -1.0 1.0
    (2.0 + 3.0j) (2.0 + 3.0j)
    (10.0 + 3.0j) (2.566370614359172953850574 - 3.0j)

The inverse cosine has two branch points: `x = \pm 1`. :func:`acos`
places the branch cuts along the line segments `(-\infty, -1)` and
`(+1, +\infty)`. In general,

.. math ::

    \cos^{-1}(x) = \frac{\pi}{2} + i \log\left(ix + \sqrt{1-x^2} \right)

where the principal-branch log and square root are implied.
"""

asin.__doc__ = r"""
Computes the inverse sine or arcsine of `x`, `\sin^{-1}(x)`.
Since `-1 \le \sin(x) \le 1` for real `x`, the inverse
sine is real-valued only for `-1 \le x \le 1`.
On this interval, it is defined to be a monotonically increasing
function assuming values between `-\pi/2` and `\pi/2`.

Basic values are::

    >>> from sympy.mpmath import *
    >>> mp.dps = 25
    >>> print asin(-1)
    -1.570796326794896619231322
    >>> print asin(0)
    0.0
    >>> print asin(1)
    1.570796326794896619231322
    >>> nprint(chop(taylor(asin, 0, 6)))
    [0.0, 1.0, 0.0, 0.166667, 0.0, 7.5e-2, 0.0]

:func:`asin` is defined so as to be a proper inverse function of
`\sin(\theta)` for `-\pi/2 < \theta < \pi/2`.
We have `\sin(\sin^{-1}(x)) = x` for all `x`, but
`\sin^{-1}(\sin(x)) = x` only for `-\pi/2 < \Re[x] < \pi/2`::

    >>> for x in [1, 10, -1, 1+3j, -2+3j]:
    ...     print chop(sin(asin(x))), asin(sin(x))
    ...
    1.0 1.0
    10.0 -0.5752220392306202846120698
    -1.0 -1.0
    (1.0 + 3.0j) (1.0 + 3.0j)
    (-2.0 + 3.0j) (-1.141592653589793238462643 - 3.0j)

The inverse sine has two branch points: `x = \pm 1`. :func:`asin`
places the branch cuts along the line segments `(-\infty, -1)` and
`(+1, +\infty)`. In general,

.. math ::

    \sin^{-1}(x) = -i \log\left(ix + \sqrt{1-x^2} \right)

where the principal-branch log and square root are implied.
"""

atan.__doc__ = r"""
Computes the inverse tangent or arctangent of `x`, `\tan^{-1}(x)`.
This is a real-valued function for all real `x`, with range
`(-\pi/2, \pi/2)`.

Basic values are::

    >>> from sympy.mpmath import *
    >>> mp.dps = 25
    >>> print atan(-inf)
    -1.570796326794896619231322
    >>> print atan(-1)
    -0.7853981633974483096156609
    >>> print atan(0)
    0.0
    >>> print atan(1)
    0.7853981633974483096156609
    >>> print atan(inf)
    1.570796326794896619231322
    >>> nprint(chop(taylor(atan, 0, 6)))
    [0.0, 1.0, 0.0, -0.333333, 0.0, 0.2, 0.0]

The inverse tangent is often used to compute angles. However,
the atan2 function is often better for this as it preserves sign
(see :func:`atan2`).

:func:`atan` is defined so as to be a proper inverse function of
`\tan(\theta)` for `-\pi/2 < \theta < \pi/2`.
We have `\tan(\tan^{-1}(x)) = x` for all `x`, but
`\tan^{-1}(\tan(x)) = x` only for `-\pi/2 < \Re[x] < \pi/2`::

    >>> mp.dps = 25
    >>> for x in [1, 10, -1, 1+3j, -2+3j]:
    ...     print tan(atan(x)), atan(tan(x))
    ...
    1.0 1.0
    10.0 0.5752220392306202846120698
    -1.0 -1.0
    (1.0 + 3.0j) (1.000000000000000000000001 + 3.0j)
    (-2.0 + 3.0j) (1.141592653589793238462644 + 3.0j)

The inverse tangent has two branch points: `x = \pm i`. :func:`atan`
places the branch cuts along the line segments `(-i \infty, -i)` and
`(+i, +i \infty)`. In general,

.. math ::

    \tan^{-1}(x) = \frac{i}{2}\left(\log(1-ix)-\log(1+ix)\right)

where the principal-branch log is implied.
"""

acot.__doc__ = r"""Computes the inverse cotangent of `x`,
`\mathrm{cot}^{-1}(x) = \tan^{-1}(1/x)`."""

asec.__doc__ = r"""Computes the inverse secant of `x`,
`\mathrm{sec}^{-1}(x) = \cos^{-1}(1/x)`."""

acsc.__doc__ = r"""Computes the inverse cosecant of `x`,
`\mathrm{csc}^{-1}(x) = \sin^{-1}(1/x)`."""

sinpi.__doc__ = """
Computes `\sin(\pi x)`, more accurately than the expression
``sin(pi*x)``::

    >>> from sympy.mpmath import *
    >>> mp.dps = 15
    >>> print sinpi(10**10), sin(pi*(10**10))
    0.0 -2.23936276195592e-6
    >>> print sinpi(10**10+0.5), sin(pi*(10**10+0.5))
    1.0 0.999999999998721
"""

cospi.__doc__ = """
Computes `\cos(\pi x)`, more accurately than the expression
``cos(pi*x)``::

    >>> from sympy.mpmath import *
    >>> mp.dps = 15
    >>> print cospi(10**10), cos(pi*(10**10))
    1.0 0.999999999997493
    >>> print cospi(10**10+0.5), cos(pi*(10**10+0.5))
    0.0 1.59960492420134e-6
"""

@funcwrapper
def sinc(x):
    r"""
    ``sinc(x)`` computes the unnormalized sinc function, defined as

    .. math ::

        \mathrm{sinc}(x) = \begin{cases}
            \sin(x)/x, & \mbox{if } x \ne 0 \\
            1,         & \mbox{if } x = 0.
        \end{cases}

    See :func:`sincpi` for the normalized sinc function.

    Simple values and limits include::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> print sinc(0)
        1.0
        >>> print sinc(1)
        0.841470984807897
        >>> print sinc(inf)
        0.0

    The integral of the sinc function is the sine integral Si::

        >>> print quad(sinc, [0, 1])
        0.946083070367183
        >>> print si(1)
        0.946083070367183
    """
    if isinf(x):
        return 1/x
    if not x:
        return x+1
    return sin(x)/x

@funcwrapper
def sincpi(x):
    r"""
    ``sincpi(x)`` computes the normalized sinc function, defined as

    .. math ::

        \mathrm{sinc}_{\pi}(x) = \begin{cases}
            \sin(\pi x)/(\pi x), & \mbox{if } x \ne 0 \\
            1,                   & \mbox{if } x = 0.
        \end{cases}

    Equivalently, we have
    `\mathrm{sinc}_{\pi}(x) = \mathrm{sinc}(\pi x)`.

    The normalization entails that the function integrates
    to unity over the entire real line::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> print quadosc(sincpi, [-inf, inf], period=2.0)
        1.0

    Like, :func:`sinpi`, :func:`sincpi` is evaluated accurately
    at its roots::

        >>> print sincpi(10)
        0.0
    """
    if isinf(x):
        return 1/x
    if not x:
        return x+1
    return sinpi(x)/(pi*x)

floor = mpfunc('floor', libmpf.mpf_floor, libmpc.mpc_floor, "")
floor.__doc__ = r"""
Computes the floor of `x`, `\lfloor x \rfloor`, defined as
the largest integer less than or equal to `x`::

    >>> from sympy.mpmath import *
    >>> print floor(3.5)
    3.0

Note: :func:`floor` returns a floating-point number, not a
Python ``int``. If `\lfloor x \rfloor` is too large to be
represented exactly at the present working precision, the
result will be rounded, not necessarily in the floor
direction."""

ceil = mpfunc('ceil', libmpf.mpf_ceil, libmpc.mpc_ceil, "")
ceil.__doc__ = r"""
Computes the ceiling of `x`, `\lceil x \rceil`, defined as
the smallest integer greater than or equal to `x`::

    >>> from sympy.mpmath import *
    >>> print ceil(3.5)
    4.0

Note: :func:`ceil` returns a floating-point number, not a
Python ``int``. If `\lceil x \rceil` is too large to be
represented exactly at the present working precision, the
result will be rounded, not necessarily in the ceiling
direction."""

@funcwrapper
def nthroot(x, n):
    r"""
    ``nthroot(x, n)`` computes the principal `n`-th root of `x`,
    `x^{1/n}`. Here `n` must be an integer, and can be negative
    (`x^{-1/n}` is `1/x^{1/n}`).

    For `n = 2` or `n = 3`, using this function is equivalent to
    calling :func:`sqrt` or :func:`cbrt`. In general,
    ``nthroot(x, n)`` is defined to compute `\exp(\log(x)/n)`.

    :func:`nthroot` is implemented to use Newton's method for small
    `n`. At high precision, this makes `x^{1/n}` not much more
    expensive than the regular exponentiation, `x^n`. For very large
    `n`, :func:`nthroot` falls back to use the exponential function.

    :func:`nthroot` is faster and more accurate than raising to a
    floating-point fraction::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> 16807 ** (mpf(1)/5)
        mpf('7.0000000000000009')
        >>> nthroot(16807, 5)
        mpf('7.0')

    """
    n = int(n)
    if isinstance(x, mpf):
        try:
            return make_mpf(libelefun.mpf_nthroot(x._mpf_, n, *prec_rounding))
        except ComplexResult:
            if mp.trap_complex:
                raise
            x = (x._mpf_, libmpf.fzero)
    else:
        x = x._mpc_
    return make_mpc(libmpc.mpc_nthroot(x, n, *prec_rounding))

def hypot(x, y):
    r"""
    Computes the Euclidean norm of the vector `(x, y)`, equal
    to `\sqrt{x^2 + y^2}`. Both `x` and `y` must be real."""
    x = mpmathify(x)
    y = mpmathify(y)
    return make_mpf(libmpf.mpf_hypot(x._mpf_, y._mpf_, *prec_rounding))

def ldexp(x, n):
    r"""
    Computes `x 2^n` efficiently. No rounding is performed.
    The argument `x` must be a real floating-point number (or
    possible to convert into one) and `n` must be a Python ``int``.

        >>> from sympy.mpmath import *
        >>> ldexp(1, 10)
        mpf('1024.0')
        >>> ldexp(1, -3)
        mpf('0.125')

    """
    x = mpmathify(x)
    return make_mpf(libmpf.mpf_shift(x._mpf_, n))

def frexp(x):
    r"""
    Given a real number `x`, returns `(y, n)` with `y \in [0.5, 1)`,
    `n` a Python integer, and such that `x = y 2^n`. No rounding is
    performed.

        >>> from sympy.mpmath import *
        >>> frexp(7.5)
        (mpf('0.9375'), 3)

    """
    x = mpmathify(x)
    y, n = libmpf.mpf_frexp(x._mpf_)
    return make_mpf(y), n

def sign(x):
    r"""
    Returns the sign of `x`, defined as `\mathrm{sign}(x) = x / |x|`
    (with the special case `\sign(0) = 0`)::

        >>> from sympy.mpmath import *
        >>> sign(10)
        mpf('1.0')
        >>> sign(-10)
        mpf('-1.0')
        >>> sign(0)
        mpf('0.0')

    Note that the sign function is also defined for complex numbers,
    for which it gives the projection onto the unit circle::

        >>> mp.dps = 15
        >>> print sign(1+j)
        (0.707106781186547 + 0.707106781186547j)

    """
    x = mpmathify(x)
    if not x or isnan(x):
        return x
    if isinstance(x, mpf):
        return mpf(cmp(x, 0))
    return x / abs(x)

@extraprec(5)
def arg(x):
    r"""
    Computes the complex argument (phase) of `x`, defined as the
    signed angle between the positive real axis and `x` in the
    complex plane::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> print arg(3)
        0.0
        >>> print arg(3+3j)
        0.785398163397448
        >>> print arg(3j)
        1.5707963267949
        >>> print arg(-3)
        3.14159265358979
        >>> print arg(-3j)
        -1.5707963267949

    The angle is defined to satisfy `-\pi < \arg(x) \le \pi` and
    with the sign convention that a nonnegative imaginary part
    results in a nonnegative argument.

    The value returned by :func:`arg` is an ``mpf`` instance.
    """
    x = mpc(x)
    return atan2(x.imag, x.real)

def fabs(x):
    r"""
    Returns the absolute value of `x`, `|x|`. Unlike :func:`abs`,
    :func:`fabs` converts non-mpmath numbers (such as ``int``)
    into mpmath numbers::

        >>> from sympy.mpmath import *
        >>> fabs(3)
        mpf('3.0')
        >>> fabs(-3)
        mpf('3.0')
        >>> fabs(3+4j)
        mpf('5.0')

    """
    return abs(mpmathify(x))

def re(x):
    r"""
    Returns the real part of `x`, `\Re(x)`. Unlike ``x.real``,
    :func:`re` converts `x` to a mpmath number::

        >>> from sympy.mpmath import *
        >>> re(3)
        mpf('3.0')
        >>> re(-1+4j)
        mpf('-1.0')

    """
    return mpmathify(x).real

def im(x):
    r"""
    Returns the imaginary part of `x`, `\Im(x)`. Unlike ``x.imag``,
    :func:`im` converts `x` to a mpmath number::

        >>> from sympy.mpmath import *
        >>> im(3)
        mpf('0.0')
        >>> im(-1+4j)
        mpf('4.0')

    """
    return mpmathify(x).imag

def conj(x):
    r"""
    Returns the complex conjugate of `x`, `\overline{x}`. Unlike
    ``x.conjugate()``, :func:`im` converts `x` to a mpmath number::

        >>> from sympy.mpmath import *
        >>> conj(3)
        mpf('3.0')
        >>> conj(-1+4j)
        mpc(real='-1.0', imag='-4.0')

    """
    return mpmathify(x).conjugate()

def log(x, b=None):
    r"""
    Computes the base-`b` logarithm of `x`, `\log_b(x)`. If `b` is
    unspecified, :func:`log` computes the natural (base `e`) logarithm
    and is equivalent to :func:`ln`. In general, the base `b` logarithm
    is defined in terms of the natural logarithm as
    `\log_b(x) = \ln(x)/\ln(b)`.

    By convention, we take `\log(0) = -\infty`.

    The natural logarithm is real if `x > 0` and complex if `x < 0` or if
    `x` is complex. The principal branch of the complex logarithm is
    used, meaning that `\Im(\ln(x)) = -\pi < \arg(x) \le \pi`.

    **Examples**

    Some basic values and limits::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> print log(1)
        0.0
        >>> print log(2)
        0.693147180559945
        >>> print log(1000,10)
        3.0
        >>> print log(4, 16)
        0.5
        >>> print log(j)
        (0.0 + 1.5707963267949j)
        >>> print log(-1)
        (0.0 + 3.14159265358979j)
        >>> print log(0)
        -inf
        >>> print log(inf)
        +inf

    The natural logarithm is the antiderivative of `1/x`::

        >>> print quad(lambda x: 1/x, [1, 5])
        1.6094379124341
        >>> print log(5)
        1.6094379124341
        >>> print diff(log, 10)
        0.1

    The Taylor series expansion of the natural logarithm around
    `x = 1` has coefficients `(-1)^{n+1}/n`::

        >>> nprint(taylor(log, 1, 7))
        [0.0, 1.0, -0.5, 0.333333, -0.25, 0.2, -0.166667, 0.142857]

    :func:`log` supports arbitrary precision evaluation::

        >>> mp.dps = 50
        >>> print log(pi)
        1.1447298858494001741434273513530587116472948129153
        >>> print log(pi, pi**3)
        0.33333333333333333333333333333333333333333333333333
        >>> mp.dps = 25
        >>> print log(3+4j)
        (1.609437912434100374600759 + 0.9272952180016122324285125j)

    """
    if b is None:
        return ln(x)
    wp = mp.prec + 20
    return ln(x, prec=wp) / ln(b, prec=wp)

def log10(x):
    r"""
    Computes the base-10 logarithm of `x`, `\log_{10}(x)`. ``log10(x)``
    is equivalent to ``log(x, 10)``.
    """
    return log(x, 10)

def power(x, y):
    r"""
    Converts `x` and `y` to mpmath numbers and evaluates
    `x^y = \exp(y \log(x))`::

        >>> from sympy.mpmath import *
        >>> mp.dps = 30
        >>> print power(2, 0.5)
        1.41421356237309504880168872421

    This shows the leading few digits of a large Mersenne prime
    (performing the exact calculation ``2**43112609-1`` and
    displaying the result in Python would be very slow)::

        >>> print power(2, 43112609)-1
        3.16470269330255923143453723949e+12978188

    """
    return mpmathify(x) ** mpmathify(y)

def modf(x,y):
    r"""
    Converts `x` and `y` to mpmath numbers and returns `x \mod y`.
    For mpmath numbers, this is equivalent to ``x % y``.

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> print modf(100, pi)
        2.61062773871641

    You can use :func:`modf` to compute fractional parts of numbers::

        >>> print modf(10.25, 1)
        0.25

    """
    x = mpmathify(x)
    y = mpmathify(y)
    return x % y

def degrees(x):
    r"""
    Converts the radian angle `x` to a degree angle::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> print degrees(pi/3)
        60.0
    """
    return x / degree

def radians(x):
    r"""
    Converts the degree angle `x` to radians::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> print radians(60)
        1.0471975511966
    """
    return x * degree

def atan2(y, x):
    r"""
    Computes the two-argument arctangent, `\mathrm{atan2}(y, x)`,
    giving the signed angle between the positive `x`-axis and the
    point `(x, y)` in the 2D plane. This function is defined for
    real `x` and `y` only.

    The two-argument arctangent essentially computes
    `\mathrm{atan}(y/x)`, but accounts for the signs of both
    `x` and `y` to give the angle for the correct quadrant. The
    following examples illustrate the difference::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> print atan2(1,1), atan(1/1.)
        0.785398163397448 0.785398163397448
        >>> print atan2(1,-1), atan(1/-1.)
        2.35619449019234 -0.785398163397448
        >>> print atan2(-1,1), atan(-1/1.)
        -0.785398163397448 -0.785398163397448
        >>> print atan2(-1,-1), atan(-1/-1.)
        -2.35619449019234 0.785398163397448

    The angle convention is the same as that used for the complex
    argument; see :func:`arg`.
    """
    x = mpmathify(x)
    y = mpmathify(y)
    return make_mpf(libelefun.mpf_atan2(y._mpf_, x._mpf_, *prec_rounding))


fib = fibonacci = mpfunc('fibonacci', libelefun.mpf_fibonacci, libmpc.mpc_fibonacci, "")

fibonacci.__doc__ = r"""
``fibonacci(n)`` computes the `n`-th Fibonacci number, `F(n)`. The
Fibonacci numbers are defined by the recurrence `F(n) = F(n-1) + F(n-2)`
with the initial values `F(0) = 0`, `F(1) = 1`. :func:`fibonacci`
extends this definition to arbitrary real and complex arguments
using the formula

.. math ::

  F(z) = \frac{\phi^z - \cos(\pi z) \phi^{-z}}{\sqrt 5}

where `\phi` is the golden ratio. :func:`fibonacci` also uses this
continuous formula to compute `F(n)` for extremely large `n`, where
calculating the exact integer would be wasteful.

For convenience, :func:`fib` is available as an alias for
:func:`fibonacci`.

**Basic examples**

Some small Fibonacci numbers are::

    >>> from sympy.mpmath import *
    >>> mp.dps = 15
    >>> for i in range(10):
    ...     print fibonacci(i),
    ...
    0.0 1.0 1.0 2.0 3.0 5.0 8.0 13.0 21.0 34.0

    >>> print fibonacci(50)
    12586269025.0

The recurrence for `F(n)` extends backwards to negative `n`::

    >>> for i in range(10):
    ...     print fibonacci(-i),
    ...
    0.0 1.0 -1.0 2.0 -3.0 5.0 -8.0 13.0 -21.0 34.0

Large Fibonacci numbers will be computed approximately unless
the precision is set high enough::

    >>> print fib(200)
    2.8057117299251e+41
    >>> mp.dps = 45
    >>> print fib(200)
    280571172992510140037611932413038677189525.0

:func:`fibonacci` can compute approximate Fibonacci numbers
of stupendous size::

    >>> mp.dps = 15
    >>> print fibonacci(10**25)
    3.49052338550226e+2089876402499787337692720

**Real and complex arguments**

The extended Fibonacci function is an analytic function. The
property `F(z) = F(z-1) + F(z-2)` holds for arbitrary `z`::

    >>> mp.dps = 15
    >>> print fib(pi)
    2.1170270579161
    >>> print fib(pi-1) + fib(pi-2)
    2.1170270579161
    >>> print fib(3+4j)
    (-5248.51130728372 - 14195.962288353j)
    >>> print fib(2+4j) + fib(1+4j)
    (-5248.51130728372 - 14195.962288353j)

The Fibonacci function has infinitely many roots on the
negative half-real axis. The first root is at 0, the second is
close to -0.18, and then there are infinitely many roots that
asymptotically approach `-n+1/2`::

    >>> print findroot(fib, -0.2)
    -0.183802359692956
    >>> print findroot(fib, -2)
    -1.57077646820395
    >>> print findroot(fib, -17)
    -16.4999999596115
    >>> print findroot(fib, -24)
    -23.5000000000479

**Mathematical relationships**

For large `n`, `F(n+1)/F(n)` approaches the golden ratio::

    >>> mp.dps = 50
    >>> print fibonacci(101)/fibonacci(100)
    1.6180339887498948482045868343656381177203127439638
    >>> print phi
    1.6180339887498948482045868343656381177203091798058

The sum of reciprocal Fibonacci numbers converges to an irrational
number for which no closed form expression is known::

    >>> mp.dps = 15
    >>> print nsum(lambda n: 1/fib(n), [1, inf])
    3.35988566624318

Amazingly, however, the sum of odd-index reciprocal Fibonacci
numbers can be expressed in terms of a Jacobi theta function::

    >>> print nsum(lambda n: 1/fib(2*n+1), [0, inf])
    1.82451515740692
    >>> print sqrt(5)*jtheta(2,0,(3-sqrt(5))/2)**2/4
    1.82451515740692

Some related sums can be done in closed form::

    >>> print nsum(lambda k: 1/(1+fib(2*k+1)), [0, inf])
    1.11803398874989
    >>> print phi - 0.5
    1.11803398874989
    >>> f = lambda k:(-1)**(k+1) / sum(fib(n)**2 for n in range(1,k+1))
    >>> print nsum(f, [1, inf])
    0.618033988749895
    >>> print phi-1
    0.618033988749895

**References**

1. http://mathworld.wolfram.com/FibonacciNumber.html
"""


zeta = mpfunc('zeta', gammazeta.mpf_zeta, gammazeta.mpc_zeta, 'Riemann zeta function')
altzeta = mpfunc('zeta', gammazeta.mpf_altzeta, gammazeta.mpc_altzeta, 'Dirichlet eta function')

zeta.__doc__ = r"""
    ``zeta(s)`` computes the Riemann zeta function, `\zeta(s)`.
    The Riemann zeta function is defined for `\Re(s) > 1` by

    .. math ::

      \zeta(s) = 1+\frac{1}{2^s}+\frac{1}{3^s}+\frac{1}{4^s}+\ldots

    and for `\Re(s) \le 1` by analytic continuation. It has a pole
    at `s = 1`.

    **Examples**

    Some exact values of the zeta function are::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> print zeta(2)
        1.64493406684823
        >>> print pi**2 / 6
        1.64493406684823
        >>> print zeta(0)
        -0.5
        >>> print zeta(-1)
        -0.0833333333333333
        >>> print zeta(-2)
        0.0

    :func:`zeta` supports arbitrary precision evaluation and
    complex arguments::

        >>> mp.dps = 50
        >>> print zeta(pi)
        1.1762417383825827588721504519380520911697389900217
        >>> print zeta(1+2j)  # doctest: +NORMALIZE_WHITESPACE
        (0.5981655697623817367034568491742186771747764868876 -
        0.35185474521784529049653859679690026505229177886045j)

    The Riemann zeta function has so-called nontrivial zeros on
    the critical line `s = 1/2 + it`::

        >>> mp.dps = 15
        >>> print findroot(zeta, 0.5+14j)
        (0.5 + 14.1347251417347j)
        >>> print findroot(zeta, 0.5+21j)
        (0.5 + 21.0220396387716j)
        >>> print findroot(zeta, 0.5+25j)
        (0.5 + 25.0108575801457j)

    For investigation of the zeta function zeros, the Riemann-Siegel
    Z-function is often more convenient than working with the Riemann
    zeta function directly (see :func:`siegelz`).

    For large positive `s`, `\zeta(s)` rapidly approaches 1::

        >>> print zeta(30)
        1.00000000093133
        >>> print zeta(100)
        1.0
        >>> print zeta(inf)
        1.0

    The following series converges and in fact has a simple
    closed form value::

        >>> print nsum(lambda k: zeta(k)-1, [2, inf])
        1.0

    **Algorithm**

    The primary algorithm is Borwein's algorithm for the Dirichlet
    eta function. Three separate implementations are used: for general
    real arguments, general complex arguments, and for integers. The
    reflection formula is applied to arguments in the negative
    half-plane. For very large real arguments, either direct
    summation or the Euler prime product is used.

    It should be noted that computation of `\zeta(s)` gets very slow
    when `s` is far away from the real axis.

    **References**

    1. http://mathworld.wolfram.com/RiemannZetaFunction.html

    2. http://www.cecm.sfu.ca/personal/pborwein/PAPERS/P155.pdf
"""

altzeta.__doc__ = r"""
    Computes the Dirichlet eta function, `\eta(s)`, also known as the
    alternating zeta function. This function is defined in analogy
    with the Riemann zeta function as providing the sum of the
    alternating series

    .. math ::

        \eta(s) = 1-\frac{1}{2^s}+\frac{1}{3^s}-\frac{1}{4^s}+\ldots

    Note that `\eta(1) = \log(2)` is the alternating harmonic series.
    The eta function unlike the Riemann zeta function is an entire
    function, having a finite value for all complex `s`.

    The alternating and non-alternating zeta functions are related
    via the simple formula

    .. math ::

        \eta(s) = (1 - 2^{1-s}) \zeta(s).

    This formula can be used to define `\eta(s)` for `\Re(s) \le 0`,
    where the series diverges.

    **Examples**

    Some special values are::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> print altzeta(1)
        0.693147180559945
        >>> print altzeta(0)
        0.5
        >>> print altzeta(-1)
        0.25
        >>> print altzeta(-2)
        0.0

    An example of a sum that can be computed more accurately and
    efficiently via :func:`altzeta` than via numerical summation::

        >>> sum(-(-1)**n / n**2.5 for n in range(1, 100))
        0.86720495150398402
        >>> print altzeta(2.5)
        0.867199889012184

    At positive even integers, the Dirichlet eta function
    evaluates to a rational multiple of a power of `\pi`::

        >>> print altzeta(2)
        0.822467033424113
        >>> print pi**2/12
        0.822467033424113

    Like the Riemann zeta function, `\eta(s)`, approaches 1
    as `s` approaches positive infinity, although it does
    so from below rather than from above::

        >>> print altzeta(30)
        0.999999999068682
        >>> print altzeta(inf)
        1.0
        >>> altzeta(1000, rounding='d')
        mpf('0.99999999999999989')
        >>> altzeta(1000, rounding='u')
        mpf('1.0')

    **References**

    1. http://mathworld.wolfram.com/DirichletEtaFunction.html

    2. http://en.wikipedia.org/wiki/Dirichlet_eta_function

"""

gamma = mpfunc('gamma', gammazeta.mpf_gamma, gammazeta.mpc_gamma, "gamma function")
factorial = mpfunc('factorial', gammazeta.mpf_factorial, gammazeta.mpc_factorial, "factorial")
fac = factorial

factorial.__doc__ = r"""
Computes the factorial, `x!`. For integers `n \ge 0`, we have
`n! = 1 \cdot 2 \cdots (n-1) \cdot n` and more generally the factorial
is defined for real or complex `x` by `x! = \Gamma(x+1)`.

**Examples**

Basic values and limits::

    >>> from sympy.mpmath import *
    >>> mp.dps = 15
    >>> for k in range(6):
    ...     print k, fac(k)
    ...
    0 1.0
    1 1.0
    2 2.0
    3 6.0
    4 24.0
    5 120.0
    >>> print fac(inf)
    +inf
    >>> print fac(0.5), sqrt(pi)/2
    0.886226925452758 0.886226925452758

For large positive `x`, `x!` can be approximated by
Stirling's formula::

    >>> x = 10**10
    >>> print fac(x)
    2.32579620567308e+95657055186
    >>> print sqrt(2*pi*x)*(x/e)**x
    2.32579597597705e+95657055186

:func:`fac` supports evaluation for astronomically large values::

    >>> print fac(10**30)
    6.22311232304258e+29565705518096748172348871081098

Reciprocal factorials appear in the Taylor series of the
exponential function (among many other contexts)::

    >>> print nsum(lambda k: 1/fac(k), [0, inf]), exp(1)
    2.71828182845905 2.71828182845905
    >>> print nsum(lambda k: pi**k/fac(k), [0, inf]), exp(pi)
    23.1406926327793 23.1406926327793

"""

gamma.__doc__ = r"""
Computes the gamma function, `\Gamma(x)`. The gamma function is a
shifted version of the ordinary factorial, satisfying
`\Gamma(n) = (n-1)!` for integers `n > 0`. More generally, it
is defined by

.. math ::

    \Gamma(x) = \int_0^{\infty} t^{x-1} e^{-t}\, dt

for any real or complex `x` with `\Re(x) > 0` and for `\Re(x) < 0`
by analytic continuation.

**Examples**

Basic values and limits::

    >>> from sympy.mpmath import *
    >>> mp.dps = 15
    >>> for k in range(1, 6):
    ...     print k, gamma(k)
    ...
    1 1.0
    2 1.0
    3 2.0
    4 6.0
    5 24.0
    >>> print gamma(inf)
    +inf
    >>> print gamma(0)
    Traceback (most recent call last):
      ...
    ValueError: gamma function pole

The gamma function of a half-integer is a rational multiple of
`\sqrt{\pi}`::

    >>> print gamma(0.5), sqrt(pi)
    1.77245385090552 1.77245385090552
    >>> print gamma(1.5), sqrt(pi)/2
    0.886226925452758 0.886226925452758

We can check the integral definition::

    >>> print gamma(3.5)
    3.32335097044784
    >>> print quad(lambda t: t**2.5*exp(-t), [0,inf])
    3.32335097044784

:func:`gamma` supports arbitrary-precision evaluation and
complex arguments::

    >>> mp.dps = 50
    >>> print gamma(sqrt(3))
    0.91510229697308632046045539308226554038315280564184
    >>> mp.dps = 25
    >>> print gamma(2j)
    (0.009902440080927490985955066 - 0.07595200133501806872408048j)

Arguments can also be large. Note that the gamma function grows
very quickly::

    >>> mp.dps = 15
    >>> print gamma(10**20)
    6.33636415517321e+1956570547910964391727

"""

def psi(m, z):
    r"""
    Gives the polygamma function of order `m` of `z`, `\psi^{(m)}(z)`.
    Special cases are known as the *digamma function* (`\psi^{(0)}(z)`),
    the *trigamma function* (`\psi^{(1)}(z)`), etc. The polygamma
    functions are defined as the logarithmic derivatives of the gamma
    function:

    .. math ::

        \psi^{(m)}(z) = \left(\frac{d}{dz}\right)^{m+1} \log \Gamma(z)

    In particular, `\psi^{(0)}(z) = \Gamma'(z)/\Gamma(z)`. In the
    present implementation of :func:`psi`, the order `m` must be a
    nonnegative integer, while the argument `z` may be an arbitrary
    complex number (with exception for the polygamma function's poles
    at `z = 0, -1, -2, \ldots`).

    **Examples**

    For various rational arguments, the polygamma function reduces to
    a combination of standard mathematical constants::

        >>> from sympy.mpmath import *
        >>> mp.dps = 25
        >>> print psi(0, 1), -euler
        -0.5772156649015328606065121 -0.5772156649015328606065121
        >>> print psi(1, '1/4'), pi**2+8*catalan
        17.19732915450711073927132 17.19732915450711073927132
        >>> print psi(2, '1/2'), -14*apery
        -16.82879664423431999559633 -16.82879664423431999559633

    The polygamma functions are derivatives of each other::

        >>> print diff(lambda x: psi(3, x), pi), psi(4, pi)
        -0.1105749312578862734526952 -0.1105749312578862734526952
        >>> print quad(lambda x: psi(4, x), [2, 3]), psi(3,3)-psi(3,2)
        -0.375 -0.375

    The digamma function diverges logarithmically as `z \to \infty`,
    while higher orders tend to zero::

        >>> print psi(0,inf), psi(1,inf), psi(2,inf)
        +inf 0.0 0.0

    Evaluation for a complex argument::

        >>> print psi(2, -1-2j)
        (0.03902435405364952654838445 + 0.1574325240413029954685366j)

    Evaluation is supported for large orders `m` and/or large
    arguments `z`::

        >>> print psi(3, 10**100)
        2.0e-300
        >>> print psi(250, 10**30+10**20*j)
        (-1.293142504363642687204865e-7010 + 3.232856260909107391513108e-7018j)

    **Application to infinite series**

    Any infinite series where the summand is a rational function of
    the index `k` can be evaluated in closed form in terms of polygamma
    functions of the roots and poles of the summand::

        >>> a = sqrt(2)
        >>> b = sqrt(3)
        >>> print nsum(lambda k: 1/((k+a)**2*(k+b)), [0, inf])
        0.4049668927517857061917531
        >>> print (psi(0,a)-psi(0,b)-a*psi(1,a)+b*psi(1,a))/(a-b)**2
        0.4049668927517857061917531

    This follows from the series representation (`m > 0`)

    .. math ::

        \psi^{(m)}(z) = (-1)^{m+1} m! \sum_{k=0}^{\infty}
            \frac{1}{(z+k)^{m+1}}.

    Since the roots of a polynomial may be complex, it is sometimes
    necessary to use the complex polygamma function to evaluate
    an entirely real-valued sum::

        >>> print nsum(lambda k: 1/(k**2-2*k+3), [0, inf])
        1.694361433907061256154665
        >>> nprint(polyroots([1,-2,3]))
        [(1.0 - 1.41421j), (1.0 + 1.41421j)]
        >>> r1 = 1-sqrt(2)*j
        >>> r2 = r1.conjugate()
        >>> print (psi(0,-r2)-psi(0,-r1))/(r1-r2)
        (1.694361433907061256154665 + 0.0j)

    """
    z = mpmathify(z)
    m = int(m)
    if isinstance(z, mpf):
        return make_mpf(gammazeta.mpf_psi(m, z._mpf_, *prec_rounding))
    else:
        return make_mpc(gammazeta.mpc_psi(m, z._mpc_, *prec_rounding))

def psi0(z):
    """Shortcut for psi(0,z) (the digamma function)"""
    return psi(0, z)

def psi1(z):
    """Shortcut for psi(1,z) (the trigamma function)"""
    return psi(1, z)

def psi2(z):
    """Shortcut for psi(2,z) (the tetragamma function)"""
    return psi(2, z)

def psi3(z):
    """Shortcut for psi(3,z) (the pentagamma function)"""
    return psi(3, z)

polygamma = psi
digamma = psi0
trigamma = psi1
tetragamma = psi2
pentagamma = psi3

harmonic = mpfunc('harmonic', gammazeta.mpf_harmonic, gammazeta.mpc_harmonic,
    "nth harmonic number")

harmonic.__doc__ = r"""
    If `n` is an integer, ``harmonic(n)`` gives a floating-point
    approximation of the `n`-th harmonic number `H(n)`, defined as

    .. math ::

        H(n) = 1 + \frac{1}{2} + \frac{1}{3} + \ldots + \frac{1}{n}

    The firrst few harmonic numbers are::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> for n in range(8):
        ...     print n, harmonic(n)
        ...
        0 0.0
        1 1.0
        2 1.5
        3 1.83333333333333
        4 2.08333333333333
        5 2.28333333333333
        6 2.45
        7 2.59285714285714

    The infinite harmonic series `1 + 1/2 + 1/3 + \ldots` diverges::

        >>> print harmonic(inf)
        +inf

    :func:`harmonic` is evaluated using the digamma function rather
    than by summing the harmonic series term by term. It can therefore
    be computed quickly for arbitrarily large `n`, and even for
    nonintegral arguments::

        >>> print harmonic(10**100)
        230.835724964306
        >>> print harmonic(0.5)
        0.613705638880109
        >>> print harmonic(3+4j)
        (2.24757548223494 + 0.850502209186044j)

    :func:`harmonic` supports arbitrary precision evaluation::

        >>> mp.dps = 50
        >>> print harmonic(11)
        3.0198773448773448773448773448773448773448773448773
        >>> print harmonic(pi)
        1.8727388590273302654363491032336134987519132374152

    The harmonic series diverges, but at a glacial pace. It is possible
    to calculate the exact number of terms required before the sum
    exceeds a given amount, say 100::

        >>> mp.dps = 50
        >>> v = 10**findroot(lambda x: harmonic(10**x) - 100, 10)
        >>> print v
        15092688622113788323693563264538101449859496.864101
        >>> v = int(ceil(v))
        >>> print v
        15092688622113788323693563264538101449859497
        >>> print harmonic(v-1)
        99.999999999999999999999999999999999999999999942747
        >>> print harmonic(v)
        100.000000000000000000000000000000000000000000009

    """

def bernoulli(n):
    r"""
    Computes the nth Bernoulli number, `B_n`, for any integer `n \ge 0`.

    The Bernoulli numbers are rational numbers, but this function
    returns a floating-point approximation. To obtain an exact
    fraction, use :func:`bernfrac` instead.

    **Examples**

    Numerical values of the first few Bernoulli numbers::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> for n in range(15):
        ...     print n, bernoulli(n)
        ...
        0 1.0
        1 -0.5
        2 0.166666666666667
        3 0.0
        4 -0.0333333333333333
        5 0.0
        6 0.0238095238095238
        7 0.0
        8 -0.0333333333333333
        9 0.0
        10 0.0757575757575758
        11 0.0
        12 -0.253113553113553
        13 0.0
        14 1.16666666666667

    Bernoulli numbers can be approximated with arbitrary precision::

        >>> mp.dps = 50
        >>> print bernoulli(100)
        -2.8382249570693706959264156336481764738284680928013e+78

    Arbitrarily large `n` are supported::

        >>> mp.dps = 15
        >>> print bernoulli(10**20 + 2)
        3.09136296657021e+1876752564973863312327

    The Bernoulli numbers are related to the Riemann zeta function
    at integer arguments::

        >>> print -bernoulli(8) * (2*pi)**8 / (2*fac(8))
        1.00407735619794
        >>> print zeta(8)
        1.00407735619794

    **Algorithm**

    For small `n` (`n < 3000`) :func:`bernoulli` uses a recurrence
    formula due to Ramanujan. All results in this range are cached,
    so sequential computation of small Bernoulli numbers is
    guaranteed to be fast.

    For larger `n`, `B_n` is evaluated in terms of the Riemann zeta
    function.
    """
    return make_mpf(gammazeta.mpf_bernoulli(int(n), *prec_rounding))

bernfrac = gammazeta.bernfrac

stieltjes_cache = {}

def stieltjes(n, a=1):
    r"""
    For a nonnegative integer `n`, ``stieltjes(n)`` computes the
    `n`-th Stieltjes constant `\gamma_n`, defined as the
    `n`-th coefficient in the Laurent series expansion of the
    Riemann zeta function around the pole at `s = 1`. That is,
    we have:

    .. math ::

      \zeta(s) = \frac{1}{s-1} \sum_{n=0}^{\infty}
          \frac{(-1)^n}{n!} \gamma_n (s-1)^n

    More generally, ``stieltjes(n, a)`` gives the corresponding
    coefficient `\gamma_n(a)` for the Hurwitz zeta function
    `\zeta(s,a)` (with `\gamma_n = \gamma_n(1)`).

    **Examples**

    The zeroth Stieltjes constant is just Euler's constant `\gamma`::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> print stieltjes(0)
        0.577215664901533

    Some more values are::

        >>> print stieltjes(1)
        -0.0728158454836767
        >>> print stieltjes(10)
        0.000205332814909065
        >>> print stieltjes(30)
        0.00355772885557316
        >>> print stieltjes(1000)
        -1.57095384420474e+486
        >>> print stieltjes(2000)
        2.680424678918e+1109
        >>> print stieltjes(1, 2.5)
        -0.23747539175716

    An alternative way to compute `\gamma_1`::

        >>> print diff(extradps(25)(lambda x: 1/(x-1) - zeta(x)), 1)
        -0.0728158454836767

    :func:`stieltjes` supports arbitrary precision evaluation::

        >>> mp.dps = 50
        >>> print stieltjes(2)
        -0.0096903631928723184845303860352125293590658061013408

    **Algorithm**

    :func:`stieltjes` numerically evaluates the integral in
    the following representation due to Ainsworth, Howell and
    Coffey [1], [2]:

    .. math ::

      \gamma_n(a) = \frac{\log^n a}{2a} - \frac{\log^{n+1}(a)}{n+1} +
          \frac{2}{a} \Re \int_0^{\infty}
          \frac{(x/a-i)\log^n(a-ix)}{(1+x^2/a^2)(e^{2\pi x}-1)} dx.

    For some reference values with `a = 1`, see e.g. [4].

    **References**

    1. O. R. Ainsworth & L. W. Howell, "An integral representation of
       the generalized Euler-Mascheroni constants", NASA Technical
       Paper 2456 (1985),
       http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19850014994_1985014994.pdf

    2. M. W. Coffey, "The Stieltjes constants, their relation to the
       `\eta_j` coefficients, and representation of the Hurwitz
       zeta function", 	arXiv:0706.0343v1 http://arxiv.org/abs/0706.0343

    3. http://mathworld.wolfram.com/StieltjesConstants.html

    4. http://pi.lacim.uqam.ca/piDATA/stieltjesgamma.txt

    """
    n = mpmathify(n)
    a = mpmathify(a)
    if n < 0:
        raise ValueError("Stieltjes constants defined for n >= 0")
    if a == 1:
        if n == 0:
            return +euler
        if n in stieltjes_cache:
            prec, s = stieltjes_cache[n]
            if prec >= mp.prec:
                return +s
    mag = 1
    def f(x):
        xa = x/a
        v = (xa-j)*log(a-j*x)**n/(1+xa**2)/(exp(2*pi*x)-1)
        return v.real / mag
    from quadrature import quad
    orig = mp.prec
    try:
        # Normalize integrand by approx. magnitude to
        # speed up quadrature (which uses absolute error)
        if n > 50:
            mp.prec = 20
            mag = quad(f, [0,inf], maxdegree=3)
        mp.prec = orig + 10 + int(n**0.5)
        s = quad(f, [0,inf], maxdegree=20)
        v = log(a)**n/(2*a) - log(a)**(n+1)/(n+1) + 2*s/a*mag
    finally:
        mp.prec = orig
    if a == 1 and isint(n):
        stieltjes_cache[n] = (mp.prec, v)
    return +v

def isnpint(x):
    if not x:
        return True
    if isinstance(x, mpf):
        sign, man, exp, bc = x._mpf_
        return sign and exp >= 0
    if isinstance(x, mpc):
        return not x.imag and isnpint(x.real)

def gammaprod(a, b):
    r"""
    Given iterables `a` and `b`, ``gammaprod(a, b)`` computes the
    product / quotient of gamma functions:

    .. math ::

        \frac{\Gamma(a_0) \Gamma(a_1) \cdots \Gamma(a_p)}
             {\Gamma(b_0) \Gamma(b_1) \cdots \Gamma(b_q)}

    Unlike direct calls to :func:`gamma`, :func:`gammaprod` considers
    the entire product as a limit and evaluates this limit properly if
    any of the numerator or denominator arguments are nonpositive
    integers such that poles of the gamma function are encountered.
    That is, :func:`gammaprod` evaluates

    .. math ::

        \lim_{\epsilon \to 0}
        \frac{\Gamma(a_0+\epsilon) \Gamma(a_1+\epsilon) \cdots
            \Gamma(a_p+\epsilon)}
             {\Gamma(b_0+\epsilon) \Gamma(b_1+\epsilon) \cdots
            \Gamma(b_q+\epsilon)}

    In particular:

    * If there are equally many poles in the numerator and the
      denominator, the limit is a rational number times the remaining,
      regular part of the product.

    * If there are more poles in the numerator, :func:`gammaprod`
      returns ``+inf``.

    * If there are more poles in the denominator, :func:`gammaprod`
      returns 0.

    **Examples**

    The reciprocal gamma function `1/\Gamma(x)` evaluated at `x = 0`::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> gammaprod([], [0])
        mpf('0.0')

    A limit::

        >>> gammaprod([-4], [-3])
        mpf('-0.25')
        >>> limit(lambda x: gamma(x-1)/gamma(x), -3, direction=1)
        mpf('-0.25')
        >>> limit(lambda x: gamma(x-1)/gamma(x), -3, direction=-1)
        mpf('-0.25')

    """
    a = [mpmathify(x) for x in a]
    b = [mpmathify(x) for x in b]
    poles_num = []
    poles_den = []
    regular_num = []
    regular_den = []
    for x in a: [regular_num, poles_num][isnpint(x)].append(x)
    for x in b: [regular_den, poles_den][isnpint(x)].append(x)
    # One more pole in numerator or denominator gives 0 or inf
    if len(poles_num) < len(poles_den): return mpf(0)
    if len(poles_num) > len(poles_den): return mpf('+inf')
    # All poles cancel
    # lim G(i)/G(j) = (-1)**(i+j) * gamma(1-j) / gamma(1-i)
    p = mpf(1)
    orig = mp.prec
    try:
        mp.prec = orig + 15
        while poles_num:
            i = poles_num.pop()
            j = poles_den.pop()
            p *= (-1)**(i+j) * gamma(1-j) / gamma(1-i)
        for x in regular_num: p *= gamma(x)
        for x in regular_den: p /= gamma(x)
    finally:
        mp.prec = orig
    return +p

def beta(x, y):
    r"""
    Computes the beta function,
    `B(x,y) = \Gamma(x) \Gamma(y) / \Gamma(x+y)`.
    The beta function is also commonly defined by the integral
    representation

    .. math ::

      B(x,y) = \int_0^1 t^{x-1} (1-t)^{y-1} \, dt

    **Examples**

    For integer and half-integer arguments where all three gamma
    functions are finite, the beta function becomes either rational
    number or a rational multiple of `\pi`::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> print beta(5, 2)
        0.0333333333333333
        >>> print beta(1.5, 2)
        0.266666666666667
        >>> print 16*beta(2.5, 1.5)
        3.14159265358979

    Where appropriate, :func:`beta` evaluates limits. A pole
    of the beta function is taken to result in ``+inf``::

        >>> print beta(-0.5, 0.5)
        0.0
        >>> print beta(-3, 3)
        -0.333333333333333
        >>> print beta(-2, 3)
        +inf
        >>> print beta(inf, 1)
        0.0
        >>> print beta(inf, 0)
        nan

    :func:`beta` supports complex numbers and arbitrary precision
    evaluation::

        >>> print beta(1, 2+j)
        (0.4 - 0.2j)
        >>> mp.dps = 25
        >>> print beta(j,0.5)
        (1.079424249270925780135675 - 1.410032405664160838288752j)
        >>> mp.dps = 50
        >>> print beta(pi, e)
        0.037890298781212201348153837138927165984170287886464

    Various integrals can be computed by means of the
    beta function::

        >>> mp.dps = 15
        >>> print quad(lambda t: t**2.5*(1-t)**2, [0, 1])
        0.0230880230880231
        >>> print beta(3.5, 3)
        0.0230880230880231
        >>> print quad(lambda t: sin(t)**4 * sqrt(cos(t)), [0, pi/2])
        0.319504062596158
        >>> print beta(2.5, 0.75)/2
        0.319504062596158

    """
    x = mpmathify(x)
    y = mpmathify(y)
    if isinf(y):
        x, y = y, x
    if isinf(x):
        if x == inf and not y.imag:
            if y == -inf:
                return nan
            if y > 0:
                return zero
            if isint(y):
                return nan
            if y < 0:
                return sign(gamma(y)) * inf
        return nan
    return gammaprod([x, y], [x+y])

def binomial(n, k):
    r"""
    Computes the binomial coefficient

    .. math ::

        {n \choose k} = \frac{n!}{k!(n-k)!}.

    The binomial coefficient gives the number of ways that `k` items
    can be chosen from a set of `n` items. More generally, the binomial
    coefficient is a well-defined function of arbitrary real or
    complex `n` and `k`, via the gamma function.

    **Examples**

    Generate Pascal's triangle::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> for n in range(5):
        ...     nprint([binomial(n,k) for k in range(n+1)])
        ...
        [1.0]
        [1.0, 1.0]
        [1.0, 2.0, 1.0]
        [1.0, 3.0, 3.0, 1.0]
        [1.0, 4.0, 6.0, 4.0, 1.0]

    There is 1 way to select 0 items from the empty set, and 0 ways to
    select 1 item from the empty set::

        >>> print binomial(0, 0)
        1.0
        >>> print binomial(0, 1)
        0.0

    :func:`binomial` supports large arguments::

        >>> print binomial(10**20, 10**20-5)
        8.33333333333333e+157
        >>> print binomial(10**20, 10**10)
        2.60784095465201e+104342944813

    Nonintegral binomial coefficients find use in series
    expansions::

        >>> nprint(taylor(lambda x: (1+x)**0.25, 0, 4))
        [1.0, 0.25, -9.375e-2, 5.46875e-2, -3.75977e-2]
        >>> nprint([binomial(0.25, k) for k in range(5)])
        [1.0, 0.25, -9.375e-2, 5.46875e-2, -3.75977e-2]

    An integral representation::

        >>> n, k = 5, 3
        >>> f = lambda t: exp(-j*k*t)*(1+exp(j*t))**n
        >>> print chop(quad(f, [-pi,pi])/(2*pi))
        10.0
        >>> print binomial(n,k)
        10.0

    """
    return gammaprod([n+1], [k+1, n-k+1])

def rf(x, n):
    r"""
    Computes the rising factorial or Pochhammer symbol,

    .. math ::

        x^{(n)} = x (x+1) \cdots (x+n-1) = \frac{\Gamma(x+n)}{\Gamma(x)}

    where the rightmost expression is valid for nonintegral `n`.

    **Examples**

    For integral `n`, the rising factorial is a polynomial::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> for n in range(5):
        ...     nprint(taylor(lambda x: rf(x,n), 0, n))
        ...
        [1.0]
        [0.0, 1.0]
        [0.0, 1.0, 1.0]
        [0.0, 2.0, 3.0, 1.0]
        [0.0, 6.0, 11.0, 6.0, 1.0]

    Evaluation is supported for arbitrary arguments::

        >>> print rf(2+3j, 5.5)
        (-7202.03920483347 - 3777.58810701527j)
    """
    return gammaprod([x+n], [x])

def ff(x, n):
    r"""
    Computes the falling factorial,

    .. math ::

        (x)_n = x (x-1) \cdots (x-n+1) = \frac{\Gamma(x+1)}{\Gamma(x-n+1)}

    where the rightmost expression is valid for nonintegral `n`.

    **Examples**

    For integral `n`, the falling factorial is a polynomial::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> for n in range(5):
        ...     nprint(taylor(lambda x: ff(x,n), 0, n))
        ...
        [1.0]
        [0.0, 1.0]
        [0.0, -1.0, 1.0]
        [0.0, 2.0, -3.0, 1.0]
        [0.0, -6.0, 11.0, -6.0, 1.0]

    Evaluation is supported for arbitrary arguments::

        >>> print ff(2+3j, 5.5)
        (-720.41085888203 + 316.101124983878j)
    """
    return gammaprod([x+1], [x-n+1])

@funcwrapper
def fac2(x):
    r"""
    Computes the double factorial `x!!`, defined for integers
    `x > 0` by

    .. math ::

        x!! = \begin{cases}
            1 \cdot 3 \cdots (x-2) \cdot x & x \;\mathrm{odd} \\
            2 \cdot 4 \cdots (x-2) \cdot x & x \;\mathrm{even}
        \end{cases}

    and more generally by [1]

    .. math ::

        x!! = 2^{x/2} \left(\frac{\pi}{2}\right)^{(\cos(\pi x)-1)/4}
              \Gamma\left(\frac{x}{2}+1\right).

    **Examples**

    The integer sequence of double factorials begins::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> nprint([fac2(n) for n in range(10)])
        [1.0, 1.0, 2.0, 3.0, 8.0, 15.0, 48.0, 105.0, 384.0, 945.0]

    For large `x`, double factorials follow a Stirling-like asymptotic
    approximation::

        >>> x = mpf(10000)
        >>> print fac2(x)
        5.97272691416282e+17830
        >>> print sqrt(pi)*x**((x+1)/2)*exp(-x/2)
        5.97262736954392e+17830

    The recurrence formula `x!! = x (x-2)!!` can be reversed to
    define the double factorial of negative odd integers (but
    not negative even integers)::

        >>> print fac2(-1), fac2(-3), fac2(-5), fac2(-7)
        1.0 -1.0 0.333333333333333 -0.0666666666666667
        >>> fac2(-2)
        Traceback (most recent call last):
          ...
        ValueError: gamma function pole

    With the exception of the poles at negative even integers,
    :func:`fac2` supports evaluation for arbitrary complex arguments.
    The recurrence formula is valid generally::

        >>> print fac2(pi+2j)
        (-1.3697207890154e-12 + 3.93665300979176e-12j)
        >>> print (pi+2j)*fac2(pi-2+2j)
        (-1.3697207890154e-12 + 3.93665300979176e-12j)

    Double factorials should not be confused with nested factorials,
    which are immensely larger::

        >>> print fac(fac(20))
        5.13805976125208e+43675043585825292774
        >>> print fac2(20)
        3715891200.0

    Double factorials appear, among other things, in series expansions
    of Gaussian functions and the error function. Infinite series
    include::

        >>> print nsum(lambda k: 1/fac2(k), [0, inf])
        3.05940740534258
        >>> print sqrt(e)*(1+sqrt(pi/2)*erf(sqrt(2)/2))
        3.05940740534258
        >>> print nsum(lambda k: 2**k/fac2(2*k-1), [1, inf])
        4.06015693855741
        >>> print e * erf(1) * sqrt(pi)
        4.06015693855741

    A beautiful Ramanujan sum::

        >>> print nsum(lambda k: (-1)**k*(fac2(2*k-1)/fac2(2*k))**3, [0,inf])
        0.90917279454693
        >>> print (gamma('9/8')/gamma('5/4')/gamma('7/8'))**2
        0.90917279454693

    **References**

    1. http://functions.wolfram.com/GammaBetaErf/Factorial2/27/01/0002/

    2. http://mathworld.wolfram.com/DoubleFactorial.html

    """
    if isinf(x):
        if x == inf:
            return x
        return nan
    return 2**(x/2)*(pi/2)**((cospi(x)-1)/4)*gamma(x/2+1)


#---------------------------------------------------------------------------#
#                                                                           #
#                          Hypergeometric functions                         #
#                                                                           #
#---------------------------------------------------------------------------#

class _mpq(tuple):
    @property
    def _mpf_(self):
        return (mpf(self[0])/self[1])._mpf_
    def __add__(self, other):
        if isinstance(other, _mpq):
            a, b = self
            c, d = other
            return _mpq((a*d+b*c, b*d))
        return NotImplemented
    def __sub__(self, other):
        if isinstance(other, _mpq):
            a, b = self
            c, d = other
            return _mpq((a*d-b*c, b*d))
        return NotImplemented

mpq_1 = _mpq((1,1))
mpq_0 = _mpq((0,1))

def parse_param(x):
    if isinstance(x, tuple):
        p, q = x
        return [[p, q]], [], []
    if isinstance(x, (int, long)):
        return [[x, 1]], [], []
    x = mpmathify(x)
    if isinstance(x, mpf):
        return [], [x._mpf_], []
    if isinstance(x, mpc):
        return [], [], [x._mpc_]

def _as_num(x):
    if isinstance(x, list):
        return _mpq(x)
    return x

def hypsum(ar, af, ac, br, bf, bc, x):
    prec, rnd = prec_rounding
    if hasattr(x, "_mpf_") and not (ac or bc):
        v = libhyper.hypsum_internal(ar, af, ac, br, bf, bc, x._mpf_, None, prec, rnd)
        return make_mpf(v)
    else:
        if hasattr(x, "_mpc_"):
            re, im = x._mpc_
        else:
            re, im = x._mpf_, libmpf.fzero
        v = libhyper.hypsum_internal(ar, af, ac, br, bf, bc, re, im, prec, rnd)
        return make_mpc(v)

def eval_hyp2f1(a,b,c,z):
    prec, rnd = prec_rounding
    ar, af, ac = parse_param(a)
    br, bf, bc = parse_param(b)
    cr, cf, cc = parse_param(c)
    absz = abs(z)
    if absz == 1:
        # TODO: determine whether it actually does, and otherwise
        # return infinity instead
        print "Warning: 2F1 might not converge for |z| = 1"
    if absz <= 1:
        # All rational
        if ar and br and cr:
            return sum_hyp2f1_rat(ar[0], br[0], cr[0], z)
        return hypsum(ar+br, af+bf, ac+bc, cr, cf, cc, z)
    # Use 1/z transformation
    a = (ar and _as_num(ar[0])) or mpmathify(a)
    b = (br and _as_num(br[0])) or mpmathify(b)
    c = (cr and _as_num(cr[0])) or mpmathify(c)
    orig = mp.prec
    try:
        mp.prec = orig + 15
        h1 = eval_hyp2f1(a, mpq_1-c+a, mpq_1-b+a, 1/z)
        h2 = eval_hyp2f1(b, mpq_1-c+b, mpq_1-a+b, 1/z)
        #s1 = G(c)*G(b-a)/G(b)/G(c-a) * (-z)**(-a) * h1
        #s2 = G(c)*G(a-b)/G(a)/G(c-b) * (-z)**(-b) * h2
        f1 = gammaprod([c,b-a],[b,c-a])
        f2 = gammaprod([c,a-b],[a,c-b])
        s1 = f1 * (-z)**(mpq_0-a) * h1
        s2 = f2 * (-z)**(mpq_0-b) * h2
        v = s1 + s2
    finally:
        mp.prec = orig
    return +v

def sum_hyp0f1_rat(a, z):
    prec, rnd = prec_rounding
    if hasattr(z, "_mpf_"):
        return make_mpf(libhyper.mpf_hyp0f1_rat(a, z._mpf_, prec, rnd))
    else:
        return make_mpc(libhyper.mpc_hyp0f1_rat(a, z._mpc_, prec, rnd))

def sum_hyp1f1_rat(a, b, z):
    prec, rnd = prec_rounding
    if hasattr(z, "_mpf_"):
        return make_mpf(libhyper.mpf_hyp1f1_rat(a, b, z._mpf_, prec, rnd))
    else:
        return make_mpc(libhyper.mpc_hyp1f1_rat(a, b, z._mpc_, prec, rnd))

def sum_hyp2f1_rat(a, b, c, z):
    prec, rnd = prec_rounding
    if hasattr(z, "_mpf_"):
        return make_mpf(libhyper.mpf_hyp2f1_rat(a, b, c, z._mpf_, prec, rnd))
    else:
        return make_mpc(libhyper.mpc_hyp2f1_rat(a, b, c, z._mpc_, prec, rnd))


#---------------------------------------------------------------------------#
#                      And now the user-friendly versions                   #
#---------------------------------------------------------------------------#

def hyper(a_s, b_s, z):
    r"""
    Evaluates the generalized hypergeometric function

    .. math ::

        \,_pF_q(a_1,\ldots,a_p; b_1,\ldots,b_q; z) =
        \sum_{n=0}^\infty \frac{(a_1)_n (a_2)_n \ldots (a_p)_n}
           {(b_1)_n(b_2)_n\ldots(b_q)_n} \frac{z^n}{n!}

    where `(x)_n` denotes the rising factorial (see :func:`rf`).

    The parameters lists ``a_s`` and ``b_s`` may contain integers,
    real numbers, complex numbers, as well as exact fractions given in
    the form of tuples `(p, q)`. :func:`hyper` is optimized to handle
    integers and fractions more efficiently than arbitrary
    floating-point parameters (since rational parameters are by
    far the most common).

    **Examples**

    We can compare the output of :func:`hyper` with :func:`nsum`::

        >>> from sympy.mpmath import *
        >>> mp.dps = 25
        >>> a,b,c,d = 2,3,4,5
        >>> x = 0.25
        >>> print hyper([a,b],[c,d],x)
        1.078903941164934876086237
        >>> fn = lambda n: rf(a,n)*rf(b,n)/rf(c,n)/rf(d,n)*x**n/fac(n)
        >>> print nsum(fn, [0, inf])
        1.078903941164934876086237

    The parameters can be any combination of integers, fractions,
    floats and complex numbers::

        >>> a, b, c, d, e = 1, (-1,2), pi, 3+4j, (2,3)
        >>> x = 0.2j
        >>> print hyper([a,b],[c,d,e],x)
        (0.9923571616434024810831887 - 0.005753848733883879742993122j)
        >>> b, e = -0.5, mpf(2)/3
        >>> fn = lambda n: rf(a,n)*rf(b,n)/rf(c,n)/rf(d,n)/rf(e,n)*x**n/fac(n)
        >>> print nsum(fn, [0, inf])
        (0.9923571616434024810831887 - 0.005753848733883879742993122j)

    """
    p = len(a_s)
    q = len(b_s)
    z = mpmathify(z)
    degree = p, q
    if degree == (0, 1):
        br, bf, bc = parse_param(b_s[0])
        if br:
            return sum_hyp0f1_rat(br[0], z)
        return hypsum([], [], [], br, bf, bc, z)
    if degree == (1, 1):
        ar, af, ac = parse_param(a_s[0])
        br, bf, bc = parse_param(b_s[0])
        if ar and br:
            a, b = ar[0], br[0]
            return sum_hyp1f1_rat(a, b, z)
        return hypsum(ar, af, ac, br, bf, bc, z)
    if degree == (2, 1):
        return eval_hyp2f1(a_s[0], a_s[1], b_s[0], z)
    ars, afs, acs, brs, bfs, bcs = [], [], [], [], [], []
    for a in a_s:
        r, f, c = parse_param(a)
        ars += r
        afs += f
        acs += c
    for b in b_s:
        r, f, c = parse_param(b)
        brs += r
        bfs += f
        bcs += c
    return hypsum(ars, afs, acs, brs, bfs, bcs, z)

def hyp0f1(a, z):
    r"""Hypergeometric function `\,_0F_1`. ``hyp0f1(a,z)`` is equivalent
    to ``hyper([],[a],z)``; see documentation for :func:`hyper` for more
    information."""
    return hyper([], [a], z)

def hyp1f1(a,b,z):
    r"""Hypergeometric function `\,_1F_1`. ``hyp1f1(a,b,z)`` is equivalent
    to ``hyper([a],[b],z)``; see documentation for :func:`hyper` for more
    information."""
    return hyper([a], [b], z)

def hyp2f1(a,b,c,z):
    r"""Hypergeometric function `\,_2F_1`. ``hyp2f1(a,b,c,z)`` is equivalent
    to ``hyper([a,b],[c],z)``; see documentation for :func:`hyper` for more
    information."""
    return hyper([a,b], [c], z)

def _lower_gamma(z, b):
    return hyp1f1(1, 1+z, b) * b**z * exp(-b) / z

def _check_pos(x):
    return isinstance(x, mpf) and x > 0

@funcwrapper
def gammainc(z, a=0, b=inf, regularized=False):
    r"""
    ``gammainc(z, a=0, b=inf)`` computes the (generalized) incomplete
    gamma function with integration limits `[a, b]`:

    .. math ::

      \Gamma(z,a,b) = \int_a^b t^{z-1} e^{-t} \, dt

    The generalized incomplete gamma function reduces to the
    following special cases when one or both endpoints are fixed:

    * `\Gamma(z,0,\infty)` is the standard ("complete")
      gamma function, `\Gamma(z)` (available directly
      as the mpmath function :func:`gamma`)
    * `\Gamma(z,a,\infty)` is the "upper" incomplete gamma
      function, `\Gamma(z,a)`
    * `\Gamma(z,0,b)` is the "lower" incomplete gamma
      function, `\gamma(z,b)`.

    Of course, we have
    `\Gamma(z,0,x) + \Gamma(z,x,\infty) = \Gamma(z)`
    for all `z` and `x`.

    Note however that some authors reverse the order of the
    arguments when defining the lower and upper incomplete
    gamma function, so one should be careful to get the correct
    definition.

    If also given the keyword argument ``regularized=True``,
    :func:`gammainc` computes the "regularized" incomplete gamma
    function

    .. math ::

      P(z,a,b) = \frac{\Gamma(z,a,b)}{\Gamma(z)}.

    **Examples**

    We can compare with numerical quadrature to verify that
    :func:`gammainc` computes the integral in the definition::

        >>> from sympy.mpmath import *
        >>> mp.dps = 20
        >>> print gammainc(2+3j, 4, 10)
        (0.009772126686277051606 - 0.077063730631298989245j)
        >>> print quad(lambda t: t**(2+3j-1) * exp(-t), [4, 10])
        (0.009772126686277051606 - 0.077063730631298989245j)

    The incomplete gamma functions satisfy simple recurrence
    relations::

        >>> mp.dps = 15
        >>> z = 3.5
        >>> a = 2
        >>> print gammainc(z+1, a), z*gammainc(z,a) + a**z*exp(-a)
        10.6013029693353 10.6013029693353
        >>> print gammainc(z+1,0,a), z*gammainc(z,0,a) - a**z*exp(-a)
        1.03042542723211 1.03042542723211

    If `z` is an integer, the recurrence reduces the incomplete gamma
    function to `P(a) \exp(-a) + Q(b) \exp(-b)` where `P` and
    `Q` are polynomials::

        >>> mp.dps = 15
        >>> print gammainc(1, 2), exp(-2)
        0.135335283236613 0.135335283236613
        >>> mp.dps = 50
        >>> identify(gammainc(6, 1, 2), ['exp(-1)', 'exp(-2)'])
        '(326*exp(-1) + (-872)*exp(-2))'

    The incomplete gamma functions reduce to functions such as
    the exponential integral Ei and the error function for special
    arguments::

        >>> mp.dps = 15
        >>> print gammainc(0, 4), -ei(-4)
        0.00377935240984891 0.00377935240984891
        >>> print gammainc(0.5, 0, 2), sqrt(pi)*erf(sqrt(2))
        1.6918067329452 1.6918067329452

    """
    if b == inf:
        if a == 0:
            v = gamma(z)
        else:
            if z == 0:
                # Reduces to exponential integral. Mind branch cuts.
                if _check_pos(a):
                    return -ei(-a)
                else:
                    return -ei(-a) + (log(-a)-log(-1/a))/2-log(a)
            # XXX: avoid poles
            v = gamma(z) - _lower_gamma(z, a)
    elif a == 0:
        v = _lower_gamma(z, b)
    else:
        if z == 0:
            # Reduces to exponential integral
            if _check_pos(a) and _check_pos(b):
                return ei(-b) - ei(-a)
            else:
                return ei(-b)-ei(-a) + \
                    (log(-a)-log(-1/a))/2-log(a) + \
                    (log(-1/b)-log(-b))/2+log(b)
        # XXX: avoid poles
        v = _lower_gamma(z, b) - _lower_gamma(z, a)
    if regularized:
        return v / gamma(z)
    else:
        return v


erf = mpfunc("erf", libhyper.mpf_erf, libhyper.mpc_erf,
    "Error function, erf(z)")

erf.__doc__ = r"""
Computes the error function, `\mathrm{erf}(x)`. The error
function is the normalized antiderivative of the Gaussian function
`\exp(-t^2)`. More precisely,

.. math::

  \mathrm{erf}(x) = \frac{2}{\sqrt \pi} \int_0^x \exp(-t^2) \,dt

**Basic examples**

Simple values and limits include::

    >>> from sympy.mpmath import *
    >>> mp.dps = 15
    >>> print erf(0)
    0.0
    >>> print erf(1)
    0.842700792949715
    >>> print erf(-1)
    -0.842700792949715
    >>> print erf(inf)
    1.0
    >>> print erf(-inf)
    -1.0

For large real `x`, `\mathrm{erf}(x)` approaches 1 very
rapidly::

    >>> print erf(3)
    0.999977909503001
    >>> print erf(5)
    0.999999999998463

The error function is an odd function::

    >>> nprint(chop(taylor(erf, 0, 5)))
    [0.0, 1.12838, 0.0, -0.376126, 0.0, 0.112838]

:func:`erf` implements arbitrary-precision evaluation and
supports complex numbers::

    >>> mp.dps = 50
    >>> print erf(0.5)
    0.52049987781304653768274665389196452873645157575796
    >>> mp.dps = 25
    >>> print erf(1+j)
    (1.316151281697947644880271 + 0.1904534692378346862841089j)

**Related functions**

See also :func:`erfc`, which is more accurate for large `x`,
and :func:`erfi` which gives the antiderivative of
`\exp(t^2)`.

The Fresnel integrals :func:`fresnels` and :func:`fresnelc`
are also related to the error function.
"""


erfc = mpfunc("erfc", libhyper.mpf_erfc, libhyper.mpc_erfc,
    "Complementary error function, erfc(z) = 1-erf(z)")

erfc.__doc__ = r"""
Computes the complementary error function,
`\mathrm{erfc}(x) = 1-\mathrm{erf}(x)`.
This function avoids cancellation that occurs when naively
computing the complementary error function as ``1-erf(x)``::

    >>> from sympy.mpmath import *
    >>> mp.dps = 15
    >>> print 1 - erf(10)
    0.0
    >>> print erfc(10)
    2.08848758376254e-45

:func:`erfc` works accurately even for ludicrously large
arguments::

    >>> print erfc(10**10)
    4.3504398860243e-43429448190325182776

"""

@funcwrapper
def erfi(z):
    r"""
    Computes the imaginary error function, `\mathrm{erfi}(x)`.
    The imaginary error function is defined in analogy with the
    error function, but with a positive sign in the integrand:

    .. math ::

      \mathrm{erfi}(x) = \frac{2}{\sqrt \pi} \int_0^x \exp(t^2) \,dt

    Whereas the error function rapidly converges to 1 as `x` grows,
    the imaginary error function rapidly diverges to infinity.
    The functions are related as
    `\mathrm{erfi}(x) = -i\,\mathrm{erf}(ix)` for all complex
    numbers `x`.

    **Examples**

    Basic values and limits::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> print erfi(0)
        0.0
        >>> print erfi(1)
        1.65042575879754
        >>> print erfi(-1)
        -1.65042575879754
        >>> print erfi(inf)
        +inf
        >>> print erfi(-inf)
        -inf

    Note the symmetry between erf and erfi::

        >>> print erfi(3j)
        (0.0 + 0.999977909503001j)
        >>> print erf(3)
        0.999977909503001
        >>> print erf(1+2j)
        (-0.536643565778565 - 5.04914370344703j)
        >>> print erfi(2+1j)
        (-5.04914370344703 - 0.536643565778565j)

    **Possible issues**

    The current implementation of :func:`erfi` is much less efficient
    and accurate than the one for erf.

    """
    return (2/sqrt(pi)*z) * sum_hyp1f1_rat((1,2),(3,2), z**2)

@funcwrapper
def erfinv(x):
    r"""
    Computes the inverse error function, satisfying

    .. math ::

        \mathrm{erf}(\mathrm{erfinv}(x)) =
        \mathrm{erfinv}(\mathrm{erf}(x)) = x.

    This function is defined only for `-1 \le x \le 1`.

    **Examples**

    Special values include::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> print erfinv(0)
        0.0
        >>> print erfinv(1)
        +inf
        >>> print erfinv(-1)
        -inf

    The domain is limited to the standard interval::

        >>> erfinv(2)
        Traceback (most recent call last):
          ...
        ValueError: erfinv(x) is defined only for -1 <= x <= 1

    It is simple to check that :func:`erfinv` computes inverse values of
    :func:`erf` as promised::

        >>> print erf(erfinv(0.75))
        0.75
        >>> print erf(erfinv(-0.995))
        -0.995

    :func:`erfinv` supports arbitrary-precision evaluation::

        >>> mp.dps = 50
        >>> erf(3)
        mpf('0.99997790950300141455862722387041767962015229291260075')
        >>> erfinv(_)
        mpf('3.0')

    A definite integral involving the inverse error function::

        >>> mp.dps = 15
        >>> print quad(erfinv, [0, 1])
        0.564189583547756
        >>> print 1/sqrt(pi)
        0.564189583547756

    The inverse error function can be used to generate random numbers
    with a Gaussian distribution (although this is a relatively
    inefficient algorithm)::

        >>> nprint([erfinv(2*rand()-1) for n in range(6)]) # doctest: +SKIP
        [-0.586747, 1.10233, -0.376796, 0.926037, -0.708142, -0.732012]

    """
    if x.imag or (x < -1) or (x > 1):
        raise ValueError("erfinv(x) is defined only for -1 <= x <= 1")
    if isnan(x): return x
    if not x: return x
    if x == 1: return inf
    if x == -1: return -inf
    if abs(x) < 0.9:
        a = 0.53728*x**3 + 0.813198*x
    else:
        # An asymptotic formula
        u = log(2/pi/(abs(x)-1)**2)
        a = sign(x) * sqrt(u - log(u))/sqrt(2)
    from optimization import findroot
    return findroot(lambda t: erf(t)-x, a)

@funcwrapper
def npdf(x, mu=0, sigma=1):
    r"""
    ``npdf(x, mu=0, sigma=1)`` evaluates the probability density
    function of a normal distribution with mean value `\mu`
    and variance `\sigma^2`.

    Elementary properties of the probability distribution can
    be verified using numerical integration::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> print quad(npdf, [-inf, inf])
        1.0
        >>> print quad(lambda x: npdf(x, 3), [3, inf])
        0.5
        >>> print quad(lambda x: npdf(x, 3, 2), [3, inf])
        0.5

    See also :func:`ncdf`, which gives the cumulative
    distribution.
    """
    sigma = mpmathify(sigma)
    return exp(-(x-mu)**2/(2*sigma**2)) / (sigma*sqrt(2*pi))

@funcwrapper
def ncdf(x, mu=0, sigma=1):
    r"""
    ``ncdf(x, mu=0, sigma=1)`` evaluates the cumulative distribution
    function of a normal distribution with mean value `\mu`
    and variance `\sigma^2`.

    See also :func:`npdf`, which gives the probability density.

    Elementary properties include::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> print ncdf(pi, mu=pi)
        0.5
        >>> print ncdf(-inf)
        0.0
        >>> print ncdf(+inf)
        1.0

    The cumulative distribution is the integral of the density
    function having identical mu and sigma::

        >>> mp.dps = 15
        >>> print diff(ncdf, 2)
        0.053990966513188
        >>> print npdf(2)
        0.053990966513188
        >>> print diff(lambda x: ncdf(x, 1, 0.5), 0)
        0.107981933026376
        >>> print npdf(0, 1, 0.5)
        0.107981933026376

    """
    a = (x-mu)/(sigma*sqrt(2))
    if a < 0:
        return erfc(-a)/2
    else:
        return (1+erf(a))/2

def ei_as(a):
    extra = 10
    mp.dps += extra
    s = k = p = 1
    while abs(p) > eps:
        p = (p*k)/a
        s += p
        k += 1
    s = (s * exp(a))/a
    mp.dps -= extra
    return s

@funcwrapper
def ei(z):
    r"""
    Computes the exponential integral or Ei-function, `\mathrm{Ei}(x)`.
    The exponential integral is defined as

    .. math ::

      \mathrm{Ei}(x) = \int_{-\infty\,}^x \frac{e^t}{t} \, dt.

    When the integration range includes `t = 0`, the exponential
    integral is interpreted as providing the Cauchy principal value.

    For real `x`, the Ei-function behaves roughly like
    `\mathrm{Ei}(x) \approx \exp(x) + \log(|x|)`.

    This function should not be confused with the family of related
    functions denoted by `E_n` which are also called "exponential
    integrals".

    **Basic examples**

    Some basic values and limits are::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> print ei(0)
        -inf
        >>> print ei(1)
        1.89511781635594
        >>> print ei(inf)
        +inf
        >>> print ei(-inf)
        0.0

    For `x < 0`, the defining integral can be evaluated
    numerically as a reference::

        >>> print ei(-4)
        -0.00377935240984891
        >>> print quad(lambda t: exp(t)/t, [-inf, -4])
        -0.00377935240984891

    :func:`ei` supports complex arguments and arbitrary
    precision evaluation::

        >>> mp.dps = 50
        >>> mp.dps = 50
        >>> print ei(pi)
        10.928374389331410348638445906907535171566338835056
        >>> mp.dps = 25
        >>> print ei(3+4j)
        (-4.154091651642689822535359 + 4.294418620024357476985535j)

    **Related functions**

    The exponential integral is closely related to the logarithmic
    integral. See :func:`li` for additional information.

    The exponential integral is related to the hyperbolic
    and trigonometric integrals (see :func:`chi`, :func:`shi`,
    :func:`ci`, :func:`si`) similarly to how the ordinary
    exponential function is related to the hyperbolic and
    trigonometric functions::

        >>> mp.dps = 15
        >>> print ei(3)
        9.93383257062542
        >>> print chi(3) + shi(3)
        9.93383257062542
        >>> print ci(3j) - j*si(3j) - pi*j/2
        (9.93383257062542 + 0.0j)

    Beware that logarithmic corrections, as in the last example
    above, are required to obtain the correct branch in general.
    For details, see [1].

    The exponential integral is also a special case of the
    hypergeometric function `\,_2F_2`::

        >>> z = 0.6
        >>> print z*hyper([1,1],[2,2],z) + (ln(z)-ln(1/z))/2 + euler
        0.769881289937359
        >>> print ei(z)
        0.769881289937359

    For x large enough use the asymptotic expansion
    ei_as(x) = exp(x)/x * Sum(k!/x^k, (k,0,inf))
    k!/x^k  goes as exp(f(k))
    f(k) = k*log(k/(x*e)) + log(k)/2, with extremal point in
    log(k/x) + 1/(2*k) = 0; therefore the smallest term of the
    asympotic series is k!/x^k ~= e^(-k - 1/2)
    requiring this to be equal to e^-prec one gets x ~= k ~= prec*log(2)
    so that one should use ei_as(x) for x > prec*log(2)

    **References**

    1. Relations between Ei and other functions:
       http://functions.wolfram.com/GammaBetaErf/ExpIntegralEi/27/01/

    2. Abramowitz & Stegun, section 5:
       http://www.math.sfu.ca/~cbm/aands/page_228.htm

    3. Asymptotic expansion for Ei:
       http://mathworld.wolfram.com/En-Function.html
    """
    if z == inf:
        return z
    if z == -inf:
        return -mpf(0)
    if not z:
        return -inf
    if abs(z) > mp.prec * 0.7 + 50:
        r = ei_as(z)
        if z.imag > 0:
            r += j*pi
        elif z.imag < 0:
            r -= j*pi
        return r
    v = z*hypsum([[1,1],[1,1]],[],[],[[2,1],[2,1]],[],[],z) + \
        (log(z)-log(1/z))/2 + euler
    if isinstance(z, mpf) and z < 0:
        return v.real
    return v

@funcwrapper
def li(z):
    r"""
    Computes the logarithmic integral or li-function
    `\mathrm{li}(x)`, defined by

    .. math ::

        \mathrm{li}(x) = \int_0^x \frac{1}{\log t} \, dt

    The logarithmic integral has a singularity at `x = 1`.

    Note that there is a second logarithmic integral, the Li
    function, defined by

    .. math ::

        \mathrm{Li}(x) = \int_2^x \frac{1}{\log t} \, dt

    This "offset logarithmic integral" can be computed via
    :func:`li` using the simple identity
    `\mathrm{Li}(x) = \mathrm{li}(x) - \mathrm{li}(2)`.

    The logarithmic integral should also not be confused with
    the polylogarithm (also denoted by Li), which is implemented
    as :func:`polylog`.

    **Examples**

    Some basic values and limits::

        >>> from sympy.mpmath import *
        >>> mp.dps = 30
        >>> print li(0)
        0.0
        >>> print li(1)
        -inf
        >>> print li(1)
        -inf
        >>> print li(2)
        1.04516378011749278484458888919
        >>> print findroot(li, 2)
        1.45136923488338105028396848589
        >>> print li(inf)
        +inf

    The logarithmic integral can be evaluated for arbitrary
    complex arguments::

        >>> mp.dps = 20
        >>> print li(3+4j)
        (3.1343755504645775265 + 2.6769247817778742392j)

    The logarithmic integral is related to the exponential integral::

        >>> print ei(log(3))
        2.1635885946671919729
        >>> print li(3)
        2.1635885946671919729

    The logarithmic integral grows like `O(x/\log(x))`::

        >>> mp.dps = 15
        >>> x = 10**100
        >>> print x/log(x)
        4.34294481903252e+97
        >>> print li(x)
        4.3619719871407e+97

    The prime number theorem states that the number of primes less
    than `x` is asymptotic to `\mathrm{li}(x)`. For example,
    it is known that there are exactly 1,925,320,391,606,803,968,923
    prime numbers less than `10^{23}` [1]. The logarithmic integral
    provides a very accurate estimate::

        >>> print li(2) + li(10**23)
        1.92532039161405e+21

    A definite integral is::

        >>> print quad(li, [0, 1])
        -0.693147180559945
        >>> print -ln(2)
        -0.693147180559945

    **References**

    1. http://mathworld.wolfram.com/PrimeCountingFunction.html

    2. http://mathworld.wolfram.com/LogarithmicIntegral.html

    """
    if not z:
        return z
    if z == 1:
        return -inf
    return ei(log(z))

ci = mpfunc('ci', libhyper.mpf_ci, libhyper.mpc_ci, '')
si = mpfunc('si', libhyper.mpf_si, libhyper.mpc_si, '')

ci.__doc__ = r"""
Computes the cosine integral,

.. math ::

    \mathrm{Ci}(x) = -\int_x^{\infty} \frac{\cos t}{t}\,dt
    = \gamma + \log x + \int_0^x \frac{\cos t - 1}{t}\,dt

**Examples**

Some values and limits::

    >>> from sympy.mpmath import *
    >>> mp.dps = 25
    >>> print ci(0)
    -inf
    >>> print ci(1)
    0.3374039229009681346626462
    >>> print ci(pi)
    0.07366791204642548599010096
    >>> print ci(inf)
    0.0
    >>> print ci(-inf)
    (0.0 + 3.141592653589793238462643j)
    >>> print ci(2+3j)
    (1.408292501520849518759125 - 2.983617742029605093121118j)

The cosine integral behaves roughly like the sinc function
(see :func:`sinc`) for large real `x`::

    >>> print ci(10**10)
    -4.875060251748226537857298e-11
    >>> print sinc(10**10)
    -4.875060250875106915277943e-11
    >>> print chop(limit(ci, inf))
    0.0

It has infinitely many roots on the positive real axis::

    >>> print findroot(ci, 1)
    0.6165054856207162337971104
    >>> print findroot(ci, 2)
    3.384180422551186426397851

We can evaluate the defining integral as a reference::

    >>> mp.dps = 15
    >>> print -quadosc(lambda t: cos(t)/t, [5, inf], omega=1)
    -0.190029749656644
    >>> print ci(5)
    -0.190029749656644

Some infinite series can be evaluated using the
cosine integral::

    >>> print nsum(lambda k: (-1)**k/(fac(2*k)*(2*k)), [1,inf])
    -0.239811742000565
    >>> print ci(1) - euler
    -0.239811742000565

"""

si.__doc__ = r"""
Computes the sine integral,

.. math ::

    \mathrm{Si}(x) = \int_0^x \frac{\sin t}{t}\,dt.

The sine integral is thus the antiderivative of the sinc
function (see :func:`sinc`).

**Examples**

Some values and limits::

    >>> from sympy.mpmath import *
    >>> mp.dps = 25
    >>> print si(0)
    0.0
    >>> print si(1)
    0.9460830703671830149413533
    >>> print si(-1)
    -0.9460830703671830149413533
    >>> print si(pi)
    1.851937051982466170361053
    >>> print si(inf)
    1.570796326794896619231322
    >>> print si(-inf)
    -1.570796326794896619231322
    >>> print si(2+3j)
    (4.547513889562289219853204 + 1.399196580646054789459839j)

The sine integral approaches `\pi/2` for large real `x`::

    >>> print si(10**10)
    1.570796326707584656968511
    >>> print pi/2
    1.570796326794896619231322

We can evaluate the defining integral as a reference::

    >>> mp.dps = 15
    >>> print quad(sinc, [0, 5])
    1.54993124494467
    >>> print si(5)
    1.54993124494467

Some infinite series can be evaluated using the
sine integral::

    >>> print nsum(lambda k: (-1)**k/(fac(2*k+1)*(2*k+1)), [0,inf])
    0.946083070367183
    >>> print si(1)
    0.946083070367183

"""

@funcwrapper
def chi(z):
    r"""
    Computes the hyperbolic cosine integral, defined
    in analogy with the cosine integral (see :func:`ci`) as

    .. math ::

        \mathrm{Chi}(x) = -\int_x^{\infty} \frac{\cosh t}{t}\,dt
        = \gamma + \log x + \int_0^x \frac{\cosh t - 1}{t}\,dt

    Some values and limits::

        >>> from sympy.mpmath import *
        >>> mp.dps = 25
        >>> print chi(0)
        -inf
        >>> print chi(1)
        0.8378669409802082408946786
        >>> print chi(inf)
        +inf
        >>> print findroot(chi, 0.5)
        0.5238225713898644064509583
        >>> print chi(2+3j)
        (-0.1683628683277204662429321 + 2.625115880451325002151688j)

    """
    if not z:
        return -inf
    z2 = (z/2)**2
    return euler + log(z) + \
        z2*hypsum([[1,1],[1,1]],[],[],[[2,1],[2,1],[3,2]],[],[],z2)

@funcwrapper
def shi(z):
    r"""
    Computes the hyperbolic sine integral, defined
    in analogy with the sine integral (see :func:`si`) as

    .. math ::

        \mathrm{Shi}(x) = \int_0^x \frac{\sinh t}{t}\,dt.

    Some values and limits::

        >>> from sympy.mpmath import *
        >>> mp.dps = 25
        >>> print shi(0)
        0.0
        >>> print shi(1)
        1.057250875375728514571842
        >>> print shi(-1)
        -1.057250875375728514571842
        >>> print shi(inf)
        +inf
        >>> print shi(2+3j)
        (-0.1931890762719198291678095 + 2.645432555362369624818525j)

    """
    z2 = (z/2)**2
    return z*hypsum([[1,2]],[],[],[[3,2],[3,2]],[],[],z2)

@funcwrapper
def fresnels(z):
    r"""
    Computes the Fresnel sine integral

    .. math ::

        S(x) = \int_0^x \sin\left(\frac{\pi t^2}{2}\right) \,dt

    Note that some sources define this function
    without the normalization factor `\pi/2`.

    **Examples**

    Some basic values and limits::

        >>> from sympy.mpmath import *
        >>> mp.dps = 25
        >>> print fresnels(0)
        0.0
        >>> print fresnels(inf)
        0.5
        >>> print fresnels(-inf)
        -0.5
        >>> print fresnels(1)
        0.4382591473903547660767567
        >>> print fresnels(1+2j)
        (36.72546488399143842838788 + 15.58775110440458732748279j)

    Comparing with the definition::

        >>> print fresnels(3)
        0.4963129989673750360976123
        >>> print quad(lambda t: sin(pi*t**2/2), [0,3])
        0.4963129989673750360976123

    """
    if z == inf:
        return mpf(0.5)
    if z == -inf:
        return mpf(-0.5)
    return pi*z**3/6*hypsum([[3,4]],[],[],[[3,2],[7,4]],[],[],-pi**2*z**4/16)

@funcwrapper
def fresnelc(z):
    r"""
    Computes the Fresnel cosine integral

    .. math ::

        C(x) = \int_0^x \cos\left(\frac{\pi t^2}{2}\right) \,dt

    Note that some sources define this function
    without the normalization factor `\pi/2`.

    **Examples**

    Some basic values and limits::

        >>> from sympy.mpmath import *
        >>> mp.dps = 25
        >>> print fresnelc(0)
        0.0
        >>> print fresnelc(inf)
        0.5
        >>> print fresnelc(-inf)
        -0.5
        >>> print fresnelc(1)
        0.7798934003768228294742064
        >>> print fresnelc(1+2j)
        (16.08787137412548041729489 - 36.22568799288165021578758j)

    Comparing with the definition::

        >>> print fresnelc(3)
        0.6057207892976856295561611
        >>> print quad(lambda t: cos(pi*t**2/2), [0,3])
        0.6057207892976856295561611

    """
    if z == inf:
        return mpf(0.5)
    if z == -inf:
        return mpf(-0.5)
    return z*hypsum([[1,4]],[],[],[[1,2],[5,4]],[],[],-pi**2*z**4/16)

@funcwrapper
def airyai(z):
    r"""
    Computes the Airy function `\mathrm{Ai}(x)`, which is
    a solution of the Airy differential equation `y''-xy=0`.
    The Ai-function behaves roughly like a slowly decaying
    sine wave for `x < 0` and like a decreasing exponential for
    `x > 0`.

    Limits and values include::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> print airyai(0), 1/(3**(2/3.)*gamma(2/3.))
        0.355028053887817 0.355028053887817
        >>> print airyai(1)
        0.135292416312881
        >>> print airyai(-1)
        0.535560883292352
        >>> print airyai(inf)
        0.0
        >>> print airyai(-inf)
        0.0

    :func:`airyai` uses a series expansion around `x = 0`,
    so it is slow for extremely large arguments. Here are
    some evaluations for moderately large arguments::

        >>> print airyai(-100)
        0.176753393239553
        >>> print airyai(100)
        2.63448215208818e-291
        >>> print airyai(50+50j)
        (-5.31790195707456e-68 - 1.16358800377071e-67j)
        >>> print airyai(-50+50j)
        (1.04124253736317e+158 + 3.3475255449236e+157j)

    The first negative root is::

        >>> print findroot(airyai, -2)
        -2.33810741045977

    We can verify the differential equation::

        >>> for x in [-3.4, 0, 2.5, 1+2j]:
        ...     print abs(diff(airyai, x, 2) - x*airyai(x)) < eps
        ...
        True
        True
        True
        True

    The Taylor series expansion around `x = 0` starts with
    the following coefficients (note that every third term
    is zero)::

        >>> nprint(chop(taylor(airyai, 0, 5)))
        [0.355028, -0.258819, 0.0, 5.91713e-2, -2.15683e-2, 0.0]

    The Airy functions are a special case of Bessel functions.
    For `x < 0`, we have::

        >>> x = 3
        >>> print airyai(-x)
        -0.378814293677658
        >>> p = 2*(x**1.5)/3
        >>> print sqrt(x)*(besselj(1/3.,p) + besselj(-1/3.,p))/3
        -0.378814293677658

    """
    if z == inf or z == -inf:
        return 1/z
    if z.real > 2:
        # cancellation: both terms are ~ 2^(z^1.5),
        # result is ~ 2^(-z^1.5), so need ~2*z^1.5 extra bits
        mp.prec += 2*int(z.real**1.5)
    z3 = z**3 / 9
    a = sum_hyp0f1_rat((2,3), z3) / (cbrt(9) * gamma(mpf(2)/3))
    b = z * sum_hyp0f1_rat((4,3), z3) / (cbrt(3) * gamma(mpf(1)/3))
    return a - b

@funcwrapper
def airybi(z):
    r"""
    Computes the Airy function `\mathrm{Bi}(x)`, which is
    a solution of the Airy differential equation `y''-xy=0`.
    The Bi-function behaves roughly like a slowly decaying
    sine wave for `x < 0` and like an increasing exponential
    for `x > 0`.

    Limits and values include::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> print airybi(0), 1/(3**(1/6.)*gamma(2/3.))
        0.614926627446001 0.614926627446001
        >>> print airybi(1)
        1.20742359495287
        >>> print airybi(-1)
        0.103997389496945
        >>> print airybi(inf)
        +inf
        >>> print airybi(-inf)
        0.0

    :func:`airyai` uses a series expansion around `x = 0`,
    so it is slow for extremely large arguments. Here are
    some evaluations for moderately large arguments::

        >>> print airybi(-100)
        0.0242738876801601
        >>> print airybi(100)
        6.0412239966702e+288
        >>> print airybi(50+50j)
        (-5.32207626732144e+63 + 1.47845029116524e+65j)
        >>> print airybi(-50+50j)
        (-3.3475255449236e+157 + 1.04124253736317e+158j)

    The first negative root is::

        >>> print findroot(airybi, -1)
        -1.17371322270913

    We can verify the differential equation::

        >>> for x in [-3.4, 0, 2.5, 1+2j]:
        ...     print abs(diff(airybi, x, 2) - x*airybi(x)) < eps
        ...
        True
        True
        True
        True

    The Taylor series expansion around `x = 0` starts with
    the following coefficients (note that every third term
    is zero)::

        >>> nprint(chop(taylor(airybi, 0, 5)))
        [0.614927, 0.448288, 0.0, 0.102488, 3.73574e-2, 0.0]

    The Airy functions are a special case of Bessel functions.
    For `x < 0`, we have::

        >>> x = 3
        >>> print airybi(-x)
        -0.198289626374927
        >>> p = 2*(x**1.5)/3
        >>> print sqrt(x/3)*(besselj(-1/3.,p) - besselj(1/3.,p))
        -0.198289626374926
    """
    if z == inf:
        return z
    if z == -inf:
        return 1/z
    z3 = z**3 / 9
    rt = nthroot(3, 6)
    a = sum_hyp0f1_rat((2,3), z3) / (rt * gamma(mpf(2)/3))
    b = z * rt * sum_hyp0f1_rat((4,3), z3) / gamma(mpf(1)/3)
    return a + b

ellipk = mpfunc('ellipk', libhyper.mpf_ellipk, libhyper.mpc_ellipk, '')
ellipe = mpfunc('ellipe', libhyper.mpf_ellipe, libhyper.mpc_ellipe, '')

ellipk.__doc__ = \
r"""
Evaluates the complete elliptic integral of the first kind,
`K(m)`, defined by

.. math ::

    K(m) = \int_0^{\pi/2} \frac{1}{\sqrt{1-m \sin^2 t}} dt.

Note that the argument is the parameter `m = k^2`,
not the modulus `k` which is sometimes used.

Alternatively, in terms of a hypergeometric function,
we have:

.. math ::

    K(m) = \frac{\pi}{2} \,_2F_1(1/2, 1/2, 1, m)

**Examples**

Values and limits include::

    >>> from sympy.mpmath import *
    >>> mp.dps = 25
    >>> print ellipk(0)
    1.570796326794896619231322
    >>> print ellipk(inf)
    (0.0 + 0.0j)
    >>> print ellipk(-inf)
    0.0
    >>> print ellipk(1)
    +inf
    >>> print ellipk(-1)
    1.31102877714605990523242
    >>> print ellipk(2)
    (1.31102877714605990523242 - 1.31102877714605990523242j)

Verifying the defining integral and hypergeometric
representation::

    >>> print ellipk(0.5)
    1.85407467730137191843385
    >>> print quad(lambda t: (1-0.5*sin(t)**2)**-0.5, [0, pi/2])
    1.85407467730137191843385
    >>> print pi/2*hyp2f1(0.5,0.5,1,0.5)
    1.85407467730137191843385

Evaluation is supported for arbitrary complex `m`::

    >>> print ellipk(3+4j)
    (0.9111955638049650086562171 + 0.6313342832413452438845091j)

A definite integral::

    >>> print quad(ellipk, [0, 1])
    2.0
"""

ellipe.__doc__ = \
r"""
Evaluates the complete elliptic integral of the second kind,
`E(m)`, defined by

.. math ::

    E(m) = \int_0^{\pi/2} \sqrt{1-m \sin^2 t} dt.

Note that the argument is the parameter `m = k^2`,
not the modulus `k` which is sometimes used.

Alternatively, in terms of a hypergeometric function,
we have:

.. math ::

    E(m) = \frac{\pi}{2} \,_2F_1(1/2, -1/2, 1, m)

**Examples**

Basic values and limits::

    >>> from sympy.mpmath import *
    >>> mp.dps = 25
    >>> print ellipe(0)
    1.570796326794896619231322
    >>> print ellipe(1)
    1.0
    >>> print ellipe(-1)
    1.910098894513856008952381
    >>> print ellipe(2)
    (0.5990701173677961037199612 + 0.5990701173677961037199612j)
    >>> print ellipe(inf)
    (0.0 + +infj)
    >>> print ellipe(-inf)
    +inf

Verifying the defining integral and hypergeometric
representation::

    >>> print ellipe(0.5)
    1.350643881047675502520175
    >>> print quad(lambda t: sqrt(1-0.5*sin(t)**2), [0, pi/2])
    1.350643881047675502520175
    >>> print pi/2*hyp2f1(0.5,-0.5,1,0.5)
    1.350643881047675502520175

Evaluation is supported for arbitrary complex `m`::

    >>> print ellipe(0.5+0.25j)
    (1.360868682163129682716687 - 0.1238733442561786843557315j)
    >>> print ellipe(3+4j)
    (1.499553520933346954333612 - 1.577879007912758274533309j)

A definite integral::

    >>> print quad(ellipe, [0,1])
    1.333333333333333333333333

"""

def agm(a, b=1):
    r"""
    ``agm(a, b)`` computes the arithmetic-geometric mean of `a` and
    `b`, defined as the limit of the following iteration:

    .. math ::

        a_0 = a

        b_0 = b

        a_{n+1} = \frac{a_n+b_n}{2}

        b_{n+1} = \sqrt{a_n b_n}

    This function can be called with a single argument, computing
    `\mathrm{agm}(a,1) = \mathrm{agm}(1,a)`.

    **Examples**

    It is a well-known theorem that the geometric mean of
    two distinct positive numbers is less than the arithmetic
    mean. It follows that the arithmetic-geometric mean lies
    between the two means::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> a = mpf(3)
        >>> b = mpf(4)
        >>> print sqrt(a*b)
        3.46410161513775
        >>> print agm(a,b)
        3.48202767635957
        >>> print (a+b)/2
        3.5

    The arithmetic-geometric mean is scale-invariant::

        >>> print agm(10*e, 10*pi)
        29.261085515723
        >>> print 10*agm(e, pi)
        29.261085515723

    As an order-of-magnitude estimate, `\mathrm{agm}(1,x) \approx x`
    for large `x`::

        >>> print agm(10**10)
        643448704.760133
        >>> print agm(10**50)
        1.34814309345871e+48

    For tiny `x`, `\mathrm{agm}(1,x) \approx -\pi/(2 \log(x/4))`::

        >>> print agm('0.01')
        0.262166887202249
        >>> print -pi/2/log('0.0025')
        0.262172347753122

    The arithmetic-geometric mean can also be computed for complex
    numbers::

        >>> print agm(3, 2+j)
        (2.51055133276184 + 0.547394054060638j)

    The AGM iteration converges very quickly (each step doubles
    the number of correct digits), so :func:`agm` supports efficient
    high-precision evaluation::

        >>> mp.dps = 10000
        >>> a = agm(1,2)
        >>> str(a)[-10:]
        '1679581912'

    **Mathematical relations**

    The arithmetic-geometric mean may be used to evaluate the
    following two parametric definite integrals:

    .. math ::

      I_1 = \int_0^{\infty}
        \frac{1}{\sqrt{(x^2+a^2)(x^2+b^2)}} \,dx

      I_2 = \int_0^{\pi/2}
        \frac{1}{\sqrt{a^2 \cos^2(x) + b^2 \sin^2(x)}} \,dx

    We have::

        >>> mp.dps = 15
        >>> a = 3
        >>> b = 4
        >>> f1 = lambda x: ((x**2+a**2)*(x**2+b**2))**-0.5
        >>> f2 = lambda x: ((a*cos(x))**2 + (b*sin(x))**2)**-0.5
        >>> print quad(f1, [0, inf])
        0.451115405388492
        >>> print quad(f2, [0, pi/2])
        0.451115405388492
        >>> print pi/(2*agm(a,b))
        0.451115405388492

    A formula for `\Gamma(1/4)`::

        >>> print gamma(0.25)
        3.62560990822191
        >>> print sqrt(2*sqrt(2*pi**3)/agm(1,sqrt(2)))
        3.62560990822191

    **Possible issues**

    The branch cut chosen for complex `a` and `b` is somewhat
    arbitrary.

    """
    if b == 1:
        return agm1(a)
    a = mpmathify(a)
    b = mpmathify(b)
    prec, rounding = prec_rounding
    if isinstance(a, mpf) and isinstance(b, mpf):
        try:
            v = libhyper.mpf_agm(a._mpf_, b._mpf_, prec, rounding)
            return make_mpf(v)
        except ComplexResult:
            pass
    if isinstance(a, mpf): a = (a._mpf_, libmpf.fzero)
    else: a = a._mpc_
    if isinstance(b, mpf): b = (b._mpf_, libmpf.fzero)
    else: b = b._mpc_
    return make_mpc(libhyper.mpc_agm(a, b, prec, rounding))

agm1 = mpfunc('agm1', libhyper.mpf_agm1, libhyper.mpc_agm1,
    'Fast alias for agm(1,a) = agm(a,1)')

@funcwrapper
def jacobi(n, a, b, x):
    r"""
    ``jacobi(n, a, b, x)`` evaluates the Jacobi polynomial
    `P_n^{(a,b)}(x)`. The Jacobi polynomials are a special
    case of the hypergeometric function `\,_2F_1` given by:

    .. math ::

        P_n^{(a,b)}(x) = {n+a \choose n}
          \,_2F_1\left(-n,1+a+b+n,a+1,\frac{1-x}{2}\right).

    Note that this definition generalizes to nonintegral values
    of `n`. When `n` is an integer, the hypergeometric series
    terminates after a finite number of terms, giving
    a polynomial in `x`.

    **Evaluation of Jacobi polynomials**

    A special evaluation is `P_n^{(a,b)}(1) = {n+a \choose n}`::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> print jacobi(4, 0.5, 0.25, 1)
        2.4609375
        >>> print binomial(4+0.5, 4)
        2.4609375

    A Jacobi polynomial of degree `n` is equal to its
    Taylor polynomial of degree `n`. The explicit
    coefficients of Jacobi polynomials can therefore
    be recovered easily using :func:`taylor`::

        >>> for n in range(5):
        ...     nprint(taylor(lambda x: jacobi(n,1,2,x), 0, n))
        ...
        [1.0]
        [-0.5, 2.5]
        [-0.75, -1.5, 5.25]
        [0.5, -3.5, -3.5, 10.5]
        [0.625, 2.5, -11.25, -7.5, 20.625]

    For nonintegral `n`, the Jacobi "polynomial" is no longer
    a polynomial::

        >>> nprint(taylor(lambda x: jacobi(0.5,1,2,x), 0, 4))
        [0.309983, 1.84119, -1.26933, 1.26699, -1.34808]

    **Orthogonality**

    The Jacobi polynomials are orthogonal on the interval
    `[-1, 1]` with respect to the weight function
    `w(x) = (1-x)^a (1+x)^b`. That is,
    `w(x) P_n^{(a,b)}(x) P_m^{(a,b)}(x)` integrates to
    zero if `m \ne n` and to a nonzero number if `m = n`.

    The orthogonality is easy to verify using numerical
    quadrature::

        >>> P = jacobi
        >>> f = lambda x: (1-x)**a * (1+x)**b * P(m,a,b,x) * P(n,a,b,x)
        >>> a = 2
        >>> b = 3
        >>> m, n = 3, 4
        >>> print chop(quad(f, [-1, 1]), 1)
        0.0
        >>> m, n = 4, 4
        >>> print quad(f, [-1, 1])
        1.9047619047619

    **Differential equation**

    The Jacobi polynomials are solutions of the differential
    equation

    .. math ::

      (1-x^2) y'' + (b-a-(a+b+2)x) y' + n (n+a+b+1) y = 0.

    We can verify that :func:`jacobi` approximately satisfies
    this equation::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> a = 2.5
        >>> b = 4
        >>> n = 3
        >>> y = lambda x: jacobi(n,a,b,x)
        >>> x = pi
        >>> A0 = n*(n+a+b+1)*y(x)
        >>> A1 = (b-a-(a+b+2)*x)*diff(y,x)
        >>> A2 = (1-x**2)*diff(y,x,2)
        >>> nprint(A2 + A1 + A0, 1)
        4.0e-12

    The difference of order `10^{-12}` is as close to zero as
    it could be at 15-digit working precision, since the terms
    are large::

        >>> print A0, A1, A2
        26560.2328981879 -21503.7641037294 -5056.46879445852

    """
    return binomial(n+a,n) * hyp2f1(-n,1+n+a+b,a+1,(1-x)/2)

@funcwrapper
def legendre(n, x):
    r"""
    ``legendre(n, x)`` evaluates the Legendre polynomial `P_n(x)`.
    The Legendre polynomials are given by the formula

    .. math ::

        P_n(x) = \frac{1}{2^n n!} \frac{d^n}{dx^n} (x^2 -1)^n.

    Alternatively, they can be computed recursively using

    .. math ::

        P_0(x) = 1

        P_1(x) = x

        (n+1) P_{n+1}(x) = (2n+1) x P_n(x) - n P_{n-1}(x).

    A third definition is in terms of the hypergeometric function
    `\,_2F_1`, whereby they can be generalized to arbitrary `n`:

    .. math ::

        P_n(x) = \,_2F_1\left(-n, n+1, 1, \frac{1-x}{2}\right)

    **Basic evaluation**

    The Legendre polynomials assume fixed values at the points
    `x = -1` and `x = 1`::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> nprint([legendre(n, 1) for n in range(6)])
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        >>> nprint([legendre(n, -1) for n in range(6)])
        [1.0, -1.0, 1.0, -1.0, 1.0, -1.0]

    The coefficients of Legendre polynomials can be recovered
    using degree-`n` Taylor expansion::

        >>> for n in range(5):
        ...     nprint(chop(taylor(lambda x: legendre(n, x), 0, n)))
        ...
        [1.0]
        [0.0, 1.0]
        [-0.5, 0.0, 1.5]
        [0.0, -1.5, 0.0, 2.5]
        [0.375, 0.0, -3.75, 0.0, 4.375]

    The roots of Legendre polynomials are located symmetrically
    on the interval `[-1, 1]`::

        >>> for n in range(5):
        ...     nprint(polyroots(taylor(lambda x: legendre(n, x), 0, n)[::-1]))
        ...
        []
        [0.0]
        [-0.57735, 0.57735]
        [-0.774597, 0.0, 0.774597]
        [-0.861136, -0.339981, 0.339981, 0.861136]

    An example of an evaluation for arbitrary `n`::

        >>> print legendre(0.75, 2+4j)
        (1.94952805264875 + 2.1071073099422j)

    **Orthogonality**

    The Legendre polynomials are orthogonal on `[-1, 1]` with respect
    to the trivial weight `w(x) = 1`. That is, `P_m(x) P_n(x)`
    integrates to zero if `m \ne n` and to `2/(2n+1)` if `m = n`::

        >>> m, n = 3, 4
        >>> print quad(lambda x: legendre(m,x)*legendre(n,x), [-1, 1])
        0.0
        >>> m, n = 4, 4
        >>> print quad(lambda x: legendre(m,x)*legendre(n,x), [-1, 1])
        0.222222222222222

    **Differential equation**

    The Legendre polynomials satisfy the differential equation

    .. math ::

        ((1-x^2) y')' + n(n+1) y' = 0.

    We can verify this numerically::

        >>> n = 3.6
        >>> x = 0.73
        >>> P = legendre
        >>> A = diff(lambda t: (1-t**2)*diff(lambda u: P(n,u), t), x)
        >>> B = n*(n+1)*P(n,x)
        >>> nprint(A+B,1)
        9.0e-16

    """
    if isint(n):
        n = int(n)
    if x == -1:
        # TODO: hyp2f1 should handle this
        if isint(n):
            return (-1)**(n + (n>=0)) * mpf(-1)
        if not int(floor(re(n))) % 2:
            return -inf
        return inf
    return hyp2f1(-n,n+1,1,(1-x)/2)

@funcwrapper
def chebyt(n, x):
    r"""
    ``chebyt(n, x)`` evaluates the Chebyshev polynomial of the first
    kind `T_n(x)`, defined by the identity

    .. math ::

        T_n(\cos x) = \cos(n x).

    The Chebyshev polynomials of the first kind are a special
    case of the Jacobi polynomials, and by extension of the
    hypergeometric function `\,_2F_1`. They can thus also be
    evaluated for nonintegral `n`.

    **Basic evaluation**

    The coefficients of the `n`-th polynomial can be recovered
    using using degree-`n` Taylor expansion::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> for n in range(5):
        ...     nprint(chop(taylor(lambda x: chebyt(n, x), 0, n)))
        ...
        [1.0]
        [0.0, 1.0]
        [-1.0, 0.0, 2.0]
        [0.0, -3.0, 0.0, 4.0]
        [1.0, 0.0, -8.0, 0.0, 8.0]

    **Orthogonality**

    The Chebyshev polynomials of the first kind are orthogonal
    on the interval `[-1, 1]` with respect to the weight
    function `w(x) = 1/\sqrt{1-x^2}`::

        >>> f = lambda x: chebyt(m,x)*chebyt(n,x)/sqrt(1-x**2)
        >>> m, n = 3, 4
        >>> nprint(quad(f, [-1, 1]),1)
        0.0
        >>> m, n = 4, 4
        >>> print quad(f, [-1, 1])
        1.57079632596448

    """
    return hyp2f1(-n,n,0.5,(1-x)/2)

@funcwrapper
def chebyu(n, x):
    r"""
    ``chebyu(n, x)`` evaluates the Chebyshev polynomial of the second
    kind `U_n(x)`, defined by the identity

    .. math ::

        U_n(\cos x) = \frac{\sin((n+1)x)}{\sin(x)}.

    The Chebyshev polynomials of the second kind are a special
    case of the Jacobi polynomials, and by extension of the
    hypergeometric function `\,_2F_1`. They can thus also be
    evaluated for nonintegral `n`.

    **Basic evaluation**

    The coefficients of the `n`-th polynomial can be recovered
    using using degree-`n` Taylor expansion::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> for n in range(5):
        ...     nprint(chop(taylor(lambda x: chebyu(n, x), 0, n)))
        ...
        [1.0]
        [0.0, 2.0]
        [-1.0, 0.0, 4.0]
        [0.0, -4.0, 0.0, 8.0]
        [1.0, 0.0, -12.0, 0.0, 16.0]

    **Orthogonality**

    The Chebyshev polynomials of the second kind are orthogonal
    on the interval `[-1, 1]` with respect to the weight
    function `w(x) = \sqrt{1-x^2}`::

        >>> f = lambda x: chebyu(m,x)*chebyu(n,x)*sqrt(1-x**2)
        >>> m, n = 3, 4
        >>> print quad(f, [-1, 1])
        0.0
        >>> m, n = 4, 4
        >>> print quad(f, [-1, 1])
        1.5707963267949

    """
    return (n+1) * hyp2f1(-n, n+2, 1.5, (1-x)/2)

@funcwrapper
def besselj(v, x):
    r"""
    ``besselj(n,x)`` computes the Bessel function of the first kind
    `J_n(x)`. Bessel functions of the first kind are defined as
    solutions of the differential equation

    .. math ::

        x^2 y'' + x y' + (x^2 - n^2) y = 0

    which appears, among other things, when solving the radial
    part of Laplace's equation in cylindrical coordinates. This
    equation has two solutions for given `n`, where the
    `J_n`-function is the solution that is nonsingular at `x = 0`.
    For positive integer `n`, `J_n(x)` behaves roughly like a sine
    (odd `n`) or cosine (even `n`) multiplied by a magnitude factor
    that decays slowly as `x \to \pm\infty`.

    Generally, `J_n` is a special case of the hypergeometric
    function `\,_0F_1`:

    .. math ::

        J_n(x) = \frac{x^n}{2^n \Gamma(n+1)}
                 \,_0F_1\left(n+1,-\frac{x^2}{4}\right)

    **Examples**

    Evaluation is supported for arbitrary arguments, and at
    arbitrary precision::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> print besselj(2, 1000)
        -0.024777229528606
        >>> print besselj(4, 0.75)
        0.000801070086542314
        >>> print besselj(2, 1000j)
        (-2.48071721019185e+432 + 0.0j)
        >>> mp.dps = 25
        >>> print besselj(0.75j, 3+4j)
        (-2.778118364828153309919653 - 1.5863603889018621585533j)
        >>> mp.dps = 50
        >>> print besselj(1, pi)
        0.28461534317975275734531059968613140570981118184947

    The Bessel functions of the first kind satisfy simple
    symmetries around `x = 0`::

        >>> mp.dps = 15
        >>> nprint([besselj(n,0) for n in range(5)])
        [1.0, 0.0, 0.0, 0.0, 0.0]
        >>> nprint([besselj(n,pi) for n in range(5)])
        [-0.304242, 0.284615, 0.485434, 0.333458, 0.151425]
        >>> nprint([besselj(n,-pi) for n in range(5)])
        [-0.304242, -0.284615, 0.485434, -0.333458, 0.151425]

    Roots of Bessel functions are often used::

        >>> nprint([findroot(j0, k) for k in [2, 5, 8, 11, 14]])
        [2.40483, 5.52008, 8.65373, 11.7915, 14.9309]
        >>> nprint([findroot(j1, k) for k in [3, 7, 10, 13, 16]])
        [3.83171, 7.01559, 10.1735, 13.3237, 16.4706]

    The roots are not periodic, but the distance between successive
    roots asymptotically approaches `2 \pi`. Bessel functions of
    the first kind have the following normalization::

        >>> print quadosc(j0, [0, inf], period=2*pi)
        1.0
        >>> print quadosc(j1, [0, inf], period=2*pi)
        1.0

    For `n = 1/2` or `n = -1/2`, the Bessel function reduces to a
    trigonometric function::

        >>> x = 10
        >>> print besselj(0.5, x), sqrt(2/(pi*x))*sin(x)
        -0.13726373575505 -0.13726373575505
        >>> print besselj(-0.5, x), sqrt(2/(pi*x))*cos(x)
        -0.211708866331398 -0.211708866331398

    """
    if isint(v):
        v = int(v)
        if isinstance(x, mpf):
            return make_mpf(libhyper.mpf_besseljn(v, x._mpf_, mp.prec))
        if isinstance(x, mpc):
            return make_mpc(libhyper.mpc_besseljn(v, x._mpc_, mp.prec))
    hx = x/2
    return hx**v * hyp0f1(v+1, -hx**2) / factorial(v)

def j0(x):
    """Computes the Bessel function `J_0(x)`. See :func:`besselj`."""
    return besselj(0, x)

def j1(x):
    """Computes the Bessel function `J_1(x)`.  See :func:`besselj`."""
    return besselj(1, x)

@funcwrapper
def bessely(n,x):
    r"""
    ``bessely(n,x)`` computes the Bessel function of the second kind,

    .. math ::

        Y_n(x) = \frac{J_n(x) \cos(\pi n) - J_{-n}(x)}{\sin(\pi n)}.

    For `n` an integer, this formula should be understood as a
    limit.

    **Examples**

    Some values of `Y_n(x)`::

        >>> from sympy.mpmath import *
        >>> mp.dps = 25
        >>> print bessely(0,0), bessely(1,0), bessely(2,0)
        -inf -inf -inf
        >>> print bessely(1, pi)
        0.3588729167767189594679827
        >>> print bessely(0.5, 3+4j)
        (9.242861436961450520325216 - 3.085042824915332562522402j)

    """
    intdist = abs(n.imag) + abs(n.real-floor(n.real+0.5))
    if not intdist:
        h = +eps
        mp.prec *= 2
        n += h
    else:
        mp.prec += -int(log(intdist, 2)+1)
    return (besselj(n,x)*cospi(n) - besselj(-n,x))/sinpi(n)

@funcwrapper
def besseli(n,x):
    r"""
    ``besseli(n,x)`` computes the modified Bessel function of the first
    kind,

    .. math ::

        I_n(x) = i^{-n} J_n(ix)

    **Examples**

    Some values of `I_n(x)`::

        >>> from sympy.mpmath import *
        >>> mp.dps = 25
        >>> print besseli(0,0)
        1.0
        >>> print besseli(1,0)
        0.0
        >>> print besseli(0,1)
        1.266065877752008335598245
        >>> print besseli(3.5, 2+3j)
        (-0.2904369752642538144289025 - 0.4469098397654815837307006j)

    For integers `n`, the following integral representation holds::

        >>> mp.dps = 15
        >>> n = 3
        >>> x = 2.3
        >>> print quad(lambda t: exp(x*cos(t))*cos(n*t), [0,pi])/pi
        0.349223221159309
        >>> print besseli(n,x)
        0.349223221159309

    """
    if isint(n):
        n = abs(int(n))
    hx = x/2
    return hx**n * hyp0f1(n+1, hx**2) / factorial(n)

@funcwrapper
def besselk(n,x):
    r"""
    ``besseli(n,x)`` computes the modified Bessel function of the
    second kind,

    .. math ::

        K_n(x) = \frac{\pi}{2} \frac{I_{-n}(x)-I_{n}(x)}{\sin(\pi n)}

    For `n` an integer, this formula should be understood as a
    limit.

    **Examples**

    Some values and limits of `K_n(x)`::

        >>> from sympy.mpmath import *
        >>> mp.dps = 25
        >>> print besselk(0,0)
        +inf
        >>> print besselk(1,0)
        +inf
        >>> print besselk(0,1)
        0.4210244382407083333356274
        >>> print besselk(3.5, 2+3j)
        (-0.02090732889633760668464128 + 0.2464022641351420167819697j)

    """
    intdist = abs(n.imag) + abs(n.real-floor(n.real+0.5))
    if not intdist:
        h = +eps
        mp.prec *= 2
        n += h
    else:
        mp.prec += -int(log(intdist, 2)+1)
    return pi*(besseli(-n,x)-besseli(n,x))/(2*sinpi(n))

def hankel1(n,x):
    r"""
    ``hankel1(n,x)`` computes the Hankel function of the first kind,
    which is the complex combination of Bessel functions given by

    .. math ::

        H_n^{(1)}(x) = J_n(x) + i Y_n(x).

    **Examples**

    The Hankel function is generally complex-valued::

        >>> from sympy.mpmath import *
        >>> mp.dps = 25
        >>> print hankel1(2, pi)
        (0.4854339326315091097054957 - 0.0999007139290278787734903j)
        >>> print hankel1(3.5, pi)
        (0.2340002029630507922628888 - 0.6419643823412927142424049j)

    """
    return besselj(n,x) + j*bessely(n,x)

def hankel2(n,x):
    r"""
    ``hankel2(n,x)`` computes the Hankel function of the second kind,
    which is the complex combination of Bessel functions given by

    .. math ::

        H_n^{(2)}(x) = J_n(x) - i Y_n(x).

    **Examples**

    The Hankel function is generally complex-valued::

        >>> from sympy.mpmath import *
        >>> mp.dps = 25
        >>> print hankel2(2, pi)
        (0.4854339326315091097054957 + 0.0999007139290278787734903j)
        >>> print hankel2(3.5, pi)
        (0.2340002029630507922628888 + 0.6419643823412927142424049j)

    """
    return besselj(n,x) - j*bessely(n,x)

@funcwrapper
def lambertw(z, k=0, approx=None):
    r"""
    The Lambert W function `W(z)` is defined as the inverse function
    of `w \exp(w)`. In other words, the value of `W(z)` is such that
    `z = W(z) \exp(W(z))` for any complex number `z`.

    The Lambert W function is a multivalued function with infinitely
    many branches. Each branch gives a separate solution of the
    equation `w \exp(w)`. All branches are supported by
    :func:`lambertw`:

    * ``lambertw(z)`` gives the principal solution (branch 0)

    * ``lambertw(z, k)`` gives the solution on branch `k`

    The Lambert W function has two partially real branches: the
    principal branch (`k = 0`) is real for real `z > -1/e`, and the
    `k = -1` branch is real for `-1/e < z < 0`. All branches except
    `k = 0` have a logarithmic singularity at `z = 0`.

    **Basic examples**

    The Lambert W function is the inverse of `w \exp(w)`::

        >>> from sympy.mpmath import *
        >>> mp.dps = 35
        >>> w = lambertw(1)
        >>> print w
        0.56714329040978387299996866221035555
        >>> print w*exp(w)
        1.0

    Any branch gives a valid inverse::

        >>> w = lambertw(1, k=3)
        >>> print w    # doctest: +NORMALIZE_WHITESPACE
        (-2.8535817554090378072068187234910812 +
          17.113535539412145912607826671159289j)
        >>> print w*exp(w)
        (1.0 + 3.5075477124212226194278700785075126e-36j)

    **Applications to equation-solving**

    The Lambert W function may be used to solve various kinds of
    equations, such as finding the value of the infinite power
    tower `z^{z^{z^{\ldots}}}`::

        >>> def tower(z, n):
        ...     if n == 0:
        ...         return z
        ...     return z ** tower(z, n-1)
        ...
        >>> tower(0.5, 100)
        0.641185744504986
        >>> mp.dps = 50
        >>> print -lambertw(-log(0.5))/log(0.5)
        0.6411857445049859844862004821148236665628209571911

    **Properties**

    The Lambert W function grows roughly like the natural logarithm
    for large arguments::

        >>> mp.dps = 15
        >>> print lambertw(1000)
        5.2496028524016
        >>> print log(1000)
        6.90775527898214
        >>> print lambertw(10**100)
        224.843106445119
        >>> print log(10**100)
        230.258509299405

    The principal branch of the Lambert W function has a rational
    Taylor series expansion around `z = 0`::

        >>> nprint(taylor(lambertw, 0, 6), 10)
        [0.0, 1.0, -1.0, 1.5, -2.666666667, 5.208333333, -10.8]

    Some special values and limits are::

        >>> mp.dps = 15
        >>> print lambertw(0)
        0.0
        >>> print lambertw(1)
        0.567143290409784
        >>> print lambertw(e)
        1.0
        >>> print lambertw(inf)
        +inf
        >>> print lambertw(0, k=-1)
        -inf
        >>> print lambertw(0, k=3)
        -inf
        >>> print lambertw(inf, k=3)
        (+inf + 18.8495559215388j)

    The `k = 0` and `k = -1` branches join at `z = -1/e` where
    `W(z) = -1` for both branches. Since `-1/e` can only be represented
    approximately with mpmath numbers, evaluating the Lambert W function
    at this point only gives `-1` approximately::

        >>> mp.dps = 25
        >>> print lambertw(-1/e, 0)
        -0.999999999999837133022867
        >>> print lambertw(-1/e, -1)
        -1.00000000000016286697718

    If `-1/e` happens to round in the negative direction, there might be
    a small imaginary part::

        >>> mp.dps = 15
        >>> print lambertw(-1/e)
        (-1.0 + 8.22007971511612e-9j)

    **Possible issues**

    The evaluation can become inaccurate very close to the branch point
    at `-1/e`. In some corner cases, :func:`lambertw` might currently
    fail to converge, or can end up on the wrong branch.

    **Algorithm**

    Halley's iteration is used to invert `w \exp(w)`, using a first-order
    asymptotic approximation (`O(\log(w))` or `O(w)`) as the initial
    estimate.

    The definition, implementation and choice of branches is based
    on Corless et al, "On the Lambert W function", Adv. Comp. Math. 5
    (1996) 329-359, available online here:
    http://www.apmaths.uwo.ca/~djeffrey/Offprints/W-adv-cm.pdf

    TODO: use a series expansion when extremely close to the branch point
    at `-1/e` and make sure that the proper branch is chosen there
    """
    if isnan(z):
        return z
    mp.prec += 20
    # We must be extremely careful near the singularities at -1/e and 0
    u = exp(-1)
    if abs(z) <= u:
        if not z:
            # w(0,0) = 0; for all other branches we hit the pole
            if not k:
                return z
            return -inf
        if not k:
            w = z
        # For small real z < 0, the -1 branch behaves roughly like log(-z)
        elif k == -1 and not z.imag and z.real < 0:
            w = log(-z)
        # Use a simple asymptotic approximation.
        else:
            w = log(z)
            # The branches are roughly logarithmic. This approximation
            # gets better for large |k|; need to check that this always
            # works for k ~= -1, 0, 1.
            if k: w += k * 2*pi*j
    elif k == 0 and z.imag and abs(z) <= 0.6:
        w = z
    else:
        if z == inf:
            if k == 0:
                return z
            else:
                return z + 2*k*pi*j
        if z == -inf:
            return (-z) + (2*k+1)*pi*j
        # Simple asymptotic approximation as above
        w = log(z)
        if k: w += k * 2*pi*j
    # Use Halley iteration to solve w*exp(w) = z
    two = mpf(2)
    weps = ldexp(eps, 15)
    for i in xrange(100):
        ew = exp(w)
        wew = w*ew
        wewz = wew-z
        wn = w - wewz/(wew+ew-(w+two)*wewz/(two*w+two))
        if abs(wn-w) < weps*abs(wn):
            return wn
        else:
            w = wn
    print "Warning: Lambert W iteration failed to converge:", z
    return wn

@funcwrapper
def barnesg(z):
    r"""
    Evaluates the Barnes G-function, which generalizes the
    superfactorial (:func:`superfac`) and by extension also the
    hyperfactorial (:func:`hyperfac`) to the complex numbers
    in an analogous way to how the gamma function generalizes
    the ordinary factorial.

    The Barnes G-function may be defined in terms of a Weierstrass
    product:

    .. math ::

        G(z+1) = (2\pi)^{z/2} e^{-[z(z+1)+\gamma z^2]/2}
        \prod_{n=1}^\infty
        \left[\left(1+\frac{z}{n}\right)^ne^{-z+z^2/(2n)}\right]

    For positive integers `n`, we have have relation to superfactorials
    `G(n) = \mathrm{sf}(n-2) = 0! \cdot 1! \cdots (n-2)!`.

    **Examples**

    Some elementary values and limits of the Barnes G-function::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> print barnesg(1), barnesg(2), barnesg(3)
        1.0 1.0 1.0
        >>> print barnesg(4)
        2.0
        >>> print barnesg(5)
        12.0
        >>> print barnesg(6)
        288.0
        >>> print barnesg(7)
        34560.0
        >>> print barnesg(8)
        24883200.0
        >>> print barnesg(inf)
        +inf
        >>> print barnesg(0), barnesg(-1), barnesg(-2)
        0.0 0.0 0.0

    Closed-form values are known for some rational arguments::

        >>> print barnesg('1/2')
        0.603244281209446
        >>> print sqrt(exp(0.25+log(2)/12)/sqrt(pi)/glaisher**3)
        0.603244281209446
        >>> print barnesg('1/4')
        0.29375596533861
        >>> print nthroot(exp('3/8')/exp(catalan/pi)/
        ...      gamma(0.25)**3/sqrt(glaisher)**9, 4)
        0.29375596533861

    The Barnes G-function satisfies the functional equation
    `G(z+1) = \Gamma(z) G(z)`::

        >>> z = pi
        >>> print barnesg(z+1)
        2.39292119327948
        >>> print gamma(z)*barnesg(z)
        2.39292119327948

    The asymptotic growth rate of the Barnes G-function is related to
    the Glaisher-Kinkelin constant::

        >>> print limit(lambda n: barnesg(n+1)/(n**(n**2/2-mpf(1)/12)*
        ...     (2*pi)**(n/2)*exp(-3*n**2/4)), inf)
        0.847536694177301
        >>> print exp('1/12')/glaisher
        0.847536694177301

    The Barnes G-function can be differentiated in closed form::

        >>> z = 3
        >>> print diff(barnesg, z)
        0.264507203401607
        >>> print barnesg(z)*((z-1)*psi(0,z)-z+(log(2*pi)+1)/2)
        0.264507203401607

    Evaluation is supported for arbitrary arguments and at arbitrary
    precision::

        >>> print barnesg(6.5)
        2548.7457695685
        >>> print barnesg(-pi)
        0.00535976768353037
        >>> print barnesg(3+4j)
        (-0.000676375932234244 - 4.42236140124728e-5j)
        >>> mp.dps = 50
        >>> print barnesg(1/sqrt(2))
        0.81305501090451340843586085064413533788206204124732
        >>> q = barnesg(10j)
        >>> print q.real
        0.000000000021852360840356557241543036724799812371995850552234
        >>> print q.imag
        -0.00000000000070035335320062304849020654215545839053210041457588

    **References**

    1. Whittaker & Watson, *A Course of Modern Analysis*,
       Cambridge University Press, 4th edition (1927), p.264
    2. http://en.wikipedia.org/wiki/Barnes_G-function
    3. http://mathworld.wolfram.com/BarnesG-Function.html

    """
    if isinf(z):
        if z == inf:
            return z
        return nan
    if isnan(z):
        return z
    if (not z.imag) and z.real <= 0 and isint(z.real):
        return z*0
    # Account for size (would not be needed if computing log(G))
    if abs(z) > 5:
        mp.dps += 2*log(abs(z),2)
    # Estimate terms for asymptotic expansion
    N = mp.dps // 2 + 5
    G = 1
    while re(z) < N:
        G /= gamma(z)
        z += 1
    z -= 1
    s = mpf(1)/12
    s -= log(glaisher)
    s += z*log(2*pi)/2
    s += (z**2/2-mpf(1)/12)*log(z)
    s -= 3*z**2/4
    z2k = z2 = z**2
    for k in xrange(1, N+1):
        t = bernoulli(2*k+2) / (4*k*(k+1)*z2k)
        if abs(t) < eps:
            #print k, N      # check how many terms were needed
            break
        z2k *= z2
        s += t
    #if k == N:
    #    print "warning: series for barnesg failed to converge"
    return G*exp(s)

def superfac(z):
    r"""
    Computes the superfactorial, defined as the product of
    consecutive factorials

    .. math ::

        \mathrm{sf}(n) = \prod_{k=1}^n k!

    For general complex `z`, `\mathrm{sf}(z)` is defined
    in terms of the Barnes G-function (see :func:`barnesg`).

    **Examples**

    The first few superfactorials are (OEIS A000178)::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> for n in range(10):
        ...     print n, superfac(n)
        ...
        0 1.0
        1 1.0
        2 2.0
        3 12.0
        4 288.0
        5 34560.0
        6 24883200.0
        7 125411328000.0
        8 5.05658474496e+15
        9 1.83493347225108e+21

    Superfactorials grow very rapidly::

        >>> print superfac(1000)
        3.24570818422368e+1177245
        >>> print superfac(10**10)
        2.61398543581249e+467427913956904067453

    Evaluation is supported for arbitrary arguments::

        >>> mp.dps = 25
        >>> print superfac(pi)
        17.20051550121297985285333
        >>> print superfac(2+3j)
        (-0.005915485633199789627466468 + 0.008156449464604044948738263j)
        >>> print diff(superfac, 1)
        0.2645072034016070205673056

    **References**

    1. http://www.research.att.com/~njas/sequences/A000178

    """
    return barnesg(z+2)

@funcwrapper
def hyperfac(z):
    r"""
    Computes the hyperfactorial, defined for integers as the product

    .. math ::

        H(n) = \prod_{k=1}^n k^k.


    The hyperfactorial satisfies the recurrence formula `H(z) = z^z H(z-1)`.
    It can be defined more generally in terms of the Barnes G-function (see
    :func:`barnesg`) and the gamma function by the formula

    .. math ::

        H(z) = \frac{\Gamma(z+1)^z}{G(z)}.

    The extension to complex numbers can also be done via
    the integral representation

    .. math ::

        H(z) = (2\pi)^{-z/2} \exp \left[
            {z+1 \choose 2} + \int_0^z \log(t!)\,dt
            \right].

    **Examples**

    The rapidly-growing sequence of hyperfactorials begins
    (OEIS A002109)::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> for n in range(10):
        ...     print n, hyperfac(n)
        ...
        0 1.0
        1 1.0
        2 4.0
        3 108.0
        4 27648.0
        5 86400000.0
        6 4031078400000.0
        7 3.3197663987712e+18
        8 5.56964379417266e+25
        9 2.15779412229419e+34

    Some even larger hyperfactorials are::

        >>> print hyperfac(1000)
        5.46458120882585e+1392926
        >>> print hyperfac(10**10)
        4.60408207642219e+489142638002418704309

    The hyperfactorial can be evaluated for arbitrary arguments::

        >>> print hyperfac(0.5)
        0.880449235173423
        >>> print diff(hyperfac, 1)
        0.581061466795327
        >>> print hyperfac(pi)
        205.211134637462
        >>> print hyperfac(-10+1j)
        (3.01144471378225e+46 - 2.45285242480185e+46j)

    The recurrence property of the hyperfactorial holds
    generally::

        >>> z = 3-4*j
        >>> print hyperfac(z)
        (-4.49795891462086e-7 - 6.33262283196162e-7j)
        >>> print z**z * hyperfac(z-1)
        (-4.49795891462086e-7 - 6.33262283196162e-7j)
        >>> z = mpf(-0.6)
        >>> print chop(z**z * hyperfac(z-1))
        1.28170142849352
        >>> print hyperfac(z)
        1.28170142849352

    The hyperfactorial may also be computed using the integral
    definition::

        >>> z = 2.5
        >>> print hyperfac(z)
        15.9842119922237
        >>> print (2*pi)**(-z/2)*exp(binomial(z+1,2) +
        ...     quad(lambda t: loggamma(t+1), [0, z]))
        15.9842119922237

    :func:`hyperfac` supports arbitrary-precision evaluation::

        >>> mp.dps = 50
        >>> print hyperfac(10)
        215779412229418562091680268288000000000000000.0
        >>> print hyperfac(1/sqrt(2))
        0.89404818005227001975423476035729076375705084390942

    **References**

    1. http://www.research.att.com/~njas/sequences/A002109
    2. http://mathworld.wolfram.com/Hyperfactorial.html

    """
    # XXX: estimate needed extra bits accurately
    if z == inf:
        return z
    if abs(z) > 5:
        extra = 4*int(log(abs(z),2))
    else:
        extra = 0
    mp.prec += extra
    if not z.imag and z.real < 0 and isint(z.real):
        n = int(re(z))
        h = hyperfac(-n-1)
        if ((n+1)//2) & 1:
            h = -h
        if isinstance(z, mpc):
            return h + 0j
        return h
    zp1 = z+1
    # Wrong branch cut
    #v = gamma(zp1)**z
    #mp.prec -= extra
    #return v / barnesg(zp1)
    v = exp(z*loggamma(zp1))
    mp.prec -= extra
    return v / barnesg(zp1)

@funcwrapper
def loggamma(z):
    r"""
    Computes the log-gamma function. Unlike `\ln(\Gamma(z))`, which
    has infinitely many complex branch cuts, the log-gamma function
    only has a single branch cut along the negative half-axis.
    The functions are identical only on (and very close to) the positive
    half-axis; elsewhere they differ by `2 n \pi i` (the real parts
    agree)::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> print loggamma(13.2), log(gamma(13.2))
        20.494004194566 20.494004194566
        >>> print loggamma(3+4j)
        (-1.75662678460378 + 4.74266443803466j)
        >>> print log(gamma(3+4j))
        (-1.75662678460378 - 1.54052086914493j)

    Note: this is a placeholder implementation. It is slower than
    :func:`gamma`, and is in particular *not* faster than :func:`gamma`
    for large arguments.
    """
    a = z.real
    b = z.imag
    if not b and a > 0:
        return log(gamma(z))
    u = arg(z)
    w = log(gamma(z))
    if b:
        gi = -b - u/2 + a*u + b*log(abs(z))
        n = floor((gi-w.imag)/(2*pi)+0.5) * (2*pi)
        return w + n*j
    elif a < 0:
        n = int(floor(a))
        w += (n-(n%2))*pi*j
    return w

@funcwrapper
def siegeltheta(t):
    r"""
    Computes the Riemann-Siegel theta function,

    .. math ::

        \theta(t) = \frac{
        \log\Gamma\left(\frac{1+2it}{4}\right) -
        \log\Gamma\left(\frac{1-2it}{4}\right)
        }{2i} - \frac{\log \pi}{2} t.

    The Riemann-Siegel theta function is important in
    providing the phase factor for the Z-function
    (see :func:`siegelz`). Evaluation is supported for real and
    complex arguments::

        >>> from sympy.mpmath import *
        >>> mp.dps = 25
        >>> print siegeltheta(0)
        0.0
        >>> print siegeltheta(inf)
        +inf
        >>> print siegeltheta(-inf)
        -inf
        >>> print siegeltheta(1)
        -1.767547952812290388302216
        >>> print siegeltheta(10+0.25j)
        (-3.068638039426838572528867 + 0.05804937947429712998395177j)

    The Riemann-Siegel theta function has odd symmetry around `t = 0`,
    two local extreme points and three real roots including 0 (located
    symmetrically)::

        >>> nprint(chop(taylor(siegeltheta, 0, 5)))
        [0.0, -2.68609, 0.0, 2.69433, 0.0, -6.40218]
        >>> print findroot(diffun(siegeltheta), 7)
        6.28983598883690277966509
        >>> print findroot(siegeltheta, 20)
        17.84559954041086081682634

    For large `t`, there is a famous asymptotic formula
    for `\theta(t)`, to first order given by::

        >>> t = mpf(10**6)
        >>> print siegeltheta(t)
        5488816.353078403444882823
        >>> print -t*log(2*pi/t)/2-t/2
        5488816.745777464310273645
    """
    if t.imag:
        # XXX: cancellation occurs
        a = loggamma(0.25+0.5j*t)
        b = loggamma(0.25-0.5j*t)
        return -log(pi)/2*t - 0.5j*(a-b)
    else:
        if isinf(t):
            return t
        return loggamma(0.25+0.5j*t).imag - log(pi)/2*t

@funcwrapper
def grampoint(n):
    r"""
    Gives the `n`-th Gram point `g_n`, defined as the solution
    to the equation `\theta(g_n) = \pi n` where `\theta(t)`
    is the Riemann-Siegel theta function (:func:`siegeltheta`).

    The first few Gram points are::

        >>> from sympy.mpmath import *
        >>> mp.dps = 25
        >>> print grampoint(0)
        17.84559954041086081682634
        >>> print grampoint(1)
        23.17028270124630927899664
        >>> print grampoint(2)
        27.67018221781633796093849
        >>> print grampoint(3)
        31.71797995476405317955149

    Checking the definition::

        >>> print siegeltheta(grampoint(3))
        9.42477796076937971538793
        >>> print 3*pi
        9.42477796076937971538793

    A large Gram point::

        >>> print grampoint(10**10)
        3293531632.728335454561153

    Gram points are useful when studying the Z-function
    (:func:`siegelz`). See the documentation of that function
    for additional examples.

    :func:`grampoint` can solve the defining equation for
    nonintegral `n`. There is a fixed point where `g(x) = x`::

        >>> print findroot(lambda x: grampoint(x) - x, 10000)
        9146.698193171459265866198

    **References**

    1. http://mathworld.wolfram.com/GramPoint.html

    """
    # Asymptotic expansion, from
    # http://mathworld.wolfram.com/GramPoint.html
    g = 2*pi*exp(1+lambertw((8*n+1)/(8*e)))
    from optimization import findroot
    return findroot(lambda t: siegeltheta(t)-pi*n, g)

@funcwrapper
def siegelz(t):
    r"""
    Computes the Z-function, also known as the Riemann-Siegel Z function,

    .. math ::

        Z(t) = e^{i \theta(t)} \zeta(1/2+it)

    where `\zeta(s)` is the Riemann zeta function (:func:`zeta`)
    and where `\theta(t)` denotes the Riemann-Siegel theta function
    (see :func:`siegeltheta`).

    Evaluation is supported for real and complex arguments::

        >>> from sympy.mpmath import *
        >>> mp.dps = 25
        >>> print siegelz(1)
        -0.7363054628673177346778998
        >>> print siegelz(3+4j)
        (-0.1852895764366314976003936 - 0.2773099198055652246992479j)

    The Z-function has a Maclaurin expansion::

        >>> nprint(chop(taylor(siegelz, 0, 4)))
        [-1.46035, 0.0, 2.73588, 0.0, -8.39357]

    The Z-function `Z(t)` is equal to `\pm |\zeta(s)|` on the
    critical strip `s = 1/2+it` (i.e. for real arguments `t`
    to `Z`).  Its zeros coincide with those of the Riemann zeta
    function::

        >>> print findroot(siegelz, 14)
        14.13472514173469379045725
        >>> print findroot(siegelz, 20)
        21.02203963877155499262848
        >>> print findroot(zeta, 0.5+14j)
        (0.5 + 14.13472514173469379045725j)
        >>> print findroot(zeta, 0.5+20j)
        (0.5 + 21.02203963877155499262848j)

    Since the Z-function is real-valued on the critical strip
    (and unlike `|\zeta(s)|` analytic), it is useful for
    investigating the zeros of the Riemann zeta function.
    For example, one can use a root-finding algorithm based
    on sign changes::

        >>> print findroot(siegelz, [100, 200], solver='bisect')
        176.4414342977104188888926

    To locate roots, Gram points `g_n` which can be computed
    by :func:`grampoint` are useful. If `(-1)^n Z(g_n)` is
    positive for two consecutive `n`, then `Z(t)` must have
    a zero between those points::

        >>> g10 = grampoint(10)
        >>> g11 = grampoint(11)
        >>> (-1)**10 * siegelz(g10) > 0
        True
        >>> (-1)**11 * siegelz(g11) > 0
        True
        >>> print findroot(siegelz, [g10, g11], solver='bisect')
        56.44624769706339480436776
        >>> print g10, g11
        54.67523744685325626632663 57.54516517954725443703014

    """
    v = exp(j*siegeltheta(t))*zeta(0.5+j*t)
    if isinstance(t, mpf):
        return v.real
    return v

@funcwrapper
def bernpoly(n, z):
    """
    Evaluates the Bernoulli polynomial `B_n(z)`.

    The first few Bernoulli polynomials are::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> for n in range(6):
        ...     nprint(chop(taylor(lambda x: bernpoly(n,x), 0, n)))
        ...
        [1.0]
        [-0.5, 1.0]
        [0.166667, -1.0, 1.0]
        [0.0, 0.5, -1.5, 1.0]
        [-3.33333e-2, 0.0, 1.0, -2.0, 1.0]
        [0.0, -0.166667, 0.0, 1.66667, -2.5, 1.0]

    At `z = 0`, the Bernoulli polynomial evaluates to a
    Bernoulli number (see :func:`bernoulli`)::

        >>> print bernpoly(12, 0), bernoulli(12)
        -0.253113553113553 -0.253113553113553
        >>> print bernpoly(13, 0), bernoulli(13)
        0.0 0.0

    """
    n = int(n)
    assert n >= 0
    # XXX: optimize
    return sum(binomial(n,k)*bernoulli(k)*z**(n-k) for k in xrange(0,n+1))

# TODO: this should be implemented low-level
def polylog_series(s, z):
    tol = +eps
    l = mpf(0)
    k = 1
    zk = z
    while 1:
        term = zk / k**s
        l += term
        if abs(term) < tol:
            break
        zk *= z
        k += 1
    return l

def polylog_continuation(n, z):
    if n < 0:
        return z*0
    a = -(2*pi*j)**n/fac(n) * bernpoly(n, log(z)/(2*pi*j))
    if isinstance(z, mpf) and z < 0:
        a = a.real
    if z.imag < 0 or (z.imag == 0 and z.real >= 1):
        a -= 2*pi*j*log(z)**(n-1)/fac(n-1)
    return a

def polylog_unitcircle(n, z):
    tol = +eps
    if n > 1:
        l = mpf(0)
        logz = log(z)
        logmz = mpf(1)
        m = 0
        while 1:
            if (n-m) != 1:
                term = zeta(n-m) * logmz / fac(m)
                if term and abs(term) < tol:
                    break
                l += term
            logmz *= logz
            m += 1
        l += log(z)**(n-1)/fac(n-1)*(harmonic(n-1)-log(-log(z)))
    elif n < 1:  # else
        l = fac(-n)*(-log(z))**(n-1)
        logz = log(z)
        logkz = mpf(1)
        k = 0
        while 1:
            b = bernoulli(k-n+1)
            if b:
                term = b*logkz/(fac(k)*(k-n+1))
                if abs(term) < tol:
                    break
                l -= term
            logkz *= logz
            k += 1
    else:
        raise ValueError
    if isinstance(z, mpf) and z < 0:
        l = l.real
    return l

@funcwrapper
def polylog(s, z):
    r"""
    Computes the polylogarithm, defined by the sum

    .. math ::

        \mathrm{Li}_s(z) = \sum_{k=1}^{\infty} \frac{z^k}{k^s}.

    This series is convergent only for `|z| < 1`, so elsewhere
    the analytic continuation is implied.

    The polylogarithm should not be confused with the logarithmic
    integral (also denoted by Li or li), which is implemented
    as :func:`li`.

    **Examples**

    The polylogarithm satisfies a huge number of functional identities.
    A sample of polylogarithm evaluations is shown below::

        >>> from sympy.mpmath import *
        >>> mp.dps = 15
        >>> print polylog(1,0.5), log(2)
        0.693147180559945 0.693147180559945
        >>> print polylog(2,0.5), (pi**2-6*log(2)**2)/12
        0.582240526465012 0.582240526465012
        >>> print polylog(2,-phi), -log(phi)**2-pi**2/10
        -1.21852526068613 -1.21852526068613
        >>> print polylog(3,0.5), 7*zeta(3)/8-pi**2*log(2)/12+log(2)**3/6
        0.53721319360804 0.53721319360804

    :func:`polylog` can evaluate the analytic continuation of the
    polylogarithm when `s` is an integer::

        >>> print polylog(2, 10)
        (0.536301287357863 - 7.23378441241546j)
        >>> print polylog(2, -10)
        -4.1982778868581
        >>> print polylog(2, 10j)
        (-3.05968879432873 + 3.71678149306807j)
        >>> print polylog(-2, 10)
        -0.150891632373114
        >>> print polylog(-2, -10)
        0.067618332081142
        >>> print polylog(-2, 10j)
        (0.0384353698579347 + 0.0912451798066779j)

    Some more examples, with arguments on the unit circle (note that
    the series definition cannot be used for computation here)::

        >>> print polylog(2,j)
        (-0.205616758356028 + 0.915965594177219j)
        >>> print j*catalan-pi**2/48
        (-0.205616758356028 + 0.915965594177219j)
        >>> print polylog(3,exp(2*pi*j/3))
        (-0.534247512515375 + 0.765587078525922j)
        >>> print -4*zeta(3)/9 + 2*j*pi**3/81
        (-0.534247512515375 + 0.765587078525921j)

    Polylogarithms of different order are related by integration
    and differentiation::

        >>> s, z = 3, 0.5
        >>> print polylog(s+1, z)
        0.517479061673899
        >>> print quad(lambda t: polylog(s,t)/t, [0, z])
        0.517479061673899
        >>> print z*diff(lambda t: polylog(s+2,t), z)
        0.517479061673899

    Taylor series expansions around `z = 0` are::

        >>> for n in range(-3, 4):
        ...     nprint(taylor(lambda x: polylog(n,x), 0, 5))
        ...
        [0.0, 1.0, 8.0, 27.0, 64.0, 125.0]
        [0.0, 1.0, 4.0, 9.0, 16.0, 25.0]
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
        [0.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        [0.0, 1.0, 0.5, 0.333333, 0.25, 0.2]
        [0.0, 1.0, 0.25, 0.111111, 6.25e-2, 4.0e-2]
        [0.0, 1.0, 0.125, 3.7037e-2, 1.5625e-2, 8.0e-3]

    The series defining the polylogarithm is simultaneously
    a Taylor series and an L-series. For certain values of `z`, the
    polylogarithm reduces to a pure zeta function::

        >>> print polylog(pi, 1), zeta(pi)
        1.17624173838258 1.17624173838258
        >>> print polylog(pi, -1), -altzeta(pi)
        -0.909670702980385 -0.909670702980385

    Evaluation for arbitrary, nonintegral `s` is supported
    for `z` within the unit circle:

        >>> print polylog(3+4j, 0.25)
        (0.24258605789446 - 0.00222938275488344j)
        >>> print nsum(lambda k: 0.25**k / k**(3+4j), [1,inf])
        (0.24258605789446 - 0.00222938275488344j)

    **References**

    1. Richard Crandall, "Note on fast polylogarithm computation"
       http://people.reed.edu/~crandall/papers/Polylog.pdf
    2. http://en.wikipedia.org/wiki/Polylogarithm
    3. http://mathworld.wolfram.com/Polylogarithm.html

    """
    if z == 1:
        return zeta(s)
    if z == -1:
        return -altzeta(s)
    if s == 0:
        return z/(1-z)
    if s == 1:
        return -log(1-z)
    if s == -1:
        return z/(1-z)**2
    if abs(z) <= 0.75 or (not isint(s) and abs(z) < 0.99):
        return polylog_series(s, z)
    if abs(z) >= 1.4 and isint(s):
        return (-1)**(s+1)*polylog_series(s, 1/z) + polylog_continuation(s, z)
    if isint(s):
        return polylog_unitcircle(int(s), z)
    raise NotImplementedError("polylog for arbitrary s and z")
    # This could perhaps be used in some cases
    #from quadrature import quad
    #return quad(lambda t: t**(s-1)/(exp(t)/z-1),[0,inf])/gamma(s)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
