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
import libintmath

from mptypes import (\
    MultiPrecisionArithmetic,
    def_mp_builtin,
    defun_wrapped,
    defun,
    defun_static,
    mp,
    constant,
    ComplexResult,
)

def_mp_builtin = def_mp_builtin
libelefun = libelefun
libmpc = libmpc
libhyper = libhyper
gammazeta = gammazeta

# The following multiprecision functions are implemented entirely with
# low-level code
sqrt = def_mp_builtin('sqrt', libelefun.mpf_sqrt, libmpc.mpc_sqrt, libmpi.mpi_sqrt, "principal square root")
cbrt = def_mp_builtin('cbrt', libelefun.mpf_cbrt, libmpc.mpc_cbrt, None, "principal cubic root")
exp = def_mp_builtin('exp', libelefun.mpf_exp, libmpc.mpc_exp, libmpi.mpi_exp, "exponential function")
ln = def_mp_builtin('ln', libelefun.mpf_log, libmpc.mpc_log, libmpi.mpi_log, "natural logarithm")
cos = def_mp_builtin('cos', libelefun.mpf_cos, libmpc.mpc_cos, libmpi.mpi_cos, "cosine")
sin = def_mp_builtin('sin', libelefun.mpf_sin, libmpc.mpc_sin, libmpi.mpi_sin, "sine")
tan = def_mp_builtin('tan', libelefun.mpf_tan, libmpc.mpc_tan, libmpi.mpi_tan, "tangent")
cosh = def_mp_builtin('cosh', libelefun.mpf_cosh, libmpc.mpc_cosh, None, "hyperbolic cosine")
sinh = def_mp_builtin('sinh', libelefun.mpf_sinh, libmpc.mpc_sinh, None, "hyperbolic sine")
tanh = def_mp_builtin('tanh', libelefun.mpf_tanh, libmpc.mpc_tanh, None, "hyperbolic tangent")
acos = def_mp_builtin('acos', libelefun.mpf_acos, libmpc.mpc_acos, None, "inverse cosine")
asin = def_mp_builtin('asin', libelefun.mpf_asin, libmpc.mpc_asin, None, "inverse sine")
atan = def_mp_builtin('atan', libelefun.mpf_atan, libmpc.mpc_atan, None, "inverse tangent")
asinh = def_mp_builtin('asinh', libelefun.mpf_asinh, libmpc.mpc_asinh, None, "inverse hyperbolic sine")
acosh = def_mp_builtin('acosh', libelefun.mpf_acosh, libmpc.mpc_acosh, None, "inverse hyperbolic cosine")
atanh = def_mp_builtin('atanh', libelefun.mpf_atanh, libmpc.mpc_atanh, None, "inverse hyperbolic tangent")
cospi = def_mp_builtin('cospi', libelefun.mpf_cos_pi, libmpc.mpc_cos_pi, None, "")
sinpi = def_mp_builtin('sinpi', libelefun.mpf_sin_pi, libmpc.mpc_sin_pi, None, "")
floor = def_mp_builtin('floor', libmpf.mpf_floor, libmpc.mpc_floor, None, "")
ceil = def_mp_builtin('ceil', libmpf.mpf_ceil, libmpc.mpc_ceil, None, "")
fibonacci = def_mp_builtin('fibonacci', libelefun.mpf_fibonacci, libmpc.mpc_fibonacci, None, "")
zeta = def_mp_builtin('zeta', gammazeta.mpf_zeta, gammazeta.mpc_zeta, None, "Riemann zeta function")
altzeta = def_mp_builtin('altzeta', gammazeta.mpf_altzeta, gammazeta.mpc_altzeta, None, "Dirichlet eta function")
gamma = def_mp_builtin('gamma', gammazeta.mpf_gamma, gammazeta.mpc_gamma, None, "gamma function")
factorial = def_mp_builtin('factorial', gammazeta.mpf_factorial, gammazeta.mpc_factorial, None, "factorial")
harmonic = def_mp_builtin('harmonic', gammazeta.mpf_harmonic, gammazeta.mpc_harmonic, None, "nth harmonic number")
erf = def_mp_builtin("erf", libhyper.mpf_erf, libhyper.mpc_erf, None, "Error function, erf(z)")
erfc = def_mp_builtin("erfc", libhyper.mpf_erfc, libhyper.mpc_erfc, None, "Complementary error function, erfc(z) = 1-erf(z)")
ci = def_mp_builtin('ci', libhyper.mpf_ci, libhyper.mpc_ci, None, "")
si = def_mp_builtin('si', libhyper.mpf_si, libhyper.mpc_si, None, "")
ellipk = def_mp_builtin('ellipk', libhyper.mpf_ellipk, libhyper.mpc_ellipk, None, "")
ellipe = def_mp_builtin('ellipe', libhyper.mpf_ellipe, libhyper.mpc_ellipe, None, "")
agm1 = def_mp_builtin('agm1', libhyper.mpf_agm1, libhyper.mpc_agm1, None, "Fast alias for agm(1,a) = agm(a,1)")

fac = MultiPrecisionArithmetic.fac = factorial
fib = MultiPrecisionArithmetic.fib = fibonacci


# The main reason why each constant is a class and not just an instance
# is that Sphinx won't show docstrings for single instances

def defconst(name, func, descr):
    MultiPrecisionArithmetic._constants.append((name, func, descr))

defconst("pi", libelefun.mpf_pi, "pi")
defconst("degree", libelefun.mpf_degree, "degree")
defconst("e", libelefun.mpf_e, "e")
defconst("ln2", libelefun.mpf_ln2, "ln(2)")
defconst("ln10", libelefun.mpf_ln10, "ln(10)")
defconst("phi", libelefun.mpf_phi, "Golden ratio (phi)")
defconst("euler", gammazeta.mpf_euler, "Euler's constant (gamma)")
defconst("catalan", gammazeta.mpf_catalan, "Catalan's constant")
defconst("khinchin", gammazeta.mpf_khinchin, "Khinchin's constant")
defconst("glaisher", gammazeta.mpf_glaisher, "Glaisher's constant")
defconst("apery", gammazeta.mpf_apery, "Apery's constant")
defconst("mertens", gammazeta.mpf_mertens, "Mertens' constant")
defconst("twinprime", gammazeta.mpf_twinprime, "Twin prime constant")

mp._create_constants(globals())

def funcwrapper(f):
    def g(*args, **kwargs):
        orig = mp.prec
        try:
            args = [mp.mpmathify(z) for z in args]
            mp.prec = orig + 10
            v = f(*args, **kwargs)
        finally:
            mp.prec = orig
        return +v
    g.__name__ = f.__name__
    g.__doc__ = f.__doc__
    return g

def altfunc(f, name, desc):
    def g(self, x):
        orig = self.prec
        try:
            self.prec = orig + 10
            return self.one/f(x)
        finally:
            self.prec = orig
    g.__name__ = name
    g.__doc__ = "Computes the %s of x, 1/%s(x)" % (desc, f.__name__)
    return defun(g)

def altinvfunc(f, name, desc):
    def g(self, x):
        orig = self.prec
        try:
            self.prec = orig + 10
            return f(self.one/x)
        finally:
            self.prec = orig
    setattr(MultiPrecisionArithmetic, name, g)
    g.__name__ = name
    g.__doc__ = "Computes the inverse %s of x, %s(1/x)" % (desc, f.__name__)
    return defun(g)

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


@defun_wrapped
def sinc(ctx, x):
    if ctx.isinf(x):
        return 1/x
    if not x:
        return x+1
    return ctx.sin(x)/x

@defun_wrapped
def sincpi(ctx, x):
    if ctx.isinf(x):
        return 1/x
    if not x:
        return x+1
    return ctx.sinpi(x)/(ctx.pi*x)

@defun
def nthroot(ctx, x, n):
    x = ctx.convert(x)
    n = int(n)
    if hasattr(x, '_mpf_'):
        try:
            return ctx.make_mpf(libelefun.mpf_nthroot(x._mpf_, n, *ctx._prec_rounding))
        except ComplexResult:
            if ctx.trap_complex:
                raise
            x = (x._mpf_, libmpf.fzero)
    else:
        x = x._mpc_
    return ctx.make_mpc(libmpc.mpc_nthroot(x, n, *ctx._prec_rounding))

@defun
def hypot(ctx, x, y):
    r"""
    Computes the Euclidean norm of the vector `(x, y)`, equal
    to `\sqrt{x^2 + y^2}`. Both `x` and `y` must be real."""
    x = ctx.convert(x)
    y = ctx.convert(y)
    return ctx.make_mpf(libmpf.mpf_hypot(x._mpf_, y._mpf_, *ctx._prec_rounding))

@defun
def ldexp(ctx, x, n):
    r"""
    Computes `x 2^n` efficiently. No rounding is performed.
    The argument `x` must be a real floating-point number (or
    possible to convert into one) and `n` must be a Python ``int``.

        >>> from mpmath import *
        >>> ldexp(1, 10)
        mpf('1024.0')
        >>> ldexp(1, -3)
        mpf('0.125')

    """
    x = ctx.convert(x)
    return ctx.make_mpf(libmpf.mpf_shift(x._mpf_, n))

@defun
def frexp(ctx, x):
    r"""
    Given a real number `x`, returns `(y, n)` with `y \in [0.5, 1)`,
    `n` a Python integer, and such that `x = y 2^n`. No rounding is
    performed.

        >>> from mpmath import *
        >>> frexp(7.5)
        (mpf('0.9375'), 3)

    """
    x = ctx.convert(x)
    y, n = libmpf.mpf_frexp(x._mpf_)
    return ctx.make_mpf(y), n

@defun
def sign(ctx, x):
    r"""
    Returns the sign of `x`, defined as `\mathrm{sign}(x) = x / |x|`
    (with the special case `\sign(0) = 0`)::

        >>> from mpmath import *
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
    x = ctx.convert(x)
    if not x or ctx.isnan(x):
        return x
    if ctx.is_real_type(x):
        return ctx.mpf(cmp(x, 0))
    return x / abs(x)

@defun
def arg(ctx, x):
    r"""
    Computes the complex argument (phase) of `x`, defined as the
    signed angle between the positive real axis and `x` in the
    complex plane::

        >>> from mpmath import *
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
    x = ctx.convert(x)
    return ctx.atan2(x.imag, x.real)

@defun
def fabs(ctx, x):
    r"""
    Returns the absolute value of `x`, `|x|`. Unlike :func:`abs`,
    :func:`fabs` converts non-mpmath numbers (such as ``int``)
    into mpmath numbers::

        >>> from mpmath import *
        >>> fabs(3)
        mpf('3.0')
        >>> fabs(-3)
        mpf('3.0')
        >>> fabs(3+4j)
        mpf('5.0')

    """
    return abs(ctx.convert(x))

@defun
def re(ctx, x):
    r"""
    Returns the real part of `x`, `\Re(x)`. Unlike ``x.real``,
    :func:`re` converts `x` to a mpmath number::

        >>> from mpmath import *
        >>> re(3)
        mpf('3.0')
        >>> re(-1+4j)
        mpf('-1.0')

    """
    return ctx.convert(x).real

@defun
def im(ctx, x):
    r"""
    Returns the imaginary part of `x`, `\Im(x)`. Unlike ``x.imag``,
    :func:`im` converts `x` to a mpmath number::

        >>> from mpmath import *
        >>> im(3)
        mpf('0.0')
        >>> im(-1+4j)
        mpf('4.0')

    """
    return ctx.convert(x).imag

@defun
def conj(ctx, x):
    r"""
    Returns the complex conjugate of `x`, `\overline{x}`. Unlike
    ``x.conjugate()``, :func:`im` converts `x` to a mpmath number::

        >>> from mpmath import *
        >>> conj(3)
        mpf('3.0')
        >>> conj(-1+4j)
        mpc(real='-1.0', imag='-4.0')

    """
    return ctx.convert(x).conjugate()


@defun
def log(ctx, x, b=None):
    if b is None:
        return ln(x)
    wp = ctx.prec + 20
    return ctx.ln(x, prec=wp) / ctx.ln(b, prec=wp)

@defun
def log10(ctx, x):
    r"""
    Computes the base-10 logarithm of `x`, `\log_{10}(x)`. ``log10(x)``
    is equivalent to ``log(x, 10)``.
    """
    return ctx.log(x, 10)

@defun
def power(ctx, x, y):
    return ctx.convert(x) ** ctx.convert(y)

@defun
def modf(ctx,x,y):
    return ctx.convert(x) % ctx.convert(y)

@defun
def degrees(ctx,x):
    return x / ctx.degree

@defun
def radians(ctx,x):
    return x * ctx.degree

@defun
def atan2(ctx, y, x):
    x = ctx.convert(x)
    y = ctx.convert(y)
    return ctx.make_mpf(libelefun.mpf_atan2(y._mpf_, x._mpf_, *ctx._prec_rounding))

@defun
def psi(ctx, m, z):
    z = ctx.convert(z)
    m = int(m)
    if ctx.is_real_type(z):
        return ctx.make_mpf(gammazeta.mpf_psi(m, z._mpf_, *ctx._prec_rounding))
    else:
        return ctx.make_mpc(gammazeta.mpc_psi(m, z._mpc_, *ctx._prec_rounding))

@defun
def psi0(ctx, z):
    """Shortcut for psi(0,z) (the digamma function)"""
    return ctx.psi(0, z)

@defun
def psi1(ctx, z):
    """Shortcut for psi(1,z) (the trigamma function)"""
    return ctx.psi(1, z)

@defun
def psi2(ctx, z):
    """Shortcut for psi(2,z) (the tetragamma function)"""
    return ctx.psi(2, z)

@defun
def psi3(ctx, z):
    """Shortcut for psi(3,z) (the pentagamma function)"""
    return ctx.psi(3, z)

polygamma = MultiPrecisionArithmetic.polygamma = psi
digamma = MultiPrecisionArithmetic.digamma = psi0
trigamma = MultiPrecisionArithmetic.trigamma = psi1
tetragamma = MultiPrecisionArithmetic.tetragamma = psi2
pentagamma = MultiPrecisionArithmetic.pentagamma = psi3

@defun
def bernoulli(ctx, n):
    return ctx.make_mpf(gammazeta.mpf_bernoulli(int(n), *ctx._prec_rounding))

bernfrac = defun_static(gammazeta.bernfrac)

@defun
def stieltjes(ctx, n, a=1):
    n = ctx.convert(n)
    a = ctx.convert(a)
    if n < 0:
        return ctx.bad_domain("Stieltjes constants defined for n >= 0")
    if hasattr(ctx, "stieltjes_cache"):
        stieltjes_cache = ctx.stieltjes_cache
    else:
        stieltjes_cache = ctx.stieltjes_cache = {}
    if a == 1:
        if n == 0:
            return +ctx.euler
        if n in stieltjes_cache:
            prec, s = stieltjes_cache[n]
            if prec >= ctx.prec:
                return +s
    mag = 1
    def f(x):
        xa = x/a
        v = (xa-ctx.j)*ctx.ln(a-ctx.j*x)**n/(1+xa**2)/(ctx.exp(2*ctx.pi*x)-1)
        return v.real / mag
    orig = ctx.prec
    try:
        # Normalize integrand by approx. magnitude to
        # speed up quadrature (which uses absolute error)
        if n > 50:
            ctx.prec = 20
            mag = ctx.quad(f, [0,ctx.inf], maxdegree=3)
        ctx.prec = orig + 10 + int(n**0.5)
        s = ctx.quad(f, [0,ctx.inf], maxdegree=20)
        v = ctx.ln(a)**n/(2*a) - ctx.ln(a)**(n+1)/(n+1) + 2*s/a*mag
    finally:
        ctx.prec = orig
    if a == 1 and ctx.isint(n):
        stieltjes_cache[n] = (ctx.prec, v)
    return +v

@defun
def gammaprod(ctx, a, b):
    a = [ctx.convert(x) for x in a]
    b = [ctx.convert(x) for x in b]
    poles_num = []
    poles_den = []
    regular_num = []
    regular_den = []
    for x in a: [regular_num, poles_num][ctx.isnpint(x)].append(x)
    for x in b: [regular_den, poles_den][ctx.isnpint(x)].append(x)
    # One more pole in numerator or denominator gives 0 or inf
    if len(poles_num) < len(poles_den): return ctx.zero
    if len(poles_num) > len(poles_den): return ctx.inf
    # All poles cancel
    # lim G(i)/G(j) = (-1)**(i+j) * gamma(1-j) / gamma(1-i)
    p = ctx.one
    orig = ctx.prec
    try:
        ctx.prec = orig + 15
        while poles_num:
            i = poles_num.pop()
            j = poles_den.pop()
            p *= (-1)**(i+j) * ctx.gamma(1-j) / ctx.gamma(1-i)
        for x in regular_num: p *= ctx.gamma(x)
        for x in regular_den: p /= ctx.gamma(x)
    finally:
        ctx.prec = orig
    return +p

@defun
def beta(ctx, x, y):
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

        >>> from mpmath import *
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
    x = ctx.convert(x)
    y = ctx.convert(y)
    if ctx.isinf(y):
        x, y = y, x
    if ctx.isinf(x):
        if x == ctx.inf and not y.imag:
            if y == ctx.ninf:
                return ctx.nan
            if y > 0:
                return ctx.zero
            if ctx.isint(y):
                return ctx.nan
            if y < 0:
                return ctx.sign(ctx.gamma(y)) * ctx.inf
        return ctx.nan
    return ctx.gammaprod([x, y], [x+y])

@defun
def binomial(ctx, n, k):
    return ctx.gammaprod([n+1], [k+1, n-k+1])

@defun
def rf(ctx, x, n):
    return ctx.gammaprod([x+n], [x])

@defun
def ff(ctx, x, n):
    return ctx.gammaprod([x+1], [x-n+1])

@defun_wrapped
def fac2(ctx, x):
    if ctx.isinf(x):
        if x == ctx.inf:
            return x
        return ctx.nan
    return 2**(x/2)*(ctx.pi/2)**((ctx.cospi(x)-1)/4)*ctx.gamma(x/2+1)


#---------------------------------------------------------------------------#
#                                                                           #
#                          Hypergeometric functions                         #
#                                                                           #
#---------------------------------------------------------------------------#

from libmpf import from_rational

class _mpq(tuple):

    def _mpmath_(self, prec, rounding):
        # XXX
        return mp.make_mpf(from_rational(self[0], self[1], prec, rounding))
        #(mpf(self[0])/self[1])._mpf_

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

@defun
def _hyp_parse_param(ctx, x):
    if isinstance(x, tuple):
        p, q = x
        return [[p, q]], [], []
    if isinstance(x, (int, long)):
        return [[x, 1]], [], []
    x = ctx.convert(x)
    if hasattr(x, '_mpf_'):
        return [], [x._mpf_], []
    if hasattr(x, '_mpc_'):
        return [], [], [x._mpc_]

def _as_num(x):
    if isinstance(x, list):
        return _mpq(x)
    return x

@defun
def hypsum(ctx, ar, af, ac, br, bf, bc, x):
    prec, rnd = ctx._prec_rounding
    if hasattr(x, '_mpf_') and not (ac or bc):
        v = libhyper.hypsum_internal(ar, af, ac, br, bf, bc, x._mpf_, None, prec, rnd)
        return ctx.make_mpf(v)
    else:
        if hasattr(x, '_mpc_'):
            re, im = x._mpc_
        else:
            re, im = x._mpf_, libmpf.fzero
        v = libhyper.hypsum_internal(ar, af, ac, br, bf, bc, re, im, prec, rnd)
        return ctx.make_mpc(v)

@defun
def eval_hyp2f1(ctx,a,b,c,z):
    prec, rnd = ctx._prec_rounding
    ar, af, ac = ctx._hyp_parse_param(a)
    br, bf, bc = ctx._hyp_parse_param(b)
    cr, cf, cc = ctx._hyp_parse_param(c)
    absz = abs(z)
    if absz == 1:
        # TODO: determine whether it actually does, and otherwise
        # return infinity instead
        print "Warning: 2F1 might not converge for |z| = 1"
    if absz <= 1:
        # All rational
        if ar and br and cr:
            return ctx.sum_hyp2f1_rat(ar[0], br[0], cr[0], z)
        return ctx.hypsum(ar+br, af+bf, ac+bc, cr, cf, cc, z)
    # Use 1/z transformation
    a = (ar and _as_num(ar[0])) or ctx.convert(a)
    b = (br and _as_num(br[0])) or ctx.convert(b)
    c = (cr and _as_num(cr[0])) or ctx.convert(c)
    orig = ctx.prec
    try:
        ctx.prec = orig + 15
        h1 = ctx.eval_hyp2f1(a, mpq_1-c+a, mpq_1-b+a, 1/z)
        h2 = ctx.eval_hyp2f1(b, mpq_1-c+b, mpq_1-a+b, 1/z)
        #s1 = G(c)*G(b-a)/G(b)/G(c-a) * (-z)**(-a) * h1
        #s2 = G(c)*G(a-b)/G(a)/G(c-b) * (-z)**(-b) * h2
        f1 = ctx.gammaprod([c,b-a],[b,c-a])
        f2 = ctx.gammaprod([c,a-b],[a,c-b])
        s1 = f1 * (-z)**(mpq_0-a) * h1
        s2 = f2 * (-z)**(mpq_0-b) * h2
        v = s1 + s2
    finally:
        ctx.prec = orig
    return +v

@defun
def sum_hyp0f1_rat(ctx, a, z):
    prec, rnd = ctx._prec_rounding
    if hasattr(z, "_mpf_"):
        return ctx.make_mpf(libhyper.mpf_hyp0f1_rat(a, z._mpf_, prec, rnd))
    else:
        return ctx.make_mpc(libhyper.mpc_hyp0f1_rat(a, z._mpc_, prec, rnd))

@defun
def sum_hyp1f1_rat(ctx, a, b, z):
    prec, rnd = ctx._prec_rounding
    if hasattr(z, "_mpf_"):
        return ctx.make_mpf(libhyper.mpf_hyp1f1_rat(a, b, z._mpf_, prec, rnd))
    else:
        return ctx.make_mpc(libhyper.mpc_hyp1f1_rat(a, b, z._mpc_, prec, rnd))

@defun
def sum_hyp2f1_rat(ctx, a, b, c, z):
    prec, rnd = ctx._prec_rounding
    if hasattr(z, "_mpf_"):
        return ctx.make_mpf(libhyper.mpf_hyp2f1_rat(a, b, c, z._mpf_, prec, rnd))
    else:
        return ctx.make_mpc(libhyper.mpc_hyp2f1_rat(a, b, c, z._mpc_, prec, rnd))


#---------------------------------------------------------------------------#
#                      And now the user-friendly versions                   #
#---------------------------------------------------------------------------#

@defun
def hyper(ctx, a_s, b_s, z):
    p = len(a_s)
    q = len(b_s)
    z = ctx.convert(z)
    degree = p, q
    if degree == (0, 1):
        br, bf, bc = ctx._hyp_parse_param(b_s[0])
        if br:
            return ctx.sum_hyp0f1_rat(br[0], z)
        return ctx.hypsum([], [], [], br, bf, bc, z)
    if degree == (1, 1):
        ar, af, ac = ctx._hyp_parse_param(a_s[0])
        br, bf, bc = ctx._hyp_parse_param(b_s[0])
        if ar and br:
            a, b = ar[0], br[0]
            return ctx.sum_hyp1f1_rat(a, b, z)
        return ctx.hypsum(ar, af, ac, br, bf, bc, z)
    if degree == (2, 1):
        return ctx.eval_hyp2f1(a_s[0], a_s[1], b_s[0], z)
    ars, afs, acs, brs, bfs, bcs = [], [], [], [], [], []
    for a in a_s:
        r, f, c = ctx._hyp_parse_param(a)
        ars += r
        afs += f
        acs += c
    for b in b_s:
        r, f, c = ctx._hyp_parse_param(b)
        brs += r
        bfs += f
        bcs += c
    return ctx.hypsum(ars, afs, acs, brs, bfs, bcs, z)

@defun
def hyp0f1(ctx, a, z):
    r"""Hypergeometric function `\,_0F_1`. ``hyp0f1(a,z)`` is equivalent
    to ``hyper([],[a],z)``; see documentation for :func:`hyper` for more
    information."""
    return ctx.hyper([], [a], z)

@defun
def hyp1f1(ctx,a,b,z):
    r"""Hypergeometric function `\,_1F_1`. ``hyp1f1(a,b,z)`` is equivalent
    to ``hyper([a],[b],z)``; see documentation for :func:`hyper` for more
    information."""
    return ctx.hyper([a], [b], z)

@defun
def hyp2f1(ctx,a,b,c,z):
    r"""Hypergeometric function `\,_2F_1`. ``hyp2f1(a,b,c,z)`` is equivalent
    to ``hyper([a,b],[c],z)``; see documentation for :func:`hyper` for more
    information."""
    return ctx.hyper([a,b], [c], z)

@defun
def _lower_gamma(ctx, z, b):
    return ctx.hyp1f1(1, 1+z, b) * b**z * ctx.exp(-b) / z

def _check_pos(x):
    try:
        return x > 0
    except TypeError:
        return False

@defun_wrapped
def gammainc(ctx, z, a=0, b=None, regularized=False):
    if b is None:
        b = ctx.inf
    ln = ctx.ln
    ei = ctx.ei
    if b == ctx.inf:
        if not a:
            v = ctx.gamma(z)
        else:
            if not z:
                # Reduces to exponential integral. Mind branch cuts.
                if _check_pos(a):
                    return -ei(-a)
                else:
                    return -ei(-a) + (ln(-a)-ln(-1/a))/2-ln(a)
            # XXX: avoid poles
            v = ctx.gamma(z) - ctx._lower_gamma(z, a)
    elif not a:
        v = ctx._lower_gamma(z, b)
    else:
        if not z:
            # Reduces to exponential integral
            if _check_pos(a) and _check_pos(b):
                return ei(-b) - ei(-a)
            else:
                return ei(-b)-ei(-a) + \
                    (ln(-a)-ln(-1/a))/2-ln(a) + \
                    (ln(-1/b)-ln(-b))/2+ln(b)
        # XXX: avoid poles
        v = ctx._lower_gamma(z, b) - ctx._lower_gamma(z, a)
    if regularized:
        return v / ctx.gamma(z)
    else:
        return v

@defun_wrapped
def erfi(ctx, z):
    return (2/ctx.sqrt(ctx.pi)*z) * ctx.sum_hyp1f1_rat((1,2),(3,2), z**2)

@defun_wrapped
def erfinv(ctx, x):
    if x.imag or (x < -1) or (x > 1):
        return ctx.bad_domain("erfinv(x) is defined only for -1 <= x <= 1")
    if ctx.isnan(x): return x
    if not x: return x
    if x == 1: return ctx.inf
    if x == -1: return ctx.ninf
    if abs(x) < 0.9:
        a = 0.53728*x**3 + 0.813198*x
    else:
        # An asymptotic formula
        u = ctx.ln(2/ctx.pi/(abs(x)-1)**2)
        a = ctx.sign(x) * ctx.sqrt(u - ctx.ln(u))/ctx.sqrt(2)
    return ctx.findroot(lambda t: ctx.erf(t)-x, a)

@defun_wrapped
def npdf(ctx, x, mu=0, sigma=1):
    sigma = ctx.convert(sigma)
    return ctx.exp(-(x-mu)**2/(2*sigma**2)) / (sigma*ctx.sqrt(2*ctx.pi))

@defun_wrapped
def ncdf(ctx, x, mu=0, sigma=1):
    a = (x-mu)/(sigma*ctx.sqrt(2))
    if a < 0:
        return ctx.erfc(-a)/2
    else:
        return (1+ctx.erf(a))/2

def ei_as(ctx, a):
    extra = 10
    ctx.dps += extra
    s = k = p = 1
    while abs(p) > ctx.eps:
        p = (p*k)/a
        s += p
        k += 1
    s = (s * ctx.exp(a))/a
    ctx.dps -= extra
    return s

@defun_wrapped
def ei(ctx, z):
    if z == ctx.inf:
        return z
    if z == ctx.ninf:
        return -ctx.zero
    if not z:
        return ctx.ninf
    if abs(z) > ctx.prec * 0.7 + 50:
        r = ei_as(ctx, z)
        if z.imag > 0:
            r += ctx.j*ctx.pi
        elif z.imag < 0:
            r -= ctx.j*ctx.pi
        return r
    v = z*hypsum([[1,1],[1,1]],[],[],[[2,1],[2,1]],[],[],z) + \
        (ctx.ln(z)-ctx.ln(1/z))/2 + ctx.euler
    if ctx.is_real_type(z) and z < 0:
        return v.real
    return v

@defun_wrapped
def li(ctx, z):
    if not z:
        return z
    if z == 1:
        return ctx.ninf
    return ctx.ei(ctx.ln(z))

@defun_wrapped
def chi(ctx, z):
    if not z:
        return ctx.ninf
    z2 = (z/2)**2
    return ctx.euler + ctx.ln(z) + \
        z2*ctx.hypsum([[1,1],[1,1]],[],[],[[2,1],[2,1],[3,2]],[],[],z2)

@defun_wrapped
def shi(ctx, z):
    z2 = (z/2)**2
    return z*ctx.hypsum([[1,2]],[],[],[[3,2],[3,2]],[],[],z2)

@defun_wrapped
def fresnels(ctx, z):
    if z == ctx.inf:
        return ctx.mpf(0.5)
    if z == ctx.ninf:
        return ctx.mpf(-0.5)
    return ctx.pi*z**3/6*ctx.hypsum([[3,4]],[],[],[[3,2],[7,4]],[],[],-ctx.pi**2*z**4/16)

@defun_wrapped
def fresnelc(ctx, z):
    if z == ctx.inf:
        return ctx.mpf(0.5)
    if z == ctx.ninf:
        return ctx.mpf(-0.5)
    return z*ctx.hypsum([[1,4]],[],[],[[1,2],[5,4]],[],[],-ctx.pi**2*z**4/16)

@defun_wrapped
def airyai(ctx, z):
    if z == ctx.inf or z == ctx.ninf:
        return 1/z
    if z.real > 2:
        # cancellation: both terms are ~ 2^(z^1.5),
        # result is ~ 2^(-z^1.5), so need ~2*z^1.5 extra bits
        ctx.prec += 2*int(z.real**1.5)
    z3 = z**3 / 9
    a = ctx.sum_hyp0f1_rat((2,3), z3) / (ctx.cbrt(9) * ctx.gamma(ctx.mpf(2)/3))
    b = z * ctx.sum_hyp0f1_rat((4,3), z3) / (ctx.cbrt(3) * ctx.gamma(ctx.mpf(1)/3))
    return a - b

@defun_wrapped
def airybi(ctx, z):
    if z == ctx.inf:
        return z
    if z == ctx.ninf:
        return 1/z
    z3 = z**3 / 9
    rt = ctx.nthroot(3, 6)
    a = ctx.sum_hyp0f1_rat((2,3), z3) / (rt * ctx.gamma(ctx.mpf(2)/3))
    b = z * rt * ctx.sum_hyp0f1_rat((4,3), z3) / ctx.gamma(ctx.mpf(1)/3)
    return a + b

@defun
def agm(ctx, a, b=1):
    if b == 1:
        return ctx.agm1(a)
    a = ctx.convert(a)
    b = ctx.convert(b)
    prec, rounding = ctx._prec_rounding
    if hasattr(a, '_mpf_') and hasattr(b, '_mpf_'):
        try:
            v = libhyper.mpf_agm(a._mpf_, b._mpf_, prec, rounding)
            return ctx.make_mpf(v)
        except ComplexResult:
            pass
    if hasattr(a, '_mpf_'): a = (a._mpf_, libmpf.fzero)
    else: a = a._mpc_
    if hasattr(b, '_mpf_'): b = (b._mpf_, libmpf.fzero)
    else: b = b._mpc_
    return ctx.make_mpc(libhyper.mpc_agm(a, b, prec, rounding))

@defun_wrapped
def jacobi(ctx, n, a, b, x):
    return ctx.binomial(n+a,n) * ctx.hyp2f1(-n,1+n+a+b,a+1,(1-x)/2)

@defun_wrapped
def legendre(ctx, n, x):
    if ctx.isint(n):
        n = int(n)
    if x == -1:
        # TODO: hyp2f1 should handle this
        if ctx.isint(n):
            return (-1)**(n + (n>=0)) * ctx.mpf(-1)
        if not int(ctx.floor(ctx.re(n))) % 2:
            return ctx.ninf
        return ctx.inf
    return ctx.hyp2f1(-n,n+1,1,(1-x)/2)

@defun_wrapped
def chebyt(ctx, n, x):
    return ctx.hyp2f1(-n,n,0.5,(1-x)/2)

@defun_wrapped
def chebyu(ctx, n, x):
    return (n+1) * ctx.hyp2f1(-n, n+2, 1.5, (1-x)/2)

@defun_wrapped
def _besselj(ctx, v, x):
    hx = x/2
    return hx**v * ctx.hyp0f1(v+1, -hx**2) / ctx.factorial(v)

@defun
def besselj(ctx, v, x):
    if ctx.isint(v):
        x = ctx.convert(x)
        v = int(v)
        prec, rounding = ctx._prec_rounding
        if hasattr(x, '_mpf_'):
            return ctx.make_mpf(libhyper.mpf_besseljn(v, x._mpf_, prec, rounding))
        if hasattr(x, '_mpc_'):
            return ctx.make_mpc(libhyper.mpc_besseljn(v, x._mpc_, prec, rounding))
    return ctx._besselj(v, x)

@defun
def j0(ctx, x):
    """Computes the Bessel function `J_0(x)`. See :func:`besselj`."""
    return ctx.besselj(0, x)

@defun
def j1(ctx, x):
    """Computes the Bessel function `J_1(x)`.  See :func:`besselj`."""
    return ctx.besselj(1, x)

@defun_wrapped
def bessely(ctx,n,x):
    intdist = abs(n.imag) + abs(n.real-ctx.floor(n.real+0.5))
    if not intdist:
        h = +ctx.eps
        ctx.prec *= 2
        n += h
    else:
        ctx.prec += -int(ctx.log(intdist, 2)+1)
    return (ctx.besselj(n,x)*ctx.cospi(n) - ctx.besselj(-n,x))/ctx.sinpi(n)

@defun_wrapped
def besseli(ctx,n,x):
    if ctx.isint(n):
        n = abs(int(n))
    hx = x/2
    return hx**n * ctx.hyp0f1(n+1, hx**2) / ctx.factorial(n)

@defun_wrapped
def besselk(ctx,n,x):
    intdist = abs(n.imag) + abs(n.real-ctx.floor(n.real+0.5))
    if not intdist:
        h = +ctx.eps
        ctx.prec *= 2
        n += h
    else:
        ctx.prec += -int(ctx.log(intdist, 2)+1)
    return ctx.pi*(ctx.besseli(-n,x)-ctx.besseli(n,x))/(2*ctx.sinpi(n))

@defun_wrapped
def hankel1(ctx,n,x):
    return ctx.besselj(n,x) + ctx.j*ctx.bessely(n,x)

@defun_wrapped
def hankel2(ctx,n,x):
    return ctx.besselj(n,x) - ctx.j*ctx.bessely(n,x)

@defun_wrapped
def lambertw(ctx, z, k=0, approx=None):
    if ctx.isnan(z):
        return z
    ctx.prec += 20
    # We must be extremely careful near the singularities at -1/e and 0
    u = ctx.exp(-1)
    if abs(z) <= u:
        if not z:
            # w(0,0) = 0; for all other branches we hit the pole
            if not k:
                return z
            return ctx.ninf
        if not k:
            w = z
        # For small real z < 0, the -1 branch behaves roughly like log(-z)
        elif k == -1 and not ctx.im(z) and ctx.re(z) < 0:
            w = ctx.ln(-z)
        # Use a simple asymptotic approximation.
        else:
            w = ctx.ln(z)
            # The branches are roughly logarithmic. This approximation
            # gets better for large |k|; need to check that this always
            # works for k ~= -1, 0, 1.
            if k: w += k * 2*ctx.pi*ctx.j
    elif k == 0 and ctx.im(z) and abs(z) <= 0.6:
        w = z
    else:
        if z == ctx.inf:
            if k == 0:
                return z
            else:
                return z + 2*k*ctx.pi*ctx.j
        if z == ctx.ninf:
            return (-z) + (2*k+1)*ctx.pi*ctx.j
        # Simple asymptotic approximation as above
        w = ctx.ln(z)
        if k: w += k * 2*ctx.pi*ctx.j
    # Use Halley iteration to solve w*exp(w) = z
    two = ctx.mpf(2)
    weps = ctx.ldexp(ctx.eps, 15)
    for i in xrange(100):
        ew = ctx.exp(w)
        wew = w*ew
        wewz = wew-z
        wn = w - wewz/(wew+ew-(w+two)*wewz/(two*w+two))
        if abs(wn-w) < weps*abs(wn):
            return wn
        else:
            w = wn
    print "Warning: Lambert W iteration failed to converge:", z
    return wn

@defun_wrapped
def barnesg(ctx, z):
    if ctx.isinf(z):
        if z == ctx.inf:
            return z
        return ctx.nan
    if ctx.isnan(z):
        return z
    if (not z.imag) and z.real <= 0 and ctx.isint(z.real):
        return z*0
    # Account for size (would not be needed if computing log(G))
    if abs(z) > 5:
        ctx.dps += 2*ctx.log(abs(z),2)
    # Estimate terms for asymptotic expansion
    N = ctx.dps // 2 + 5
    G = 1
    while ctx.re(z) < N:
        G /= ctx.gamma(z)
        z += 1
    z -= 1
    s = ctx.mpf(1)/12
    s -= ctx.log(ctx.glaisher)
    s += z*ctx.log(2*ctx.pi)/2
    s += (z**2/2-ctx.mpf(1)/12)*ctx.log(z)
    s -= 3*z**2/4
    z2k = z2 = z**2
    for k in xrange(1, N+1):
        t = ctx.bernoulli(2*k+2) / (4*k*(k+1)*z2k)
        if abs(t) < ctx.eps:
            #print k, N      # check how many terms were needed
            break
        z2k *= z2
        s += t
    #if k == N:
    #    print "warning: series for barnesg failed to converge"
    return G*ctx.exp(s)

@defun
def superfac(ctx, z):
    return ctx.barnesg(z+2)

@defun_wrapped
def hyperfac(ctx, z):
    # XXX: estimate needed extra bits accurately
    if z == ctx.inf:
        return z
    if abs(z) > 5:
        extra = 4*int(ctx.log(abs(z),2))
    else:
        extra = 0
    ctx.prec += extra
    if not z.imag and z.real < 0 and ctx.isint(z.real):
        n = int(ctx.re(z))
        h = ctx.hyperfac(-n-1)
        if ((n+1)//2) & 1:
            h = -h
        if ctx.is_complex_type(z):
            return h + 0j
        return h
    zp1 = z+1
    # Wrong branch cut
    #v = ctx.gamma(zp1)**z
    #ctx.prec -= extra
    #return v / ctx.barnesg(zp1)
    v = ctx.exp(z*ctx.loggamma(zp1))
    ctx.prec -= extra
    return v / ctx.barnesg(zp1)

@defun_wrapped
def loggamma(ctx, z):
    a = z.real
    b = z.imag
    if not b and a > 0:
        return ctx.ln(ctx.gamma(z))
    u = ctx.arg(z)
    w = ctx.ln(ctx.gamma(z))
    if b:
        gi = -b - u/2 + a*u + b*ctx.ln(abs(z))
        n = ctx.floor((gi-w.imag)/(2*ctx.pi)+0.5) * (2*ctx.pi)
        return w + n*ctx.j
    elif a < 0:
        n = int(ctx.floor(a))
        w += (n-(n%2))*ctx.pi*ctx.j
    return w

@defun_wrapped
def siegeltheta(ctx, t):
    if t.imag:
        # XXX: cancellation occurs
        a = ctx.loggamma(0.25+0.5j*t)
        b = ctx.loggamma(0.25-0.5j*t)
        return -ctx.ln(ctx.pi)/2*t - 0.5j*(a-b)
    else:
        if ctx.isinf(t):
            return t
        return ctx.loggamma(0.25+0.5j*t).imag - ctx.ln(ctx.pi)/2*t

@defun_wrapped
def grampoint(ctx, n):
    # ctxsymptotic expansion, from
    # http://mathworld.wolfram.com/GramPoint.html
    g = 2*ctx.pi*ctx.exp(1+ctx.lambertw((8*n+1)/(8*ctx.e)))
    return ctx.findroot(lambda t: ctx.siegeltheta(t)-ctx.pi*n, g)

@defun_wrapped
def siegelz(ctx, t):
    v = ctx.exp(ctx.j*ctx.siegeltheta(t))*ctx.zeta(0.5+ctx.j*t)
    if ctx.is_real_type(t):
        return v.real
    return v

_zeta_zeros = [
14.134725142,21.022039639,25.010857580,30.424876126,32.935061588,
37.586178159,40.918719012,43.327073281,48.005150881,49.773832478,
52.970321478,56.446247697,59.347044003,60.831778525,65.112544048,
67.079810529,69.546401711,72.067157674,75.704690699,77.144840069,
79.337375020,82.910380854,84.735492981,87.425274613,88.809111208,
92.491899271,94.651344041,95.870634228,98.831194218,101.317851006,
103.725538040,105.446623052,107.168611184,111.029535543,111.874659177,
114.320220915,116.226680321,118.790782866,121.370125002,122.946829294,
124.256818554,127.516683880,129.578704200,131.087688531,133.497737203,
134.756509753,138.116042055,139.736208952,141.123707404,143.111845808,
146.000982487,147.422765343,150.053520421,150.925257612,153.024693811,
156.112909294,157.597591818,158.849988171,161.188964138,163.030709687,
165.537069188,167.184439978,169.094515416,169.911976479,173.411536520,
174.754191523,176.441434298,178.377407776,179.916484020,182.207078484,
184.874467848,185.598783678,187.228922584,189.416158656,192.026656361,
193.079726604,195.265396680,196.876481841,198.015309676,201.264751944,
202.493594514,204.189671803,205.394697202,207.906258888,209.576509717,
211.690862595,213.347919360,214.547044783,216.169538508,219.067596349,
220.714918839,221.430705555,224.007000255,224.983324670,227.421444280,
229.337413306,231.250188700,231.987235253,233.693404179,236.524229666,
]

def _load_zeta_zeros(url):
    import urllib
    d = urllib.urlopen(url)
    L = [float(x) for x in d.readlines()]
    # Sanity check
    assert round(L[0]) == 14
    _zeta_zeros[:] = L

@defun
def zetazero(ctx, n, url='http://www.dtc.umn.edu/~odlyzko/zeta_tables/zeros1'):
    n = int(n)
    if n < 0:
        return zetazero(-n).conjugate()
    if n == 0:
        raise ValueError("n must be nonzero")
    if n > len(_zeta_zeros) and n <= 100000:
        _load_zeta_zeros(url)
    if n > len(_zeta_zeros):
        raise NotImplementedError("n too large for zetazeros")
    return ctx.mpc(0.5, ctx.findroot(ctx.siegelz, _zeta_zeros[n-1]))

@defun_wrapped
def riemannr(ctx, x):
    if x == 0:
        return ctx.zero
    # Check if a simple asymptotic estimate is accurate enough
    if abs(x) > 1000:
        a = ctx.li(x)
        b = 0.5*ctx.li(ctx.sqrt(x))
        if abs(b) < abs(a)*ctx.eps:
            return a
    if abs(x) < 0.01:
        # XXX
        ctx.prec += int(-ctx.log(abs(x),2))
    # Sum Gram's series
    s = t = ctx.one
    u = ctx.ln(x)
    k = 1
    while abs(t) > abs(s)*ctx.eps:
        t = t * u / k
        s += t / (k * ctx.zeta(k+1))
        k += 1
    return s

@defun_static
def primepi(x):
    x = int(x)
    if x < 2:
        return 0
    from gammazeta import list_primes
    return len(list_primes(x))

@defun_wrapped
def primepi2(ctx, x):
    x = int(x)
    if x < 2:
        return ctx.mpi(0,0)
    if x < 2657:
        return ctx.mpi(ctx.primepi(x))
    mid = ctx.li(x)
    # Schoenfeld's estimate for x >= 2657, assuming RH
    err = ctx.sqrt(x,rounding='u')*ctx.ln(x,rounding='u')/8/ctx.pi(rounding='d')
    a = ctx.floor((ctx.mpi(mid)-err).a, rounding='d')
    b = ctx.ceil((ctx.mpi(mid)+err).b, rounding='u')
    return ctx.mpi(a, b)

@defun_wrapped
def primezeta(ctx, s):
    if ctx.isnan(s):
        return s
    if ctx.re(s) <= 0:
        raise ValueError("prime zeta function defined only for re(s) > 0")
    if s == 1:
        return ctx.inf
    if s == 0.5:
        return ctx.mpc(ctx.ninf, ctx.pi)
    r = ctx.re(s)
    if r > ctx.prec:
        return 0.5**s
    else:
        wp = ctx.prec + int(r)
        def terms():
            orig = ctx.prec
            # zeta ~ 1+eps; need to set precision
            # to get logarithm accurately
            k = 0
            while 1:
                k += 1
                u = libintmath.moebius(k)
                if not u:
                    continue
                ctx.prec = wp
                t = u*ctx.ln(ctx.zeta(k*s))/k
                if not t:
                    return
                #print ctx.prec, ctx.nstr(t)
                ctx.prec = orig
                yield t
    return sum_accurately(ctx, terms)

@defun_wrapped
def bernpoly(ctx, n, z):
    n = int(n)
    assert n >= 0
    # XXX: optimize
    return sum(ctx.binomial(n,k)*ctx.bernoulli(k)*z**(n-k) for k in xrange(0,n+1))

# TODO: this should be implemented low-level
def polylog_series(ctx, s, z):
    tol = +ctx.eps
    l = ctx.zero
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

def polylog_continuation(ctx, n, z):
    if n < 0:
        return z*0
    twopij = 2j * ctx.pi
    a = -twopij**n/ctx.fac(n) * ctx.bernpoly(n, ctx.ln(z)/twopij)
    if ctx.is_real_type(z) and z < 0:
        a = a.real
    if z.imag < 0 or (z.imag == 0 and z.real >= 1):
        a -= twopij*ctx.ln(z)**(n-1)/ctx.fac(n-1)
    return a

def polylog_unitcircle(ctx, n, z):
    tol = +ctx.eps
    if n > 1:
        l = ctx.zero
        logz = ctx.ln(z)
        logmz = ctx.one
        m = 0
        while 1:
            if (n-m) != 1:
                term = ctx.zeta(n-m) * logmz / ctx.fac(m)
                if term and abs(term) < tol:
                    break
                l += term
            logmz *= logz
            m += 1
        l += ctx.ln(z)**(n-1)/ctx.fac(n-1)*(ctx.harmonic(n-1)-ctx.ln(-ctx.ln(z)))
    elif n < 1:  # else
        l = ctx.fac(-n)*(-ctx.ln(z))**(n-1)
        logz = ctx.ln(z)
        logkz = ctx.one
        k = 0
        while 1:
            b = ctx.bernoulli(k-n+1)
            if b:
                term = b*logkz/(ctx.fac(k)*(k-n+1))
                if abs(term) < tol:
                    break
                l -= term
            logkz *= logz
            k += 1
    else:
        raise ValueError
    if ctx.is_real_type(z) and z < 0:
        l = l.real
    return l

@defun_wrapped
def polylog(ctx, s, z):
    if z == 1:
        return ctx.zeta(s)
    if z == -1:
        return -ctx.altzeta(s)
    if s == 0:
        return z/(1-z)
    if s == 1:
        return -ctx.ln(1-z)
    if s == -1:
        return z/(1-z)**2
    if abs(z) <= 0.75 or (not ctx.isint(s) and abs(z) < 0.99):
        return polylog_series(ctx, s, z)
    if abs(z) >= 1.4 and ctx.isint(s):
        return (-1)**(s+1)*polylog_series(ctx, s, 1/z) + polylog_continuation(ctx, s, z)
    if ctx.isint(s):
        return polylog_unitcircle(ctx, int(s), z)
    raise NotImplementedError("polylog for arbitrary s and z")
    # This could perhaps be used in some cases
    #from quadrature import quad
    #return quad(lambda t: t**(s-1)/(exp(t)/z-1),[0,inf])/gamma(s)

# Experimental code; could be used elsewhere
def sum_accurately(ctx, terms, check_step=1):
    orig = ctx.prec
    extra = 10
    while 1:
        ctx.prec = orig + extra
        max_term = ctx.ninf
        s = 0
        k = 0
        for term in terms():
            s += term
            if not k % check_step and term:
                abs_term = abs(term)
                abs_sum = abs(s)
                max_term = max(max_term, abs_term)
                if abs_term <= ctx.eps*abs_sum:
                    break
            k += 1
        if abs_sum:
            cancellation = int(max(0,ctx.log(max_term/abs_sum,2)))
        else:
            cancellation = ctx.prec
        if cancellation < extra:
            break
        else:
            extra += cancellation
    return s

@defun_wrapped
def bell(ctx, n, x=1):
    x = ctx.convert(x)
    if not n:
        if ctx.isnan(x):
            return x
        return type(x)(1)
    if ctx.isinf(x) or ctx.isinf(n) or ctx.isnan(x) or ctx.isnan(n):
        return x**n
    if n == 1: return x
    if n == 2: return x*(x+1)
    if x == 0: return ctx.sincpi(n)
    return _polyexp(ctx, n, x, True) / ctx.exp(x)

def _polyexp(ctx, n, x, extra=False):
    def _terms():
        if extra:
            yield ctx.sincpi(n)
        t = x
        k = 1
        while 1:
            yield k**n * t
            k += 1
            t = t*x/k
    return sum_accurately(ctx, _terms, check_step=4)

@defun_wrapped
def polyexp(ctx, s, z):
    if ctx.isinf(z) or ctx.isinf(s) or ctx.isnan(z) or ctx.isnan(s):
        return z**s
    if z == 0: return z*s
    if s == 0: return ctx.expm1(z)
    if s == 1: return ctx.exp(z)*z
    if s == 2: return ctx.exp(z)*z*(z+1)
    return _polyexp(ctx, s, z)

# TODO: tests; improve implementation
@defun_wrapped
def expm1(ctx, x):
    """
    Accurately computes exp(x)-1.
    """
    if not x:
        return type(x)(1)
    return sum_accurately(ctx, lambda: iter([ctx.exp(x),-1]),1)

# hack
if __name__ == '__main__':
    #import doctest
    #doctest.testmod()
    import doctest
    globs = globals().copy()
    for obj in globs: #sorted(globs.keys()):
        print obj
        doctest.run_docstring_examples(globs[obj], {})
