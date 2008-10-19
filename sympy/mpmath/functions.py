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

from mptypes import (\
    mpnumeric, convert_lossless,
    mpf, make_mpf,
    mpc, make_mpc,
    mpi, make_mpi,
    constant,
    prec_rounding, mp,
    extraprec,
    zero, one, inf, ninf, nan, j, isnan, isinf, isint, eps,
    ComplexResult,
)

# Mathematical constants
pi = constant(libelefun.mpf_pi, "pi")
degree = constant(libelefun.mpf_degree, "degree")
e = constant(libelefun.mpf_e, "e")
ln2 = constant(libelefun.mpf_ln2, "ln(2)")
ln10 = constant(libelefun.mpf_ln10, "ln(10)")
phi = constant(libelefun.mpf_phi, "Golden ratio (phi)")
euler = constant(gammazeta.mpf_euler, "Euler's constant (gamma)")
catalan = constant(gammazeta.mpf_catalan, "Catalan's constant")
khinchin = constant(gammazeta.mpf_khinchin, "Khinchin's constant")
glaisher = constant(gammazeta.mpf_glaisher, "Glaisher's constant")
apery = constant(gammazeta.mpf_apery, "Apery's constant")

def funcwrapper(f):
    def g(*args, **kwargs):
        orig = mp.prec
        try:
            args = [convert_lossless(z) for z in args]
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
            x = convert_lossless(x)
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
    f.__doc__ = "Returns the %s of x" % doc
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
    g.__doc__ = "Returns the %s of x, 1/%s(x)" % (desc, f.__name__)
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
    g.__doc__ = "Returns the inverse %s of x, %s(1/x)" % (desc, f.__name__)
    return g

sqrt = mpfunc('sqrt', libelefun.mpf_sqrt, libmpc.mpc_sqrt, "principal square root", libmpi.mpi_sqrt)
cbrt = mpfunc('cbrt', libelefun.mpf_cbrt, libmpc.mpc_cbrt, "principal cubic root")
exp = mpfunc('exp', libelefun.mpf_exp, libmpc.mpc_exp, "exponential function", libmpi.mpi_exp)
ln = mpfunc('ln', libelefun.mpf_log, libmpc.mpc_log, "natural logarithm", libmpi.mpi_log)

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

floor = mpfunc('floor', libmpf.mpf_floor, libmpc.mpc_floor, "")
floor.__doc__ = """Computes the floor function of x. Note: returns a floating-point
number, not a Python int. If x is larger than the precision, it will be rounded,
not necessarily in the floor direction."""

ceil = mpfunc('ceil', libmpf.mpf_ceil, libmpc.mpc_ceil, "")
ceil.__doc__ = """Computes the ceiling function of x. Note: returns a floating-point
number, not a Python int. If x is larger than the precision, it will be rounded,
not necessarily in the floor direction."""

@funcwrapper
def nthroot(x, n):
    """principal n-th root"""
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
    """Returns the Euclidean distance sqrt(x*x + y*y). Both x and y
    must be real."""
    x = convert_lossless(x)
    y = convert_lossless(y)
    return make_mpf(libmpf.mpf_hypot(x._mpf_, y._mpf_, *prec_rounding))

def ldexp(x, n):
    """Calculate mpf(x) * 2**n efficiently. No rounding is performed."""
    x = convert_lossless(x)
    return make_mpf(libmpf.mpf_shift(x._mpf_, n))

def frexp(x):
    """Convert x to a scaled number y in the range [0.5, 1). Returns
    (y, n) such that x = y * 2**n. No rounding is performed."""
    x = convert_lossless(x)
    y, n = libmpf.mpf_frexp(x._mpf_)
    return make_mpf(y), n

def sign(x):
    """Return sign(x), defined as x/abs(x), or 0 for x = 0."""
    x = convert_lossless(x)
    if not x or isnan(x):
        return x
    if isinstance(x, mpf):
        return cmp(x, 0)
    return x / abs(x)

@extraprec(5)
def arg(x):
    """Returns the complex argument (phase) of x. The returned value is
    an mpf instance. The argument is here defined to satisfy
    -pi < arg(x) <= pi. On the negative real half-axis, it is taken to
    be +pi."""
    x = mpc(x)
    return atan2(x.imag, x.real)

def log(x, b=None):
    """Returns the base-b logarithm of x. If b is unspecified, return
    the natural (base-e) logarithm. log(x, b) is defined as
    log(x)/log(b). log(0) raises ValueError.

    The natural logarithm is real if x > 0 and complex if x < 0 or if x
    is complex. The principal branch of the complex logarithm is chosen,
    for which Im(log(x)) = -pi < arg(x) <= pi. """
    if b is None:
        return ln(x)
    wp = mp.prec + 20
    return ln(x, prec=wp) / ln(b, prec=wp)

def log10(x):
    """Base-10 logarithm. Equivalent to log(x,10)."""
    return log(x, 10)

def power(x, y):
    """Converts x and y to mpf or mpc and returns x**y = exp(y*log(x))."""
    return convert_lossless(x) ** convert_lossless(y)

def modf(x,y):
    """Converts x and y to mpf or mpc and returns x % y"""
    x = convert_lossless(x)
    y = convert_lossless(y)
    return x % y

def degrees(x):
    """Convert x given in radians to degrees"""
    return x / degree

def radians(x):
    """Convert x given in degrees to radians"""
    return x * degree

def atan2(y,x):
    """atan2(y, x) has the same magnitude as atan(y/x) but accounts for
    the signs of y and x. (Defined for real x and y only.)"""
    x = convert_lossless(x)
    y = convert_lossless(y)
    return make_mpf(libelefun.mpf_atan2(y._mpf_, x._mpf_, *prec_rounding))


cospi = mpfunc('cospi', gammazeta.mpf_cos_pi, gammazeta.mpc_cos_pi, 'computes cos(pi*x) accurately')
sinpi = mpfunc('sinpi', gammazeta.mpf_sin_pi, gammazeta.mpc_sin_pi, 'computes sin(pi*x) accurately')

zeta = mpfunc('zeta', gammazeta.mpf_zeta, gammazeta.mpc_zeta, 'Riemann zeta function')
gamma = mpfunc('gamma', gammazeta.mpf_gamma, gammazeta.mpc_gamma, "gamma function")
factorial = mpfunc('factorial', gammazeta.mpf_factorial, gammazeta.mpc_factorial, "factorial")
fac = factorial

def psi(m, z):
    """
    Gives the polygamma function of order m of z, psi^(m)(z). Special
    cases are the digamma function (psi0), trigamma function (psi1),
    tetragamma (psi2) and pentagamma (psi4) functions.

    The parameter m should be a nonnegative integer.
    """
    z = convert_lossless(z)
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

def bernoulli(n):
    """nth Bernoulli number, B_n"""
    return make_mpf(gammazeta.mpf_bernoulli(int(n), *prec_rounding))

bernfrac = gammazeta.bernfrac

stieltjes_cache = {}

def stieltjes(n):
    """Computes the nth Stieltjes constant."""
    n = int(n)
    if n == 0:
        return +euler
    if n < 0:
        raise ValueError("Stieltjes constants defined for n >= 0")
    if n in stieltjes_cache:
        prec, s = stieltjes_cache[n]
        if prec >= mp.prec:
            return +s
    from quadrature import quadgl
    def f(x):
        r = exp(pi*j*x)
        return (zeta(r+1) / r**n).real
    orig = mp.prec
    try:
        p = int(log(factorial(n), 2) + 35)
        mp.prec += p
        u = quadgl(f, [-1, 1])
        v = mpf(-1)**n * factorial(n) * u / 2
    finally:
        mp.prec = orig
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
    """
    Computes the product / quotient of gamma functions

        G(a_0) G(a_1) ... G(a_p)
        ------------------------
        G(b_0) G(b_1) ... G(a_q)

    with proper cancellation of poles (interpreting the expression as a
    limit). Returns +inf if the limit diverges.
    """
    a = [convert_lossless(x) for x in a]
    b = [convert_lossless(x) for x in b]
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

def binomial(n, k):
    """Binomial coefficient, C(n,k) = n!/(k!*(n-k)!)."""
    return gammaprod([n+1], [k+1, n-k+1])

def rf(x, n):
    """Rising factorial (Pochhammer symbol), x^(n)"""
    return gammaprod([x+n], [x])

def ff(x, n):
    """Falling factorial, x_(n)"""
    return gammaprod([x+1], [x-n+1])


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
    x = convert_lossless(x)
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
    a = (ar and _as_num(ar[0])) or convert_lossless(a)
    b = (br and _as_num(br[0])) or convert_lossless(b)
    c = (cr and _as_num(cr[0])) or convert_lossless(c)
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
    """
    Hypergeometric function pFq,

          [ a_1, a_2, ..., a_p |    ]
      pFq [                    |  z ]
          [ b_1, b_2, ..., b_q |    ]

    The parameter lists a_s and b_s may contain real or complex numbers.
    Exact rational parameters can be given as tuples (p, q).
    """
    p = len(a_s)
    q = len(b_s)
    z = convert_lossless(z)
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
    """Hypergeometric function 0F1. hyp0f1(a,z) is equivalent
    to hyper([], [a], z); see documentation for hyper() for more
    information."""
    return hyper([], [a], z)

def hyp1f1(a,b,z):
    """Hypergeometric function 1F1. hyp1f1(a,b,z) is equivalent
    to hyper([a], [b], z); see documentation for hyper() for more
    information."""
    return hyper([a], [b], z)

def hyp2f1(a,b,c,z):
    """Hypergeometric function 2F1. hyp2f1(a,b,c,z) is equivalent
    to hyper([a,b], [c], z); see documentation for hyper() for more
    information."""
    return hyper([a,b], [c], z)

@funcwrapper
def lower_gamma(a,z):
    """Lower incomplete gamma function gamma(a, z)"""
    # XXX: may need more precision
    return hyp1f1(1, 1+a, z) * z**a * exp(-z) / a

@funcwrapper
def upper_gamma(a,z):
    """Upper incomplete gamma function Gamma(a, z)"""
    return gamma(a) - lower_gamma(a, z)

erf = mpfunc("erf", libhyper.mpf_erf, libhyper.mpc_erf,
    "Error function, erf(z)")
erfc = mpfunc("erfc", libhyper.mpf_erfc, libhyper.mpc_erfc,
    "Complementary error function, erfc(z) = 1-erf(z)")

@funcwrapper
def erfi(z):
    """Imaginary error function, erfi(z)"""
    return (2/sqrt(pi)*z) * sum_hyp1f1_rat((1,2),(3,2), z**2)

@funcwrapper
def npdf(x, mu=0, sigma=1):
    """
    npdf(x, mu=0, sigma=1) -- probability density function of a
    normal distribution with mean value mu and variance sigma^2.
    """
    sigma = convert_lossless(sigma)
    return exp(-(x-mu)**2/(2*sigma**2)) / (sigma*sqrt(2*pi))

@funcwrapper
def ncdf(x, mu=0, sigma=1):
    """
    ncdf(x, mu=0, sigma=1) -- cumulative distribution function of
    a normal distribution with mean value mu and variance sigma^2.
    """
    a = (x-mu)/(sigma*sqrt(2))
    if a < 0:
        return erfc(-a)/2
    else:
        return (1+erf(a))/2

@funcwrapper
def ei(z):
    """Exponential integral, Ei(z)"""
    if z == inf:
        return z
    if z == -inf:
        return -mpf(0)
    v = z*hypsum([[1,1],[1,1]],[],[],[[2,1],[2,1]],[],[],z) + \
        (log(z)-log(1/z))/2 + euler
    if isinstance(z, mpf) and z < 0:
        return v.real
    return v

@funcwrapper
def li(z):
    """Logarithmic integral, li(z)"""
    if not z:
        return z
    if z == 1:
        return -inf
    return ei(log(z))

@funcwrapper
def ci(z):
    """Cosine integral, Ci(z)"""
    if z == inf:
        return 1/z
    if not z:
        return -inf
    z2 = -(z/2)**2
    return euler + log(z) + \
        z2*hypsum([[1,1],[1,1]],[],[],[[2,1],[2,1],[3,2]],[],[],z2)

@funcwrapper
def si(z):
    """Sine integral, Si(z)"""
    if z == inf:
        return pi/2
    if z == -inf:
        return -pi/2
    z2 = -(z/2)**2
    return z*hypsum([[1,2]],[],[],[[3,2],[3,2]],[],[],z2)

@funcwrapper
def chi(z):
    """Hyperbolic cosine integral, Chi(z)"""
    if not z:
        return -inf
    z2 = (z/2)**2
    return euler + log(z) + \
        z2*hypsum([[1,1],[1,1]],[],[],[[2,1],[2,1],[3,2]],[],[],z2)

@funcwrapper
def shi(z):
    """Hyperbolic sine integral, Shi(z)"""
    z2 = (z/2)**2
    return z*hypsum([[1,2]],[],[],[[3,2],[3,2]],[],[],z2)

@funcwrapper
def fresnels(z):
    """Fresnel integral S, S(z)"""
    if z == inf:
        return mpf(0.5)
    if z == -inf:
        return mpf(-0.5)
    return pi*z**3/6*hypsum([[3,4]],[],[],[[3,2],[7,4]],[],[],-pi**2*z**4/16)

@funcwrapper
def fresnelc(z):
    """Fresnel integral C, C(z)"""
    if z == inf:
        return mpf(0.5)
    if z == -inf:
        return mpf(-0.5)
    return z*hypsum([[1,4]],[],[],[[1,2],[5,4]],[],[],-pi**2*z**4/16)

@funcwrapper
def airyai(z):
    """Airy function, Ai(z)"""
    if z == inf:
        return 1/z
    if z == -inf:
        return mpf(0)
    z3 = z**3 / 9
    a = sum_hyp0f1_rat((2,3), z3) / (cbrt(9) * gamma(mpf(2)/3))
    b = z * sum_hyp0f1_rat((4,3), z3) / (cbrt(3) * gamma(mpf(1)/3))
    return a - b

@funcwrapper
def airybi(z):
    """Airy function, Bi(z)"""
    if z == inf:
        return z
    if z == -inf:
        return mpf(0)
    z3 = z**3 / 9
    rt = nthroot(3, 6)
    a = sum_hyp0f1_rat((2,3), z3) / (rt * gamma(mpf(2)/3))
    b = z * rt * sum_hyp0f1_rat((4,3), z3) / gamma(mpf(1)/3)
    return a + b

@funcwrapper
def ellipe(m):
    """Complete elliptic integral of the second kind, E(m). Note that
    the argument is the parameter m = k^2, not the modulus k."""
    if m == 1:
        return m
    return pi/2 * sum_hyp2f1_rat((1,2),(-1,2),(1,1), m)

@funcwrapper
def ellipk(m):
    """Complete elliptic integral of the first kind, K(m). Note that
    the argument is the parameter m = k^2, not the modulus k."""
    # Poor implementation:
    # return pi/2 * sum_hyp2f1_rat((1,2),(1,2),(1,1), m)
    if m == 1:
        return inf
    if isnan(m):
        return m
    if isinf(m):
        return 1/m
    s = sqrt(m)
    a = (1-s)/(1+s)
    v = pi/4*(1+a)/agm(1,a)
    if isinstance(m, mpf) and m < 1:
        return v.real
    return v

# TODO: for complex a, b handle the branch cut correctly
@funcwrapper
def agm(a, b=1):
    """Arithmetic-geometric mean of a and b. Can be called with
    a single argument, computing agm(a,1) = agm(1,a)."""
    if not a or not b:
        return a*b
    weps = eps * 16
    half = mpf(0.5)
    while abs(a-b) > weps:
        a, b = (a+b)*half, (a*b)**half
    return a

@funcwrapper
def jacobi(n, a, b, x):
    """Jacobi polynomial P_n^(a,b)(x)."""
    return binomial(n+a,n) * hyp2f1(-n,1+n+a+b,a+1,(1-x)/2)

@funcwrapper
def legendre(n, x):
    """Legendre polynomial P_n(x)."""
    if isint(n):
        n = int(n)
    if x == -1:
        # TODO: hyp2f1 should handle this
        if x == int(x):
            return (-1)**(n + (n>=0)) * mpf(-1)
        return inf
    return hyp2f1(-n,n+1,1,(1-x)/2)

@funcwrapper
def chebyt(n, x):
    """Chebyshev polynomial of the first kind T_n(x)."""
    return hyp2f1(-n,n,0.5,(1-x)/2)

@funcwrapper
def chebyu(n, x):
    """Chebyshev polynomial of the second kind U_n(x)."""
    return (n+1) * hyp2f1(-n, n+2, 1.5, (1-x)/2)

@funcwrapper
def jv(v, x):
    """Bessel function J_v(x)."""
    if isint(v):
        if isinstance(x, mpf):
            return make_mpf(libhyper.mpf_besseljn(int(v), x._mpf_, mp.prec))
        if isinstance(x, mpc):
            return make_mpc(libhyper.mpc_besseljn(int(v), x._mpc_, mp.prec))
    hx = x/2
    return hx**v * hyp0f1(v+1, -hx**2) / factorial(v)

jn = jv

def j0(x):
    """Bessel function J_0(x)."""
    return jv(0, x)

def j1(x):
    """Bessel function J_1(x)."""
    return jv(1, x)

@funcwrapper
def lambertw(z, k=0, approx=None):
    """
    lambertw(z,k) gives the kth branch of the Lambert W function W(z),
    defined as the kth solution of z = W(z)*exp(W(z)).

    lambertw(z) == lambertw(z, k=0) gives the principal branch
    value (0th branch solution), which is real for z > -1/e .

    The k = -1 branch is real for -1/e < z < 0. All branches except
    k = 0 have a logarithmic singularity at 0.

    The definition, implementation and choice of branches is based
    on Corless et al, "On the Lambert W function", Adv. Comp. Math. 5
    (1996) 329-359, available online here:
    http://www.apmaths.uwo.ca/~djeffrey/Offprints/W-adv-cm.pdf

    TODO: use a series expansion when extremely close to the branch point
    at -1/e and make sure that the proper branch is chosen there
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
        if z == inf: return z
        if z == -inf: return nan
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
