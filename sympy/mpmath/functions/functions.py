class SpecialFunctions(object):
    """
    This class implements special functions using high-level code.

    Elementary and some other functions (e.g. gamma function, basecase
    hypergeometric series) are assumed to be predefined by the context as
    "builtins" or "low-level" functions.
    """
    defined_functions = {}

    # The series for the Jacobi theta functions converge for |q| < 1;
    # in the current implementation they throw a ValueError for
    # abs(q) > THETA_Q_LIM
    THETA_Q_LIM = 1 - 10**-7

    def __init__(self):
        cls = self.__class__
        for name in cls.defined_functions:
            f, wrap = cls.defined_functions[name]
            cls._wrap_specfun(name, f, wrap)

        self.mpq_1 = self._mpq((1,1))
        self.mpq_0 = self._mpq((0,1))
        self.mpq_1_2 = self._mpq((1,2))
        self.mpq_3_2 = self._mpq((3,2))
        self.mpq_1_4 = self._mpq((1,4))
        self.mpq_1_16 = self._mpq((1,16))
        self.mpq_3_16 = self._mpq((3,16))
        self.mpq_5_2 = self._mpq((5,2))
        self.mpq_3_4 = self._mpq((3,4))
        self.mpq_7_4 = self._mpq((7,4))
        self.mpq_5_4 = self._mpq((5,4))

        self._aliases.update({
            'phase' : 'arg',
            'conjugate' : 'conj',
            'nthroot' : 'root',
            'polygamma' : 'psi',
            'hurwitz' : 'zeta',
            #'digamma' : 'psi0',
            #'trigamma' : 'psi1',
            #'tetragamma' : 'psi2',
            #'pentagamma' : 'psi3',
            'fibonacci' : 'fib',
            'factorial' : 'fac',
        })

    # Default -- do nothing
    @classmethod
    def _wrap_specfun(cls, name, f, wrap):
        setattr(cls, name, f)

    # Optional fast versions of common functions in common cases.
    # If not overridden, default (generic hypergeometric series)
    # implementations will be used
    def _besselj(ctx, n, z): raise NotImplementedError
    def _erf(ctx, z): raise NotImplementedError
    def _erfc(ctx, z): raise NotImplementedError
    def _gamma_upper_int(ctx, z, a): raise NotImplementedError
    def _expint_int(ctx, n, z): raise NotImplementedError
    def _zeta(ctx, s): raise NotImplementedError
    def _zetasum_fast(ctx, s, a, n, derivatives, reflect): raise NotImplementedError
    def _ei(ctx, z): raise NotImplementedError
    def _e1(ctx, z): raise NotImplementedError
    def _ci(ctx, z): raise NotImplementedError
    def _si(ctx, z): raise NotImplementedError
    def _altzeta(ctx, s): raise NotImplementedError

def defun_wrapped(f):
    SpecialFunctions.defined_functions[f.__name__] = f, True

def defun(f):
    SpecialFunctions.defined_functions[f.__name__] = f, False

def defun_static(f):
    setattr(SpecialFunctions, f.__name__, f)

@defun_wrapped
def cot(ctx, z): return ctx.one / ctx.tan(z)

@defun_wrapped
def sec(ctx, z): return ctx.one / ctx.cos(z)

@defun_wrapped
def csc(ctx, z): return ctx.one / ctx.sin(z)

@defun_wrapped
def coth(ctx, z): return ctx.one / ctx.tanh(z)

@defun_wrapped
def sech(ctx, z): return ctx.one / ctx.cosh(z)

@defun_wrapped
def csch(ctx, z): return ctx.one / ctx.sinh(z)

@defun_wrapped
def acot(ctx, z): return ctx.atan(ctx.one / z)

@defun_wrapped
def asec(ctx, z): return ctx.acos(ctx.one / z)

@defun_wrapped
def acsc(ctx, z): return ctx.asin(ctx.one / z)

@defun_wrapped
def acoth(ctx, z): return ctx.atanh(ctx.one / z)

@defun_wrapped
def asech(ctx, z): return ctx.acosh(ctx.one / z)

@defun_wrapped
def acsch(ctx, z): return ctx.asinh(ctx.one / z)

@defun
def sign(ctx, x):
    x = ctx.convert(x)
    if not x or ctx.isnan(x):
        return x
    if ctx._is_real_type(x):
        return ctx.mpf(cmp(x, 0))
    return x / abs(x)

@defun
def agm(ctx, a, b=1):
    if b == 1:
        return ctx.agm1(a)
    a = ctx.convert(a)
    b = ctx.convert(b)
    return ctx._agm(a, b)

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

# TODO: tests; improve implementation
@defun_wrapped
def expm1(ctx, x):
    if not x:
        return ctx.zero
    # exp(x) - 1 ~ x
    if ctx.mag(x) < -ctx.prec:
        return x + 0.5*x**2
    # TODO: accurately eval the smaller of the real/imag parts
    return ctx.sum_accurately(lambda: iter([ctx.exp(x),-1]),1)

@defun_wrapped
def powm1(ctx, x, y):
    mag = ctx.mag
    one = ctx.one
    w = x**y - one
    M = mag(w)
    # Only moderate cancellation
    if M > -8:
        return w
    # Check for the only possible exact cases
    if not w:
        if (not y) or (x in (1, -1, 1j, -1j) and ctx.isint(y)):
            return w
    x1 = x - one
    magy = mag(y)
    lnx = ctx.ln(x)
    # Small y: x^y - 1 ~ log(x)*y + O(log(x)^2 * y^2)
    if magy + mag(lnx) < -ctx.prec:
        return lnx*y + (lnx*y)**2/2
    # TODO: accurately eval the smaller of the real/imag part
    return ctx.sum_accurately(lambda: iter([x**y, -1]), 1)

@defun
def _rootof1(ctx, k, n):
    k = int(k)
    n = int(n)
    k %= n
    if not k:
        return ctx.one
    elif 2*k == n:
        return -ctx.one
    elif 4*k == n:
        return ctx.j
    elif 4*k == 3*n:
        return -ctx.j
    return ctx.expjpi(2*ctx.mpf(k)/n)

@defun
def root(ctx, x, n, k=0):
    n = int(n)
    x = ctx.convert(x)
    if k:
        # Special case: there is an exact real root
        if (n & 1 and 2*k == n-1) and (not ctx.im(x)) and (ctx.re(x) < 0):
            return -ctx.root(-x, n)
        # Multiply by root of unity
        prec = ctx.prec
        try:
            ctx.prec += 10
            v = ctx.root(x, n, 0) * ctx._rootof1(k, n)
        finally:
            ctx.prec = prec
        return +v
    return ctx._nthroot(x, n)

@defun
def unitroots(ctx, n, primitive=False):
    gcd = ctx._gcd
    prec = ctx.prec
    try:
        ctx.prec += 10
        if primitive:
            v = [ctx._rootof1(k,n) for k in range(n) if gcd(k,n) == 1]
        else:
            # TODO: this can be done *much* faster
            v = [ctx._rootof1(k,n) for k in range(n)]
    finally:
        ctx.prec = prec
    return [+x for x in v]

@defun
def arg(ctx, x):
    x = ctx.convert(x)
    re = ctx._re(x)
    im = ctx._im(x)
    return ctx.atan2(im, re)

@defun
def fabs(ctx, x):
    return abs(ctx.convert(x))

@defun
def re(ctx, x):
    x = ctx.convert(x)
    if hasattr(x, "real"):    # py2.5 doesn't have .real/.imag for all numbers
        return x.real
    return x

@defun
def im(ctx, x):
    x = ctx.convert(x)
    if hasattr(x, "imag"):    # py2.5 doesn't have .real/.imag for all numbers
        return x.imag
    return ctx.zero

@defun
def conj(ctx, x):
    return ctx.convert(x).conjugate()

@defun
def polar(ctx, z):
    return (ctx.fabs(z), ctx.arg(z))

@defun_wrapped
def rect(ctx, r, phi):
    return r * ctx.mpc(*ctx.cos_sin(phi))

@defun
def log(ctx, x, b=None):
    if b is None:
        return ctx.ln(x)
    wp = ctx.prec + 20
    return ctx.ln(x, prec=wp) / ctx.ln(b, prec=wp)

@defun
def log10(ctx, x):
    return ctx.log(x, 10)

@defun
def modf(ctx, x, y):
    return ctx.convert(x) % ctx.convert(y)

@defun
def degrees(ctx, x):
    return x / ctx.degree

@defun
def radians(ctx, x):
    return x * ctx.degree

@defun_wrapped
def lambertw(ctx, z, k=0):
    k = int(k)
    if ctx.isnan(z):
        return z
    ctx.prec += 20
    mag = ctx.mag(z)
    # Start from fp approximation
    if ctx is ctx._mp and abs(mag) < 900 and abs(k) < 10000 and \
        abs(z+0.36787944117144) > 0.01:
        w = ctx._fp.lambertw(z, k)
    else:
        absz = abs(z)
        # We must be extremely careful near the singularities at -1/e and 0
        u = ctx.exp(-1)
        if absz <= u:
            if not z:
                # w(0,0) = 0; for all other branches we hit the pole
                if not k:
                    return z
                return ctx.ninf
            if not k:
                w = z
            # For small real z < 0, the -1 branch aves roughly like log(-z)
            elif k == -1 and not ctx.im(z) and ctx.re(z) < 0:
                w = ctx.ln(-z)
            # Use a simple asymptotic approximation.
            else:
                w = ctx.ln(z)
                # The branches are roughly logarithmic. This approximation
                # gets better for large |k|; need to check that this always
                # works for k ~= -1, 0, 1.
                if k: w += k * 2*ctx.pi*ctx.j
        elif k == 0 and ctx.im(z) and absz <= 0.7:
            # Both the W(z) ~= z and W(z) ~= ln(z) approximations break
            # down around z ~= -0.5 (converging to the wrong branch), so patch
            # with a constant approximation (adjusted for sign)
            if abs(z+0.5) < 0.1:
                if ctx.im(z) > 0:
                    w = ctx.mpc(0.7+0.7j)
                else:
                    w = ctx.mpc(0.7-0.7j)
            else:
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
            if k:
                w += k * 2*ctx.pi*ctx.j
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
    ctx.warn("Lambert W iteration failed to converge for %s" % z)
    return wn

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
    return ctx.sum_accurately(_terms, check_step=4)

@defun_wrapped
def polyexp(ctx, s, z):
    if ctx.isinf(z) or ctx.isinf(s) or ctx.isnan(z) or ctx.isnan(s):
        return z**s
    if z == 0: return z*s
    if s == 0: return ctx.expm1(z)
    if s == 1: return ctx.exp(z)*z
    if s == 2: return ctx.exp(z)*z*(z+1)
    return _polyexp(ctx, s, z)

@defun_wrapped
def cyclotomic(ctx, n, z):
    n = int(n)
    assert n >= 0
    p = ctx.one
    if n == 0:
        return p
    if n == 1:
        return z - p
    if n == 2:
        return z + p
    # Use divisor product representation. Unfortunately, this sometimes
    # includes singularities for roots of unity, which we have to cancel out.
    # Matching zeros/poles pairwise, we have (1-z^a)/(1-z^b) ~ a/b + O(z-1).
    a_prod = 1
    b_prod = 1
    num_zeros = 0
    num_poles = 0
    for d in range(1,n+1):
        if not n % d:
            w = ctx.moebius(n//d)
            # Use powm1 because it is important that we get 0 only
            # if it really is exactly 0
            b = -ctx.powm1(z, d)
            if b:
                p *= b**w
            else:
                if w == 1:
                    a_prod *= d
                    num_zeros += 1
                elif w == -1:
                    b_prod *= d
                    num_poles += 1
    #print n, num_zeros, num_poles
    if num_zeros:
        if num_zeros > num_poles:
            p *= 0
        else:
            p *= a_prod
            p /= b_prod
    return p
