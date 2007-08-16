"""
Routines for arbitrary-precision numerical quadrature
"""

from sympy import oo

#----------------------------------------------------------------------------#
#                        General classes / utilities                         #
#----------------------------------------------------------------------------#

from float_ import Float
from constants import *
from functions import *
from evalf_ import polyfunc

import math


class Quadrature:
    """
    This class provides a standard interface for quadrature rules. A
    rule is applied through a function call (e.g. to integrate sin(x)
    from x=3 to x=8, SomeQuadRule(<options>)(lambda x: sin(x), 3, 8))

    Each Quadrature subclass should define an '_eval' method that
    returns a tuple containing two values: (y, err) where y is the
    calculated value of the integral, and err is an estimate of
    the approximation error (or None if the error cannot be estimated
    from the rule).

    By default, the quadrature rule is assumed to be defined on
    [-1, 1] or (-1, 1). __call__ automatically transforms other
    (finite) intervals to this default interval. This behavior can
    be modified by overriding the 'transform' method.
    """

    def transform(self, f, a, b):
        """Given an integrand function f and boundaries a, b,
        return an equivalent integrand defined on [-1, 1]"""
        a = Float(a)
        b = Float(b)
        if (a, b) == (-1, 1):
            return f
        else:
            C = (b-a)/2
            D = (b+a)/2
            def g(x):
                return C * f(D + C*x)
            return g

    def __call__(self, f, a, b, extraprec=3, verbose=False):
        Float.store()
        Float.extraprec(extraprec)
        f = self.transform(f, a, b)
        try:
            s, err = self._eval(f, verbose)
        except Exception, e:
            Float.revert()
            raise e
        Float.revert()
        return +s, err

    def adaptive(self, f, a, b, eps, steps=0, maxsteps=5000, verbose=False):
        """Apply rule adaptively (must support error estimation)"""
        s, err = self(f, a, b, verbose=verbose)
        if err <= eps or steps >= maxsteps:
            return s, err, steps+1
        if verbose:
            print steps, a, b
        mid = (a+b)/2
        s1, e1, steps = self.adaptive(f, a, mid, eps, steps+1, maxsteps, verbose)
        s2, e2, steps = self.adaptive(f, mid, b, eps, steps, maxsteps, verbose)
        return s1 + s2, e1 + e2, steps


#----------------------------------------------------------------------------#
#                       Trivial quadrature rules                             #
#----------------------------------------------------------------------------#

# Mainly for demonstration purposes...

class Midpoint(Quadrature):
    def _eval(self, f, verbose=False):
        return 2*f(Float(0)), None

class Trapezoid(Quadrature):
    def _eval(self, f, verbose=False):
        return f(Float(-1)) + f(Float(1)), None


#----------------------------------------------------------------------------#
#                          Gaussian quadrature                               #
#----------------------------------------------------------------------------#

def _iterfixpoint(h, x, eps):
    while 1:
        new = h(x)
        if x.ae(new, eps):
            return new
        x = new

def _gaussnode(pf, k, n):
    import math
    x = Float(-math.cos(math.pi*(k+1-0.25)/(n+0.5)))
    def h(r):
        p, q = pf(r)
        return r - p/q
    eps = Float((1, -Float._prec//2))
    x = _iterfixpoint(h, x, eps)
    w = 2 / (1-x**2) / (pf(x)[1])**2
    return x, w

_gausscache = {}


class GaussLegendre(Quadrature):
    """
    Gauss-Legendre quadrature is highly efficient for smooth integrands
    on finite intervals.
    """

    def __init__(self, n=None, verbose=False):
        self.n = n
        prec = dps = Float.getdps()
        self.prec = prec

        # Reuse old nodes
        if n in _gausscache:
            cacheddps, xs, ws = _gausscache[n]
            if cacheddps >= dps:
                self.x = [x for x in xs]
                self.w = [w for w in ws]
                return

        if verbose:
            print ("calculating nodes for degree-%i Gauss-Legendre "
                "quadrature..." % n)

        Float.store()
        Float.setdps(2*prec + 5)

        self.x = [None] * n
        self.w = [None] * n

        from sympy.specfun import legendre
        pf = polyfunc(legendre(n, 'x'), True)

        for k in xrange(n//2 + 1):
            if verbose and k % 4 == 0:
                print "  node", k, "of", n//2
            x, w = _gaussnode(pf, k, n)
            self.x[k] = x
            self.x[n-k-1] = -x
            self.w[k] = self.w[n-k-1] = w

        _gausscache[n] = (dps, self.x, self.w)

        Float.revert()

    def _eval(self, f, verbose=False):
        s = Float(0)
        for i in xrange(self.n):
            s += self.w[i] * f(self.x[i])
        return s, None


class CompositeGaussLegendre(Quadrature):
    """
    Gauss-Legendre quadrature with error estimation. Two Gauss-Legendre
    rules of different degree are used, and the error is estimated as
    the difference between the two results.
    """

    def __init__(self, deg1, deg2, verbose=False):
        self.rule1 = GaussLegendre(deg1, verbose)
        self.rule2 = GaussLegendre(deg2, verbose)

    def _eval(self, f, verbose=False):
        s = self.rule1._eval(f, verbose)[0]
        t = self.rule2._eval(f, verbose)[0]
        return s, abs(s-t)


#----------------------------------------------------------------------------#
#                            Tanh-sinh quadrature                            #
#----------------------------------------------------------------------------#

def _tsnode(k, h):
    # x[k] = tanh(pi/2 * sinh(k*h))
    # w[k] = pi/2 * cosh(k*h) / cosh(pi/2 sinh(k*h))**2
    t = Float(k) * h
    a = exp(t)
    pi = pi_float()
    ar = Float(1) / a
    sinht = (a - ar) / 2
    cosht = (a + ar) / 2
    b = exp((pi * sinht) / 2)
    br = Float(1) / b
    sinhb = (b - br) / 2
    coshb = (b + br) / 2
    return sinhb/coshb, (pi/2)*cosht/coshb**2

_tscache = {}

class TanhSinh(Quadrature):
    """
    The tanh-sinh quadrature rule (also known as the doubly-exponential
    quadrature rule) is well suited for integrals with singularities
    and/or extremely high precision levels (tens or hundreds of digits).

    The level m is a small integer between, say, 5 for low precision
    levels and 10 for extremely high precision (it must be set higher
    if the integrand is ill-behaved).

    This implementation of the tanh-sinh algorithm is based on the
    description given in Borwein, Bailey & Girgensohn, "Experimentation
    in Mathematics - Computational Paths to Discovery", A K Peters,
    2003, pages 312-313. It also described in various documents
    available online.
    """

    def __init__(self, eps, m, verbose=False):
        # Compute abscissas and weights

        prec = Float.getdps()

        self.prec = prec
        self.eps = eps
        self.m = m
        self.h = h = Float(1) / 2**m
        self.x = []
        self.w = []

        if (eps, m) in _tscache:
            self.x, self.w = _tscache[(eps, m)]
            return

        if verbose:
            print ("calculating nodes for tanh-sinh quadrature with "
                "epsilon %s and degree %i" % (eps, m))

        Float.store()
        Float.setdps(prec + 10)

        for k in xrange(20 * 2**m + 1):
            x, w = _tsnode(k, h)
            self.x.append(x)
            self.w.append(w)
            diff = abs(self.x[-1] - Float(1))
            if diff <= eps:
                break
            if verbose and m > 6 and k % 300 == 299:
                Float.store(); Float.setdps(5)
                print "  progress", -log(diff, 10) / prec
                Float.revert()

        _tscache[(eps, m)] = self.x, self.w

        Float.revert()

    def _estimate_error(self, res):
        try:
            D1 = log(abs(res[-1]-res[-2]), 10)
            D2 = log(abs(res[-1]-res[-3]), 10)
        except ValueError:
            return self.eps
        D3 = -self.prec
        D4 = min(0, max(D1**2/D2, 2*D1, D3))
        return Float('0.1') ** -int(D4)

    def _eval(self, f, verbose=False):
        S = Float(0)
        h = Float(1)
        m = self.m
        res = []
        for k in xrange(1, m+1):
            h = h / 2
            for i in xrange(0, len(self.x), 2**(m-k)):
                if i % (2**(m-k+1)) != 0 or k == 1:
                    if i == 0:
                        S = S + self.w[0]*f(Float(0))
                    else:
                        S = S + (self.w[i])*(f(-self.x[i]) + f(self.x[i]))
            res.append(h*S)
            if k > 2:
                err = self._estimate_error(res)
                if err <= self.eps:
                    break
        return +res[-1], self._estimate_error(res)


class AdaptiveTanhSinh(Quadrature):
    """
    Adaptive tanh-sinh quadrature algorithm.
    """

    def __init__(self, initial):
        self.initial = initial

    def adaptive(self, f, a, b, eps, steps=0, maxsteps=5000, verbose=False):
        prec = Float.getprec()
        for m in xrange(self.initial, 50):
            if verbose:
                print "using tanh-sinh rule with m =", m
            ts = TanhSinh(eps, m, verbose=verbose)
            s, err = ts(f, a, b, verbose=verbose)
            steps = 2*len(ts.x)
            if err <= eps or steps >= maxsteps:
                return s, err, steps


#----------------------------------------------------------------------------#
#                          User-friendly functions                           #
#----------------------------------------------------------------------------#


def nintegrate(f, a, b, method=0, maxsteps=5000, verbose=False):
    """
    Basic usage
    ===========

    nintegrate(f, a, b) numerically evaluates the integral

           - b
          |    f(x) dx.
         -   a


    The integrand f should be a callable that accepts a Float as input
    and outputs a Float or ComplexFloat; a and b should be Floats,
    numbers that can be converted to Floats, or +/- infinity.

    A simple example:

        >>> Float.setdps(15)
        >>> print nintegrate(lambda x: 2*x**2 - 4*x, 2, 3)
        2.66666666666667

    Calculating the area of a unit circle, with 30 digits:

        >>> Float.setdps(30)
        >>> print nintegrate(lambda x: 4*sqrt(1-x**2), 0, 1)
        3.14159265358979323846264338328

    The integration interval can be infinite or semi-infinite:
    
        >>> Float.setdps(15)
        >>> print nintegrate(lambda x: exp(-x)*sin(x), 0, oo)
        0.5

    Integration methods and accuracy
    ================================

    Nintegrate attempts to obtain a value that is fully accurate within
    the current working precision (i.e., correct to 15 decimal places
    at the default precision level). If nintegrate fails to reach full
    accuracy after a certain number of steps, it prints a warning
    message.

    This message signifies either that the integral is either divergent
    or, if convergent, ill-behaved. It may still be possible to
    evaluate an ill-behaved integral by increasing the 'maxsteps'
    setting, changing the integration method, and/or manually
    transforming the integrand.

    Nintegrate currently supports the following integration methods:

        method = 0  :  Gaussian quadrature (default)
        method = 1  :  tanh-sinh quadrature

    Gaussian quadrature is generally very efficient if the integration
    interval is finite and the integrand is smooth on the entire range
    (including the endpoints). It may fail if the integrand has many
    discontinuities, is highly oscillatory, or possesses integrable
    singularities.

    The tanh-sinh algorithm is often better if the integration interval
    is infinite or if singularities are present at the endpoints;
    especially at very high precision levels. It does not perform well
    if there are singularities between the endpoints or the integrand
    is bumpy or oscillatory.

    It may help to manually transform the integrand, e.g. changing
    variables to remove singularities or breaking up the integration
    interval so that singularities appear only at the endpoints.

    The 'verbose' flag can be set to track the computation's progress.
    """
    dps = Float.getdps()
    Float.store()
    Float.setdps(dps + 3)
    prec = Float.getprec()

    # Transform infinite or semi-infinite interval to a finite interval
    if a == -oo or b == oo:
        g = f
        if a == -oo and b == oo:
            def f(x):
                # make adaptive quadrature work from the left
                x = 1 - x
                return g(x) + g(Float(1)/x)/x**2 + g(-x) + g(Float(-1)/x)/(-x)**2
        elif b == oo:
            aa = Float(a)
            def f(x):
                x = 1 - x
                return g(x + aa) + g(Float(1)/x + aa)/x**2
        elif a == -oo:
            bb = Float(b)
            def f(x):
                x = 1 - x
                return g(-x + bb) + g(Float(-1)/x + bb)/(-x)**2
        a, b = Float(0), Float(1)
    else:
        a, b = Float(a), Float(b)

    eps = Float((1, -prec+4))

    if method == 0:
        degree = int(5 + dps**0.8)
        rule = CompositeGaussLegendre(degree, degree//2, verbose)
    elif method == 1:
        rule = AdaptiveTanhSinh(initial = int(3 + max(0, math.log(prec/30.0, 2))))

    else:
        Float.revert()
        raise ValueError("unknown method")

    s, err, steps = rule.adaptive(f, Float(a), Float(b), eps, steps=0, maxsteps=maxsteps, verbose=verbose)

    Float.revert()

    if not err.ae(0):
        Float.store()
        Float.setdps(1)
        print ("Warning: failed to reach full accuracy. "
            "Estimated magnitude of error:", str(err))
        Float.revert()

    return +s
