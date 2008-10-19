from mptypes import (mp, mpf, convert_lossless, inf,
   eps, nstr, make_mpf, AS_POINTS)

from functions import pi, exp, log, ldexp

from libmpf import mpf_neg

import math

def NEG(x):
    return make_mpf(mpf_neg(x._mpf_))

def transform(f, a, b):
    """
    Given an integrand f defined over the interval [a, b], return an
    equivalent integrand g defined on the standard interval [-1, 1].

    If a and b are finite, this is achived by means of a linear change
    of variables. If at least one point is infinite, the substitution
    t = 1/x is used.
    """
    a = convert_lossless(a)
    b = convert_lossless(b)
    if (a, b) == (-1, 1):
        return f
    one = mpf(1)
    half = mpf(0.5)
    # The transformation 1/x sends [1, inf] to [0, 1], which in turn
    # can be transformed to [-1, 1] the usual way. For a double
    # infinite interval, we simply evaluate the function symmetrically
    if (a, b) == (-inf, inf):
        # return transform(lambda x: (f(-1/x+1)+f(1/x-1))/x**2, 0, 1)
        def y(x):
            u = 2/(x+one)
            w = one - u
            return half * (f(w)+f(-w)) * u**2
        return y
    if a == -inf:
        # return transform(lambda x: f(-1/x+b+1)/x**2, 0, 1)
        b1 = b+1
        def y(x):
            u = 2/(x+one)
            return half * f(b1-u) * u**2
        return y
    if b == inf:
        # return transform(lambda x: f(1/x+a-1)/x**2, 0, 1)
        a1 = a-1
        def y(x):
            u = 2/(x+one)
            return half * f(a1+u) * u**2
        return y
    # Simple linear change of variables
    C = (b-a)/2
    D = (b+a)/2
    def g(x):
        return C * f(D + C*x)
    return g


class Quadrature:

    @classmethod
    def get_nodes(cls, prec, level, verbose=False):
        if (prec, level) in cls.cached_nodes:
            return cls.cached_nodes[prec, level]
        orig = mp.prec
        try:
            nodes = cls.calc_nodes(prec, level, verbose)
        finally:
            mp.prec = orig
        cls.cached_nodes[prec, level] = nodes
        return nodes

    @classmethod
    def guess_level(cls, prec):
        """
        Estimate a reasonable default maximum level for quadrature that
        doubles the accuracy when the level is incremented.
        (This is the case for tanh-sinh quadrature as implemented here.)

            50 bits -> 6
            100 bits -> 7
            500 bits -> 10
            3000 bits -> 12

        These numbers are based purely on a limited amount of
        experimentation and will sometimes be wrong.
        """
        # Expected level
        g = int(4 + max(0, log(prec/30.0, 2)))
        # Reasonable "worst case"
        g += 2
        return g

    @classmethod
    def estimate_error(cls, results, prec, epsilon):
        """
        Estimate error of the calculation at the present level by
        comparing it to the results from two previous levels. The
        algorithm is given by Borwein, Bailey & Girgensohn.
        """
        try:
            if results[-1] == results[-2] == results[-3]:
                return mpf(0)
            D1 = log(abs(results[-1]-results[-2]), 10)
            D2 = log(abs(results[-1]-results[-3]), 10)
        except ValueError:
            return epsilon
        D3 = -prec
        D4 = min(0, max(D1**2/D2, 2*D1, D3))
        return mpf(10) ** int(D4)

    @classmethod
    def summation(cls, f, points, prec, epsilon, max_level, verbose=False):
        """
        Main summation function
        """
        I = err = mpf(0)
        for i in xrange(len(points)-1):
            a, b = points[i], points[i+1]
            if a == b:
                continue
            g = transform(f, a, b)
            results = []
            for level in xrange(1, max_level+1):
                if verbose:
                    print "Integrating from %s to %s (level %s of %s)" % \
                        (nstr(a), nstr(b), level, max_level)
                results.append(cls.sum_next(prec, level, results, g, verbose))
                if level > 2:
                    err = cls.estimate_error(results, prec, epsilon)
                    if err <= epsilon:
                        break
                    if verbose:
                        print "Estimated error:", nstr(err)
            I += results[-1]
        if err > epsilon:
            if verbose:
                print "Failed to reach full accuracy. Estimated error:", nstr(err)
        return I, err

class TanhSinh(Quadrature):
    """
    This class implements "tanh-sinh" or "doubly exponential"
    quadrature. This quadrature rule is based on the Euler-Maclaurin
    integral formula. By performing a change of variables involving
    nested exponentials / hyperbolic functions (hence the name), the
    derivatives at the endpoints vanish rapidly.

    Since the error term in the Euler-Maclaurin formula depends
    on the derivatives at the endpoints, a simple step sum becomes
    extremely accurate. In practice, this means that doubling the
    number of evaluation points roughly doubles the number of accurate
    digits.

    Comparison to Gauss-Legendre:
      * Initial computation of nodes is usually faster
      * Handles endpoint singularities better
      * Handles infinite integration intervals better
      * Is slower for smooth integrands once nodes have been computed

    The implementation of the tanh-sinh algorithm is based on the
    description given in Borwein, Bailey & Girgensohn, "Experimentation
    in Mathematics - Computational Paths to Discovery", A K Peters,
    2003, pages 312-313.

    A few improvements have been made:
      * A more efficient scheme is used to compute nodes
      * The nodes are computed successively instead of all at once

    Various documents describing the algorithm are available online, e.g.:

      * http://crd.lbl.gov/~dhbailey/dhbpapers/dhb-tanh-sinh.pdf
      * http://users.cs.dal.ca/~jborwein/tanh-sinh.pdf
    """

    cached_nodes = {}

    @classmethod
    def sum_next(cls, prec, level, previous, f, verbose=False):
        h = mpf(2)**(-level)
        # Abscissas overlap, so reusing saves half of the time
        if previous:
            S = previous[-1]/(h*2)
        else:
            S = mpf(0)
        for x, w in cls.get_nodes(prec, level, verbose=False):
            S += w*(f(NEG(x)) + f(x))
        return h*S

    @classmethod
    def calc_nodes(cls, prec, level, verbose=False):
        """
        The abscissas and weights for tanh-sinh quadrature are given by

            x[k] = tanh(pi/2 * sinh(t))
            w[k] = pi/2 * cosh(t) / cosh(pi/2 sinh(t))**2

        Here t varies uniformly with k: t0, t0+h, t0+2*h, ...

        The list of nodes is actually infinite, but the weights
        die off so rapidly that only a few are needed.
        """
        nodes = []

        extra = 20
        mp.prec += extra
        eps = ldexp(1, -prec-10)
        pi4 = pi/4

        # For simplicity, we work in steps h = 1/2^n, with the first point
        # offset so that we can reuse the sum from the previous level

        # We define level 1 to include the "level 0" steps, including
        # the point x = 0. (It doesn't work well otherwise; not sure why.)
        t0 = ldexp(1, -level)
        if level == 1:
            nodes.append((mpf(0), pi4))
            h = t0
        else:
            h = t0*2

        # Since h is fixed, we can compute the next exponential
        # by simply multiplying by exp(h)
        expt0 = exp(t0)
        a = pi4 * expt0
        b = pi4 / expt0
        udelta = exp(h)
        urdelta = 1/udelta

        for k in xrange(0, 20*2**level+1):
            # Reference implementation:
            # t = t0 + k*h
            # x = tanh(pi/2 * sinh(t))
            # w = pi/2 * cosh(t) / cosh(pi/2 * sinh(t))**2

            # Fast implementation. Note that c = exp(pi/2 * sinh(t))
            c = exp(a-b)
            d = 1/c
            co = (c+d)/2
            si = (c-d)/2
            x = si / co
            w = (a+b) / co**2
            diff = abs(x-1)
            if diff <= eps:
                break

            nodes.append((x, w))
            a *= udelta
            b *= urdelta

            if verbose and k % 300 == 150:
                # Note: the number displayed is rather arbitrary. Should
                # figure out how to print something that looks more like a
                # percentage
                print "Calculating nodes:", nstr(-log(diff, 10) / prec)

        mp.prec -= extra
        return nodes


class GaussLegendre(Quadrature):
    """
    This class implements Gauss-Legendre quadrature, which is
    exceptionally efficient for polynomials and polynomial-like (i.e.
    very smooth) integrands.

    The abscissas and weights are given by roots and values of
    Legendre polynomials, which are the orthogonal polynomials
    on [-1, 1] with respect to the unit weight.

    In this implementation, we take the "level" m of the quadrature
    to denote a Gauss-Legendre rule of degree 3*2^m (following Borwein,
    Bailey & Girgensohn). This ensures quadratic convergence.

    Comparison to tanh-sinh quadrature:
      * Is faster for smooth integrands once nodes have been computed
      * Initial computation of nodes is usually slower
      * Handles endpoint singularities worse
      * Handles infinite integration intervals worse

    """

    cached_nodes = {}

    @classmethod
    def calc_nodes(cls, prec, level, verbose=False):
        # It is important that the epsilon is set lower than the
        # "real" epsilon
        epsilon = ldexp(1, -prec-8)
        # Fairly high precision might be required for accurate
        # evaluation of the roots
        orig = mp.prec
        mp.prec = int(prec*1.5)
        nodes = []
        n = 3*2**(level-1)
        upto = n//2 + 1
        for j in xrange(1, upto):
            # Asymptotic formula for the roots
            r = mpf(math.cos(math.pi*(j-0.25)/(n+0.5)))
            # Newton iteration
            while 1:
                t1, t2 = 1, 0
                # Evaluates the Legendre polynomial using its defining
                # recurrence relation
                for j1 in xrange(1,n+1):
                    t3, t2, t1 = t2, t1, ((2*j1-1)*r*t1 - (j1-1)*t2)/j1
                t4 = n*(r*t1- t2)/(r**2-1)
                t5 = r
                a = t1/t4
                r = r - a
                if abs(a) < epsilon:
                    break
            x = r
            w = 2/((1-r**2)*t4**2)
            if verbose  and j % 30 == 15:
                print "Computing nodes (%i of %i)" % (j, upto)
            nodes.append((x, w))
        mp.prec = orig
        return nodes

    @classmethod
    def sum_next(cls, prec, level, previous, f, verbose=False):
        s = mpf(0)
        for x, w in cls.get_nodes(prec, level, verbose):
            s += w * (f(NEG(x)) + f(x))
        return s


def quad(f, *points, **kwargs):
    """
    Computes a single, double or triple integral over a given
    interval, rectangle, or box.

        quad(f(x), [x1, x2])
        quad(f(x,y), [x1, x2], [y1, y2])
        quad(f(x,y,z), [x1, x2], [y1, y2], [z1, z2])

    By default, tanh-sinh quadrature is used. A custom method
    can be specified via the 'method' keyword argument. Basic examples:

        >>> from mpmath import *
        >>> print quad(lambda x: exp(-x**2), [0, inf])
        0.886226925452758
        >>> f = lambda x, y: exp(x*sin(y))
        >>> print quad(f, [0, 1], [0, pi])
        4.47046663466179

    The functions quadgl(...) and quadts(...) act as shortcuts for
    quad(..., method='gauss-legendre') and quad(..., method='tanh-sinh').

    An interval may contain more than two points. In this case, the
    integration is split into subintervals, between each pair of
    consecutive points. This is useful for dealing with
    mid-interval discontinuities, or integrating over large
    intervals where the function is irregular or oscillates:

        >>> print quadgl(lambda x: abs(sin(x)), [0, pi, 2*pi])
        4.0
        >>> print quadgl(sin, arange(0, 1000+1, 10))
        0.437620923709297
        >>> print cos(0) - cos(1000)
        0.437620923709297

    Additional keyword options:

        verbose=True -- print details about progress

        error=True   -- return (value, err) where err is the estimated
                        error

    """
    rule = kwargs.get('method', TanhSinh)
    if type(rule) is str:
        rule = {'tanh-sinh':TanhSinh, 'gauss-legendre':GaussLegendre}[rule]
    verbose = kwargs.get('verbose')
    dim = len(points)
    orig = prec = mp.prec
    epsilon = eps/8
    m = kwargs.get('maxlevel') or rule.guess_level(prec)
    points = [AS_POINTS(p) for p in points]
    try:
        mp.prec += 20
        if dim == 1:
            v, err = rule.summation(f, points[0], prec, epsilon, m, verbose)
        elif dim == 2:
            v, err = rule.summation(lambda x: \
                    rule.summation(lambda y: f(x,y), \
                    points[1], prec, epsilon, m)[0],
                points[0], prec, epsilon, m, verbose)
        elif dim == 3:
            v, err = rule.summation(lambda x: \
                    rule.summation(lambda y: \
                        rule.summation(lambda z: f(x,y,z), \
                        points[2], prec, epsilon, m)[0],
                    points[1], prec, epsilon, m)[0],
                points[0], prec, epsilon, m, verbose)
        else:
            raise NotImplementedError("quadrature must have dim 1, 2 or 3")
    finally:
        mp.prec = orig
    if kwargs.get("error"):
        return +v, err
    return +v

def quadts(*args, **kwargs):
    """
    Performs tanh-sinh quadrature. The call

        quadts(func, *points, ...)

    is simply a shortcut for:

        quad(func, *points, ..., method=TanhSinh)

    For example, a single integral and a double integral:

        quadts(lambda x: exp(cos(x)), [0, 1])
        quadts(lambda x, y: exp(cos(x+y)), [0, 1], [0, 1])

    See the documentation for quad for information about how points
    arguments and keyword arguments are parsed.

    See documentation for TanhSinh for algorithmic information about
    tanh-sinh quadrature.
    """
    kwargs['method'] = TanhSinh
    return quad(*args, **kwargs)

def quadgl(*args, **kwargs):
    """
    Performs Gauss-Legendre quadrature. The call

        quadgl(func, *points, ...)

    is simply a shortcut for:

        quad(func, *points, ..., method=TanhSinh)

    For example, a single integral and a double integral:

        quadgl(lambda x: exp(cos(x)), [0, 1])
        quadgl(lambda x, y: exp(cos(x+y)), [0, 1], [0, 1])

    See the documentation for quad for information about how points
    arguments and keyword arguments are parsed.

    See documentation for TanhSinh for algorithmic information about
    tanh-sinh quadrature.
    """
    kwargs['method'] = GaussLegendre
    return quad(*args, **kwargs)

