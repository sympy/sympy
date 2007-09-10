from sympy import Polynomial, Symbol
from float_ import Float, ComplexFloat
from evalf_ import polyfunc
from utils_ import bitcount
from functions import *

import random
rng = random.Random()

class ConvergenceError(Exception):
    pass


def polyroots(poly, maxsteps=20):
    """
    Numerically locate all (complex) roots of a given polynomial 'poly'
    using the Durand-Kerner method (http://en.wikipedia.org/wiki/
    Durand-Kerner_method)

    This function returns a tuple (roots, err) where roots is a list of
    ComplexFloats sorted by absolute value, and err is an estimate of
    the maximum error. The 'poly' argument should be a SymPy expression
    representing a univariate polynomial.

    Example:
        >>> Float.setdps(15)
        >>> x = Symbol('x')
        >>> r, e = polyroots((x-3)*(x-2))
        >>> r[0]
        ComplexFloat(real='1.9999999999999993', imag='-4.5917748078995606E-41')
        >>> r[1]
        ComplexFloat(real='3', imag='-4.6222318665293660E-33')
        >>> e
        Float('2.2204460492503131E-16')

    Polyroots attempts to achieve a maximum error less than the epsilon
    of the current working precision, but may fail to do so if the
    polynomial is badly conditioned. Usually the error can be reduced
    by increasing the 'maxsteps' parameter, although this will of
    course also reduce the speed proprtionally. It may also be
    necessary to increase the working precision.

    In particular, polyroots is easily fooled by multiple roots and
    roots that are very close together. In general, if n-multiple roots
    are present, polyroots will typically only locate their values to
    within 1/n of the working precision. For example, with the standard
    precision of ~15 decimal digits, the expected accuracy for a double
    root is 7-8 decimals:

        >>> Float.setdps(15)
        >>> r, e = polyroots((x-2)**2)
        >>> r[0]
        ComplexFloat(real='1.9999999832018203', imag='-2.5469994476319085E-8')
        >>> r[1]
        ComplexFloat(real='2.0000000110216827', imag='1.4784776969781909E-8')
        >>> e
        Float('5.4076851354766810E-8')

    An effective cure is to multiply the working precision n times (in
    the case of a double root, doubling it):

        >>> Float.setdps(30)
        >>> r, e = polyroots((x-2)**2, 40)
        >>> r[0].real
        Float('1.9999999999999999075008441778756')
        >>> e
        Float('1.9055980809840526328094566161453E-16')

    The result is good to 15 decimals as expected.
    """
    # Must be monic
    poly = Polynomial(poly)
    poly = Polynomial(poly / poly.coeffs[0][0])
    deg = poly.coeffs[0][1]
    f = polyfunc(poly)
    roots = [ComplexFloat(0.4+0.9j)**n for n in range(deg)]
    error = [Float(1) for n in range(deg)]
    for step in range(maxsteps):
        if max(error).ae(0):
            break
        for i in range(deg):
            if not error[i].ae(0):
                p = roots[i]
                x = f(p)
                for j in range(deg):
                    if i != j:
                        try:
                            x /= (p-roots[j])
                        except ZeroDivisionError:
                            continue
                roots[i] = p - x
                error[i] = abs(x)
    roots.sort(key=abs)
    err = max(error)
    err = max(err, Float((1, -Float._prec+1)))
    return roots, err


def sign(x):
    if x < 0:
        return -1
    return 1


def bisect(f, a, b, eps=None, maxsteps=None, verbose=False):
    """
    Numerical root-finding using the bisection method.

    Given a real-valued function f and an interval [a, b] (not
    necessarily real) such that f(a) and f(b) have opposite signs,
    narrow the interval where f crosses the x axis to a relative
    width at most equal to 'eps' through repeated bisections. If
    not specified, eps is set to give a full precision estimate
    of the root.

    A tuple for the narrowed interval is returned.

    If f is continuous, bisect is guaranteed to find a root. If f
    jumps discontinuously from positive negative values, bisect will
    locate the point of the discontinuity. bisect quits early and
    returns an interval of width zero if it encounters a point x
    where f(x) = 0 exactly.

    Optionally, perform no more than 'maxsteps' bisections (by
    default perform as many as needed); if the 'verbose' flag is set,
    print status at each step.

    Examples
    ========

    Find a bounding interval for pi/2 (which is a root of cos):

        >>> Float.setdps(15)
        >>> a, b = bisect(cos, 1, 2, 1e-2, verbose=True)
          bisect step  1: a=1.000000  b=2.000000  delta=1
          bisect step  2: a=1.500000  b=2.000000  delta=0.5
          bisect step  3: a=1.500000  b=1.750000  delta=0.25
          bisect step  4: a=1.500000  b=1.625000  delta=0.125
          bisect step  5: a=1.562500  b=1.625000  delta=0.0625
          bisect step  6: a=1.562500  b=1.593750  delta=0.03125
        >>> print a, b
        1.5625 1.578125

    Calculate the same value to full precision:

        >>> a, b = bisect(cos, 1, 2)
        >>> a, b
        (Float('1.5707963267948966'), Float('1.5707963267948968'))
        >>> print cos(a)
        6.12320117570134E-17
        >>> print cos(b)
        -1.60812593168018E-16

    Although generally reliable, the bisection method has a slow rate
    of convergence. The previous computation required about 50 steps,
    which can be compared to the secant method which only requires
    5-6 steps to obtain a full precision value from an initial value
    near 1.5.
    """

    a, b = a*Float(1), b*Float(1)   # Convert to Float or ComplexFloat
    fa0, fb0 = f(a), f(b)

    # Sanity check
    if not fa0 * fb0 <= 0:
        raise ValueError("bisect: f(a) and f(b) have the same sign at a=%r "
            " and b=%r" %  (a, b))

    fa0sign = sign(fa0)

    # Default eps / set eps to something sane if zero
    mineps = Float((1, -Float.getprec()+1))
    if not eps or eps < mineps:
        eps = mineps

    if maxsteps is None:
        # Use maxsteps as a safeguard should the loop condition fail
        abdelta = abs(a-b)
        maxsteps = bitcount(abdelta.man) - abdelta.exp + Float.getprec()

    i = 1
    while not a.ae(b, eps) and i <= maxsteps:
        if verbose:
            print "  bisect step %2i: a=%8f  b=%8f  delta=%g" % (i,a,b,b-a)
        mid = (a+b)*0.5
        fmid = f(mid)
        if fmid == 0:
            return (mid, mid)
        if sign(fmid) == fa0sign:
            a = mid
        else:
            b = mid
        i += 1

    return a, b


def perturb(x, eps):
    #rng.seed(x)
    if x == 0:
        return eps
    else:
        #return x * (1 + rng.random()*eps)
        return x * (1 + eps)


def secant(f, x0, x1=None, eps=None, maxsteps=50, bracket=None,
    on_error=1, verbose=False):
    """
    secant(f, x0) uses the secant method to find a numerical root of f,
    starting from the initial approximation/guess x0.

    This algorithm is very efficient for finding regular roots of
    smooth functions; nearly as efficient as Newton's method. Each
    iteration roughly multiplies the number of correct digits in the
    answer by 1.6 (in practice, a 15-digit value can usually be
    found in less than 10 steps).

    The secant method can be inefficient or may fail to converge in
    the presence of high-order roots or if the initial approximation
    is far from the true root. Generally, it works best if f is
    monotonic and steeply sloped everywhere in the surrounding of
    the root (and x0 is in that surrounding).

    Advanced settings
    =================

    The secant method actually requires two initial values x0 and x1.
    By default, an x1 value is generated by perturbing x0 slightly;
    in some cases this perturbation may be too small and it may then
    help to specify a custom value for x1.

    Several other advanced options are supported as guards against
    poor convergence:

      eps
          Only attempt to locate the root to within this epsilon.
          By default, eps is set to give a full-precision value.

      maxsteps
          Break when this many iterations have been performed.

      bracket
          If bracket = [a, b], break if the iteration ends up somewhere
          outside the interval [a, b]. This is useful to protect from
          the iterate x_n "shooting off to space" in the presence of
          a nearly horizontal function value. If the bracket is a
          single number, break if abs(x_0-x_n) exceeds this value.

    The 'on_error' parameter can be set to control what happens in the
    case of a convergence failure:

       on_error = 0   print warning if in verbose mode; return value
       on_error = 1   print warning; return value (default)
       on_error = 2   raise an exception

    Examples
    ========

    A simple example, calculating pi:

        >>> Float.setdps(15)
        >>> print secant(sin, 3)
        3.14159265358979

    If we try the same for cosine, starting from x=0, we get the
    "wrong" root because the cosine is nearly horizontal around x=0,
    sending the initial iterate far away (verbose=True shows what
    happens in more detail):

        >>> Float.setdps(15)
        >>> print secant(cos, 0)
        391.128285371929

    This can be helped by providing a second point that better matches
    the geometry of the function's graph or by providing a closer
    initial estimate:

        >>> print secant(cos, x0=0, x1=1)
        1.57079632679490
        >>> print secant(cos, 1.5)
        1.57079632679490

    As another example, a high-precision calculation of log(3):

        >>> Float.setdps(50)
        >>> print secant(lambda x: exp(x)-3, 1)
        1.0986122886681096913952452369225257046474905578227
        >>> Float.revert()
    """

    prec = Float.getprec()
    if eps is None:
        eps = Float((1, -prec+2))
    else:
        eps = Float(eps)

    x = x0 = x0 * Float(1)

    if x1 is None:
        # Perturb the initial value to generate a second point as
        # needed to start the secant method iteration
        xprev = perturb(x, 0.01)
    else:
        xprev = x1 * Float(1)

    if isinstance(bracket, (list, tuple)):
        brackets = Float(bracket[0]), Float(bracket[1])
    else:
        brackets = None

    exception = None
    bracketmsg = "Failed to converge (outside bracket)"
    stepsmsg = "Failed to converge (maxsteps exceeded)"
    derivmsg = "Failed to converge (inft./zero derivative or " \
        "loss of precision in derivative)"

    for i in xrange(1, maxsteps+1):
        if verbose:
            try:
                print "  secant: x_%i=%8f  delta=%g" % (i,x,x-xprev)
            except TypeError:
                print "  secant: x_%i=%s  delta=%s" % (i,x,x-xprev)

        # Calculate function values
        fx = f(x)
        fdiff = fx - f(xprev)
        xdiff = x - xprev

        if xdiff.ae(0, eps):
            break

        # Try to calculate the finite difference approximation for the
        # derivative and update accordingly. In floating-point
        # arithmetic, this may cause a division by zero when f(x) is
        # extremely close to 0, which we have to watch out for
        try:
            deriv = xdiff/fdiff
            x, xprev = x - fx*deriv, x
        except ZeroDivisionError:
            exception = derivmsg
            break

        # Check if within reasonable bounds
        if brackets:
            if not brackets[0] <= x <= brackets[1]:
                exception = bracketmsg
                break
        elif bracket is not None:
            print "heh", x, x0, x-x0, bracket
            if abs(x-x0) >= bracket:
                exception = bracketmsg
                break

    if i == maxsteps:
        exception = stepsmsg

    if exception:
        if on_error == 2: raise ConvergenceError(exception)
        if on_error == 1 or verbose: print exception

    return x
