"""
Transcendental functions for real numbers:
* exp
* log
* sin/cos/tan
* sinh/cosh/tanh

"""

from util import *
from floatop import *
from squareroot import *
from constants import *
from convert import *


"""
The exponential function has a rapidly convergent Maclaurin series:

    exp(x) = 1 + x + x**2/2! + x**3/3! + x**4/4! + ...

The series can be summed very easily using fixed-point arithmetic.
The convergence can be improved further, using a trick due to
Richard P. Brent: instead of computing exp(x) directly, we choose a
small integer r (say, r=10) and compute exp(x/2**r)**(2**r).

The optimal value for r depends on the Python platform, the magnitude
of x and the target precision, and has to be estimated from
experimental timings. One test with x ~= 0.3 showed that
r = 2.2*prec**0.42 gave a good fit to the optimal values for r for
prec between 1 and 10000 bits, on one particular machine.

This optimization makes the summation about twice as fast at
low precision levels and much faster at high precision
(roughly five times faster at 1000 decimal digits).

If |x| is very large, we first rewrite it as t + n*log(2) with the
integer n chosen such that |t| <= log(2), and then calculate
exp(x) as exp(t)*(2**n), using the Maclaurin series for exp(t)
(the multiplication by 2**n just amounts to shifting the exponent).
"""

def exp_series(x, prec):
    r = int(2.2 * prec ** 0.42)
    # XXX: more careful calculation of guard bits
    guards = r + 3
    if prec > 60:
        guards += int(math.log(prec))
    prec2 = prec + guards
    x = rshift_quick(x, r - guards)
    s = (1 << prec2) + x
    a = x
    k = 2
    # Sum exp(x/2**r)
    while 1:
        a = ((a*x) >> prec2) // k
        if not a: break
        s += a
        k += 1
    # Calculate s**(2**r) by repeated squaring
    for j in range(r):
        s = (s*s) >> prec2
    return s >> guards

def fexp(x, prec, rounding):
    man, exp, bc = x
    # extra precision needs to be similar in magnitude to log_2(|x|)
    prec2 = prec + 6 + max(0, bc+exp)
    t = make_fixed(x, prec2)
    # abs(x) > 1?
    if exp+bc > 1:  #fcmp(fabs(x), fone) > 0:
        lg2 = log2_fixed(prec2)
        n, t = divmod(t, lg2)
    else:
        n = 0
    return normalize(exp_series(t, prec2), -prec2+n, prec, rounding)


"""
The basic strategy for computing log(x) is to set r = log(x) and use
Newton's method to solve the equation exp(r) = x. We set the initial
value r_0 to math.log(x) and then iterate r_{n+1} = r_n + exp(-r_n) - 1
until convergence. As with square roots, we increase the working
precision dynamically during the process so that only one full-precision
evaluation of exp is required.

log(x) is small for most inputs, so the r values can safely be
computed using fixed-point arithmetic. However, when x has a very
large or small exponent, we can improve performance through the
normalization log(t * 2**n) = log(t) + n*log(2), choosing n such
that 0.5 <= t <= 1 (for example).

There are some caveats: if x is extremely close to 1, the working
precision must be increased to maintain high relative precision in the
output (alternatively, the series approximation for log(1+x) could
be used in that case).
"""

# This function performs the Newton iteration using fixed-point
# arithmetic. x is assumed to have magnitude ~= 1
def _log_newton(x, prec):
    extra = 8
    # 50-bit approximation
    #r = int(_clog(Float((x, -prec), 64)) * 2.0**50)
    fx = math.log(to_float((x, -prec, bitcount(x))))
    r = int(fx * 2.0**50)
    prevp = 50
    for p in giant_steps(50, prec+extra):
        rb = lshift_quick(r, p-prevp)
        e = exp_series(-rb, p)
        r = rb + ((rshift_quick(x, prec-p)*e)>>p) - (1 << p)
        prevp = p
    return r >> extra

def flog(x, prec, rounding):
    if x == fzero: raise ValueError, "logarithm of 0"
    if x == fone:  return fzero
    man, exp, bc = x
    if man < 0: raise ValueError, "logarithm of a negative number"
    # Estimated precision needed for log(t) + n*log(2)
    prec2 = prec + int(math.log(1+abs(bc+exp), 2)) + 10
    # Watch out for the case when x is very close to 1
    if -1 < bc + exp < 2:
        near_one = fabs(fsub(x, fone, STANDARD_PREC, ROUND_FLOOR), STANDARD_PREC, ROUND_FLOOR)
        if near_one == 0:
            return fzero
        # estimate how close
        prec2 += -(near_one[1]) - bitcount(near_one[0])
    # Separate mantissa and exponent, calculate, join parts
    t = rshift_quick(man, bc-prec2)
    l = _log_newton(t, prec2)
    a = (exp + bc) * log2_fixed(prec2)
    return normalize(l+a, -prec2, prec, rounding)



"""
We compute sin(x) around 0 from its Taylor series, and cos(x) around 0
from sqrt(1-sin(x)**2). This way we can simultaneously compute sin and
cos, which are often needed together (e.g. for the tangent function or
the complex exponential), with little extra cost compared to computing
just one of them. The main reason for computing sin first (and not sin
from cos) is to obtain high relative accuracy for x extremely close to
0, where the operation sqrt(1-cos(x)**2) can cause huge cancellations.

For any value of x, we can reduce it to the interval A = [-pi/4, pi/4]
(where the Taylor series converges quickly) by translations, changing
signs, and switching the roles of cos and sin:

   A : sin(x) = sin(x)           cos(x) = cos(x)
   B : sin(x) = cos(x-pi/2)      cos(x) = -sin(x-pi/2)
   C : sin(x) = -sin(x-pi)       cos(x) = -cos(x-pi)
   D : sin(x) = -cos(x-3*pi/2)   cos(x) = sin(x-3*pi/2)

|     A      |      B     |      C     |     D     |
v            v            v            v           v

   1 |  ____   ..........                            ____
     |      _..          ..                        __
     |      . __           .                     __
     |    ..    _           ..                  _
     |   .       __           .               __
-----| -.----------_-----------.-------------_-----------
     | .            _           ..          _           .
     |               __           .       __           .
     |                 _           ..    _           ..
     |                  __           . __           .
     |                    __         _..          ..
  -1 |                      _________   ..........
      0                       pi                     2*pi


TODO: could use cos series too when extremely close to 0
"""

def _sin_series(x, prec):
    x2 = (x*x) >> prec
    s = a = x
    k = 3
    while a:
        a = ((a * x2) >> prec) // (-k*(k-1))
        s += a
        k += 2
    return s

def _trig_reduce(x, prec):
    pi_ = pi_fixed(prec)
    pi4 = pi_ >> 2
    pi2 = pi_ >> 1
    n, rem = divmod(x + pi4, pi2)
    rem -= pi4
    return n, rem

def cos_sin(x, prec, rounding):
    """Simultaneously compute (cos(x), sin(x)) for real x."""
    man, exp, bc = x
    bits_from_unit = abs(bc + exp)
    prec2 = prec + bits_from_unit + 15
    xf = make_fixed(x, prec2)
    n, rx = _trig_reduce(xf, prec2)
    case = n % 4
    one = 1 << prec2
    if case == 0:
        s = _sin_series(rx, prec2)
        c = sqrt_fixed(one - ((s*s)>>prec2), prec2)
    elif case == 1:
        c = -_sin_series(rx, prec2)
        s = sqrt_fixed(one - ((c*c)>>prec2), prec2)
    elif case == 2:
        s = -_sin_series(rx, prec2)
        c = -sqrt_fixed(one - ((s*s)>>prec2), prec2)
    elif case == 3:
        c = _sin_series(rx, prec2)
        s = -sqrt_fixed(one - ((c*c)>>prec2), prec2)
    c = normalize(c, -prec2, prec, rounding)
    s = normalize(s, -prec2, prec, rounding)
    return c, s

def fcos(x, prec, rounding):
    return cos_sin(x, prec, rounding)[0]

def fsin(x, prec, rounding):
    return cos_sin(x, prec, rounding)[1]

def ftan(x, prec, rounding):
    c, s = cos_sin(x, prec+6, ROUND_FLOOR)
    return fdiv(s, c, prec, rounding)


#----------------------------------------------------------------------
# Hyperbolic functions
#

def _sinh_series(x, prec):
    x2 = (x*x) >> prec
    s = a = x
    k = 3
    while a:
        a = ((a * x2) >> prec) // (k*(k-1))
        s += a
        k += 2
    return s

def cosh_sinh(x, prec, rounding):
    """Simultaneously compute (cosh(x), sinh(x)) for real x"""

    man, exp, bc = x
    high_bit = exp + bc
    prec2 = prec + 6

    if high_bit < -3:
        # Extremely close to 0, sinh(x) ~= x and cosh(x) ~= 1
        # TODO: support directed rounding
        if high_bit < -prec-2:
            return (fone, fpos(x, prec, rounding))

        # Avoid cancellation when computing sinh
        # TODO: might be faster to use sinh series directly
        prec2 += (-high_bit) + 4

    # In the general case, we use
    #    cosh(x) = (exp(x) + exp(-x))/2
    #    sinh(x) = (exp(x) - exp(-x))/2
    # and note that the exponential only needs to be computed once.
    ep = fexp(x, prec2, ROUND_FLOOR)
    em = fdiv(fone, ep, prec2, ROUND_FLOOR)
    ch = fshift_exact(fadd(ep, em, prec, rounding), -1)
    sh = fshift_exact(fsub(ep, em, prec, rounding), -1)
    return ch, sh

def fcosh(x, prec, rounding):
    """Compute cosh(x) for a real argument x"""
    return cosh_sinh(x, prec, rounding)[0]

def fsinh(x, prec, rounding):
    """Compute sinh(x) for a real argument x"""
    return cosh_sinh(x, prec, rounding)[1]

def ftanh(x, prec, rounding):
    """Compute tanh(x) for a real argument x"""
    ch, sh = cosh_sinh(x, prec+6, ROUND_FLOOR)
    return fdiv(sh, ch, prec, rounding)


#----------------------------------------------------------------------
# Inverse tangent
#

"""
Near x = 0, use atan(x) = x - x**3/3 + x**5/5 - ...
Near x = 1, use atan(x) = y/x * (1 + 2/3*y + 2*4/3/5*y**2 + ...)
where y = x**2/(1+x**2).

TODO: these series are not impressively fast. It is probably better
to calculate atan from tan, using Newton's method or even the
secant method.
"""

def _atan_series_1(x, prec, rounding):
    man, exp, bc = x
    # Increase absolute precision when extremely close to 0
    bc = bitcount(man)
    diff = -(bc + exp)
    prec2 = prec
    if diff > 10:
        if 3*diff - 4 > prec:  # x**3 term vanishes; atan(x) ~x
            return normalize(man, exp, prec, rounding)
        prec2 = prec + diff
    prec2 += 15  # XXX: better estimate for number of guard bits
    x = make_fixed(x, prec2)
    x2 = (x*x)>>prec2; one = 1<<prec2; s=a=x
    for n in xrange(1, 1000000):
        a = (a*x2) >> prec2
        s += a // ((-1)**n * (n+n+1))
        if -100 < a < 100:
            break
    return normalize(s, -prec2, prec, rounding)

def _atan_series_2(x, prec, rounding):
    prec2 = prec + 15
    x = make_fixed(x, prec2)
    one = 1<<prec2; x2 = (x*x)>>prec2; y=(x2<<prec2)//(one+x2)
    s = a = one
    for n in xrange(1, 1000000):
        a = ((a*y)>>prec2) * (2*n) // (2*n+1)
        if a < 100:
            break
        s += a
    return normalize(y*s//x, -prec2, prec, rounding)

_cutoff_1 = (5, -3, 3)   # ~0.6
_cutoff_2 = (3, -1, 2)   # 1.5

def fatan(x, prec, rounding):
    if x[0] < 0:
        t = fatan(fneg_exact(x), prec+4, ROUND_FLOOR)
        return normalize(-t[0], t[1], prec, rounding)
    if fcmp(x, _cutoff_1) < 0:
        return _atan_series_1(x, prec, rounding)
    if fcmp(x, _cutoff_2) < 0:
        return _atan_series_2(x, prec, rounding)
    # For large x, use atan(x) = pi/2 - atan(1/x)
    if x[1] > 10*prec:
        pi = fpi(prec, rounding)
        pihalf = pi[0], pi[1]-1, pi[2]
    else:
        pi = fpi(prec+4, ROUND_FLOOR)
        pihalf = pi[0], pi[1]-1, pi[2]
        t = fatan(fdiv(fone, x, prec+4, ROUND_FLOOR), prec+4, ROUND_FLOOR)
        return fsub(pihalf, t, prec, rounding)