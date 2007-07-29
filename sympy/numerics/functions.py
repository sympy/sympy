"""
This module implements standard mathematical functions
(exp, log, sin, ...) for real arguments. Complex versions will
be provied in a complementary file cfunctions.py.
"""

from float_ import Float
from constants import pi_float, log2_fixed, log2_float, log10_float
from utils_ import bitcount, make_fixed
from math import log as _clog

def _quadratic_steps(start, target):
    L = [target]
    while L[-1] > start*2:
        L = L + [L[-1]//2 + 1]
    return L[::-1]

def _rshift(x, n):
    if n >= 0: return x >> n
    else:      return x << (-n)

def _lshift(x, n):
    if n >= 0: return x << n
    else:      return x >> (-n)

def _with_extraprec(n, f):
    Float._prec += n
    v = f()
    Float._prec -= n
    return +v


#----------------------------------------------------------------------
# Square root
#

"""
Square roots are most efficiently computed with Newton's method.
Two functions are implemented: _sqrt_fixed and _sqrt_fixed2.

  _sqrt_fixed uses the iteration r_{n+1} = (r_n + y/r_n)/2,
  which is just Newton's method applied to the equation r**2 = y.

  _sqrt_fixed2 uses the iteration r_{n+1} = r_n*(3 - y*r_n**2)
  to calculate 1/sqrt(y), and then multiplies by y to obtain
  sqrt(y).

The first iteration is slightly faster at low precision levels, since it
essentially just requires one division at each step, compared
to the three multiplications in the second formula. However, the second
iteration is much better at extremely high precision levels. This is
likely due to the fact that Python uses the Karatsuba algorithm
for integer multiplication, which is asymptotically faster than
its division algorithm.

For optimal speed, we exploit the "self-correcting" nature of
Newton's method to perform subcomputations at as low a precision level
as possible. Starting from a 50-bit floating-point estimate, the
first step can be computed using 100-bit precision, the second
at 200-bit precision, and so on; full precision is only needed for
the final step.

Both functions use fixed-point arithmetic and assume that the input y
is a big integer, i.e. given the integer y and precision prec,
they return floor(sqrt(x) * 2**prec) where y = floor(x * 2**prec).

The functions currently assume that x ~= 1. (TODO: make the code
work for x of arbitrary magnitude.) The main sqrt() function
fiddles with the exponent of the input to reduce it to unit
magnitude before passing it to _sqrt_fixed or _sqrt_fixed2.

"""

def _sqrt_fixed(y, prec):
    # get 50-bit initial guess from regular float math
    if prec < 200:
        r = int(y**0.5 * 2.0**(50-prec*0.5))
    else:
        r = int((y >> (prec-100))**0.5)
    prevp = 50
    for p in _quadratic_steps(50, prec+8):
        # Newton iteration: r_{n+1} = (r_{n} + y/r_{n})/2
        # print "sqrt", p
        r = _lshift(r, p-prevp-1) + (_rshift(y, prec-p-prevp+1)//r)
        prevp = p
    return r >> 8

def _sqrt_fixed2(y, prec):
    r = float(Float((y, -prec), 64))**-0.5
    r = int(r * 2**50)
    prevp = 50
    for p in _quadratic_steps(50, prec+8):
        # print "sqrt", p
        r2 = _rshift(r*r, 2*prevp - p)
        A = _lshift(r, p-prevp)
        T = _rshift(y, prec-p)
        S = (T*r2) >> p
        B = (3 << p) - S
        r = (A*B)>>(p+1)
        prevp = p
    r = (r * y) >> prec
    return r >> 8

def sqrt(x):
    """
    sqrt(x) returns the square root of x as a Float, rounded to
    the current working precision. This function only handles real-
    valued square roots (it raises ValueError if x is negative).
    """
    if not isinstance(x, Float):
        x = Float(x)
    if x == 0:
        return Float(0)
    if x < 0:
        raise ValueError
    prec = Float._prec + 4
    # Convert to a fixed-point number with prec bits. Adjust
    # exponents to be even so that they can be divided in half
    if prec & 1: prec += 1
    man = x.man
    exp = x.exp
    if exp & 1:
        exp -= 1
        man <<= 1
    shift = bitcount(man) - prec
    shift -= shift & 1
    man = _rshift(man, shift)
    if prec < 65000:
        man = _sqrt_fixed(man, prec)
    else:
        man = _sqrt_fixed2(man, prec)
    return Float((man, (exp+shift-prec)//2))

def hypot(x, y):
    """hypot(x, y) computes the distance function sqrt(x**2+y**2)"""
    if not isinstance(x, Float): x = Float(x)
    if not isinstance(y, Float): y = Float(y)
    if not y: return abs(x)
    if not x: return abs(y)
    return _with_extraprec(4, lambda: sqrt(x*x + y*y))


#----------------------------------------------------------------------
# Exponential function
#

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

def _exp_series(x, prec):
    r = int(2.2 * prec ** 0.42)
    # XXX: more careful calculation of guard bits
    guards = r + 3
    if prec > 60:
        guards += int(_clog(prec))
    prec2 = prec + guards
    x = _rshift(x, r - guards)
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

def exp(x):
    """
    exp(x) -- compute the exponential function of the real number x
    """
    if not isinstance(x, Float):
        x = Float(x)
    # extra precision needs to be similar in magnitude to log_2(|x|)
    prec = Float._prec + 4 + max(0, bitcount(x.man) + x.exp)
    t = make_fixed(x, prec)
    if abs(x) > 1:
        lg2 = log2_fixed(prec)
        n, t = divmod(t, lg2)
    else:
        n = 0
    y = _exp_series(t, prec)
    return Float((y, -prec+n))


#----------------------------------------------------------------------
# Natural logarithm
#

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
    # 50-bit approximation
    r = int(_clog(Float((x, -prec), 64)) * 2.0**50)
    prevp = 50
    for p in _quadratic_steps(50, prec+8):
        rb = _lshift(r, p-prevp)
        e = _exp_series(-rb, p)
        r = rb + ((_rshift(x, prec-p)*e)>>p) - (1 << p)
        prevp = p
    return r >> 8

def log(x, b=None):
    """
    log(x)    -> the natural (base e) logarithm of x
    log(x, b) -> the base b logarithm of x

    Both x and b must be positive real numbers (with b != 1).
    """
    # Basic input management
    if b is not None:
        Float._prec += 3
        if   b == 2:  blog = log2_float()
        elif b == 10: blog = log10_float()
        else: blog = log(b)
        l = log(x) / blog
        Float._prec -= 3
        return l
    if not isinstance(x, Float): x = Float(x)
    if x == 1: return Float((0, 0))
    if x <= 0: raise ValueError, "complex logarithm not implemented"

    bc = bitcount(x.man)
    # Estimated precision needed for log(t) + n*log(2)
    prec = Float._prec + int(_clog(1+abs(bc+x.exp), 2)) + 10

    # Watch out for the case when x is very close to 1
    if -1 < bc+x.exp < 2:
        near_one = abs(x-1)
        if near_one == 0:
            return Float((0, 0))
        prec += -(near_one.exp) - bitcount(near_one.man)

    # Separate mantissa and exponent, calculate fixed-point
    # approximation and put it all together
    t = _rshift(x.man, bc-prec)
    l = _log_newton(t, prec)
    a = (x.exp+bc) * log2_fixed(prec)
    return Float((l+a, -prec))


#----------------------------------------------------------------------
# Trigonometric functions
#

"""
We compute sin(x) around 0 from its Taylor series, and cos(x) around 0
from sqrt(1-sin(x)**2). This way we can simultaneously compute sin and
cos, which are often needed together (e.g. for the tangent function or
the complex exponential), with little extra cost compared to computing
just one of them. The main reason for computing sin first (and not cos
from sin) is to obtain high relative accuracy for x extremely close to
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

TODO: check if the reduced x is very close to 0 and in that case repeat
the mod operation at a higher level of precision.
"""

from constants import pi_fixed

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

def cos_sin(x):
    """
    cos_sin(x) calculates both the cosine and the sine of x rounded
    to the nearest Float value, and returns the tuple (cos(x), sin(x)).
    """
    if not isinstance(x, Float):
        x = Float(x)
    bits_from_unit = abs(bitcount(x.man) + x.exp)
    prec = Float._prec + bits_from_unit + 15
    xf = make_fixed(x, prec)
    n, rx = _trig_reduce(xf, prec)
    case = n % 4
    one = 1<<prec
    if case == 0:
        s = _sin_series(rx, prec)
        c = _sqrt_fixed(one - ((s*s)>>prec), prec)
    elif case == 1:
        c = -_sin_series(rx, prec)
        s = _sqrt_fixed(one - ((c*c)>>prec), prec)
    elif case == 2:
        s = -_sin_series(rx, prec)
        c = -_sqrt_fixed(one - ((s*s)>>prec), prec)
    elif case == 3:
        c = _sin_series(rx, prec)
        s = -_sqrt_fixed(one - ((c*c)>>prec), prec)
    return Float((c, -prec)), Float((s, -prec))

def cos(x):
    return cos_sin(x)[0]

def sin(x):
    return cos_sin(x)[1]

def tan(x):
    Float._prec += 2
    c, s = cos_sin(x)
    t = s / c
    Float._prec -= 2
    return t


#----------------------------------------------------------------------
# Inverse tangent
#

"""
Near x = 0, use atan(x) = x - x**3/3 + x**5/5 - ...
Near x = 1, use atan(x) = y/x * (1 + 2/3*y + 2*4/3/5*y**2 + ...)
where y = x**2/(1+x**2).

TODO: these series are not impressively fast. It is probably better
to calculate atan from tan, using Newton's method.

"""

def _atan_series_1(x):
    prec = Float._prec
    # Increase absolute precision when extremely close to 0
    bc = bitcount(x.man)
    diff = -(bc + x.exp)
    if diff > 10:
        if 3*diff - 4 > prec: # x**3 term vanishes; atan(x) ~x
            return +x
        prec = prec + diff
    prec += 15 # XXX: better estimate for number of guard bits
    x = make_fixed(x, prec)
    x2 = (x*x)>>prec; one = 1<<prec; s=a=x
    for n in xrange(1, 1000000):
        a = (a*x2) >> prec
        s += a // ((-1)**n * (n+n+1))
        if -100 < a < 100:
            break
    return Float((s, -prec))

def _atan_series_2(x):
    prec = Float._prec
    prec = prec + 15
    x = make_fixed(x, prec)
    one = 1<<prec; x2 = (x*x)>>prec; y=(x2<<prec)//(one+x2)
    s = a = one
    for n in xrange(1, 1000000):
        a = ((a*y)>>prec) * (2*n) // (2*n+1)
        if a < 100: break
        s += a
    return Float(((y*s)//x, -prec))

def atan(x):
    """Compute atan(x) for a real number x"""
    if not isinstance(x, Float): x = Float(x)
    if x < -0.6: return _with_extraprec(2, -atan(-x))
    if x < 0.6: return _atan_series_1(x)
    if x < 1.5: return _atan_series_2(x)
    # For large x, use atan(x) = pi/2 - atan(1/x)
    Float._prec += 4
    prec = Float._prec
    if x.exp > 10*prec:
        t = pi_float()/2   # XXX
    else:
        t = pi_float()/2 - atan(Float(1)/x)
    Float._prec -= 4
    return +t

def atan2(y,x):
    """atan2(y,x) has the same magnitude as atan(y/x) but
    accounts for the signs of y and x"""
    if y < 0: return -atan2(-y, x)
    if not x and not y: return Float(0)
    if y > 0 and x == 0:
        Float._prec += 2; t = pi_float()/2; Float._prec -= 2
        return t
    Float._prec += 2
    if x > 0: a = atan(y/x)
    else:     a = pi_float() - atan(-y/x)
    Float._prec -= 2
    return +a
