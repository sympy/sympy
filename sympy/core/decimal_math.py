"""
Decimal math functions.

Heavily copied from http://code.google.com/p/dmath/, generalized pow to non-integer exponents.

"""

__all__ = ['acos', 'asin', 'atan', 'atan2', 'ceiling', 'cos', 'cosh', 'degrees',
           'e', 'exp', 'floor', 'golden_ratio', 'hypot', 'log', 'log10', 'pi',
           'pow', 'radians', 'sign', 'sin', 'sinh', 'sqrt', 'tan', 'tanh',
           'cot','coth']

import math, decimal

from decimal import getcontext
from decimal import Decimal as D

def pi():
    """Compute Pi to the current precision."""
    getcontext().prec += 2
    lasts, t, s, n, na, d, da = 0, D(3), 3, 1, 0, 0, 24
    while s != lasts:
        lasts = s
        n, na = n + na, na + 8
        d, da = d + da, da + 32
        t = (t * n) / d
        s += t
    getcontext().prec -= 2
    return +s

def e():
    """Compute the base of the natural logarithm to the current precision."""
    return exp(D(1))

def golden_ratio():
    """Calculate the golden ratio to the current precision."""
    return  +((1 + D(5).sqrt()) / 2)


def exp(x):
    """Return e raised to the power of x.
    """
    getcontext().prec += 2
    i, lasts, s, fact, num = 0, 0, 1, 1, 1
    while s != lasts:
        lasts = s    
        i += 1
        fact *= i
        num *= x     
        s += num / fact   
    getcontext().prec -= 2
    return +s

def cos(x):
    """Return the cosine of x as measured in radians.
    """
    getcontext().prec += 2
    i, lasts, s, fact, num, sign = 0, 0, 1, 1, 1, 1
    while s != lasts:
        lasts = s    
        i += 2
        fact *= i * (i - 1)
        num *= x * x
        sign *= -1
        s += num / fact * sign 
    getcontext().prec -= 2        
    return +s

def sin(x):
    """Return the sine of x as measured in radians.
    """
    getcontext().prec += 2
    i, lasts, s, fact, num, sign = 1, 0, x, 1, x, 1
    while s != lasts:
        lasts = s    
        i += 2
        fact *= i * (i - 1)
        num *= x * x
        sign *= -1
        s += num / fact * sign
    getcontext().prec -= 2
    return +s

def cosh(x):
    """Return the hyperbolic cosine of Decimal x."""
    if x == 0:
        return D(1)
    
    getcontext().prec += 2
    i, lasts, s, fact, num = 0, 0, 1, 1, 1
    while s != lasts:
        lasts = s
        i += 2
        num *= x * x
        fact *= i * (i - 1)
        s += num / fact
    getcontext().prec -= 2
    return +s

def sinh(x):
    """Return the hyperbolic sine of Decimal x."""
    if x == 0:
        return D(0)
    
    getcontext().prec += 2
    i, lasts, s, fact, num = 1, 0, x, 1, x
    while s != lasts:
        lasts = s
        i += 2
        num *= x * x
        fact *= i * (i - 1)
        s += num / fact
    getcontext().prec -= 2
    return +s

def asin(x):
    """Return the arc sine (measured in radians) of Decimal x."""
    if abs(x) > 1:
        raise ValueError("Domain error: asin accepts -1 <= x <= 1")
    
    if x == -1:
        return pi() / -2
    elif x == 0:
        return D(0)
    elif x == 1:
        return pi() / 2
    
    return atan2(x, D.sqrt(1 - x ** 2))

def acos(x):
    """Return the arc cosine (measured in radians) of Decimal x."""
    if abs(x) > 1:
        raise ValueError("Domain error: acos accepts -1 <= x <= 1")

    if x == -1:
        return pi()
    elif x == 0:
        return pi() / 2
    elif x == 1:
        return D(0)
    
    return pi() / 2 - atan2(x, D.sqrt(1 - x ** 2))

def tan(x):
    """Return the tangent of Decimal x (measured in radians)."""
    return +(sin(x) / cos(x))

def tanh(x):
    """Return the hyperbolic tangent of Decimal x."""
    return +(sinh(x) / cosh(x))

def cot(x):
    """Return the cotangent of Decimal x (measured in radians)."""
    return +(cos(x) / sin(x))

def coth(x):
    """Return the hyperbolic cotangent of Decimal x."""
    return +(cosh(x) / sinh(x))

def atan(x):
    """Return the arc tangent (measured in radians) of Decimal x."""
    if x == D('-Inf'):
        return pi() / -2
    elif x == 0:
        return D(0)
    elif x == D('Inf'):
        return pi() / 2
    
    if x < -1:
        c = pi() / -2
        x = 1 / x
    elif x > 1:
        c = pi() / 2
        x = 1 / x
    else:
        c = 0
    
    getcontext().prec += 2
    x_squared = x ** 2
    y = x_squared / (1 + x_squared)
    y_over_x = y / x
    i, lasts, s, coeff, num = D(0), 0, y_over_x, 1, y_over_x
    while s != lasts:
        lasts = s 
        i += 2
        coeff *= i / (i + 1)
        num *= y
        s += coeff * num
    if c:
        s = c - s
    getcontext().prec -= 2
    return +s

def sign(x):
    """Return -1 for negative numbers and 1 for positive numbers."""
    return 2 * D(x >= 0) - 1

def atan2(y, x):
    """Return the arc tangent (measured in radians) of y/x.
    Unlike atan(y/x), the signs of both x and y are considered.
    """
    abs_y = abs(y)
    abs_x = abs(x)
    y_is_real = abs_y != D('Inf')
    
    if x:
        if y_is_real:
            a = y and atan(y / x) or D(0)
            if x < 0:
                a += sign(y) * pi()
            return a
        elif abs_y == abs_x:
            x = sign(x)
            y = sign(y)
            return pi() * (D(2) * abs(x) - x) / (D(4) * y)
    if y:
        return atan(sign(y) * D('Inf'))
    elif x < 0:
        return sign(y) * pi()
    else:
        return D(0)

def log(x, base=None):
    """log(x[, base]) -> the logarithm of Decimal x to the given Decimal base.
    If the base not specified, returns the natural logarithm (base e) of x.
    """
    if x < 0:
        return D('NaN')
    elif base == 1:
        raise ValueError("Base was 1!")
    elif x == base:
        return D(1)
    elif x == 0:
        return D('-Inf')

    getcontext().prec += 2

    if x < 1 or x > 10:
        # Consider x = a / 10^b, with |b| > 1 and 1 < a < 10. As long as our
        # log calculation is fast for 1 <= x <= 10, then we can consider x
        # in this form and then log(x) = log(a/10^b) = log(a) - b*log(10)
        r1 = list(x.as_tuple())
        r1[2] = 1 - len(r1[1])
        a = D(tuple(r1))

        r2 = list(x.as_tuple())
        b = 1 - r2[2] - len(r2[1])

        ret = log(a, base) - b*log(D(10), base)
        getcontext().prec -= 2
        return +ret

    if base is None:
        log_base = 1
        approx = math.log(x)
    else:
        log_base = log(base)
        approx = math.log(x, base)

    lasts, s = 0, D(repr(approx))
    while lasts != s:
        lasts = s
        s = s - 1 + x / exp(s)
    s /= log_base

    getcontext().prec -= 2
    return +s

def log10(x):
    """log10(x) -> the base 10 logarithm of Decimal x."""
    return log(x, D(10))

def pow(x, y):
    """pow(x,y) -> x ** y."""
    if isinstance(y,(int,long)):
        return D.__pow__(x, y)
    if isinstance(y,D) and y==D('0.5'):
        return sqrt(x)
    return exp(y * log(x)) # this is slow

sqrt = D.sqrt

def degrees(x):
    """degrees(x) -> converts Decimal angle x from radians to degrees"""
    return +(x * 180 / pi())

def radians(x):
    """radians(x) -> converts Decimal angle x from degrees to radians"""
    return +(x * pi() / 180)

def ceiling(x):
    """Return the smallest integral value >= x."""
    return x.to_integral(rounding=decimal.ROUND_CEILING)

def floor(x):
    """Return the largest integral value <= x."""
    return x.to_integral(rounding=decimal.ROUND_FLOOR)

def hypot(x, y):
    """Return the Euclidean distance, sqrt(x*x + y*y)."""
    return sqrt(x * x + y * y)

# EOF
