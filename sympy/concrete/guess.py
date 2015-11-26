"""Various algorithms for helping identifying numbers and sequences."""
from __future__ import print_function, division

from sympy.utilities import public

from sympy.core import Function, var
from sympy import sympify, floor
from mpmath import pslq, sqrt, mp

@public
def find_simple_recurrence(v, A=Function('a'), N=var('n'), maxcoeff=1000):
    """
    Detects and returns a recurrence relation from a sequence of several integer
    (or rational) terms. The name of the function in the returned expression is
    'a' by default; the main variable is 'n' by default. The smallest index in
    the returned expression is always n (and never n-1, n-2, etc.).

    Examples
    ========

    >>> from sympy.concrete.guess import find_simple_recurrence
    >>> from mpmath import fib
    >>> find_simple_recurrence( [ fib(k) for k in range(12) ] )
    -a(n) - a(n + 1) + a(n + 2)

    >>> from sympy import Function, var
    >>> a = [1, 1, 1]
    >>> for k in range(15): a.append(5*a[-1]-3*a[-2]+8*a[-3])
    >>> find_simple_recurrence(a, A=Function('f'), N=var('i'))
    -8*f(i) + 3*f(i + 1) - 5*f(i + 2) + f(i + 3)

    """
    l = len(v)>>1

    previous = mp.prec # save current precision
    mp.prec = 128
    b = [ sum( sqrt((l>>1)**2 + k)*v[-1-k-i] for k in range(l) )
          for i in range(l) ]
    p = pslq(   b,
                maxcoeff = maxcoeff,
                maxsteps = 128 + 4*l
            )
    mp.prec = previous # restore current precision

    first, last = 0, l-1
    while p[first]==0: first +=1
    while p[last]==0: last -=1

    rel = 0
    for k in range(first, last+1):
        rel += A(N+last-k)*p[k]

    return rel


@public
def rationalize(x, maxcoeff=10000):
    """
    Helps identifying a rational number from a float (or mpmath.mpf) value by
    using a continued fraction. The algorithm stops as soon as a large partial
    quotient is detected (greater than 10000 by default).

    Examples
    ========

    >>> from sympy.concrete.guess import rationalize
    >>> from mpmath import cos, pi
    >>> rationalize( cos(pi/3) )
    1/2

    >>> from mpmath import mpf
    >>> rationalize( mpf("0.333333333333333") )
    1/3

    While the function is rather intended to help 'identifying' rational
    values, it may be used in some cases for approximating real numbers.
    (Though other functions may be more relevant in that case.)

    >>> rationalize( pi, maxcoeff = 250 )
    355/113

    See also
    ========
    Several other methods can approximate a real number as a rational, like:

      * fractions.Fraction.from_decimal
      * fractions.Fraction.from_float
      * mpmath.pslq by using the following syntax: mpmath.pslq([x, 1])
      * sympy.simplify.nsimplify (which is a more general function)

    The main difference between the current function and all these variants is
    that control focuses on magnitude of partial quotients here rather than on
    global precision of the approximation. If the real is "known to be" a
    rational number, the current function should be able to detect it correctly
    with the default settings even when denominator is great (unless its
    expansion contains unusually big partial quotients) which may occur
    when studying sequences of increasing numbers. If the user cares more
    on getting simple fractions, other methods may be more convenient.

    """
    p0, p1 = 0, 1
    q0, q1 = 1, 0
    a = floor(x)
    while a < maxcoeff:
        p = a*p1 + p0
        q = a*q1 + q0
        p0, p1 = p1, p
        q0, q1 = q1, q
        if x==a: break
        x = 1/(x-a)
        a = floor(x)
    return sympify(p) / q


#
# TODO
# ====
# Roadmap: add a function for guessing a generating function by using
# Pade approximants
