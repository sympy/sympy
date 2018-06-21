"""
Recurrence
"""
from __future__ import print_function, division

from sympy.core import S, Symbol, sympify
from sympy.core.compatibility import as_int, range, iterable

def linrec(coeffs, init, n):
    """
    Linear recurrence evaluation of homogeneous type, having coefficients
    independent of the variable on which the recurrence is defined.

    Parameters
    ==========

    coeffs : iterable
        coefficients of recurrence
    init : iterable
        initial values of recurrence
    n : Integer
        point of evaluation of recurrence

    Explanation
    ===========

    Let, y(n) be the recurrence of given type, c be the sequence
    of coefficients, b be the sequence of intial/base values of the
    recurrence and k (equal to len(c)) be the order of recurrence.
    Then,

    if n < k,
        y(n) = b[n]
    otherwise,
        y(n) = c[0]*y(n - 1) + c[1]*y(n - 2) + ... + c[k - 1]*y(n - k)

    Examples
    ========

    >>> from sympy import linrec
    >>> from sympy.abc import x, y, z

    >>> linrec(coeffs=[1, 1], init=[0, 1], n=10)
    55
    >>> linrec(coeffs=[1, 1], init=[x, y], n=10)
    34*x + 55*y
    >>> linrec(coeffs=[x, y], init=[0, 1], n=5)
    x**2*y + x*(x**3 + 2*x*y) + y**2
    >>> linrec(coeffs=[1, 2, 3, 0, 0, 4], init=[x, y, z], n=16)
    13576*x + 5676*y + 2356*z

    """

    if not coeffs:
        return 0

    if not iterable(coeffs):
        raise TypeError("Expected a sequence of coefficients for"
                        " the recurrence")

    if not iterable(init):
        raise TypeError("Expected a sequence of values for the initialization"
                        " of the recurrence")

    n = as_int(n)
    if n < 0:
        raise ValueError("Point of evaluation of recurrence must be a "
                        "non-negative integer")

    c = [sympify(arg) for arg in coeffs]
    b = [sympify(arg) for arg in init]
    k = len(c)

    if len(b) > k:
        raise TypeError("Count of initial values should not exceed the "
                        "order of the recurrence")
    else:
        b += [S.Zero]*(k - len(b)) # remaining initial values default to zero

    def _square_and_reduce(u, offset):
        # squares (u[0] + u[1] * x + u[2] * x**2 + ... + u[k - 1] * x**(k - 1))
        # (multiplies by x if offset is 1) and reduces the above result of length
        # upto 2*k to k using the characterstic equation of the recurrence,
        # which is, x**k = c[0] * x**(k - 1) + c[1] * x**(k - 2) + ... + c[k - 1]

        w = [S.Zero]*(2*len(u) - 1 + offset)
        for i, p in enumerate(u):
            for j, q in enumerate(u):
                w[offset + i + j] += p*q

        for j in range(len(w) - 1, k - 1, -1):
            for i in range(k):
                w[j - i - 1] += w[j]*c[i]

        return w[:k]

    def _final_coeffs(n):
        # computes the final coefficient list - `cf` corresponding to the
        # point at which recurrence is to be evalauted - `n`, such that,
        # y(n) = cf[0]*y(k - 1) + cf[1]*y(k - 2) + ... + cf[k - 1]*y(0)

        if n < k:
            return [S.Zero]*n + [S.One] + [S.Zero]*(k - n - 1)
        else:
            return _square_and_reduce(_final_coeffs(n // 2), offset=n%2)

    return b[n] if n < k else sum(u*v for u, v in zip(_final_coeffs(n), b))
