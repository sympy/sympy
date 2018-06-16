"""
Recurrence evaluation
"""
from __future__ import print_function, division

from sympy.core import S, sympify
from sympy.core.compatibility import as_int, range, iterable

def reval_lhcc(coeffs, init, n):
    """
    Recurrence evaluation for type - linear, homogeneous, and
    with constant coefficients.

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

    Let, y be the recurrence of given type and c be the sequence of
    coeffs, and k (equal to len(c)) be the order of recurrence. Then,

    if n < k,
        y(n) = init(n)
    otherwise,
        y(n) = c(0)*y(n - 1) + c(1)*y(n - 2) + ... + c(k - 1)*y(n - k)

    """

    if not coeffs:
        return 0

    if not iterable(coeffs):
        raise TypeError("Expected a sequence of coefficients for the "
                        "recurrence")

    if not iterable(init):
        raise TypeError("Expected a sequence of values for the initialization"
                        "of the recurrence")

    c = [sympify(arg) for arg in coeffs] # coefficients of recurrence
    b = [sympify(arg) for arg in init] # intial/base values of recurrence
    n, k = as_int(n), len(c)

    if n < 0:
        raise ValueError("Point of evaluation of recurrence "
                        "must be a non-negative integer")

    if len(b) > k:
        raise TypeError("Count of initial values should not exceed "
                        "the count of coefficients of the recurrence")
    else:
        # remaining initial values are assumed to be 0
        b += [S.Zero]*(k - len(b))

    def _square_and_reduce_cht(u, offset):
        # Squares the coefficient vector u of length at
        # most k, to get coefficient vector of length at most
        # 2*k and reduces it to length at most k using
        # Cayley Hamilton Theorem (CHT)
        w = [S.Zero]*(2*len(u) - 1 + offset)
        for i, p in enumerate(u):
            for j, q in enumerate(u):
                w[offset + i + j] += p*q

        for j in range(len(w) - 1, k - 1, -1):
            for i in range(k):
                w[j - i - 1] += w[j]*c[i]

        return w[:k]

    def _final_coeffs(n):
        # Computes the final coefficient vector corresponding
        # to given point of evaluation of recurrence
        if n < k:
            a = [S.Zero]*k
            a[n] = S.One
            return a

        a = _final_coeffs(n // 2)

        return _square_and_reduce_cht(a, offset=n%2)

    return b[n] if n < k else sum(u*v for u, v in zip(_final_coeffs(n), b))
