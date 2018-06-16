"""
Recurrence evaluation
"""
from __future__ import print_function, division

from sympy.core import S, Symbol, sympify
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

    Let, y(n) be the recurrence of given type, c be the sequence
    of coefficients, b be the sequence of intial/base values of the
    reccurence and k (equal to len(c)) be the order of recurrence.
    Then,

    if n < k,
        y(n) = b[n]
    otherwise,
        y(n) = c[0]*y(n - 1) + c[1]*y(n - 2) + ... + c[k - 1]*y(n - k)

    """

    if not coeffs:
        return 0

    if not iterable(coeffs):
        raise TypeError("Expected a sequence of constant coefficients for"
                        " the recurrence")

    if not iterable(init):
        raise TypeError("Expected a sequence of values for the initialization"
                        " of the recurrence")

    c = [sympify(arg) for arg in coeffs]
    b = [sympify(arg) for arg in init]
    n, k = as_int(n), len(c)

    if n < 0:
        raise ValueError("Point of evaluation of recurrence must be a "
                        "non-negative integer")

    if any(x.has(Symbol) for x in c):
        raise ValueError("Expected a sequence of constant coefficients for"
                        " the recurrence")

    if len(b) > k:
        raise TypeError("Count of initial values should not exceed the "
                        "order of the recurrence")
    else:
        b += [S.Zero]*(k - len(b)) # remaining initial values default to zero

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
            return [S.Zero]*n + [S.One] + [S.Zero]*(k - n - 1)
        else:
            return _square_and_reduce_cht(_final_coeffs(n // 2), offset=n%2)

    return b[n] if n < k else sum(u*v for u, v in zip(_final_coeffs(n), b))
