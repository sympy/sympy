"""
Discrete Fourier Transform, Number Theoretic Transform,
Walsh Hadamard Transform, Zeta Transform, Mobius Transform
"""
from __future__ import print_function, division

from sympy.core import S, Symbol, sympify
from sympy.core.compatibility import as_int, range, iterable
from sympy.core.function import expand, expand_mul
from sympy.core.numbers import pi, I
from sympy.functions.elementary.exponential import exp
from sympy.functions.elementary.trigonometric import sin, cos
from sympy.ntheory import isprime, primitive_root
from sympy.utilities.iterables import ibin


#----------------------------------------------------------------------------#
#                                                                            #
#                         Discrete Fourier Transform                         #
#                                                                            #
#----------------------------------------------------------------------------#

def _fourier_transform(seq, symbolic, inverse=False):
    """Utility function for the Discrete Fourier Transform (DFT)"""

    if not iterable(seq):
        raise TypeError("Expected a sequence of numeric coefficients " +
                        "for Fourier Transform")

    a = [sympify(arg) for arg in seq]
    if any(x.has(Symbol) for x in a):
        raise ValueError("Expected non-symbolic coefficients")

    n = len(a)
    if n < 2:
        return a

    b = n.bit_length() - 1
    if n&(n - 1): # not a power of 2
        b += 1
        n = 2**b

    a += [S.Zero]*(n - len(a))
    for i in range(1, n):
        j = int(ibin(i, b, str=True)[::-1], 2)
        if i < j:
            a[i], a[j] = a[j], a[i]

    ang = -2*pi/n if inverse else 2*pi/n

    if not symbolic:
        ang = ang.evalf()

    w = [cos(ang*i) + I*sin(ang*i) for i in range(n // 2)]

    h = 2
    while h <= n:
        hf, ut = h // 2, n // h
        for i in range(0, n, h):
            for j in range(hf):
                u, v = a[i + j], expand_mul(a[i + j + hf]*w[ut * j])
                a[i + j], a[i + j + hf] = u + v, u - v
        h *= 2

    if inverse:
        for i in range(n):
            a[i] /= n

    return a


def fft(seq, symbolic=True):
    r"""
    Performs the Discrete Fourier Transform (DFT) in the complex domain.

    The sequence is automatically padded to the right with zeros, as the
    radix 2 FFT requires the number of sample points to be a power of 2.

    Parameters
    ==========

    seq : iterable
        The sequence on which DFT is to be applied.
    symbolic : bool
        Determines if DFT is to be performed using symbolic values or
        numerical values.

    Examples
    ========

    >>> from sympy import fft, ifft

    >>> fft([1, 2, 3, 4])
    [10, -2 - 2*I, -2, -2 + 2*I]
    >>> ifft(_)
    [1, 2, 3, 4]

    >>> ifft([1, 2, 3, 4])
    [5/2, -1/2 + I/2, -1/2, -1/2 - I/2]
    >>> fft(_)
    [1, 2, 3, 4]

    >>> fft([5])
    [5]
    >>> ifft([7])
    [7]

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
    .. [2] http://mathworld.wolfram.com/FastFourierTransform.html

    """

    return _fourier_transform(seq, symbolic=symbolic)


def ifft(seq, symbolic=True):
    return _fourier_transform(seq, symbolic=symbolic, inverse=True)

ifft.__doc__ = fft.__doc__


#----------------------------------------------------------------------------#
#                                                                            #
#                         Number Theoretic Transform                         #
#                                                                            #
#----------------------------------------------------------------------------#

def _number_theoretic_transform(seq, q, inverse=False):
    """Utility function for the Number Theoretic transform (NTT)"""

    if not iterable(seq):
        raise TypeError("Expected a sequence of integer coefficients " +
                        "for Number Theoretic Transform")

    q = as_int(q)
    if isprime(q) == False:
        raise ValueError("Expected prime modulo for " +
                        "Number Theoretic Transform")

    a = [as_int(x) for x in seq]

    n = len(a)
    if n < 1:
        return a

    b = n.bit_length() - 1
    if n&(n - 1):
        b += 1
        n = 2**b

    if (q - 1) % n:
        raise ValueError("Expected prime modulo of the form (m*2**k + 1)")

    a += [0]*(n - len(a))
    for i in range(1, n):
        j = int(ibin(i, b, str=True)[::-1], 2)
        if i < j:
            a[i], a[j] = a[j] % q, a[i] % q

    pr = primitive_root(q)

    rt = pow(pr, (q - 1) // n, q)
    if inverse:
        rt = pow(rt, q - 2, q)

    w = [1]*(n // 2)
    for i in range(1, n // 2):
        w[i] = w[i - 1] * rt % q

    h = 2
    while h <= n:
        hf, ut = h // 2, n // h
        for i in range(0, n, h):
            for j in range(hf):
                u, v = a[i + j], a[i + j + hf]*w[ut * j]
                a[i + j], a[i + j + hf] = (u + v) % q, (u - v) % q
        h *= 2

    if inverse:
        rv = pow(n, q - 2, q)
        for i in range(n):
            a[i] = a[i]*rv % q

    return a


def ntt(seq, q):
    r"""
    Performs the Number Theoretic Transform (NTT), which specializes the
    Discrete Fourier Transform (DFT) over quotient ring Z/pZ for prime p
    instead of complex numbers C.


    The sequence is automatically padded to the right with zeros, as the
    radix 2 NTT requires the number of sample points to be a power of 2.

    Examples
    ========

    >>> from sympy import ntt, intt
    >>> ntt([1, 2, 3, 4], 3*2**8 + 1)
    [10, 643, 767, 122]
    >>> intt(_, 3*2**8 + 1)
    [1, 2, 3, 4]
    >>> intt([1, 2, 3, 4], 3*2**8 + 1)
    [387, 415, 384, 353]
    >>> ntt(_, 3*2**8 + 1)
    [1, 2, 3, 4]

    References
    ==========

    .. [1] http://www.apfloat.org/ntt.html
    .. [2] http://mathworld.wolfram.com/NumberTheoreticTransform.html

    """

    return _number_theoretic_transform(seq, q)


def intt(seq, q):
    return _number_theoretic_transform(seq, q, inverse=True)

intt.__doc__ = ntt.__doc__
