"""
Discrete Fourier Transform, Walsh Hadamard Transform,
Number Theoretic Transform, Zeta Transform, Mobius Transform
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
    w = None
    if symbolic:
        w = [cos(ang*i) + I*sin(ang*i) for i in range(n // 2)]
    else:
        w = [exp(I*ang*i).evalf() for i in range(n // 2)]

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
