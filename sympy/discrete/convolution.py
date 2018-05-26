"""
Linear Convolution, Arithmetic Convolution (bitwise XOR, AND, OR),
Subset Convolution, Superset Convolution
"""
from __future__ import print_function, division

from sympy.core import S, sympify
from sympy.core.compatibility import range, iterable
from sympy.core.function import expand_mul
from sympy.discrete.transforms import fft, ifft, ntt, intt


def convolution(a, b, **hints):
    """
    Performs convolution by determining the type of desired
    convolution using hints.

    If no hints are given, linear convolution is performed using
    FFT.

    Parameters
    ==========

    a, b : iterables
        The sequences for which convolution is performed.
    hints : dict
        Specifies the type of convolution to be performed.

        dps : Integer
            Specifies the number of decimal digits for precision for
            performing FFT on the sequence.

        p : Integer
            Prime modulus of the form (m*2**k + 1) to be used for
            performing NTT on the sequence.

    Examples
    ========

    >>> from sympy import S, I, convolution

    >>> convolution([1 + 2*I, 4 + 3*I], [S(5)/4, 6], dps=3)
    [1.25 + 2.5*I, 11.0 + 15.8*I, 24.0 + 18.0*I]

    >>> convolution([111, 777], [888, 444], p=19*2**10 + 1)
    [1283, 19351, 14219]
    """

    p = hints.pop('p', None)
    dps = hints.pop('dps', None)

    if sum(x is not None for x in (p, dps)) > 1:
        raise TypeError("Ambiguity in determining the convolution type")

    if p is not None:
        return convolution_ntt(a, b, p)
    elif hints.pop('ntt', False):
        raise TypeError("Prime modulus must be specified for performing NTT")

    return convolution_fft(a, b, dps)


#----------------------------------------------------------------------------#
#                                                                            #
#                       Convolution for Complex domain                       #
#                                                                            #
#----------------------------------------------------------------------------#

def convolution_fft(a, b, dps=None):
    """
    Performs linear convolution using Fast Fourier Transform.

    Parameters
    ==========

    a, b : iterables
        The sequences for which convolution is performed.
    dps : Integer
        Specifies the number of decimal digits for precision.

    Examples
    ========

    >>> from sympy import S, I, convolution_fft

    >>> convolution_fft([2, 3], [4, 5])
    [8, 22, 15]
    >>> convolution_fft([2, 5], [6, 7, 3])
    [12, 44, 41, 15]

    >>> convolution_fft([1 + 2*I, 4 + 3*I], [S(5)/4, 6])
    [5/4 + 5*I/2, 11 + 63*I/4, 24 + 18*I]

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Convolution_theorem
    .. [1] https://en.wikipedia.org/wiki/Discrete_Fourier_transform_(general)

    """
    a, b = a[:], b[:]
    n = m = len(a) + len(b) - 1 # convolution size

    nb = n.bit_length()
    if n&(n - 1): # not a power of 2
        n = 2**nb

    # padding with zeros
    a += [S.Zero]*(n - len(a))
    b += [S.Zero]*(n - len(b))

    a, b = fft(a, dps), fft(b, dps)
    for i in range(n):
        a[i] = expand_mul(a[i]*b[i])

    a = ifft(a, dps)[:m]

    return a


#----------------------------------------------------------------------------#
#                                                                            #
#                           Convolution for GF(p)                            #
#                                                                            #
#----------------------------------------------------------------------------#

def convolution_ntt(a, b, p):
    """
    Performs linear convolution using Number Theoretic Transform.

    Parameters
    ==========

    a, b : iterables
        The sequences for which convolution is performed.
    p : Integer
        Prime modulus of the form (m*2**k + 1) to be used for performing
        NTT on the sequence.

    Examples
    ========

    >>> from sympy import convolution_ntt

    >>> convolution_ntt([2, 3], [4, 5], 19*2**10 + 1)
    [8, 22, 15]
    >>> convolution_ntt([2, 5], [6, 7, 3], 19*2**10 + 1)
    [12, 44, 41, 15]

    >>> convolution_ntt([333, 555], [222, 666], 19*2**10 + 1)
    [15555, 14219, 19404]

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Convolution_theorem
    .. [2] https://en.wikipedia.org/wiki/Discrete_Fourier_transform_(general)

    """

    a, b = a[:], b[:]
    n = m = len(a) + len(b) - 1 # convolution size

    nb = n.bit_length()
    if n&(n - 1): # not a power of 2
        n = 2**nb

    # padding with zeros
    a += [0]*(n - len(a))
    b += [0]*(n - len(b))

    a, b = ntt(a, p), ntt(b, p)
    for i in range(n):
        a[i] = a[i]*b[i] % p

    a = intt(a, p)[:m]

    return a
