"""
Convolution (using FFT, NTT, FWHT), Subset Convolution,
Covering Product, Intersecting Product
"""
from __future__ import print_function, division

from sympy.core import S
from sympy.core.compatibility import range, as_int
from sympy.core.function import expand_mul
from sympy.discrete.transforms import (
    fft, ifft, ntt, intt, fwht, ifwht)


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
        The following hints can be given as keyword arguments.
        dps : Integer
            Specifies the number of decimal digits for precision for
            performing FFT on the sequence.
        prime : Integer
            Prime modulus of the form (m*2**k + 1) to be used for
            performing NTT on the sequence.
        cycle : Integer
            Specifies the length for doing cyclic convolution.
        dyadic : bool
            Identifies the convolution type as dyadic (XOR)
            convolution, which is performed using FWHT.

    Examples
    ========

    >>> from sympy import convolution, symbols, S, I

    >>> convolution([1 + 2*I, 4 + 3*I], [S(5)/4, 6], dps=3)
    [1.25 + 2.5*I, 11.0 + 15.8*I, 24.0 + 18.0*I]

    >>> convolution([1, 2, 3], [4, 5, 6], cycle=3)
    [31, 31, 28]

    >>> convolution([111, 777], [888, 444], prime=19*2**10 + 1)
    [1283, 19351, 14219]

    >>> convolution([111, 777], [888, 444], prime=19*2**10 + 1, cycle=2)
    [15502, 19351]

    >>> u, v, x, y, z = symbols('u v x y z')
    >>> convolution([u, v], [x, y, z], dyadic=True)
    [u*x + v*y, u*y + v*x, u*z, v*z]

    """

    fft = hints.pop('fft', None)
    dps = hints.pop('dps', None)
    p = hints.pop('prime', None)
    c = as_int(hints.pop('cycle', 0))
    dyadic = hints.pop('dyadic', None)

    if c < 0:
        raise ValueError("The length for cyclic convolution must be non-negative")

    fft = True if fft else None
    dyadic = True if dyadic else None
    if sum(x is not None for x in (p, dps, dyadic)) > 1 or \
            sum(x is not None for x in (fft, dyadic)) > 1:
        raise TypeError("Ambiguity in determining the convolution type")

    if p is not None:
        ls = convolution_ntt(a, b, prime=p)
        return ls if not c else [sum(ls[i::c]) % p for i in range(c)]

    elif hints.pop('ntt', False):
        raise TypeError("Prime modulus must be specified for performing NTT")

    if dyadic:
        ls = convolution_fwht(a, b)
    else:
        ls = convolution_fft(a, b, dps=dps)

    return ls if not c else [sum(ls[i::c]) for i in range(c)]


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

    >>> from sympy import S, I
    >>> from sympy.discrete.convolution import convolution_fft

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

    if n > 0 and n&(n - 1): # not a power of 2
        n = 2**n.bit_length()

    # padding with zeros
    a += [S.Zero]*(n - len(a))
    b += [S.Zero]*(n - len(b))

    a, b = fft(a, dps), fft(b, dps)
    a = [expand_mul(x*y) for x, y in zip(a, b)]
    a = ifft(a, dps)[:m]

    return a


#----------------------------------------------------------------------------#
#                                                                            #
#                           Convolution for GF(p)                            #
#                                                                            #
#----------------------------------------------------------------------------#

def convolution_ntt(a, b, prime):
    """
    Performs linear convolution using Number Theoretic Transform.

    Parameters
    ==========

    a, b : iterables
        The sequences for which convolution is performed.
    prime : Integer
        Prime modulus of the form (m*2**k + 1) to be used for performing
        NTT on the sequence.

    Examples
    ========

    >>> from sympy.discrete.convolution import convolution_ntt

    >>> convolution_ntt([2, 3], [4, 5], prime=19*2**10 + 1)
    [8, 22, 15]
    >>> convolution_ntt([2, 5], [6, 7, 3], prime=19*2**10 + 1)
    [12, 44, 41, 15]
    >>> convolution_ntt([333, 555], [222, 666], prime=19*2**10 + 1)
    [15555, 14219, 19404]

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Convolution_theorem
    .. [2] https://en.wikipedia.org/wiki/Discrete_Fourier_transform_(general)

    """

    a, b, p = a[:], b[:], as_int(prime)
    n = m = len(a) + len(b) - 1 # convolution size

    if n > 0 and n&(n - 1): # not a power of 2
        n = 2**n.bit_length()

    # padding with zeros
    a += [0]*(n - len(a))
    b += [0]*(n - len(b))

    a, b = ntt(a, p), ntt(b, p)
    a = [x*y % p for x, y in zip(a, b)]
    a = intt(a, p)[:m]

    return a


#----------------------------------------------------------------------------#
#                                                                            #
#                         Convolution for 2**n-group                         #
#                                                                            #
#----------------------------------------------------------------------------#

def convolution_fwht(a, b):
    """
    Performs dyadic (XOR) convolution using Fast Walsh Hadamard Transform.

    The convolution is automatically padded to the right with zeros, as the
    radix 2 FWHT requires the number of sample points to be a power of 2.

    Parameters
    ==========

    a, b : iterables
        The sequences for which convolution is performed.

    Examples
    ========

    >>> from sympy import symbols, S, I
    >>> from sympy.discrete.convolution import convolution_fwht

    >>> u, v, x, y = symbols('u v x y')
    >>> convolution_fwht([u, v], [x, y])
    [u*x + v*y, u*y + v*x]

    >>> convolution_fwht([2, 3], [4, 5])
    [23, 22]
    >>> convolution_fwht([2, 5 + 4*I, 7], [6*I, 7, 3 + 4*I])
    [56 + 68*I, -10 + 30*I, 6 + 50*I, 48 + 32*I]

    >>> convolution_fwht([S(33)/7, S(55)/6, S(7)/4], [S(2)/3, 5])
    [2057/42, 1870/63, 7/6, 35/4]

    References
    ==========

    .. [1] https://researchgate.net/publication/26511536_Walsh_-_Hadamard_Transformation_of_a_Convolution
    .. [2] https://en.wikipedia.org/wiki/Hadamard_transform

    """

    if not a or not b:
        return []

    a, b = a[:], b[:]
    n = max(len(a), len(b))

    if n&(n - 1): # not a power of 2
        n = 2**n.bit_length()

    # padding with zeros
    a += [S.Zero]*(n - len(a))
    b += [S.Zero]*(n - len(b))

    a, b = fwht(a), fwht(b)
    a = [expand_mul(x*y) for x, y in zip(a, b)]
    a = ifwht(a)

    return a
