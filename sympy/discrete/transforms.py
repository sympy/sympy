"""
Discrete Fourier Transform, Walsh Hadamard Transform,
Number Theoretic Transform, Zeta Transform, Mobius Transform
"""
from __future__ import print_function, division

from sympy.core import S, Symbol, sympify
from sympy.core.compatibility import as_int, range
from sympy.core.function import expand, expand_mul
from sympy.core.numbers import pi, I
from sympy.functions.elementary.exponential import exp
from sympy.ntheory import isprime, primitive_root
from sympy.utilities.iterables import ibin

def _fourier_transform(seq, inverse=False):
    """
    Performs the DFT in complex field, by padding zeros as
    radix 2 FFT, requires the number of sample points to be
    a power of 2.
    """

    if not hasattr(seq, '__iter__'):
        raise TypeError("Expected a sequence of numeric coefficients " +
                        "for Fourier Transform")

    a = sympify(seq)
    if any(x.has(Symbol) for x in a):
        raise ValueError("Expected non-symbolic coefficients")

    n = len(a)
    if n < 2:
        return

    while n&(n - 1):
        n += n&-n

    a += [S.Zero]*(n - len(a))
    b = n.bit_length() - 1
    for i in range(1, n):
        j = int(ibin(i, b, str=True)[::-1], 2)
        if i < j:
            a[i], a[j] = a[j], a[i]

    ang = -2*pi/n if inverse else 2*pi/n
    w = [exp(I*ang*i).expand(complex=True) for i in range(n // 2)]

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


def fft(seq):
    return _fourier_transform(seq)


def ifft(seq):
    return _fourier_transform(seq, inverse=True)
