# -*- coding: utf-8 -*-

"""
Discrete Fourier Transform, Number Theoretic Transform,
Walsh Hadamard Transform, Mobius Transform
"""
from __future__ import print_function, division, unicode_literals

from sympy.core import S, Symbol, sympify
from sympy.core.compatibility import as_int, range, iterable
from sympy.core.function import expand_mul
from sympy.core.numbers import pi, I
from sympy.functions.elementary.trigonometric import sin, cos
from sympy.ntheory import isprime, primitive_root
from sympy.utilities.iterables import ibin


#----------------------------------------------------------------------------#
#                                                                            #
#                         Discrete Fourier Transform                         #
#                                                                            #
#----------------------------------------------------------------------------#

def _fourier_transform(seq, dps, inverse=False, fourier_params=(1, 1)):
    """Utility function for the Discrete Fourier Transform"""

    if not iterable(seq):
        raise TypeError("Expected a sequence of numeric coefficients "
                        "for Fourier Transform")

    a = [sympify(arg) for arg in seq]

    n = len(a)
    if n < 2:
        return a

    (alpha, beta) = sympify(fourier_params)

    b = n.bit_length() - 1
    if n&(n - 1): # not a power of 2
        b += 1
        n = 2**b

    a += [S.Zero]*(n - len(a))
    for i in range(1, n):
        j = int(ibin(i, b, str=True)[::-1], 2)
        if i < j:
            a[i], a[j] = a[j], a[i]

    ang = -2*pi*beta/n if inverse else 2*pi*beta/n

    if dps is not None:
        ang = ang.evalf(dps + 2)

    w = [cos(ang*i) + I*sin(ang*i) for i in range(n // 2)]

    h = 2
    while h <= n:
        hf, ut = h // 2, n // h
        for i in range(0, n, h):
            for j in range(hf):
                u, v = a[i + j], expand_mul(a[i + j + hf]*w[ut * j])
                a[i + j], a[i + j + hf] = u + v, u - v
        h *= 2

    n = sympify(n)
    denom = n**((S(1)+alpha)/S(2)) if inverse else n**((S(1)-alpha)/S(2))
    a = [(x/denom).evalf(dps) for x in a] if dps else [x/denom for x in a]

    return a


def fft(seq, dps=None, fourier_params=(1, 1)):
    r"""Performs the Fast Fourier Transform (**FFT**).

    The sequence is automatically padded to the right with zeros, as the
    *radix-2 FFT* requires the number of sample points to be a power of 2.

    Parameters
    ==========

    seq : iterable
        The sequence on which **DFT** is to be applied.

    dps : Integer
        Specifies the number of decimal digits for precision.

    fourier_params : (Expr, Expr)
        Specifies the formula used for discrete Fourier transformation.

        Different libraries or CAS may use different conventions which
        should be taken some precaution into.
        The conventional settings are :

        SymPy : (1, 1)

        NumPy, MATLAB : (1, -1)

        Mathematica : (0, 1)

        Maple : (0, -1)

        You may also pass the symbols to get the full symbolic result
        for the FFT.
        However, we also do not block you from passing the numbers which
        may be invalid to be consider a fourier transform, or out of
        convention.

        Explicit formula is described in the notes section below.

    Notes
    =====

    The formula for computing Fourier transformation may vary across libraries,
    and caution should be taken.

    SymPy's default implementation uses the formula

    `A[k] = \sum_{n=0}^{N-1} a[n] \omega_{N}^{k n}` for ``fft``

    and `a[k] =\frac{1}{N} \sum_{n=0}^{N-1} A[n] \omega_{N}^{k n}`
    for ``ifft``

    where `\omega_{N} = e^{2 \pi i/N}`.

    While it can be generalized as

    `A[k] = \frac{1}{\sqrt{N^{1-\alpha}}}
    \sum_{n=0}^{N-1} a[n] \omega_{N}^{k n}`
    for ``fft``

    and
    `a[k] =\frac{1}{\sqrt{N^{1+\alpha}}}
    \sum_{n=0}^{N-1} A[n] \omega_{N}^{k n}`
    for ``ifft``

    where `\omega_{N} = e^{2 \pi i \beta/N}`

    while `(\alpha,\beta)` are ``fourier_params``.

    However, there may be certain restrictions which can be placed onto
    `\alpha` and `\beta`, if you want to make your FFT scheme valid for
    most of the context it is being used.
    `\alpha` can work fine with any real numbers because it may only
    scale the output result.
    However, specifying `\beta` into zero or any number which is not
    relatively prime to the length of the list may make the fourier
    matrix singular, which may result in no theoritical way of inversing
    the transform by another transform.
    Specifying `\beta` into non-integer may result not make the fourier
    matrix singular in the loose definition, but cooley-tukey scheme
    may break and the algorithm may not yield the same result as the
    corresponding DFT matrix with the same parameter.
    These behaviors may easily be observed by experimenting with our
    library, and the conventional reason that most of the textbooks or
    numerical libraries are setting `\beta` into `1` or `-1` is that
    it is the most practical way of specifying the relatively-prime
    integer without having to worry about cases where the length of the
    the list to transform may vary.

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

    >>> ifft([1, 7, 3, 4], dps=15)
    [3.75, -0.5 - 0.75*I, -1.75, -0.5 + 0.75*I]
    >>> fft(_)
    [1.0, 7.0, 3.0, 4.0]

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
    .. [2] http://mathworld.wolfram.com/FastFourierTransform.html

    """

    return _fourier_transform(seq, dps=dps,
        fourier_params=fourier_params)


def ifft(seq, dps=None, fourier_params=(1, 1)):
    return _fourier_transform(seq, dps=dps, inverse=True,
        fourier_params=fourier_params)

ifft.__doc__ = fft.__doc__


#----------------------------------------------------------------------------#
#                                                                            #
#                         Number Theoretic Transform                         #
#                                                                            #
#----------------------------------------------------------------------------#

def _number_theoretic_transform(seq, prime, inverse=False):
    """Utility function for the Number Theoretic Transform"""

    if not iterable(seq):
        raise TypeError("Expected a sequence of integer coefficients "
                        "for Number Theoretic Transform")

    p = as_int(prime)
    if not isprime(p):
        raise ValueError("Expected prime modulus for "
                        "Number Theoretic Transform")

    a = [as_int(x) % p for x in seq]

    n = len(a)
    if n < 1:
        return a

    b = n.bit_length() - 1
    if n&(n - 1):
        b += 1
        n = 2**b

    if (p - 1) % n:
        raise ValueError("Expected prime modulus of the form (m*2**k + 1)")

    a += [0]*(n - len(a))
    for i in range(1, n):
        j = int(ibin(i, b, str=True)[::-1], 2)
        if i < j:
            a[i], a[j] = a[j], a[i]

    pr = primitive_root(p)

    rt = pow(pr, (p - 1) // n, p)
    if inverse:
        rt = pow(rt, p - 2, p)

    w = [1]*(n // 2)
    for i in range(1, n // 2):
        w[i] = w[i - 1]*rt % p

    h = 2
    while h <= n:
        hf, ut = h // 2, n // h
        for i in range(0, n, h):
            for j in range(hf):
                u, v = a[i + j], a[i + j + hf]*w[ut * j]
                a[i + j], a[i + j + hf] = (u + v) % p, (u - v) % p
        h *= 2

    if inverse:
        rv = pow(n, p - 2, p)
        a = [x*rv % p for x in a]

    return a


def ntt(seq, prime):
    r"""
    Performs the Number Theoretic Transform (**NTT**), which specializes the
    Discrete Fourier Transform (**DFT**) over quotient ring `Z/pZ` for prime
    `p` instead of complex numbers `C`.

    The sequence is automatically padded to the right with zeros, as the
    *radix-2 NTT* requires the number of sample points to be a power of 2.

    Parameters
    ==========

    seq : iterable
        The sequence on which **DFT** is to be applied.
    prime : Integer
        Prime modulus of the form `(m 2^k + 1)` to be used for performing
        **NTT** on the sequence.

    Examples
    ========

    >>> from sympy import ntt, intt
    >>> ntt([1, 2, 3, 4], prime=3*2**8 + 1)
    [10, 643, 767, 122]
    >>> intt(_, 3*2**8 + 1)
    [1, 2, 3, 4]
    >>> intt([1, 2, 3, 4], prime=3*2**8 + 1)
    [387, 415, 384, 353]
    >>> ntt(_, prime=3*2**8 + 1)
    [1, 2, 3, 4]

    References
    ==========

    .. [1] http://www.apfloat.org/ntt.html
    .. [2] http://mathworld.wolfram.com/NumberTheoreticTransform.html
    .. [3] https://en.wikipedia.org/wiki/Discrete_Fourier_transform_(general%29

    """

    return _number_theoretic_transform(seq, prime=prime)


def intt(seq, prime):
    return _number_theoretic_transform(seq, prime=prime, inverse=True)

intt.__doc__ = ntt.__doc__


#----------------------------------------------------------------------------#
#                                                                            #
#                          Walsh Hadamard Transform                          #
#                                                                            #
#----------------------------------------------------------------------------#

def _walsh_hadamard_transform(seq, inverse=False):
    """Utility function for the Walsh Hadamard Transform"""

    if not iterable(seq):
        raise TypeError("Expected a sequence of coefficients "
                        "for Walsh Hadamard Transform")

    a = [sympify(arg) for arg in seq]
    n = len(a)
    if n < 2:
        return a

    if n&(n - 1):
        n = 2**n.bit_length()

    a += [S.Zero]*(n - len(a))
    h = 2
    while h <= n:
        hf, ut = h // 2, n // h
        for i in range(0, n, h):
            for j in range(hf):
                u, v = a[i + j], a[i + j + hf]
                a[i + j], a[i + j + hf] = u + v, u - v
        h *= 2

    if inverse:
        a = [x/n for x in a]

    return a


def fwht(seq):
    r"""
    Performs the Walsh Hadamard Transform (**WHT**), and uses Hadamard
    ordering for the sequence.

    The sequence is automatically padded to the right with zeros, as the
    *radix-2 FWHT* requires the number of sample points to be a power of 2.

    Parameters
    ==========

    seq : iterable
        The sequence on which WHT is to be applied.

    Examples
    ========

    >>> from sympy import fwht, ifwht
    >>> fwht([4, 2, 2, 0, 0, 2, -2, 0])
    [8, 0, 8, 0, 8, 8, 0, 0]
    >>> ifwht(_)
    [4, 2, 2, 0, 0, 2, -2, 0]

    >>> ifwht([19, -1, 11, -9, -7, 13, -15, 5])
    [2, 0, 4, 0, 3, 10, 0, 0]
    >>> fwht(_)
    [19, -1, 11, -9, -7, 13, -15, 5]

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Hadamard_transform
    .. [2] https://en.wikipedia.org/wiki/Fast_Walsh%E2%80%93Hadamard_transform

    """

    return _walsh_hadamard_transform(seq)


def ifwht(seq):
    return _walsh_hadamard_transform(seq, inverse=True)

ifwht.__doc__ = fwht.__doc__


#----------------------------------------------------------------------------#
#                                                                            #
#                    Möbius Transform for Subset Lattice                     #
#                                                                            #
#----------------------------------------------------------------------------#

def _mobius_transform(seq, sgn, subset):
    r"""Utility function for performing Möbius Transform using
    Yate's Dynamic Programming method"""

    if not iterable(seq):
        raise TypeError("Expected a sequence of coefficients")

    a = [sympify(arg) for arg in seq]

    n = len(a)
    if n < 2:
        return a

    if n&(n - 1):
        n = 2**n.bit_length()

    a += [S.Zero]*(n - len(a))

    if subset:
        i = 1
        while i < n:
            for j in range(n):
                if j & i:
                    a[j] += sgn*a[j ^ i]
            i *= 2

    else:
        i = 1
        while i < n:
            for j in range(n):
                if j & i:
                    continue
                a[j] += sgn*a[j ^ i]
            i *= 2

    return a


def mobius_transform(seq, subset=True):
    r"""
    Performs the Möbius Transform for subset lattice with indices of
    sequence as bitmasks.

    The indices of each argument, considered as bit strings, correspond
    to subsets of a finite set.

    The sequence is automatically padded to the right with zeros, as the
    definition of subset/superset based on bitmasks (indices) requires
    the size of sequence to be a power of 2.

    Parameters
    ==========

    seq : iterable
        The sequence on which Möbius Transform is to be applied.
    subset : bool
        Specifies if Möbius Transform is applied by enumerating subsets
        or supersets of the given set.

    Examples
    ========

    >>> from sympy import symbols
    >>> from sympy import mobius_transform, inverse_mobius_transform
    >>> x, y, z = symbols('x y z')

    >>> mobius_transform([x, y, z])
    [x, x + y, x + z, x + y + z]
    >>> inverse_mobius_transform(_)
    [x, y, z, 0]

    >>> mobius_transform([x, y, z], subset=False)
    [x + y + z, y, z, 0]
    >>> inverse_mobius_transform(_, subset=False)
    [x, y, z, 0]

    >>> mobius_transform([1, 2, 3, 4])
    [1, 3, 4, 10]
    >>> inverse_mobius_transform(_)
    [1, 2, 3, 4]
    >>> mobius_transform([1, 2, 3, 4], subset=False)
    [10, 6, 7, 4]
    >>> inverse_mobius_transform(_, subset=False)
    [1, 2, 3, 4]

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/M%C3%B6bius_inversion_formula
    .. [2] https://people.csail.mit.edu/rrw/presentations/subset-conv.pdf
    .. [3] https://arxiv.org/pdf/1211.0189.pdf

    """

    return _mobius_transform(seq, sgn=+1, subset=subset)

def inverse_mobius_transform(seq, subset=True):
    return _mobius_transform(seq, sgn=-1, subset=subset)

inverse_mobius_transform.__doc__ = mobius_transform.__doc__
