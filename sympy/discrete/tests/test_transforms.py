from __future__ import print_function, division

from sympy import sqrt
from sympy.core import S, Symbol, I
from sympy.core.compatibility import range
from sympy.discrete import fft, ifft
from sympy.utilities.pytest import raises


def test_fft_ifft():
    assert all(tf(ls) == ls for tf in (fft, ifft) \
                        for ls in ([], [S(5)/3]))

    ls = list(range(6))
    fls = [15, -7*sqrt(2)/2 - 4 - sqrt(2)*I/2 + 2*I, 2 + 3*I,
             -4 + 7*sqrt(2)/2 - 2*I - sqrt(2)*I/2, -3,
             -4 + 7*sqrt(2)/2 + sqrt(2)*I/2 + 2*I,
              2 - 3*I, -7*sqrt(2)/2 - 4 - 2*I + sqrt(2)*I/2]

    assert fft(ls) == fls
    assert ifft(fls) == ls + [S.Zero]*2

    ls = [1 + 2*I, 3 + 4*I, 5 + 6*I]
    ifls = [S(9)/4 + 3*I, -7*I/4, S(3)/4 + I, -2 - I/4]

    assert ifft(ls) == ifls
    assert fft(ifls) == ls + [S.Zero]

    x = Symbol('x', real=True)
    raises(TypeError, lambda: fft(x))
    raises(ValueError, lambda: ifft([x, 2*x, 3*x**2, 4*x**3]))
