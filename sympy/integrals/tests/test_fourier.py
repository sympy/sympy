from __future__ import annotations
from sympy.integrals.fourier import (
    fourier_transform, inverse_fourier_transform,
    FourierTransform, InverseFourierTransform)
from sympy.core.function import Function
from sympy.core.numbers import pi
from sympy.core.symbol import symbols
from sympy.functions.elementary.complexes import Abs
from sympy.functions.elementary.exponential import exp
from sympy.functions.elementary.miscellaneous import sqrt


def test_fourier_transform():
    FT = fourier_transform
    IFT = inverse_fourier_transform
    a = symbols('a', positive=True)
    x, k = symbols('x k')
    f = Function('f')

    # Test unevaluated form
    assert FT(f(x), x, k) == FourierTransform(f(x), x, k)
    assert IFT(f(k), k, x) == InverseFourierTransform(f(k), k, x)

    # Classical rules
    assert FT(exp(-x**2), x, k) == sqrt(pi)*exp(-pi**2*k**2)
    assert FT(exp(-x**2), x, k, noconds=False) == \
        (sqrt(pi)*exp(-pi**2*k**2), True)
    assert FT(exp(-a*x**2), x, k) == \
        sqrt(pi)*exp(-pi**2*k**2/a)/sqrt(a)
    assert FT(exp(-a*x**2), x, k, noconds=False) == \
        (sqrt(pi)*exp(-pi**2*k**2/a)/sqrt(a), True)
    assert FT(exp(-a*Abs(x)), x, k) == 2*a/(a**2 + 4*pi**2*k**2)
    assert FT(1/(a**2 + x**2), x, k) == pi*exp(-2*pi*a*Abs(k))/a
    assert FT(1/(1 + x**2), x, k) == pi*exp(-2*pi*Abs(k))

    # Linearity and constant coefficients
    assert FT(2*exp(-x**2), x, k) == 2*sqrt(pi)*exp(-pi**2*k**2)
    assert (FT(exp(-x**2) + exp(-a*Abs(x)), x, k) ==
            sqrt(pi)*exp(-pi**2*k**2) + 2*a/(a**2 + 4*pi**2*k**2))

    # doit on unevaluated objects
    assert (2*FourierTransform(exp(-x**2), x, k)).doit() == \
        2*sqrt(pi)*exp(-pi**2*k**2)
    assert FourierTransform(f(x), x, k).doit(noconds=False) == \
        (FourierTransform(f(x), x, k), True)

    # The inverse transform is unchanged so far (integration backend)
    assert IFT(sqrt(pi/a)*exp(-(pi*k)**2/a), k, x) == exp(-a*x**2)
