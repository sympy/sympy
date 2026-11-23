"""Fourier Transforms"""
from __future__ import annotations
from sympy.core import S
from sympy.integrals.transforms import FourierTypeTransform


class FourierTransform(FourierTypeTransform):
    """
    Class representing unevaluated Fourier transforms.

    For usage of this class, see the :class:`IntegralTransform` docstring.

    For how to compute Fourier transforms, see the :func:`fourier_transform`
    docstring.
    """

    _name = 'Fourier'

    def a(self):
        return 1

    def b(self):
        return -2*S.Pi


def fourier_transform(f, x, k, **hints):
    r"""
    Compute the unitary, ordinary-frequency Fourier transform of ``f``, defined
    as

    .. math:: F(k) = \int_{-\infty}^\infty f(x) e^{-2\pi i x k} \mathrm{d} x.

    Explanation
    ===========

    If the transform cannot be computed in closed form, this
    function returns an unevaluated :class:`FourierTransform` object.

    For other Fourier transform conventions, see the function
    :func:`sympy.integrals.transforms._fourier_type_transform`.

    For a description of possible hints, refer to the docstring of
    :func:`sympy.integrals.transforms.IntegralTransform.doit`.
    Note that for this transform, by default ``noconds=True``.

    Examples
    ========

    >>> from sympy import fourier_transform, exp
    >>> from sympy.abc import x, k
    >>> fourier_transform(exp(-x**2), x, k)
    sqrt(pi)*exp(-pi**2*k**2)
    >>> fourier_transform(exp(-x**2), x, k, noconds=False)
    (sqrt(pi)*exp(-pi**2*k**2), True)

    See Also
    ========

    inverse_fourier_transform
    sine_transform, inverse_sine_transform
    cosine_transform, inverse_cosine_transform
    hankel_transform, inverse_hankel_transform
    mellin_transform, laplace_transform
    """
    return FourierTransform(f, x, k).doit(**hints)


class InverseFourierTransform(FourierTypeTransform):
    """
    Class representing unevaluated inverse Fourier transforms.

    For usage of this class, see the :class:`IntegralTransform` docstring.

    For how to compute inverse Fourier transforms, see the
    :func:`inverse_fourier_transform` docstring.
    """

    _name = 'Inverse Fourier'

    def a(self):
        return 1

    def b(self):
        return 2*S.Pi


def inverse_fourier_transform(F, k, x, **hints):
    r"""
    Compute the unitary, ordinary-frequency inverse Fourier transform of `F`,
    defined as

    .. math:: f(x) = \int_{-\infty}^\infty F(k) e^{2\pi i x k} \mathrm{d} k.

    Explanation
    ===========

    If the transform cannot be computed in closed form, this
    function returns an unevaluated :class:`InverseFourierTransform` object.

    For other Fourier transform conventions, see the function
    :func:`sympy.integrals.transforms._fourier_type_transform`.

    For a description of possible hints, refer to the docstring of
    :func:`sympy.integrals.transforms.IntegralTransform.doit`.
    Note that for this transform, by default ``noconds=True``.

    Examples
    ========

    >>> from sympy import inverse_fourier_transform, exp, sqrt, pi
    >>> from sympy.abc import x, k
    >>> inverse_fourier_transform(sqrt(pi)*exp(-(pi*k)**2), k, x)
    exp(-x**2)
    >>> inverse_fourier_transform(sqrt(pi)*exp(-(pi*k)**2), k, x, noconds=False)
    (exp(-x**2), True)

    See Also
    ========

    fourier_transform
    sine_transform, inverse_sine_transform
    cosine_transform, inverse_cosine_transform
    hankel_transform, inverse_hankel_transform
    mellin_transform, laplace_transform
    """
    return InverseFourierTransform(F, k, x).doit(**hints)
