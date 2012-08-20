from sympy.core import S, Dummy
from sympy.polys.orthopolys import legendre_poly, laguerre_poly
from sympy.polys.rootoftools import RootOf

def gauss_legendre(n, n_digits):
    r"""
    Computes the Gauss-Legendre quadrature [1] points and weights.

    Parameters
    ==========

    n : the order of quadrature

    n_digits : number of significant digits of the points and weights to
               return

    Returns
    =======

    (x, w) : the ``x`` and ``w`` are lists of points and weights as Floats

    The Gauss-Legendre quadrature approximates the integral:

    .. math::

        \int_{-1}^1 f(x)\,dx \approx \sum_{i=1}^n w_i f(x_i)

    The points `x_i` and weights `w_i` are returned as ``(x, w)`` tuple of
    lists.

    Examples
    ========

    >>> from sympy.integrals.quadrature import gauss_legendre
    >>> x, w = gauss_legendre(3, 5)
    >>> x
    [-0.7746, 0, 0.7746]
    >>> w
    [0.55556, 0.88889, 0.55556]
    >>> x, w = gauss_legendre(4, 5)
    >>> x
    [-0.86114, -0.33998, 0.33998, 0.86114]
    >>> w
    [0.34786, 0.65215, 0.65215, 0.34786]

    [1] http://en.wikipedia.org/wiki/Gaussian_quadrature
    """
    x = Dummy("x")
    p = legendre_poly(n, x, polys=True)
    pd = p.diff(x)
    xi = []
    w  = []
    for r in p.real_roots():
        if isinstance(r, RootOf):
            r = r.eval_rational(S(1)/10**(n_digits+2))
        xi.append(r.n(n_digits))
        w.append((2/((1-r**2) * pd.subs(x, r)**2)).n(n_digits))
    return xi, w

def gauss_laguerre(n, n_digits):
    r"""
    Computes the Gauss-Laguerre quadrature [1] points and weights.

    Parameters
    ==========

    n : the order of quadrature

    n_digits : number of significant digits of the points and weights to
               return

    Returns
    =======

    (x, w) : the ``x`` and ``w`` are lists of points and weights as Floats

    The Gauss-Laguerre quadrature approximates the integral:

    .. math::

        \int_0^{\infty} e^{-x} f(x)\,dx \approx \sum_{i=1}^n w_i f(x_i)

    The points `x_i` and weights `w_i` are returned as ``(x, w)`` tuple of
    lists.

    Examples
    ========

    >>> from sympy.integrals.quadrature import gauss_laguerre
    >>> x, w = gauss_laguerre(3, 5)
    >>> x
    [0.41577, 2.2943, 6.2899]
    >>> w
    [0.71109, 0.27852, 0.010389]
    >>> x, w = gauss_laguerre(6, 5)
    >>> x
    [0.22285, 1.1889, 2.9927, 5.7751, 9.8375, 15.983]
    >>> w
    [0.45896, 0.417, 0.11337, 0.010399, 0.00026102, 8.9855e-7]

    [1] http://en.wikipedia.org/wiki/Gauss%E2%80%93Laguerre_quadrature
    """
    x = Dummy("x")
    p  = laguerre_poly(n, x, polys=True)
    p1 = laguerre_poly(n+1, x, polys=True)
    xi = []
    w  = []
    for r in p.real_roots():
        if isinstance(r, RootOf):
            r = r.eval_rational(S(1)/10**(n_digits+2))
        xi.append(r.n(n_digits))
        w.append((r/((n+1)**2 * p1.subs(x, r)**2)).n(n_digits))
    return xi, w
