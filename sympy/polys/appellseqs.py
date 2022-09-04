r"""
Efficient functions for generating Appell sequences.

An Appell sequence is a zero-indexed sequence of polynomials `p_i(x)`
satisfying `p_{i+1}'(x)=(i+1)p_i(x)` for all `i`. This definition leads
to the following iterative algorithm:

.. math :: p_i(x) = i \int_0^x p_{i-1}(t)\,dt + c_i

where the new constant terms `c_i` are usually determined from the
just-evaluated integral and `i`.

Appell sequences satisfy the following identity from umbral calculus:

.. math :: p_n(x+y) = \sum_{k=0}^n \binom{n}{k} p_k(x) y^{n-k}

References
==========

.. [1] https://en.wikipedia.org/wiki/Appell_sequence
.. [2] Peter Luschny, "An introduction to the Bernoulli function",
       https://arxiv.org/abs/2009.06743
"""
from sympy.core.symbol import Dummy
from sympy.polys.densearith import dup_mul_ground, dup_sub_ground
from sympy.polys.densetools import dup_eval, dup_integrate
from sympy.polys.domains import ZZ, QQ
from sympy.polys.polyclasses import DMP
from sympy.polys.polytools import Poly, PurePoly
from sympy.utilities import public

def dup_bernoulli(n, K):
    """Low-level implementation of Bernoulli polynomials."""
    seq = [[1], [1, K(-1,2)]]
    for i in range(2, n+1):
        q = dup_integrate(dup_mul_ground(seq[-1], i, K), 1, K)
        if i&1 ^ 1:
            q = dup_sub_ground(q, dup_eval(q, K(1,2), K) * K(1<<(i-1), (1<<i)-1), K)
        seq.append(q)
    return seq[n]

@public
def bernoulli_poly(n, x=None, polys=False):
    r"""Generates the Bernoulli polynomial `\operatorname{B}_n(x)`.

    `\operatorname{B}_n(x)` is the unique polynomial satisfying

    .. math :: \int_{x}^{x+1} B_n(t) \,dt = x^n.

    Based on this, we have for nonnegative integer `s` and integer
    `a` and `b`

    .. math :: \sum_{k=a}^{b} k^s = \frac{\operatorname{B}_{s+1}(b+1) -
            \operatorname{B}_{s+1}(a)}{s+1}

    which is related to Jakob Bernoulli's original motivation for introducing
    the Bernoulli numbers, the values of these polynomials at `x = 1`.

    Examples
    ========

    >>> from sympy import summation
    >>> from sympy.abc import x
    >>> from sympy.polys import bernoulli_poly
    >>> bernoulli_poly(5, x)
    x**5 - 5*x**4/2 + 5*x**3/3 - x/6

    >>> def psum(p, a, b):
    ...     return (bernoulli_poly(p+1,b+1) - bernoulli_poly(p+1,a)) / (p+1)
    >>> psum(4, -6, 27)
    3144337
    >>> summation(x**4, (x, -6, 27))
    3144337

    >>> psum(1, 1, x).factor()
    x*(x + 1)/2
    >>> psum(2, 1, x).factor()
    x*(x + 1)*(2*x + 1)/6
    >>> psum(3, 1, x).factor()
    x**2*(x + 1)**2/4

    Parameters
    ==========

    n : int
        Degree of the polynomial.
    x : optional
    polys : bool, optional
        If True, return a Poly, otherwise (default) return an expression.

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Bernoulli_polynomials
    """
    if n < 0:
        raise ValueError("Cannot generate Bernoulli polynomial of degree %s" % n)
    poly = DMP(dup_bernoulli(int(n), QQ), QQ)
    if x is not None:
        poly = Poly.new(poly, x)
    else:
        poly = PurePoly.new(poly, Dummy('x'))
    return poly if polys else poly.as_expr()


def dup_bernoulli_c(n, K):
    """Low-level implementation of central Bernoulli polynomials."""
    seq = [[1], [1, 0]]
    for i in range(2, n+1):
        q = dup_integrate(dup_mul_ground(seq[-1], i, K), 1, K)
        if i&1 ^ 1:
            q = dup_sub_ground(q, dup_eval(q, 1, K) * K((1<<(i-1))-1, (1<<i)-1), K)
        seq.append(q)
    return seq[n]

@public
def bernoulli_c_poly(n, x=None, polys=False):
    r"""Generates the central Bernoulli polynomial `\operatorname{B}_n^c(x)`.

    These are scaled and shifted versions of the plain Bernoulli polynomials,
    done in such a way that `operatorname{B}_n^c(x)` is an even or odd function
    for even or odd `n` respectively:

    .. math :: \operatorname{B}_n^c(x) = 2^n \operatorname{B}_n
            \left(\frac{x+1}{2}\right)

    Parameters
    ==========

    n : int
        Degree of the polynomial.
    x : optional
    polys : bool, optional
        If True, return a Poly, otherwise (default) return an expression.
    """
    if n < 0:
        raise ValueError("Cannot generate central Bernoulli polynomial of degree %s" % n)
    poly = DMP(dup_bernoulli_c(int(n), QQ), QQ)
    if x is not None:
        poly = Poly.new(poly, x)
    else:
        poly = PurePoly.new(poly, Dummy('x'))
    return poly if polys else poly.as_expr()


def dup_genocchi(n, K):
    """Low-level implementation of Genocchi polynomials."""
    seq = [[0], [-1]]
    for i in range(2, n+1):
        q = dup_integrate(dup_mul_ground(seq[-1], i, K), 1, K)
        if i&1 ^ 1:
            q = dup_sub_ground(q, dup_eval(q, 1, K) >> 1, K)
        seq.append(q)
    return seq[n]

@public
def genocchi_poly(n, x=None, polys=False):
    """Generates the Genocchi polynomial `\operatorname{G}_n(x)`.

    `\operatorname{G}_n(x)` is twice the difference between the plain and
    central Bernoulli polynomials, so has degree `n-1`:

    .. math :: \operatorname{G}_n(x) = 2 (\operatorname{B}_n(x) -
            \operatorname{B}_n^c(x))

    The factor of 2 in the definition endows `\operatorname{G}_n(x)` with
    integer coefficients.

    Parameters
    ==========

    n : int
        Degree of the polynomial plus one.
    x : optional
    polys : bool, optional
        If True, return a Poly, otherwise (default) return an expression.
    """
    if n < 0:
        raise ValueError("Cannot generate Genocchi polynomial of degree %s" % (n-1))
    poly = DMP(dup_genocchi(int(n), ZZ), ZZ)
    if x is not None:
        poly = Poly.new(poly, x)
    else:
        poly = PurePoly.new(poly, Dummy('x'))
    return poly if polys else poly.as_expr()


def dup_euler(n, K):
    """Low-level implementation of Euler polynomials."""
    seq = [[1], [1, K(-1,2)]]
    for i in range(2, n+1):
        q = dup_integrate(dup_mul_ground(seq[-1], i, K), 1, K)
        if i&1:
            q = dup_sub_ground(q, dup_eval(q, 1, K) / 2, K)
        seq.append(q)
    return seq[n]

@public
def euler_poly(n, x=None, polys=False):
    """Generates the Euler polynomial `\operatorname{E}_n(x)`.

    Parameters
    ==========

    n : int
        Degree of the polynomial.
    x : optional
    polys : bool, optional
        If True, return a Poly, otherwise (default) return an expression.
    """
    if n < 0:
        raise ValueError("Cannot generate Euler polynomial of degree %s" % n)
    poly = DMP(dup_euler(int(n), QQ), QQ)
    if x is not None:
        poly = Poly.new(poly, x)
    else:
        poly = PurePoly.new(poly, Dummy('x'))
    return poly if polys else poly.as_expr()


def dup_andre(n, K):
    """Low-level implementation of Andre polynomials."""
    seq = [[1], [1, 0]]
    for i in range(2, n+1):
        q = dup_integrate(dup_mul_ground(seq[-1], i, K), 1, K)
        if i&1 ^ 1:
            q = dup_sub_ground(q, dup_eval(q, 1, K), K)
        seq.append(q)
    return seq[n]

@public
def andre_poly(n, x=None, polys=False):
    """Generates the Andre polynomial of degree `n` in `x`.

    Luschny calls these the *Swiss-knife polynomials*.

    Parameters
    ==========

    n : int
        Degree of the polynomial.
    x : optional
    polys : bool, optional
        If True, return a Poly, otherwise (default) return an expression.

    References
    ==========

    .. [1] Peter Luschny, "An introduction to the Bernoulli function",
           https://arxiv.org/abs/2009.06743
    """
    if n < 0:
        raise ValueError("Cannot generate Andre polynomial of degree %s" % n)
    poly = DMP(dup_andre(int(n), ZZ), ZZ)
    if x is not None:
        poly = Poly.new(poly, x)
    else:
        poly = PurePoly.new(poly, Dummy('x'))
    return poly if polys else poly.as_expr()
