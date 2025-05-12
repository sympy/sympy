from __future__ import annotations

from sympy.core.singleton import S
from sympy.core.symbol import Symbol
from sympy.polys.polytools import lcm
from sympy.utilities import public
from sympy import Expr, Poly
from typing import Generator
from sympy.polys.polytools import gcdex_steps

@public
def approximants(l, X=Symbol('x'), simplify=False):
    """
    Return a generator for consecutive Pade approximants for a series.
    It can also be used for computing the rational generating function of a
    series when possible, since the last approximant returned by the generator
    will be the generating function (if any).

    Explanation
    ===========

    The input list can contain more complex expressions than integer or
    rational numbers; symbols may also be involved in the computation. An
    example below show how to compute the generating function of the whole
    Pascal triangle.

    The generator can be asked to apply the sympy.simplify function on each
    generated term, which will make the computation slower; however it may be
    useful when symbols are involved in the expressions.

    Examples
    ========

    >>> from sympy.series import approximants
    >>> from sympy import lucas, fibonacci, symbols, binomial
    >>> g = [lucas(k) for k in range(16)]
    >>> [e for e in approximants(g)]
    [2, -4/(x - 2), (5*x - 2)/(3*x - 1), (x - 2)/(x**2 + x - 1)]

    >>> h = [fibonacci(k) for k in range(16)]
    >>> [e for e in approximants(h)]
    [x, -x/(x - 1), (x**2 - x)/(2*x - 1), -x/(x**2 + x - 1)]

    >>> x, t = symbols("x,t")
    >>> p=[sum(binomial(k,i)*x**i for i in range(k+1)) for k in range(16)]
    >>> y = approximants(p, t)
    >>> for k in range(3): print(next(y))
    1
    (x + 1)/((-x - 1)*(t*(x + 1) + (x + 1)/(-x - 1)))
    nan

    >>> y = approximants(p, t, simplify=True)
    >>> for k in range(3): print(next(y))
    1
    -1/(t*(x + 1) - 1)
    nan

    See Also
    ========

    sympy.concrete.guess.guess_generating_function_rational
    sympy.series.approximants.pade_approximant
    sympy.series.approximants.pade_approximants
    """
    from sympy.simplify import simplify as simp
    from sympy.simplify.radsimp import denom
    p1, q1 = [S.One], [S.Zero]
    p2, q2 = [S.Zero], [S.One]
    while len(l):
        b = 0
        while l[b]==0:
            b += 1
            if b == len(l):
                return
        m = [S.One/l[b]]
        for k in range(b+1, len(l)):
            s = 0
            for j in range(b, k):
                s -= l[j+1] * m[b-j-1]
            m.append(s/l[b])
        l = m
        a, l[0] = l[0], 0
        p = [0] * max(len(p2), b+len(p1))
        q = [0] * max(len(q2), b+len(q1))
        for k in range(len(p2)):
            p[k] = a*p2[k]
        for k in range(b, b+len(p1)):
            p[k] += p1[k-b]
        for k in range(len(q2)):
            q[k] = a*q2[k]
        for k in range(b, b+len(q1)):
            q[k] += q1[k-b]
        while p[-1]==0: p.pop()
        while q[-1]==0: q.pop()
        p1, p2 = p2, p
        q1, q2 = q2, q

        # yield result
        c = 1
        for x in p:
            c = lcm(c, denom(x))
        for x in q:
            c = lcm(c, denom(x))
        out = ( sum(c*e*X**k for k, e in enumerate(p))
              / sum(c*e*X**k for k, e in enumerate(q)) )
        if simplify:
            yield(simp(out))
        else:
            yield out
    return


@public
def pade_approximants(
    f: Poly, order: int
) -> Generator[tuple[Poly, Poly], None, None]:
    """
    Returns all pade approximants of the expression `f` of the desired order.

    Description
    ===========

    Computes all pade approximants `f` of order `[m/n]` with `m + n =` order.

    Parameters
    ==========

    f : Poly
        The polynomial to approximate, must be univariate
    order : int
        The desired order of the pade approximants

    Returns
    =======

    pade approximants : Generator[tuple[Poly, Poly], None, None]
        Generator for (numerator, denominator) pairs of polynomials
        representing the numerator and denominators of the approximants.

    Examples
    ========

    >>> from sympy.abc import x
    >>> from sympy import Poly
    >>> from sympy.series.approximants import pade_approximants

    >>> exp_taylor_poly = Poly(1 + x + x**2/2 + x**3/6 + x**4/24, x, domain='QQ')
    >>> pade_exp = pade_approximants(exp_taylor_poly, 2)
    >>> for num, denom in pade_exp: print(num, denom)
    Poly(1/2*x**2 + x + 1, x, domain='QQ') Poly(1, x, domain='QQ')
    Poly(2*x + 4, x, domain='QQ') Poly(-2*x + 4, x, domain='QQ')
    Poly(1, x, domain='QQ') Poly(1/2*x**2 - x + 1, x, domain='QQ')

    To get pade approximants, divide the numerators by the denominators

    >>> sin_taylor_poly = Poly(x - x**3/6 + x**5/120, domain='QQ')
    >>> pade_sin = pade_approximants(sin_taylor_poly, 5)
    >>> for num, denom in pade_sin: print(num/denom)
    x**5/120 - x**3/6 + x
    -(20*x**4 - 120*x**2)/(120*x)
    (-7*x**3/60 + x)/(x**2/20 + 1)
    360*x**2/(7*(60*x**3/7 + 360*x/7))
    x/(7*x**4/360 + x**2/6 + 1)

    See also
    ========

    sympy.series.approximants.approximants
    sympy.series.approximants.pade_approximant_gcdex
    sympy.series.approximants.pade_approximant
    """

    if order < 0:
        raise ValueError("'order' must be a non-negative integer")

    if not f.is_univariate:
        raise ValueError("f must be a univariate polynomial")

    # truncate f so f.degree() <= order
    if f.degree() > order:
        f = f.slice(0, order+1)

    x = f.gens[0]

    remainder_monomial = Poly(x ** (order + 1), x)

    yield f, f.one

    euclidean_algorithm_steps = gcdex_steps(remainder_monomial, f)

    for _, t, r in euclidean_algorithm_steps:
        yield r, t


@public
def pade_approximant_gcdex(
    f: Poly, m: int, n: int | None = None
    ) -> tuple[Poly, Poly] | tuple[None, None]:
    """
    `[m/n]` pade approximant of `f` around `x=0`. Computed using the extended euclidean algorithm.

    Description
    ===========

    For a polynomial `f` and integers `m` and `n`, the `[m/n]` pade approximant
    of `f` is a rational function `p(x)/q(x)`, where `p(x)` and `q(x)` are
    polynomials of degree at most `m` and `n` respectively. The pade
    approximation, when it exists, is such that the taylor series of `p(x)/q(x)`
    around `x=0` matches `f` up to order `(m + n)` in `x`.

    Parameters
    ==========

    f : Poly
        The polynomial to approximate, must be univariate
    m : int
        The degree of the numerator polynomial
    n : int | None
        The degree of the denominator polynomial. If None, `n = m`

    Returns
    =======

    pade approximant : tuple[Poly, Poly]
        returns `p(x), q(x)`, where `p(x)/q(x)` is the pade approximant

    Examples
    ========

    >>> import sympy as sp
    >>> from sympy.series.approximants import pade_approximant_gcdex
    >>> x = sp.symbols('x')

    >>> exp_taylor_poly = sp.Poly(x**4/24 + x**3/6 + x**2/2 + x + 1, x, domain='QQ')
    >>> num_exp, denom_exp = pade_approximant_gcdex(exp_taylor_poly, 2)
    >>> print(num_exp, denom_exp)
    Poly(1/4*x**2 + 3/2*x + 3, x, domain='QQ') Poly(1/4*x**2 - 3/2*x + 3, x, domain='QQ')

    The numerators and denominators of an `[m/n]` pade approximant do not
    necessarily have order `m` and `n` respectively.

    >>> sin_taylor_poly = sp.Poly(x**5/120 - x**3/6 + x, x, domain='QQ')
    >>> num_sin, denom_sin = pade_approximant_gcdex(sin_taylor_poly, 1, 3)
    >>> print(num_sin, denom_sin)
    Poly(36*x, x, domain='QQ') Poly(6*x**2 + 36, x, domain='QQ')

    The `[m/0]` pade approximant is the same a truncating to order `m`

    >>> poly = sp.Poly(4*x**4 + 3*x**3 + 2*x**2 + x + 1, x, domain='QQ')
    >>> num_cos, denom_cos = pade_approximant_gcdex(poly, 3, 0)
    >>> print(num_cos, denom_cos)
    Poly(3*x**3 + 2*x**2 + x + 1, x, domain='QQ') Poly(1, x, domain='QQ')

    The pade approximant does not always exist, for example, the `[1/1]` approximant for `cos(x)`

    >>> cos_taylor_poly = sp.Poly(1 - x**2/2, x, domain='QQ')
    >>> num_cos11, denom_cos11 = pade_approximant_gcdex(cos_taylor_poly, 1)
    >>> print(num_cos11, denom_cos11)
    Poly(2*x, x, domain='QQ') Poly(2*x, x, domain='QQ')

    `num_cos11/denom_cos11 == 1`, and `1` does not match the taylor expansion of `cos(x)`
    to second order in `x`.

    See also
    ========

    sympy.series.approximants.approximants
    sympy.series.approximants.pade_approximant
    sympy.series.approximants.pade_approximants
    """
    if n is None:
        n = m

    if m < 0 or n < 0:
        raise ValueError("m and n must be non-negative integers")

    if not f.is_univariate:
        raise ValueError("f must be a univariate polynomial")

    min_degree = f.EM()[0]
    if min_degree > m:
        raise ValueError(f"polynomial has zero of order {min_degree}, " \
        f"which is greater than {m}, the requested order of the numerator.")

    approximants = pade_approximants(f, m + n)

    numerator, denominator = next(approximants)
    for next_numerator, next_denominator in approximants:
        if next_numerator.degree() < m:
            return numerator, denominator

        numerator, denominator = next_numerator, next_denominator

    if numerator.degree() == m:
        return numerator, denominator

    # if it hasn't returned yet, there is no [m/n] pade approximant
    return None, None


@public
def pade_approximant(
    f: Expr, x: Expr, m: int, n: int | None = None
    ) -> Expr:
    """
    `[m/n]` pade approximant of `f` around `x=0`.

    Description
    ===========

    For a function `f` and integers `m` and `n`, the `[m/n]` pade approximant
    of `f` is the rational function `p(x)/q(x)`, where `p(x)` and `q(x)` are
    polynomials of degree at most `m` and `n` respectively. The pade
    approximation, when it exists, is such that the taylor series of `p(x)/q(x)`
    around `x=0`matches the taylor series of `f(x)` up to at least order `(m + n)`
    in `x`.

    Parameters
    ==========

    f : Expr
        The function to approximate
    x : Expr
        The variable of the function
    m : int
        The degree of the numerator polynomial
    n : int | None
        The degree of the denominator polynomial. If None, `n = m`

    Returns
    =======

    pade approximant : tuple[Poly, Poly]
        returns p(x)

    Examples
    ========

    >>> import sympy as sp
    >>> from sympy.series.approximants import pade_approximant
    >>> x = sp.symbols('x')

    >>> pade_exp = pade_approximant(sp.exp(x), x, 2)
    >>> pade_exp
    (x**2/4 + 3*x/2 + 3)/(x**2/4 - 3*x/2 + 3)

    The numerators and denominators of an `[m/n]` pade approximant do not
    necessarily have order `m` and `n` respectively.

    >>> pade_sin = pade_approximant(sp.sin(x), x, 1, 3)
    >>> pade_sin
    36*x/(6*x**2 + 36)

    The `[m/0]` pade approximant is the mth order taylor polynomial of `f`

    >>> pade_cos = pade_approximant(sp.cos(x), x, 4, 0)
    >>> pade_cos
    x**4/24 - x**2/2 + 1

    The pade approximant does not always exist, for example, the `[1/1]` approximant for `cos(x)`

    >>> pade_cos11 = pade_approximant(sp.cos(x), x, 1)
    >>> pade_cos11
    1

    `1` does not match the taylor expansion of `cos(x)` to second order in `x`.

    See also
    ========

    sympy.series.approximants.approximants
    sympy.series.approximants.pade_approximant_gcdex
    sympy.series.approximants.pade_approximants
    """
    if n is None:
        n = m

    if m < 0 or n < 0:
        raise ValueError("m and n must be non-negative integers")

    f_taylor_poly = f.series(x, n=m + n + 1).removeO().as_poly(x, field=True)

    # TODO: if other methods for computing pade approximants are implemented,
    # this function should be modified to automatically select the best method.
    numerator, denominator = pade_approximant_gcdex(f_taylor_poly, m, n)

    return numerator/denominator
