from sympy.simplify.radsimp import numer, denom
from sympy.polys.polytools import Poly, lcm
from sympy.solvers.solvers import solve
from sympy.core.symbol import symbols
from sympy.core.mul import prod
from sympy.concrete.gosper import gosper_normal


def _zeil(f, k, n, N, order):
    k_quotient = (f.subs(k, k + 1) / f).combsimp()
    n_quotients = [(f.subs(n, n + i) / f).combsimp() for i in range(order)]
    numers = [numer(q) for q in n_quotients]
    denoms = [denom(q) for q in n_quotients]

    D = prod(denoms)
    lead_quotient = (k_quotient * D / D.subs(k, k + 1)).simplify()
    r, s, p = gosper_normal(numer(lead_quotient), denom(lead_quotient), k)

    var_polys = symbols("p:{}".format(order), dummy=True)
    poly_part = p * sum(vp * a * D / b for vp, a, b in zip(var_polys, numers, denoms))
    poly_part = poly_part.simplify()

    for x_deg in range(6):
        cs = symbols("c:{}".format(x_deg+1), dummy=True)
        x = Poly(sum(c * k**i for i, c in enumerate(cs)), k)
        eqn = x.shift(1) * r - x * s.shift(-1) - poly_part
        soln = solve(eqn.all_coeffs(), set(var_polys) | set(cs))

        if soln and x.subs(soln) != 0:
            break

    free_vars = {v: 1 for v in set(var_polys) | set(cs) if v not in soln}

    operator = sum(vp * N**i for i, vp in enumerate(var_polys))
    pre_cert = x * s.shift(-1) / poly_part

    operator = Poly(operator.subs(soln).subs(free_vars), N)
    pre_cert = pre_cert.subs(soln).subs(free_vars)

    adjust = sum(quot * coeff for quot, coeff in zip(n_quotients[::-1], operator.all_coeffs()))
    cert = pre_cert * adjust
    cert = cert.factor()

    if operator != 0:
        clear = 1
        for c in operator.all_coeffs():
            clear = lcm(clear, denom(c), n)

        operator = (operator * clear).simplify()
        cert = (cert * clear).simplify()

    return operator.expr.simplify(), cert

def zeil(f, k, n, N, maxorder=6):
    r"""
    Execute Zeilberger's algorithm on ``f(n, k)``.

    Explanation
    ===========

    If ``f(n, k)`` is hypergeometric in both n and k, then it satisfies an
    identity of the form

    .. math::
        \sum_{i = 0}^d p_i(n) f(n + i, k) = \Delta_k G(n, k)

    where ``d`` is a positive integer, the ``p_i(n)`` are polynomials in ``n``,
    ``G(n, k)`` is a hypergeometric term in ``n`` and ``k``, and ``\Delta_k``
    is the difference operator in ``k``. This is called a "creative
    telescoping" recurrence. Zeilberger's algorithm computes this identity for
    any hypergeometric term.

    The above identity can be written as

    .. math::
        \left(\sum_{i = 0}^d p_i(n) N^i\right) f(n, k) = \Delta_k R(n, k) f(n, k)

    where ``N`` is the shift operator in ``n`` and ``R(n, k)`` is a rational
    function called the "certificate."

    This procedure returns ``(operator, certificate)`` where ``operator`` is
    the above polynomial in ``N`` and ``n`` ``certificate`` is the certifying
    rational function ``R(n, k)``.

    Examples
    ========

    >>> from sympy.concrete.zeil import zeil
    >>> from sympy.abc import k, n, N
    >>> from sympy import binomial

    >>> zeil(binomial(n, k), k, n, N)
    (N - 2, k/(k - n - 1))
    >>> zeil(binomial(n, k) / 2**n, k, n, N)
    (N - 1, k/(2*(k - n - 1)))
    >>> zeil(binomial(n, k) * (-1)**(n - k), k, n, N)
    (1, -k/n)
    >>> zeil(binomial(n, k)**2, k, n, N)
    (N*(n + 1) - 4*n - 2, k**2*(2*k - 3*n - 3)/(-k + n + 1)**2)

    """
    for order in range(1, maxorder+1):
        res = _zeil(f, k, n, N, order)
        if res and res[0] != 0:
            return res

    return None
