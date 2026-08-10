from __future__ import annotations

from typing import TYPE_CHECKING

from sympy.polys.polyerrors import ExactQuotientFailed
from sympy.polys.sparsetools import (
    _smp_imul_ground,
    smp_add,
    smp_coeff_wrt,
    smp_degree,
    smp_div_list,
    smp_mul,
    smp_neg,
    smp_pow_generic,
    smp_sub,
)

if TYPE_CHECKING:
    from sympy.polys.domains.domain import Domain, Er
    from sympy.polys.sparsetools import smp


def smp_prem(f: smp[Er], g: smp[Er], i: int, n: int, dom: Domain[Er]) -> smp[Er]:

    df = smp_degree(f, i, n, dom)
    dg = smp_degree(g, i, n, dom)
    if dg < 0:
        raise ZeroDivisionError("polynomial division")

    r, dr = f, df

    if df < dg:
        return r

    N = df - dg + 1

    lc_g = smp_coeff_wrt(g, i, dg, n, dom)

    while True:
        lc_r = smp_coeff_wrt(r, i, dr, n, dom)
        j, N = dr - dg, N - 1

        R = smp_mul(r, lc_g, dom, n)

        xp = {(0,) * i + (j,) + (0,) * (n - i - 1): dom.one}

        G = smp_mul(smp_mul(g, lc_r, dom, n), xp, dom, n)

        r = smp_sub(R, G, dom, n)

        dr = smp_degree(r, i, n, dom)
        # dr = r._degree_int(x)

        if dr < dg:
            break

    # c = lc_g**N
    c = smp_pow_generic(lc_g, N, dom, n)

    return smp_mul(c, r, dom, n)


def smp_pdiv(
    f: smp[Er], g: smp[Er], i: int, n: int, dom: Domain[Er]
) -> tuple[smp[Er], smp[Er]]:
    df = smp_degree(f, i, n, dom)
    dg = smp_degree(g, i, n, dom)

    if dg < 0:
        raise ZeroDivisionError("polynomial division")

    q: smp[Er] = {}
    r, dr = f, df

    if df < dg:
        return q, r

    N = df - dg + 1
    lc_g = smp_coeff_wrt(g, i, dg, n, dom)

    while True:
        lc_r = smp_coeff_wrt(r, i, dr, n, dom)
        j, N = dr - dg, N - 1

        Q = smp_mul(q, lc_g, dom, n)

        xp = {(0,) * i + (j,) + (0,) * (n - i - 1): dom.one}

        q = smp_add(Q, smp_mul(lc_r, xp, dom, n), dom, n)

        R = smp_mul(r, lc_g, dom, n)

        G = smp_mul(smp_mul(g, lc_r, dom, n), xp, dom, n)

        r = smp_sub(R, G, dom, n)

        dr = smp_degree(r, i, n, dom)

        if dr < dg:
            break

    c = smp_pow_generic(lc_g, N, dom, n)

    q = smp_mul(q, c, dom, n)
    r = smp_mul(r, c, dom, n)

    return q, r


def smp_pquo(f: smp[Er], g: smp[Er], i: int, n: int, dom: Domain[Er]) -> smp[Er]:
    return smp_pdiv(f, g, i, n, dom)[0]


def smp_pexquo(f: smp[Er], g: smp[Er], i: int, n: int, dom: Domain[Er]) -> smp[Er]:
    q, r = smp_pdiv(f, g, i, n, dom)

    if r == {}:
        return q
    else:
        raise ExactQuotientFailed(f, g)


def smp_subresultants(f: smp[Er], g: smp[Er], i: int, n: int, dom: Domain[Er]) -> list[smp[Er]]:

        l = smp_degree(f, i, n, dom)
        m = smp_degree(g, i, n, dom)

        if l < m:
            f, g = g, f
            l, m = m, l

        if not f:
            return [{}, {}]

        if not g:
            zm = (0,) * n
            return [f, {zm: dom.one}]

        R = [f, g]

        d = l - m
        b = dom((-1) ** (d + 1))

        # Compute the pseudo-remainder for f and g
        h = smp_prem(f, g, i, n, dom)

        _smp_imul_ground(h, b, n, dom)

        # Compute the coefficient of g with respect to x**m
        lc = smp_coeff_wrt(g, i, m, n, dom)

        c = smp_pow_generic(lc, d, dom, n)

        c = smp_neg(c, n, dom)
        while h:

            k = smp_degree(h, i, n, dom)

            R.append(h)
            f, g, m, d = g, h, k, m - k

            b = smp_mul(smp_neg(lc, n, dom), smp_pow_generic(c, d, dom, n), dom, n)

            h = smp_prem(f, g, i, n, dom)

            [h], _ = smp_div_list(h, [b], n, dom)
            lc = smp_coeff_wrt(g, i, k, n, dom)

            if d > 1:

                p = smp_pow_generic(smp_neg(lc, n, dom), d, dom, n)

                q = smp_pow_generic(c, (d - 1), dom, n)

                [c], _ = smp_div_list(p, [q], n, dom)
            else:

                c = smp_neg(lc, n, dom)

        return R
