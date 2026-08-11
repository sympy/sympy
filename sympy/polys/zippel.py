"""Low-level Zippel algorithm functions on raw dict representations.

These functions deliberately do not depend on PolyRing or PolyElement.
They operate on dicts mapping exponent tuples to coefficients and use
only the provided domain and number of variables.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from sympy.polys.densebasic import dup_from_dict, dup_to_dict
from sympy.polys.galoistools import gf_gcd, gf_quo, gf_from_dict
from sympy.ntheory.modular import crt
from sympy.polys.sparsetools import (
    smp_is_ground,
    smp_LC,
    smp_quo_ground,
    smp_trunc_ground,
)

if TYPE_CHECKING:
    from sympy.external.gmpy import MPZ
    from sympy.polys.densebasic import dup
    from sympy.polys.domains.domain import Domain
    from sympy.polys.monomials import monom
    from sympy.polys.sparsetools import smp


def smp_gf_gcd(fp: smp[MPZ], gp: smp[MPZ], p: MPZ, dom: Domain[MPZ]) -> smp[MPZ]:
    r"""
    Compute the GCD of two raw univariate polynomial dictionaries in
    `\mathbb{Z}_p[x]`.
    """
    f_list = dup_from_dict(fp, dom)
    g_list = dup_from_dict(gp, dom)

    gcd_list = gf_gcd(f_list, g_list, p, dom)
    gcd = dup_to_dict(gcd_list, dom)
    return smp_trunc_ground(gcd, p, 1, dom)


def smp_trivial_gcd(
    f: smp[MPZ], g: smp[MPZ], n: int, dom: Domain[MPZ]
) -> tuple[smp[MPZ], smp[MPZ], smp[MPZ]] | None:
    """
    Compute the GCD and cofactors for trivial cases.
    """
    zm = (0,) * n
    zero: smp[MPZ] = {}
    one: smp[MPZ] = {zm: dom.one}

    def neg(poly: smp[MPZ]) -> smp[MPZ]:
        return {mon: -coeff for mon, coeff in poly.items() if coeff}

    if not (f or g):
        return zero, zero, zero
    elif not f:
        if smp_LC(g, n, dom) < dom.zero:
            return neg(g), zero, neg(one)
        else:
            return g.copy(), zero, one
    elif not g:
        if smp_LC(f, n, dom) < dom.zero:
            return neg(f), neg(one), zero
        else:
            return f.copy(), one, zero
    elif smp_is_ground(f, n, dom) and smp_is_ground(g, n, dom):
        c = dom.gcd(smp_LC(f, n, dom), smp_LC(g, n, dom))
        h: smp[MPZ] = {zm: c}
        return h, smp_quo_ground(f, c, n, dom), smp_quo_ground(g, c, n, dom)
    return None


def smp_primitive_wrt_last(
    f: smp[MPZ], n: int, dom: Domain[MPZ], p: MPZ
) -> tuple[dup[MPZ], smp[MPZ]]:
    """
    Computes the content and primitive part of the poly
    with respect to the first n-1 variables: therefore content is in the last
    variable and prim. part is in the other variables.
    """
    coeffs: dict[monom, dict[int, MPZ]] = {}
    for mon, coeff in f.items():
        if mon[:-1] not in coeffs:
            coeffs[mon[:-1]] = {}
        coeffs[mon[:-1]][mon[-1]] = coeff

    cont: dup[MPZ] = []
    dense_coeffs: dict[monom, dup[MPZ]] = {}

    for mon, coeff_dict in coeffs.items():
        dense_coeff = gf_from_dict(coeff_dict, p, dom)
        dense_coeffs[mon] = dense_coeff
        cont = gf_gcd(cont, dense_coeff, p, dom)

    prim: smp[MPZ] = {}
    for mon, dense_coeff in dense_coeffs.items():
        quotient = gf_quo(dense_coeff, cont, p, dom)
        deg = len(quotient)
        for i, el in enumerate(quotient):
            if el != 0:
                prim[mon + (deg - 1 - i,)] = el

    return cont, smp_trunc_ground(prim, p, n, dom)


def smp_deg_wrt_last(f: smp[MPZ], n: int) -> monom:
    """
    Computes the highest degree monomial with respect to the first n-1 variables.
    """
    degf: monom = (0,) * (n - 1)
    for mon in f:
        if mon[:-1] > degf:
            degf = mon[:-1]
    return degf


def smp_LC_wrt_last(f: smp[MPZ], n: int, dom: Domain[MPZ]) -> smp[MPZ]:
    """
    Computes the LC with respect to the first n-1 variables.
    """
    degf = smp_deg_wrt_last(f, n)
    lcf: smp[MPZ] = {}

    for mon, coeff in f.items():
        if mon[:-1] == degf:
            lcf[(mon[-1],)] = coeff

    return lcf


def smp_chinese_remainder_reconstruction_multivariate(
    hp: smp[MPZ],
    hq: smp[MPZ],
    p: MPZ,
    q: MPZ,
    dom: Domain[MPZ],
    n: int,
) -> smp[MPZ]:
    r"""
    Construct a polynomial `h_{pq}` in
    `\mathbb{Z}_{p q}[x_0, \ldots, x_{k-1}]` such that

    .. math ::

        h_{pq} = h_p \; \mathrm{mod} \, p

        h_{pq} = h_q \; \mathrm{mod} \, q

    for relatively prime integers `p` and `q` and polynomials
    `h_p` and `h_q` in `\mathbb{Z}_p[x_0, \ldots, x_{k-1}]` and
    `\mathbb{Z}_q[x_0, \ldots, x_{k-1}]` respectively.

    The coefficients of the polynomial `h_{pq}` are computed with the
    Chinese Remainder Theorem. The symmetric representation in
    `\mathbb{Z}_p[x_0, \ldots, x_{k-1}]`,
    `\mathbb{Z}_q[x_0, \ldots, x_{k-1}]` and
    `\mathbb{Z}_{p q}[x_0, \ldots, x_{k-1}]` is used.

    Parameters
    ==========

    hp : dict
        raw multivariate integer polynomial with coefficients in `\mathbb{Z}_p`
    hq : dict
        raw multivariate integer polynomial with coefficients in `\mathbb{Z}_q`
    p : Integer
        modulus of `h_p`, relatively prime to `q`
    q : Integer
        modulus of `h_q`, relatively prime to `p`

    """
    hpmonoms = set(hp.keys())
    hqmonoms = set(hq.keys())
    monoms = hpmonoms.intersection(hqmonoms)
    hpmonoms.difference_update(monoms)
    hqmonoms.difference_update(monoms)

    zero = dom.zero

    hpq: smp[MPZ] = {}

    def crt_scalar(cp, cq, p, q):
        return dom(crt([p, q], [cp, cq], symmetric=True)[0])

    for mon in monoms:
        coeff = crt_scalar(hp[mon], hq[mon], p, q)
        if coeff:
            hpq[mon] = coeff
    for mon in hpmonoms:
        coeff = crt_scalar(hp[mon], zero, p, q)
        if coeff:
            hpq[mon] = coeff
    for mon in hqmonoms:
        coeff = crt_scalar(zero, hq[mon], p, q)
        if coeff:
            hpq[mon] = coeff

    return hpq
