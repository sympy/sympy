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
from sympy.polys.sparse_tools import (
    dict_is_ground,
    dict_LC,
    dict_quo_ground,
    dict_trunc_ground,
)

if TYPE_CHECKING:
    from sympy.external.gmpy import MPZ
    from sympy.polys.densebasic import dup
    from sympy.polys.domains.domain import Domain

Monom = tuple[int, ...]


def dict_gf_gcd(fp: dict[Monom, MPZ], gp: dict[Monom, MPZ], p: MPZ,
    dom: Domain[MPZ]) -> dict[Monom, MPZ]:
    r"""
    Compute the GCD of two raw univariate polynomial dictionaries in
    `\mathbb{Z}_p[x]`.
    """
    f_list = dup_from_dict(fp, dom)
    g_list = dup_from_dict(gp, dom)

    gcd_list = gf_gcd(f_list, g_list, p, dom)
    gcd = dup_to_dict(gcd_list, dom)
    return dict_trunc_ground(gcd, p, 1, dom)


def dict_trivial_gcd(f: dict[Monom, MPZ], g: dict[Monom, MPZ],
    n: int, dom: Domain[MPZ]
    ) -> tuple[dict[Monom, MPZ], dict[Monom, MPZ],
               dict[Monom, MPZ]] | None:
    """
    Compute the GCD and cofactors for trivial cases.
    """
    zero_monom = (0,)*n
    zero: dict[Monom, MPZ] = {}
    one: dict[Monom, MPZ] = {zero_monom: dom.one}

    def neg(poly: dict[Monom, MPZ]) -> dict[Monom, MPZ]:
        return {mon: -coeff for mon, coeff in poly.items() if coeff}

    if not (f or g):
        return zero, zero, zero
    elif not f:
        if dict_LC(g, n, dom) < dom.zero:
            return neg(g), zero, neg(one)
        else:
            return g.copy(), zero, one
    elif not g:
        if dict_LC(f, n, dom) < dom.zero:
            return neg(f), neg(one), zero
        else:
            return f.copy(), one, zero
    elif dict_is_ground(f, n, dom) and dict_is_ground(g, n, dom):
        c = dom.gcd(dict_LC(f, n, dom), dict_LC(g, n, dom))
        h: dict[Monom, MPZ] = {zero_monom: c}
        return h, dict_quo_ground(f, c, n, dom), dict_quo_ground(g, c, n, dom)
    return None


def dict_primitive_wrt_last(f: dict[Monom, MPZ], n: int, dom: Domain[MPZ],
    p: MPZ) -> tuple[dup[MPZ], dict[Monom, MPZ]]:
    """
    Computes the content and primitive part of the poly
    with respect to the first n-1 variables: therefore content is in the last
    variable and prim. part is in the other variables.
    """
    coeffs: dict[Monom, dict[int, MPZ]] = {}
    for monom, coeff in f.items():
        if monom[:-1] not in coeffs:
            coeffs[monom[:-1]] = {}
        coeffs[monom[:-1]][monom[-1]] = coeff

    cont: dup[MPZ] = []
    dense_coeffs: dict[Monom, dup[MPZ]] = {}

    for mon, coeff_dict in coeffs.items():
        dense_coeff = gf_from_dict(coeff_dict, p, dom)
        dense_coeffs[mon] = dense_coeff
        cont = gf_gcd(cont, dense_coeff, p, dom)

    prim: dict[Monom, MPZ] = {}
    for mon, dense_coeff in dense_coeffs.items():
        quotient = gf_quo(dense_coeff, cont, p, dom)
        deg = len(quotient)
        for i, el in enumerate(quotient):
            if el != 0:
                prim[mon + (deg-1-i,)] = el

    return cont, dict_trunc_ground(prim, p, n, dom)


def dict_deg_wrt_last(f: dict[Monom, MPZ], n: int) -> Monom:
    """
    Computes the highest degree Monom with respect to the first n-1 variables.
    """
    degf: Monom = (0,) * (n-1)
    for monom in f:
        if monom[:-1] > degf:
            degf = monom[:-1]
    return degf


def dict_LC_wrt_last(f: dict[Monom, MPZ], n: int,
    dom: Domain[MPZ]) -> dict[Monom, MPZ]:
    """
    Computes the LC with respect to the first n-1 variables.
    """
    degf = dict_deg_wrt_last(f, n)
    lcf: dict[Monom, MPZ] = {}

    for monom, coeff in f.items():
        if monom[:-1] == degf:
            lcf[(monom[-1],)] = coeff

    return lcf


def dict_chinese_remainder_reconstruction_multivariate(
    hp: dict[Monom, MPZ], hq: dict[Monom, MPZ], p: MPZ, q: MPZ, dom: Domain[MPZ], n: int,
) -> dict[Monom, MPZ]:
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

    hpq: dict[Monom, MPZ] = {}

    def crt_scalar(cp, cq, p, q):
        return dom(crt([p, q], [cp, cq], symmetric=True)[0])

    for monom in monoms:
        coeff = crt_scalar(hp[monom], hq[monom], p, q)
        if coeff:
            hpq[monom] = coeff
    for monom in hpmonoms:
        coeff = crt_scalar(hp[monom], zero, p, q)
        if coeff:
            hpq[monom] = coeff
    for monom in hqmonoms:
        coeff = crt_scalar(zero, hq[monom], p, q)
        if coeff:
            hpq[monom] = coeff

    return hpq
