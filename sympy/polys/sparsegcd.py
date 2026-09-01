from __future__ import annotations

from sympy.polys.domains.integerring import ZZ
from sympy.polys.heuristicgcd import smp_heugcd
from typing import TYPE_CHECKING

from sympy.polys.monomials import monom, monomial_ldiv
from sympy.polys.monomials import monomial_gcd
from sympy.polys.orderings import MonomialOrder, lex
from sympy.polys.polyerrors import ExactQuotientFailed, HeuristicGCDFailed
from sympy.polys.sparseprs import smp_subresultants
from sympy.polys.sparsetools import  (smp_LC, smp_clear_denoms, smp_coeff_wrt, smp_deflate, smp_degrees, smp_div_list, smp_domain_convert, smp_inflate, smp_is_one,
    smp_mul, smp_mul_ground, smp_neg, smp_quo_ground, smp_quo_term, smp_reorg_poly, smp_swap_var)
from sympy.polys.zippel import _smp_zippel_gcd, _smp_zippel_gcd_mod

if TYPE_CHECKING:
    from sympy.external.gmpy import MPZ
    from sympy.polys.sparsetools import smp
    from sympy.polys.domains.domain import Domain, Er


def smp_cofactors_norm(
    gcd: smp[Er],
    cff: smp[Er],
    cfg: smp[Er],
    n: int,
    domain: Domain[Er],
    order: MonomialOrder = lex,
) -> tuple[smp[Er], smp[Er], smp[Er]]:
    # normalizes the gcd and cofactors.
    if not gcd:
        return gcd, cff, cfg

    if domain.is_Field:
        u = domain.revert(smp_LC(gcd, n, domain, order))

    else:
        u = domain.canonical_unit(smp_LC(gcd, n, domain, order))

    gcd = smp_mul_ground(gcd, u, n, domain)
    cff = smp_quo_ground(cff, u, n, domain)
    cfg = smp_quo_ground(cfg, u, n, domain)

    return gcd, cff, cfg


def smp_primitive_wrt_ZZ(g: smp[MPZ], i: int, n: int
) -> tuple[smp[MPZ], smp[MPZ]]:
    """
    Returns the content and primitive part of a polynomial with respect to a
    specified variable.
    """
    g_r = smp_reorg_poly(g, [i], n, ZZ)
    cont = smp_gcd_list_ZZ(list(g_r.values()), n-1)
    map = [j for j in range(n) if j != i]
    cont_ = smp_elevate(cont, map, n, ZZ)
    [prim], _ = smp_div_list(g, [cont_], n, ZZ)
    return cont, prim


def smp_zippel_gcd(
    A: smp[MPZ],
    B: smp[MPZ],
    n: int,
) -> tuple[smp[MPZ], smp[MPZ], smp[MPZ]]:

    #Check whether the GCD is monic in any variable.

    #If no suitable variable is found, extract the content with respect
    # to the first variable and compute the GCD of the primitive parts.

    for i in range(n):
        m_A = max(mon[i] for mon in A)
        m_B = max(mon[i] for mon in B)

        lc_A = smp_coeff_wrt(A, i, m_A, n, ZZ)
        lc_B = smp_coeff_wrt(B, i, m_B, n, ZZ)

        if len(lc_A) == 1 or len(lc_B) == 1:
            if i == 0:
                return _smp_zippel_gcd(A, B, n)

            A_ = smp_swap_var(A, i, n, ZZ)
            B_ = smp_swap_var(B, i, n, ZZ)

            gcd, cfa, cfb = _smp_zippel_gcd(A_, B_, n)

            return (
                smp_swap_var(gcd, i, n, ZZ),
                smp_swap_var(cfa, i, n, ZZ),
                smp_swap_var(cfb, i, n, ZZ),
            )

    A_cont, A_prim = smp_primitive_wrt_ZZ(A, 0, n)
    B_cont, B_prim = smp_primitive_wrt_ZZ(B, 0, n)

    cont_gcd, cont_cfa, cont_cfb = smp_cofactors_ZZ(
        A_cont,
        B_cont,
        n - 1,
    )

    gcd, cfa, cfb = _smp_zippel_gcd(A_prim, B_prim, n)

    indices = list(range(1, n))

    cont_gcd = smp_elevate(cont_gcd, indices, n, ZZ)
    cont_cfa = smp_elevate(cont_cfa, indices, n, ZZ)
    cont_cfb = smp_elevate(cont_cfb, indices, n, ZZ)

    return (
        smp_mul(cont_gcd, gcd, ZZ, n),
        smp_mul(cont_cfa, cfa, ZZ, n),
        smp_mul(cont_cfb, cfb, ZZ, n),
    )


def smp_elevate(f, map: list[int], n: int, domain: Domain[Er]) -> smp[Er]:
    # adds variables to a polynomial
    h = {}
    for mon, coeff in f.items():
        m = [0]*n
        for i, exp in zip(map, mon):
            m[i] = exp
        h[tuple(m)] = coeff
    return h


def smp_monomial_extract_ZZ(A: smp[MPZ], B: smp[MPZ], n: int
    ) -> tuple[smp[MPZ], smp[MPZ], smp[MPZ]]:
    # Extracts any common monomial from the polynomials.
    zm = (0,) * n
    monoms: list[monom] = []
    monoms.extend(A.keys())
    monoms.extend(B.keys())

    # Check for the presence of zero_monom in any polynomial
    if zm in monoms:
        return A, B, {zm: ZZ.one}

    mgcd = monoms[0]
    for mon in monoms[1:]:
        mgcd = monomial_gcd(mgcd, mon)

        # If gcd becomes the zero monomial, return early
        if mgcd == zm:
            return A, B, {zm: ZZ.one}

    gcd_term = (mgcd, ZZ.one)
    gcd: smp[MPZ] = {mgcd: ZZ.one}
    A = smp_quo_term(A, gcd_term, n, ZZ)
    B = smp_quo_term(B, gcd_term, n, ZZ)
    return A, B, gcd


def smp_cofactors(A: smp[Er], B: smp[Er], n: int, domain: Domain[Er], order: MonomialOrder = lex
) -> tuple[smp[Er], smp[Er], smp[Er]]:
    # Entry point of the gcd API:

    result = _smp_cofactors_trivial(A, B, n, domain, order)

    if result is not None:
        return result

    if result is None:
        J, (A, B) = smp_deflate((A, B), n, domain)

        #changes the domain of the inputs in one managable by the sympy algorithms
        if not domain.is_Exact:
            domain_ex = domain.get_exact()
            A_ex = smp_domain_convert(A, domain, domain_ex, n)
            B_ex = smp_domain_convert(B, domain, domain_ex, n)

            h_ex, cff_ex, cfg_ex = smp_cofactors(
                A_ex, B_ex, n, domain_ex, order
            )

            h = smp_domain_convert(h_ex, domain_ex, domain, n)
            cff = smp_domain_convert(cff_ex, domain_ex, domain, n)
            cfg = smp_domain_convert(cfg_ex, domain_ex, domain, n)
        else:
            h, cff, cfg = _smp_cofactors(A, B, n, domain)

        result = (
            smp_inflate(h, J, n, domain),
            smp_inflate(cff, J, n, domain),
            smp_inflate(cfg, J, n, domain),
        )

    return smp_cofactors_norm(*result, n, domain, order)


def _smp_cofactors(A: smp[Er], B: smp[Er], n: int, domain: Domain[Er]
) -> tuple[smp[Er], smp[Er], smp[Er]]:
    if domain.is_QQ:
        return smp_cofactors_QQ(A, B, n, domain)

    elif domain.is_ZZ:
        return smp_cofactors_ZZ(A, B, n) # type: ignore

    elif domain.is_FiniteField:
        return smp_cofactors_FF(A, B, n, domain)

    else:
        return smp_cofactors_prs(A, B, n, domain)


def _smp_cofactors_trivial(
    A: smp[Er],
    B: smp[Er],
    n: int,
    domain: Domain[Er],
    order: MonomialOrder = lex,
) -> tuple[smp[Er], smp[Er], smp[Er]] | None:
    # Returns the gcd and cofactors of the polynmials in
    # trivial cases.
    if not A and not B:
        return {}, {}, {}

    if not A:
        return smp_gcd_zero(A, B, n, domain, order)

    if not B:
        h, cfg, cff = smp_gcd_zero(B, A, n, domain, order)
        return h, cff, cfg

    if len(A) == 1:
        return smp_gcd_monom(A, B, n, domain)

    if len(B) == 1:
        h, cfg, cff = smp_gcd_monom(B, A, n, domain)
        return h, cff, cfg

    return None


def smp_cofactors_prs(A: smp[Er], B: smp[Er], n: int, domain: Domain[Er]
) -> tuple[smp[Er], smp[Er], smp[Er]]:
    # Computes the gcd and cofactors using the PRS.
    gcd = smp_gcd_prs(A, B, n, domain)

    [cff], r = smp_div_list(A, [gcd], n, domain)
    if r:
        raise ExactQuotientFailed(A, gcd)

    [cfg], r = smp_div_list(B, [gcd], n, domain)
    if r:
        raise ExactQuotientFailed(B, gcd)

    return gcd, cff, cfg


def smp_zippel_gcd_mod(A: smp[Er], B: smp[Er], n: int, domain: Domain[Er],
) -> smp[Er] | None:
    """
    Check whether the GCD is monic in any variable.

    If no suitable variable is found, extract the content with respect
    to the first variable and compute the GCD of the primitive parts.
    Then it converts the inputs coefficients to Z to perform the gcd computation
    with _smp_zippel_gcd_mod, and if the computation was successful it converts the
    coefficients back to the finite field. If it wasn't it returns None.
    """

    p = ZZ(domain.characteristic())

    for i in range(n):
        m_A = max(el[i] for el in A.keys())
        m_B = max(el[i] for el in B.keys())
        lc_A = smp_coeff_wrt(A, i, m_A, n, domain)
        lc_B = smp_coeff_wrt(B, i, m_B, n, domain)
        if len(lc_A) == 1 or len(lc_B) == 1:
            A_z = smp_domain_convert(A, domain, ZZ, n)
            B_z = smp_domain_convert(B, domain, ZZ, n)
            if i == 0:
                gcd = _smp_zippel_gcd_mod(A_z, B_z, p, n)
                if gcd is None:
                    return None
                return smp_domain_convert(gcd, ZZ, domain, n)
            A_z_ = smp_swap_var(A_z, i, n, ZZ)
            B_z_ = smp_swap_var(B_z, i, n, ZZ)
            gcd = _smp_zippel_gcd_mod(A_z_, B_z_, p, n)
            if gcd is None:
                return None
            gcd = smp_swap_var(gcd, i, n, ZZ)
            return smp_domain_convert(gcd, ZZ, domain, n)

    Ac, Ap = smp_primitive_wrt(A, 0, n, domain)
    Bc, Bp = smp_primitive_wrt(B, 0, n, domain)
    Ap_z = smp_domain_convert(Ap, domain, ZZ, n)
    Bp_z = smp_domain_convert(Bp, domain, ZZ, n)
    cont_gcd, _, _ = smp_cofactors_FF(Ac, Bc, n-1, domain)
    gcdp = _smp_zippel_gcd_mod(Ap_z, Bp_z, p, n)
    if gcdp is None:
        return None
    cont_gcd = smp_elevate(cont_gcd, list(range(1, n)), n, domain)
    gcdp_dom: smp[Er] = smp_domain_convert(gcdp, ZZ, domain, n)
    return smp_mul(cont_gcd, gcdp_dom, domain, n)


def smp_primitive_wrt_FF(g: smp[Er], i: int, n: int, domain
) -> tuple[smp[Er], smp[Er]]:
    """
    Returns the content and primitive part of a polynomial with respect to a
    specified variable.
    """
    g_r = smp_reorg_poly(g, [i], n, domain)
    cont = smp_gcd_list_FF(list(g_r.values()), n-1, domain)
    map = [j for j in range(n) if j != i]
    cont_ = smp_elevate(cont, map, n, domain)
    [prim], _ = smp_div_list(g, [cont_], n, domain)
    return cont, prim


def smp_primitive_wrt(g: smp[Er], i: int, n: int, domain
) -> tuple[smp[Er], smp[Er]]:
    """
    Returns the content and primitive part of a polynomial with respect to a
    specified variable.
    """
    g_r = smp_reorg_poly(g, [i], n, domain)
    cont = smp_gcd_list(list(g_r.values()), n-1, domain)
    map = [j for j in range(n) if j != i]
    cont_ = smp_elevate(cont, map, n, domain)
    [prim], _ = smp_div_list(g, [cont_], n, domain)
    return cont, prim


def smp_gcd_list(polys: list[smp[Er]], n: int, domain) -> smp[Er]:
    # Returns the gcd of a list of polynomials.
    if not polys:
        return {}
    gcd = polys[0]
    for pol in polys[1:]:
        result = _smp_cofactors_trivial(gcd, pol, n, domain)

        if result is None:
            gcd, _, _ = smp_cofactors(gcd, pol, n, domain)
        else:
            gcd = result[0]
        if smp_is_one(gcd, n, domain):
            break
    return gcd


def smp_cofactors_FF(A: smp[Er], B: smp[Er], n: int, domain: Domain[Er],
) -> tuple[smp[Er], smp[Er], smp[Er]]:
    # Returns the gcd and cofactors when domain is a finite field.
    p = domain.characteristic()
    # 10**5 was chosen because in a low characteristic field Zippel's algorithm could
    # run out of evaluation points, and the algorithm could run into an infinite loop.
    # The choice is under the assumption that the number of var and the degree of the polynomials
    # will be considerably smaller than 10**5. The choice could be improved, allowing smaller fields.
    if p > 10**5:
        gcd = smp_zippel_gcd_mod(A, B, n, domain)
        # smp_zippel_gcd_mod could be rewritten so that when it performs test
        # division it also returns the cofactors
        if gcd:
            [cff], rff = smp_div_list(A, [gcd], n, domain)
            [cfg], rfg = smp_div_list(B, [gcd], n, domain)

            if not rff and not rfg:
                return gcd, cff, cfg

    return smp_cofactors_prs(A, B, n, domain)


def smp_gcd_prs(A: smp[Er], B: smp[Er], n: int, domain: Domain[Er]
) -> smp[Er]:
    """
    Compute a greatest common divisor using a subresultant polynomial
    remainder sequence (PRS).

    The result is not normalized and may differ from a canonical GCD by a
    unit of the coefficient domain. Normalization is performed by
    ``smp_cofactors``.

    Examples
    ========

    >>> from sympy.polys.sparsegcd import smp_gcd_prs
    >>> from sympy import ZZ, ring
    >>> R, x, y = ring("x, y", ZZ)
    >>> c = x**2 + x*y + y**2 - 1
    >>> f = (x + 1)*c
    >>> g = (x - 1)*c
    >>> h = R.from_dict(smp_gcd_prs(f, g, 2, ZZ))
    >>> h == c or h == -c
    True
    """

    if n == 0:
        a = next(iter(A.values()), domain.zero)
        b = next(iter(B.values()), domain.zero)
        gc = domain.gcd(a, b)
        return {(): gc} if gc else {}

    map = list(range(1,n))
    c1, pp1 = smp_primitive_wrt(A, 0, n, domain)
    c2, pp2 = smp_primitive_wrt(B, 0, n, domain)

    cont = smp_gcd_prs(c1, c2, n-1, domain)

    h = smp_subresultants(pp1, pp2, 0, n, domain)[-1]
    _, h = smp_primitive_wrt(h, 0, n, domain)
    map = list(range(1,n))
    h = smp_mul(h, smp_elevate(cont, map, n, domain), domain, n)

    return h


def smp_cofactors_QQ(
    A: smp[Er],
    B: smp[Er],
    n: int,
    domain: Domain[Er],
) -> tuple[smp[Er], smp[Er], smp[Er]]:
    # Returns gcd and cofactors when domain is Q.
    cf, A_q = smp_clear_denoms(A, n, domain)
    cg, B_q = smp_clear_denoms(B, n, domain)

    A_z = smp_domain_convert(A_q, domain, ZZ, n)
    B_z = smp_domain_convert(B_q, domain, ZZ, n)

    h_z, cff_z, cfg_z = smp_cofactors_ZZ(A_z, B_z, n)

    h = smp_domain_convert(h_z, ZZ, domain, n)
    cff = smp_domain_convert(cff_z, ZZ, domain, n)
    cfg = smp_domain_convert(cfg_z, ZZ, domain, n)

    cff = smp_quo_ground(cff, cf, n, domain)
    cfg = smp_quo_ground(cfg, cg, n, domain)

    return h, cff, cfg


def smp_cofactors_ZZ(A: smp[MPZ], B: smp[MPZ], n: int
) -> tuple[smp[MPZ], smp[MPZ], smp[MPZ]]:
    # Dispatches the gcd to either zippel or heugcd and returns gcd and cofactors.
    # The dispatching choice is based on a formula fitted on benchmarking data.

    zm: monom = (0,) * n
    one: smp[MPZ] = {zm: ZZ.one}

    A, B, m = smp_monomial_extract_ZZ(A, B, n)

    terms = min(len(A), len(B))

    deg_A = max(map(sum, A))
    deg_B = max(map(sum, B))
    deg = min(deg_A, deg_B)

    if n <= 3 or deg < 4:
        try:
            h, cfa, cfb = smp_heugcd(A, B, n)
        except HeuristicGCDFailed:
            h, cfa, cfb = smp_zippel_gcd(A, B, n)
        if m != one:
            h = smp_mul(h, m, ZZ, n)

        return h, cfa, cfb

    max_deg_var = min(max(smp_degrees(A, n, ZZ)), max(smp_degrees(B, n, ZZ)))

    # Total degree can overstate the cost when it is spread over many
    # variables, so use a lower effective degree for widely spread inputs.
    deg -= max(0, deg // max_deg_var - 1)
    score = (n - 3) * (deg - 3) + deg // 2

    if score >= 150 or terms * max(score - 25, 1) > 500:
        h, cfa, cfb = smp_zippel_gcd(A, B, n)
    else:
        try:
            h, cfa, cfb = smp_heugcd(A, B, n)
        except HeuristicGCDFailed:
            h, cfa, cfb = smp_zippel_gcd(A, B, n)

    if m != one:
        h = smp_mul(h, m, ZZ, n)

    return h, cfa, cfb


def smp_gcd_list_ZZ(polys: list[smp[MPZ]], n: int) -> smp[MPZ]:
    # returns the gcd of a list of polynomials.
    if not polys:
        return {}
    gcd = polys[0]
    for pol in polys[1:]:
        result = _smp_cofactors_trivial(gcd, pol, n, ZZ)

        if result is None:
            gcd, _, _ = smp_cofactors_ZZ(gcd, pol, n)
        else:
            gcd = result[0]
        if smp_is_one(gcd, n, ZZ):
            break
    return gcd


def smp_gcd_list_FF(polys: list[smp[Er]], n: int, domain) -> smp[Er]:
    # returns the gcd of a list of polynomials.
    if not polys:
        return {}
    gcd = polys[0]
    for pol in polys[1:]:
        result = _smp_cofactors_trivial(gcd, pol, n, domain)

        if result is None:
            gcd, _, _ = smp_cofactors_FF(gcd, pol, n, domain)
        else:
            gcd = result[0]
        if smp_is_one(gcd, n, domain):
            break
    return gcd


def smp_gcd_zero(
    A: smp[Er],
    B: smp[Er],
    n: int,
    domain: Domain[Er],
    order: MonomialOrder = lex,
) -> tuple[smp[Er], smp[Er], smp[Er]]:
    # Return gcd and cofactors when A is zero and B is nonzero.
    zm: monom = (0,) * n
    zero: smp[Er] = {}
    one: smp[Er] = {zm: domain.one}

    if domain.is_nonnegative(smp_LC(B, n, domain, order)):
        return B, zero, one

    return smp_neg(B, n, domain), zero, smp_neg(one, n, domain)


def smp_gcd_monom(
    A: smp[Er],
    B: smp[Er],
    n: int,
    domain: Domain[Er],
) -> tuple[smp[Er], smp[Er], smp[Er]]:
    # Return gcd and cofactors when A consists of exactly one term.
    mf, cf = next(iter(A.items()))

    mgcd = mf
    cgcd = cf

    for mg, cg in B.items():
        mgcd = monomial_gcd(mgcd, mg)
        cgcd = domain.gcd(cgcd, cg)

    h = {mgcd: cgcd}

    cff = {monomial_ldiv(mf, mgcd): domain.quo(cf, cgcd)}

    cfg = {monomial_ldiv(mg, mgcd): domain.quo(cg, cgcd) for mg, cg in B.items()}

    return h, cff, cfg
