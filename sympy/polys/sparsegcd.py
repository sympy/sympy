from __future__ import annotations

import random

from sympy.polys.domains.integerring import ZZ
from sympy.polys.domains.rationalfield import QQ
from sympy.polys.heuristicgcd import smp_heugcd
from sympy.polys.monomials import monomial_ngcd
from typing import TYPE_CHECKING, Callable, Iterable

from sympy.polys.orderings import MonomialOrder, lex
from sympy.polys.polyerrors import HeuristicGCDFailed
from sympy.polys.sparseprs import smp_subresultants
from sympy.polys.sparsetools import  (_smp_imul_ground, smp_LC, smp_active_var, smp_add, smp_clear_denoms, smp_coeff_wrt, smp_div_list,
    smp_domain_convert, smp_is_one, smp_monomial_extract, smp_mul, smp_mul_ground, smp_rem_list, smp_reorg_poly, smp_swap_var, smp_unique_polys)
from sympy.polys.zippel import _smp_zippel_gcd, _smp_zippel_gcd_mod

if TYPE_CHECKING:
    from sympy.external.gmpy import MPZ
    from sympy.polys.sparsetools import smp
    from sympy.polys.domains.domain import Domain, Er


def gcd_terms(polynomials: Iterable[smp[Er]], domain: Domain[Er]) -> smp[Er]:
    """
    Returns the greatest common divisor (GCD) of all terms in a list of
    polynomials p with respect to a given ring and domain.

    Examples
    ========

    >>> from sympy.polys.rings import gcd_terms
    >>> from sympy import ZZ, ring
    >>> R, x, y = ring("x, y", ZZ)

    >>> p1 = 2*x**3*y**2 - 2*x**2*y**3
    >>> p2 = 4*x**2*y**3 - 4*x**3*y**2
    >>> polynomials = [p1, p2]
    >>> ring = p1.ring
    >>> domain = ring.domain
    >>> R(gcd_terms(polynomials, domain))
    2*x**2*y**2
    >>> p1.gcd(p2) # Shows the difference between the gcd_terms and gcd
    2*x**3*y**2 - 2*x**2*y**3

    """
    monomials = set()
    coeffs = set()

    for pi in polynomials:
        for monomial, coeff in pi.items():
            monomials.add(monomial)
            coeffs.add(coeff)

    monom_gcd = monomial_ngcd(list(monomials))
    coeff_gcd = domain.gcdn(coeffs)

    return {monom_gcd: coeff_gcd}


def smp_gcd_prs(A: smp[Er], B: smp[Er], n: int, domain: Domain[Er],
data: tuple[Callable[..., smp[Er]], bool, int] | None = None) -> smp[Er]:
    """
    Returns the greatest common divisor (GCD) of two polynomials using the
    Polynomial Resultant Sequences (PRS) method.

    Examples
    ========

    >>> from sympy.polys.rings import gcd_prs
    >>> from sympy import ZZ, ring
    >>> R, x, y = ring("x, y", ZZ)

    >>> f = 4*x**2 + 8*x + 10*y
    >>> g = 2*x**3 + 4*x**2*y + 2*x*y**2
    >>> gcd_prs(f, g)
    2

    """
    if data:
        f, _, _ = data
    map = [i for i in range(1,n)]
    c1, pp1 = smp_primitive_wrt(A, 0, n, domain, data)
    c2, pp2 = smp_primitive_wrt(B, 0, n, domain, data)

    c = _gcd_list([c1, c2], n-1, domain, data)

    # h = pp1.subresultants(pp2, x)[-1]
    h = smp_subresultants(pp1, pp2, 0, n, domain)[-1]
    _, h = smp_primitive_wrt(h, 0, n, domain, data)

    if not domain.is_Field:
        k = domain.canonical_unit(smp_LC(h, n, domain))
        _smp_imul_ground(c, k, n-1, domain)
    map = [i for i in range(1,n)]
    h = smp_mul(h, smp_elevate(c, map, n, domain), domain, n)

    if domain.is_Field:
        # h = h.monic()
        k = smp_LC(h, n, domain)
        _smp_imul_ground(h, domain.revert(k), n, domain)

    return h


def gcd_list(polys: list[smp[Er]], n: int, domain: Domain[Er],
order: MonomialOrder = lex) -> smp[Er]:
    """
    changes the domain of the inputs in one managable by the sympy algorithms,
    calls _gcd_list for the actual gcd computation,
    restores the original domain and normalizes the output
    """
    if not domain.is_Exact:
        domain_ex = domain.get_exact()
        polys_ex = [smp_domain_convert(pol, domain, domain_ex, n) for pol in polys]
        gcd_ex = gcd_list(polys_ex, n, domain_ex, order)
        gcd = smp_domain_convert(gcd_ex, domain_ex, domain, n)

    elif domain.is_QQ:
        polys = [smp_clear_denoms(pol, n, domain)[1] for pol in polys]
        polys_z = [smp_domain_convert(pol, domain, ZZ, n) for pol in polys]
        gcd_z = _gcd_list(polys_z, n, ZZ)
        gcd = smp_domain_convert(gcd_z, ZZ, domain, n)

    else:
        gcd = _gcd_list(polys, n, domain)

    return smp_normalize_gcd(gcd, n, domain, order)


def smp_normalize_gcd(gcd: smp[Er], n: int, domain: Domain[Er],
order: MonomialOrder = lex) -> smp[Er]:
    # normalizes the gcd
    if not gcd:
        return {}

    if domain.is_Field:
        k = smp_LC(gcd, n, domain, order)
        gcd = smp_mul_ground(gcd, domain.revert(k), n, domain)

    else:
        k = domain.canonical_unit(smp_LC(gcd, n, domain, order))
        gcd = smp_mul_ground(gcd, k, n, domain)
    return gcd


def _gcd_list(polys: list[smp[Er]], n: int, domain: Domain[Er],
data: tuple[Callable[..., smp[Er]], bool, int] | None = None) -> smp[Er]:
    """
    Returns the greatest common divisor (GCD) of a list of polynomials.

    Examples
    ========

    >>> from sympy.polys.rings import gcd_extract
    >>> from sympy import ZZ, ring
    >>> R, x, y = ring("x, y", ZZ)

    >>> polynomials = x - y, x - y
    >>> gcd_extract(polynomials, "heugcd")
    x - y

    """
    # l'idea è di rendere questa privata, e creare una funzione gcd_list
    # che renda esatti i domini, e che delega la scelta dell'algoritmo a algorithm_selector
    polys = [polynomial for polynomial in polys if polynomial]

    if not polys:
        return {}

    if len(polys) == 1:
        return polys[0]

    if any(len(pi) == 1 for pi in polys):
        return gcd_terms(polys, domain)

    polys, mon_gcd = smp_monomial_extract(polys, n, domain)

    polys, common_var = _gcd_preprocess_polys(polys, n, domain)

    # quick exits
    if common_var is None:
        gcd = polys[0]
        if mon_gcd is not None:
            return smp_mul(gcd, mon_gcd, domain, n)
        return gcd

    m = len(common_var)

    if len(polys) == 1:
        gcd = polys[0]
        gcd = smp_elevate(gcd, common_var, n, domain)
        if mon_gcd is not None:
            return smp_mul(gcd, mon_gcd, domain, n)
        return gcd

    if data is None:
        f_, comb, p = gcd_algorithm_selector(polys, m, domain)
    else:
        f_, comb, p = data

    # forse ha più senso fare restituire a gcd_algorithm_slector un bool
    # sulla base del quale decidere se usare le combinazioni lineari oppure l'iterazione
    if comb:
        gcd = smp_gcd_list_lin_comb(polys, p, m, domain, f_, data)

    else:
        gcd = smp_gcd_list_iterative(polys, f_, m, domain, data)

    gcd = smp_elevate(gcd, common_var, n, domain)
    if mon_gcd is not None:
        return smp_mul(gcd, mon_gcd, domain, n)
    return gcd


def smp_gcd_list_iterative(polys: list[smp[Er]], f: Callable[..., smp[Er]],
n: int, domain: Domain[Er], data: tuple[Callable[..., smp[Er]], bool, int] | None = None
) -> smp[Er]:
    gcd = polys[0]
    for pol in polys[1:]:
        gcd = f(gcd, pol, n, domain, data)
        if smp_is_one(gcd, n, domain):
            break
    return gcd


def gcd_algorithm_selector(polys: list[smp[Er]], n: int, domain: Domain[Er]
) -> tuple[Callable[..., smp[Er]], bool, int]:
    ...


def adapter_heugcd(A: smp[MPZ], B: smp[MPZ], n: int, domain: Domain[MPZ],
data: tuple[Callable[..., smp[MPZ]], bool, int] | None = None) -> smp[MPZ]:
    try:
        gcd, _, _ =  smp_heugcd(A, B, n)
    except HeuristicGCDFailed:
        gcd, _, _ = smp_zippel_gcd(A, B, n, data)
    return gcd


def adapter_zippel(A: smp[MPZ], B: smp[MPZ], n: int, domain: Domain[MPZ],
data: tuple[Callable[..., smp[MPZ]], bool, int] | None = None) -> smp[MPZ]:
    gcd, _, _ = smp_zippel_gcd(A, B, n, data)
    return gcd


def adapter_zippel_mod(A: smp[Er], B: smp[Er], n: int, domain: Domain[Er],
data: tuple[Callable[..., smp[Er]], bool, int] | None = None) -> smp[Er]:
    gcd = smp_zippel_gcd_mod(A, B, n, domain)
    if gcd is not None:
        return gcd
    return smp_gcd_prs(A, B, n, domain, data)


def smp_compress_polys(polys: list[smp[Er]], var: list[int], n: int, dom: Domain[Er]
) -> list[smp[Er]]:
    var.sort()
    compr_polys = []
    for pol in polys:
        compr_pol = {}
        for mon, coeff in pol.items():
            compr_mon = tuple(mon[i] for i in var)
            compr_pol[compr_mon] = coeff
        compr_polys.append(compr_pol)
    return compr_polys


def _gcd_preprocess_polys(polynomials: list[smp[Er]], n: int, dom: Domain[Er]
) -> tuple[list[smp[Er]], list[int] | None]:
    """
    Simplify a list of polynomials whose gcd is wanted. Returns a possibly
    longer list of simpler polynomials having the same gcd as the input.

    This is done by eliminating symbols that can not be part of the gcd
    because they do not appear in each item of the input. In the output
    list, all items have exactly the same symbols. The set of those symbols
    is also returned.

    Examples
    ========

    >>> from sympy.core import ordered
    >>> from sympy.polys.rings import _gcd_preprocess_polys
    >>> from sympy import ring, ZZ
    >>> R, x, y = ring("x, y", ZZ)

    >>> f = x**2 - y**2
    >>> g = x**2 - 2*x*y + y**2
    >>> h = x - y
    >>> polynomials = [f, g, h]
    >>> result = _gcd_preprocess_polys(polynomials, 2, dom)
    >>> list(ordered(result[0])), result[1] # Set ordering is random
    ([x - y, x**2 - y**2, x**2 - 2*x*y + y**2], {0, 1})

    """
    all_polys = polynomials
    common_global = list(range(n))
    m = n
    while True:
        m = len(common_global)
        # Quick exits are most efficient if we start from the simplest polys
        polynomials = sorted(smp_unique_polys(all_polys, m, dom), key=len)

        # Find the intersection of symbols for each poly:
        common = smp_active_var(polynomials[0], m, dom)
        allsame = True
        for pol in polynomials[1:]:

            syms = smp_active_var(pol, m, dom)
            if allsame and syms != common:
                allsame = False
            common &= syms

            # Quick exit
            if not common:
                gcd = dom.gcdn(coeff for polynomial in polynomials for coeff in polynomial.values())
                return [{(0,)*n: gcd}], None

        common_global = [var for i, var in enumerate(common_global) if i in common]
        # The loop is complete if they all have the same symbols.
        if allsame:
            polynomials = smp_compress_polys(polynomials, list(common), m, dom)
            return polynomials, common_global

        # Extract coefficients as polys containing only the common symbols.
        all_polys = []
        rem_var = list(set(range(m)) - common)
        for i, pol in enumerate(polynomials):
            r_pol = smp_reorg_poly(pol, rem_var, m, dom)
            all_polys.extend(r_pol.values())

        # Quick exit
        if any(len(c) == 1 for c in all_polys):
            gcd = gcd_terms(all_polys, dom)
            return [gcd], common_global


def smp_primitive_wrt(g: smp[Er], i: int, n: int, domain: Domain[Er],
data: tuple[Callable[..., smp[Er]], bool, int] | None = None) -> tuple[smp[Er], smp[Er]]:
    """
    Returns the content and primitive part of a polynomial with respect to a
    specified variable.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x, y = ring("x, y", ZZ)
    >>> p = 6*x**2*y - 9*x**3*y**2 + 3*x*y
    >>> p.primitive_wrt(x)
    (3*y, -3*x**3*y + 2*x**2 + x)

    >>> p.primitive() # Distinguishing primitive and primitive_wrt outcomes
    (3, -3*x**3*y**2 + 2*x**2*y + x*y)
    """
    g_r = smp_reorg_poly(g, [i], n, domain)
    cont = _gcd_list(list(g_r.values()), n-1, domain, data)
    map = [j for j in range(n) if j != i]
    cont_ = smp_elevate(cont, map, n, domain)
    [prim], _ = smp_div_list(g, [cont_], n, domain)
    return cont, prim


def smp_zippel_gcd(A: smp[MPZ], B: smp[MPZ], n: int,
data: tuple[Callable[..., smp[MPZ]], bool, int] | None = None
) -> tuple[smp[MPZ], smp[MPZ], smp[MPZ]]:
    """
    checks if in any variable the gcd is monic
    if it doesn't find any variable with monic gcd
    it extracts content wrt the first variable
    and proceeds to compute the gcd
    maybe it is better to call the generic sympy gcd rather than
    modgcd_multivariate once the whole dispatching heuristic is implemented
    """
    for i in range(n):
        m_A = max(el[i] for el in A.keys())
        m_B = max(el[i] for el in B.keys())
        lc_A = smp_coeff_wrt(A, i, m_A, n, ZZ)
        lc_B = smp_coeff_wrt(B, i, m_B, n, ZZ)
        if len(lc_A) == 1 or len(lc_B) == 1:
            if i == 0:
                return _smp_zippel_gcd(A, B, n)
            A_ = smp_swap_var(A, i, n, ZZ) # swap non scambia le variabili, crea nuova funzione
            B_ = smp_swap_var(B, i, n, ZZ)
            gcd, cff, cfg = _smp_zippel_gcd(A_, B_, n)
            return smp_swap_var(gcd, i, n, ZZ), smp_swap_var(cff, i, n, ZZ), smp_swap_var(cfg, i, n, ZZ)

    A_cont, A_prim = smp_primitive_wrt(A, 0, n, ZZ, data)
    B_cont, B_prim = smp_primitive_wrt(B, 0, n, ZZ, data)
    cont_gcd = _gcd_list([A_cont, B_cont], n-1, ZZ, data)
    gcd, cff, cfg = _smp_zippel_gcd(A_prim, B_prim, n)
    cont_gcd = smp_elevate(cont_gcd, [el for el in range(1, n)], n, ZZ)
    return smp_mul(cont_gcd, gcd, ZZ, n), cff, cfg # wrong cofactors, multiply by content/content_gcd


def smp_zippel_gcd_mod(A: smp[Er], B: smp[Er], n: int, domain: Domain[Er],
data: tuple[Callable[..., smp[Er]], bool, int] | None = None) -> smp[Er] | None:
    """
    checks if in any variable the gcd is monic
    if it doesn't find any variable with monic gcd
    it extracts content wrt the first variable
    and proceeds to compute the gcd
    maybe it is better to call the generic sympy gcd rather than
    modgcd_multivariate once the whole dispatching heuristic is implemented
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
                if gcd == None:
                    return None
                return smp_domain_convert(gcd, ZZ, domain, n)
            A_z_ = smp_swap_var(A_z, i, n, ZZ)
            B_z_ = smp_swap_var(B_z, i, n, ZZ)
            gcd = _smp_zippel_gcd_mod(A_z_, B_z_, p, n)
            if gcd == None:
                return None
            gcd = smp_swap_var(gcd, i, n, ZZ)
            return smp_domain_convert(gcd, ZZ, domain, n)

    Ac, Ap = smp_primitive_wrt(A, 0, n, domain, data)
    Bc, Bp = smp_primitive_wrt(B, 0, n, domain, data)
    Ap_z = smp_domain_convert(Ap, domain, ZZ, n)
    Bp_z = smp_domain_convert(Bp, domain, ZZ, n)
    cont_gcd = _gcd_list([Ac, Bc], n-1, domain, data)
    gcdp = _smp_zippel_gcd_mod(Ap_z, Bp_z, p, n)
    if gcdp == None:
        return None
    cont_gcd = smp_elevate(cont_gcd, [el for el in range(1, n)], n, ZZ)
    gcdp = smp_domain_convert(gcdp, ZZ, domain, n)
    return smp_mul(cont_gcd, gcdp, domain, n)


def smp_elevate(f, map: list[int], n, domain):
    # adds variables to a polynomial
    h = {}
    for mon, coeff in f.items():
        m = [0]*n
        for i, exp in zip(map, mon):
            m[i] = exp
        h[tuple(m)] = coeff
    return h


def smp_gcd_list_lin_comb(polys: list[smp[Er]], t: int,
n: int, domain: Domain[Er], f: Callable[..., smp[Er]],
data: tuple[Callable[..., smp[Er]], bool, int] | None = None
) -> smp[Er]:
    # faster than computing the gcd iteratively
    # when gcd is computed with Zippel

    if len(polys) == 1:
        return polys[0]
    if len(polys) == 2:
        gcd = f(polys[0], polys[1], n, domain, data)
        return gcd
    v_sor = sorted(polys, key=lambda p: len(p))
    while True:
        B = v_sor[1].copy()
        for pol in v_sor[2:]:
            k = domain(random.randint(1, t))
            smp_mul_ground(pol, k, n, domain)
            smp_add(B, pol, domain, n)
        gcd = f(v_sor[0], B, n, domain, data)
        restart = False
        for pol in polys:
            if smp_rem_list(pol, [gcd], n, domain) != {}:
                restart = True
                break
        if not restart:
            break
    return gcd
