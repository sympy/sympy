"""Low-level Zippel algorithm functions on raw dict representations.

These functions deliberately do not depend on PolyRing or PolyElement.
They operate on dicts mapping exponent tuples to coefficients and use
only the provided domain and number of variables.
"""

from __future__ import annotations

from collections import defaultdict
from random import randint
from typing import TYPE_CHECKING

from sympy.ntheory.generate import nextprime
from sympy.polys.densebasic import dup_from_dict, dup_to_dict
from sympy.polys.domains.finitefield import FF
from sympy.polys.domains.integerring import ZZ
from sympy.polys.galoistools import gf_add, gf_div, gf_eval, gf_gcd, gf_mul, gf_quo, gf_from_dict
from sympy.ntheory.modular import crt
from sympy.polys.matrices.domainmatrix import DomainMatrix
from sympy.polys.matrices.exceptions import DMNonInvertibleMatrixError
from sympy.polys.sparsetools import (
    _smp_itrunc_ground,
    smp_coeff_wrt,
    smp_degree,
    smp_degrees,
    smp_div_list,
    smp_evaluate,
    smp_is_ground,
    smp_LC,
    smp_is_one,
    smp_leading_expv,
    smp_mul,
    smp_mul_ground,
    smp_primitive,
    smp_quo_ground,
    smp_rem_list,
    smp_subs_drop,
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


def smp_swap(f, i, n, domain):
    """
    Make the variable `x_i` the leading one in a multivariate polynomial `f`.
    """
    fswap = {}
    for monom, coeff in f.items():
        monomswap = (monom[i],) + monom[:i] + monom[i+1:]
        fswap[monomswap] = coeff
    return fswap


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


def lag_basis(evalpoints: list[MPZ], p: MPZ) -> list[dup[MPZ]]:
    r"""
    Computes the Lagrange basis associated to the given
    list of evaluation points over `\mathbb{Z}_p`.

    Example
    =======

    >>> from sympy.polys.zippel import lag_basis
    >>> evalpoints = [1, 2, 5]
    >>> p = 11
    >>> lag_basis(evalpoints, p)
    [[3, 1, 8], [7, 2, 2], [1, 8, 2]]

    References
    ==========

    .. [1] [Yang09]_

    """
    master_pol = [ZZ.one]
    for k in evalpoints:
        master_pol = gf_mul(master_pol, [ZZ.one, -k], p, ZZ)

    v = [gf_div(master_pol, [ZZ.one, -k], p, ZZ)[0] for k in evalpoints]
    for i, poly in enumerate(v):
        c = gf_eval(poly, evalpoints[i], p, ZZ)
        c = ZZ.invert(c, p)
        v[i] = [(c * el) % p for el in v[i]]
    return v


def vandermonde_interp(
    basis: list[dup[MPZ]], values: list[MPZ], p: MPZ, trans: bool = True
) -> list[MPZ]:
    r"""
    Solves for `x` the linear system `A^T x = \mathrm{values}` in
    `\mathbb{Z}_p`, where `A` is the Vandermonde matrix such that
    `a_{i,j} = k_i^j`.
    If `\mathrm{trans} = \mathrm{False}`, it solves the linear system
    `A x = \mathrm{values}`.
    Both systems are solved in `O(n^2)` time.

    The parameter basis can be computed using
    `\operatorname{lag\_basis}(k, p)`.

    Note that solving a linear system with associated Vandermonde matrix is
    equivalent to interpolating an univariate polynomial.
    Note also that the interpolation of a multivariate polynomial in
    `\mathbb{Z}_p[x_1, \ldots, x_n]`, can be made equivalent to the resolution
    of a linear system with the transposed of a Vandermonde matrix as the
    associated matrix.
    It is sufficient to choose a random tuple `(t_1, \ldots, t_n)`,
    and to use `(t_1^j, \ldots, t_n^j)` for
    `j \in \{0, 1, \ldots\}` as evaluation points.

    Example
    =======

    >>> from sympy.polys.zippel import vandermonde_interp, lag_basis
    >>> from sympy import Matrix

    >>> k = [1, 2, 5]
    >>> p = 9973
    >>> values = [4, 9, 36]
    >>> A = Matrix([[1, 1, 1], [1, 2, 4], [1, 5, 25]])

    >>> bas = lag_basis(k, p)

    >>> x = vandermonde_interp(bas, values, p, trans = False)

    >>> A * Matrix(x) == Matrix(values)
    True

    References
    ==========

    .. [1] [Yang09]_

    """
    deg = len(values) - 1
    if trans:
        sol = []
        for poly in basis:
            sol.append(sum(val * poly[deg - i] for i, val in enumerate(values)) % p)
    else:
        sol = []
        for i in range(len(values)):
            sol.append(
                sum(values[j] * basis[j][deg - i] for j in range(len(values))) % p
            )

    return sol


def incremental_newton_interp(
    x: list[MPZ], v: list[MPZ], xk: MPZ, uk: MPZ, p: MPZ
) -> MPZ:
    r"""
    Computes the next Newton interpolation coefficient `v_k` over
    `\mathbb{Z}_p`.

    Given a list of `k` evaluation points
    `x = [x_0, \ldots, x_{k-1}]` and the corresponding coefficients
    `v = [v_0, \ldots, v_{k-1}]` of the Newton interpolation polynomial
    `P_{k-1}(x)`, this function calculates the new coefficient `v_k` required
    to interpolate a new point `(x_k, u_k)`.

    References
    ==========

    .. [1] [Geddes92]_ p. 186.

    Examples
    ========

    >>> from sympy.polys.zippel import incremental_newton_interp
    >>> x = [0, 1]
    >>> v = [67, 47]
    >>> xk = 2
    >>> uk = 66
    >>> p = 97
    >>> v2 = incremental_newton_interp(x, v, xk, uk, p)
    >>> v2
    1

    """
    s = v[0]
    c = ZZ.one
    for i, el in enumerate(v[1:]):
        c = (c * (xk - x[i])) % p
        s = (s + el * c) % p
    inv = ZZ.invert(c * (xk - x[len(v) - 1]), p)
    return ((uk - s) * inv) % p


def from_newt_to_poly(x: list[MPZ], v: list[MPZ], p: MPZ) -> dup[MPZ]:
    r"""
    Given `k` evaluation points `x = [x_0, \ldots, x_{k-1}]` and the
    corresponding coefficients `v = [v_0, \ldots, v_{k-1}]` of the Newton
    interpolation polynomial `P_{k-1}(x)`, this function calculates explicitly
    and returns the interpolation polynomial `P_{k-1}(x)` over
    `\mathbb{Z}_p`.

    References
    ==========

    .. [1] [Geddes92]_ p. 186.

    Examples
    ========

    >>> from sympy.polys.zippel import from_newt_to_poly
    >>> x = [0, 1, 2]
    >>> v = [67, 47, 1]
    >>> p = 97
    >>> from_newt_to_poly(x, v, p)
    [1, 46, 67]

    """
    pol = [v[-1]]
    for i in range(len(v) - 2, -1, -1):
        binomial = [ZZ.one, -x[i]]
        pol = gf_mul(pol, binomial, p, ZZ)
        pol = gf_add(pol, [v[i]], p, ZZ)

    return pol


def skeleton_sorter(
    G: dict[monom, MPZ],
) -> tuple[dict[int, list[tuple[monom, list[tuple[int, int]]]]], list[list[MPZ]], bool, bool]:
    r"""
    Reorganizes the skeleton of a sparse polynomial for multivariate interpolation.

    This function extracts the monomials of a sparse polynomial.
    It groups the monomials by the degree of the first variable,
    sorts these groups by the number of monomials they contain
    (in ascending order), and builds a compact representation of the non-zero
    degrees for the remaining variables, useful to evaluate quickly the monomials.

    It also saves the coefficients of each monomial in a list of lists, and
    2 booleans needed by the function ``zippel_interp()`` to choose the right routine.

    Examples
    ========

    If we take the sparse polynomial
    `G = 5x^2y^2 + 7xy^5z^3 + 8xz^4` in `\mathbb{Z}[x, y, z]`,
    ``skeleton_sorter(G)`` returns the following outputs:

    .. code-block:: python

        S = {
            2: [((2, 2, 0), [(0, 2)])],
            1: [((1, 5, 3), [(0, 5), (1, 3)]), ((1, 0, 4), [(1, 4)])]
        }

        h = [[5], [7], [8]]

        monic = True
        pseudomonic = True

    """
    S_ = defaultdict(list)
    for mon, _ in G.items():
        S_[mon[0]].append((mon, []))
    S = {deg: S_[deg] for deg in sorted(S_.keys(), key=lambda x: len(S_[x]))}
    h = []
    lc = S[max(S.keys())]

    monic = False
    pseudomonic = False
    if len(lc) == 1:
        monic = True
        for el in lc[0][0][1:]:
            if el != 0:
                pseudomonic = True

    for _, mons in S.items():
        for mon_tuple in mons:
            h.append([G[mon_tuple[0]]])
            for i, el in enumerate(mon_tuple[0][1:]):
                if el != 0:
                    mon_tuple[1].append((i, el))

    return S, h, monic, pseudomonic


def zippel_gcd(f: smp[MPZ], g: smp[MPZ], n: int) -> tuple[smp[MPZ], smp[MPZ], smp[MPZ]]:
    r"""
    Compute the GCD of two polynomials in `\mathbb{Z}[x_0, \ldots, x_{k-1}]`
    using the Zippel algorithm.

    The algorithm computes the GCD of two multivariate integer polynomials
    `f` and `g` by calculating the GCD in
    `\mathbb{Z}_p[x_0, \ldots, x_{k-1}]` for suitable primes `p` and then
    reconstructing the coefficients with the Chinese Remainder Theorem. To
    compute the multivariate GCD over `\mathbb{Z}_p` the recursive
    subroutine :func:`sparse_gcd` is used. To verify the result in
    `\mathbb{Z}[x_0, \ldots, x_{k-1}]`, trial division is done, but only for
    candidates which are very likely the desired GCD.

    Parameters
    ==========

    f : PolyElement
        multivariate integer polynomial
    g : PolyElement
        multivariate integer polynomial

    Returns
    =======

    h : PolyElement
        GCD of the polynomials `f` and `g`
    cff : PolyElement
        cofactor of `f`, i.e. `\frac{f}{h}`
    cfg : PolyElement
        cofactor of `g`, i.e. `\frac{g}{h}`

    Examples
    ========

    >>> from sympy.polys.zippel import zippel_gcd
    >>> from sympy.polys import ring, ZZ

    >>> R, x, y = ring("x, y", ZZ)

    >>> f = x**2 - y**2
    >>> g = x**2 + 2*x*y + y**2

    >>> h, cff, cfg = zippel_gcd(f, g, 2)
    >>> h, cff, cfg
    (x + y, x - y, x + y)

    >>> cff * h == f
    True
    >>> cfg * h == g
    True

    >>> R, x, y, z = ring("x, y, z", ZZ)

    >>> f = x*z**2 - y*z**2
    >>> g = x**2*z + z

    >>> h, cff, cfg = zippel_gcd(f, g, 3)
    >>> h, cff, cfg
    (z, x*z - y*z, x**2 + 1)

    >>> cff * h == f
    True
    >>> cfg * h == g
    True

    References
    ==========

    References
    ==========

    .. [1] [Geddes92]_
    .. [2] [Yang09]_

    See also
    ========

    sparse_gcd

    """

    result = smp_trivial_gcd(f, g, n, ZZ)
    if result is not None:
        return result

    # divide out integer content
    cf, f = smp_primitive(f, n, ZZ)
    cg, g = smp_primitive(g, n, ZZ)
    ch = ZZ.gcd(cf, cg)


    gamma = ZZ.gcd(smp_LC(f, n, ZZ), smp_LC(g, n, ZZ))

    badprimes = ZZ.one
    for i in range(n):
        sf = smp_swap(f, i, n, ZZ)
        sg = smp_swap(g, i, n, ZZ)
        badprimes *= ZZ.gcd(smp_LC(sf, n, ZZ), smp_LC(sg, n, ZZ))

    # empirical evidence shows that brown is faster using one
    # while Zippel needs large primes to not run out of evaluation points
    m = ZZ.one

    p = 10**9

    hlastm = {}
    while True:
        p = ZZ(nextprime(p))
        while badprimes % p == 0:
            p = ZZ(nextprime(p))

        fp = smp_trunc_ground(f, p, n, ZZ)
        gp =  smp_trunc_ground(g, p, n, ZZ)

        hp = sparse_gcd(fp, gp, p, n)
        if hp is not None:
            if hlastm != {}:
                if smp_degrees(hp, n, ZZ) > smp_degrees(hlastm, n, ZZ):
                    continue
                elif smp_degrees(hp, n, ZZ) < smp_degrees(hlastm, n, ZZ):
                    hlastm = {}
                    m = ZZ.one

            hp_k = ZZ.invert(smp_LC(hp, n, ZZ), p)
            for mon, coeff in hp.items():
                hp[mon] = coeff * hp_k
                if not hp[mon]:
                    del hp[mon]
        else:
            continue

        hp = smp_mul_ground(hp, gamma, n, ZZ)
        hp = smp_trunc_ground(hp, p, n, ZZ)
        if m == 1:
            m = p
            hlastm = hp
            continue

        hm = smp_chinese_remainder_reconstruction_multivariate(hp, hlastm, p, m, ZZ, n)
        m *= p

        if not hm == hlastm:
            hlastm = hm
            continue

        _, h = smp_primitive(hm, n, ZZ)
        fquo, frem = smp_div_list(f, [h], n, ZZ)
        gquo, grem = smp_div_list(g, [h], n, ZZ)
        if not frem and not grem:
            if smp_LC(h, n, ZZ) < 0:
                ch = -ch
            h = smp_mul_ground(h, ch, n, ZZ)

            cff = smp_mul_ground(fquo[0], cf // ch, n, ZZ)
            cfg = smp_mul_ground(gquo[0], cg // ch, n, ZZ)
            return h, cff, cfg


def sparse_gcd(A: smp[MPZ], B: smp[MPZ], p: MPZ, n: int
    ) -> smp[MPZ] | None:
    r"""
    Compute the GCD of two multivariate polynomials over a finite field.

    The polynomials ``A`` and ``B`` are represented by sparse dictionaries with
    coefficients in `\mathbb{Z}` reduced modulo the prime `p`. The computation is
    performed recursively by evaluating the last variable, computing lower
    dimensional GCDs, and recovering the missing variable by Newton
    interpolation.

    Parameters
    ==========

    A, B :
        Polynomials over `\mathbb{F}_p`.
    p :
        Prime modulus.
    n :
        Number of variables.

    Returns
    =======

    dict or None
        smp representation of the GCD of ``A`` and ``B`` modulo `p`.
        A return value of ``None`` indicates a failed probabilistic computation for
        the current prime and evaluation points

    References
    ==========

    .. [1] [Geddes92]_
    .. [2] [Yang09]_

    """
    if n == 1:
        G = smp_gf_gcd(A, B, p, ZZ)
        return G
    a, A = smp_primitive_wrt_last(A, n, ZZ, p)
    b, B = smp_primitive_wrt_last(B, n, ZZ, p)

    c = gf_gcd(a, b, p, ZZ)   # here c is a list poly
    # here c is elevated to a dict poly with all variables
    c = {(0,)*(n-1) + (i,): coeff for i, coeff in enumerate(c[::-1]) if coeff != 0}

    g = smp_gf_gcd(smp_LC_wrt_last(A, n, ZZ), smp_LC_wrt_last(B, n, ZZ), p, ZZ)
    M = []
    h = []
    skippable = set()
    pt_chances = 0
    while True:
        t = ZZ(randint(1, int(p - 1)))
        gk = smp_evaluate(g, {0:t}, 1, ZZ) %p
        if t in M:
            continue
        if gk == 0:
            pt_chances +=1
            if pt_chances < 3:
                continue
            else:
                return None
        # here tuples of the dicts are reduced by one
        A_ = smp_subs_drop(A, {n - 1: t}, n, ZZ)
        _smp_itrunc_ground(A_, p, n - 1, ZZ)

        B_ = smp_subs_drop(B, {n - 1: t}, n, ZZ)
        _smp_itrunc_ground(B_, p, n - 1, ZZ)
        G = sparse_gcd(A_, B_, p, n-1)
        if G == None:
            M = []
            h = []
            skippable = set()
            pt_chances += 1
            if pt_chances < 3:
                continue
            else:
                return None
        else:
            if smp_is_one(G, n - 1, ZZ):
                return c
            G_k = ZZ.invert(smp_LC(G, n-1, ZZ), p)

            M.append(t)
            # normalization
            for mon in G:
                G[mon] = ((G[mon] * G_k) * gk)
            G = smp_trunc_ground(G, p, n-1, ZZ)

            G_s, h, monic, pseudomonic = skeleton_sorter(G)
            initial_M = M.copy()
            initial_h = [coeffs.copy() for coeffs in h]
        skeleton_chances = 0
        while True:
            t = ZZ(randint(1, int(p - 1)))
            gk = smp_evaluate(g, {0:t}, 1, ZZ) %p
            if t in M or gk == 0:
                continue

            A_ = smp_subs_drop(A, {n-1: t}, n, ZZ)
            B_ = smp_subs_drop(B, {n-1: t}, n, ZZ)
            _smp_itrunc_ground(A_, p, n-1, ZZ)
            _smp_itrunc_ground(B_, p, n-1, ZZ)
            G_ = zippel_interp(A_, B_, G_s, p, monic, pseudomonic, n-1)
            if G_ == None:
                skeleton_chances += 1
                if skeleton_chances < 3:
                    continue
                else:
                    return None

            G_k = ZZ.invert(smp_LC(G_, n-1, ZZ), p)# inverse of leading coeff
            for mon in G_:
                G_[mon] = ((G_[mon] * G_k) * gk) %p
            repeat = False

            for i, (_, coeff) in enumerate(G_.items()):
                # if a coefficient is zero it means the poly has reached its deg:
                # every coefficient from that point on will (very likely) be zero
                if i in skippable:
                    continue
                vk = incremental_newton_interp(M[:len(h[i])], h[i], t, coeff, p)
                if vk != 0:
                    repeat = True
                    h[i].append(vk)
                else:
                    skippable.add(i)
            M.append(t)
            if not repeat:
                gcd = {}
                for i, mon in enumerate(G_.keys()):
                    pol = from_newt_to_poly(M[:len(h[i])], h[i], p)[::-1]
                    for j, el in enumerate(pol):
                        gcd[mon + (j,)] = el

                _, gcd = smp_primitive_wrt_last(gcd, n, ZZ, p)

                # modular conversion to perform test divisions
                A_mod = {}
                for mon, coeff in A.items():
                    if coeff % p != 0:
                        A_mod[mon] = FF(p)(coeff)
                B_mod = {}
                for mon, coeff in B.items():
                    if coeff % p != 0:
                        B_mod[mon] = FF(p)(coeff)
                gcd_mod = {}
                for mon, coeff in gcd.items():
                    if coeff % p != 0:
                        gcd_mod[mon] = FF(p)(coeff)

                if smp_rem_list(A_mod, [gcd_mod], n, FF(p)) == {} and smp_rem_list(B_mod, [gcd_mod], n, FF(p)) == {}:
                    gcd = smp_mul(gcd, c, ZZ, n)
                    _smp_itrunc_ground(gcd, p, n, ZZ)
                    return gcd
                else:
                    skeleton_chances += 1
                    if skeleton_chances < 3:
                        M = initial_M.copy()
                        h = [coeffs.copy() for coeffs in initial_h]
                        skippable = set()
                    else:
                        return None


def zippel_interp(A: smp[MPZ], B: smp[MPZ],
    G: dict[int, list[tuple[monom, list[tuple[int, int]]]]], p: MPZ, monic: bool,
    pseudomonic: bool, n: int) -> smp[MPZ] | None:

    r"""
    Recover a sparse modular GCD from a prescribed skeleton.

    Given two polynomials over
    `\mathbb{F}_p[x_0, \ldots, x_{n-1}]`, this function reconstructs the
    coefficients of their GCD using sparse interpolation. The expected
    monomial support is supplied by ``G`` and is verified against the specialized
    univariate GCDs computed during the interpolation.

    Parameters
    ==========

    A, B :
        Polynomials over `\mathbb{F}_p`.
    G :
        Skeleton of the expected GCD, grouped by degree in `x_0`.
    p :
        Prime modulus.
    monic :
        Whether the highest power of `x_0` appears only in
        one monomial.
    pseudomonic :
        Whether the monomial with the highest power of `x_0`
        depends on variables other than `x_0`.
    n :
        Number of variables in ``A`` and ``B``.

    Returns
    =======

    dict or None
        The interpolated polynomial with support prescribed by ``G``, or
        ``None`` if the skeleton is incompatible with the computed
        specializations or the interpolation cannot be completed.

    References
    ==========

    .. [1] [Zippel79]
    .. [2] [Yang09]_

    """
    deg_a = smp_degree(A, 0, n, ZZ)
    lc_A = smp_coeff_wrt(A, 0, deg_a, n, ZZ)
    lengths = list(len(el) for el in G.values())
    nt = lengths[-1]
    # base case: for one variable its sufficient to take the gcd's coefficients
    # the function runs also without this base case, but it allows for an early exit
    if n == 1:
        C = {}
        gcd = smp_gf_gcd(A, B, p, ZZ)
        if smp_degree(gcd, 0, 1, ZZ) != max(G):
                return None
        for mon in gcd:
            if mon[0] not in G:
                return None
        for deg, monoms in G.items():
            C[monoms[0][0]] = gcd.get((deg,), 0)
        return C

    for _ in range(3):
        while True:
            t = tuple(ZZ(randint(int(1), int(p-1))) for _ in range(n-1))
            is_bad_tuple = False
            lc_A_ev = {}
            all_vand_basis = []
            # checking that the leading coeff of A isn't zero in the chosen tuple
            for mon, coeff in lc_A.items():
                j = ZZ.one
                for i, k in enumerate(mon[1:]):
                    j = (j *pow(t[i], k, p)) %p
                lc_A_ev[j] = (lc_A_ev.get(j, 0) + coeff) %p
            for i in range(nt):
                j = 0
                for val, coeff in lc_A_ev.items():
                    j = (j + pow(val, i, p) * coeff) %p
                if j == 0:
                    is_bad_tuple = True
                    break
            if is_bad_tuple:
                continue

            for deg, el in G.items():
                vand_bas = []
                for mon in el:
                    j = ZZ.one
                    for i in mon[1]:
                        j =  (j * pow(t[i[0]], i[1], p)) %p
                    if j in vand_bas:
                        is_bad_tuple = True
                        break

                    vand_bas.append(j)
                if is_bad_tuple:
                    break
                all_vand_basis.append(vand_bas)
            if is_bad_tuple:
                continue
            else:
                break

        deg_b = smp_degree(B, 0, n, ZZ)
        A_flat = list({} for _ in range(deg_a +1))
        B_flat = list({} for _ in range(deg_b +1))

        for mon, coeff in A.items():
            j = 1
            for i, k in enumerate(mon[1:]):
                if k != 0:
                    j = (j * pow(t[i], k, p)) %p
            A_flat[mon[0]][j] = (A_flat[mon[0]].get(j, 0) + coeff) %p

        for mon, coeff in B.items():
            j = 1
            for i, k in enumerate(mon[1:]):
                if k != 0:
                    j = (j * pow(t[i], k, p)) %p
            B_flat[mon[0]][j] = (B_flat[mon[0]].get(j, 0) + coeff) %p

        # represented A, B as lists (list index is corresponding power of x1)
        # containing dicts with as keys the evaluated monomials
        # and as values the related coeffs. example: {t1: c1,..., tn: cn}
        # so that to evaluate in powers of the tuple I can do a linear combination t1**k*c1+...+tn**k*cn

        eval_points = {deg:[] for deg in G.keys()}
        if monic:
            z = nt
        else:
            y = len(lengths)
            z = max(nt, (sum(lengths) + y - 3) // (y - 1))

        lc_ev = ZZ.one
        if pseudomonic:
            lc = G[max(G.keys())][0][1]
            for el in lc:
                lc_ev *= pow(t[el[0]], el[1], p) % p

        for i in range(z):
            A_ev = []
            B_ev = []
            for pol in A_flat:
                if pol == {}:
                    A_ev.append(0)
                else:
                    h = 0
                    for j, k in pol.items():
                        # this i+1 is to prevent
                        # the first evaluation point from always being (1,..,1),
                        # we start the powers of the tuple from 1 and not from 0
                        h = (h + pow(j, (i+1), p) *k) %p
                    A_ev.append(h)

            for pol in B_flat:
                if pol == {}:
                    B_ev.append(0)
                else:
                    h = 0
                    for j, k in pol.items():
                        h = (h+ pow(j, (i+1), p) *k) %p
                    B_ev.append(h)

            G_ev = gf_gcd(A_ev[::-1], B_ev[::-1], p, ZZ)
            G_ev.reverse()
            # checking that the skeleton wasn't wrongly determined
            max_deg = max(G)
            if len(G_ev) != max_deg + 1:
                return None
            for j, el in enumerate(G_ev):
                if el != 0 and j not in G.keys():
                    return None
            for deg in G.keys():
                if pseudomonic:
                    eval_points[deg].append((G_ev[deg] * pow(lc_ev, i, p)) % p)
                else:
                    eval_points[deg].append(G_ev[deg])

        C = {}
        if monic:
            for i, (deg, monoms) in enumerate(G.items()):
                v = lag_basis(all_vand_basis[i], p)
                c = vandermonde_interp(v, eval_points[deg][:len(v)], p, trans=True)
                for j, mon in enumerate(monoms):

                    C[mon[0]] = (c[j] * pow(all_vand_basis[i][j], -1, p) ) %p
            return C
        else:
            # building the matrix to determine scaling coeffs
            mat = []
            zeros = [0] * (z-1)
            b = []
            for j, evalp in enumerate(eval_points.values()):
                for r in range(len(mat)):
                    mat[r] = mat[r][:-(z-1)] + [0]*lengths[j] + mat[r][-(z-1):]
                row = [0] * sum(lengths[:j])
                for i in range(z):
                    new_row = row + [pow(el, i+1, p) for el in all_vand_basis[j]]
                    if i != 0:
                        zeros[i-1] = - evalp[i]
                        new_row += zeros
                        zeros[i-1] = 0
                    else:
                        new_row += zeros
                    mat.append(new_row)

                b += [[evalp[0]]] + [[el] for el in zeros]

                # the first matrix will always be underdetermined
                if j == 0:
                    continue

                A_m = DomainMatrix.from_list(mat, FF(p))
                B_m = DomainMatrix.from_list(b, FF(p))

                try:
                    xnum, xden = A_m.solve_den_rref(B_m)
                except DMNonInvertibleMatrixError:
                    continue

                xden_inv = ZZ.invert(ZZ(xden), p)
                already_sol = j
                sol = [(int(el) * xden_inv) % p for el in xnum.to_list_flat()]
                scal_coeffs = [1] + sol[-(z - 1):]

                for i, (deg, monoms) in enumerate(G.items()):
                    if i <= already_sol:
                        for j, mon in enumerate(monoms):
                            C[mon[0]] = sol[sum(lengths[:i]) + j]

                    else:
                        v = lag_basis(all_vand_basis[i], p)
                        scaled_evalp = []
                        for x, y in zip(eval_points[deg][:len(v)], scal_coeffs[:len(v)]):
                            scaled_evalp.append((x * y) % p)

                        c = vandermonde_interp(v, scaled_evalp, p, trans=True)
                        for j, mon in enumerate(monoms):
                            C[mon[0]] = (c[j] * ZZ.invert(all_vand_basis[i][j], p) ) %p
                return C
    return None
