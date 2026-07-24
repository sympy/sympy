from __future__ import annotations
from typing import TYPE_CHECKING
from collections import defaultdict
from sympy.polys.domains.integerring import ZZ
from sympy.polys.galoistools import gf_add, gf_gcd, gf_mul, gf_quo, gf_from_dict
from sympy.ntheory.modular import crt
from sympy.polys.domains import PolynomialRing

if TYPE_CHECKING:
    from sympy.external.gmpy import MPZ
    from sympy.polys.densebasic import dup
    from sympy.polys.monomials import monom


def _gf_gcd(fp, gp, p):
    r"""
    Compute the GCD of two univariate polynomials in `\mathbb{Z}_p[x]`.

    Parameters
    ==========

    fp : PolyElement
        Univariate integer polynomial modulo `p`.
    gp : PolyElement
        Univariate integer polynomial modulo `p`.
    p : int
        Prime modulus.

    Returns
    =======

    gcd : PolyElement
        The monic GCD of `fp` and `gp` in `\mathbb{Z}_p[x]`.

    Examples
    ========

    >>> from sympy.polys.zippel import _gf_gcd
    >>> from sympy.polys import ring, ZZ

    >>> R, x = ring("x", ZZ)
    >>> p = 5

    >>> f = (x**2 - 1).trunc_ground(p)
    >>> g = (x**2 - 2*x + 1).trunc_ground(p)

    >>> _gf_gcd(f, g, p)
    x - 1

    """
    dom = fp.ring.domain
    f_list = fp.to_dense()
    g_list = gp.to_dense()

    gcd_list = gf_gcd(f_list, g_list, p, dom)

    return fp.ring.from_dense(gcd_list).trunc_ground(p)


def _trivial_gcd(f, g):
    """
    Compute the GCD of two polynomials in trivial cases, i.e. when one
    or both polynomials are zero, or when both are constant.

    Parameters
    ==========

    f : PolyElement
        A multivariate polynomial.
    g : PolyElement
        A multivariate polynomial.

    Returns
    =======

    A tuple (gcd, cff, cfg) with the gcd and cofactors,
    or None if the GCD cannot be computed trivially.

    Examples
    ========

    >>> from sympy.polys.zippel import _trivial_gcd
    >>> from sympy.polys import ring, ZZ

    >>> R, x, y = ring("x, y", ZZ)

    >>> _trivial_gcd(R.zero, R.zero)
    (0, 0, 0)

    >>> _trivial_gcd(x**2 + y, R.zero)
    (x**2 + y, 1, 0)

    >>> f = R.ground_new(15)
    >>> g = R.ground_new(10)
    >>> _trivial_gcd(f, g)
    (5, 3, 2)

    >>> _trivial_gcd(x + 1, y - 1) is None
    True

    """
    ring = f.ring

    if not (f or g):
        return ring.zero, ring.zero, ring.zero
    elif not f:
        if g.LC < ring.domain.zero:
            return -g, ring.zero, -ring.one
        else:
            return g, ring.zero, ring.one
    elif not g:
        if f.LC < ring.domain.zero:
            return -f, -ring.one, ring.zero
        else:
            return f, ring.one, ring.zero
    elif f.is_ground and g.is_ground:
        c = ring.domain.gcd(f.LC, g.LC)
        h = ring.ground_new(c)
        return h, f.quo_ground(c), g.quo_ground(c)
    return None


def _primitive(f, p):
    r"""
    Compute the content and the primitive part of a polynomial in
    `\mathbb{Z}_p[x_0, \ldots, x_{k-2}, y] \cong \mathbb{Z}_p[y][x_0, \ldots, x_{k-2}]`.

    Parameters
    ==========

    f : PolyElement
        integer polynomial in `\mathbb{Z}_p[x0, \ldots, x{k-2}, y]`
    p : Integer
        modulus of `f`

    Returns
    =======

    contf : PolyElement
        integer polynomial in `\mathbb{Z}_p[y]`, content of `f`
    ppf : PolyElement
        primitive part of `f`, i.e. `\frac{f}{contf}`

    Examples
    ========

    >>> from sympy.polys.zippel import _primitive
    >>> from sympy.polys import ring, ZZ

    >>> R, x, y = ring("x, y", ZZ)
    >>> p = 3

    >>> f = x**2*y**2 + x**2*y - y**2 - y
    >>> _primitive(f, p)
    (y**2 + y, x**2 - 1)

    >>> R, x, y, z = ring("x, y, z", ZZ)

    >>> f = x*y*z - y**2*z**2
    >>> _primitive(f, p)
    (z, x*y - y**2*z)

    """
    ring = f.ring
    dom = ring.domain
    k = ring.ngens

    coeffs = {}
    for monom, coeff in f.iterterms():
        if monom[:-1] not in coeffs:
            coeffs[monom[:-1]] = {}
        coeffs[monom[:-1]][monom[-1]] = coeff

    cont = []

    for mon, coeff in coeffs.items():
        coeff = gf_from_dict(coeff, p, dom)
        coeffs[mon] = coeff
        cont = gf_gcd(cont, coeff, p, dom)

    prim = {}
    for mon, coeff in coeffs.items():
        coeff = gf_quo(coeff, cont, p, dom)
        deg = len(coeff)
        for i, el in enumerate(coeff):
            if el != 0:
                prim[mon + (deg-1-i,)] = el

    yring = ring.clone(symbols=[ring.symbols[k-1]])
    contf = yring.from_dense(cont).trunc_ground(p)

    return contf, ring.from_dict(prim).trunc_ground(p)


def _LC(f):
    r"""
    Compute the leading coefficient of a multivariate polynomial
    `f \in K[x_0, \ldots, x_{k-2}, y] \cong K[y][x_0, \ldots, x_{k-2}]`.

    Parameters
    ==========

    f : PolyElement
        polynomial in `K[x_0, \ldots, x_{k-2}, y]`

    Returns
    =======

    lcf : PolyElement
        polynomial in `K[y]`, leading coefficient of `f`

    Examples
    ========

    >>> from sympy.polys.zippel import _LC
    >>> from sympy.polys import ring, ZZ

    >>> R, x, y = ring("x, y", ZZ)

    >>> f = x**2*y**2 + x**2*y - 1
    >>> _LC(f)
    y**2 + y

    >>> R, x, y, z = ring("x, y, z", ZZ)

    >>> f = x**2*y**2 + x**2*y - 1
    >>> _LC(f)
    1

    >>> f = x*y*z - y**2*z**2
    >>> _LC(f)
    z

    """
    ring = f.ring
    k = ring.ngens
    yring = ring.clone(symbols=ring.symbols[k-1])
    y = yring.gens[0]
    degf = _deg(f)

    lcf = yring.zero
    for monom, coeff in f.iterterms():
        if monom[:-1] == degf:
            lcf += coeff*y**monom[-1]
    return lcf


def _deg(f):
    r"""
    Compute the degree of a multivariate polynomial
    `f \in K[x_0, \ldots, x_{k-2}, y] \cong K[y][x_0, \ldots, x_{k-2}]`.

    Parameters
    ==========

    f : PolyElement
        polynomial in `K[x_0, \ldots, x_{k-2}, y]`

    Returns
    =======

    degf : Integer tuple
        degree of `f` in `x_0, \ldots, x_{k-2}`

    Examples
    ========

    >>> from sympy.polys.zippel import _deg
    >>> from sympy.polys import ring, ZZ

    >>> R, x, y = ring("x, y", ZZ)

    >>> f = x**2*y**2 + x**2*y - 1
    >>> _deg(f)
    (2,)

    >>> R, x, y, z = ring("x, y, z", ZZ)

    >>> f = x**2*y**2 + x**2*y - 1
    >>> _deg(f)
    (2, 2)

    >>> f = x*y*z - y**2*z**2
    >>> _deg(f)
    (1, 1)

    """
    k = f.ring.ngens
    degf = (0,) * (k-1)
    for monom in f.itermonoms():
        if monom[:-1] > degf:
            degf = monom[:-1]
    return degf


def _chinese_remainder_reconstruction_multivariate(hp, hq, p, q):
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

    hp : PolyElement
        multivariate integer polynomial with coefficients in `\mathbb{Z}_p`
    hq : PolyElement
        multivariate integer polynomial with coefficients in `\mathbb{Z}_q`
    p : Integer
        modulus of `h_p`, relatively prime to `q`
    q : Integer
        modulus of `h_q`, relatively prime to `p`

    Examples
    ========

    >>> from sympy.polys.zippel import _chinese_remainder_reconstruction_multivariate
    >>> from sympy.polys import ring, ZZ

    >>> R, x, y = ring("x, y", ZZ)
    >>> p = 3
    >>> q = 5

    >>> hp = x**3*y - x**2 - 1
    >>> hq = -x**3*y - 2*x*y**2 + 2

    >>> hpq = _chinese_remainder_reconstruction_multivariate(hp, hq, p, q)
    >>> hpq
    4*x**3*y + 5*x**2 + 3*x*y**2 + 2

    >>> hpq.trunc_ground(p) == hp
    True
    >>> hpq.trunc_ground(q) == hq
    True

    >>> R, x, y, z = ring("x, y, z", ZZ)
    >>> p = 6
    >>> q = 5

    >>> hp = 3*x**4 - y**3*z + z
    >>> hq = -2*x**4 + z

    >>> hpq = _chinese_remainder_reconstruction_multivariate(hp, hq, p, q)
    >>> hpq
    3*x**4 + 5*y**3*z + z

    >>> hpq.trunc_ground(p) == hp
    True
    >>> hpq.trunc_ground(q) == hq
    True

    """
    hpmonoms = set(hp.monoms())
    hqmonoms = set(hq.monoms())
    monoms = hpmonoms.intersection(hqmonoms)
    hpmonoms.difference_update(monoms)
    hqmonoms.difference_update(monoms)

    domain = hp.ring.domain
    zero = domain.zero

    hpq = hp.ring.zero

    if isinstance(hp.ring.domain, PolynomialRing):
        crt_ = _chinese_remainder_reconstruction_multivariate
    else:
        def crt_(cp, cq, p, q):
            return domain(crt([p, q], [cp, cq], symmetric=True)[0])

    for monom in monoms:
        hpq[monom] = crt_(hp[monom], hq[monom], p, q)
    for monom in hpmonoms:
        hpq[monom] = crt_(hp[monom], zero, p, q)
    for monom in hqmonoms:
        hpq[monom] = crt_(zero, hq[monom], p, q)

    return hpq


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
) -> tuple[dict[int, list[list[monom]]], list[list[MPZ]], bool, bool]:
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
            2: [[(2, 2, 0), (0, 2)]],
            1: [[(1, 5, 3), (0, 5), (1, 3)], [(1, 0, 4), (1, 4)]]
        }

        h = [[5], [7], [8]]

        monic = True
        pseudomonic = True

    """
    S_ = defaultdict(list)
    for mon, _ in G.items():
        S_[mon[0]].append([mon])
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
        for mon_list in mons:
            h.append([G[mon_list[0]]])
            for i, el in enumerate(mon_list[0][1:]):
                if el != 0:
                    mon_list.append((i, el))

    return S, h, monic, pseudomonic
