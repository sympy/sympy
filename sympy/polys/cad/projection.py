"""Projection operators for cylindrical algebraic decomposition.

The projection phase of CAD takes a set of polynomials in
`x_1, \\ldots, x_k` and produces a set of polynomials in
`x_1, \\ldots, x_{k-1}` such that, over every cell of a decomposition of
`\\mathbb{R}^{k-1}` which is sign-invariant for the projected set, the real
roots of the original polynomials in `x_k` are delineable: they are given by
finitely many continuous functions, which do not cross, and the polynomials
have constant signs between them.

Two operators are implemented:

* :func:`mccallum_projection` is the projection of McCallum [1]_ applied to a
  squarefree basis: coefficients, discriminants and pairwise resultants. It
  gives small projection sets but it is only guaranteed to work if no
  polynomial of the projection factor sets vanishes identically on a cell of
  positive dimension (the set is *well-oriented*). The lifting phase checks
  this and falls back to the second operator if the check fails.

* :func:`hong_projection` is the projection of Hong [2]_, an improvement of
  the original operator of Collins, based on the reducta of the polynomials
  and on their principal subresultant coefficients. It is always correct but
  produces many more polynomials.

References
==========

.. [1] S. McCallum, An improved projection operation for cylindrical
       algebraic decomposition, in: Quantifier Elimination and Cylindrical
       Algebraic Decomposition, Springer, 1998, pp. 242-268.
.. [2] H. Hong, An improvement of the projection operator in cylindrical
       algebraic decomposition, ISSAC 1990, pp. 261-264.
"""
from __future__ import annotations

from itertools import combinations

from sympy.polys.densebasic import dmp_strip, dmp_zero_p
from sympy.polys.domains import ZZ
from sympy.polys.euclidtools import dmp_psc
from sympy.polys.polyclasses import DMP
from sympy.polys.polyerrors import PolynomialError
from sympy.polys.polytools import Poly


def _to_polys(polys, gens):
    """Convert ``polys`` to primitive polynomials over ZZ in ``gens``."""
    result = []
    for p in polys:
        if isinstance(p, Poly):
            p = p.as_expr()
        p = Poly(p, *gens)
        if not p.domain.is_ZZ:
            if not (p.domain.is_QQ or p.domain.is_ZZ):
                raise PolynomialError(
                    "polynomial with rational coefficients expected, got %s" % p)
            _, p = p.clear_denoms()
            p = p.set_domain(ZZ)
        result.append(p)
    return result


def _level(f, gens):
    """Index of the last generator of ``gens`` that ``f`` depends on."""
    for k in range(len(gens), 0, -1):
        if gens[k - 1] in f.gens and f.degree(gens[k - 1]) > 0:
            return k
    return 0


def _factors(f):
    """Irreducible non-constant factors of ``f`` with positive leading
    coefficient."""
    factors = []
    for g, _ in f.factor_list()[1]:
        if g.is_ground:
            continue
        if g.LC() < 0:
            g = -g
        factors.append(g)
    return factors


def squarefree_basis(polys, gens):
    """Finest squarefree basis of ``polys``, split by level.

    Returns a list ``[B_1, ..., B_n]`` where ``n = len(gens)`` and ``B_k``
    contains the irreducible factors (over ZZ, with positive leading
    coefficient) of the polynomials that depend on ``gens[k-1]`` but not on
    any later generator. Constants are dropped.

    Examples
    ========

    >>> from sympy.abc import x, y
    >>> from sympy.polys.cad import squarefree_basis
    >>> squarefree_basis([(x**2 - 1)*(x + y)**2, 6], [x, y])
    [[Poly(x - 1, x, y, domain='ZZ'), Poly(x + 1, x, y, domain='ZZ')], [Poly(x + y, x, y, domain='ZZ')]]

    The polynomials in ``B_k`` are pairwise coprime and squarefree, so that
    the sign of each input polynomial is determined by the signs of the basis
    elements.
    """
    levels = [[] for _ in gens]
    for f in _to_polys(polys, gens):
        for g in _factors(f):
            k = _level(g, gens)
            if k and g not in levels[k - 1]:
                levels[k - 1].append(g)
    return levels


def _main_first(f, x):
    """Reorder the generators of ``f`` so that ``x`` comes first."""
    gens = [x] + [g for g in f.gens if g != x]
    return f.reorder(*gens)


def _from_dense(rep, f):
    """Polynomial in the generators of ``f`` after the first one from the
    dense representation ``rep``. Univariate ``f`` give domain elements."""
    rest = f.gens[1:]
    if not rest:
        return rep
    return Poly.new(DMP(rep, f.domain, len(rest) - 1), *rest)


def _coefficients(f):
    """Coefficients of ``f`` with respect to its first generator, from the
    leading one down, as polynomials in the other generators."""
    return [_from_dense(c, f) for c in f.rep.to_list()]


def _psc(f, g):
    """Principal subresultant coefficients of ``f`` and ``g`` with respect
    to their first generator."""
    lev = len(f.gens) - 1
    psc = dmp_psc(f.rep.to_list(), g.rep.to_list(), lev, f.domain)
    return [_from_dense(c, f) for c in psc]


def _reducta(f):
    """The nonzero reducta of ``f`` with respect to its first generator:
    ``f``, ``f`` minus its leading term, and so on."""
    rep = f.rep.to_list()
    lev = len(f.gens) - 1
    reducta = []
    while not dmp_zero_p(rep, lev):
        reducta.append(Poly.new(DMP(rep, f.domain, lev), *f.gens))
        rep = dmp_strip(rep[1:], lev)
    return reducta


def _is_constant(c):
    return not isinstance(c, Poly) or c.is_ground


def _is_zero(c):
    return c.is_zero if isinstance(c, Poly) else not c


def _collect(items):
    """Irreducible factors of the non-constant polynomials in ``items``."""
    result = []
    for c in items:
        if _is_constant(c):
            continue
        for g in _factors(c):
            if g not in result:
                result.append(g)
    return result


def mccallum_projection(polys, x):
    """McCallum's projection of ``polys`` with respect to ``x``.

    ``polys`` must be a squarefree basis (pairwise coprime squarefree
    polynomials) of positive degree in ``x``, see :func:`squarefree_basis`.
    The projection consists of the irreducible factors of the coefficients
    with respect to ``x``, of the discriminants and of the pairwise
    resultants.

    Only the coefficients from the leading one down to the first nonzero
    constant coefficient are used: they are all that is needed for the
    degree to be invariant on each cell, and if some coefficient is a nonzero
    constant the polynomial cannot vanish identically anywhere.

    Examples
    ========

    >>> from sympy import Poly
    >>> from sympy.abc import x, y
    >>> from sympy.polys.cad import mccallum_projection
    >>> mccallum_projection([Poly(x**2 + y**2 - 1, x, y)], y)
    [Poly(x - 1, x, domain='ZZ'), Poly(x + 1, x, domain='ZZ')]
    >>> mccallum_projection([Poly(x*y - 1, x, y), Poly(y - x, x, y)], y)
    [Poly(x, x, domain='ZZ'), Poly(x - 1, x, domain='ZZ'), Poly(x + 1, x, domain='ZZ')]
    """
    polys = [_main_first(f, x) for f in polys]
    items = []
    for f in polys:
        for c in _coefficients(f):
            if not _is_constant(c):
                items.append(c)
            elif not _is_zero(c):
                break
        if f.degree() >= 2:
            items.append(f.discriminant())
    for f, g in combinations(polys, 2):
        items.append(f.resultant(g))
    return _collect(items)


def hong_projection(polys, x):
    """Hong's projection of ``polys`` with respect to ``x``.

    ``polys`` must be a squarefree basis of positive degree in ``x``. For
    every reductum `r` of every polynomial the projection contains the
    leading coefficient of `r` and the principal subresultant coefficients of
    `r` and its derivative; for every pair `f, g` it contains the principal
    subresultant coefficients of every reductum of `f` with `g`.

    Examples
    ========

    >>> from sympy import Poly
    >>> from sympy.abc import x, y
    >>> from sympy.polys.cad import hong_projection
    >>> hong_projection([Poly(x**2 + y**2 - 1, x, y)], y)
    [Poly(x - 1, x, domain='ZZ'), Poly(x + 1, x, domain='ZZ')]
    >>> hong_projection([Poly(x*y - 1, x, y), Poly(y - x, x, y)], y)
    [Poly(x, x, domain='ZZ'), Poly(x - 1, x, domain='ZZ'), Poly(x + 1, x, domain='ZZ')]
    """
    polys = [_main_first(f, x) for f in polys]
    items = []
    for f in polys:
        for r in _reducta(f):
            items.append(_coefficients(r)[0])
            if r.degree() >= 2:
                items.extend(_psc(r, r.diff(r.gens[0])))
    for f, g in combinations(polys, 2):
        for r in _reducta(f):
            if r.degree() == 0:
                items.append(_coefficients(r)[0])
            else:
                items.extend(_psc(r, g))
    return _collect(items)


_PROJECTIONS = {
    'mccallum': mccallum_projection,
    'hong': hong_projection,
}


def projection_sets(polys, gens, method='mccallum'):
    """Projection factor sets of ``polys`` for the variables ``gens``.

    The last generator is projected first. Returns a list
    ``[P_1, ..., P_n]`` where ``P_k`` is a squarefree basis of polynomials in
    ``gens[:k]`` of positive degree in ``gens[k-1]``: ``P_n`` is the basis of
    the input and each ``P_k`` for ``k < n`` contains the projection of
    ``P_{k+1}`` together with the factors of the input (and of the
    projections of higher levels) that only depend on ``gens[:k]``.

    ``method`` is ``'mccallum'`` (default) or ``'hong'``.

    Examples
    ========

    >>> from sympy.abc import x, y
    >>> from sympy.polys.cad import projection_sets
    >>> projection_sets([x**2 + y**2 - 1, x - y], [x, y])
    [[Poly(x - 1, x, domain='ZZ'), Poly(x + 1, x, domain='ZZ'), Poly(2*x**2 - 1, x, domain='ZZ')], [Poly(x**2 + y**2 - 1, x, y, domain='ZZ'), Poly(x - y, x, y, domain='ZZ')]]
    """
    try:
        project = _PROJECTIONS[method]
    except KeyError:
        raise ValueError("unknown projection method %r" % (method,))

    gens = list(gens)
    if not gens:
        raise ValueError("at least one generator is needed")

    levels = squarefree_basis(polys, gens)
    n = len(gens)

    result = [None]*n
    for k in range(n, 0, -1):
        current = [Poly(f.as_expr(), *gens[:k]) for f in levels[k - 1]]
        result[k - 1] = current
        if k == 1:
            break
        for g in project(current, gens[k - 1]):
            g = Poly(g.as_expr(), *gens)
            j = _level(g, gens)
            if j and g not in levels[j - 1]:
                levels[j - 1].append(g)
    return result
