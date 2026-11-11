"""Heuristic polynomial GCD algorithm (HEUGCD)."""

from __future__ import annotations

from typing import TYPE_CHECKING
from sympy.polys.domains.integerring import ZZ
from sympy.polys.orderings import MonomialOrder, lex
from sympy.polys.sparsetools import (
    _smp_imul_ground,
    _smp_iprimitive,
    _smp_iquo_ground,
    _smp_isub,
    smp_LC,
    smp_content,
    smp_div_list,
    smp_evaluate,
    smp_neg,
    smp_quo_ground,
    smp_subs_drop,
    smp_trunc_ground,
)

from .polyerrors import HeuristicGCDFailed

if TYPE_CHECKING:
    from sympy.polys.sparsetools import smp
    from sympy.external.gmpy import MPZ
    from sympy.polys.rings import PolyElement

HEU_GCD_MAX = 6


def heugcd(
    f: PolyElement[MPZ], g: PolyElement[MPZ]
) -> tuple[PolyElement[MPZ], PolyElement[MPZ], PolyElement[MPZ]]:
    """
    Heuristic polynomial GCD in ``Z[X]``.

    Given univariate polynomials ``f`` and ``g`` in ``Z[X]``, returns
    their GCD and cofactors, i.e. polynomials ``h``, ``cff`` and ``cfg``
    such that::

          h = gcd(f, g), cff = quo(f, h) and cfg = quo(g, h)

    The algorithm is purely heuristic which means it may fail to compute
    the GCD. This will be signaled by raising an exception. In this case
    you will need to switch to another GCD method.

    The algorithm computes the polynomial GCD by evaluating polynomials
    ``f`` and ``g`` at certain points and computing (fast) integer GCD
    of those evaluations. The polynomial GCD is recovered from the integer
    image by interpolation. The evaluation process reduces f and g variable
    by variable into a large integer. The final step is to verify if the
    interpolated polynomial is the correct GCD. This gives cofactors of
    the input polynomials as a side effect.

    Examples
    ========

    >>> from sympy.polys.heuristicgcd import heugcd
    >>> from sympy.polys import ring, ZZ

    >>> R, x,y, = ring("x,y", ZZ)

    >>> f = x**2 + 2*x*y + y**2
    >>> g = x**2 + x*y

    >>> h, cff, cfg = heugcd(f, g)
    >>> h, cff, cfg
    (x + y, x + y, x)

    >>> cff*h == f
    True
    >>> cfg*h == g
    True

    References
    ==========

    .. [1] [Liao95]_

    """
    ring = f.ring
    assert ring == g.ring and ring.domain.is_ZZ
    h, cff, cfg = smp_heugcd(f, g, ring.ngens, ring.order)
    return f.new(h), f.new(cff), f.new(cfg)


def smp_heugcd(
    f: smp[MPZ], g: smp[MPZ], n: int, order: MonomialOrder = lex
) -> tuple[smp[MPZ], smp[MPZ], smp[MPZ]]:

    fc = smp_content(f, n, ZZ)
    gc = smp_content(g, n, ZZ)
    gcd = ZZ.gcd(fc, gc)

    f = smp_quo_ground(f, gcd, n, ZZ)
    g = smp_quo_ground(g, gcd, n, ZZ)

    abs = ZZ.abs

    f_norm = max(abs(el) for el in f.values())
    g_norm = max(abs(el) for el in g.values())
    B = ZZ(2 * min(f_norm, g_norm) + 29)

    x = max(
        min(B, 99 * ZZ.sqrt(B)),
        2
        * min(
            f_norm // abs(smp_LC(f, n, ZZ, order)),
            g_norm // abs(smp_LC(g, n, ZZ, order)),
        )
        + 4,
    )

    h0: MPZ | smp[MPZ]
    cff0: MPZ | smp[MPZ]
    cfg0: MPZ | smp[MPZ]

    for _ in range(0, HEU_GCD_MAX):
        if n == 1:
            ff = smp_evaluate(f, {0: x}, n, ZZ)
            gg = smp_evaluate(g, {0: x}, n, ZZ)
            if ff and gg:
                h0, cff0, cfg0 = ZZ.cofactors(ff, gg)
            else:
                x = 73794 * x * ZZ.sqrt(ZZ.sqrt(x)) // 27011
                continue
        else:
            ff = smp_subs_drop(f, {0: x}, n, ZZ)
            gg = smp_subs_drop(g, {0: x}, n, ZZ)
            if ff and gg:
                h0, cff0, cfg0 = smp_heugcd(ff, gg, n - 1, order)
            else:
                x = 73794 * x * ZZ.sqrt(ZZ.sqrt(x)) // 27011
                continue

        h = _gcd_interpolate(h0, x, n, order)
        _smp_iprimitive(h, n, ZZ)

        [cff_], r = smp_div_list(f, [h], n, ZZ, order)

        if not r:
            [cfg_], r = smp_div_list(g, [h], n, ZZ, order)
            if not r:
                _smp_imul_ground(h, gcd, n, ZZ)
                return h, cff_, cfg_

        cff = _gcd_interpolate(cff0, x, n, order)

        [h], r = smp_div_list(f, [cff], n, ZZ, order)

        if not r:
            [cfg_], r = smp_div_list(g, [h], n, ZZ, order)

            if not r:
                _smp_imul_ground(h, gcd, n, ZZ)
                return h, cff, cfg_

        cfg = _gcd_interpolate(cfg0, x, n, order)

        [h], r = smp_div_list(g, [cfg], n, ZZ, order)
        if not r:
            [cff_], r = smp_div_list(f, [h], n, ZZ, order)
            if not r:
                _smp_imul_ground(h, gcd, n, ZZ)
                return h, cff_, cfg

        x = 73794 * x * ZZ.sqrt(ZZ.sqrt(x)) // 27011
    raise HeuristicGCDFailed("no luck")


def _gcd_interpolate(
    h: MPZ | smp[MPZ], x: MPZ, n: int, order: MonomialOrder
) -> smp[MPZ]:
    """Interpolate polynomial GCD from integer GCD."""
    f: smp[MPZ] = {}
    i = 0

    # TODO: don't expose poly repr implementation details
    if isinstance(h, dict):  # h is a polynomial
        while h:
            g = smp_trunc_ground(h, x, n - 1, ZZ)
            _smp_isub(h, g, n - 1, ZZ)
            _smp_iquo_ground(h, x, n - 1, ZZ)

            # f += X**i*g
            if g:
                for monom, coeff in g.items():
                    f[(i,) + monom] = coeff
            i += 1

    else:  # h is a scalar
        while h:
            c = h % x
            if c > x // 2:
                c -= x
            h = (h - c) // x

            # f += X**i*c
            if c:
                f[(i,)] = c
            i += 1

    if smp_LC(f, n, ZZ, order) < 0:
        return smp_neg(f, n, ZZ)
    else:
        return f
