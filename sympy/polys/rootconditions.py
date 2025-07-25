from __future__ import annotations
from sympy.polys.domains import EXRAW
from sympy.polys.densebasic import dup_convert
from sympy.polys.densetools import dup_clear_denoms
from sympy.external.gmpy import MPQ
from sympy.core.mul import prod

def routh_hurwitz(p: Poly):
    """
    Return the conditions for the polynomial to have all roots with negative
    real part.

    """
    coeffs = p.rep.to_list()
    domain = p.domain

    conds = routh_hurwitz_dom(coeffs, domain)

    # return And(*[domain.to_sympy(c) > 0 for c in conds])
    return [domain.to_sympy(c) for c in conds]


def routh_hurwitz_dom(p: list[Er], K: Domain[Er]) -> list[Er]:

    if K.is_QQ:
        return _routh_hurwitz_qq(p, K)

    elif K.is_ZZ or K.is_RR:
        pq = dup_convert(p, K, QQ)
        conds = _routh_hurwitz_qq(pq, QQ)
        return dup_convert(conds, QQ, K)

    elif K.is_PolynomialRing:
        return _routh_hurwitz_poly(p, K)

    elif K.is_FractionField:
        _, pp = dup_clear_denoms(p, K, convert=True)
        conds = _routh_hurwitz_poly(pp, K)
        return dup_convert(conds, K.get_ring(), K)

    else:
        raise NotImplementedError


def _routh_hurwitz_qq(p: list[MPQ], K: RationalField) -> list[MPQ]:
    if len(p) < 2:
        return [K.one]
    elif p[0] * p[1] <= 0:
        return [-K.one]
    elif len(p) == 2:
        return [K.one]

    qs = p.copy()
    for i in range(1, len(p), 2):
        qs[i - 1] -= p[i] * p[0] / p[1]
    qs = qs[1:]

    return _routh_hurwitz_qq(qs, K)


def _routh_hurwitz_poly(p: list[PolyElement[Er]], K: PolynomialRing[Er]):
    return _rec_routh_hurwitz_poly(p, [], K)


def _rec_routh_hurwitz_poly(p: list[PolyElement[Er]],
                            previous_cond: list[PolyElement[Er]],
                            K: PolynomialRing[Er]):
    if len(p) < 2:
        return [K.one]

    if len(p) == 2:
        return [_clear_cond(p[0] * p[1], previous_cond, K)]

    qs = [p[1] * qi for qi in p]
    for i in range(1, len(p), 2):
        qs[i - 1] -= p[i] * p[0]
    qs = qs[1:]

    cond = p[0] * p[1]
    cond = _clear_cond(cond, previous_cond, K)

    return [cond] + _rec_routh_hurwitz_poly(qs, previous_cond, K)

def _clear_cond(cond: PolyElement[Er], previous_cond, K):
    # Divide out factors known to be positive
    for c in previous_cond:
        cond_quo, r = K.div(cond, c)
        while not r:
            cond = cond_quo
            cond_quo, r = K.div(cond, c)

    # Remove factors of even degree and reduce odd degree factors
    _, facs_m = cond.sqf_list()
    cond = prod([fac for fac, m in facs_m if m % 2 == 1])

    previous_cond += [cond]

    return cond

# TODO Implement conditions for discrete time systems
# Possuble ways are:
# - Map z = (1+s)/(1-s) and use routh hurwitz
# - Use Jury's test
