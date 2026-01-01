from __future__ import annotations

from sympy.polys.densebasic import dup
from sympy.polys.domains.domain import Er, Domain
from sympy.core.numbers import I

from sympy.polys.densearith import dup_exquo_ground
from sympy.polys.domains import EXRAW, EX, QQ, ZZ
from sympy.polys.densebasic import dup_convert
from sympy.polys.densetools import dup_clear_denoms, dup_transform, dup_eval
from sympy.core.exprtools import factor_terms
from sympy.core.expr import Expr
from sympy.simplify.simplify import signsimp
from sympy.core.mul import Mul


def dup_routh_hurwitz(f: dup[Er], K: Domain[Er]) -> list[Er]:
    """
    Computes the Routh Hurwitz criteria of ``f``.

    Explanation
    ===========

    Returns a list of conditions ``[c1, c2, ...]`` such that the leading
    coefficient is nonzero and all roots of the polynomial have strictly
    negative real part if and only if ``ci > 0`` for all ``i``.

    References
    ==========

    .. [1] G. Meinsma: Elementary proof of the Routh-Hurwitz test.
           Systems & Control Letters, Volume 25, Issue 4, 1995, Pages 237-242,
           https://courses.washington.edu/mengr471/resources/Routh_Hurwitz_Proof.pdf

    """
    if not f:
        raise ValueError("zero polynomial")

    excluded_domains = [EX, EXRAW]
    if not K in excluded_domains and I in K:
        raise NotImplementedError(
            "Routh-Hurwitz is not implemented for complex domains"
        )
    if K.is_RR:
        pq = dup_convert(f, K, QQ)
        _, pz = dup_clear_denoms(pq, QQ, convert=True)
        conds: list = _dup_routh_hurwitz_fraction_free(pz, ZZ)
        return [K.convert_from(c, ZZ) for c in conds]

    if K.is_QQ:
        _, pz = dup_clear_denoms(f, K, convert=True)
        conds = _dup_routh_hurwitz_fraction_free(pz, ZZ)
        return [K.convert_from(c, ZZ) for c in conds]

    elif K.is_ZZ:
        return _dup_routh_hurwitz_fraction_free(f, K)

    elif K.is_PolynomialRing:
        return _dup_routh_hurwitz_fraction_free(f, K)

    elif K.is_FractionField:
        _, pp = dup_clear_denoms(f, K, convert=True)
        conds = _dup_routh_hurwitz_fraction_free(pp, K)
        R = K.get_ring()
        return [K.convert_from(c, R) for c in conds]

    else:
        pe = dup_convert(f, K, EXRAW)
        conds = _dup_routh_hurwitz_exraw(pe)
        return [K.convert_from(c, EXRAW) for c in conds]


def _dup_routh_hurwitz_fraction_free(p: dup[Er], K: Domain[Er]) -> list[Er]:
    if len(p) == 1:
        return []
    elif len(p) == 2:
        return [p[0] * p[1]]
    elif len(p) == 3:
        return [p[0] * p[1], p[0] * p[2]]

    LC = p[0]
    TC = p[-1]
    monic = K.is_one(LC)

    p1s = [p[1]]

    while len(p) > 3:
        qs = [p[1] * qi for qi in p[1:]]
        for i in range(1, len(qs) - 1, 2):
            qs[i] = qs[i] - p[i + 2] * p[0]
        p = qs

        p1 = p[1]
        if K.is_zero(p1):
            return [-K.one]

        if len(p1s) >= 2:
            p1 = K.exquo(p1, p1s[-2])
        if len(p1s) >= 3:
            p1 = K.exquo(p1, p1s[-3])
            p1 = K.exquo(p1, p1s[-3])
            p = dup_exquo_ground(p, p1s[-3], K)
            p = dup_exquo_ground(p, p1s[-3], K)

        p1s.append(p1)

    if not monic:
        for i in range(0, len(p1s), 2):
            p1s[i] *= LC

    p1s.append(LC * TC)

    return p1s


def _dup_routh_hurwitz_exraw(p: list[Expr]) -> list[Expr]:
    if all(c.is_Number for c in p):
        return _dup_routh_hurwitz_exraw_div(p)

    return _dup_routh_hurwitz_exraw_no_div(p)


def _dup_routh_hurwitz_exraw_div(p: list[Expr]) -> list[Expr]:
    """
    Stability check with divisions, used for numeric cases.

    """
    if len(p) < 2:
        return []

    if (p[0] * p[1]).is_nonpositive:
        return [-EXRAW.one]

    qs = p.copy()
    for i in range(1, len(p), 2):
        qs[i - 1] -= p[i] * p[0] / p[1]
    qs = qs[1:]

    return [p[0] * p[1]] + _dup_routh_hurwitz_exraw_div(qs)


def _dup_routh_hurwitz_exraw_no_div(p: list[Expr]) -> list[Expr]:
    """
    Stability check without divisions, used for EXRAW.

    """

    return _rec_dup_routh_hurwitz_exraw_no_div(p, [], set())


def _rec_dup_routh_hurwitz_exraw_no_div(
    p: list[Expr],
    previous_cond: list[tuple[int, set[Expr]]],
    nonzeros: set[Expr],
) -> list[Expr]:
    if len(p) < 2:
        return [z**2 for z in nonzeros]

    qs = [p[1] * qi for qi in p]
    for i in range(1, len(p), 2):
        qs[i - 1] -= p[i] * p[0]
    qs = qs[1:]

    cond, previous_cond, nonzeros = _clear_cond_exraw(
        p[0] * p[1], previous_cond, nonzeros
    )
    return [cond] + _rec_dup_routh_hurwitz_exraw_no_div(
        qs, previous_cond, nonzeros
    )


def _clear_cond_exraw(
    cond: Expr, previous_cond: list[tuple[int, set[Expr]]], nonzeros: set[Expr]
) -> tuple[Expr, list[tuple[int, set[Expr]]], set[Expr]]:
    """
    Clear the condition by removing even powers and simplifying odd powers.
    Also, remove factors that are already present in previous conditions.

    This function returns a simplified product of factors and the
    `previous_cond` list with the new condition.

    """
    factors_repr, new_nonzeros1 = _build_simplified_factors(cond, previous_cond)

    nonzeros_updated = nonzeros | new_nonzeros1
    first_factor_iteration = factors_repr

    # Second iteration to simplify the factors further using `factor_terms`.
    # To reach a better simplification, it could be better to use `factor`, but
    # it could be slower.
    # Because the point of the algorithm for EXRAW is to be fast, we prefer to
    # use `factor_terms`.
    new_cond = Mul(factors_repr[0], *factors_repr[1])
    factors_repr, new_nonzeros2 = _build_simplified_factors(
        factor_terms(new_cond), previous_cond
    )

    nonzeros_updated.update(new_nonzeros2)

    # It's important to keep the first factor iteration, because it follows the
    # pattern of the unsimplified coefficients, which appears in the next
    # iteration of the Routh-Hurwitz algorithm.
    # Adding the last factor form, probably lead to a worse simplification.
    previous_cond = previous_cond + [first_factor_iteration]

    clear_cond = Mul(factors_repr[0], *factors_repr[1])
    return clear_cond, previous_cond, nonzeros_updated


def _build_simplified_factors(
    cond: Expr, previous_cond: list[tuple[int, set[Expr]]]
) -> tuple[tuple[int, set[Expr]], set[Expr]]:
    """
    Build a set of simplified factors from the expression rem oving even powers
    and simplifying odd powers.
    This function also removes factors that are already present in previous
    conditions.

    """
    cond, sign = signsimp(cond, evaluate=False), 1

    if cond.could_extract_minus_sign():
        cond, sign = -cond, -1

    # Build a set of factors without even powers
    factors = set()
    nonzeros = set()

    for b, e in cond.as_powers_dict().items():
        if not e.is_Integer:
            factors.add(b**e)
        elif e % 2 != 0:
            factors.add(b)
        else:
            nonzeros.add(b)

    # Remove factors that are already present in previous conditions.
    for sign_prev, factors_prev in previous_cond:
        if factors.issuperset(factors_prev):
            sign *= sign_prev
            factors -= factors_prev

        nonzeros -= factors_prev

    return (sign, factors), nonzeros


def dup_schur_conditions(f: dup[Er], K: Domain[Er]) -> list[Er]:
    """
    Computes the schur stability conditions of ``f``.

    Explanation
    ===========

    Returns expressions ``[e1, e2, ...]`` such that all roots of the
    polynomial lie inside the unit circle if and only if ``ei > 0``
    for all ``i``.

    References
    ==========

    .. [1] https://faculty.washington.edu/chx/teaching/me547/2_1_stability.pdf#:~:text=2.6%20Routh,plane%20Real

    """
    if not f:
        raise ValueError("zero polynomial")

    excluded_domains = [EX, EXRAW]
    if not K in excluded_domains and I in K:
        raise NotImplementedError(
            "Schur conditions is not implemented for complex domains"
        )

    if K.is_Exact:
        return _dup_schur_conditions(f, K)
    else:
        K_exact = K.get_exact()
        pe = dup_convert(f, K, K_exact)

        conds = _dup_schur_conditions(pe, K_exact)
        return [K.convert_from(c, K_exact) for c in conds]


def _dup_schur_conditions(f: dup[Er], K: Domain[Er]) -> list[Er]:
    """
    Computes the Schur stability conditions, mapping the unit circle to the left
    half plane, transforming the problem into a Routh-Hurwitz stability problem.

    The transformation is done via ``z = (1 + s)/(1 - s)`` to the polynomial
    ``f``.

    """
    # Check if -1 is a root, since the transformation is not defined in this case
    if K.is_zero(dup_eval(f, -K.one, K)):
        return [-K.one]

    f_t = dup_transform(f, [K.one, K.one], [-K.one, K.one], K)
    return dup_routh_hurwitz(f_t, K)
