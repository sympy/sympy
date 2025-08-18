from __future__ import annotations

from sympy.polys.densebasic import dup
from sympy.polys.domains.domain import Er, Domain
from sympy.external.gmpy import MPQ

from sympy.polys.domains import EXRAW, QQ, PolynomialRing, RationalField
from sympy.polys.rings import PolyElement
from sympy.polys.densebasic import dup_convert
from sympy.polys.densetools import dup_clear_denoms
from sympy.core.mul import prod
from sympy.core.exprtools import factor_terms
from sympy.core.expr import Expr
from sympy.simplify.simplify import signsimp
from sympy.core.mul import Mul

def dup_routh_hurwitz(f: dup[Er], K: Domain[Er]) -> list[Er]:
    """
    Computes the Routh Hurwitz criteria of ``f``.

    The conditions, in general, are represented by a list of multivariate
    polynomials which must be positive to ensure stability.

    Note: This method assumes that the leading coefficient is non-zero.
    In the opposite case, additional verification is required.

    Depending on the domain, a different approach is used.
    In non-numeric cases, the algorithm is modified to avoid divisions.

    References
    ==========

    .. [1] G. Meinsma: Elementary proof of the Routh-Hurwitz test.
           Systems & Control Letters, Volume 25, Issue 4, 1995, Pages 237-242,
           https://courses.washington.edu/mengr471/resources/Routh_Hurwitz_Proof.pdf

    """
    if K.is_QQ:
        return _dup_routh_hurwitz_qq(f, K) # type: ignore

    elif K.is_ZZ or K.is_RR:
        pq = dup_convert(f, K, QQ)
        conds = _dup_routh_hurwitz_qq(pq, QQ)
        return dup_convert(conds, QQ, K)

    elif K.is_PolynomialRing:
        return _dup_routh_hurwitz_poly(f, K) # type: ignore

    elif K.is_FractionField:
        _, pp = dup_clear_denoms(f, K, convert=True)
        conds = _dup_routh_hurwitz_poly(pp, K) # type: ignore
        return dup_convert(conds, K.get_ring(), K)

    else:
        return _dup_routh_hurwitz_exraw(dup_convert(f, K, EXRAW)) # type: ignore


def _dup_routh_hurwitz_qq(p: list[MPQ], K: RationalField) -> list[MPQ]:
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

    return _dup_routh_hurwitz_qq(qs, K)


def _dup_routh_hurwitz_poly(p: list[PolyElement[Er]], K: PolynomialRing[Er]):
    return _rec_dup_routh_hurwitz_poly(p, [], K)


def _rec_dup_routh_hurwitz_poly(p: list[PolyElement[Er]],
                            previous_cond: list[PolyElement[Er]],
                            K: PolynomialRing[Er]):
    if len(p) < 2:
        return [K.one]

    if len(p) == 2:
        return previous_cond + [_clear_cond_poly(p[0] * p[1], previous_cond, K)]

    qs = [p[1] * qi for qi in p]
    for i in range(1, len(p), 2):
        qs[i - 1] -= p[i] * p[0]
    qs = qs[1:]

    cond = p[0] * p[1]

    # This check stops the recursion immediately and avoid a possible infinite
    # while loop in _clear_cond_poly.
    if cond.is_zero:
        return [K(-1)]

    cond = _clear_cond_poly(cond, previous_cond, K)

    # We avoid adding 1 or -1 to the previous conditions, since it will lead to
    # an infinite loop in the next _clear_cond_poly.
    if cond == K(-1):
        return [K(-1)]

    if not cond.is_one:
        previous_cond.append(cond)

    return _rec_dup_routh_hurwitz_poly(qs, previous_cond, K)


def _clear_cond_poly(cond: PolyElement[Er], previous_cond, K):
    # Divide out factors known to be positive from previous conditions.

    # There are not controls in that functions on cond and previous_cond,
    # we assume that at this point, there are no zeroes, ones and negative ones
    # in previous_cond.
    # If there are, there will be an infinite while loop.
    for c in previous_cond:
        cond_quo, r = K.div(cond, c)
        while not r:
            cond = cond_quo
            cond_quo, r = K.div(cond, c)

    # Remove factors of even degree and reduce odd degree factors
    common_factor, facs_m = cond.sqf_list()

    cond = prod([fac for fac, m in facs_m if m % 2 == 1], K.one)
    if common_factor < 0:
        cond = -K.one* cond

    return cond


def _dup_routh_hurwitz_exraw(p: list[Expr]) -> list[Expr]:
    if all(c.is_Number for c in p):
        return _dup_routh_hurwitz_exraw_div(p)

    return _dup_routh_hurwitz_exraw_no_div(p)


def _dup_routh_hurwitz_exraw_div(p: list[Expr]) -> list[Expr]:
    """
    Stability check with divisions, used for numeric cases.

    """
    if len(p) < 2:
        return [EXRAW.one]

    if (p[0]*p[1]).is_nonpositive:
        return [-EXRAW.one]

    if len(p) == 2:
        return [p[0]*p[1]]

    qs = p.copy()
    for i in range(1, len(p), 2):
        qs[i - 1] -= p[i] * p[0] / p[1]
    qs = qs[1:]

    return [p[0] * p[1]] + _dup_routh_hurwitz_exraw_div(qs)


def _dup_routh_hurwitz_exraw_no_div(p: list[Expr]) -> list[Expr]:
    """
    Stability check without divisions, used for EXRAW.

    """
    return _rec_dup_routh_hurwitz_exraw_no_div(p, [])


def _rec_dup_routh_hurwitz_exraw_no_div(p: list[Expr],
                                previous_cond: list[tuple[int, set[Expr]]]) -> list[Expr]:
    if len(p) < 2:
        return [EXRAW.one]

    if len(p) == 2:
        return [_clear_cond_exraw(p[0]*p[1], previous_cond)[0]]

    qs = [p[1] * qi for qi in p]
    for i in range(1, len(p), 2):
        qs[i - 1] -= p[i] * p[0]
    qs = qs[1:]

    cond, previous_cond = _clear_cond_exraw(p[0] * p[1], previous_cond)
    return [cond] + _rec_dup_routh_hurwitz_exraw_no_div(qs, previous_cond)


def _clear_cond_exraw(cond: Expr, previous_cond: list[tuple[int, set]]) -> tuple[Expr, list[tuple[int, set[Expr]]]]:
    """
    Clear the condition by removing even powers and simplifying odd powers.
    Also, remove factors that are already present in previous conditions.

    This function returns a simplified product of factors and the
    `previous_cond` list with the new condition.

    """
    factors_repr: tuple[int, set[Expr]] = _build_simplified_factors(cond, previous_cond)

    first_factor_iteration = factors_repr

    # Second iteration to simplify the factors further using `factor_terms`.
    # To reach a better simplification, it could be better to use `factor`, but
    # it could be slower.
    # Because the point of the algorithm for EXRAW is to be fast, we prefer to
    # use `factor_terms`.
    new_cond = Mul(factors_repr[0], *factors_repr[1])
    factors_repr = _build_simplified_factors(factor_terms(new_cond), previous_cond)

    # It's important to keep the first factor iteration, because it follows the
    # pattern of the unsimplified coefficients, which appears in the next
    # iteration of the Routh-Hurwitz algorithm.
    # Adding the last factor form, probably lead to a worse simplification.
    previous_cond = previous_cond + [first_factor_iteration]

    return Mul(factors_repr[0], *factors_repr[1]), previous_cond


def _build_simplified_factors(cond: Expr,
                              previous_cond: list[tuple[int, set[Expr]]]) -> tuple[int, set[Expr]]:
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
    powers_dict = dict(cond.as_powers_dict())

    for factor in powers_dict:
        if powers_dict[factor] % 2 != 0:
                factors.add(factor)

    factors_rep = (sign, factors)

    # Remove factors that are already present in previous conditions.
    for prev_factors in previous_cond:
        cond_quo = _div(factors_rep, prev_factors)
        if cond_quo is not None:
            factors_rep = cond_quo

    return factors_rep


def _div(n: tuple[int, set[Expr]],
        d: tuple[int, set[Expr]]) -> tuple[int, set[Expr]] | None:
    """Divide n by d if possible else return None."""
    nsign, nfactors = n
    dsign, dfactors = d
    if nfactors.issuperset(dfactors):
        return nsign*dsign, nfactors - dfactors
    else:
        return None

# TODO Implement conditions for discrete time systems
# Possible ways are:
# - Map z = (1+s)/(1-s) and use routh-hurwitz
# - Use Jury's test
