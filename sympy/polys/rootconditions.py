from __future__ import annotations


from typing import TYPE_CHECKING, TypeAlias, Union
from sympy.core.relational import StrictGreaterThan
from sympy.logic.boolalg import Boolean

from sympy.polys.densebasic import dup, dmp, dmp_tup, _T
from sympy.polys.domains.domain import Er, Domain
from sympy.external.gmpy import MPQ

from sympy.polys.domains import EXRAW, QQ, PolynomialRing, RationalField
from sympy.polys.rings import PolyElement
from sympy.polys.densebasic import dup_convert
from sympy.polys.densetools import dup_clear_denoms
from sympy.core.mul import prod

from sympy.core.exprtools import factor_terms
from sympy.simplify.simplify import signsimp
from sympy.core.mul import Mul

conditions: TypeAlias = "list[Union[StrictGreaterThan, Boolean]]"

# The _dup and _dmp functions do not do anything but are needed so that a type
# checker can understand the conversion between the two types.
#
# A dup is a list of domain elements. A dmp is a list of lists of domain
# elements of arbitrary depth.

if TYPE_CHECKING:
    def _dup(p: dmp[_T], /) -> dup[_T]: ...
    def _dmp(p: dup[_T], /) -> dmp[_T]: ...
    def _dmp_tup(p: tuple[_T, ...], /) -> dmp_tup[_T]: ...
    def _idup(ps: tuple[dmp[_T], ...], /) -> tuple[dup[_T], ...]: ...
    def _idmp(ps: tuple[dup[_T], ...], /) -> tuple[dmp[_T], ...]: ...
else:

    def _dup(p, /):
        return p

    def _dmp(p, /):
        return p

    def _dmp_tup(p, /):
        return p

    def _idup(ps, /):
        return ps

    def _idmp(ps, /):
        return ps

def dup_routh_hurwitz(f: dup[Er], K: Domain[Er]) -> list[conditions]:
    """
    Return the conditions for the polynomial to have all roots with negative
    real part.

    Note: This method assumes that the leading coefficient is non-zero.
    In the opposite case, additional verification is required.

    Algorithm from:
    https://courses.washington.edu/mengr471/resources/Routh_Hurwitz_Proof.pdf
    In non-numeric cases, the algorithm is modified to avoid divisions.

    """
    conds = dup_routh_hurwitz_dom(f, K)

    # return [K.to_sympy(c) for c in conds]
    return [K.to_sympy(c) > 0 for c in conds]


def dup_routh_hurwitz_dom(p: list[Er], K: Domain[Er]) -> list[Er]:
    """
    Calculate the Routh-Hurwitz conditions for the polynomial with a different
    algorithm depending on the domain.

    """

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
        # For everything else, we use the EXRAW domain
        return _routh_hurwitz_exraw(p)


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
        return previous_cond + [_clear_cond_poly(p[0] * p[1], previous_cond, K)]

    qs = [p[1] * qi for qi in p]
    for i in range(1, len(p), 2):
        qs[i - 1] -= p[i] * p[0]
    qs = qs[1:]

    cond = p[0] * p[1]

    # This check is really important. It stops the recursion immediately and
    # avoid a possible infinite while loop in _clear_cond_poly.
    if cond.is_zero:
        return [K.zero]

    cond = _clear_cond_poly(cond, previous_cond, K)

    # These checks are also important, because adding 1 or -1 to the previous
    # conditions will lead to an infinite loop in the next _clear_cond_poly.
    if cond == K(-1):
        return [K(-1)]

    if not cond.is_one:
        previous_cond.append(cond)

    return _rec_routh_hurwitz_poly(qs, previous_cond, K)

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

    cond = K(prod([fac for fac, m in facs_m if m % 2 == 1]))
    if common_factor < 0:
        cond = cond * K(-1)

    return cond


def _routh_hurwitz_exraw(p: list[Er]) -> list[Er]:
    if all(c.is_number for c in p):
        return _calc_conditions_div(p)

    return _calc_conditions_no_div(p)

def _calc_conditions_div(p: list[Er]) -> list[Er]:
    """
    Stability check with divisions, used for numeric cases.

    """
    if len(p) < 2:
        return [EXRAW.one]

    if (p[0]*p[1]).is_nonpositive:
        return [EXRAW(-1)]

    if len(p) == 2:
        return [p[0]*p[1]]

    qs = p.copy()
    for i in range(1, len(p), 2):
        qs[i - 1] -= p[i] * p[0] / p[1]
    qs = qs[1:]

    return [p[0] * p[1]] + _calc_conditions_div(qs)

def _calc_conditions_no_div(p: list[Er]) -> list[Er]:
    """
    Stability check without divisions, used for EXRAW.

    """
    return _rec_calc_conditions_no_div(p, [])

def _rec_calc_conditions_no_div(p: list[Er],
                                previous_cond: list[set]) -> list[Er]:
    if len(p) < 2:
        return [EXRAW.one]

    if len(p) == 2:
        return [_clear_cond_exraw(p[0]*p[1], previous_cond)]

    qs = [p[1] * qi for qi in p]
    for i in range(1, len(p), 2):
        qs[i - 1] -= p[i] * p[0]
    qs = qs[1:]

    cond: Mul = _clear_cond_exraw(p[0] * p[1], previous_cond)
    return [cond] + _rec_calc_conditions_no_div(qs, previous_cond)

def _clear_cond_exraw(cond: Mul, previous_cond: list[set]) -> Mul:
    """
    Clear the condition by removing even powers and simplifying odd powers.
    Also, remove factors that are already present in previous conditions.

    This function returns a simplified product of factors and UPDATES the
    `previous_cond` list with the new condition.

    """
    factors: set[Er] = _build_simplified_factors(cond, previous_cond)

    first_factor_iteration = factors.copy()

    # Second iteration to simplify the factors further using factor_terms.
    # To reach a better simplification, it could be better to use `factor`, but
    # it could be slower.
    # Because the point of the algorithm for EXRAW is to be fast, we prefer to
    # use `factor_terms`.
    factors = _build_simplified_factors(factor_terms(Mul(*factors)),
                                        previous_cond)

    # It's important to keep the first factor iteration, because it follows the
    # pattern of the unsimplified coefficients, which appears in the next
    # iteration of the Routh-Hurwitz algorithm.
    # Adding the last factor form, probably lead to a worse simplification.
    previous_cond.append(first_factor_iteration)
    return Mul(*factors)

def _build_simplified_factors(cond: Mul,
                              previous_cond: list[set[Er]]) -> set[Er]:
    """
    Build a set of simplified factors from the expression removing even powers
    and simplifying odd powers.
    This function also removes factors that are already present in previous
    conditions.

    """
    factors: set[Er] = _extract_odd_factors(cond)
    _remove_previous_conditions(factors, previous_cond)
    return factors

def _extract_odd_factors(cond: Mul) -> set[Er]:
    """
    Extract factors with odd powers from a Mul expression.
    Even powers are ignored, as they do not contribute to the condition.
    Odd powers are simplified to the exponent of 1.

    """
    powers_dict = dict(cond.as_powers_dict()) # Dict of factors with their powers
    factors = set()
    n_minus_one = 0
    for factor in powers_dict:
        if powers_dict[factor] % 2 != 0:
            signsimp_factor = signsimp(factor, evaluate=False)
            if signsimp_factor.could_extract_minus_sign():
                factors.add(-signsimp_factor)
                n_minus_one += 1
            else:
                factors.add(signsimp_factor)

    if n_minus_one % 2 == 1:
        factors.add(-1)

    return factors

def _remove_previous_conditions(factors: set, previous_cond: list[set[Er]]):
    """
    Remove factors that are already present in previous conditions.

    """
    for prev_factors in previous_cond:
        if prev_factors.issubset(factors):
            factors -= prev_factors
        elif -1 in prev_factors:
            if (prev_factors - {-1}).issubset(factors):
                factors -= (prev_factors - {-1})
                factors.add(-1)

# TODO Implement conditions for discrete time systems
# Possible ways are:
# - Map z = (1+s)/(1-s) and use routh-hurwitz
# - Use Jury's test
