"""
Module for equation rearrangement with assumptions
"""

from sympy import Add, ask, log, Mul, Pow, Q, refine
from sympy.multipledispatch import dispatch
from sympy.relation.binrel import BinaryRelation
from sympy.relation.equality import Equal
from sympy.relation.inequality import InEqual


def rearrange(eqn, assumptions=True):
    """
    Rearrange the equation with given assumptions.

    This does not swap the equation or rewrite the equation with respect
    to certain variable.

    Note that this function is experimental. It may return incorrect
    result, or be changed at any time. Use with discretion.
    """
    new_eqn = _try_rearrange(eqn, assumptions)
    while eqn != new_eqn:
        eqn = new_eqn
        new_eqn = _try_rearrange(eqn, assumptions)
        lhs, rhs = refine(new_eqn.lhs, assumptions), refine(new_eqn.rhs, assumptions)
        lhs, rhs = lhs.simplify(), rhs.simplify()
        new_eqn = new_eqn.function(lhs, rhs)

    return new_eqn


def _try_rearrange(eqn, assumptions):
    try:
        result = _rearrange_equation(eqn.function, eqn.lhs, eqn.rhs, assumptions=assumptions)
    except NotImplementedError:
        result = eqn
    return result


@dispatch(BinaryRelation, Add, Add)
def _rearrange_equation(rel, lhs, rhs, assumptions=True): # noqa:F811
    lhs_coeff, lhs_terms = lhs.as_coeff_add()
    _, rhs_terms = rhs.as_coeff_add()
    commonterms = _commonelem_commut(lhs_terms, rhs_terms)
    # XXX Perhaps need to use SymPy's newly-introduced add() function
    commonterm = Add(*commonterms, evaluate=True)
    if lhs_coeff != 0:
        commonterm += lhs_coeff

    eqn = rel(lhs, rhs)
    if ask(Q.infinite(commonterm), assumptions) is False:
        eqn -= commonterm
    return eqn


@dispatch(Equal, Mul, Mul)
def _rearrange_equation(rel, lhs, rhs, assumptions=True): # noqa:F811
    lhs_coeff, lhs_terms = lhs.as_coeff_mul()
    _, rhs_terms = rhs.as_coeff_mul()
    commonterms = _commonelem_commut(lhs_terms, rhs_terms)
    # XXX Perhaps need to use SymPy's newly-introduced mul() function
    commonterm = Mul(*commonterms, evaluate=True)
    if lhs_coeff != 1:
        commonterm *= lhs_coeff

    eqn = rel(lhs, rhs)
    if ask(Q.nonzero(commonterm), assumptions):
        eqn /= commonterm
    return eqn


@dispatch(InEqual, Mul, Mul)
def _rearrange_equation(rel, lhs, rhs, assumptions=True): # noqa:F811
    lhs_coeff, lhs_terms = lhs.as_coeff_mul()
    _, rhs_terms = rhs.as_coeff_mul()
    commonterms = _commonelem_commut(lhs_terms, rhs_terms)
    # XXX Perhaps need to use SymPy's newly-introduced mul() function
    commonterm = Mul(*commonterms, evaluate=True)
    if lhs_coeff != 1:
        commonterm *= lhs_coeff

    eqn = rel(lhs, rhs)

    is_pos = ask(Q.positive(commonterm), assumptions)
    is_neg = ask(Q.negative(commonterm), assumptions)

    if (is_pos or is_neg):
        eqn /= commonterm
    return eqn


@dispatch(BinaryRelation, Pow, Pow)
def _rearrange_equation(rel, lhs, rhs, assumptions=True): # noqa:F811
    lhs_b, lhs_e = lhs.as_base_exp()
    rhs_b, rhs_e = rhs.as_base_exp()

    eqn = rel(lhs, rhs)
    if not ask(Q.eq(lhs_e, 1) & Q.eq(rhs_e, 1), assumptions):
        if ask(Q.positive(lhs_b) & Q.positive(rhs_b), assumptions):
            lhs_logb, rhs_logb = log(lhs_b), log(rhs_b)
            eqn = rel(lhs_logb*lhs_e, rhs_logb*rhs_e)

    return eqn


@dispatch(BinaryRelation, log, log)
def _rearrange_equation(rel, lhs, rhs, assumptions=True): # noqa:F811
    if ask(Q.complex(lhs) & Q.complex(rhs), assumptions):
        return rel(lhs.args[0], lhs.args[1])
    return rel(lhs, rhs)

def _commonelem_commut(iter1, iter2):
    # get the list of common elements in two iterables, with duplicates
    list1, list2 = list(iter1), list(iter2)

    nomore = set()
    result = []
    while list1:

        if not list2:
            break

        elem = list1.pop()

        if elem in nomore:
            continue

        try:
            index = list2.index(elem)
        except ValueError:
            nomore.add(elem)
            continue

        list2.pop(index)
        result.append(elem)
    return result
