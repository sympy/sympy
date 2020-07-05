from sympy.core.relational import Eq, is_eq
from sympy.core.logic import fuzzy_and, fuzzy_bool
from sympy.logic.boolalg import And
from sympy.multipledispatch import dispatch
from sympy.sets.sets import tfn, ProductSet, Interval, FiniteSet


@dispatch(Interval, FiniteSet)
def _eval_is_eq(lhs, rhs): # noqa: F811
    return False


@dispatch(FiniteSet, Interval)
def _eval_is_eq(lhs, rhs): # noqa: F811
    return False


@dispatch(Interval, Interval)
def _eval_is_eq(lhs, rhs): # noqa: F811
    return And(Eq(lhs.left, rhs.left),
               Eq(lhs.right, rhs.right),
               lhs.left_open == rhs.left_open,
               lhs.right_open == rhs.right_open)


@dispatch(FiniteSet, Interval)
def _eval_is_eq(lhs, rhs): # noqa: F811
    return False


@dispatch(FiniteSet, FiniteSet)
def _eval_is_eq(lhs, rhs): # noqa: F811
    def all_in_both():
        s_set = set(lhs.args)
        o_set = set(rhs.args)
        yield fuzzy_and(lhs._contains(e) for e in o_set - s_set)
        yield fuzzy_and(rhs._contains(e) for e in s_set - o_set)

    return tfn[fuzzy_and(all_in_both())]


@dispatch(ProductSet, ProductSet)
def _eval_is_eq(lhs, rhs): # noqa: F811
    if len(lhs.sets) != len(rhs.sets):
        return False

    eqs = (is_eq(x, y) for x, y in zip(lhs.sets, rhs.sets))
    return tfn[fuzzy_and(map(fuzzy_bool, eqs))]
