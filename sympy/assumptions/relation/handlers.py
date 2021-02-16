"""
Handlers for equalities and inequalities
"""
from sympy.assumptions import ask, Q
from sympy.combinatorics import Permutation
from sympy.core import Basic, Expr, Symbol, Tuple
from sympy.core.function import AppliedUndef
from sympy.core.logic import fuzzy_bool, fuzzy_and, fuzzy_not
from sympy.core.numbers import Number, Rational
from sympy.core.relational import is_le
from sympy.core.sympify import sympify
from sympy.functions import ceiling, floor, frac
from sympy.logic.boolalg import And, Boolean
from sympy.matrices import ImmutableDenseMatrix, MatrixExpr
from sympy.multipledispatch import MDNotImplementedError
from sympy.polys import ComplexRootOf
from sympy.sets.sets import tfn, ProductSet, Interval, FiniteSet, Set

from .equality import EqualityPredicate, UnequalityPredicate


### EqualityPredicate ###

# If type A is superclass of type B, dispatch as (A, B) if possible
# to avoid unambiguity.

@EqualityPredicate.register(Basic, Basic)
def _(lhs, rhs, assumptions):
    return None

@EqualityPredicate.register(Boolean, Boolean)
def _(lhs, rhs, assumptions):
    # Only boolean can be equal to boolean, except Symbol which can
    # represent arbitrary boolean object.
    if (lhs.is_Symbol or rhs.is_Symbol):
        raise MDNotImplementedError  # Delegate to (Expr, Expr)
    lhs_eval, rhs_eval = ask(lhs, assumptions), ask(rhs, assumptions)
    if lhs_eval is None or rhs_eval is None:
        return None
    return ask(lhs, assumptions) == ask(rhs, assumptions)

@EqualityPredicate.register(Basic, Boolean)
def _(lhs, rhs, assumptions):
    # Only boolean can be equal to boolean, except Symbol which can
    # represent arbitrary boolean object.
    if (lhs.is_Symbol or rhs.is_Symbol):
        raise MDNotImplementedError  # Delegate to (Expr, Expr)
    return False

@EqualityPredicate.register(Expr, Tuple)
@EqualityPredicate.register(Tuple, Number)
def _(lhs, rhs, assumptions):
    return False

@EqualityPredicate.register(AppliedUndef, Tuple)
@EqualityPredicate.register(Symbol, Tuple)
def _(lhs, rhs, assumptions):
    return None

@EqualityPredicate.register(Tuple, Tuple)
def _(lhs, rhs, assumptions):
    if len(lhs) != len(rhs):
        return False
    return fuzzy_and(fuzzy_bool(ask(Q.eq(s, o), assumptions)) for s, o in zip(lhs, rhs))

@EqualityPredicate.register(Expr, floor)
def _(lhs, rhs, assumptions):
    return Q.eq.eval((lhs, rhs.rewrite(ceiling)), assumptions) or \
        Q.eq.eval((lhs, rhs.rewrite(frac)), assumptions)

@EqualityPredicate.register(Basic, ceiling)
def _(lhs, rhs, assumptions):
    return Q.eq.eval((lhs, rhs.rewrite(floor)), assumptions) or \
        Q.eq.eval((lhs, rhs.rewrite(frac)), assumptions)

@EqualityPredicate.register(Basic, frac)
def _(lhs, rhs, assumptions):
    if (lhs == rhs.rewrite(floor)) or \
        (lhs == rhs.rewrite(ceiling)):
        return True
    # Check if other < 0
    if lhs.is_extended_negative:
        return False
    # Check if other >= 1
    res = rhs._value_one_or_more(lhs)
    if res is not None:
        return False

@EqualityPredicate.register(Permutation, Permutation)
def _(lhs, rhs, assumptions):
    if lhs._size != rhs._size:
        return None
    return lhs._array_form == rhs._array_form

@EqualityPredicate.register(ImmutableDenseMatrix, ImmutableDenseMatrix)
def _(lhs, rhs, assumptions):
    if lhs.shape != rhs.shape:
        return False
    return (lhs - rhs).is_zero_matrix

@EqualityPredicate.register(MatrixExpr, MatrixExpr)
def _(lhs, rhs, assumptions):
    if lhs.shape != rhs.shape:
        return False
    if (lhs - rhs).is_ZeroMatrix:
        return True

@EqualityPredicate.register(Expr, MatrixExpr)
def _(lhs, rhs, assumptions):
    return False

@EqualityPredicate.register(ComplexRootOf, ComplexRootOf)
def _(lhs, rhs, assumptions):
    return lhs == rhs

@EqualityPredicate.register(Basic, ComplexRootOf)
def _(lhs, rhs, assumptions):
    # CRootOf represents a Root, so if lhs is that root, it should set
    # the expression to zero *and* it should be in the interval of the
    # CRootOf instance. It must also be a number that agrees with the
    # is_real value of the CRootOf instance.
    if not lhs.is_number:
        return None
    if not lhs.is_finite:
        return False
    z = rhs.expr.subs(rhs.expr.free_symbols.pop(), lhs).is_zero
    if z is False:  # all roots will make z True but we don't know
        # whether this is the right root if z is True
        return False
    o = lhs.is_real, lhs.is_imaginary
    s = rhs.is_real, rhs.is_imaginary
    assert None not in s  # this is part of initial refinement
    if o != s and None not in o:
        return False
    re, im = lhs.as_real_imag()
    if rhs.is_real:
        if im:
            return False
        i = rhs._get_interval()
        a, b = [Rational(str(_)) for _ in (i.a, i.b)]
        return sympify(a <= lhs and lhs <= b)
    i = rhs._get_interval()
    r1, r2, i1, i2 = [Rational(str(j)) for j in (
        i.ax, i.bx, i.ay, i.by)]
    return is_le(r1, re) and is_le(re,r2) and is_le(i1,im) and is_le(im,i2)

@EqualityPredicate.register(Set, Set) # type:ignore
def _(lhs, rhs, assumptions):
    return tfn[fuzzy_and(a.is_subset(b) for a, b in [(lhs, rhs), (rhs, lhs)])]

@EqualityPredicate.register(Basic, Set) # type:ignore
def _(lhs, rhs, assumptions):
    return False

@EqualityPredicate.register(Interval, FiniteSet) # type:ignore
@EqualityPredicate.register(FiniteSet, Interval)
def _(lhs, rhs, assumptions):
    return False

@EqualityPredicate.register(Interval, Interval) # type:ignore
def _(lhs, rhs, assumptions):
    lhs_left = ask(Q.eq(lhs.left, rhs.left), assumptions)
    lhs_right = ask(Q.eq(lhs.right, rhs.right), assumptions)

    if lhs_left is None or lhs_right is None:
        return None

    return And(lhs_left,
               lhs_right,
               lhs.left_open == rhs.left_open,
               lhs.right_open == rhs.right_open)

@EqualityPredicate.register(FiniteSet, FiniteSet) # type:ignore
def _(lhs, rhs, assumptions):
    def all_in_both():
        s_set = set(lhs.args)
        o_set = set(rhs.args)
        yield fuzzy_and(lhs._contains(e) for e in o_set - s_set)
        yield fuzzy_and(rhs._contains(e) for e in s_set - o_set)

    return tfn[fuzzy_and(all_in_both())]

@EqualityPredicate.register(ProductSet, ProductSet) # type:ignore
def _(lhs, rhs, assumptions):
    if len(lhs.sets) != len(rhs.sets):
        return False

    eqs = (ask(Q.eq(x, y), assumptions) for x, y in zip(lhs.sets, rhs.sets))
    return tfn[fuzzy_and(map(fuzzy_bool, eqs))]


### UnequalityPredicate ###

@UnequalityPredicate.register(object, object)
def _(lhs, rhs, assumptions):
    return fuzzy_not(ask(Q.eq(lhs, rhs), assumptions))
