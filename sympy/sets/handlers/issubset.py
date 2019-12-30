from sympy import S
from sympy.core.logic import fuzzy_and, fuzzy_bool, fuzzy_not, fuzzy_or
from sympy.core.relational import Eq
from sympy.sets.sets import FiniteSet, Interval, Set, Union
from sympy.sets.fancysets import Complexes, Reals, Range, Rationals
from sympy.multipledispatch import dispatch


_inf_sets = [S.Naturals, S.Naturals0, S.Integers, S.Rationals, S.Reals, S.Complexes]

@dispatch(Set, Set)
def is_subset_sets(a, b): # noqa:F811
    return None

@dispatch(Interval, Interval)
def is_subset_sets(a, b): # noqa:F811
    # This is correct but can be made more comprehensive...
    if fuzzy_bool(a.start < b.start):
        return False
    if fuzzy_bool(a.end > b.end):
        return False
    if (b.left_open and not a.left_open and fuzzy_bool(Eq(a.start, b.start))):
        return False
    if (b.right_open and not a.right_open and fuzzy_bool(Eq(a.end, b.end))):
        return False

@dispatch(Interval, FiniteSet)
def is_subset_sets(a_interval, b_fs): # noqa:F811
    # An Interval can only be a subset of a finite set if it is finite
    # which can only happen if it has zero measure.
    if fuzzy_not(a_interval.measure.is_zero):
        return False

@dispatch(Interval, Union)
def is_subset_sets(a_interval, b_u): # noqa:F811
    if all(isinstance(s, (Interval, FiniteSet)) for s in b_u.args):
        intervals = [s for s in b_u.args if isinstance(s, Interval)]
        if all(fuzzy_bool(a_interval.start < s.start) for s in intervals):
            return False
        if all(fuzzy_bool(a_interval.end > s.end) for s in intervals):
            return False
        if a_interval.measure.is_nonzero:
            no_overlap = lambda s1, s2: fuzzy_or([
                    fuzzy_bool(s1.end <= s2.start),
                    fuzzy_bool(s1.start >= s2.end),
                    ])
            if all(no_overlap(s, a_interval) for s in intervals):
                return False

@dispatch(Range, Range)
def is_subset_sets(a, b): # noqa:F811
    if a.step == b.step == 1:
        return fuzzy_and([fuzzy_bool(a.start >= b.start),
                          fuzzy_bool(a.stop <= b.stop)])

@dispatch(Range, Interval)
def is_subset_sets(a_range, b_interval): # noqa:F811
    if a_range.step.is_positive:
        if b_interval.left_open and a_range.inf.is_finite:
            cond_left = a_range.inf > b_interval.left
        else:
            cond_left = a_range.inf >= b_interval.left
        if b_interval.right_open and a_range.sup.is_finite:
            cond_right = a_range.sup < b_interval.right
        else:
            cond_right = a_range.sup <= b_interval.right
        return fuzzy_and([cond_left, cond_right])

@dispatch(Interval, Range)
def is_subset_sets(a_interval, b_range): # noqa:F811
    if a_interval.measure.is_extended_nonzero:
        return False

@dispatch(Interval, Rationals)
def is_subset_sets(a_interval, b_rationals): # noqa:F811
    if a_interval.measure.is_extended_nonzero:
        return False

@dispatch(Range, Complexes)
def is_subset_sets(a, b): # noqa:F811
    return True

@dispatch(Complexes, Interval)
def is_subset_sets(a, b): # noqa:F811
    return False

@dispatch(Complexes, Range)
def is_subset_sets(a, b): # noqa:F811
    return False

@dispatch(Complexes, Rationals)
def is_subset_sets(a, b): # noqa:F811
    return False

@dispatch(Rationals, Reals)
def is_subset_sets(a, b): # noqa:F811
    return True

@dispatch(Rationals, Range)
def is_subset_sets(a, b): # noqa:F811
    return False
