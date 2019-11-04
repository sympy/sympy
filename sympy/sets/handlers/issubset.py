from sympy import S
from sympy.core.logic import fuzzy_and, fuzzy_bool, fuzzy_not, fuzzy_or
from sympy.core.relational import Eq
from sympy.sets.sets import FiniteSet, Interval, Set, Union
from sympy.sets.fancysets import Reals
from sympy.multipledispatch import dispatch


_inf_sets = [S.Naturals, S.Naturals0, S.Integers, S.Rationals, S.Reals, S.Complexes]

@dispatch(Set, Set)
def is_subset_sets(a, b):
    if a in _inf_sets:
        if b in _inf_sets:
            return _inf_sets.index(a) <= _inf_sets.index(b)
        if isinstance(b, FiniteSet):
            return False

@dispatch(Interval, Interval)
def is_subset_sets(a, b):
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
def is_subset_sets(a_interval, b_fs):
    # An Interval can only be a subset of a finite set if it is finite
    # which can only happen if it has zero measure.
    if fuzzy_not(a_interval.measure.is_zero):
        return False

@dispatch(Interval, Union)
def is_subset_sets(a_interval, b_u):
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

@dispatch(FiniteSet, Set)
def is_subset_sets(a_fs, b):
    return fuzzy_and(b._contains(e) for e in a_fs.args)

@dispatch(Union, Set)
def is_subset_sets(a_u, b):
    return fuzzy_and(s.is_subset(b) for s in a_u.args)
