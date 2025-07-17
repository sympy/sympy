from sympy.core.singleton import S
from sympy.core.symbol import Symbol
from sympy.core.logic import fuzzy_and, fuzzy_bool, fuzzy_not, fuzzy_or
from sympy.sets.sets import FiniteSet, Interval, Set, Union, ProductSet
from sympy.sets.fancysets import Complexes, Reals, Range, Rationals, Integers, ImageSet
from sympy.multipledispatch import Dispatcher


_inf_sets = [S.Naturals, S.Naturals0, S.Integers, S.Rationals, S.Reals, S.Complexes]


is_subset_sets = Dispatcher('is_subset_sets')


@is_subset_sets.register(Set, Set)
def _(a, b):
    return None

@is_subset_sets.register(Interval, Interval)
def _(a, b):

    if b.left_open and not a.left_open:
        # [a1, a2] <= (b1, b2]
        left_in = fuzzy_bool(a.start > b.start)
    else:
        # (a1, a2] <= (b1, b2]
        # (a1, a2] <= [b1, b2]
        # [a1, a2] <= [b1, b2]
        left_in = fuzzy_bool(a.start >= b.start)

    if left_in is False:
        return False

    if b.right_open and not a.right_open:
        # [a1, a2] <= [b1, b2)
        right_in = fuzzy_bool(a.end < b.end)
    else:
        # [a1, a2) <= [b1, b2)
        # [a1, a2) <= [b1, b2]
        # [a1, a2] <= [b1, b2]
        right_in = fuzzy_bool(a.end <= b.end)

    return fuzzy_and([left_in, right_in])

@is_subset_sets.register(Interval, FiniteSet)
def _(a_interval, b_fs):
    # An Interval can only be a subset of a finite set if it is finite
    # which can only happen if it has zero measure.
    if fuzzy_not(a_interval.measure.is_zero):
        return False

@is_subset_sets.register(Interval, Union)
def _(a_interval, b_u):
    if fuzzy_or(a_interval.is_subset(s) for s in b_u.args):
        return True
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

@is_subset_sets.register(Range, Range)
def _(a, b):
    if a.step == b.step == 1:
        return fuzzy_and([fuzzy_bool(a.start >= b.start),
                          fuzzy_bool(a.stop <= b.stop)])

@is_subset_sets.register(Range, Interval)
def _(a_range, b_interval):
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

@is_subset_sets.register(Range, FiniteSet)
def _(a_range, b_finiteset):
    try:
        a_size = a_range.size
    except ValueError:
        # symbolic Range of unknown size
        return None
    if a_size > len(b_finiteset):
        return False
    elif any(arg.has(Symbol) for arg in a_range.args):
        return fuzzy_and(b_finiteset.contains(x) for x in a_range)
    else:
        # Checking A \ B == EmptySet is more efficient than repeated naive
        # membership checks on an arbitrary FiniteSet.
        a_set = set(a_range)
        b_remaining = len(b_finiteset)
        # Symbolic expressions and numbers of unknown type (integer or not) are
        # all counted as "candidates", i.e. *potentially* matching some a in
        # a_range.
        cnt_candidate = 0
        for b in b_finiteset:
            if b.is_Integer:
                a_set.discard(b)
            elif fuzzy_not(b.is_integer):
                pass
            else:
                cnt_candidate += 1
            b_remaining -= 1
            if len(a_set) > b_remaining + cnt_candidate:
                return False
            if len(a_set) == 0:
                return True
        return None

@is_subset_sets.register(ImageSet, Reals)
def _(a_imageset, b_reals):
    f = a_imageset.lamda
    base = a_imageset.base_set
    var = f.variables
    if base.is_subset(b_reals):
        if f.expr.as_poly(*var, domain='QQ') is not None:
            return True

@is_subset_sets.register(Interval, Range)
def _(a_interval, b_range):
    if a_interval.measure.is_extended_nonzero:
        return False

@is_subset_sets.register(Interval, Rationals)
def _(a_interval, b_rationals):
    if a_interval.measure.is_extended_nonzero:
        return False

@is_subset_sets.register(Range, Complexes)
def _(a, b):
    return True

@is_subset_sets.register(Complexes, Interval)
def _(a, b):
    return False

@is_subset_sets.register(Complexes, Range)
def _(a, b):
    return False

@is_subset_sets.register(Complexes, Rationals)
def _(a, b):
    return False

@is_subset_sets.register(Rationals, Complexes)
def _(a, b):
    return True

@is_subset_sets.register(Reals, Complexes)
def _(a, b):
    return True

@is_subset_sets.register(Interval, Reals)
def _(a, b):
    return True

@is_subset_sets.register(Integers, Rationals)
def _(a, b):
    return True

@is_subset_sets.register(Rationals, Reals)
def _(a, b):
    return True

@is_subset_sets.register(Rationals, Range)
def _(a, b):
    return False

@is_subset_sets.register(Range, Rationals)
def _(a, b):
    return True

@is_subset_sets.register(ProductSet, ProductSet)
def _(a, b):
    if len(a.sets) != len(b.sets):
        return False
    return fuzzy_and(a_i.is_subset(b_i) for a_i, b_i in zip(a.sets, b.sets))

@is_subset_sets.register(ProductSet, FiniteSet)
def _(a_ps, b_fs):
    return fuzzy_and(b_fs.contains(x) for x in a_ps)
