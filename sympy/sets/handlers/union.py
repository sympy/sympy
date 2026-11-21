from __future__ import annotations
from sympy.core.numbers import oo
from sympy.core.singleton import S
from sympy.core.sympify import sympify
from sympy.functions.elementary.miscellaneous import Min, Max
from sympy.sets.sets import (EmptySet, FiniteSet, Intersection,
    Interval, ProductSet, Set, Union, UniversalSet)
from sympy.sets.fancysets import (ComplexRegion, Naturals, Naturals0,
    Integers, Rationals, Reals, Range)
from sympy.multipledispatch import Dispatcher


union_sets = Dispatcher('union_sets')


@union_sets.register(Naturals0, Naturals)
def _(a, b):
    return a

@union_sets.register(Rationals, Naturals)
def _(a, b):
    return a

@union_sets.register(Rationals, Naturals0)
def _(a, b):
    return a

@union_sets.register(Reals, Naturals)
def _(a, b):
    return a

@union_sets.register(Reals, Naturals0)
def _(a, b):
    return a

@union_sets.register(Reals, Rationals)
def _(a, b):
    return a

@union_sets.register(Integers, Set)
def _(a, b):
    intersect = Intersection(a, b)
    if intersect == a:
        return b
    elif intersect == b:
        return a

@union_sets.register(ComplexRegion, Set)
def _(a, b):
    if b.is_subset(S.Reals):
        # treat a subset of reals as a complex region
        b = ComplexRegion.from_real(b)

    if b.is_ComplexRegion:
        # a in rectangular form
        if (not a.polar) and (not b.polar):
            return ComplexRegion(Union(a.sets, b.sets))
        # a in polar form
        elif a.polar and b.polar:
            return ComplexRegion(Union(a.sets, b.sets), polar=True)
    return None

@union_sets.register(EmptySet, Set)
def _(a, b):
    return b


@union_sets.register(UniversalSet, Set)
def _(a, b):
    return a

@union_sets.register(ProductSet, ProductSet)
def _(a, b):
    if b.is_subset(a):
        return a
    if len(b.sets) != len(a.sets):
        return None
    if len(a.sets) == 2:
        a1, a2 = a.sets
        b1, b2 = b.sets
        if a1 == b1:
            return a1 * Union(a2, b2)
        if a2 == b2:
            return Union(a1, b1) * a2
    return None

@union_sets.register(ProductSet, Set)
def _(a, b):
    if b.is_subset(a):
        return a
    return None

@union_sets.register(Interval, Interval)
def _(a, b):
    if a._is_comparable(b):
        # Non-overlapping intervals
        end = Min(a.end, b.end)
        start = Max(a.start, b.start)
        if (end < start or
           (end == start and (end not in a and end not in b))):
            return None
        else:
            start = Min(a.start, b.start)
            end = Max(a.end, b.end)

            left_open = ((a.start != start or a.left_open) and
                         (b.start != start or b.left_open))
            right_open = ((a.end != end or a.right_open) and
                          (b.end != end or b.right_open))
            return Interval(start, end, left_open, right_open)

@union_sets.register(Interval, UniversalSet)
def _(a, b):
    return S.UniversalSet

@union_sets.register(Interval, Set)
def _(a, b):
    # If I have open end points and these endpoints are contained in b
    # But only in case, when endpoints are finite. Because
    # interval does not contain oo or -oo.
    open_left_in_b_and_finite = (a.left_open and
                                     sympify(b.contains(a.start)) is S.true and
                                     a.start.is_finite)
    open_right_in_b_and_finite = (a.right_open and
                                      sympify(b.contains(a.end)) is S.true and
                                      a.end.is_finite)
    if open_left_in_b_and_finite or open_right_in_b_and_finite:
        # Fill in my end points and return
        open_left = a.left_open and not open_left_in_b_and_finite
        open_right = a.right_open and not open_right_in_b_and_finite
        new_a = Interval(a.start, a.end, open_left, open_right)
        return {new_a, b}
    return None

@union_sets.register(FiniteSet, FiniteSet)
def _(a, b):
    return FiniteSet(*(a._elements | b._elements))

@union_sets.register(FiniteSet, Set)
def _(a, b):
    # If `b` set contains one of my elements, remove it from `a`
    if any(b.contains(x) == True for x in a):
        return {
            FiniteSet(*[x for x in a if b.contains(x) != True]), b}
    return None

@union_sets.register(Range, Range)
def _(a, b):
    if a.is_empty:
        return b
    if b.is_empty:
        return a
    if a == b:
        return a
    if a.is_subset(b):
        return b
    if b.is_subset(a):
        return a
    if not all(all(v.is_number for v in r.args) for r in [a, b]):
        return None

    if a.step < 0 and a.size.is_finite:
        a = a.reversed
    if b.step < 0 and b.size.is_finite:
        b = b.reversed

    if a.step == b.step:
        st = a.step
        if st > 0:
            elem_a = a.start if a.start.is_finite else a.sup
            elem_b = b.start if b.start.is_finite else b.sup
            if (elem_b - elem_a) % st != 0:
                return None

            if a.inf > b.inf:
                a, b = b, a
            if b.inf <= a.sup + st:
                new_inf = a.inf
                new_sup = Max(a.sup, b.sup)
                new_stop = oo if new_sup is oo else new_sup + st
                return Range(new_inf, new_stop, st)

        elif st < 0:
            elem_a = a.start if a.start.is_finite else a.inf
            elem_b = b.start if b.start.is_finite else b.inf
            if (elem_b - elem_a) % st != 0:
                return None

            if a.sup < b.sup:
                a, b = b, a
            if b.sup >= a.inf + st:
                new_sup = a.sup
                new_inf = Min(a.inf, b.inf)
                new_stop = -oo if new_inf is -oo else new_inf + st
                return Range(new_sup, new_stop, st)

    return None

@union_sets.register(Set, Set)
def _(a, b):
    return None
