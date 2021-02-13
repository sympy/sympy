from sympy import (Interval, Intersection, Set, EmptySet, S, sympify,
                   FiniteSet, Union, ComplexRegion, ProductSet)
from sympy.multipledispatch import dispatch
from sympy.sets.fancysets import (Naturals, Naturals0, Integers, Rationals,
                                  Reals)
from sympy.sets.sets import UniversalSet


@dispatch(Naturals0, Naturals)  # type: ignore # noqa:F811
def union_sets(a, b): # noqa:F811
    return a

@dispatch(Rationals, Naturals)  # type: ignore # noqa:F811
def union_sets(a, b): # noqa:F811
    return a

@dispatch(Rationals, Naturals0)  # type: ignore # noqa:F811
def union_sets(a, b): # noqa:F811
    return a

@dispatch(Reals, Naturals)  # type: ignore # noqa:F811
def union_sets(a, b): # noqa:F811
    return a

@dispatch(Reals, Naturals0)  # type: ignore # noqa:F811
def union_sets(a, b): # noqa:F811
    return a

@dispatch(Reals, Rationals)  # type: ignore # noqa:F811
def union_sets(a, b): # noqa:F811
    return a

@dispatch(Integers, Set)  # type: ignore # noqa:F811
def union_sets(a, b): # noqa:F811
    intersect = Intersection(a, b)
    if intersect == a:
        return b
    elif intersect == b:
        return a

@dispatch(ComplexRegion, Set)  # type: ignore # noqa:F811
def union_sets(a, b): # noqa:F811
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

@dispatch(type(EmptySet), Set)  # type: ignore # noqa:F811
def union_sets(a, b): # noqa:F811
    return b


@dispatch(UniversalSet, Set)  # type: ignore # noqa:F811
def union_sets(a, b): # noqa:F811
    return a

@dispatch(ProductSet, ProductSet)  # type: ignore # noqa:F811
def union_sets(a, b): # noqa:F811
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

@dispatch(ProductSet, Set)  # type: ignore # noqa:F811
def union_sets(a, b): # noqa:F811
    if b.is_subset(a):
        return a
    return None

@dispatch(Interval, Interval)  # type: ignore # noqa:F811
def union_sets(a, b): # noqa:F811
    if a._is_comparable(b):
        from sympy.functions.elementary.miscellaneous import Min, Max
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

@dispatch(Interval, UniversalSet)  # type: ignore # noqa:F811
def union_sets(a, b): # noqa:F811
    return S.UniversalSet

@dispatch(Interval, Set)  # type: ignore # noqa:F811
def union_sets(a, b): # noqa:F811
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
        open_left = a.left_open and a.start not in b
        open_right = a.right_open and a.end not in b
        new_a = Interval(a.start, a.end, open_left, open_right)
        return {new_a, b}
    return None

@dispatch(FiniteSet, FiniteSet)  # type: ignore # noqa:F811
def union_sets(a, b): # noqa:F811
    return FiniteSet(*(a._elements | b._elements))

@dispatch(FiniteSet, Set)  # type: ignore # noqa:F811
def union_sets(a, b): # noqa:F811
    # If `b` set contains one of my elements, remove it from `a`
    if any(b.contains(x) == True for x in a):
        return {
            FiniteSet(*[x for x in a if b.contains(x) != True]), b}
    return None

@dispatch(Set, Set)  # type: ignore # noqa:F811
def union_sets(a, b): # noqa:F811
    return None
