from sympy import (Interval, Intersection, Set, EmptySet, S, sympify,
                   FiniteSet, Union, ComplexRegion, ProductSet)
from sympy.core.expr import Expr
from sympy.core.function import Lambda
from sympy.core.symbol import Dummy
from sympy.multipledispatch import dispatch
from sympy.polys.polytools import cancel
from sympy.sets.fancysets import (Naturals, Naturals0, Integers, Rationals,
                                  Reals, ImageSet)
from sympy.sets.sets import UniversalSet, imageset


@dispatch(Naturals0, Naturals)  # type: ignore
def union_sets(a, b): # noqa:F811
    return a

@dispatch(Rationals, Naturals)  # type: ignore
def union_sets(a, b): # noqa:F811
    return a

@dispatch(Rationals, Naturals0)  # type: ignore
def union_sets(a, b): # noqa:F811
    return a

@dispatch(Reals, Naturals)  # type: ignore
def union_sets(a, b): # noqa:F811
    return a

@dispatch(Reals, Naturals0)  # type: ignore
def union_sets(a, b): # noqa:F811
    return a

@dispatch(Reals, Rationals)  # type: ignore
def union_sets(a, b): # noqa:F811
    return a

@dispatch(Integers, Set)  # type: ignore
def union_sets(a, b): # noqa:F811
    intersect = Intersection(a, b)
    if intersect == a:
        return b
    elif intersect == b:
        return a

@dispatch(ComplexRegion, Set)  # type: ignore
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

@dispatch(type(EmptySet), Set)  # type: ignore
def union_sets(a, b): # noqa:F811
    return b


@dispatch(UniversalSet, Set)  # type: ignore
def union_sets(a, b): # noqa:F811
    return a

@dispatch(ProductSet, ProductSet)  # type: ignore
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

@dispatch(ProductSet, Set)  # type: ignore
def union_sets(a, b): # noqa:F811
    if b.is_subset(a):
        return a
    return None

@dispatch(Interval, Interval)  # type: ignore
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

@dispatch(Interval, UniversalSet)  # type: ignore
def union_sets(a, b): # noqa:F811
    return S.UniversalSet

@dispatch(Interval, Set)  # type: ignore
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
        return set((new_a, b))
    return None

@dispatch(FiniteSet, FiniteSet)  # type: ignore
def union_sets(a, b): # noqa:F811
    return FiniteSet(*(a._elements | b._elements))

@dispatch(FiniteSet, Set)  # type: ignore
def union_sets(a, b): # noqa:F811
    # If `b` set contains one of my elements, remove it from `a`
    if any(b.contains(x) == True for x in a):
        return set((
            FiniteSet(*[x for x in a if b.contains(x) != True]), b))
    return None

@dispatch(ImageSet, ImageSet)  # type: ignore
def union_sets(a, b): # noqa:F811
    """Simplify union of two arithmetic sequences if possible.

    Given two ``ImageSet``s *a* and *b* representing arithmetic sequences
    respectively defined by $n \mapsto p n + q$ and $n \mapsto r n + s$, with
    $p, q, r, s$ complex numbers, ``union_sets()`` checks whether it is
    possible to reduce the two arithmetic sequences to a single sequence and,
    if possible, returns a new ``ImageSet`` based on this combined arithmetic
    sequence.

    No other type of ``ImageSet`` unions is currently supported.

    """
    if (a.base_set != S.Integers or b.base_set != S.Integers or
        not all(isinstance(imgset.lamda.expr, Expr) for imgset in (a, b))):
            return None

    def extract_linear_coeffs(imgset):
        """Returns $p$, $q$ in $n \mapsto p n + q$; else None."""
        expr = imgset.lamda.expr
        var = imgset.lamda.variables
        if len(var) > 1:
            return None
        var = var[0]
        poly = expr.as_poly(var)
        if poly and poly.is_linear:
            return poly.LC(), poly.TC()

    pq = extract_linear_coeffs(a)
    if not pq:
        return
    rs = extract_linear_coeffs(b)
    if not rs:
        return
    p, q, r, s = *pq, *rs

    def in_sequence(x, p, q):
        # x == p*n + q iff n == (x - q)/p
        return cancel((x - q)/p).is_integer == True
    def is_subset(p, q, r, s):
        # checks whether {p*n + q} is a subset of {r*m + s}
        if not cancel(p/r).is_integer:
            return False
        return in_sequence(q, r, s)

    # 1st possibility: `a` is a subset of `b` or vice-versa
    if is_subset(p, q, r, s):
        return b
    elif is_subset(r, s, p, q):
        return a
    else:
        # 2nd possibility: Both sequences are alternating, i.e. they have the
        # same period and the midpoint of `a`'s period interval must be
        # contained in the second sequence.
        midpoint = q + p/2
        if in_sequence(midpoint, r, s):
            u, v = p/2, q
            n = Dummy('n')
            return imageset(Lambda(n, u*n + v), S.Integers)
    return None

@dispatch(Set, Set)  # type: ignore
def union_sets(a, b): # noqa:F811
    return None
