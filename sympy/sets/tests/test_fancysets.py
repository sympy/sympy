from sympy.core.compatibility import range
from sympy.sets.fancysets import (ImageSet, Range, normalize_theta_set,
                                  ComplexRegion)
from sympy.sets.sets import (FiniteSet, Interval, imageset, EmptySet, Union,
                             Intersection)
from sympy import (S, Symbol, Lambda, symbols, cos, sin, pi, oo, Basic,
                   Rational, sqrt, tan, log, Abs, I)
from sympy.utilities.pytest import XFAIL, raises
import itertools

x = Symbol('x')


def test_naturals():
    N = S.Naturals
    assert 5 in N
    assert -5 not in N
    assert 5.5 not in N
    ni = iter(N)
    a, b, c, d = next(ni), next(ni), next(ni), next(ni)
    assert (a, b, c, d) == (1, 2, 3, 4)
    assert isinstance(a, Basic)

    assert N.intersect(Interval(-5, 5)) == Range(1, 6)
    assert N.intersect(Interval(-5, 5, True, True)) == Range(1, 5)

    assert N.boundary == N

    assert N.inf == 1
    assert N.sup == oo

def test_naturals0():
    N = S.Naturals0
    assert 0 in N
    assert -1 not in N
    assert next(iter(N)) == 0

def test_integers():
    Z = S.Integers
    assert 5 in Z
    assert -5 in Z
    assert 5.5 not in Z
    zi = iter(Z)
    a, b, c, d = next(zi), next(zi), next(zi), next(zi)
    assert (a, b, c, d) == (0, 1, -1, 2)
    assert isinstance(a, Basic)

    assert Z.intersect(Interval(-5, 5)) == Range(-5, 6)
    assert Z.intersect(Interval(-5, 5, True, True)) == Range(-4, 5)

    assert Z.inf == -oo
    assert Z.sup == oo

    assert Z.boundary == Z


def test_ImageSet():
    squares = ImageSet(Lambda(x, x**2), S.Naturals)
    assert 4 in squares
    assert 5 not in squares
    assert FiniteSet(*range(10)).intersect(squares) == FiniteSet(1, 4, 9)

    assert 16 not in squares.intersect(Interval(0, 10))

    si = iter(squares)
    a, b, c, d = next(si), next(si), next(si), next(si)
    assert (a, b, c, d) == (1, 4, 9, 16)

    harmonics = ImageSet(Lambda(x, 1/x), S.Naturals)
    assert Rational(1, 5) in harmonics
    assert Rational(.25) in harmonics
    assert 0.25 not in harmonics
    assert Rational(.3) not in harmonics

    assert harmonics.is_iterable

def test_image_is_ImageSet():
    assert isinstance(imageset(x, sqrt(sin(x)), Range(5)), ImageSet)


@XFAIL
def test_halfcircle():
    # This test sometimes works and sometimes doesn't.
    # It may be an issue with solve? Maybe with using Lambdas/dummys?
    # I believe the code within fancysets is correct
    r, th = symbols('r, theta', real=True)
    L = Lambda((r, th), (r*cos(th), r*sin(th)))
    halfcircle = ImageSet(L, Interval(0, 1)*Interval(0, pi))

    assert (1, 0) in halfcircle
    assert (0, -1) not in halfcircle
    assert (0, 0) in halfcircle

    assert not halfcircle.is_iterable


def test_ImageSet_iterator_not_injetive():
    L = Lambda(x, x - x % 2)  # produces 0, 2, 2, 4, 4, 6, 6, ...
    evens = ImageSet(L, S.Naturals)
    i = iter(evens)
    # No repeats here
    assert (next(i), next(i), next(i), next(i)) == (0, 2, 4, 6)


def test_Range():
    assert Range(5) == Range(0, 5) == Range(0, 5, 1)

    r = Range(10, 20, 2)
    assert 12 in r
    assert 8 not in r
    assert 11 not in r
    assert 30 not in r

    assert list(Range(0, 5)) == list(range(5))
    assert list(Range(5, 0, -1)) == list(range(1, 6))

    assert Range(0, 10, -1) == S.EmptySet

    assert Range(5, 15).sup == 14
    assert Range(5, 15).inf == 5
    assert Range(15, 5, -1).sup == 15
    assert Range(15, 5, -1).inf == 6
    assert Range(10, 67, 10).sup == 60
    assert Range(60, 7, -10).inf == 10

    assert len(Range(10, 38, 10)) == 3
    assert Range(0, 0, 5) == S.EmptySet

    assert Range(1, 1) == S.EmptySet
    raises(ValueError, lambda: Range(0, oo, oo))
    raises(ValueError, lambda: Range(-oo, oo))
    raises(ValueError, lambda: Range(-oo, oo, 2))
    raises(ValueError, lambda: Range(0, pi, 1))

    assert 5 in Range(0, oo, 5)
    assert -5 in Range(-oo, 0, 5)

    assert Range(0, oo)
    assert Range(-oo, 0)
    assert Range(0, -oo, -1)

    assert Range(0, oo, 2)._last_element is oo
    assert Range(-oo, 1, 1)._last_element is S.Zero

    it = iter(Range(-oo, 0, 2))
    assert (next(it), next(it)) == (-2, -4)

    assert Range(-1, 10, 1).intersect(S.Integers) == Range(-1, 10, 1)
    assert Range(-1, 10, 1).intersect(S.Naturals) == Range(1, 10, 1)

    assert Range(1, 10, 1)._ith_element(5) == 6 # the index starts from zero
    assert Range(1, 10, 1)._last_element == 9

    assert Range(1, 10, 1).boundary == Range(1, 10, 1)


def test_range_interval_intersection():
    # Intersection with intervals
    assert FiniteSet(*Range(0, 10, 1).intersect(Interval(2, 6))) == \
        FiniteSet(2, 3, 4, 5, 6)

    # Open Intervals are removed
    assert (FiniteSet(*Range(0, 10, 1).intersect(Interval(2, 6, True, True)))
            == FiniteSet(3, 4, 5))

    # Try this with large steps
    assert (FiniteSet(*Range(0, 100, 10).intersect(Interval(15, 55))) ==
            FiniteSet(20, 30, 40, 50))

    # Going backwards
    assert FiniteSet(*Range(10, -9, -3).intersect(Interval(-5, 6))) == \
        FiniteSet(-5, -2, 1, 4)
    assert FiniteSet(*Range(10, -9, -3).intersect(Interval(-5, 6, True))) == \
        FiniteSet(-2, 1, 4)


def test_fun():
    assert (FiniteSet(*ImageSet(Lambda(x, sin(pi*x/4)),
        Range(-10, 11))) == FiniteSet(-1, -sqrt(2)/2, 0, sqrt(2)/2, 1))


def test_Reals():
    assert 5 in S.Reals
    assert S.Pi in S.Reals
    assert -sqrt(2) in S.Reals
    assert (2, 5) not in S.Reals
    assert sqrt(-1) not in S.Reals
    assert S.Reals == Interval(-oo, oo)
    assert S.Reals != Interval(0, oo)


def test_Complex():
    assert 5 in S.Complexes
    assert 5 + 4*I in S.Complexes
    assert S.Pi in S.Complexes
    assert -sqrt(2) in S.Complexes
    assert -I in S.Complexes
    assert sqrt(-1) in S.Complexes
    assert S.Complexes.intersect(S.Reals) == S.Reals
    # assert S.Complexes.union(S.Reals) == S.Complexes
    assert S.Complexes == ComplexRegion(S.Reals*S.Reals)


def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(itertools.islice(iterable, n))


def test_intersections():
    assert S.Integers.intersect(S.Reals) == S.Integers
    assert 5 in S.Integers.intersect(S.Reals)
    assert 5 in S.Integers.intersect(S.Reals)
    assert -5 not in S.Naturals.intersect(S.Reals)
    assert 5.5 not in S.Integers.intersect(S.Reals)
    assert 5 in S.Integers.intersect(Interval(3, oo))
    assert -5 in S.Integers.intersect(Interval(-oo, 3))
    assert all(x.is_Integer
            for x in take(10, S.Integers.intersect(Interval(3, oo)) ))


def test_infinitely_indexed_set_1():
    from sympy.abc import n, m, t
    assert imageset(Lambda(n, n), S.Integers) == imageset(Lambda(m, m), S.Integers)

    assert imageset(Lambda(n, 2*n), S.Integers).intersect(imageset(Lambda(m, 2*m + 1), S.Integers)) == \
            EmptySet()

    assert imageset(Lambda(n, 2*n), S.Integers).intersect(imageset(Lambda(n, 2*n + 1), S.Integers)) == \
            EmptySet()

    assert imageset(Lambda(m, 2*m), S.Integers).intersect(imageset(Lambda(n, 3*n), S.Integers)) == \
            ImageSet(Lambda(t, 6*t), S.Integers)


def test_infinitely_indexed_set_2():
    from sympy import exp
    from sympy.abc import n
    a = Symbol('a', integer=True)
    assert imageset(Lambda(n, n), S.Integers) == imageset(Lambda(n, n + a), S.Integers)
    assert imageset(Lambda(n, n), S.Integers) == imageset(Lambda(n, -n + a), S.Integers)
    assert imageset(Lambda(n, -6*n), S.Integers) == ImageSet(Lambda(n, 6*n), S.Integers)
    assert imageset(Lambda(n, 2*n + pi), S.Integers) == ImageSet(Lambda(n, 2*n + pi), S.Integers)
    assert imageset(Lambda(n, pi*n + pi), S.Integers) == ImageSet(Lambda(n, pi*n + pi), S.Integers)
    assert imageset(Lambda(n, exp(n)), S.Integers) != imageset(Lambda(n, n), S.Integers)


def test_imageset_intersect_real():
    from sympy import I
    from sympy.abc import n
    assert imageset(Lambda(n, n + (n - 1)*(n + 1)*I), S.Integers).intersect(S.Reals) == \
            FiniteSet(-1, 1)

    s = ImageSet(Lambda(n, -I*(I*(2*pi*n - pi/4) + log(Abs(sqrt(-I))))), S.Integers)
    assert s.intersect(S.Reals) == imageset(Lambda(n, 2*n*pi - pi/4), S.Integers)


@XFAIL
def test_infinitely_indexed_failed_diophantine():
    from sympy.abc import n, m, t
    assert imageset(Lambda(m, 2*pi*m), S.Integers).intersect(imageset(Lambda(n, 3*pi*n), S.Integers)) == \
            ImageSet(Lambda(t, -6*pi*t), S.Integers)


@XFAIL
def test_infinitely_indexed_set_3():
    from sympy.abc import n
    assert imageset(Lambda(n, 2*n + 1), S.Integers) == imageset(Lambda(n, 2*n - 1), S.Integers)
    assert imageset(Lambda(n, 3*n + 2), S.Integers) == imageset(Lambda(n, 3*n - 1), S.Integers)


def test_ImageSet_simplification():
    from sympy.abc import n, m
    assert imageset(Lambda(n, n), S.Integers) == S.Integers
    assert imageset(Lambda(n, sin(n)),
                    imageset(Lambda(m, tan(m)), S.Integers)) == \
            imageset(Lambda(m, sin(tan(m))), S.Integers)


def test_ComplexRegion_contains():

    # contains in ComplexRegion
    a = Interval(2, 3)
    b = Interval(4, 6)
    c = Interval(7, 9)
    c1 = ComplexRegion(a*b)
    c2 = ComplexRegion(Union(a*b, c*a))
    assert 2.5 + 4.5*I in c1
    assert 2 + 4*I in c1
    assert 3 + 4*I in c1
    assert 8 + 2.5*I in c2
    assert 2.5 + 6.1*I not in c1
    assert 4.5 + 3.2*I not in c1

    r1 = Interval(0, 1)
    theta1 = Interval(0, 2*S.Pi)
    c3 = ComplexRegion(r1*theta1, polar=True)
    assert 0.5 + 0.6*I in c3
    assert I in c3
    assert 1 in c3
    assert 0 in c3
    assert 1 + I not in c3
    assert 1 - I not in c3


def test_ComplexRegion_intersect():

    # Polar form
    X_axis = ComplexRegion(Interval(0, oo)*FiniteSet(0, S.Pi), polar=True)

    unit_disk = ComplexRegion(Interval(0, 1)*Interval(0, 2*S.Pi), polar=True)
    upper_half_unit_disk = ComplexRegion(Interval(0, 1)*Interval(0, S.Pi), polar=True)
    upper_half_disk = ComplexRegion(Interval(0, oo)*Interval(0, S.Pi), polar=True)
    lower_half_disk = ComplexRegion(Interval(0, oo)*Interval(S.Pi, 2*S.Pi), polar=True)
    right_half_disk = ComplexRegion(Interval(0, oo)*Interval(-S.Pi/2, S.Pi/2), polar=True)
    first_quad_disk = ComplexRegion(Interval(0, oo)*Interval(0, S.Pi/2), polar=True)

    assert upper_half_disk.intersect(unit_disk) == upper_half_unit_disk
    assert right_half_disk.intersect(first_quad_disk) == first_quad_disk
    assert upper_half_disk.intersect(right_half_disk) == first_quad_disk
    assert upper_half_disk.intersect(lower_half_disk) == X_axis

    c1 = ComplexRegion(Interval(0, 4)*Interval(0, 2*S.Pi), polar=True)
    assert c1.intersect(Interval(1, 5)) == Interval(1, 4)
    assert c1.intersect(Interval(4, 9)) == FiniteSet(4)
    assert c1.intersect(Interval(5, 12)) == EmptySet()

    # Rectangular form
    X_axis = ComplexRegion(Interval(-oo, oo)*FiniteSet(0))

    unit_square = ComplexRegion(Interval(-1, 1)*Interval(-1, 1))
    upper_half_unit_square = ComplexRegion(Interval(-1, 1)*Interval(0, 1))
    upper_half_plane = ComplexRegion(Interval(-oo, oo)*Interval(0, oo))
    lower_half_plane = ComplexRegion(Interval(-oo, oo)*Interval(-oo, 0))
    right_half_plane = ComplexRegion(Interval(0, oo)*Interval(-oo, oo))
    first_quad_plane = ComplexRegion(Interval(0, oo)*Interval(0, oo))

    assert upper_half_plane.intersect(unit_square) == upper_half_unit_square
    assert right_half_plane.intersect(first_quad_plane) == first_quad_plane
    assert upper_half_plane.intersect(right_half_plane) == first_quad_plane
    assert upper_half_plane.intersect(lower_half_plane) == X_axis

    c1 = ComplexRegion(Interval(-5, 5)*Interval(-10, 10))
    assert c1.intersect(Interval(2, 7)) == Interval(2, 5)
    assert c1.intersect(Interval(5, 7)) == FiniteSet(5)
    assert c1.intersect(Interval(6, 9)) == EmptySet()

    # unevaluated object
    C1 = ComplexRegion(Interval(0, 1)*Interval(0, 2*S.Pi), polar=True)
    C2 = ComplexRegion(Interval(-1, 1)*Interval(-1, 1))
    assert C1.intersect(C2) == Intersection(C1, C2, evaluate=False)


def test_ComplexRegion_union():

    # Polar form
    c1 = ComplexRegion(Interval(0, 1)*Interval(0, 2*S.Pi), polar=True)
    c2 = ComplexRegion(Interval(0, 1)*Interval(0, S.Pi), polar=True)
    c3 = ComplexRegion(Interval(0, oo)*Interval(0, S.Pi), polar=True)
    c4 = ComplexRegion(Interval(0, oo)*Interval(S.Pi, 2*S.Pi), polar=True)

    p1 = Union(Interval(0, 1)*Interval(0, 2*S.Pi), Interval(0, 1)*Interval(0, S.Pi))
    p2 = Union(Interval(0, oo)*Interval(0, S.Pi), Interval(0, oo)*Interval(S.Pi, 2*S.Pi))

    assert c1.union(c2) == ComplexRegion(p1, polar=True)
    assert c3.union(c4) == ComplexRegion(p2, polar=True)

    # Rectangular form
    c5 = ComplexRegion(Interval(2, 5)*Interval(6, 9))
    c6 = ComplexRegion(Interval(4, 6)*Interval(10, 12))
    c7 = ComplexRegion(Interval(0, 10)*Interval(-10, 0))
    c8 = ComplexRegion(Interval(12, 16)*Interval(14, 20))

    p3 = Union(Interval(2, 5)*Interval(6, 9), Interval(4, 6)*Interval(10, 12))
    p4 = Union(Interval(0, 10)*Interval(-10, 0), Interval(12, 16)*Interval(14, 20))

    assert c5.union(c6) == ComplexRegion(p3)
    assert c7.union(c8) == ComplexRegion(p4)

    assert c1.union(Interval(2, 4)) == Union(c1, Interval(2, 4), evaluate=False)
    assert c5.union(Interval(2, 4)) == Union(c5, Interval(2, 4), evaluate=False)


def test_ComplexRegion_measure():
    a, b = Interval(2, 5), Interval(4, 8)
    theta1, theta2 = Interval(0, 2*S.Pi), Interval(0, S.Pi)
    c1 = ComplexRegion(a*b)
    c2 = ComplexRegion(Union(a*theta1, b*theta2), polar=True)

    assert c1.measure == 12
    assert c2.measure == 9*pi


def test_normalize_theta_set():

    # Interval
    assert normalize_theta_set(Interval(pi, 2*pi)) == \
        Union(FiniteSet(0), Interval(pi, 2*pi, False, True))
    assert normalize_theta_set(Interval(9*pi/2, 5*pi)) == Interval(pi/2, pi)
    assert normalize_theta_set(Interval(-3*pi/2, pi/2)) == \
        Interval(0, 2*pi, False, True)
    assert normalize_theta_set(Interval(-3*pi/2, pi/2, True, True)) == \
        Union(Interval(0, pi/2, False, True), Interval(pi/2, 2*pi, True, True))
    assert normalize_theta_set(Interval(-7*pi/2, -3*pi/2, True, True)) == \
        Union(Interval(0, pi/2, False, True), Interval(pi/2, 2*pi, True, True))
    assert normalize_theta_set(Interval(-pi/2, pi/2)) == \
        Union(Interval(0, pi/2), Interval(3*pi/2, 2*pi, False, True))
    assert normalize_theta_set(Interval(-pi/2, pi/2, True, True)) == \
        Union(Interval(0, pi/2, False, True), Interval(3*pi/2, 2*pi, True, True))
    assert normalize_theta_set(Interval(-4*pi, 3*pi)) == \
        Interval(0, 2*pi, False, True)
    assert normalize_theta_set(Interval(-3*pi/2, -pi/2)) == \
        Interval(pi/2, 3*pi/2)
    assert normalize_theta_set(Interval(0, 2*pi, True, True)) == \
        Interval(0, 2*pi, True, True)
    assert normalize_theta_set(Interval(-pi/2, pi/2, False, True)) == \
        Union(Interval(0, pi/2, False, True), Interval(3*pi/2, 2*pi, False, True))
    assert normalize_theta_set(Interval(-pi/2, pi/2, True, False)) == \
        Union(Interval(0, pi/2), Interval(3*pi/2, 2*pi, True, True))
    assert normalize_theta_set(Interval(-pi/2, pi/2, False, False)) == \
        Union(Interval(0, pi/2), Interval(3*pi/2, 2*pi, False, True))
    assert normalize_theta_set(Interval(4*pi, 9*pi/2, True, True)) == \
        Interval(0, pi/2, True, True)
    assert normalize_theta_set(Interval(4*pi, 9*pi/2, True, False)) == \
        Interval(0, pi/2, True, False)
    assert normalize_theta_set(Interval(4*pi, 9*pi/2,False, True)) == \
        Interval(0, pi/2, False, True)
    assert normalize_theta_set(Interval(3*pi, 5*pi, True, True)) == \
        Union(Interval(0, pi, False, True), Interval(pi, 2*pi, True, True))

    # FiniteSet
    assert normalize_theta_set(FiniteSet(0, pi, 3*pi)) == FiniteSet(0, pi)
    assert normalize_theta_set(FiniteSet(0, pi/2, pi, 2*pi)) == \
        FiniteSet(0, pi/2, pi)
    assert normalize_theta_set(FiniteSet(0, -pi/2, -pi, -2*pi)) == \
        FiniteSet(0, pi, 3*pi/2)
    assert normalize_theta_set(FiniteSet(-3*pi/2, pi/2)) == \
        FiniteSet(pi/2)
    assert normalize_theta_set(FiniteSet(2*pi)) == FiniteSet(0)

    # Unions
    assert normalize_theta_set(Union(Interval(0, pi/3), Interval(pi/2, pi))) == \
        Union(Interval(0, pi/3), Interval(pi/2, pi))
    assert normalize_theta_set(Union(Interval(0, pi), Interval(2*pi, 7*pi/3))) == \
        Interval(0, pi)

    # ValueError for non-real sets
    raises(ValueError, lambda: normalize_theta_set(S.Complexes))


def test_ComplexRegion_FiniteSet():
    x, y, z, a, b, c = symbols('x y z a b c')

    # Issue #9669
    assert ComplexRegion(FiniteSet(a, b, c)*FiniteSet(x, y, z)) == \
        FiniteSet(a + I*x, a + I*y, a + I*z, b + I*x, b + I*y,
                  b + I*z, c + I*x, c + I*y, c + I*z)
    assert ComplexRegion(FiniteSet(2)*FiniteSet(3)) == FiniteSet(2 + 3*I)


def test_union_RealSubSet():
    assert (S.Complexes).union(Interval(1, 2)) == S.Complexes
    assert (S.Complexes).union(S.Integers) == S.Complexes
