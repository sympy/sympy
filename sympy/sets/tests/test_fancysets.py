from sympy.core.compatibility import range
from sympy.sets.fancysets import ImageSet, Range
from sympy.sets.sets import FiniteSet, Interval, imageset, EmptySet
from sympy import (S, Symbol, Lambda, symbols, cos, sin, pi, oo, Basic,
        Rational, sqrt, tan, log, Abs)
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


def test_ImageSet_intersection():
    n_pi = imageset(x, x*pi, S.Integers)

    assert n_pi.intersect(Range(10)) == FiniteSet(0)
    assert n_pi.intersect(Interval(0, 10)) == \
        FiniteSet(0, pi, 2*pi, 3*pi)


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


def test_reals():
    assert 5 in S.Reals
    assert S.Pi in S.Reals
    assert -sqrt(2) in S.Reals
    assert (2, 5) not in S.Reals
    assert sqrt(-1) not in S.Reals


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
