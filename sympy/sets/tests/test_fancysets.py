from sympy.sets.fancysets import ImageSet, Range
from sympy.core.sets import FiniteSet, Interval, imageset
from sympy import (S, Symbol, Lambda, symbols, cos, sin, pi, oo, Basic,
        Rational, sqrt)
from sympy.utilities.pytest import XFAIL
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

    assert N.inf == 1
    assert N.sup == oo

def test_naturals0():
    N = S.Naturals0
    assert 0 in N
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


def test_ImageSet():
    squares = ImageSet(Lambda(x, x**2), S.Naturals)
    assert 4 in squares
    assert 5 not in squares
    assert FiniteSet(range(10)).intersect(squares) == FiniteSet(1, 4, 9)

    assert 16 not in squares.intersect(Interval(0, 10))

    si = iter(squares)
    a, b, c, d = next(si), next(si), next(si), next(si)
    assert (a, b, c, d) == (1, 4, 9, 16)

    harmonics = ImageSet(Lambda(x, 1/x), S.Naturals)
    assert Rational(1, 5) in harmonics
    assert .25 in harmonics
    assert .3 not in harmonics

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

    assert Range(5, 15).sup == 14
    assert Range(5, 15).inf == 5
    assert Range(15, 5, -1).sup == 15
    assert Range(15, 5, -1).inf == 6
    assert Range(10, 67, 10).sup == 60
    assert Range(60, 7, -10).inf == 10

    assert len(Range(10, 38, 10)) == 3

    assert Range(0, 0, 5) == S.EmptySet


def test_range_interval_intersection():
    # Intersection with intervals
    assert FiniteSet(Range(0, 10, 1).intersect(Interval(2, 6))) == \
        FiniteSet(2, 3, 4, 5, 6)

    # Open Intervals are removed
    assert (FiniteSet(Range(0, 10, 1).intersect(Interval(2, 6, True, True)))
            == FiniteSet(3, 4, 5))

    # Try this with large steps
    assert (FiniteSet(Range(0, 100, 10).intersect(Interval(15, 55))) ==
            FiniteSet(20, 30, 40, 50))

    # Going backwards
    assert FiniteSet(Range(10, -9, -3).intersect(Interval(-5, 6))) == \
        FiniteSet(-5, -2, 1, 4)
    assert FiniteSet(Range(10, -9, -3).intersect(Interval(-5, 6, True))) == \
        FiniteSet(-2, 1, 4)


def test_fun():
    assert (FiniteSet(ImageSet(Lambda(x, sin(pi*x/4)),
        Range(-10, 11))) == FiniteSet(-1, -sqrt(2)/2, 0, sqrt(2)/2, 1))


def test_reals():
    assert 5 in S.Reals
    assert S.Pi in S.Reals
    assert -sqrt(2) in S.Reals
    assert (2, 5) not in S.Reals


@XFAIL  # this is because contains is now very strict
def test_reals_fail():
    assert sqrt(-1) not in S.Reals


def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(itertools.islice(iterable, n))


def test_intersections():
    assert 5 in S.Integers.intersect(S.Reals)
    assert 5 in S.Integers.intersect(S.Reals)
    assert -5 not in S.Naturals.intersect(S.Reals)
    assert 5.5 not in S.Integers.intersect(S.Reals)
    assert 5 in S.Integers.intersect(Interval(3, oo))
    assert -5 in S.Integers.intersect(Interval(-oo, 3))
    assert all(x.is_Integer
            for x in take(10, S.Integers.intersect(Interval(3, oo)) ))
