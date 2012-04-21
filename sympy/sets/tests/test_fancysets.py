from sympy.sets.fancysets import TransformationSet, RangeSet
from sympy.core.sets import FiniteSet, Interval
from sympy import (S, Symbol, Lambda, symbols, cos, sin, pi, oo, Basic,
        Rational, sqrt)
from sympy.utilities.pytest import XFAIL

x = Symbol('x')

def test_naturals():
    N = S.Naturals
    assert 5 in N
    assert -5 not in N
    assert 5.5 not in N
    ni = iter(N)
    a,b,c,d = ni.next(), ni.next(), ni.next(), ni.next()
    assert (a,b,c,d) == (1,2,3,4)
    assert isinstance(a, Basic)

    assert N.intersect(Interval(-5, 5)) == FiniteSet(1, 2, 3, 4, 5)
    assert N.intersect(Interval(-5, 5, True, True)) == FiniteSet(1, 2, 3, 4)

    assert N.inf == 1
    assert N.sup == oo

def test_integers():
    Z = S.Integers
    assert 5 in Z
    assert -5 in Z
    assert 5.5 not in Z
    zi = iter(Z)
    a,b,c,d = zi.next(), zi.next(), zi.next(), zi.next()
    assert (a,b,c,d) == (0, 1, -1, 2)
    assert isinstance(a, Basic)

    assert Z.intersect(Interval(-5, 5)) == FiniteSet(range(-5, 6))
    assert Z.intersect(Interval(-5, 5, True, True)) == FiniteSet(range(-4,5))

    assert Z.inf == -oo
    assert Z.sup == oo

def test_TransformationSet():
    squares = TransformationSet(Lambda(x, x**2), S.Naturals)
    assert 4 in squares
    assert 5 not in squares
    assert FiniteSet(range(10)).intersect(squares) == FiniteSet(1, 4, 9)

    assert 16 not in squares.intersect(Interval(0, 10))

    si = iter(squares)
    a,b,c,d = si.next(), si.next(), si.next(), si.next()
    assert (a,b,c,d) == (1, 4, 9, 16)

    harmonics = TransformationSet(Lambda(x, 1/x), S.Naturals)
    assert Rational(1,5) in harmonics
    assert .25 in harmonics
    assert .3 not in harmonics

    assert harmonics.is_iterable

@XFAIL
def test_halfcircle():
    # This test sometimes works and sometimes doesn't.
    # It may be an issue with solve? Maybe with using Lambdas/dummys?
    # I believe the code within fancysets is correct
    r, th = symbols('r, theta', real=True)
    L = Lambda((r, th), (r*cos(th), r*sin(th)))
    halfcircle = TransformationSet(L, Interval(0, 1)*Interval(0, pi))

    assert (1, 0) in halfcircle
    assert (0, -1) not in halfcircle
    assert (0, 0) in halfcircle

    assert not halfcircle.is_iterable

def test_transformation_iterator_not_injetive():
    L = Lambda(x, x - x%2) # produces 0, 2, 2, 4, 4, 6, 6, ...
    evens = TransformationSet(L, S.Naturals)
    i = iter(evens)
    # No repeats here
    assert (i.next(), i.next(), i.next(), i.next()) == (0, 2, 4, 6)

def test_RangeSet():
    assert RangeSet(5) == RangeSet(0, 5) == RangeSet(0, 5, 1)

    r = RangeSet(10, 20, 2)
    assert 12 in r
    assert 8 not in r
    assert 11 not in r
    assert 30 not in r

    assert list(RangeSet(0, 5)) == range(5)
    assert list(RangeSet(5, 0, -1)) == range(1, 6)

    assert RangeSet(5, 15).sup == 14
    assert RangeSet(5, 15).inf == 5
    assert RangeSet(15, 5, -1).sup == 15
    assert RangeSet(15, 5, -1).inf == 6
    assert RangeSet(10, 67, 10).sup == 60
    assert RangeSet(60, 7, -10).inf == 10

    assert len(RangeSet(10, 38, 10)) == 3

    assert RangeSet(0, 0, 5) == S.EmptySet


def test_range_interval_intersection():
    # Intersection with intervals
    assert FiniteSet(RangeSet(0, 10, 1).intersect(Interval(2, 6))) == \
            FiniteSet(2, 3, 4, 5, 6)

    # Open Intervals are removed
    assert (FiniteSet(RangeSet(0, 10, 1).intersect(Interval(2, 6, True, True)))
            == FiniteSet(3, 4, 5))

    # Try this with large steps
    assert (FiniteSet(RangeSet(0, 100, 10).intersect(Interval(15, 55))) ==
            FiniteSet(20, 30, 40, 50))

    # Going backwards
    assert FiniteSet(RangeSet(10, -9, -3).intersect(Interval(-5, 6))) == \
            FiniteSet(-5, -2, 1, 4)
    assert FiniteSet(RangeSet(10, -9, -3).intersect(Interval(-5, 6, True))) == \
            FiniteSet(-2, 1, 4)

def test_fun():
    assert (FiniteSet(TransformationSet(Lambda(x, sin(pi*x/4)),
        RangeSet(-10,11))) == FiniteSet(-1, -sqrt(2)/2, 0, sqrt(2)/2, 1))
