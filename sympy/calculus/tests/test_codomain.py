from sympy import Symbol, S, oo, sqrt
from sympy.calculus.codomain import codomain, not_empty_in
from sympy.sets.sets import Interval, FiniteSet, Complement, Union
from sympy.utilities.pytest import XFAIL, raises


def test_codomain():
    x = Symbol('x', real=True)
    assert codomain(x, Interval(-1, 1), x) == Interval(-1, 1)
    assert codomain(x, Interval(0, 1, True, True), x) == \
            Interval(0, 1, True, True)
    assert codomain(x, Interval(1, 2, True, False), x) == Interval(1, 2, True, False)
    assert codomain(x, Interval(1, 2, False, True), x) == Interval(1, 2, False, True)
    assert codomain(x**2, Interval(-1, 1), x) == Interval(0, 1)
    assert codomain(x**3, Interval(0, 1), x) == Interval(0, 1)
    assert codomain(x/(x**2 - 4), Interval(3, 4), x) == Interval(S(1)/3, S(3)/5)
    assert codomain(1, Interval(-1, 4), x) == FiniteSet(1)
    assert codomain(x, Interval(-oo, oo), x) == S.Reals

    assert codomain(1/x**2, FiniteSet(1, 2, -1, 0), x) == FiniteSet(1, S(1)/4)
    assert codomain(x, FiniteSet(1, -1, 3, 5), x) == FiniteSet(-1, 1, 3, 5)
    assert codomain(x**2 - x, FiniteSet(1, -1, 3, 5, -oo), x) == \
            FiniteSet(0, 2, 6, 20, oo)
    assert codomain(x**2/(x - 4), FiniteSet(4), x) == S.EmptySet
    assert codomain(x**2 - x, FiniteSet(S(1)/2, -oo, oo, 2), x) == \
            FiniteSet(S(-1)/4, 2, oo)

    assert codomain(x**2, Interval(-1, 1, True, True), x) == Interval(0, 1, False, True)
    assert codomain(x**2, Interval(-1, 1, False, True), x) == Interval(0, 1)
    assert codomain(x**2, Interval(-1, 1, True, False), x) == Interval(0, 1)

    assert codomain(1/x, Interval(0, 1), x) == Interval(1, oo)
    assert codomain(1/x, Interval(-1, 1), x) == Union(Interval(-oo, -1), Interval(1, oo))
    assert codomain(1/x**2, Interval(-1, 1), x) == Interval(1, oo)
    assert codomain(1/x**2, Interval(-1, 1, True, False), x) == Interval(1, oo)
    assert codomain(1/x**2, Interval(-1, 1, True, True), x) == \
            Interval(1, oo, True, True)
    assert codomain(1/x**2, Interval(-1, 1, False, True), x) == Interval(1, oo)
    assert codomain(1/x, Interval(1, 2), x) == Interval(S(1)/2, 1)
    assert codomain(1/x**2, Interval(-2, -1, True, True), x) == \
            Interval(S(1)/4, 1, True, True)
    assert codomain(x**2/(x - 4), Interval(-oo, oo), x) == \
            Complement(S.Reals, Interval(0, 16, True, True))
    assert codomain(x**2/(x - 4), Interval(3, 4), x) == Interval(-oo, -9)
    assert codomain(-x**2/(x - 4), Interval(3, 4), x) == Interval(9, oo)
    assert codomain((x**2 - x)/(x**3 - 1), S.Reals, x) == Interval(-1, S(1)/3, False, True)
    assert codomain(-x**2 + 1/x, S.Reals, x) == S.Reals
    assert codomain(x**2 - 1/x, S.Reals, x) == S.Reals

    assert codomain(x**2, Union(Interval(1, 2), FiniteSet(3)), x) == \
            Union(Interval(1, 4), FiniteSet(9))
    assert codomain(x/(x**2 - 4), Union(Interval(-oo, 1), Interval(0, oo)), x) == S.Reals
    assert codomain(x, Union(Interval(-1, 1), FiniteSet(-oo)), x) == \
            Union(Interval(-1, 1), FiniteSet(-oo))
    assert codomain(x**2 - x, Interval(1, oo), x) == Interval(0, oo)
    raises(ValueError, lambda: codomain(sqrt(x), Interval(-1, 2), x))


def test_not_empty_in():
    from sympy.abc import x
    a = Symbol('a', real=True)
    assert not_empty_in(FiniteSet(x, 2*x).intersect(Interval(1, 2)), x) == Interval(S(1)/2, 2)
    assert not_empty_in(FiniteSet(x, x**2).intersect(Interval(1, 2)), x) == \
        Union(Interval(-sqrt(2), -1), Interval(1, 2))
    assert not_empty_in(FiniteSet(x**2 + x, x).intersect(Interval(2, 4)), x) == \
        Union(Interval(-sqrt(17)/2 - S(1)/2, -2), Interval(1, -S(1)/2 + sqrt(17)/2), Interval(2, 4))
    assert not_empty_in(FiniteSet(x/(x - 1)).intersect(S.Reals), x) == Complement(S.Reals, FiniteSet(1))
    assert not_empty_in(FiniteSet(a/(a - 1)).intersect(S.Reals), a) == Complement(S.Reals, FiniteSet(1))
    assert not_empty_in(FiniteSet((x**2 - 3*x + 2)/(x - 1)).intersect(S.Reals), x) == \
        Complement(S.Reals, FiniteSet(1))
    assert not_empty_in(FiniteSet(3, 4, x/(x - 1)).intersect(Interval(2, 3)), x) == \
        Union(Interval(S(3)/2, 2), FiniteSet(3))


@XFAIL
def test_failing_not_empty_in():
    assert not_empty_in(FiniteSet(x/(x**2 - 1)).intersect(S.Reals), x) == \
        Complement(S.Reals, FiniteSet(-1, 1))
