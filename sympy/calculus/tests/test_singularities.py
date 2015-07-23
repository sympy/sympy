from sympy import Symbol, exp, log, S, oo
from sympy.calculus.singularities import singularities, range_func
from sympy.sets.sets import Interval, FiniteSet, Complement, Union

from sympy.utilities.pytest import XFAIL


def test_singularities():
    x = Symbol('x', real=True)

    assert singularities(x**2, x) == ()
    assert singularities(x/(x**2 + 3*x + 2), x) == (-2, -1)


@XFAIL
def test_singularities_non_rational():
    x = Symbol('x', real=True)

    assert singularities(exp(1/x), x) == (0)
    assert singularities(log((x - 2)**2), x) == (2)


def test_range_func():
    x = Symbol('x', real=True)
    assert range_func(x, Interval(-1, 1)) == Interval(-1, 1)
    assert range_func(x, Interval(0, 1, True, True)) == Interval(0, 1, True, True)
    assert range_func(x, Interval(1, 2, True, False)) == Interval(1, 2, True, False)
    assert range_func(x**2, Interval(1, 5)) == Interval(1, 25)
    assert range_func(x**3, Interval(0, 1)) == Interval(0, 1)
    assert range_func(x/(x - 1), Interval(2, 3)) == Interval(S(3)/2, 2)
    assert range_func(x/(x**2 - 4), Interval(3, 4)) == Interval(S(1)/3, S(3)/5)
    assert range_func(1/x**2, FiniteSet(1, 2, -1, 0)) == FiniteSet(1, S(1)/4)
    assert range_func(1, Interval(-1, 4)) == FiniteSet(1)
    assert range_func(x, Interval(1, 2, False, True)) == Interval(1, 2, False, True)
    assert range_func(x**2, Interval(-1, 1)) == Interval(0, 1)
    assert range_func(x**2, Interval(-1, 1, True, True)) == Interval(0, 1, False, True)
    assert range_func(x**2, Interval(-1, 1, False, True)) == Interval(0, 1)
    assert range_func(x**2, Interval(-1, 1, True, False)) == Interval(0, 1)
    assert range_func(x, Interval(-1, 1)) == Interval(-1, 1)
    assert range_func(x, FiniteSet(1, -1, 3, 5)) == FiniteSet(-1, 1, 3, 5)
    assert range_func(x**2 - x, FiniteSet(1, -1, 3, 5, -oo)) == FiniteSet(0, 2, 6, 20, oo)
    assert range_func(1/x, Interval(0, 1)) == Interval(1, oo)
    assert range_func(1/x, Interval(-1, 1)) == Union(Interval(-oo, -1), Interval(1, oo))
    assert range_func(x**2 - x, FiniteSet(1, -1, 3, 5, -oo)) == FiniteSet(0, 2, 6, 20, oo)
    assert range_func(1/x**2, Interval(-1, 1)) == Interval(1, oo)
    assert range_func(1/x**2, Interval(-1, 1, True, False)) == Interval(1, oo)
    assert range_func(1/x**2, Interval(-1, 1, True, True)) == Interval(1, oo, True, True)
    assert range_func(1/x**2, Interval(-1, 1, False, True)) == Interval(1, oo)
    assert range_func(1/(x - 4), Interval(0, 5)) == Union(Interval(-oo, S(-1)/4), Interval(1, oo))
    assert range_func(1/x, Interval(1, 2, False, True)) == Interval(S(1)/2, 1, True, False)
    assert range_func(1/x, Interval(1, 2, True, False)) == Interval(S(1)/2, 1, False, True)
    assert range_func(1/x, Interval(1, 2)) == Interval(S(1)/2, 1)
    assert range_func(1/x**2, Interval(-2, -1, True, True)) == Interval(S(1)/4, 1, True, True)
    assert range_func(x, Interval(-oo, oo)) == Interval(-oo, oo)
    assert range_func(x**2/(x - 4), Interval(-oo, oo)) == Complement(S.Reals, Interval(0, 16, True, True))
    assert range_func(x**2/(x - 4), FiniteSet(4)) == S.EmptySet
    assert range_func(x**2/(x - 4), Interval(3, 4)) == Interval(-oo, -9)
    assert range_func(-x**2/(x - 4), Interval(3, 4)) == Interval(9, oo)
