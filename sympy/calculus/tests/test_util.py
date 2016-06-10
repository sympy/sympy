from sympy import (Symbol, S, exp, log, sqrt, oo, E, zoo, tan,
        sin, pi)
from sympy.calculus.util import not_empty_in, AccumBounds
from sympy.core import Add, Mul, Pow
from sympy.sets.sets import Interval, FiniteSet, Complement, Union
from sympy.utilities.pytest import raises
from sympy.abc import x

a = Symbol('a', real=True)


def test_not_empty_in():
    assert not_empty_in(FiniteSet(x, 2*x).intersect(Interval(1, 2, True, False)), x) == \
        Interval(S(1)/2, 2, True, False)
    assert not_empty_in(FiniteSet(x, x**2).intersect(Interval(1, 2)), x) == \
        Union(Interval(-sqrt(2), -1), Interval(1, 2))
    assert not_empty_in(FiniteSet(x**2 + x, x).intersect(Interval(2, 4)), x) == \
        Union(Interval(-sqrt(17)/2 - S(1)/2, -2),
              Interval(1, -S(1)/2 + sqrt(17)/2), Interval(2, 4))
    assert not_empty_in(FiniteSet(x/(x - 1)).intersect(S.Reals), x) == \
        Complement(S.Reals, FiniteSet(1))
    assert not_empty_in(FiniteSet(a/(a - 1)).intersect(S.Reals), a) == \
        Complement(S.Reals, FiniteSet(1))
    assert not_empty_in(FiniteSet((x**2 - 3*x + 2)/(x - 1)).intersect(S.Reals), x) == \
        Complement(S.Reals, FiniteSet(1))
    assert not_empty_in(FiniteSet(3, 4, x/(x - 1)).intersect(Interval(2, 3)), x) == \
        Union(Interval(S(3)/2, 2), FiniteSet(3))
    assert not_empty_in(FiniteSet(x/(x**2 - 1)).intersect(S.Reals), x) == \
        Complement(S.Reals, FiniteSet(-1, 1))
    assert not_empty_in(FiniteSet(x, x**2).intersect(Union(Interval(1, 3, True, True),
                                                           Interval(4, 5))), x) == \
        Union(Interval(-sqrt(5), -2), Interval(-sqrt(3), -1, True, True),
              Interval(1, 3, True, True), Interval(4, 5))
    assert not_empty_in(FiniteSet(1).intersect(Interval(3, 4)), x) == S.EmptySet
    assert not_empty_in(FiniteSet(x**2/(x + 2)).intersect(Interval(1, oo)), x) == \
        Union(Interval(-2, -1, True, False), Interval(2, oo))


def test_AccumBounds():
    assert AccumBounds(1, 2).args == (1, 2)
    assert AccumBounds(1, 2).delta == S(1)
    assert AccumBounds(1, 2).mid == S(3)/2
    assert AccumBounds(1, 3).is_real == True

    assert AccumBounds(1, 1) == S(1)

    assert AccumBounds(1, 2) + 1 == AccumBounds(2, 3)
    assert 1 + AccumBounds(1, 2) == AccumBounds(2, 3)
    assert AccumBounds(1, 2) + AccumBounds(2, 3) == AccumBounds(3, 5)

    assert -AccumBounds(1, 2) == AccumBounds(-2, -1)

    assert AccumBounds(1, 2) - 1 == AccumBounds(0, 1)
    assert 1 - AccumBounds(1, 2) == AccumBounds(-1, 0)
    assert AccumBounds(2, 3) - AccumBounds(1, 2) == AccumBounds(0, 2)

    assert x + AccumBounds(1, 2) == Add(AccumBounds(1, 2), x)
    assert a + AccumBounds(1, 2) == AccumBounds(1 + a, 2 + a)
    assert AccumBounds(1, 2) - x == Add(AccumBounds(1, 2), -x)

    assert AccumBounds(-oo, 1) + oo == AccumBounds(-oo, oo)
    assert AccumBounds(1, oo) + oo == oo
    assert AccumBounds(1, oo) - oo == AccumBounds(-oo, oo)
    assert (-oo - AccumBounds(-1, oo)) == -oo
    assert AccumBounds(-oo, 1) - oo == -oo

    assert AccumBounds(1, oo) - oo == AccumBounds(-oo, oo)
    assert AccumBounds(-oo, 1) - (-oo) == AccumBounds(-oo, oo)
    assert (oo - AccumBounds(1, oo)) == AccumBounds(-oo, oo)
    assert (-oo - AccumBounds(1, oo)) == -oo

    assert AccumBounds(1, 2)/2 == AccumBounds(S(1)/2, 1)
    assert 2/AccumBounds(2, 3) == AccumBounds(S(2)/3, 1)
    assert 1/AccumBounds(-1, 1) == AccumBounds(-oo, oo)

    assert abs(AccumBounds(1, 2)) == AccumBounds(1, 2)
    assert abs(AccumBounds(-2, -1)) == AccumBounds(1, 2)
    assert abs(AccumBounds(-2, 1)) == AccumBounds(0, 2)
    assert abs(AccumBounds(-1, 2)) == AccumBounds(0, 2)


def test_AccumBounds_mul():
    assert AccumBounds(1, 2)*2 == AccumBounds(2, 4)
    assert 2*AccumBounds(1, 2) == AccumBounds(2, 4)
    assert AccumBounds(1, 2)*AccumBounds(2, 3) == AccumBounds(2, 6)

    assert AccumBounds(1, 2)*0 == 0
    assert AccumBounds(1, oo)*0 == AccumBounds(0, oo)
    assert AccumBounds(-oo, 1)*0 == AccumBounds(-oo, 0)
    assert AccumBounds(-oo, oo)*0 == AccumBounds(-oo, oo)

    assert AccumBounds(1, 2)*x == Mul(AccumBounds(1, 2), x, evaluate=False)

    assert AccumBounds(0, 2)*oo == AccumBounds(0, oo)
    assert AccumBounds(-2, 0)*oo == AccumBounds(-oo, 0)
    assert AccumBounds(0, 2)*(-oo) == AccumBounds(-oo, 0)
    assert AccumBounds(-2, 0)*(-oo) == AccumBounds(0, oo)
    assert AccumBounds(-1, 1)*oo == AccumBounds(-oo, oo)
    assert AccumBounds(-1, 1)*(-oo) == AccumBounds(-oo, oo)
    assert AccumBounds(-oo, oo)*oo == AccumBounds(-oo, oo)


def test_AccumBounds_div():
    assert AccumBounds(-1, 3)/AccumBounds(3, 4) == AccumBounds(-S(1)/3, 1)
    assert AccumBounds(-2, 4)/AccumBounds(-3, 4) == AccumBounds(-oo, oo)
    assert AccumBounds(-3, -2)/AccumBounds(-4, 0) == AccumBounds(S(1)/2, oo)

    # these two tests can have a better answer
    # after Union of AccumBounds is improved
    assert AccumBounds(-3, -2)/AccumBounds(-2, 1) == AccumBounds(-oo, oo)
    assert AccumBounds(2, 3)/AccumBounds(-2, 2) == AccumBounds(-oo, oo)

    assert AccumBounds(-3, -2)/AccumBounds(0, 4) == AccumBounds(-oo, -S(1)/2)
    assert AccumBounds(2, 4)/AccumBounds(-3, 0) == AccumBounds(-oo, -S(2)/3)
    assert AccumBounds(2, 4)/AccumBounds(0, 3) == AccumBounds(S(2)/3, oo)

    assert AccumBounds(0, 1)/AccumBounds(0, 1) == AccumBounds(0, oo)
    assert AccumBounds(-1, 0)/AccumBounds(0, 1) == AccumBounds(-oo, 0)
    assert AccumBounds(-1, 2)/AccumBounds(-2, 2) == AccumBounds(-oo, oo)

    assert 1/AccumBounds(-1, 2) == AccumBounds(-oo, oo)
    assert 1/AccumBounds(0, 2) == AccumBounds(S(1)/2, oo)
    assert (-1)/AccumBounds(0, 2) == AccumBounds(-oo, -S(1)/2)
    assert 1/AccumBounds(-oo, 0) == AccumBounds(-oo, 0)
    assert 1/AccumBounds(-1, 0) == AccumBounds(-oo, -1)
    assert (-2)/AccumBounds(-oo, 0) == AccumBounds(0, oo)
    assert 1/AccumBounds(-oo, -1) == AccumBounds(-1, 0)

    assert AccumBounds(1, 2)/a == Mul(AccumBounds(1, 2), 1/a, evaluate=False)

    assert AccumBounds(1, 2)/0 == AccumBounds(1, 2)*zoo
    assert AccumBounds(1, oo)/oo == AccumBounds(0, oo)
    assert AccumBounds(1, oo)/(-oo) == AccumBounds(-oo, 0)
    assert AccumBounds(-oo, -1)/oo == AccumBounds(-oo, 0)
    assert AccumBounds(-oo, -1)/(-oo) == AccumBounds(0, oo)
    assert AccumBounds(-oo, oo)/oo == AccumBounds(-oo, oo)
    assert AccumBounds(-oo, oo)/(-oo) == AccumBounds(-oo, oo)
    assert AccumBounds(-1, oo)/oo == AccumBounds(0, oo)
    assert AccumBounds(-1, oo)/(-oo) == AccumBounds(-oo, 0)
    assert AccumBounds(-oo, 1)/oo == AccumBounds(-oo, 0)
    assert AccumBounds(-oo, 1)/(-oo) == AccumBounds(0, oo)


def test_AccumBounds_func():
    assert (x**2 + 2*x + 1).subs(x, AccumBounds(-1, 1)) == AccumBounds(-1, 4)
    assert exp(AccumBounds(0, 1)) == AccumBounds(1, E)
    assert exp(AccumBounds(-oo, oo)) == AccumBounds(0, oo)
    assert log(AccumBounds(3, 6)) == AccumBounds(log(3), log(6))


def test_AccumBounds_pow():
    assert AccumBounds(0, 2)**2 == AccumBounds(0, 4)
    assert AccumBounds(-1, 1)**2 == AccumBounds(0, 1)
    assert AccumBounds(1, 2)**2 == AccumBounds(1, 4)
    assert AccumBounds(-1, 2)**3 == AccumBounds(-1, 8)
    assert AccumBounds(-1, 1)**0 == 1

    assert AccumBounds(1, 2)**(S(5)/2) == AccumBounds(1, 4*sqrt(2))
    assert AccumBounds(-1, 2)**(S(1)/3) == AccumBounds(-1, 2**(S(1)/3))
    assert AccumBounds(0, 2)**(S(1)/2) == AccumBounds(0, sqrt(2))

    assert AccumBounds(-4, 2)**(S(2)/3) == AccumBounds(0, 2*2**(S(1)/3))

    assert AccumBounds(-1, 5)**(S(1)/2) == AccumBounds(0, sqrt(5))
    assert AccumBounds(-oo, 2)**(S(1)/2) == AccumBounds(0, sqrt(2))
    assert AccumBounds(-2, 3)**(S(-1)/4) == AccumBounds(0, oo)

    assert AccumBounds(1, 5)**(-2) == AccumBounds(S(1)/25, 1)
    assert AccumBounds(-1, 3)**(-2) == AccumBounds(0, oo)
    assert AccumBounds(0, 2)**(-2) == AccumBounds(S(1)/4, oo)
    assert AccumBounds(-1, 2)**(-3) == AccumBounds(-oo, oo)
    assert AccumBounds(-3, -2)**(-3) == AccumBounds(S(-1)/8, -S(1)/27)
    assert AccumBounds(-3, -2)**(-2) == AccumBounds(S(1)/9, S(1)/4)
    assert AccumBounds(0, oo)**(S(1)/2) == AccumBounds(0, oo)
    assert AccumBounds(-oo, -1)**(S(1)/3) == AccumBounds(-oo, -1)
    assert AccumBounds(-2, 3)**(-S(1)/3) == AccumBounds(-oo, oo)
    assert AccumBounds(-oo, 0)**(-2) == AccumBounds(0, oo)
    assert AccumBounds(-2, 0)**(-2) == AccumBounds(S(1)/4, oo)

    assert AccumBounds(S(1)/3, S(1)/2)**oo == S(0)
    assert AccumBounds(0, S(1)/2)**oo == S(0)
    assert AccumBounds(S(1)/2, 1)**oo == AccumBounds(0, oo)
    assert AccumBounds(0, 1)**oo == AccumBounds(0, oo)
    assert AccumBounds(2, 3)**oo == oo
    assert AccumBounds(1, 2)**oo == AccumBounds(0, oo)
    assert AccumBounds(S(1)/2, 3)**oo == AccumBounds(0, oo)
    assert AccumBounds(-S(1)/3, -S(1)/4)**oo == S(0)
    assert AccumBounds(-1, -S(1)/2)**oo == AccumBounds(-oo, oo)
    assert AccumBounds(-3, -2)**oo == FiniteSet(-oo, oo)
    assert AccumBounds(-2, -1)**oo == AccumBounds(-oo, oo)
    assert AccumBounds(-2, -S(1)/2)**oo == AccumBounds(-oo, oo)
    assert AccumBounds(-S(1)/2, S(1)/2)**oo == S(0)
    assert AccumBounds(-S(1)/2, 1)**oo == AccumBounds(0, oo)
    assert AccumBounds(-S(2)/3, 2)**oo == AccumBounds(0, oo)
    assert AccumBounds(-1, 1)**oo == AccumBounds(-oo, oo)
    assert AccumBounds(-1, S(1)/2)**oo == AccumBounds(-oo, oo)
    assert AccumBounds(-1, 2)**oo == AccumBounds(-oo, oo)
    assert AccumBounds(-2, S(1)/2)**oo == AccumBounds(-oo, oo)

    assert AccumBounds(1, 2)**x == Pow(AccumBounds(1, 2), x, evaluate=False)

    assert AccumBounds(2, 3)**(-oo) == S(0)
    assert AccumBounds(0, 2)**(-oo) == AccumBounds(0, oo)
    assert AccumBounds(-1, 2)**(-oo) == AccumBounds(-oo, oo)

    assert (tan(x)**sin(2*x)).subs(x, AccumBounds(0, pi/2)) == \
        Pow(AccumBounds(-oo, oo), AccumBounds(0, 1), evaluate=False)


def test_comparison_AccumBounds():
    assert (AccumBounds(1, 3) < 4) == S.true
    assert (AccumBounds(1, 3) < -1) == S.false
    assert (AccumBounds(1, 3) < 2) is None

    assert (AccumBounds(1, 3) > 4) == S.false
    assert (AccumBounds(1, 3) > -1) == S.true
    assert (AccumBounds(1, 3) > 2) is None

    assert (AccumBounds(1, 3) < AccumBounds(4, 6)) == S.true
    assert (AccumBounds(1, 3) < AccumBounds(2, 4)) is None
    assert (AccumBounds(1, 3) < AccumBounds(-2, 0)) == S.false


def test_contains_AccumBounds():
    assert (1 in AccumBounds(1, 2)) == S.true
    raises(TypeError, lambda: a in AccumBounds(1, 2))
    assert (-oo in AccumBounds(1, oo)) == S.true
    assert (oo in AccumBounds(-oo, 0)) == S.true
