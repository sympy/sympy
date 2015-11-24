from sympy import Symbol, S, exp, log, sqrt, oo, E, zoo
from sympy.calculus.util import not_empty_in, Limits
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


def test_Limits():
    assert Limits(1, 2).delta == S(1)
    assert Limits(1, 2).mid == S(3)/2

    assert Limits(1, 1) == S(1)

    assert Limits(1, 2) + 1 == Limits(2, 3)
    assert 1 + Limits(1, 2) == Limits(2, 3)
    assert Limits(1, 2) + Limits(2, 3) == Limits(3, 5)

    assert -Limits(1, 2) == Limits(-2, -1)

    assert Limits(1, 2) - 1 == Limits(0, 1)
    assert 1 - Limits(1, 2) == Limits(-1, 0)
    assert Limits(2, 3) - Limits(1, 2) == Limits(0, 2)

    assert x + Limits(1, 2) == Add(Limits(1, 2), x)
    assert a + Limits(1, 2) == Limits(1 + a, 2 + a)
    assert Limits(1, 2) - x == Add(Limits(1, 2), -x)

    assert Limits(-oo, 1) + oo == Limits(-oo, oo)
    assert Limits(1, oo) + oo == oo
    assert Limits(1, oo) - oo == Limits(-oo, oo)
    assert (-oo - Limits(-1, oo)) == -oo
    assert Limits(-oo, 1) - oo == -oo

    assert Limits(1, oo) - oo == Limits(-oo, oo)
    assert Limits(-oo, 1) - (-oo) == Limits(-oo, oo)
    assert (oo - Limits(1, oo)) == Limits(-oo, oo)
    assert (-oo - Limits(1, oo)) == -oo

    assert Limits(1, 2)/2 == Limits(S(1)/2, 1)
    assert 2/Limits(2, 3) == Limits(S(2)/3, 1)
    assert 1/Limits(-1, 1) == Limits(-oo, oo)

    assert abs(Limits(1, 2)) == Limits(1, 2)
    assert abs(Limits(-2, -1)) == Limits(1, 2)
    assert abs(Limits(-2, 1)) == Limits(0, 2)
    assert abs(Limits(-1, 2)) == Limits(0, 2)


def test_Limits_mul():
    assert Limits(1, 2)*2 == Limits(2, 4)
    assert 2*Limits(1, 2) == Limits(2, 4)
    assert Limits(1, 2)*Limits(2, 3) == Limits(2, 6)

    assert Limits(1, 2)*0 == 0
    assert Limits(1, oo)*0 == Limits(0, oo)
    assert Limits(-oo, 1)*0 == Limits(-oo, 0)
    assert Limits(-oo, oo)*0 == Limits(-oo, oo)

    assert Limits(1, 2)*x == Mul(Limits(1, 2), x, evaluate=False)

    assert Limits(0, 2)*oo == Limits(0, oo)
    assert Limits(-2, 0)*oo == Limits(-oo, 0)
    assert Limits(0, 2)*(-oo) == Limits(-oo, 0)
    assert Limits(-2, 0)*(-oo) == Limits(0, oo)
    assert Limits(-1, 1)*oo == Limits(-oo, oo)
    assert Limits(-1, 1)*(-oo) == Limits(-oo, oo)
    assert Limits(-oo, oo) == Limits(-oo, oo)


def test_Limits_div():
    assert Limits(-1, 3)/Limits(3, 4) == Limits(-S(1)/3, 1)
    assert Limits(-2, 4)/Limits(-3, 4) == Limits(-oo, oo)
    assert Limits(-3, -2)/Limits(-4, 0) == Limits(S(1)/2, oo)
    assert Limits(-3, -2)/Limits(-2, 1) == Union(Limits(-oo, -2), Limits(1, oo))
    assert Limits(-3, -2)/Limits(0, 4) == Limits(-oo, -S(1)/2)
    assert Limits(2, 4)/Limits(-3, 0) == Limits(-oo, -S(2)/3)
    assert Limits(2, 3)/Limits(-2, 2) == Union(Limits(-oo, -1), Limits(1, oo))
    assert Limits(2, 4)/Limits(0, 3) == Limits(S(2)/3, oo)

    assert Limits(0, 1)/Limits(0, 1) == Limits(0, oo)
    assert Limits(-1, 0)/Limits(0, 1) == Limits(-oo, 0)
    assert Limits(-1, 2)/Limits(-2, 2) == Limits(-oo, oo)

    assert 1/Limits(-1, 2) == Limits(-oo, oo)
    assert 1/Limits(0, 2) == Limits(S(1)/2, oo)
    assert (-1)/Limits(0, 2) == Limits(-oo, -S(1)/2)

    assert Limits(1, 2)/a == Mul(Limits(1, 2), 1/a, evaluate=False)

    assert Limits(1, 2)/0 == Limits(1, 2)*zoo
    assert Limits(1, oo)/oo == Limits(0, oo)
    assert Limits(1, oo)/(-oo) == Limits(-oo, 0)
    assert Limits(-oo, -1)/oo == Limits(-oo, 0)
    assert Limits(-oo, -1)/(-oo) == Limits(0, oo)
    assert Limits(-oo, oo)/oo == Limits(-oo, oo)
    assert Limits(-oo, oo)/(-oo) == Limits(-oo, oo)
    assert Limits(-1, oo)/oo == Limits(0, oo)
    assert Limits(-1, oo)/(-oo) == Limits(-oo, 0)
    assert Limits(-oo, 1)/oo == Limits(-oo, 0)
    assert Limits(-oo, 1)/(-oo) == Limits(0, oo)


def test_Limits_subs():
    assert (x**2 + 2*x + 1).subs(x, Limits(-1, 1)) == Limits(-1, 4)


def test_Limits_pow():
    assert Limits(-1, 1)**2 == Limits(0, 1)
    assert Limits(1, 2)**2 == Limits(1, 4)
    assert Limits(-1, 2)**3 == Limits(-1, 8)
    assert Limits(-1, 1)**0 == 1
    assert exp(Limits(0, 1)) == Limits(1, E)
    assert exp(Limits(-oo, oo)) == Limits(0, oo)
    assert log(Limits(3, 6)) == Limits(log(3), log(6))

    assert Limits(1, 2)**(S(5)/2) == Limits(1, 4*sqrt(2))
    assert Limits(-1, 2)**(S(1)/3) == Limits(-1, 2**(S(1)/3))
    assert Limits(0, 2)**(S(1)/2) == Limits(0, sqrt(2))
    assert Limits(0, 2)**(S(1)/2) == Limits(0, sqrt(2))

    assert Limits(-4, 2)**(S(2)/3) == Limits(0, 2*2**(S(1)/3))

    assert Limits(-1, 5)**(S(1)/2) == Limits(0, sqrt(5))
    assert Limits(-oo, 2)**(S(1)/2) == Limits(0, sqrt(2))
    assert Limits(-2, 3)**(S(-1)/4) == Limits(0, oo)

    assert Limits(1, 5)**(-2) == Limits(S(1)/25, 1)
    assert Limits(-1, 3)**(-2) == Limits(0, oo)
    assert Limits(0, 2)**(-2) == Limits(S(1)/4, oo)
    assert Limits(-1, 2)**(-3) == Limits(-oo, oo)
    assert Limits(-3, -2)**(-3) == Limits(S(-1)/8, -S(1)/27)
    assert Limits(-3, -2)**(-2) == Limits(S(1)/9, S(1)/4)
    assert Limits(0, oo)**(S(1)/2) == Limits(0, oo)
    assert Limits(-oo, -1)**(S(1)/3) == Limits(-oo, -1)
    assert Limits(-2, 3)**(-S(1)/3) == Limits(-oo, oo)
    assert Limits(-oo, 0)**(-2) == Limits(0, oo)
    assert Limits(-2, 0)**(-2) == Limits(S(1)/4, oo)

    assert Limits(1, 2)**x == Pow(Limits(1, 2), x, evaluate=False)
