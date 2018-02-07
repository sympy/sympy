from sympy.sets.setexpr import SetExpr
from sympy.utilities.pytest import XFAIL
from sympy.sets import Interval, FiniteSet
from sympy import Expr, Set, exp, log, sin, cos, Symbol, Min, Max, S, oo

I = Interval(0, 2)


def test_setexpr():
    se = SetExpr(Interval(0, 1))
    assert isinstance(se.set, Set)
    assert isinstance(se, Expr)


def test_scalar_funcs():
    assert SetExpr(Interval(0, 1)).set == Interval(0, 1)
    a, b = Symbol('a', real=True), Symbol('b', real=True)
    a, b = 1, 2
    # TODO: add support for more functions in the future:
    for f in [exp, log]:
        input_se = f(SetExpr(Interval(a, b)))
        output = input_se.set
        expected = Interval(Min(f(a), f(b)), Max(f(a), f(b)))
        assert output == expected


def test_Add_Mul():
    assert (SetExpr(Interval(0, 1)) + 1).set == Interval(1, 2)
    assert (SetExpr(Interval(0, 1)) * 2).set == Interval(0, 2)


def test_Pow():
    assert (SetExpr(Interval(0, 2))**2).set == Interval(0, 4)


def test_compound():
    assert (exp(SetExpr(Interval(0, 1)) * 2 + 1)).set == \
        Interval(exp(1), exp(3))


def test_Interval_Interval():
    assert (SetExpr(Interval(1, 2)) + SetExpr(Interval(10, 20))).set == \
        Interval(11, 22)
    assert (SetExpr(Interval(1, 2)) * SetExpr(Interval(10, 20))).set == \
        Interval(10, 40)


def test_FiniteSet_FiniteSet():
    assert (SetExpr(FiniteSet(1, 2, 3)) + SetExpr(FiniteSet(1, 2))).set ==\
        FiniteSet(2, 3, 4, 5)
    assert (SetExpr(FiniteSet(1, 2, 3)) * SetExpr(FiniteSet(1, 2))).set ==\
        FiniteSet(1, 2, 3, 4, 6)


def test_Interval_FiniteSet():
    assert (SetExpr(FiniteSet(1, 2)) + SetExpr(Interval(0, 10))).set == \
        Interval(1, 12)


def test_Many_Sets():
    assert (SetExpr(Interval(0, 1)) +
                    SetExpr(Interval(2, 3)) +
                    SetExpr(FiniteSet(10, 11, 12))).set == Interval(12, 16)


def test_same_setexprs_are_not_identical():
    a = SetExpr(FiniteSet(0, 1))
    b = SetExpr(FiniteSet(0, 1))
    assert (a + b).set == FiniteSet(0, 1, 2)

    # Cannont detect the set being the same:
    # assert (a + a).set == FiniteSet(0, 2)


def test_Interval_arithmetic():
    i12cc = SetExpr(Interval(1, 2))
    i12lo = SetExpr(Interval.Lopen(1, 2))
    i12ro = SetExpr(Interval.Ropen(1, 2))
    i12o = SetExpr(Interval.open(1, 2))

    n23cc = SetExpr(Interval(-2, 3))
    n23lo = SetExpr(Interval.Lopen(-2, 3))
    n23ro = SetExpr(Interval.Ropen(-2, 3))
    n23o = SetExpr(Interval.open(-2, 3))

    n3n2cc = SetExpr(Interval(-3, -2))

    assert i12cc + i12cc == SetExpr(Interval(2, 4))
    assert i12cc - i12cc == SetExpr(Interval(-1, 1))
    assert i12cc * i12cc == SetExpr(Interval(1, 4))
    assert i12cc / i12cc == SetExpr(Interval(S.Half, 2))
    assert i12cc ** 2 == SetExpr(Interval(1, 4))
    assert i12cc ** 3 == SetExpr(Interval(1, 8))

    assert i12lo + i12ro == SetExpr(Interval.open(2, 4))
    assert i12lo - i12ro == SetExpr(Interval.Lopen(-1, 1))
    assert i12lo * i12ro == SetExpr(Interval.open(1, 4))
    assert i12lo / i12ro == SetExpr(Interval.Lopen(S.Half, 2))
    assert i12lo + i12lo == SetExpr(Interval.Lopen(2, 4))
    assert i12lo - i12lo == SetExpr(Interval.open(-1, 1))
    assert i12lo * i12lo == SetExpr(Interval.Lopen(1, 4))
    assert i12lo / i12lo == SetExpr(Interval.open(S.Half, 2))
    assert i12lo + i12cc == SetExpr(Interval.Lopen(2, 4))
    assert i12lo - i12cc == SetExpr(Interval.Lopen(-1, 1))
    assert i12lo * i12cc == SetExpr(Interval.Lopen(1, 4))
    assert i12lo / i12cc == SetExpr(Interval.Lopen(S.Half, 2))
    assert i12lo + i12o == SetExpr(Interval.open(2, 4))
    assert i12lo - i12o == SetExpr(Interval.open(-1, 1))
    assert i12lo * i12o == SetExpr(Interval.open(1, 4))
    assert i12lo / i12o == SetExpr(Interval.open(S.Half, 2))
    assert i12lo ** 2 == SetExpr(Interval.Lopen(1, 4))
    assert i12lo ** 3 == SetExpr(Interval.Lopen(1, 8))

    assert i12ro + i12ro == SetExpr(Interval.Ropen(2, 4))
    assert i12ro - i12ro == SetExpr(Interval.open(-1, 1))
    assert i12ro * i12ro == SetExpr(Interval.Ropen(1, 4))
    assert i12ro / i12ro == SetExpr(Interval.open(S.Half, 2))
    assert i12ro + i12cc == SetExpr(Interval.Ropen(2, 4))
    assert i12ro - i12cc == SetExpr(Interval.Ropen(-1, 1))
    assert i12ro * i12cc == SetExpr(Interval.Ropen(1, 4))
    assert i12ro / i12cc == SetExpr(Interval.Ropen(S.Half, 2))
    assert i12ro + i12o == SetExpr(Interval.open(2, 4))
    assert i12ro - i12o == SetExpr(Interval.open(-1, 1))
    assert i12ro * i12o == SetExpr(Interval.open(1, 4))
    assert i12ro / i12o == SetExpr(Interval.open(S.Half, 2))
    assert i12ro ** 2 == SetExpr(Interval.Ropen(1, 4))
    assert i12ro ** 3 == SetExpr(Interval.Ropen(1, 8))

    assert i12o + i12lo == SetExpr(Interval.open(2, 4))
    assert i12o - i12lo == SetExpr(Interval.open(-1, 1))
    assert i12o * i12lo == SetExpr(Interval.open(1, 4))
    assert i12o / i12lo == SetExpr(Interval.open(S.Half, 2))
    assert i12o + i12ro == SetExpr(Interval.open(2, 4))
    assert i12o - i12ro == SetExpr(Interval.open(-1, 1))
    assert i12o * i12ro == SetExpr(Interval.open(1, 4))
    assert i12o / i12ro == SetExpr(Interval.open(S.Half, 2))
    assert i12o + i12cc == SetExpr(Interval.open(2, 4))
    assert i12o - i12cc == SetExpr(Interval.open(-1, 1))
    assert i12o * i12cc == SetExpr(Interval.open(1, 4))
    assert i12o / i12cc == SetExpr(Interval.open(S.Half, 2))
    assert i12o ** 2 == SetExpr(Interval.open(1, 4))
    assert i12o ** 3 == SetExpr(Interval.open(1, 8))

    assert n23cc + n23cc == SetExpr(Interval(-4, 6))
    assert n23cc - n23cc == SetExpr(Interval(-5, 5))
    assert n23cc * n23cc == SetExpr(Interval(-6, 9))
    assert n23cc / n23cc == SetExpr(Interval.open(-oo, oo))
    assert n23cc + n23ro == SetExpr(Interval.Ropen(-4, 6))
    assert n23cc - n23ro == SetExpr(Interval.Lopen(-5, 5))
    assert n23cc * n23ro == SetExpr(Interval.Ropen(-6, 9))
    assert n23cc / n23ro == SetExpr(Interval.Lopen(-oo, oo))
    assert n23cc + n23lo == SetExpr(Interval.Lopen(-4, 6))
    assert n23cc - n23lo == SetExpr(Interval.Ropen(-5, 5))
    assert n23cc * n23lo == SetExpr(Interval(-6, 9))
    assert n23cc / n23lo == SetExpr(Interval.open(-oo, oo))
    assert n23cc + n23o == SetExpr(Interval.open(-4, 6))
    assert n23cc - n23o == SetExpr(Interval.open(-5, 5))
    assert n23cc * n23o == SetExpr(Interval.open(-6, 9))
    assert n23cc / n23o == SetExpr(Interval.open(-oo, oo))
    assert n23cc ** 2 == SetExpr(Interval(0, 9))
    assert n23cc ** 3 == SetExpr(Interval(-8, 27))

    n32cc = SetExpr(Interval(-3, 2))
    n32lo = SetExpr(Interval.Lopen(-3, 2))
    n32ro = SetExpr(Interval.Ropen(-3, 2))
    assert n32cc * n32lo == SetExpr(Interval.Ropen(-6, 9))
    assert n32cc * n32cc == SetExpr(Interval(-6, 9))
    assert n32lo * n32cc == SetExpr(Interval.Ropen(-6, 9))
    assert n32cc * n32ro == SetExpr(Interval(-6, 9))
    assert n32lo * n32ro == SetExpr(Interval.Ropen(-6, 9))
    assert n32cc / n32lo == SetExpr(Interval.Ropen(-oo, oo))
    assert i12cc / n32lo == SetExpr(Interval.Ropen(-oo, oo))

    assert n3n2cc ** 2 == SetExpr(Interval(4, 9))
    assert n3n2cc ** 3 == SetExpr(Interval(-27, -8))

    assert n23cc + i12cc == SetExpr(Interval(-1, 5))
    assert n23cc - i12cc == SetExpr(Interval(-4, 2))
    assert n23cc * i12cc == SetExpr(Interval(-4, 6))
    assert n23cc / i12cc == SetExpr(Interval(-2, 3))
