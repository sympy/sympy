from sympy.sets.setexpr import SetExpr
from sympy.core.sets import Interval, FiniteSet
from sympy import Expr, Set, exp, log, sin, cos, Symbol, Min, Max

I = Interval(0, 2)


def test_setexpr():
    se = SetExpr(Interval(0, 1))
    assert isinstance(se.set, Set)
    assert isinstance(se, Expr)


def test_scalar_funcs():
    assert SetExpr(Interval(0, 1)).set == Interval(0, 1)
    a, b = Symbol('a', real=True), Symbol('b', real=True)
    a, b = 1, 2
    for f in [exp, log, sin, cos]:
        input = f(SetExpr(Interval(a, b)))
        output = input.set
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

    assert (a + a).set == FiniteSet(0, 2)

