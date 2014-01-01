from sympy.sets.setexpr import SetExpr, simplify
from sympy.core.sets import Interval, FiniteSet
from sympy import Expr, Set, exp, log, sin, cos, Symbol, Min, Max

I = Interval(0, 2)


def test_setexpr():
    se = SetExpr(Interval(0, 1))
    assert isinstance(se.set, Set)
    assert isinstance(se, Expr)


def test_scalar_funcs():
    assert simplify(SetExpr(Interval(0, 1))) == SetExpr(Interval(0, 1))
    a, b = Symbol('a', real=True), Symbol('b', real=True)
    a, b = 1, 2
    for f in [exp, log, sin, cos]:
        input = f(SetExpr(Interval(a, b)))
        output = simplify(input)
        expected = SetExpr(Interval(Min(f(a), f(b)), Max(f(a), f(b))))
        assert output == expected


def test_Add_Mul():
    assert simplify(SetExpr(Interval(0, 1)) + 1) == SetExpr(Interval(1, 2))
    assert simplify(SetExpr(Interval(0, 1)) * 2) == SetExpr(Interval(0, 2))


def test_Pow():
    assert simplify(SetExpr(Interval(0, 2))**2) == SetExpr(Interval(0, 4))


def test_compound():
    assert simplify(exp(SetExpr(Interval(0, 1)) * 2 + 1)) == \
        SetExpr(Interval(exp(1), exp(3)))


def test_Interval_Interval():
    assert simplify(SetExpr(Interval(1, 2)) + SetExpr(Interval(10, 20))) == \
        SetExpr(Interval(11, 22))
    assert simplify(SetExpr(Interval(1, 2)) * SetExpr(Interval(10, 20))) == \
        SetExpr(Interval(10, 40))


def test_FiniteSet_FiniteSet():
    assert simplify(SetExpr(FiniteSet(1, 2, 3)) + SetExpr(FiniteSet(1, 2))) ==\
        SetExpr(FiniteSet(2, 3, 4, 5))
    assert simplify(SetExpr(FiniteSet(1, 2, 3)) * SetExpr(FiniteSet(1, 2))) ==\
        SetExpr(FiniteSet(1, 2, 3, 4, 6))


def test_Interval_FiniteSet():
    assert simplify(SetExpr(FiniteSet(1, 2)) + SetExpr(Interval(0, 10))) == \
        SetExpr(Interval(1, 12))


def test_Many_Sets():
    assert simplify(SetExpr(Interval(0, 1)) +
                    SetExpr(Interval(2, 3)) +
                    SetExpr(FiniteSet(10, 11, 12))) == SetExpr(Interval(12, 16))
