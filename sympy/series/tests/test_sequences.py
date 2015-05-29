from sympy import (S, Tuple, symbols, Interval, EmptySequence, oo, SeqPer\
                   , SeqFormula, SeqFunc, Lambda, sequence)
from sympy.series.sequences import SeqExpr, SeqExprOp
from sympy.utilities.pytest import raises

x, y, z = symbols('x y z')
n, m = symbols('n m')


def test_EmptySequence():
    assert isinstance(S.EmptySequence, EmptySequence)

    assert S.EmptySequence.interval is S.EmptySet
    assert S.EmptySequence.length is S.Zero

    assert list(S.EmptySequence) == []


def test_SeqExpr():
    s = SeqExpr((1, 2, 3), (0, 10))

    assert isinstance(s, SeqExpr)
    assert s.gen == Tuple(1, 2, 3)
    assert s.interval == Interval(0, 10)
    assert s.start == 0
    assert s.stop == 10
    assert s.length == 11

    assert SeqExpr((1, 2, 3), (0, 10, 2)).length == 6
    assert SeqExpr((1, 2, 3), (0, oo)).length is oo

    assert SeqExpr((1, 2, 3), (oo, -oo)) is S.EmptySequence

    raises(ValueError, lambda: SeqExpr((1, 2, 3), (0, 1, 2, 3)))
    raises(ValueError, lambda: SeqExpr((1, 2, 3), (-oo, oo)))
    raises(ValueError, lambda: SeqExpr((1, 2, 3), (0, oo, oo)))


def test_SeqPer():
    s = SeqPer((1, 2, 3), (0, 5))

    assert isinstance(s, SeqPer)
    assert s.periodical == Tuple(1, 2, 3)
    assert s.period == 3
    assert s.coeff(3) == 1

    assert list(s) == [1, 2, 3, 1, 2, 3]
    assert s[:] == [1, 2, 3, 1, 2, 3]
    assert SeqPer((1, 2, 3), (0, 5, 2))[:] == [1, 3, 2]
    assert SeqPer((1, 2, 3), (-oo, 0))[0:6] == [1, 2, 3, 1, 2, 3]


def test_SeqFormula():
    s = SeqFormula((n**2, n), (0, 5))

    assert isinstance(s, SeqFormula)
    assert s.formula == n**2
    assert s.coeff(3) == 9

    assert list(s) == [i**2 for i in range(6)]
    assert s[:] == [i**2 for i in range(6)]
    assert SeqFormula((n**2, n), (0, 5, 2))[:] == [0, 4, 16]
    assert SeqFormula((n**2, n), (-oo, 0))[0:6] == [i**2 for i in range(6)]

    assert SeqFormula(n**2, (0, oo)) == SeqFormula((n**2, n), (0, oo))

    assert SeqFormula(n**2, (0, m)).subs(m, x) == SeqFormula(n**2, (0, x))
    assert SeqFormula((m*n**2, n), (0, oo)).subs(m, x) == \
        SeqFormula((x*n**2, n), (0, oo))


def test_SeqFunc():
    s = SeqFunc(Lambda(n, n**2), (0, 5))

    assert isinstance(s, SeqFunc)
    assert s.function == Lambda(n, n**2)
    assert s.coeff(3) == 9

    assert list(s) == [i**2 for i in range(6)]
    assert s[:] == [i**2 for i in range(6)]
    assert SeqFunc(Lambda(n, n**2), (0, 5, 2))[:] == [0, 4, 16]
    assert SeqFunc(Lambda(n, n**2), (-oo, 0))[0:6] == [i**2 for i in range(6)]


def test_sequence():
    form = SeqFormula((n**2, n), (0, 5))
    per = SeqPer((1, 2, 3), (0, 5))
    func = SeqFunc(Lambda(n, n**2), (0, 5))
    inter = SeqFormula((n**2, n))

    assert sequence(formula=(n**2, n), interval=(0, 5)) == form
    assert sequence(periodical=(1, 2, 3), interval=(0, 5)) == per
    assert sequence(func=Lambda(n, n**2), interval=(0, 5)) == func
    assert sequence(formula=(n**2, n)) == inter


def test_SeqExprOp():
    form = SeqFormula((n**2, n), (0, 10))
    per = SeqPer((1, 2, 3), (5, 10))
    func = SeqFunc(Lambda(m, m**2), (0, 10))

    s = SeqExprOp(form, per, func)
    assert s.gen == ((n**2, n), (1, 2, 3), Lambda(m, m**2))
    assert s.interval == Interval(5, 10)
    assert s.start == 5
    assert s.stop == 10
    assert s.length == 6
    assert s.variables == (n, m)
