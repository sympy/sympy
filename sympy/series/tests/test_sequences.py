from sympy import (S, Tuple, symbols, Interval, EmptySequence, oo, SeqPer\
                   , SeqFormula)
from sympy.series.sequences import SeqExpr
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

    assert SeqFormula((m*n**2, n), (0, oo)).subs(m, x) == \
        SeqFormula((x*n**2, n), (0, oo))
    assert SeqFormula((n**2, n), (0, m)).subs(m, x) == \
        SeqFormula((n**2, n), (0, x))
