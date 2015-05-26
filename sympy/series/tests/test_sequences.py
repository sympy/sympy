from sympy import (S, Tuple, Interval, EmptySequence, oo, SeqPer)
from sympy.series.sequences import SeqExpr
from sympy.utilities.pytest import raises


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
