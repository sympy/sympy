from sympy import (S, Tuple, Interval)
from sympy.series.sequences import (SeqExpr, EmptySequence)
from sympy.utilities.pytest import raises


def test_EmptySequence():
    assert isinstance(S.EmptySequence, EmptySequence)

    assert S.EmptySequence.interval is S.EmptySet
    assert S.EmptySequence.length is S.Zero


def test_SeqExpr():
    s = SeqExpr((1, 2, 3), (0, 10))

    assert isinstance(s, SeqExpr)
    assert s.gen == Tuple(1, 2, 3)
    assert s.interval == Interval(0, 10)
    assert s.start == 0
    assert s.end == 10
    assert s.length == 11

    assert SeqExpr((1, 2, 3), (0, S.Infinity)).length is S.Infinity

    assert SeqExpr((1, 2, 3), (0, -1)) is S.EmptySequence

    raises(ValueError, lambda: SeqExpr((1, 2, 3), (1, None)))
    raises(ValueError, lambda: SeqExpr((1, 2, 3), (1, 2, 3)))
