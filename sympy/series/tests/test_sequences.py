from sympy import S, EmptySequence


def test_EmptySequence():
    assert isinstance(S.EmptySequence, EmptySequence)

    assert S.EmptySequence.interval is S.EmptySet
    assert S.EmptySequence.length is S.Zero
