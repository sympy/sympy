from sympy.sets import Ordinal, Ordinals, OrdinalPow

def test_string_ordinals():
    assert str(Ordinals.w) == 'w'
    assert str(Ordinal([OrdinalPow(5, 3), OrdinalPow(3, 2)])) == '{w**5}*3 + {w**3}*2'
    assert str(Ordinal([OrdinalPow(5, 3), OrdinalPow(0, 5)])) == '{w**5}*3 + 5'

def test_addition_with_integers():
    assert Ordinal([OrdinalPow(5, 3)])+3 == Ordinal([OrdinalPow(5, 3), OrdinalPow(0, 3)])
    assert Ordinal([OrdinalPow(5, 3), OrdinalPow(0, 2)])+3 == Ordinal([OrdinalPow(5, 3), OrdinalPow(0, 5)])

def test_addition_with_ordinals():
    assert Ordinal([OrdinalPow(5, 3), OrdinalPow(3, 2)]) + Ordinal([OrdinalPow(3, 3)]) == \
        Ordinal([OrdinalPow(5, 3), OrdinalPow(3, 5)])
    assert Ordinal([OrdinalPow(5, 3), OrdinalPow(3, 2)]) + Ordinal([OrdinalPow(4, 2)]) == \
        Ordinal([OrdinalPow(5, 3), OrdinalPow(4, 2)])
    assert Ordinal([OrdinalPow(Ordinals.w, 2), OrdinalPow(3, 2)]) + Ordinal([OrdinalPow(4, 2)]) == \
        Ordinal([OrdinalPow(Ordinals.w, 2), OrdinalPow(4, 2)])

def test_comparison():
    assert Ordinal([OrdinalPow(5, 3)]) > Ordinal([OrdinalPow(4, 3), OrdinalPow(2, 0)])
    assert Ordinal([OrdinalPow(5, 3), OrdinalPow(3, 2)]) < Ordinal([OrdinalPow(5, 4)])
    assert Ordinal([OrdinalPow(5, 4)]) < Ordinal([OrdinalPow(5, 5), OrdinalPow(4, 1)])
    assert Ordinal([OrdinalPow(5, 3), OrdinalPow(3, 2)]) == Ordinal([OrdinalPow(5, 3), OrdinalPow(3, 2)])
    assert not Ordinal([OrdinalPow(5, 3), OrdinalPow(3, 2)]) == Ordinal([OrdinalPow(5, 3)])
    assert Ordinal([OrdinalPow(Ordinals.w, 3)]) > Ordinal([OrdinalPow(5, 3)])
