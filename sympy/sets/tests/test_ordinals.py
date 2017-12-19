from sympy.sets import Ordinal, Ordinals, OmegaPower

def test_string_ordinals():
    assert str(Ordinals.w) == 'w'
    assert str(Ordinal(OmegaPower(5, 3), OmegaPower(3, 2))) == '{w**5}*3 + {w**3}*2'
    assert str(Ordinal(OmegaPower(5, 3), OmegaPower(0, 5))) == '{w**5}*3 + 5'
    assert str(Ordinal(OmegaPower(1, 3), OmegaPower(0, 5))) == 'w*3 + 5'

def test_addition_with_integers():
    assert 3 + Ordinal(OmegaPower(5, 3)) == Ordinal(OmegaPower(5, 3))
    assert Ordinal(OmegaPower(5, 3))+3 == Ordinal(OmegaPower(5, 3), OmegaPower(0, 3))
    assert Ordinal(OmegaPower(5, 3), OmegaPower(0, 2))+3 == Ordinal(OmegaPower(5, 3), OmegaPower(0, 5))

def test_addition_with_ordinals():
    assert Ordinal(OmegaPower(5, 3), OmegaPower(3, 2)) + Ordinal(OmegaPower(3, 3)) == \
        Ordinal(OmegaPower(5, 3), OmegaPower(3, 5))
    assert Ordinal(OmegaPower(5, 3), OmegaPower(3, 2)) + Ordinal(OmegaPower(4, 2)) == \
        Ordinal(OmegaPower(5, 3), OmegaPower(4, 2))
    assert Ordinal(OmegaPower(Ordinals.w, 2), OmegaPower(3, 2)) + Ordinal(OmegaPower(4, 2)) == \
        Ordinal(OmegaPower(Ordinals.w, 2), OmegaPower(4, 2))

def test_comparison():
    assert Ordinal(OmegaPower(5, 3)) > Ordinal(OmegaPower(4, 3), OmegaPower(2, 0))
    assert Ordinal(OmegaPower(5, 3), OmegaPower(3, 2)) < Ordinal(OmegaPower(5, 4))
    assert Ordinal(OmegaPower(5, 4)) < Ordinal(OmegaPower(5, 5), OmegaPower(4, 1))
    assert Ordinal(OmegaPower(5, 3), OmegaPower(3, 2)) == Ordinal(OmegaPower(5, 3), OmegaPower(3, 2))
    assert not Ordinal(OmegaPower(5, 3), OmegaPower(3, 2)) == Ordinal(OmegaPower(5, 3))
    assert Ordinal(OmegaPower(Ordinals.w, 3)) > Ordinal(OmegaPower(5, 3))
