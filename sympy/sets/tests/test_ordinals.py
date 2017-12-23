from sympy.sets import Ordinal, OmegaPower, Ord0, OrdOmega

def test_string_ordinals():
    assert str(OrdOmega) == 'w'
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
    assert Ordinal(OmegaPower(OrdOmega, 2), OmegaPower(3, 2)) + Ordinal(OmegaPower(4, 2)) == \
        Ordinal(OmegaPower(OrdOmega, 2), OmegaPower(4, 2))

def test_comparison():
    assert Ordinal(OmegaPower(5, 3)) > Ordinal(OmegaPower(4, 3), OmegaPower(2, 1))
    assert Ordinal(OmegaPower(5, 3), OmegaPower(3, 2)) < Ordinal(OmegaPower(5, 4))
    assert Ordinal(OmegaPower(5, 4)) < Ordinal(OmegaPower(5, 5), OmegaPower(4, 1))
    assert Ordinal(OmegaPower(5, 3), OmegaPower(3, 2)) == Ordinal(OmegaPower(5, 3), OmegaPower(3, 2))
    assert not Ordinal(OmegaPower(5, 3), OmegaPower(3, 2)) == Ordinal(OmegaPower(5, 3))
    assert Ordinal(OmegaPower(OrdOmega, 3)) > Ordinal(OmegaPower(5, 3))

def test_multiplication_with_integers():
    w = OrdOmega
    assert 3*w == w
    assert w*9 == Ordinal(OmegaPower(1,9))

def test_multiplication():
    w = OrdOmega
    w*(w + 1) == w*w + w
    (w + 1)*(w + 1) ==  w*w + w + 1
    w*1 == w
    1*w == w
    w*Ord0 == Ord0
    Ord0*w == Ord0
