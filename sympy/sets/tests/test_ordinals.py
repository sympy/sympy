from sympy.sets import Ordinal, Ordinals

def test_string_ordinals():
    assert str(Ordinals.w) == 'w'
    assert str(Ordinal([[5, 3], [3, 2]])) == '{w**5}*3 + {w**3}*2'
    assert str(Ordinal([[5, 3], [0, 5]])) == '{w**5}*3 + 5'

def test_addition_with_integers():
    assert Ordinal([[5, 3]])+3 == Ordinal([[5, 3], [0, 3]])
    assert Ordinal([[5, 3], [0, 2]])+3 == Ordinal([[5, 3], [0, 5]])

def test_addition_with_ordinals():
    assert Ordinal([[5, 3], [3, 2]]) + Ordinal([[3, 3]]) == \
        Ordinal([[5, 3], [3, 5]])
    assert Ordinal([[5, 3], [3, 2]]) + Ordinal([[4, 2]]) == \
        Ordinal([[5, 3], [4, 2]])
    assert Ordinal([[Ordinals.w, 2], [3, 2]]) + Ordinal([[4, 2]]) == \
        Ordinal([[Ordinals.w, 2], [4, 2]])

def test_comparison():
    assert Ordinal([[5, 3]]) > Ordinal([[4, 3], [2, 0]])
    assert Ordinal([[5, 3], [3, 2]]) < Ordinal([[5, 4]])
    assert Ordinal([[5, 4]]) < Ordinal([[5, 5], [4, 1]])
    assert Ordinal([[5,3], [3,2]]) == Ordinal([[5,3], [3,2]])
    assert not Ordinal([[5,3], [3,2]]) == Ordinal([[5,3]])
    assert Ordinal([[Ordinals.w, 3]]) > Ordinal([[5, 3]])
