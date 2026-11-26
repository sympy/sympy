from sympy.core.singleton import S
from sympy.core.symbol import Symbol
from sympy.sets.subset import Subset
from sympy.sets.sets import FiniteSet
from sympy.testing.pytest import raises


def test_subset_basic():
    raises(TypeError, lambda: Subset(S.Integers, 1))
    raises(TypeError, lambda: Subset(2, S.Naturals))
    assert Subset(S.Naturals, S.Integers) is S.true
    assert Subset(S.Integers, S.Naturals) is S.false

    i = FiniteSet(Symbol('i'))
    assert Subset(i, S.Naturals) == Subset(i, S.Naturals, evaluate=False)


def test_type_error():
    # Pass in a parameter not of type "set"
    raises(TypeError, lambda: Subset(2, S.Integers))
    raises(TypeError, lambda: Subset(S.Integers, 2))
