from sympy import (
    S, Symbol, Map, Set
)
from sympy.algebras.abstract.handlers import find_structure

x = Symbol('x')
f1, g1 = Map('f', Set('X'), S.RealsField), Map('g1', Set('X'), S.RealsField)

def test_Expr():

    assert find_structure(S.One) == S.ComplexesField
    assert find_structure(x) == S.ComplexesField

def test_Map():

    assert find_structure(f1) == find_structure(g1)
