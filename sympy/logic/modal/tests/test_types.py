from __future__ import annotations
import pytest
from sympy.logic.modal.types import Universe, BoolType, FunctionType, PredicateVariable, ModalPredicate
from sympy.core.symbol import Symbol

def test_universe():
    u0 = Universe(0)
    u1 = Universe(1)
    assert u0.level == 0
    assert u1.level == 1
    assert str(u0) == "Type_0"

    with pytest.raises(ValueError):
        Universe(-1)

def test_function_type():
    u0 = Universe(0)
    b = BoolType()
    f = FunctionType(u0, b)
    assert f.domain == u0
    assert f.codomain == b

def test_predicate_variable():
    P = PredicateVariable('P', type=FunctionType(Universe(0), BoolType()))
    assert P.name == 'P'
    assert isinstance(P.type, FunctionType)

    x = Symbol('x')
    px = P(x)
    assert isinstance(px, ModalPredicate)
    assert px.name == 'P'
    assert px.args[0] == x
