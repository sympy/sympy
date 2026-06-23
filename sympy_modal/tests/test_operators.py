import pytest
from sympy.core.symbol import Symbol
from sympy_modal.types import PredicateVariable, FunctionType, Universe, BoolType
from sympy_modal.operators import Box, Diamond, ProvabilityBox, ForAllPredicates, ExistsPredicates

def test_box_diamond():
    p = Symbol('p')
    b = Box(p)
    d = Diamond(p)
    assert b.args[0] == p
    assert d.args[0] == p
    assert b.is_well_typed()
    assert d.is_well_typed()

def test_named_boxes():
    p = Symbol('p')
    pb = ProvabilityBox(p)
    assert pb.modality == "provability"
    assert isinstance(pb, Box)

def test_quantifiers():
    P = PredicateVariable('P', type=FunctionType(Universe(0), BoolType()))
    x = Symbol('x')
    px = P(x)

    fa = ForAllPredicates(P, px)
    assert fa.variable == P
    assert fa.formula == px
    assert fa.is_well_typed()

    ex = ExistsPredicates(P, px)
    assert ex.variable == P
    assert ex.formula == px
    assert ex.is_well_typed()
