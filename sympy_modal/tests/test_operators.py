import pytest
from sympy.core.symbol import Symbol
from sympy_modal.types import PredicateVariable, FunctionType, Universe, BoolType
from sympy_modal.operators import (
    Box, Diamond, ProvabilityBox, ForAllPredicates, ExistsPredicates,
    AgentBox, CommonKnowledge, Next, Until
)

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

def test_expressiveness_operators():
    p = Symbol('p')
    q = Symbol('q')

    ab = AgentBox('Alice', p)
    assert ab.agent == 'Alice'
    assert ab.modality == 'epistemic_Alice'

    ck = CommonKnowledge('GroupA', p)
    assert ck.group == 'GroupA'
    assert ck.modality == 'common_knowledge_GroupA'

    nx = Next(p)
    assert nx.modality == 'temporal_next'
    assert nx.is_well_typed()

    ut = Until(p, q)
    assert ut.left == p
    assert ut.right == q
    assert ut.is_well_typed()
def test_agent_box_subs():
    from sympy.core.symbol import Symbol
    from sympy_modal.operators import AgentBox
    p = Symbol('p')
    q = Symbol('q')
    ab = AgentBox('Alice', p)
    # This should not crash and should correctly reconstruct the AST
    ab_sub = ab.subs(p, q)
    assert isinstance(ab_sub, AgentBox)
    assert ab_sub.agent == 'Alice'
    assert ab_sub.formula == q
