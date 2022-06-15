from multiprocessing.context import assert_spawning
from sympy.core.symbol import (Symbol, symbols)
from sympy.physics.continuum_mechanics.truss import Truss

def test_truss():
    A = Symbol('A')
    B = Symbol('B')
    C = Symbol('C')
    AB, BC, AC = symbols('AB, BC, AC')
    P = Symbol('P')
    
    t = Truss()
    assert t.nodes == []
    assert t.node_labels == []
    assert t.node_positions == []
    assert t.members == []
    assert t.member_labels == []
    assert t.loads == {}
    assert t.supports == {}

    # testing the add_node method
    t.add_node(A, 0, 0)
    t.add_node(B, 2, 2)
    t.add_node(C, 3, 0)
    assert t.nodes == [(A, 0, 0), (B, 2, 2), (C, 3, 0)]
    assert t.node_labels == [A, B, C]
    assert t.node_positions == [(0, 0), (2, 2), (3, 0)]
    assert t.loads == {A: [[0, 90]], B: [[0, 90]], C: [[0, 90]]}
    assert t.supports == {A: 'none', B: 'none', C: 'none'}

    # testing the remove_node method
    t.remove_node(C)
    assert t.nodes == [(A, 0, 0), (B, 2, 2)]
    assert t.node_labels == [A, B]
    assert t.node_positions == [(0, 0), (2, 2)]
    assert t.loads == {A: [[0, 90]], B: [[0, 90]]}
    assert t.supports == {A: 'none', B: 'none'}

    t.add_node(C, 3, 0)

    # testing the add_member method
    t.add_member(AB, A, B)
    t.add_member(BC, B, C)
    t.add_member(AC, A, C)
    assert t.members == [(AB, A, B), (BC, B, C), (AC, A, C)]
    assert t.member_labels == [AB, BC, AC]

    # testing the remove_member method
    t.remove_member(BC)
    assert t.members == [(AB, A, B), (AC, A, C)]
    assert t.member_labels == [AB, AC]

    t.add_member(BC, B, C)

    # testing the apply_load method
    t.apply_load(A, P, 90)
    t.apply_load(A, P/4, 90)
    t.apply_load(A, 2*P,45)
    t.apply_load(B, P/2, 90)
    assert t.loads == {A: [[5*P/4, 90], [2*P, 45]], B: [[P/2, 90]], C: [[0, 90]]}
    assert t.loads[A] == [[5*P/4, 90], [2*P, 45]]

    # testing the add_support method
    t.add_support(A, "pinned")
    t.add_support(B, "roller")
    assert t.supports == {A: 'pinned', B: 'roller', C: 'none'}

    # testing the remove_support method
    t.remove_support(A)
    assert t.supports == {A: 'none', B: 'roller', C: 'none'}

