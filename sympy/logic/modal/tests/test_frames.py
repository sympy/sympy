from __future__ import annotations
from sympy.core.symbol import Symbol
from sympy.logic.boolalg import Implies
from sympy.logic.modal.frames import KripkeFrame
from sympy.logic.modal.operators import Box

def test_kripke_frame_s4():
    s4 = KripkeFrame.S4()
    p = Symbol('p')
    # T axiom: Box(P) -> P
    assert s4.validates(Implies(Box(p), p))
    # 4 axiom: Box(P) -> Box(Box(P))
    assert s4.validates(Implies(Box(p), Box(Box(p))))
    # Löb axiom is not generally valid in S4
    lob = Implies(Box(Implies(Box(p), p)), Box(p))
    assert not s4.validates(lob)

def test_kripke_frame_gl():
    gl = KripkeFrame.GL()
    p = Symbol('p')

    # T axiom is NOT valid in GL (irreflexive)
    assert not gl.validates(Implies(Box(p), p))

    # Löb axiom is valid in GL
    lob = Implies(Box(Implies(Box(p), p)), Box(p))
    assert gl.validates(lob)

def test_valid_inference():
    s4 = KripkeFrame.S4()
    p = Symbol('p')

    assert s4.is_valid_inference([Box(p)], p)
