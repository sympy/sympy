"""For more tests on satisfiability, see test_dimacs"""

from sympy import symbols
from sympy.logic.boolalg import Equivalent
from sympy.logic.inference import pl_true, satisfiable, find_pure_symbol, \
        find_unit_clause
from sympy.logic.kb import PropKB
from sympy.logic.algorithms.dpll import dpll, dpll_satisfiable

def test_find_pure_symbol():
    A, B, C = symbols('ABC')
    assert find_pure_symbol([A], [A]) == (A, True)
    assert find_pure_symbol([A, B, C], [ A | ~B, ~B | ~C, C | A]) == (A, True)
    assert find_pure_symbol([A, B, C], [~A |  B,  B | ~C, C | A]) == (B, True)
    assert find_pure_symbol([A, B, C], [~A | ~B, ~B | ~C, C | A]) == (B, False)
    assert find_pure_symbol([A, B, C], [~A | B, ~B | ~C, C | A]) == (None, None)

def test_unit_clause():
    A, B, C = symbols('ABC')
    assert find_unit_clause([A], {}) == (A, True)
    assert find_unit_clause([A | B], {A: True}) == (B, True)
    assert find_unit_clause([A | B], {B: True}) == (A, True)
    assert find_unit_clause([A | B | C, B | ~C, A | ~B], {A:True}) == (B, False)
    assert find_unit_clause([A | B | C, B | ~C, A | B], {A:True})  == (B, True)
    assert find_unit_clause([A | B | C, B | ~C, A ], {}) == (A, True)

def test_dpll():
    """This is also tested in test_dimacs"""
    A, B, C = symbols('ABC')
    assert dpll([A | B], [A, B], {A: True, B: True}) == {A: True, B: True}

def test_dpll_satisfiable():
    A, B, C = symbols('ABC')
    assert dpll_satisfiable( A & ~A )     == False
    assert dpll_satisfiable( A & ~B )     == {A: True, B: False}
    assert dpll_satisfiable( A & B & C  ) == {A: True, B: True, C: True}
    assert dpll_satisfiable( A | B )      in ({A: True}, {B: True})
    assert dpll_satisfiable( (A | B) & (A >> B)) == {B: True}
    assert dpll_satisfiable(Equivalent(A, B) & A) == {A: True, B: True}
    assert dpll_satisfiable(Equivalent(A, B) & ~A) == {A: False, B: False}

def test_satisfiable_1():
    """We prove expr entails alpha proving expr & ~B is unsatisfiable"""
    A, B, C = symbols('ABC')
    assert satisfiable(A & (A >> B) & ~B) == False

def test_pl_true():
    A, B, C = symbols('ABC')
    assert pl_true( A & B, {A : True, B : True}) == True
    assert pl_true( A | B, {A : True}) == True
    assert pl_true( A | B, {B : True}) == True
    assert pl_true( A | B | ~C, {A: False, B: True, C: True}) == True

    # test for false
    assert pl_true ( A & B, {A: False, B: False}) == False
    assert pl_true ( A & B, {A: False}) == False
    assert pl_true ( A & B, {B: False}) == False
    assert pl_true ( A | B, {A: False, B: False}) == False

def test_PropKB():
    A, B, C = symbols('ABC')
    kb = PropKB()
    kb.tell(A >> B)
    kb.tell(B >> C)
    assert kb.ask(A) == True
    assert kb.ask(B) == True
    assert kb.ask(C) == True
    assert kb.ask(~A) == True
    assert kb.ask(~B) == True
    assert kb.ask(~C) == True
    kb.tell(A)
    assert kb.ask(A) == True
    assert kb.ask(B) == True
    assert kb.ask(C) == True

    kb2 = PropKB()
    kb2.tell(Equivalent(A, B))
    assert kb2.ask(A) == True
    assert kb2.ask(B) == True
    kb2.tell(A)
    assert kb2.ask(A) == True

    kb3 = PropKB()
    kb3.tell(A)

def test_propKB_tolerant():
    """"tolerant to bad input"""
    kb = PropKB()
    A, B, C = symbols('ABC')
    assert kb.ask(B) == False
