"""For more tests on satisfiability, see test_dimacs"""

from sympy import symbols
from sympy.logic.boolalg import Equivalent, Implies
from sympy.logic.inference import pl_true, satisfiable, PropKB
from sympy.logic.algorithms.dpll import dpll, dpll_satisfiable, \
    find_pure_symbol, find_unit_clause, unit_propagate, \
    find_pure_symbol_int_repr, find_unit_clause_int_repr, \
    unit_propagate_int_repr
from sympy.utilities.pytest import raises, XFAIL

def test_find_pure_symbol():
    A, B, C = symbols('ABC')
    assert find_pure_symbol([A], [A]) == (A, True)
    assert find_pure_symbol([A, B], [~A | B, ~B | A]) == (None, None)
    assert find_pure_symbol([A, B, C], [ A | ~B, ~B | ~C, C | A]) == (A, True)
    assert find_pure_symbol([A, B, C], [~A |  B,  B | ~C, C | A]) == (B, True)
    assert find_pure_symbol([A, B, C], [~A | ~B, ~B | ~C, C | A]) == (B, False)
    assert find_pure_symbol([A, B, C], [~A | B, ~B | ~C, C | A]) == (None, None)

def test_find_pure_symbol_int_repr():
    assert find_pure_symbol_int_repr([1], [[1]]) == (1, True)
    assert find_pure_symbol_int_repr([1, 2], [[-1, 2], [-2, 1]]) == (None, None)
    assert find_pure_symbol_int_repr([1, 2, 3], [[1, -2], [-2, -3], [3, 1]]) == \
                                         (1, True)
    assert find_pure_symbol_int_repr([1, 2, 3], [[-1, 2], [2, -3], [3, 1]]) == \
        (2, True)
    assert find_pure_symbol_int_repr([1, 2, 3], [[-1, -2], [-2, -3], [3, 1]]) == \
        (2, False)
    assert find_pure_symbol_int_repr([1, 2, 3], [[-1, 2], [-2, -3], [3, 1]]) == \
        (None, None)

def test_unit_clause():
    A, B, C = symbols('ABC')
    assert find_unit_clause([A], {}) == (A, True)
    assert find_unit_clause([A, ~A], {}) == (A, True) ### Wrong ??
    assert find_unit_clause([A | B], {A: True}) == (B, True)
    assert find_unit_clause([A | B], {B: True}) == (A, True)
    assert find_unit_clause([A | B | C, B | ~C, A | ~B], {A:True}) == (B, False)
    assert find_unit_clause([A | B | C, B | ~C, A | B], {A:True})  == (B, True)
    assert find_unit_clause([A | B | C, B | ~C, A ], {}) == (A, True)

def test_unit_clause_int_repr():
    assert find_unit_clause_int_repr([[1]], {}) == (1, True)
    assert find_unit_clause_int_repr([[1], [-1]], {}) == (1, True)
    assert find_unit_clause_int_repr([[1,2]], {1: True}) == (2, True)
    assert find_unit_clause_int_repr([[1,2]], {2: True}) == (1, True)
    assert find_unit_clause_int_repr([[1,2,3], [2, -3], [1, -2]], {1: True}) == \
        (2, False)
    assert find_unit_clause_int_repr([[1, 2, 3], [3, -3], [1, 2]], {1: True}) == \
        (2, True)
#    assert find_unit_clause([A | B | C, B | ~C, A ], {}) == (A, True)

def test_unit_propagate():
    A, B, C = symbols('ABC')
    assert unit_propagate([A | B], A) == []
    assert unit_propagate([A | B, ~A | C, ~C | B, A], A) == [C, ~C | B, A]

def test_unit_propagate_int_repr():
    assert unit_propagate_int_repr([[1, 2]], 1) == []
    assert unit_propagate_int_repr([[1, 2], [-1, 3], [-3, 2], [1]], 1) == \
        [[3], [-3, 2], [1]]

def test_dpll():
    """This is also tested in test_dimacs"""
    A, B, C = symbols('ABC')
    assert dpll([A | B], [A, B], {A: True, B: True}) == {A: True, B: True}

def test_dpll_satisfiable():
    A, B, C = symbols('ABC')
    assert dpll_satisfiable( A & ~A ) == False
    assert dpll_satisfiable( A & ~B ) == {A: True, B: False}
    assert dpll_satisfiable( A | B ) in ({A: True}, {B: True}, {A: True, B: True})
    assert dpll_satisfiable( (~A | B) & (~B | A) ) == {A: True, B: True}
    assert dpll_satisfiable( (A | B) & (~B | C) ) == {A: True, B: False}
    assert dpll_satisfiable( A & B & C  ) == {A: True, B: True, C: True}
    assert dpll_satisfiable( (A | B) & (A >> B) ) == {B: True}
    assert dpll_satisfiable( Equivalent(A, B) & A ) == {A: True, B: True}
    assert dpll_satisfiable( Equivalent(A, B) & ~A ) == {A: False, B: False}

def test_satisfiable():
    A, B, C = symbols('ABC')
    assert satisfiable(A & (A >> B) & ~B) == False

def test_pl_true():
    A, B, C = symbols('ABC')
    assert pl_true(True) == True
    assert pl_true( A & B, {A : True, B : True}) == True
    assert pl_true( A | B, {A : True}) == True
    assert pl_true( A | B, {B : True}) == True
    assert pl_true( A | B, {A: None, B: True}) == True
    assert pl_true( A >> B, {A: False}) == True
    assert pl_true( A | B | ~C, {A: False, B: True, C: True}) == True
    assert pl_true(Equivalent(A, B), {A:False, B:False}) == True

    # test for false
    assert pl_true(False) == False
    assert pl_true( A & B, {A: False, B: False}) == False
    assert pl_true( A & B, {A: False}) == False
    assert pl_true( A & B, {B: False}) == False
    assert pl_true( A | B, {A: False, B: False}) == False

    #test for None
    assert pl_true(B, {B: None}) is None
    assert pl_true( A & B, {A: True, B: None}) is None
    assert pl_true( A >> B, {A: True, B: None}) is None
    assert pl_true(Equivalent(A, B), {A:None}) is None
    assert pl_true(Equivalent(A, B), {A:True, B:None}) is None

def test_pl_true_wrong_input():
    from sympy import pi
    raises(ValueError, "pl_true('John Cleese')")
    raises(ValueError, "pl_true(42+pi+pi**2)")
    #raises(ValueError, "pl_true(42)")  #returns None, but should it?

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
    assert kb.ask(~C) == False
    kb.retract(A)
    assert kb.ask(~C) == True

    kb2 = PropKB(Equivalent(A, B))
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
