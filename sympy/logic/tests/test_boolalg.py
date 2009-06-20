from sympy.logic.boolalg import And, Or, Not, Implies, Equivalent, to_cnf, \
    eliminate_implications, distribute_and_over_or
from sympy import symbols
from sympy.utilities.pytest import raises, XFAIL

def test_overloading():
    """Test that |, & are overloaded as expected"""
    A, B, C = symbols('ABC')
    assert A & B == And(A, B)
    assert A | B == Or(A, B)
    assert (A & B) | C == Or(And(A, B), C)
    assert A >> B == Implies(A, B)
    assert A << B == Implies(B, A)
    assert ~A == Not(A)

def test_And():
    A, B, C = symbols('ABC')
    assert And() == True
    assert And(A)== A
    assert And(True) == True
    assert And(False) == False
    assert And(True,  True ) == True
    assert And(True,  False) == False
    assert And(False, False) == False

def test_Or():
    A, B, C = symbols('ABC')
    assert Or() == False
    assert Or(A) == A
    assert Or(True) == True
    assert Or(False) == False
    assert Or(True,  True ) == True
    assert Or(True,  False) == True
    assert Or(False, False) == False

def test_Not():
    assert Not(True, True ) == [False, False]
    assert Not(True, False) == [False, True ]
    assert Not(False,False) == [True,  True ]

def test_Implies():
    A, B, C = symbols('ABC')
    raises(ValueError, "Implies()")
    raises(ValueError, "Implies(A)")
    raises(ValueError, "Implies(A, B, C)")
    assert A >> B == B << A

@XFAIL
def test_Equivalent():
    A, B, C = symbols('ABC')
    assert Equivalent(A, B) == Equivalent(B, A)

def test_bool_symbol():
    """Test that mixing symbols with boolean values
    works as expected"""
    A, B, C = symbols('ABC')
    assert And(A, True)  == A
    assert And(A, True, True) == A
    assert And(A, False) == False
    assert And(A, True, False) == False
    assert Or(A, True)   == True
    assert Or(A, False)  == A

def test_subs():
    A, B, C = symbols('ABC')
    assert (A & B).subs(A, True) == B
    assert (A & B).subs(A, False) == False
    assert (A & B).subs(B, True) == A
    assert (A & B).subs(B, False) == False
    assert (A & B).subs({A: True, B:True}) == True
    assert (A | B).subs(A, True) == True
    assert (A | B).subs(A, False) == B
    assert (A | B).subs(B, True) == True
    assert (A | B).subs(B, False) == A

@XFAIL
def test_subs_xfail():
    A, B, C = symbols('ABC')
    assert (A | B).subs({A: True, B:True}) == True

"""
we test for axioms of boolean algebra
see http://en.wikipedia.org/wiki/Boolean_algebra_(structure)
"""

def test_commutative():
    """Test for commutativity of And and Not"""
    A, B = symbols('AB')
    assert A & B == B & A
    assert A | B == B | A

def test_and_associativity():
    """Test for associativity of And"""
    A, B, C = symbols('ABC')
    assert (A & B) & C == A & (B & C)

def test_or_assicativity():
    A, B, C = symbols('ABC')
    assert ((A | B) | C) == (A | (B | C))

def test_double_negation():
    a = symbols('a')
    assert ~(~a) == a

def test_De_Morgan():
    A, B, C = symbols('ABC')
    assert ~(A & B) == (~A) | (~B)
    assert ~(A | B) == (~A) & (~B)
    assert ~(A | B | C) == ~A & ~B  & ~C

# test methods
def test_eliminate_implications():
    A, B, C = symbols('ABC')
    assert eliminate_implications( A >> B) == (~A) | B
    assert eliminate_implications(A >> (C >>Not(B))) \
        == Or(Or(Not(B), Not(C)), Not(A))

def test_distribute():
    A, B, C = symbols('ABC')
    assert distribute_and_over_or(Or(And(A, B), C)) == And(Or(A, C), Or(B, C))

def test_to_cnf():
    A, B, C = symbols('ABC')
    assert to_cnf(~(B | C)) == And(Not(B), Not(C))
    assert to_cnf((A & B) | C) == And(Or(A, C), Or(B, C))
    assert to_cnf(A >> B) == (~A) | B
    assert to_cnf(A >> (B & C)) == (~A | B) & (~A | C)

    assert to_cnf(Equivalent(A, B)) == And(Or(A, Not(B)), Or(B, Not(A)))
    assert to_cnf(Equivalent(A, B & C)) == (~A | B) & (~A | C) & (~B | ~C | A)
    assert to_cnf(Equivalent(A, B | C)) == \
    And(Or(Not(B), A), Or(Not(C), A), Or(B, C, Not(A)))
