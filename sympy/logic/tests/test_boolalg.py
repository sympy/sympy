from sympy.logic.boolalg import to_cnf, eliminate_implications, distribute_and_over_or, \
    compile_rule, conjuncts, disjuncts, to_int_repr, fuzzy_not
from sympy import symbols, And, Or, Xor, Not, Nand, Nor, Implies, Equivalent
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
    assert And(A) == A
    assert And(True) == True
    assert And(False) == False
    assert And(True,  True ) == True
    assert And(True,  False) == False
    assert And(False, False) == False
    assert And(True,  A) == A
    assert And(False, A) == False
    assert And(True, True, True) == True
    assert And(True, True , A) == A
    assert And(True, False, A) == False

def test_Or():
    A, B, C = symbols('ABC')
    assert Or() == False
    assert Or(A) == A
    assert Or(True) == True
    assert Or(False) == False
    assert Or(True,  True ) == True
    assert Or(True,  False) == True
    assert Or(False, False) == False
    assert Or(True, A) == True
    assert Or(False, A) == A
    assert Or(True, False, False) == True
    assert Or(True, False, A) == True
    assert Or(False, False, A) == A

def test_Xor():
    A, B, C = symbols('ABC')
    assert Xor() == False
    assert Xor(A) == A
    assert Xor(True) == True
    assert Xor(False) == False
    assert Xor(True,  True ) == False
    assert Xor(True,  False) == True
    assert Xor(False, False) == False
    assert Xor(True, A) == ~A
    assert Xor(False, A) == A
    assert Xor(True, False, False) == True
    assert Xor(True, False, A) == ~A
    assert Xor(False, False, A) == A

def test_Not():
    assert Not(True) == False
    assert Not(False) == True
    assert Not(True, True ) == [False, False]
    assert Not(True, False) == [False, True ]
    assert Not(False,False) == [True,  True ]

def test_Nand():
    A, B, C = symbols('ABC')
    assert Nand() == False
    assert Nand(A) == ~A
    assert Nand(True) == False
    assert Nand(False) == True
    assert Nand(True,  True ) == False
    assert Nand(True,  False) == True
    assert Nand(False, False) == True
    assert Nand(True,  A) == ~A
    assert Nand(False, A) == True
    assert Nand(True, True, True) == False
    assert Nand(True, True , A) == ~A
    assert Nand(True, False, A) == True

def test_Nor():
    A, B, C = symbols('ABC')
    assert Nor() == False
    assert Nor(A) == ~A
    assert Nor(True) == False
    assert Nor(False) == True
    assert Nor(True,  True ) == False
    assert Nor(True,  False) == False
    assert Nor(False, False) == True
    assert Nor(True,  A) == False
    assert Nor(False, A) == ~A
    assert Nor(True, True, True) == False
    assert Nor(True, True , A) == False
    assert Nor(True, False, A) == False

def test_Implies():
    A, B, C = symbols('ABC')
    Implies(True, True) == True
    Implies(False, False) == False
    assert A >> B == B << A

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

def test_conjuncts():
    A, B, C = symbols('ABC')
    assert set(conjuncts(A & B & C)) == set([A, B, C])
    assert set(conjuncts((A | B) & C)) == set([A | B, C])
    assert conjuncts(A) == [A]
    assert conjuncts(True) == [True]
    assert conjuncts(False) == [False]

def test_disjuncts():
    A, B, C = symbols('ABC')
    assert disjuncts(A | B | C) == [A, B, C]
    assert disjuncts((A | B) & C) == [(A | B) & C]
    assert disjuncts(A) == [A]
    assert disjuncts(True) == [True]
    assert disjuncts(False) == [False]

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

def test_compile_rule():
    from sympy import sympify
    assert compile_rule("A & B") == sympify("A & B")

def test_to_int_repr():
    x, y, z = symbols('x y z')
    def sorted_recursive(arg):
        try:
            return sorted(sorted_recursive(x) for x in arg)
        except TypeError:   #arg is not a sequence
            return arg

    assert sorted_recursive(to_int_repr([x | y, z | x], [x, y, z])) == \
                                            sorted_recursive([[1, 2], [1, 3]])
    assert sorted_recursive(to_int_repr([x | y, z | ~x], [x, y, z])) == \
                                            sorted_recursive([[1, 2], [3, -1]])

def test_fuzzy_not():
    assert fuzzy_not(False) == True
    assert fuzzy_not(True) == False
    assert fuzzy_not(None) == None
