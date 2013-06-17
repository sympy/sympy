from sympy import symbols, sympify, Dummy, simplify
from sympy.logic.boolalg import (
    And, Boolean, Equivalent, ITE, Implies, Nand, Nor, Not, Or, POSform,
    SOPform, Xor, conjuncts, disjuncts, distribute_or_over_and,
    distribute_and_over_or, eliminate_implications, is_cnf, is_dnf,
    simplify_logic, to_cnf, to_dnf, to_int_repr, bool_equal
)
from sympy.utilities.pytest import raises


A, B, C = symbols('A,B,C')


def test_overloading():
    """Test that |, & are overloaded as expected"""

    assert A & B == And(A, B)
    assert A | B == Or(A, B)
    assert (A & B) | C == Or(And(A, B), C)
    assert A >> B == Implies(A, B)
    assert A << B == Implies(B, A)
    assert ~A == Not(A)
    assert A ^ B == Xor(A, B)


def test_And():

    assert And() is True
    assert And(A) == A
    assert And(True) is True
    assert And(False) is False
    assert And(True, True ) is True
    assert And(True, False) is False
    assert And(False, False) is False
    assert And(True, A) == A
    assert And(False, A) is False
    assert And(True, True, True) is True
    assert And(True, True, A) == A
    assert And(True, False, A) is False
    assert And(2, A) == A
    assert And(2, 3) is True


def test_Or():

    assert Or() is False
    assert Or(A) == A
    assert Or(True) is True
    assert Or(False) is False
    assert Or(True, True ) is True
    assert Or(True, False) is True
    assert Or(False, False) is False
    assert Or(True, A) is True
    assert Or(False, A) == A
    assert Or(True, False, False) is True
    assert Or(True, False, A) is True
    assert Or(False, False, A) == A
    assert Or(2, A) is True


def test_Xor():

    assert Xor() is False
    assert Xor(A) == A
    assert Xor(True) is True
    assert Xor(False) is False
    assert Xor(True, True ) is False
    assert Xor(True, False) is True
    assert Xor(False, False) is False
    assert Xor(True, A) == ~A
    assert Xor(False, A) == A
    assert Xor(True, False, False) is True
    assert Xor(True, False, A) == ~A
    assert Xor(False, False, A) == A


def test_Not():

    raises(TypeError, lambda: Not(True, False))
    assert Not(True) is False
    assert Not(False) is True
    assert Not(0) is True
    assert Not(1) is False
    assert Not(2) is False


def test_Nand():

    assert Nand() is False
    assert Nand(A) == ~A
    assert Nand(True) is False
    assert Nand(False) is True
    assert Nand(True, True ) is False
    assert Nand(True, False) is True
    assert Nand(False, False) is True
    assert Nand(True, A) == ~A
    assert Nand(False, A) is True
    assert Nand(True, True, True) is False
    assert Nand(True, True, A) == ~A
    assert Nand(True, False, A) is True


def test_Nor():

    assert Nor() is True
    assert Nor(A) == ~A
    assert Nor(True) is False
    assert Nor(False) is True
    assert Nor(True, True ) is False
    assert Nor(True, False) is False
    assert Nor(False, False) is True
    assert Nor(True, A) is False
    assert Nor(False, A) == ~A
    assert Nor(True, True, True) is False
    assert Nor(True, True, A) is False
    assert Nor(True, False, A) is False


def test_Implies():

    raises(ValueError, lambda: Implies(A, B, C))
    assert Implies(True, True) is True
    assert Implies(True, False) is False
    assert Implies(False, True) is True
    assert Implies(False, False) is True
    assert Implies(0, A) is True
    assert Implies(1, 1) is True
    assert Implies(1, 0) is False
    assert A >> B == B << A


def test_Equivalent():

    assert Equivalent(A, B) == Equivalent(B, A) == Equivalent(A, B, A)
    assert Equivalent() is True
    assert Equivalent(A, A) == Equivalent(A) is True
    assert Equivalent(True, True) == Equivalent(False, False) is True
    assert Equivalent(True, False) == Equivalent(False, True) is False
    assert Equivalent(A, True) == A
    assert Equivalent(A, False) == Not(A)
    assert Equivalent(A, B, True) == A & B
    assert Equivalent(A, B, False) == ~A & ~B
    assert Equivalent(1, A) == A
    assert Equivalent(0, A) == Not(A)


def test_simplification():
    """
    Test working of simplification methods.
    """
    set1 = [[0, 0, 1], [0, 1, 1], [1, 0, 0], [1, 1, 0]]
    set2 = [[0, 0, 0], [0, 1, 0], [1, 0, 1], [1, 1, 1]]
    from sympy.abc import w, x, y, z
    assert SOPform('xyz', set1) == Or(And(Not(x), z), And(Not(z), x))
    assert Not(SOPform('xyz', set2)) == And(Or(Not(x), Not(z)), Or(x, z))
    assert POSform('xyz', set1 + set2) is True
    assert SOPform('xyz', set1 + set2) is True
    assert SOPform([Dummy(), Dummy(), Dummy()], set1 + set2) is True

    minterms = [[0, 0, 0, 1], [0, 0, 1, 1], [0, 1, 1, 1], [1, 0, 1, 1],
        [1, 1, 1, 1]]
    dontcares = [[0, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 1]]
    assert (
        SOPform('wxyz', minterms, dontcares) ==
        Or(And(Not(w), z), And(y, z)))
    assert POSform('wxyz', minterms, dontcares) == And(Or(Not(w), y), z)

    # test simplification
    ans = And(A, Or(B, C))
    assert simplify_logic('A & (B | C)') == ans
    assert simplify_logic('(A & B) | (A & C)') == ans
    assert simplify_logic(Implies(A, B)) == Or(Not(A), B)
    assert simplify_logic(Equivalent(A, B)) == \
           Or(And(A, B), And(Not(A), Not(B)))

    # check input
    ans = SOPform('xy', [[1, 0]])
    assert SOPform([x, y], [[1, 0]]) == ans
    assert POSform(['x', 'y'], [[1, 0]]) == ans

    raises(ValueError, lambda: SOPform('x', [[1]], [[1]]))
    assert SOPform('x', [[1]], [[0]]) is True
    assert SOPform('x', [[0]], [[1]]) is True
    assert SOPform('x', [], []) is False

    raises(ValueError, lambda: POSform('x', [[1]], [[1]]))
    assert POSform('x', [[1]], [[0]]) is True
    assert POSform('x', [[0]], [[1]]) is True
    assert POSform('x', [], []) is False

    #check working of simplify
    assert simplify('(A & B) | (A & C)') == sympify('And(A, Or(B, C))')
    assert simplify(And(x, Not(x))) == False
    assert simplify(Or(x, Not(x))) == True


def test_bool_equal():
    """
    Test working of bool_equal function.
    """

    minterms = [[0, 0, 0, 1], [0, 0, 1, 1], [0, 1, 1, 1], [1, 0, 1, 1],
        [1, 1, 1, 1]]
    from sympy.abc import a, b, c, x, y, z
    assert bool_equal(Not(Not(a)), a)
    assert bool_equal(SOPform(['w', 'x', 'y', 'z'], minterms),
        POSform(['w', 'x', 'y', 'z'], minterms))
    assert bool_equal(SOPform(['x', 'z', 'y'],[[1, 0, 1]]),
        SOPform(['a', 'b', 'c'],[[1, 0, 1]])) != False
    function1 = SOPform(['x','z','y'],[[1, 0, 1], [0, 0, 1]])
    function2 = SOPform(['a','b','c'],[[1, 0, 1], [1, 0, 0]])
    assert bool_equal(function1, function2, info=True) == \
        (function1, {y: a, z: b})


def test_bool_symbol():
    """Test that mixing symbols with boolean values
    works as expected"""

    assert And(A, True) == A
    assert And(A, True, True) == A
    assert And(A, False) is False
    assert And(A, True, False) is False
    assert Or(A, True) is True
    assert Or(A, False) == A


def test_subs():

    assert (A & B).subs(A, True) == B
    assert (A & B).subs(A, False) is False
    assert (A & B).subs(B, True) == A
    assert (A & B).subs(B, False) is False
    assert (A & B).subs({A: True, B: True}) is True
    assert (A | B).subs(A, True) is True
    assert (A | B).subs(A, False) == B
    assert (A | B).subs(B, True) is True
    assert (A | B).subs(B, False) == A
    assert (A | B).subs({A: True, B: True}) is True

"""
we test for axioms of boolean algebra
see http://en.wikipedia.org/wiki/Boolean_algebra_(structure)
"""


def test_commutative():
    """Test for commutativity of And and Or"""
    A, B = map(Boolean, symbols('A,B'))

    assert A & B == B & A
    assert A | B == B | A


def test_and_associativity():
    """Test for associativity of And"""

    assert (A & B) & C == A & (B & C)


def test_or_assicativity():

    assert ((A | B) | C) == (A | (B | C))


def test_double_negation():
    a = Boolean()
    assert ~(~a) == a


def test_De_Morgan():

    assert ~(A & B) == (~A) | (~B)
    assert ~(A | B) == (~A) & (~B)
    assert ~(A | B | C) == ~A & ~B & ~C

# test methods


def test_eliminate_implications():

    assert eliminate_implications(Implies(A, B, evaluate=False)) == (~A) | B
    assert eliminate_implications(
        A >> (C >> Not(B))) == Or(Or(Not(B), Not(C)), Not(A))


def test_conjuncts():
    assert conjuncts(A & B & C) == set([A, B, C])
    assert conjuncts((A | B) & C) == set([A | B, C])
    assert conjuncts(A) == set([A])
    assert conjuncts(True) == set([True])
    assert conjuncts(False) == set([False])


def test_disjuncts():
    assert disjuncts(A | B | C) == set([A, B, C])
    assert disjuncts((A | B) & C) == set([(A | B) & C])
    assert disjuncts(A) == set([A])
    assert disjuncts(True) == set([True])
    assert disjuncts(False) == set([False])


def test_distribute():

    assert distribute_and_over_or(Or(And(A, B), C)) == And(Or(A, C), Or(B, C))
    assert distribute_or_over_and(And(A, Or(B, C))) == Or(And(A, B), And(A, C))


def test_to_cnf():

    assert to_cnf(~(B | C)) == And(Not(B), Not(C))
    assert to_cnf((A & B) | C) == And(Or(A, C), Or(B, C))
    assert to_cnf(A >> B) == (~A) | B
    assert to_cnf(A >> (B & C)) == (~A | B) & (~A | C)
    assert to_cnf(A & (B | C) | ~A & (B | C), True) == B | C

    assert to_cnf(Equivalent(A, B)) == And(Or(A, Not(B)), Or(B, Not(A)))
    assert to_cnf(Equivalent(A, B & C)) == \
           (~A | B) & (~A | C) & (~B | ~C | A)
    assert to_cnf(Equivalent(A, B | C), True) == \
        And(Or(Not(B), A), Or(Not(C), A), Or(B, C, Not(A)))


def test_to_dnf():

    assert to_dnf(~(B | C)) == And(Not(B), Not(C))
    assert to_dnf(A & (B | C)) == Or(And(A, B), And(A, C))
    assert to_dnf(A >> B) == (~A) | B
    assert to_dnf(A >> (B & C)) == (~A) | (B & C)

    assert to_dnf(Equivalent(A, B), True) == \
           Or(And(A, B), And(Not(A), Not(B)))
    assert to_dnf(Equivalent(A, B & C), True) == \
           Or(And(A, B, C), And(Not(A), Not(B)), And(Not(A), Not(C)))


def test_to_int_repr():
    x, y, z = map(Boolean, symbols('x,y,z'))

    def sorted_recursive(arg):
        try:
            return sorted(sorted_recursive(x) for x in arg)
        except TypeError:  # arg is not a sequence
            return arg

    assert sorted_recursive(to_int_repr([x | y, z | x], [x, y, z])) == \
        sorted_recursive([[1, 2], [1, 3]])
    assert sorted_recursive(to_int_repr([x | y, z | ~x], [x, y, z])) == \
        sorted_recursive([[1, 2], [3, -1]])


def test_is_cnf():
    x, y, z = symbols('x,y,z')
    assert is_cnf(x) is True
    assert is_cnf(x | y | z) is True
    assert is_cnf(x & y & z) is True
    assert is_cnf((x | y) & z) is True
    assert is_cnf((x & y) | z) is False


def test_is_dnf():
    x, y, z = symbols('x,y,z')
    assert is_dnf(x) is True
    assert is_dnf(x | y | z) is True
    assert is_dnf(x & y & z) is True
    assert is_dnf((x & y) | z) is True
    assert is_dnf((x | y) & z) is False


def test_ITE():
    A, B, C = map(Boolean, symbols('A,B,C'))

    assert ITE(True, False, True) is False
    assert ITE(True, True, False) is True
    assert ITE(False, True, False) is False
    assert ITE(False, False, True) is True

    A = True
    assert ITE(A, B, C) == B
    A = False
    assert ITE(A, B, C) == C
    B = True
    assert ITE(And(A, B), B, C) == C
    assert ITE(Or(A, False), And(B, True), False) is False
