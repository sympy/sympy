from sympy import (symbols, sympify, Dummy, simplify, Equality, S, Interval,
                   oo, EmptySet)
from sympy.logic.boolalg import (
    And, Boolean, Equivalent, ITE, Implies, Nand, Nor, Not, Or, POSform,
    SOPform, Xor, conjuncts, disjuncts, distribute_or_over_and,
    distribute_and_over_or, eliminate_implications, is_cnf, is_dnf,
    simplify_logic, to_cnf, to_dnf, to_int_repr, bool_map, true, false,
    BooleanAtom
)
from sympy.utilities.pytest import raises, XFAIL
from sympy.utilities import cartes


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

    assert And() is true
    assert And(A) == A
    assert And(True) is true
    assert And(False) is false
    assert And(True, True ) is true
    assert And(True, False) is false
    assert And(False, False) is false
    assert And(True, A) == A
    assert And(False, A) is false
    assert And(True, True, True) is true
    assert And(True, True, A) == A
    assert And(True, False, A) is false
    assert And(2, A) == A
    assert And(2, 3) is true


def test_Or():

    assert Or() is false
    assert Or(A) == A
    assert Or(True) is true
    assert Or(False) is false
    assert Or(True, True ) is true
    assert Or(True, False) is true
    assert Or(False, False) is false
    assert Or(True, A) is true
    assert Or(False, A) == A
    assert Or(True, False, False) is true
    assert Or(True, False, A) is true
    assert Or(False, False, A) == A
    assert Or(2, A) is true


def test_Xor():

    assert Xor() is false
    assert Xor(A) == A
    assert Xor(A, A) is false
    assert Xor(True, A, A) is true
    assert Xor(A, A, A, A, A) == A
    assert Xor(True, False, False, A, B) == And(Or(A, Not(B)), Or(B, Not(A)))
    assert Xor(True) is true
    assert Xor(False) is false
    assert Xor(True, True ) is false
    assert Xor(True, False) is true
    assert Xor(False, False) is false
    assert Xor(True, A) == ~A
    assert Xor(False, A) == A
    assert Xor(True, False, False) is true
    assert Xor(True, False, A) == ~A
    assert Xor(False, False, A) == A


def test_Not():

    raises(TypeError, lambda: Not(True, False))
    assert Not(True) is false
    assert Not(False) is true
    assert Not(0) is true
    assert Not(1) is false
    assert Not(2) is false


def test_Nand():

    assert Nand() is false
    assert Nand(A) == ~A
    assert Nand(True) is false
    assert Nand(False) is true
    assert Nand(True, True ) is false
    assert Nand(True, False) is true
    assert Nand(False, False) is true
    assert Nand(True, A) == ~A
    assert Nand(False, A) is true
    assert Nand(True, True, True) is false
    assert Nand(True, True, A) == ~A
    assert Nand(True, False, A) is true


def test_Nor():

    assert Nor() is true
    assert Nor(A) == ~A
    assert Nor(True) is false
    assert Nor(False) is true
    assert Nor(True, True ) is false
    assert Nor(True, False) is false
    assert Nor(False, False) is true
    assert Nor(True, A) is false
    assert Nor(False, A) == ~A
    assert Nor(True, True, True) is false
    assert Nor(True, True, A) is false
    assert Nor(True, False, A) is false


def test_Implies():

    raises(ValueError, lambda: Implies(A, B, C))
    assert Implies(True, True) is true
    assert Implies(True, False) is false
    assert Implies(False, True) is true
    assert Implies(False, False) is true
    assert Implies(0, A) is true
    assert Implies(1, 1) is true
    assert Implies(1, 0) is false
    assert A >> B == B << A


def test_Equivalent():

    assert Equivalent(A, B) == Equivalent(B, A) == Equivalent(A, B, A)
    assert Equivalent() is true
    assert Equivalent(A, A) == Equivalent(A) is true
    assert Equivalent(True, True) == Equivalent(False, False) is true
    assert Equivalent(True, False) == Equivalent(False, True) is false
    assert Equivalent(A, True) == A
    assert Equivalent(A, False) == Not(A)
    assert Equivalent(A, B, True) == A & B
    assert Equivalent(A, B, False) == ~A & ~B
    assert Equivalent(1, A) == A
    assert Equivalent(0, A) == Not(A)
    assert Equivalent(A, Equivalent(B, C)) != Equivalent(Equivalent(A, B), C)


def test_simplification():
    """
    Test working of simplification methods.
    """
    set1 = [[0, 0, 1], [0, 1, 1], [1, 0, 0], [1, 1, 0]]
    set2 = [[0, 0, 0], [0, 1, 0], [1, 0, 1], [1, 1, 1]]
    from sympy.abc import w, x, y, z
    assert SOPform('xyz', set1) == Or(And(Not(x), z), And(Not(z), x))
    assert Not(SOPform('xyz', set2)) == And(Or(Not(x), Not(z)), Or(x, z))
    assert POSform('xyz', set1 + set2) is true
    assert SOPform('xyz', set1 + set2) is true
    assert SOPform([Dummy(), Dummy(), Dummy()], set1 + set2) is true

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
    assert simplify_logic(And(Equality(A, 2), C)) == And(Equality(A, 2), C)
    assert simplify_logic(And(Equality(A, 2), A)) == And(Equality(A, 2), A)
    assert simplify_logic(And(Equality(A, B), C)) == And(Equality(A, B), C)
    assert simplify_logic(Or(And(Equality(A, 3), B), And(Equality(A, 3), C))) \
           == And(Equality(A, 3), Or(B, C))
    assert simplify_logic(And(A, x**2-x)) == And(A, x*(x-1))
    assert simplify_logic(And(A, x**2-x), deep=False) == And(A, x**2-x)

    # check input
    ans = SOPform('xy', [[1, 0]])
    assert SOPform([x, y], [[1, 0]]) == ans
    assert POSform(['x', 'y'], [[1, 0]]) == ans

    raises(ValueError, lambda: SOPform('x', [[1]], [[1]]))
    assert SOPform('x', [[1]], [[0]]) is true
    assert SOPform('x', [[0]], [[1]]) is true
    assert SOPform('x', [], []) is false

    raises(ValueError, lambda: POSform('x', [[1]], [[1]]))
    assert POSform('x', [[1]], [[0]]) is true
    assert POSform('x', [[0]], [[1]]) is true
    assert POSform('x', [], []) is false

    #check working of simplify
    assert simplify('(A & B) | (A & C)') == sympify('And(A, Or(B, C))')
    assert simplify(And(x, Not(x))) == False
    assert simplify(Or(x, Not(x))) == True


def test_bool_map():
    """
    Test working of bool_map function.
    """

    minterms = [[0, 0, 0, 1], [0, 0, 1, 1], [0, 1, 1, 1], [1, 0, 1, 1],
        [1, 1, 1, 1]]
    from sympy.abc import a, b, c, w, x, y, z
    assert bool_map(Not(Not(a)), a) == (a, {a: a})
    assert bool_map(SOPform(['w', 'x', 'y', 'z'], minterms),
        POSform(['w', 'x', 'y', 'z'], minterms)) == \
        (And(Or(Not(w), y), Or(Not(x), y), z), {x: x, w: w, z: z, y: y})
    assert bool_map(SOPform(['x', 'z', 'y'],[[1, 0, 1]]),
        SOPform(['a', 'b', 'c'],[[1, 0, 1]])) != False
    function1 = SOPform(['x','z','y'],[[1, 0, 1], [0, 0, 1]])
    function2 = SOPform(['a','b','c'],[[1, 0, 1], [1, 0, 0]])
    assert bool_map(function1, function2) == \
        (function1, {y: a, z: b})


def test_bool_symbol():
    """Test that mixing symbols with boolean values
    works as expected"""

    assert And(A, True) == A
    assert And(A, True, True) == A
    assert And(A, False) is false
    assert And(A, True, False) is false
    assert Or(A, True) is true
    assert Or(A, False) == A


def test_subs():

    assert (A & B).subs(A, True) == B
    assert (A & B).subs(A, False) is false
    assert (A & B).subs(B, True) == A
    assert (A & B).subs(B, False) is false
    assert (A & B).subs({A: True, B: True}) is true
    assert (A | B).subs(A, True) is true
    assert (A | B).subs(A, False) == B
    assert (A | B).subs(B, True) is true
    assert (A | B).subs(B, False) == A
    assert (A | B).subs({A: True, B: True}) is true

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
    from sympy.abc import A, B, C, D
    assert eliminate_implications(Implies(A, B, evaluate=False)) == (~A) | B
    assert eliminate_implications(
        A >> (C >> Not(B))) == Or(Or(Not(B), Not(C)), Not(A))
    assert eliminate_implications(Equivalent(A, B, C, D)) == \
        (~A | B) & (~B | C) & (~C | D) & (~D | A)


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

    assert ITE(True, False, True) is false
    assert ITE(True, True, False) is true
    assert ITE(False, True, False) is false
    assert ITE(False, False, True) is true

    A = True
    assert ITE(A, B, C) == B
    A = False
    assert ITE(A, B, C) == C
    B = True
    assert ITE(And(A, B), B, C) == C
    assert ITE(Or(A, False), And(B, True), False) is false

def test_operators():
    # Mostly test __and__, __rand__, and so on
    assert True & A == A & True == A
    assert False & A == A & False == False
    assert A & B == And(A, B)
    assert True | A == A | True == True
    assert False | A == A | False == A
    assert A | B == Or(A, B)
    assert ~A == Not(A)
    assert True >> A == A << True == A
    assert False >> A == A << False == True
    assert A >> True == True << A == True
    assert A >> False == False << A == ~A
    assert A >> B == B << A == Implies(A, B)
    assert True ^ A == A ^ True == ~A
    assert False ^ A == A ^ False == A
    assert A ^ B == Xor(A, B)

def test_true_false():
    x = symbols('x')

    assert true is S.true
    assert false is S.false
    assert true is not True
    assert false is not False
    assert true
    assert not false
    assert true == True
    assert false == False
    assert not (true == False)
    assert not (false == True)
    assert not (true == false)

    assert hash(true) == hash(True)
    assert hash(false) == hash(False)
    assert len(set([true, True])) == len(set([false, False])) == 1

    assert isinstance(true, BooleanAtom)
    assert isinstance(false, BooleanAtom)
    # We don't want to subclass from bool, because bool subclasses from
    # int. But operators like &, |, ^, <<, >>, and ~ act differently on 0 and
    # 1 then we want them to on true and false.  See the docstrings of the
    # various And, Or, etc. functions for examples.
    assert not isinstance(true, bool)
    assert not isinstance(false, bool)

    # Note: using 'is' comparison is important here. We want these to return
    # true and false, not True and False

    assert Not(true) is false
    assert Not(True) is false
    assert Not(false) is true
    assert Not(False) is true
    assert ~true is false
    assert ~false is true

    for T, F in cartes([True, true], [False, false]):
        assert And(T, F) is false
        assert And(F, T) is false
        assert And(F, F) is false
        assert And(T, T) is true
        assert And(T, x) == x
        assert And(F, x) is false
        if not (T is True and F is False):
            assert T & F is false
            assert F & T is false
        if not F is False:
            assert F & F is false
        if not T is True:
            assert T & T is true

        assert Or(T, F) is true
        assert Or(F, T) is true
        assert Or(F, F) is false
        assert Or(T, T) is true
        assert Or(T, x) is true
        assert Or(F, x) == x
        if not (T is True and F is False):
            assert T | F is true
            assert F | T is true
        if not F is False:
            assert F | F is false
        if not T is True:
            assert T | T is true

        assert Xor(T, F) is true
        assert Xor(F, T) is true
        assert Xor(F, F) is false
        assert Xor(T, T) is false
        assert Xor(T, x) == ~x
        assert Xor(F, x) == x
        if not (T is True and F is False):
            assert T ^ F is true
            assert F ^ T is true
        if not F is False:
            assert F ^ F is false
        if not T is True:
            assert T ^ T is false

        assert Nand(T, F) is true
        assert Nand(F, T) is true
        assert Nand(F, F) is true
        assert Nand(T, T) is false
        assert Nand(T, x) == ~x
        assert Nand(F, x) is true

        assert Nor(T, F) is false
        assert Nor(F, T) is false
        assert Nor(F, F) is true
        assert Nor(T, T) is false
        assert Nor(T, x) is false
        assert Nor(F, x) == ~x

        assert Implies(T, F) is false
        assert Implies(F, T) is true
        assert Implies(F, F) is true
        assert Implies(T, T) is true
        assert Implies(T, x) == x
        assert Implies(F, x) is true
        assert Implies(x, T) is true
        assert Implies(x, F) == ~x
        if not (T is True and F is False):
            assert T >> F is false
            assert F << T is false
            assert F >> T is true
            assert T << F is true
        if not F is False:
            assert F >> F is true
            assert F << F is true
        if not T is True:
            assert T >> T is true
            assert T << T is true

        assert Equivalent(T, F) is false
        assert Equivalent(F, T) is false
        assert Equivalent(F, F) is true
        assert Equivalent(T, T) is true
        assert Equivalent(T, x) == x
        assert Equivalent(F, x) == ~x
        assert Equivalent(x, T) == x
        assert Equivalent(x, F) == ~x

        assert ITE(T, T, T) is true
        assert ITE(T, T, F) is true
        assert ITE(T, F, T) is false
        assert ITE(T, F, F) is false
        assert ITE(F, T, T) is true
        assert ITE(F, T, F) is false
        assert ITE(F, F, T) is true
        assert ITE(F, F, F) is false


def test_bool_as_set():
    x = symbols('x')

    assert And(x <= 2, x >= -2).as_set() == Interval(-2, 2)
    assert Or(x >= 2, x <= -2).as_set() == Interval(-oo, -2) + Interval(2, oo)
    assert Not(x > 2).as_set() == Interval(-oo, 2)
    assert true.as_set() == S.UniversalSet
    assert false.as_set() == EmptySet()


@XFAIL
def test_multivariate_bool_as_set():
    x, y = symbols('x,y')

    assert And(x >= 0, y >= 0).as_set() == Interval(0, oo)*Interval(0, oo)
    assert Or(x >= 0, y >= 0).as_set() == S.UniversalSet - \
        Interval(-oo, 0, True, True)*Interval(-oo, 0, True, True)
