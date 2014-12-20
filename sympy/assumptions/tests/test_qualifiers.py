from sympy import sympify, Q, S, Symbol
from sympy.assumptions.qualifiers import Qualified, ContradictionError
from sympy.utilities.pytest import raises


def test_equals():
    a = Symbol("a")
    x = Symbol("x")
    assert Qualified(x, a) == Qualified(x, a)


def test_nontrivial_axioms():
    e = sympify("x")
    p1 = Qualified(e, Q.positive(S("x")))
    p2 = Qualified(e, Q.nonnegative(S("x")))
    assert p1 == p1
    assert p2 == p2
    assert p1 != p2


def test_two():
    two_a = Qualified(2, sympify("2==2"))
    two_b = Qualified(2, sympify("2!=3"))
    c = two_a + two_b
    assert (two_a + two_b).unqualified == 4


def test_large_2():
    raises(ContradictionError, lambda: Qualified("2+2", sympify("2+2==5")))


def test_int_cant_be_odd_and_even():
    x = Symbol("x")

    theory = Q.odd(x) ^ Q.even(x)
    even_x = Qualified(x, Q.even(x))
    odd_x = Qualified(x, Q.odd(x))

    raises(ContradictionError, lambda: Qualified(even_x + odd_x, theory))
    raises(ContradictionError, lambda: Qualified(even_x * odd_x, theory))

