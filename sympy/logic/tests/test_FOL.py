""" Test file for SymPy First Order Logic """

from sympy.abc import X, Y, Z
from sympy.utilities.pytest import raises

from sympy.logic.boolalg import (And, Implies, Or, Xor)
from sympy.logic.FOL import (AppliedFunction, AppliedPredicate, Exists,
    ForAll, Function, Predicate)


def test_Predicate():
    A = Predicate('A')
    B = Predicate('B')
    P = A(X, Y)
    Q = B(Z)
    assert isinstance(A, Predicate)
    assert isinstance(P, AppliedPredicate)
    assert P.name == 'A'
    assert P.func == A
    raises(ValueError, lambda: A())
    assert not A(X, Y) == A(Y, X)
    assert not A(X, Y) == A(X, Y, Z)
    assert P & Q == And(A(X, Y), B(Z))
    assert Or(P) == A(X, Y)


def test_Function():
    f = Function('f')
    g = Function('g')
    fx = f(X)
    gyz = g(Y, Z)
    assert isinstance(f, Function)
    assert isinstance(fx, AppliedFunction)
    assert fx.name == 'f'
    assert fx.func == f
    raises(ValueError, lambda: f())
    assert not f(X, Y) == f(Y, X)
    assert not f(X, Y) == f(X, Y, Z)
    assert fx | gyz == Or(f(X), g(Y, Z))
    assert And(fx) == f(X)


def test_ForAll():
    A = Predicate('A')
    B = Predicate('B')
    assert ForAll(X, ForAll(Y, A(X) & B(Y))) == ForAll((X, Y), A(X) & B(Y))
    raises(ValueError, lambda: ForAll(X, A(X) | ForAll(X, B(X))))
    ForAll((), A(X)) == A(X)
    ForAll(Y, A(X)) == A(X)
    ForAll((X, Y), A(X)) == ForAll(X, A(X))
    ForAll((X, Y, Z), A(X) >> B(Y, Z)).vars == (X, Y, Z)
    ForAll((X, Y, Z), A(X) >> B(Y, Z)).expr == Implies(A(X), B(Y, Z))
    ForAll((X, Y), ForAll(Z, A(X, Y) ^ B(Y, Z))).vars == (X, Y, Z)
    ForAll((X, Y), ForAll(Z, A(X, Y) ^ B(Y, Z))).expr == Xor(A(X, Y), B(Y, Z))


def test_Exists():
    A = Predicate('A')
    B = Predicate('B')
    assert Exists(X, Exists(Y, A(X) & B(Y))) == Exists((X, Y), A(X) & B(Y))
    raises(ValueError, lambda: Exists(X, A(X) | Exists(X, B(X))))
    Exists((), A(X)) == A(X)
    Exists(Y, A(X)) == A(X)
    Exists((X, Y), A(X)) == Exists(X, A(X))
    Exists((X, Y, Z), A(X) >> B(Y, Z)).vars == (X, Y, Z)
    Exists((X, Y, Z), A(X) >> B(Y, Z)).expr == Implies(A(X), B(Y, Z))
    Exists((X, Y), Exists(Z, A(X, Y) ^ B(Y, Z))).vars == (X, Y, Z)
    Exists((X, Y), Exists(Z, A(X, Y) ^ B(Y, Z))).expr == Xor(A(X, Y), B(Y, Z))


def test_to_pnf():
    pass