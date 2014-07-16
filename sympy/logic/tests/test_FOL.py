""" Test file for SymPy First Order Logic """

from sympy.utilities.pytest import raises

from sympy.logic.boolalg import (And, Implies, Or, Xor, true)
from sympy.logic.FOL import (AppliedFunction, AppliedPredicate, Exists,
    ForAll, Function, mgu, Predicate, standardize, to_pnf, to_snf)

from sympy.abc import X, Y, Z

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


def test_standardize():
    from sympy.core.symbol import Symbol
    X0 = Symbol('X0')
    P = Predicate('P')
    Q = Predicate('Q')
    assert standardize(ForAll(X, P(X)) >> ForAll(X, Q(X))) == \
                ForAll(X, P(X)) >> ForAll(X0, Q(X0))


def test_to_pnf():
    P = Predicate('P')
    Q = Predicate('Q')
    assert to_pnf(P(X, Y)) == P(X, Y)
    assert to_pnf(Exists(X, ForAll(Y, P(X, Y)))) == \
                Exists(X, ForAll(Y, P(X, Y)))
    assert to_pnf(ForAll((X, Y), ~(P(X) >> Q(Y)))) == \
                ForAll((X, Y), P(X) & ~Q(Y))
    assert to_pnf(Exists(X, P(X)) & Exists(Y, Q(Y))) == \
                Exists((X, Y), P(X) & Q(Y))
    assert to_pnf(ForAll(X, P(X)) >> ForAll(Y, Q(Y))) == \
                ForAll(Y, Exists(X, ~P(X) | Q(Y)))


def test_to_snf():
    from sympy.abc import W
    P = Predicate('P')
    Q = Predicate('Q')
    f0 = Function('f0')
    f1 = Function('f1')
    assert to_snf(ForAll((X, Y), P(X) >> Q(Y))) == \
                ForAll((X, Y), ~P(X) | Q(Y))
    assert to_snf(Exists(X, ForAll(Y, P(X) | Q(Y)))) == \
                ForAll(Y, P(X) | Q(Y))
    assert to_snf(ForAll(X, Exists(Y, P(X) & Q(Y)))) == \
                ForAll(X, P(X) & Q(f0(X)))
    assert to_snf(ForAll(W, Exists(X, ForAll(Y, Exists(Z, P(W, X) >>
        Q(Y, Z)))))) == ForAll((W, Y), ~P(W, f0(W)) | Q(Y, f1(W, Y)))


def test_mgu():
    P = Predicate('A')
    Q = Predicate('B')
    f = Function('f')
    g = Function('g')
    assert mgu(P(X), Q(X)) is False
    assert mgu(P(X), P(X, Y)) is False
    assert mgu(P(X), P(X)) == {true: true}
    assert mgu(P(X, X), P(Y, f(Y))) is False
    assert mgu(P('a', X, f(g(Z))), P(Z, f(Y), f(Y))) == {X: f(Y), Z: 'a', Y: g('a')}


def test_resolution():
    pass
