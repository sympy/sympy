""" Test file for SymPy First Order Logic """

from sympy.utilities.pytest import raises

from sympy.logic.boolalg import (And, Implies, Or, Xor, true)
from sympy.logic.FOL import (AppliedFunction, AppliedPredicate, Constant,
    entails, Exists, fol_true, ForAll, Function, mgu, Predicate, resolve,
    standardize, to_cnf, to_dnf, to_pnf, to_snf)

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


def test_fol_true():
    P = Predicate('P')
    f = Function('f')
    _P = {0: False, 'default': True}
    _f = {0: 1, 'default': 0}
    assert fol_true(P(0)) is None
    assert fol_true(P(f(0)), {P:_P, f:_f}) is True
    assert fol_true(P(f(X)), {X:1, P:_P, f:_f}) is False
    assert fol_true(ForAll(X, P(X)), {P:_P}) is None
    assert fol_true(ForAll(X, Exists(Y, P(X) & P(Y))),
        {X:[0, 1], P:_P}) is None

    from operator import gt, lt, eq
    GT = Predicate('GT')
    LT = Predicate('LT')
    EQ = Predicate('EQ')
    add1 = Function('add1')
    _add1 = lambda x: (x + 1) % 6
    domain = range(6)
    assert fol_true(LT(add1(X), X), {X:1, LT:lt, add1:_add1}) is False
    assert fol_true(GT(X, Y) & GT(Y, Z) >> GT(X, Z),
        {X:3, Y:2, Z:1, GT:gt}) is True
    assert fol_true(ForAll(X, GT(add1(X), X)),
        {X:domain, GT:gt, add1:_add1}) is False
    assert fol_true(Exists(X, EQ(_add1(X), X),),
        {X:domain, EQ:eq, add1:_add1}) is False
    assert fol_true(ForAll((X, Y, Z), (GT(X, Y) & GT(Y, Z)) >> GT(X, Z)),
        {X:domain, Y:domain, Z:domain, GT:gt}) is True


def test_standardize():
    from sympy.core.symbol import Symbol
    X0 = Symbol('X0')
    Y0 = Symbol('Y0')
    P = Predicate('P')
    Q = Predicate('Q')
    R = Predicate('R')
    assert standardize(ForAll(X, P(X)) >> ForAll(X, Q(X))) == \
                ForAll(X, P(X)) >> ForAll(X0, Q(X0))
    assert standardize(ForAll(X, P(X)) & ForAll(X, Q(X))) == \
                ForAll(X, P(X)) & ForAll(X, Q(X))
    assert standardize(Exists(X, P(X)) | Exists(X, Q(X))) == \
                Exists(X, P(X)) | Exists(X, Q(X))
    assert standardize(ForAll(X, P(X) & Q(X)) | ForAll(X, P(X) >> Q(X))) == \
                ForAll(X, P(X) & Q(X)) | ForAll(X0, P(X0) >> Q(X0))
    assert standardize(Exists(X, P(X) & Q(X)) & Exists(X, P(X) >> Q(X))) == \
                Exists(X, P(X) & Q(X)) & Exists(X0, P(X0) >> Q(X0))
    assert standardize(Exists(X, P(X) >> ForAll(Y, R(Y))) |
                Exists((X, Y), P(X, Y))) == Exists(X, P(X) >>
                ForAll(Y0, R(Y0))) | Exists((X, Y), P(X, Y))


def test_to_pnf():
    P = Predicate('P')
    Q = Predicate('Q')
    assert to_pnf(True) is true
    assert to_pnf(P(X, Y)) == P(X, Y)
    assert to_pnf(~ForAll(X, P(X))) == Exists(X, ~P(X))
    assert to_pnf(~Exists(X, P(X))) == ForAll(X, ~P(X))
    assert to_pnf(Exists(X, ForAll(Y, P(X, Y)))) == \
                Exists(X, ForAll(Y, P(X, Y)))
    assert to_pnf(ForAll((X, Y), ~(P(X) >> Q(Y)))) == \
                ForAll((X, Y), P(X) & ~Q(Y))
    assert to_pnf(Exists(X, P(X)) & Exists(Y, Q(Y))) == \
                Exists((X, Y), P(X) & Q(Y))
    assert to_pnf(ForAll(X, P(X)) >> ForAll(Y, Q(Y))) == \
                ForAll(Y, Exists(X, ~P(X) | Q(Y)))
    assert to_pnf(ForAll(X, P(X)) & ForAll(X, Q(X))) == \
                ForAll(X, P(X) & Q(X))
    assert to_pnf(Exists(X, P(X)) | Exists(X, Q(X))) == \
                Exists(X, P(X) | Q(X))


def test_to_snf():
    from sympy.abc import W
    P = Predicate('P')
    Q = Predicate('Q')
    f0 = Function('f0')
    f1 = Function('f1')
    c0 = Constant('c0')
    assert to_snf(ForAll((X, Y), P(X) >> Q(Y))) == \
                ForAll((X, Y), ~P(X) | Q(Y))
    assert to_snf(Exists(X, ForAll(Y, P(X) | Q(Y)))) == \
                ForAll(Y, P(c0) | Q(Y))
    assert to_snf(ForAll(X, Exists(Y, P(X) & Q(Y)))) == \
                ForAll(X, P(X) & Q(f0(X)))
    assert to_snf(ForAll(W, Exists(X, ForAll(Y, Exists(Z, P(W, X) >>
        Q(Y, Z)))))) == ForAll((W, Y), ~P(W, f0(W)) | Q(Y, f1(W, Y)))


def test_to_cnf():
    P = Predicate('P')
    Q = Predicate('Q')
    assert to_cnf(to_cnf((P(X) & Q(X)) | (P(Y) & Q(Y)))) == \
        (P(X) | P(Y)) & (P(X) | Q(Y)) & (Q(X) | P(Y)) & (Q(X) | Q(Y))
    assert to_cnf(ForAll(X, P(X) | ForAll(Y, Q(Y)))) == Or(P(X), Q(Y))


def test_to_dnf():
    P = Predicate('P')
    Q = Predicate('Q')
    assert to_dnf((P(X) | Q(X)) & (P(Y) | Q(Y))) == \
        (P(X) & P(Y)) | (P(X) & Q(Y)) | (Q(X) & P(Y)) | (Q(X) & Q(Y))
    assert to_dnf(ForAll(X, P(X) | ForAll(Y, Q(Y)))) == Or(P(X), Q(Y))


def test_mgu():
    P = Predicate('P')
    Q = Predicate('Q')
    f = Function('f')
    g = Function('g')
    a = Constant('a')
    b = Constant('b')
    assert mgu(P(X), Q(X)) is False
    assert mgu(P(X), P(X, Y)) is False
    assert mgu(P(X), P(X)) == {true: true}
    assert mgu(P(X, X), P(Y, f(Y))) is False
    assert mgu(P(X, Y), P(f(X), Z)) is False
    assert mgu(P(f(a)), P(f(b))) is False
    assert mgu(P(X, f(a)), P(Y, X)) == {X: f(a), Y: f(a)}
    assert mgu(P(a, X, f(g(Z))), P(Z, f(Y), f(Y))) == {X: f(g(a)), Z: a, Y: g(a)}


def test_resolution():
    P = Predicate('P')
    Q = Predicate('Q')
    a = Constant('a')
    assert resolve(And(Or(P(X), Q(X)), Or(~P(X), Q(X)))) is True
    assert resolve(And(P(X) >> Q(X), P(a), ~Q(a))) is False


def test_entails():
    P = Predicate('P')
    a = Constant('a')
    b = Constant('b')
    c = Constant('c')
    formula_set = [
        ForAll(X, ForAll(Y, (P(X, Y) & P(Y, Z)) >> P(X, Z))),
        P(a, b),
        P(b, c)
    ]
    assert entails(P(a, c), formula_set) is True
