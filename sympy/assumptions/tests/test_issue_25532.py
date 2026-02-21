from sympy import Q, MatrixSymbol, Symbol, ask
from sympy.assumptions.satask import satask
from sympy.core.facts import InconsistentAssumptions
from sympy.testing.pytest import raises

def test_matrix_scalar_assumptions():
    X = MatrixSymbol('X', 2, 2)
    a = Symbol('a')

    assert ask(Q.positive(X)) is False
    assert ask(Q.negative(X)) is False
    assert ask(Q.prime(X)) is False
    assert ask(Q.composite(X)) is False
    assert ask(Q.even(X)) is False
    assert ask(Q.odd(X)) is False

    raises(InconsistentAssumptions, lambda: satask(Q.is_true(a), Q.is_true(a) & Q.positive(X)))
    raises(InconsistentAssumptions, lambda: satask(Q.is_true(a), Q.is_true(a) & Q.even(X)))

def test_satask_preserves_normal_behavior():
    x = Symbol('x')

    assert satask(Q.even(x), Q.integer(x) & Q.even(x)) is True
    assert satask(Q.odd(x), Q.even(x)) is False
    assert satask(Q.positive(x)) is None