from sympy import (S, symbols, FiniteSet, Eq, Matrix, MatrixSymbol, Float, And)
from sympy.stats import DiscreteMarkovChain, P, TransitionMatrixOf
from sympy.stats.rv import RandomIndexedSymbol
from sympy.stats.symbolic_probability import Probability
from sympy.utilities.pytest import raises

def test_DiscreteMarkovChain():

    # pass only the name
    X = DiscreteMarkovChain("X")
    assert X.state_space == S.Reals
    assert X.index_set == S.Naturals0
    assert X.transition_probabilities == None
    t = symbols('t', positive=True, integer=True)
    assert isinstance(X[t], RandomIndexedSymbol)

    # pass name and state_space
    Y = DiscreteMarkovChain("Y", [1, 2, 3])
    assert Y.transition_probabilities == None
    assert Y.state_space == FiniteSet(1, 2, 3)
    assert P(Eq(Y[2], 1), Eq(Y[0], 2)) == Probability(Eq(Y[2], 1), Eq(Y[0], 2))

    # pass name, state_space and transition_probabilities
    T = Matrix([[0.5, 0.2, 0.3],[0.2, 0.5, 0.3],[0.2, 0.3, 0.5]])
    TS = MatrixSymbol('T', 3, 3)
    Y = DiscreteMarkovChain("Y", [0, 1, 2], T)
    YS = DiscreteMarkovChain("Y", [0, 1, 2], TS)
    assert P(Eq(Y[3], 2), Eq(Y[1], 1)).round(2) == Float(0.36, 2)
    assert str(P(Eq(YS[3], 2), Eq(YS[1], 1))) == \
        "T[0, 2]*T[1, 0] + T[1, 1]*T[1, 2] + T[1, 2]*T[2, 2]"
    TO = Matrix([[0.25, 0.75, 0],[0, 0.25, 0.75],[0.75, 0, 0.25]])
    assert P(Eq(Y[3], 2), Eq(Y[1], 1) & TransitionMatrixOf(Y, TO)).round(3) == Float(0.375, 3)
    TSO = MatrixSymbol('T', 4, 4)
    raises(ValueError, lambda: str(P(Eq(YS[3], 2), Eq(YS[1], 1) & TransitionMatrixOf(YS, TSO))))

    # extended tests for probability queries
    TO1 = Matrix([[S(1)/4, S(3)/4, 0],[S(1)/3, S(1)/3, S(1)/3],[0, S(1)/4, S(3)/4]])
    assert P(And(Eq(Y[2], 1), Eq(Y[1], 1), Eq(Y[0], 0)),
            Eq(Probability(Eq(Y[0], 0)), S(1)/4) & TransitionMatrixOf(Y, TO1)) == S(1)/16
    assert P(And(Eq(Y[2], 1), Eq(Y[1], 1), Eq(Y[0], 0)), TransitionMatrixOf(Y, TO1)) == \
            Probability(Eq(Y[0], 0))/4
