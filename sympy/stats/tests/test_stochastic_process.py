from sympy import (S, symbols, FiniteSet, Eq, Matrix, MatrixSymbol, Float, And,
                   ImmutableMatrix, Ne, Lt, Gt, exp, Not, Rational)
from sympy.stats import (DiscreteMarkovChain, P, TransitionMatrixOf, E,
                         StochasticStateSpaceOf, variance, ContinuousMarkovChain)
from sympy.stats.joint_rv import JointDistribution
from sympy.stats.rv import RandomIndexedSymbol
from sympy.stats.symbolic_probability import Probability, Expectation
from sympy.utilities.pytest import raises


def test_DiscreteMarkovChain():

    # pass only the name
    X = DiscreteMarkovChain("X")
    assert X.state_space == S.Reals
    assert X.index_set == S.Naturals0
    assert X.transition_probabilities == None
    t = symbols('t', positive=True, integer=True)
    assert isinstance(X[t], RandomIndexedSymbol)
    assert E(X[0]) == Expectation(X[0])
    raises(TypeError, lambda: DiscreteMarkovChain(1))
    raises(NotImplementedError, lambda: X(t))

    # pass name and state_space
    Y = DiscreteMarkovChain("Y", [1, 2, 3])
    assert Y.transition_probabilities == None
    assert Y.state_space == FiniteSet(1, 2, 3)
    assert P(Eq(Y[2], 1), Eq(Y[0], 2)) == Probability(Eq(Y[2], 1), Eq(Y[0], 2))
    assert E(X[0]) == Expectation(X[0])
    raises(TypeError, lambda: DiscreteMarkovChain("Y", dict((1, 1))))

    # pass name, state_space and transition_probabilities
    T = Matrix([[0.5, 0.2, 0.3],[0.2, 0.5, 0.3],[0.2, 0.3, 0.5]])
    TS = MatrixSymbol('T', 3, 3)
    Y = DiscreteMarkovChain("Y", [0, 1, 2], T)
    YS = DiscreteMarkovChain("Y", [0, 1, 2], TS)
    assert YS._transient2transient() == None
    assert YS._transient2absorbing() == None
    assert Y.joint_distribution(1, Y[2], 3) == JointDistribution(Y[1], Y[2], Y[3])
    raises(ValueError, lambda: Y.joint_distribution(Y[1].symbol, Y[2].symbol))
    assert P(Eq(Y[3], 2), Eq(Y[1], 1)).round(2) == Float(0.36, 2)
    assert str(P(Eq(YS[3], 2), Eq(YS[1], 1))) == \
        "T[0, 2]*T[1, 0] + T[1, 1]*T[1, 2] + T[1, 2]*T[2, 2]"
    assert P(Eq(YS[1], 1), Eq(YS[2], 2)) == Probability(Eq(YS[1], 1))
    assert P(Eq(YS[3], 3), Eq(YS[1], 1)) is S.Zero
    TO = Matrix([[0.25, 0.75, 0],[0, 0.25, 0.75],[0.75, 0, 0.25]])
    assert P(Eq(Y[3], 2), Eq(Y[1], 1) & TransitionMatrixOf(Y, TO)).round(3) == Float(0.375, 3)
    assert E(Y[3], evaluate=False) == Expectation(Y[3])
    assert E(Y[3], Eq(Y[2], 1)).round(2) == Float(1.1, 3)
    TSO = MatrixSymbol('T', 4, 4)
    raises(ValueError, lambda: str(P(Eq(YS[3], 2), Eq(YS[1], 1) & TransitionMatrixOf(YS, TSO))))
    raises(TypeError, lambda: DiscreteMarkovChain("Z", [0, 1, 2], symbols('M')))
    raises(ValueError, lambda: DiscreteMarkovChain("Z", [0, 1, 2], MatrixSymbol('T', 3, 4)))
    raises(ValueError, lambda: E(Y[3], Eq(Y[2], 6)))
    raises(ValueError, lambda: E(Y[2], Eq(Y[3], 1)))


    # extended tests for probability queries
    TO1 = Matrix([[Rational(1, 4), Rational(3, 4), 0],[Rational(1, 3), Rational(1, 3), Rational(1, 3)],[0, Rational(1, 4), Rational(3, 4)]])
    assert P(And(Eq(Y[2], 1), Eq(Y[1], 1), Eq(Y[0], 0)),
            Eq(Probability(Eq(Y[0], 0)), Rational(1, 4)) & TransitionMatrixOf(Y, TO1)) == Rational(1, 16)
    assert P(And(Eq(Y[2], 1), Eq(Y[1], 1), Eq(Y[0], 0)), TransitionMatrixOf(Y, TO1)) == \
            Probability(Eq(Y[0], 0))/4
    assert P(Lt(X[1], 2) & Gt(X[1], 0), Eq(X[0], 2) &
        StochasticStateSpaceOf(X, [0, 1, 2]) & TransitionMatrixOf(X, TO1)) == Rational(1, 4)
    assert P(Ne(X[1], 2) & Ne(X[1], 1), Eq(X[0], 2) &
        StochasticStateSpaceOf(X, [0, 1, 2]) & TransitionMatrixOf(X, TO1)) is S.Zero
    assert P(And(Eq(Y[2], 1), Eq(Y[1], 1), Eq(Y[0], 0)), Eq(Y[1], 1)) == 0.1*Probability(Eq(Y[0], 0))

    # testing properties of Markov chain
    TO2 = Matrix([[S.One, 0, 0],[Rational(1, 3), Rational(1, 3), Rational(1, 3)],[0, Rational(1, 4), Rational(3, 4)]])
    TO3 = Matrix([[Rational(1, 4), Rational(3, 4), 0],[Rational(1, 3), Rational(1, 3), Rational(1, 3)],[0, Rational(1, 4), Rational(3, 4)]])
    Y2 = DiscreteMarkovChain('Y', trans_probs=TO2)
    Y3 = DiscreteMarkovChain('Y', trans_probs=TO3)
    assert Y3._transient2absorbing() == None
    raises (ValueError, lambda: Y3.fundamental_matrix())
    assert Y2.is_absorbing_chain() == True
    assert Y3.is_absorbing_chain() == False
    TO4 = Matrix([[Rational(1, 5), Rational(2, 5), Rational(2, 5)], [Rational(1, 10), S.Half, Rational(2, 5)], [Rational(3, 5), Rational(3, 10), Rational(1, 10)]])
    Y4 = DiscreteMarkovChain('Y', trans_probs=TO4)
    w = ImmutableMatrix([[Rational(11, 39), Rational(16, 39), Rational(4, 13)]])
    assert Y4.limiting_distribution == w
    assert Y4.is_regular() == True
    TS1 = MatrixSymbol('T', 3, 3)
    Y5 = DiscreteMarkovChain('Y', trans_probs=TS1)
    assert Y5.limiting_distribution(w, TO4).doit() == True
    TO6 = Matrix([[S.One, 0, 0, 0, 0],[S.Half, 0, S.Half, 0, 0],[0, S.Half, 0, S.Half, 0], [0, 0, S.Half, 0, S.Half], [0, 0, 0, 0, 1]])
    Y6 = DiscreteMarkovChain('Y', trans_probs=TO6)
    assert Y6._transient2absorbing() == ImmutableMatrix([[S.Half, 0], [0, 0], [0, S.Half]])
    assert Y6._transient2transient() == ImmutableMatrix([[0, S.Half, 0], [S.Half, 0, S.Half], [0, S.Half, 0]])
    assert Y6.fundamental_matrix() == ImmutableMatrix([[Rational(3, 2), S.One, S.Half], [S.One, S(2), S.One], [S.Half, S.One, Rational(3, 2)]])
    assert Y6.absorbing_probabilites() == ImmutableMatrix([[Rational(3, 4), Rational(1, 4)], [S.Half, S.Half], [Rational(1, 4), Rational(3, 4)]])

    # testing miscellaneous queries
    T = Matrix([[S.Half, Rational(1, 4), Rational(1, 4)],
                [Rational(1, 3), 0, Rational(2, 3)],
                [S.Half, S.Half, 0]])
    X = DiscreteMarkovChain('X', [0, 1, 2], T)
    assert P(Eq(X[1], 2) & Eq(X[2], 1) & Eq(X[3], 0),
    Eq(P(Eq(X[1], 0)), Rational(1, 4)) & Eq(P(Eq(X[1], 1)), Rational(1, 4))) == Rational(1, 12)
    assert P(Eq(X[2], 1) | Eq(X[2], 2), Eq(X[1], 1)) == Rational(2, 3)
    assert P(Eq(X[2], 1) & Eq(X[2], 2), Eq(X[1], 1)) is S.Zero
    assert P(Ne(X[2], 2), Eq(X[1], 1)) == Rational(1, 3)
    assert E(X[1]**2, Eq(X[0], 1)) == Rational(8, 3)
    assert variance(X[1], Eq(X[0], 1)) == Rational(8, 9)
    raises(ValueError, lambda: E(X[1], Eq(X[2], 1)))

def test_ContinuousMarkovChain():
    T1 = Matrix([[S(-2), S(2), S.Zero],
                 [S.Zero, S.NegativeOne, S.One],
                 [Rational(3, 2), Rational(3, 2), S(-3)]])
    C1 = ContinuousMarkovChain('C', [0, 1, 2], T1)
    assert C1.limiting_distribution() == ImmutableMatrix([[Rational(3, 19), Rational(12, 19), Rational(4, 19)]])

    T2 = Matrix([[-S.One, S.One, S.Zero], [S.One, -S.One, S.Zero], [S.Zero, S.One, -S.One]])
    C2 = ContinuousMarkovChain('C', [0, 1, 2], T2)
    A, t = C2.generator_matrix, symbols('t', positive=True)
    assert C2.transition_probabilities(A)(t) == Matrix([[S.Half + exp(-2*t)/2, S.Half - exp(-2*t)/2, 0],
                                                       [S.Half - exp(-2*t)/2, S.Half + exp(-2*t)/2, 0],
                                                       [S.Half - exp(-t) + exp(-2*t)/2, S.Half - exp(-2*t)/2, exp(-t)]])
    assert P(Eq(C2(1), 1), Eq(C2(0), 1), evaluate=False) == Probability(Eq(C2(1), 1))
    assert P(Eq(C2(1), 1), Eq(C2(0), 1)) == exp(-2)/2 + S.Half
    assert P(Eq(C2(1), 0) & Eq(C2(2), 1) & Eq(C2(3), 1),
                Eq(P(Eq(C2(1), 0)), S.Half)) == (Rational(1, 4) - exp(-2)/4)*(exp(-2)/2 + S.Half)
    assert P(Not(Eq(C2(1), 0) & Eq(C2(2), 1) & Eq(C2(3), 2)) |
                (Eq(C2(1), 0) & Eq(C2(2), 1) & Eq(C2(3), 2)),
                Eq(P(Eq(C2(1), 0)), Rational(1, 4)) & Eq(P(Eq(C2(1), 1)), Rational(1, 4))) is S.One
    assert E(C2(Rational(3, 2)), Eq(C2(0), 2)) == -exp(-3)/2 + 2*exp(Rational(-3, 2)) + S.Half
    assert variance(C2(Rational(3, 2)), Eq(C2(0), 1)) == ((S.Half - exp(-3)/2)**2*(exp(-3)/2 + S.Half)
                                                    + (Rational(-1, 2) - exp(-3)/2)**2*(S.Half - exp(-3)/2))
    raises(KeyError, lambda: P(Eq(C2(1), 0), Eq(P(Eq(C2(1), 1)), S.Half)))
    assert P(Eq(C2(1), 0), Eq(P(Eq(C2(5), 1)), S.Half)) == Probability(Eq(C2(1), 0))
    TS1 = MatrixSymbol('G', 3, 3)
    CS1 = ContinuousMarkovChain('C', [0, 1, 2], TS1)
    A = CS1.generator_matrix
    assert CS1.transition_probabilities(A)(t) == exp(t*A)
