from sympy import (exp, S, ProductSet, sqrt, pi, symbols,
                Product, gamma, Dummy)
from sympy.matrices import Determinant, Matrix, Trace, MatrixSymbol
from sympy.stats import density
from sympy.stats.matrix_distributions import (MatrixGammaDistribution,
                MatrixGamma, MatrixPSpace)
from sympy.testing.pytest import raises


def test_MatrixPSpace():
    M = MatrixGammaDistribution(1, 2, [[2, 1], [1, 2]])
    MP = MatrixPSpace('M', M, 2, 2)
    assert MP.distribution == M
    raises(ValueError, lambda: MatrixPSpace('M', M, 1.2, 2))

def test_MatrixGamma():
    M = MatrixGamma('M', 1, 2, [[1, 0], [0, 1]])
    assert M.pspace.distribution.set == ProductSet(S.Reals, S.Reals)
    assert isinstance(density(M), MatrixGammaDistribution)
    X = MatrixSymbol('X', 2, 2)
    num = exp(Trace(Matrix([[-S(1)/2, 0], [0, -S(1)/2]])*X))
    assert density(M)(X).doit() == num/(4*pi*sqrt(Determinant(X)))
    assert density(M)([[2, 1], [1, 2]]).doit() == sqrt(3)*exp(-2)/(12*pi)
    X = MatrixSymbol('X', 1, 2)
    Y = MatrixSymbol('Y', 1, 2)
    assert density(M)([X, Y]).doit() == exp(-X[0, 0]/2 - Y[0, 1]/2)/(4*pi*sqrt(
                                X[0, 0]*Y[0, 1] - X[0, 1]*Y[0, 0]))
    # symbolic
    a, b = symbols('a b', positive=True)
    d = symbols('d', positive=True, integer=True)
    Y = MatrixSymbol('Y', d, d)
    Z = MatrixSymbol('Z', 2, 2)
    SM = MatrixSymbol('SM', d, d)
    M2 = MatrixGamma('M2', a, b, SM)
    M3 = MatrixGamma('M3', 2, 3, [[2, 1], [1, 2]])
    k = Dummy('k')
    exprd = pi**(-d*(d - 1)/4)*b**(-a*d)*exp(Trace((-1/b)*SM**(-1)*Y)
        )*Determinant(SM)**(-a)*Determinant(Y)**(a - d/2 - S(1)/2)/Product(
        gamma(-k/2 + a + S(1)/2), (k, 1, d))
    assert density(M2)(Y).dummy_eq(exprd)
    raises(NotImplementedError, lambda: density(M3 + M)(Z))
    raises(ValueError, lambda: density(M)(1))
    raises(ValueError, lambda: MatrixGamma('M', -1, 2, [[1, 0], [0, 1]]))
    raises(ValueError, lambda: MatrixGamma('M', -1, -2, [[1, 0], [0, 1]]))
    raises(ValueError, lambda: MatrixGamma('M', -1, 2, [[1, 0], [2, 1]]))
    raises(ValueError, lambda: MatrixGamma('M', -1, 2, [[1, 0], [0]]))
