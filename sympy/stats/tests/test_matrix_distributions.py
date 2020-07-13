from sympy import exp, MatrixSymbol, S, ProductSet, sqrt, pi
from sympy.matrices import Determinant, Matrix, Trace
from sympy.stats import density
from sympy.stats.matrix_distributions import MatrixGammaDistribution, MatrixGamma
from sympy.testing.pytest import raises

def test_MatrixGamma():
    M = MatrixGamma('M', 1, 2, [[1, 0], [0, 1]])
    assert M.pspace.model.set == ProductSet(S.Reals, S.Reals)
    assert isinstance(density(M), MatrixGammaDistribution)
    X = MatrixSymbol('X', 2, 2)
    num = exp(Trace(Matrix([[-S(1)/2, 0], [0, -S(1)/2]])*X))
    assert density(M)(X).doit() == num/(4*pi*sqrt(Determinant(X)))
    assert density(M)([[2, 1], [1, 2]]).doit() == sqrt(3)*exp(-2)/(12*pi)
    raises(ValueError, lambda: MatrixGamma('M', -1, 2, [[1, 0], [0, 1]]))
    raises(ValueError, lambda: MatrixGamma('M', -1, -2, [[1, 0], [0, 1]]))
    raises(ValueError, lambda: MatrixGamma('M', -1, 2, [[1, 0], [2, 1]]))
    raises(ValueError, lambda: MatrixGamma('M', -1, 2, [[1, 0], [0]]))
