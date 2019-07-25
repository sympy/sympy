from sympy import sqrt, exp, Trace, pi, S, Integral, MatrixSymbol
from sympy.stats import (GaussianUnitaryEnsemble as GUE, density,
                         GaussianOrthogonalEnsemble as GOE,
                         GaussianSymplecticEnsemble as GSE)
from sympy.stats.rv import RandomMatrixSymbol, Density
from sympy.stats.random_matrix_models import GaussianEnsemble
from sympy.utilities.pytest import raises

def test_GaussianEnsemble():
    G = GaussianEnsemble('G', 3)
    assert density(G) == Density(G)
    raises(ValueError, lambda: GaussianEnsemble('G', 3.5))

def test_GaussianUnitaryEnsemble():
    H = RandomMatrixSymbol('H', 3, 3)
    G = GUE('U', 3)
    assert density(G)(H) == sqrt(2)*exp(-3*Trace(H**2)/2)/(4*pi**(S(9)/2))

def test_GaussianOrthogonalEnsemble():
    H = RandomMatrixSymbol('H', 3, 3)
    _H = MatrixSymbol('_H', 3, 3)
    G = GOE('O', 3)
    assert density(G)(H) == exp(-3*Trace(H**2)/4)/Integral(exp(-3*Trace(_H**2)/4), _H)

def test_GaussianSymplecticEnsemble():
    H = RandomMatrixSymbol('H', 3, 3)
    _H = MatrixSymbol('_H', 3, 3)
    G = GSE('O', 3)
    assert density(G)(H) == exp(-3*Trace(H**2))/Integral(exp(-3*Trace(_H**2)), _H)
