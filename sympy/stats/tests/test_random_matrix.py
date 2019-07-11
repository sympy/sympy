from sympy import sqrt, exp, Trace, pi, S
from sympy.stats import GUE, density
from sympy.stats.rv import RandomMatrixSymbol, Density
from sympy.stats.random_matrix_models import GaussianEnsemble
from sympy.utilities.pytest import raises

def test_GaussianEnsemble():
    G = GaussianEnsemble('G', 3)
    assert density(G) == Density(G)

def test_GaussianUnitaryEnsemble():
    H = RandomMatrixSymbol('H', 3, 3)
    G = GUE('U', 3)
    raises(ValueError, lambda: GUE('U', 3.5))
    assert density(G)(H) == sqrt(2)*exp(3*Trace(H**2)/2)/(4*pi**(S(9)/2))
