from __future__ import print_function, division

from sympy import Basic, exp, pi, Lambda, Trace, S
from sympy.core.sympify import _sympify
from sympy.stats.rv import _symbol_converter, Density, RandomMatrixSymbol
from sympy.stats.random_matrix import RandomMatrixPSpace

__all__ = [
    'GaussianUnitaryEnsemble',
]

class RandomMatrixEnsemble(Basic):
    """
    Abstract class for random matrix ensembles.
    It acts as an umbrella for all the ensembles
    defined in sympy.stats.
    """
    pass

class GaussianEnsemble(RandomMatrixEnsemble):
    """
    Abstract class for Gaussian ensembles.
    Contains the properties common to all the
    gaussian ensembles.

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Random_matrix#Gaussian_ensembles
    .. [2] https://arxiv.org/pdf/1712.07903.pdf
    """
    def __new__(cls, sym, dim=None):
        sym, dim = _symbol_converter(sym), _sympify(dim)
        if dim.is_integer == False:
            raise ValueError("Dimension of the random matrices must be "
                                "integers, received %s instead."%(dim))
        self = Basic.__new__(cls, sym, dim)
        rmp = RandomMatrixPSpace(sym, model=self)
        return RandomMatrixSymbol(sym, dim, dim, pspace=rmp)

    symbol = property(lambda self: self.args[0])
    dimension = property(lambda self: self.args[1])

    def density(self, expr):
        return Density(expr)

class GaussianUnitaryEnsemble(GaussianEnsemble):
    """
    Represents Gaussian Unitary Ensembles.

    Examples
    ========

    >>> from sympy.stats import GaussianUnitaryEnsemble as GUE, density
    >>> G = GUE('U', 2)
    >>> density(G)
    Lambda(H, exp(Trace(H**2))/(2*pi**2))
    """
    @property
    def normalization_constant(self):
        n = self.dimension
        return 2**(S(n)/2) * pi**(S(n**2)/2)

    def density(self, expr):
        n, ZGUE = self.dimension, self.normalization_constant
        h_pspace = RandomMatrixPSpace('P', model=self)
        H = RandomMatrixSymbol('H', n, n, pspace=h_pspace)
        return Lambda(H, exp(S(n)/2 * Trace(H**2))/ZGUE)
