from __future__ import print_function, division

from sympy import Basic, exp, pi, Lambda, Trace, S, MatrixSymbol, Integral
from sympy.core.sympify import _sympify
from sympy.stats.rv import _symbol_converter, Density, RandomMatrixSymbol
from sympy.stats.random_matrix import RandomMatrixPSpace

__all__ = [
    'GaussianUnitaryEnsemble',
    'GaussianOrthogonalEnsemble',
    'GaussianSymplecticEnsemble'
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
    Lambda(H, exp(-Trace(H**2))/(2*pi**2))
    """
    @property
    def normalization_constant(self):
        n = self.dimension
        return 2**(S(n)/2) * pi**(S(n**2)/2)

    def density(self, expr):
        n, ZGUE = self.dimension, self.normalization_constant
        h_pspace = RandomMatrixPSpace('P', model=self)
        H = RandomMatrixSymbol('H', n, n, pspace=h_pspace)
        return Lambda(H, exp(-S(n)/2 * Trace(H**2))/ZGUE)

class GaussianOrthogonalEnsemble(GaussianEnsemble):
    """
    Represents Gaussian Orthogonal Ensembles.

    Examples
    ========

    >>> from sympy.stats import GaussianOrthogonalEnsemble as GOE, density
    >>> G = GOE('U', 2)
    >>> density(G)
    Lambda(H, exp(-Trace(H**2)/2)/Integral(exp(-Trace(_H**2)/2), _H))
    """
    @property
    def normalization_constant(self):
        n = self.dimension
        _H = MatrixSymbol('_H', n, n)
        return Integral(exp(-S(n)/4 * Trace(_H**2)))

    def density(self, expr):
        n, ZGOE = self.dimension, self.normalization_constant
        h_pspace = RandomMatrixPSpace('P', model=self)
        H = RandomMatrixSymbol('H', n, n, pspace=h_pspace)
        return Lambda(H, exp(-S(n)/4 * Trace(H**2))/ZGOE)

class GaussianSymplecticEnsemble(GaussianEnsemble):
    """
    Represents Gaussian Symplectic Ensembles.

    Examples
    ========

    >>> from sympy.stats import GaussianSymplecticEnsemble as GSE, density
    >>> G = GSE('U', 2)
    >>> density(G)
    Lambda(H, exp(-2*Trace(H**2))/Integral(exp(-2*Trace(_H**2)), _H))
    """
    @property
    def normalization_constant(self):
        n = self.dimension
        _H = MatrixSymbol('_H', n, n)
        return Integral(exp(-S(n) * Trace(_H**2)))

    def density(self, expr):
        n, ZGSE = self.dimension, self.normalization_constant
        h_pspace = RandomMatrixPSpace('P', model=self)
        H = RandomMatrixSymbol('H', n, n, pspace=h_pspace)
        return Lambda(H, exp(-S(n) * Trace(H**2))/ZGSE)
