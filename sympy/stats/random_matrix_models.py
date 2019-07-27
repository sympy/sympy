from __future__ import print_function, division

from sympy import (Basic, exp, pi, Lambda, Trace, S, MatrixSymbol, Integral,
                   gamma, Product, Dummy, Sum, Abs, IndexedBase)
from sympy.core.sympify import _sympify
from sympy.multipledispatch import dispatch
from sympy.stats.rv import _symbol_converter, Density, RandomMatrixSymbol
from sympy.stats.random_matrix import RandomMatrixPSpace
from sympy.tensor.array import ArrayComprehension

__all__ = [
    'GaussianUnitaryEnsemble',
    'GaussianOrthogonalEnsemble',
    'GaussianSymplecticEnsemble',
    'joint_eigen_distribution',
    'level_spacing_distribution'
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

    def _compute_normalization_constant(self, beta, n):
        """
        Helper function for computing normalization
        constant for joint probability density of eigen
        values of Gaussian ensembles.

        References
        ==========

        .. [1] https://en.wikipedia.org/wiki/Selberg_integral#Mehta's_integral
        """
        n = S(n)
        prod_term = lambda j: gamma(1 + beta*S(j)/2)/gamma(S(1) + beta/S(2))
        j = Dummy('j', integer=True, positive=True)
        term1 = Product(prod_term(j), (j, 1, n)).doit()
        term2 = (2/(beta*n))**(beta*n*(n - 1)/4 + n/2)
        term3 = (2*pi)**(n/2)
        return term1 * term2 * term3

    def _compute_joint_eigen_dsitribution(self, beta):
        n = self.dimension
        Zbn = self._compute_normalization_constant(beta, n)
        l = IndexedBase('l')
        i = Dummy('i', integer=True, positive=True)
        j = Dummy('j', integer=True, positive=True)
        k = Dummy('k', integer=True, positive=True)
        term1 = exp((-S(n)/2) * Sum(l[k]**2, (k, 1, n)).doit())
        sub_term = Lambda(i, Product(Abs(l[j] - l[i])**beta, (j, i + 1, n)))
        term2 = Product(sub_term(i).doit(), (i, 1, n - 1)).doit()
        syms = ArrayComprehension(l[k], (k, 1, n)).doit()
        return Lambda(syms, (term1 * term2)/Zbn)

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

    def joint_eigen_distribution(self):
        return self._compute_joint_eigen_dsitribution(2)

    def level_spacing_distribution(self):
        s = Dummy('s')
        f = (32/pi**2)*(s**2)*exp((-4/pi)*s**2)
        return Lambda(s, f)

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

    def joint_eigen_distribution(self):
        return self._compute_joint_eigen_dsitribution(1)

    def level_spacing_distribution(self):
        s = Dummy('s')
        f = (pi/2)*s*exp((-pi/4)*s**2)
        return Lambda(s, f)

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

    def joint_eigen_distribution(self):
        return self._compute_joint_eigen_dsitribution(4)

    def level_spacing_distribution(self):
        s = Dummy('s')
        f = ((S(2)**18)/((S(3)**6)*(pi**3)))*(s**4)*exp((-64/(9*pi))*s**2)
        return Lambda(s, f)

@dispatch(RandomMatrixSymbol)
def joint_eigen_distribution(mat):
    return mat.pspace.model.joint_eigen_distribution()

def level_spacing_distribution(mat):
    return mat.pspace.model.level_spacing_distribution()
