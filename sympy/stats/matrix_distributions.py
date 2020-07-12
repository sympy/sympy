from sympy import S, Lambda, Basic, exp, multigamma
from sympy.core.sympify import sympify
from sympy.matrices import ImmutableMatrix, Inverse, Trace, Determinant, MatrixSymbol
from sympy.stats.random_matrix import RandomMatrixPSpace
from sympy.stats.rv import _symbol_converter, _value_check, RandomMatrixSymbol, NamedArgsMixin


def rv(symbol, dim, cls, args):
    args = list(map(sympify, args))
    dist = cls(*args)
    dist.check(*args)
    pspace = RandomMatrixPSpace(symbol, dist)
    return RandomMatrixSymbol(symbol, dim, dim, pspace)


class MatrixDistribution(Basic, NamedArgsMixin):
    def __new__(cls, *args):
        args = list(map(sympify, args))
        return Basic.__new__(cls, *args)

    @staticmethod
    def check(*args):
        pass

class MatrixGammaDistribution(MatrixDistribution):

    _argnames = ('alpha', 'beta', 'scale_matrix')

    @staticmethod
    def check(alpha, beta, scale_matrix):
        if not isinstance(scale_matrix , MatrixSymbol):
            _value_check(scale_matrix.is_positive_definite, "The shape \
                matrix must be positive definite.")
        _value_check(scale_matrix.shape[0] == scale_matrix.shape[1], "Should \
        be square matrix")
        _value_check(alpha.is_positive, "Shape parameter should be positive.")
        _value_check(beta.is_positive, "Scale parameter should be positive.")

    @property
    def set(self):
        k = self.args[3].shape[0]
        return S.Reals ** k

    def density(self, expr):
        alpha , beta , scale_matrix = self.alpha, self.beta, self.scale_matrix
        p = scale_matrix.shape[0]
        h_pspace = RandomMatrixPSpace('P', model=self)
        H = RandomMatrixSymbol('H', p, p, pspace=h_pspace)
        sigma_inv_x = - Inverse(scale_matrix)*H / beta
        term1 = exp(Trace(sigma_inv_x))/((beta**(p*alpha)) * multigamma(alpha, p))
        term2 = (Determinant(scale_matrix))**(-alpha)
        term3 = (Determinant(H))**(alpha - (p + 1)/2)
        return Lambda(H, term1 * term2 * term3)

def MatrixGamma(symbol, alpha, beta, scale_matrix):
    if isinstance(scale_matrix, list):
        scale_matrix = ImmutableMatrix(scale_matrix)
    dim = scale_matrix.shape[0]
    return rv(symbol, dim, MatrixGammaDistribution, (alpha, beta, scale_matrix))
