from sympy import S, Basic, exp, multigamma
from sympy.core.sympify import sympify, _sympify
from sympy.matrices import (ImmutableMatrix, Inverse, Trace, Determinant,
                            MatrixSymbol, MatrixBase)
from sympy.stats.rv import (_value_check, RandomMatrixSymbol, NamedArgsMixin, PSpace,
                            _symbol_converter)


class MatrixPSpace(PSpace):
    """
    Represents probability space for
    Matrix Distributions
    """
    def __new__(cls, sym, distribution, dim_n, dim_m):
        sym = _symbol_converter(sym)
        dim_n, dim_m = _sympify(dim_n), _sympify(dim_m)
        if not (dim_n.is_integer and dim_m.is_integer):
            raise ValueError("Dimensions should be integers")
        return Basic.__new__(cls, sym, distribution, dim_n, dim_m)

    distribution = property(lambda self: self.args[1])
    symbol = property(lambda self: self.args[0])

    @property
    def value(self):
        return RandomMatrixSymbol(self.symbol, self.args[2], self.args[3], self)

    def compute_density(self, expr, *args):
        rms = expr.atoms(RandomMatrixSymbol)
        if len(rms) > 2 or (not isinstance(expr, RandomMatrixSymbol)):
            raise NotImplementedError("Currently, no algorithm has been "
                    "implemented to handle general expressions containing "
                    "multiple matrix distributions.")
        return self.distribution.pdf(expr)


def rv(symbol, cls, args):
    args = list(map(sympify, args))
    dist = cls(*args)
    dist.check(*args)
    dim = dist.dimension
    pspace = MatrixPSpace(symbol, dist, dim[0], dim[1])
    return pspace.value


class MatrixDistribution(Basic, NamedArgsMixin):
    def __new__(cls, *args):
        args = list(map(sympify, args))
        return Basic.__new__(cls, *args)

    @staticmethod
    def check(*args):
        pass

    def __call__(self, expr):
        if isinstance(expr, list):
            expr = ImmutableMatrix(expr)
        return self.pdf(expr)

class MatrixGammaDistribution(MatrixDistribution):

    _argnames = ('alpha', 'beta', 'scale_matrix')

    @staticmethod
    def check(alpha, beta, scale_matrix):
        if not isinstance(scale_matrix , MatrixSymbol):
            _value_check(scale_matrix.is_positive_definite, "The shape "
                "matrix must be positive definite.")
        _value_check(scale_matrix.shape[0] == scale_matrix.shape[1], "Should "
        "be square matrix")
        _value_check(alpha.is_positive, "Shape parameter should be positive.")
        _value_check(beta.is_positive, "Scale parameter should be positive.")

    @property
    def set(self):
        k = self.scale_matrix.shape[0]
        return S.Reals ** k

    @property
    def dimension(self):
        return self.scale_matrix.shape

    def pdf(self, x):
        alpha , beta , scale_matrix = self.alpha, self.beta, self.scale_matrix
        p = scale_matrix.shape[0]
        if isinstance(x, list):
            x = ImmutableMatrix(x)
        if not isinstance(x, (MatrixBase, MatrixSymbol)):
            raise ValueError("%s should be an isinstance of Matrix "
                    "or MatrixSymbol" % str(x))
        sigma_inv_x = - Inverse(scale_matrix)*x / beta
        term1 = exp(Trace(sigma_inv_x))/((beta**(p*alpha)) * multigamma(alpha, p))
        term2 = (Determinant(scale_matrix))**(-alpha)
        term3 = (Determinant(x))**(alpha - S(p + 1)/2)
        return term1 * term2 * term3


def MatrixGamma(symbol, alpha, beta, scale_matrix):
    if isinstance(scale_matrix, list):
        scale_matrix = ImmutableMatrix(scale_matrix)
    return rv(symbol, MatrixGammaDistribution, (alpha, beta, scale_matrix))
