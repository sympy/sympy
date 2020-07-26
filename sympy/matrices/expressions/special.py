from sympy.core.singleton import S
from sympy.core.sympify import _sympify
from sympy.matrices.common import NonInvertibleMatrixError
from .matexpr import MatrixExpr


class ZeroMatrix(MatrixExpr):
    """The Matrix Zero 0 - additive identity

    Examples
    ========

    >>> from sympy import MatrixSymbol, ZeroMatrix
    >>> A = MatrixSymbol('A', 3, 5)
    >>> Z = ZeroMatrix(3, 5)
    >>> A + Z
    A
    >>> Z*A.T
    0
    """
    is_ZeroMatrix = True

    def __new__(cls, m, n):
        m, n = _sympify(m), _sympify(n)
        cls._check_dim(m)
        cls._check_dim(n)

        return super().__new__(cls, m, n)

    @property
    def shape(self):
        return (self.args[0], self.args[1])

    def _eval_power(self, exp):
        # exp = -1, 0, 1 are already handled at this stage
        if (exp < 0) == True:
            raise NonInvertibleMatrixError("Matrix det == 0; not invertible")
        return self

    def _eval_transpose(self):
        return ZeroMatrix(self.cols, self.rows)

    def _eval_trace(self):
        return S.Zero

    def _eval_determinant(self):
        return S.Zero

    def _eval_inverse(self):
        raise NonInvertibleMatrixError("Matrix det == 0; not invertible.")

    def conjugate(self):
        return self

    def _entry(self, i, j, **kwargs):
        return S.Zero


class GenericZeroMatrix(ZeroMatrix):
    """
    A zero matrix without a specified shape

    This exists primarily so MatAdd() with no arguments can return something
    meaningful.
    """
    def __new__(cls):
        # super(ZeroMatrix, cls) instead of super(GenericZeroMatrix, cls)
        # because ZeroMatrix.__new__ doesn't have the same signature
        return super(ZeroMatrix, cls).__new__(cls)

    @property
    def rows(self):
        raise TypeError("GenericZeroMatrix does not have a specified shape")

    @property
    def cols(self):
        raise TypeError("GenericZeroMatrix does not have a specified shape")

    @property
    def shape(self):
        raise TypeError("GenericZeroMatrix does not have a specified shape")

    # Avoid Matrix.__eq__ which might call .shape
    def __eq__(self, other):
        return isinstance(other, GenericZeroMatrix)

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return super().__hash__()
