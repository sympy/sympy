import mpmath

from sympy.core.basic import Basic
from sympy.core.sympify import _sympify
from sympy.core.numbers import Integer

from .matexpr import MatrixExpr, MatrixElement

class WrappedMpmathMatrix(Basic):
    """"""
    is_Atom = True

    def _hashable_content(self):
        return tuple(self._mat._toliststr())

    def __new__(cls, arg, copy=True):
        if not isinstance(arg, mpmath.matrix):
            raise ValueError('Only mpmath.matrix type is supported.')

        if copy:
            obj = Basic.__new__(cls, arg.copy())
        else:
            obj = Basic.__new__(cls, arg)
        return obj

    @property
    def _mat(self):
        return self.args[0]

class MpmathMatrix(MatrixExpr):
    is_Atom = True

    def __new__(cls, arg, copy=True):
        if not isinstance(arg, mpmath.matrix):
            raise ValueError('Only mpmath.matrix type is supported.')

        obj = MatrixExpr.__new__(cls, arg, copy=copy)
        return obj

    @property
    def _wrapped(self):
        return self.args[0]

    def doit(self, **kwargs):
        return self

    def _entry(self, i, j, **kwargs):
        return MatrixElement(self, i, j)

    def as_explicit(self):
        from sympy.matrices.immutable import ImmutableDenseMatrix
        return ImmutableDenseMatrix(self._wrapped._mat.tolist())

    @property
    def rows(self):
        return Integer(self._wrapped._mat.rows)

    @property
    def cols(self):
        return Integer(self._wrapped._mat.cols)

    @property
    def shape(self):
        return self.rows, self.cols

    def _eval_pow(self, exp):
        # XXX Mpmath only supports integer power
        return MpmathMatrix(self._wrapped._mat ** int(exp), copy=False)

    def _eval_inverse(self):
        return MpmathMatrix(self._wrapped._mat ** -1, copy=False)

    def _eval_transpose(self):
        return MpmathMatrix(self._wrapped._mat.transpose(), copy=False)

    def _eval_conjugate(self):
        return MpmathMatrix(self._wrapped._mat.conjugate(), copy=False)

    def _eval_adjoint(self):
        return MpmathMatrix(self._wrapped._mat.transpose_conj(), copy=False)

    def _eval_determinant(self):
        return _sympify(mpmath.det(self._wrapped._mat))

    def _eval_rref(self, tol=10**-10):
        # Experimental custom implementation
        fabs = mpmath.fabs
        mat = self._wrapped._mat.copy()

        def is_zero(val):
            return fabs(val) < tol

        def find_pivot(vec):
            max_index = None
            max_value = 0
            for index, value in enumerate(vec):
                if fabs(value) > fabs(max_value) and not is_zero(value):
                    value = max_value
                    max_index = index

            return max_index

        row = 0
        for col in range(mat.cols):
            if is_zero(mat[row, col]):
                pivot_pos = find_pivot(mat[row:, col])
                if pivot_pos == None:
                    continue
                mat[row, col:], mat[row + pivot_pos, col:] = \
                    mat[row + pivot_pos, col:], mat[row, col:]

            mat[row, col:] /= mat[row, col]

            for i in range(mat.rows):
                if i == row:
                    continue
                mat[i, col:] -=  mat[i, col] * mat[row, col:]

            row += 1

        return MpmathMatrix(mat, copy=False)
