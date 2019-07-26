from sympy.core.basic import Basic
from sympy.core.sympify import _sympify, converter
from sympy.core.numbers import Integer
from sympy.external import import_module

from .matexpr import MatrixExpr, MatrixElement

numpy = import_module('numpy')

def sympify_numpy_array(a):
    from sympy.tensor.array import Array

    if len(a.shape) == 2:
        return WrappedNumPyMatrix(a)
    else:
        return Array(a.flat, a.shape)

if numpy:
    converter[numpy.ndarray] = sympify_numpy_array

class WrappedNumPyMatrix(Basic):
    """"""
    is_Atom = True

    def _hashable_content(self):
        return tuple(self._mat.tostring())

    def __new__(cls, arg, copy=True):
        if not isinstance(arg, numpy.ndarray):
            raise ValueError('Only Numpy array type is supported.')

        if not len(arg.shape) == 2:
            raise ValueError('Only 2-Dimensional array is supported')

        if copy:
            obj = Basic.__new__(cls, arg.copy())
        else:
            obj = Basic.__new__(cls, arg)
        return obj

    @property
    def _mat(self):
        return self.args[0]

class NumPyMatrix(MatrixExpr):
    is_Atom = True

    def __new__(cls, arg, copy=True):
        if not isinstance(arg, numpy.ndarray):
            raise ValueError('Only Numpy array type is supported.')

        if not len(arg.shape) == 2:
            raise ValueError('Only 2-Dimensional array is supported')

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
        return Integer(self._wrapped._mat.shape[0])

    @property
    def cols(self):
        return Integer(self._wrapped._mat.shape[1])

    @property
    def shape(self):
        return self.rows, self.cols

    def _eval_pow(self, exp):
        mat = self._wrapped._mat
        return NumPyMatrix(
            numpy.linalg.matrix_power(mat, int(exp)), copy=False)

    def _eval_inverse(self):
        mat = self._wrapped._mat
        return NumPyMatrix(numpy.linalg.inv(mat), copy=False)

    def _eval_transpose(self):
        mat = self._wrapped._mat
        return NumPyMatrix(numpy.transpose(mat), copy=False)

    def _eval_conjugate(self):
        mat = self._wrapped._mat
        return NumPyMatrix(numpy.conj(mat), copy=False)

    def _eval_adjoint(self):
        mat = self._wrapped._mat
        return NumPyMatrix(numpy.conj(numpy.transpose(mat)), copy=False)

    def _eval_determinant(self):
        mat = self._wrapped._mat
        return _sympify(numpy.linalg.det(mat))

    def _eval_rref(self, tol=10**-8):
        # Experimental custom implementation
        fabs = numpy.absolute
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

        rows, cols = mat.shape

        row = 0
        for col in range(cols):
            if is_zero(mat[row, col]):
                pivot_pos = find_pivot(mat[row:, col])
                if pivot_pos == None:
                    continue
                mat[row, col:], mat[row + pivot_pos, col:] = \
                    mat[row + pivot_pos, col:], mat[row, col:]

            mat[row, col:] /= mat[row, col]

            for i in range(rows):
                if i == row:
                    continue
                mat[i, col:] -=  mat[i, col] * mat[row, col:]

            row += 1

        return NumPyMatrix(mat, copy=False)
