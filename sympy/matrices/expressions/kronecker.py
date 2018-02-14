"""Implementation of the Kronecker product"""

from __future__ import print_function, division


from sympy.core import Add, Mul, Pow, sympify, prod
from sympy.core.compatibility import range
from sympy.strategies import unpack, flatten, condition, exhaust, do_one

from sympy.functions import adjoint
from sympy.matrices.expressions.transpose import transpose

from sympy.matrices.expressions.matexpr import MatrixExpr, ShapeError
from sympy.matrices.matrices import MatrixBase


def kronecker_product(*matrices):
    """
    The Kronecker product of two or more arguments.

    Returns a symbolic ``KroneckerProduct`` object.

    Examples
    ========

    >>> from sympy import pprint
    >>> from sympy.matrices import kronecker_product, MatrixSymbol
    >>> A = MatrixSymbol('A', 2, 2)
    >>> B = MatrixSymbol('B', 2, 2)
    >>> kronecker_product(A)
    A
    >>> kronecker_product(A, B)
    AxB
    >>> kronecker_product(A, B)[0, 1]
    A[0, 0]*B[0, 1]
    >>> pprint(kronecker_product(A, B).as_explicit(), use_unicode=True)
    ⎡A₀₀⋅B₀₀  A₀₀⋅B₀₁  A₀₁⋅B₀₀  A₀₁⋅B₀₁⎤
    ⎢                                  ⎥
    ⎢A₀₀⋅B₁₀  A₀₀⋅B₁₁  A₀₁⋅B₁₀  A₀₁⋅B₁₁⎥
    ⎢                                  ⎥
    ⎢A₁₀⋅B₀₀  A₁₀⋅B₀₁  A₁₁⋅B₀₀  A₁₁⋅B₀₁⎥
    ⎢                                  ⎥
    ⎣A₁₀⋅B₁₀  A₁₀⋅B₁₁  A₁₁⋅B₁₀  A₁₁⋅B₁₁⎦
    """
    if not matrices:
        raise TypeError("Empty Kronecker product is undefined")
    validate(*matrices)
    if len(matrices) == 1:
        return matrices[0]
    else:
        return KroneckerProduct(*matrices).doit()


class KroneckerProduct(MatrixExpr):
    """
    The Kronecker product is a non-commutative product of matrices.
    Given two matrices of dimension (m, n) and (s, t) it produces a matrix
    of of dimension (m s, n t).

    This is a symbolic object that simply stores its argument without
    evaluating it. To actually compute the product, use the function
    ``kronecker_product()`` or call the the ``.doit()`` or  ``.as_explicit()``
    methods.

    >>> from sympy.matrices import kronecker_product, KroneckerProduct, MatrixSymbol
    >>> A = MatrixSymbol('A', 5, 5)
    >>> B = MatrixSymbol('B', 5, 5)
    >>> isinstance(kronecker_product(A, B), KroneckerProduct)
    True
    """
    is_KroneckerProduct = True

    def __new__(cls, *args, **kwargs):
        args = list(map(sympify, args))
        check = kwargs.get('check', True)
        if check:
            validate(*args)
        return super(KroneckerProduct, cls).__new__(cls, *args)

    @property
    def shape(self):
        rows, cols = self.args[0].shape
        for mat in self.args[1:]:
            rows *= mat.rows
            cols *= mat.cols
        return (rows, cols)

    def _entry(self, i, j):
        result = 1
        for mat in reversed(self.args):
            i, m = divmod(i, mat.rows)
            j, n = divmod(j, mat.cols)
            result *= mat[m, n]
        return result

    def _eval_adjoint(self):
        return KroneckerProduct(*list(map(adjoint, self.args))).doit()

    def _eval_conjugate(self):
        return KroneckerProduct(*[a.conjugate() for a in self.args]).doit()

    def _eval_transpose(self):
        return KroneckerProduct(*list(map(transpose, self.args))).doit()

    def _eval_trace(self):
        from .trace import trace
        return prod(trace(a) for a in self.args)

    def _eval_determinant(self):
        from .determinant import det, Determinant
        if not all(a.is_square for a in self.args):
            return Determinant(self)

        m = self.rows
        return prod(det(a)**(m/a.rows) for a in self.args)

    def _eval_inverse(self):
        try:
            return KroneckerProduct(*[a.inverse() for a in self.args])
        except ShapeError:
            from sympy.matrices.expressions.inverse import Inverse
            return Inverse(self)

    def doit(self, **kwargs):
        deep = kwargs.get('deep', True)
        if deep:
            args = [arg.doit(**kwargs) for arg in self.args]
        else:
            args = self.args
        return canonicalize(KroneckerProduct(*args))

def validate(*args):
    if not all(arg.is_Matrix for arg in args):
        raise TypeError("Mix of Matrix and Scalar symbols")


def extract_commutative(kron):
    c_part = []
    nc_part = []
    for arg in kron.args:
        c, nc = arg.args_cnc()
        c_part.extend(c)
        nc_part.append(Mul._from_args(nc))

    c_part = Mul(*c_part)
    if c_part != 1:
        return c_part*KroneckerProduct(*nc_part)
    return kron


def matrix_kronecker_product(*matrices):

    # Make sure we have a sequence of Matrices
    if not all(isinstance(m, MatrixBase) for m in matrices):
        raise TypeError(
            'Sequence of Matrices expected, got: %s' % repr(matrices)
        )

    # Pull out the first element in the product.
    matrix_expansion = matrices[-1]
    # Do the kronecker product working from right to left.
    for mat in reversed(matrices[:-1]):
        rows = mat.rows
        cols = mat.cols
        # Go through each row appending kronecker product to.
        # running matrix_expansion.
        for i in range(rows):
            start = matrix_expansion*mat[i*cols]
            # Go through each column joining each item
            for j in range(cols - 1):
                start = start.row_join(
                    matrix_expansion*mat[i*cols + j + 1]
                )
            # If this is the first element, make it the start of the
            # new row.
            if i == 0:
                next = start
            else:
                next = next.col_join(start)
        matrix_expansion = next

    MatrixClass = max(matrices, key=lambda M: M._class_priority).__class__
    if isinstance(matrix_expansion, MatrixClass):
        return matrix_expansion
    else:
        return MatrixClass(matrix_expansion)


def explicit_kronecker_product(kron):
    # Make sure we have a sequence of Matrices
    if not all(isinstance(m, MatrixBase) for m in kron.args):
        return kron

    return matrix_kronecker_product(*kron.args)


rules = (unpack,
         explicit_kronecker_product,
         flatten,
         extract_commutative)

canonicalize = exhaust(condition(lambda x: isinstance(x, KroneckerProduct),
                                 do_one(*rules)))
