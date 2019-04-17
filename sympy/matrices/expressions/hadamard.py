from __future__ import print_function, division

from sympy.core import Mul, sympify
from sympy.matrices.expressions.matexpr import MatrixExpr, ShapeError, Identity
from sympy.strategies import unpack, flatten, condition, exhaust, do_one


def hadamard_product(*matrices):
    """
    Return the elementwise (aka Hadamard) product of matrices.

    Examples
    ========

    >>> from sympy.matrices import hadamard_product, MatrixSymbol
    >>> A = MatrixSymbol('A', 2, 3)
    >>> B = MatrixSymbol('B', 2, 3)
    >>> hadamard_product(A)
    A
    >>> hadamard_product(A, B)
    A.*B
    >>> hadamard_product(A, B)[0, 1]
    A[0, 1]*B[0, 1]
    """
    if not matrices:
        raise TypeError("Empty Hadamard product is undefined")
    validate(*matrices)
    if len(matrices) == 1:
        return matrices[0]
    else:
        matrices = [i for i in matrices if not i.is_Identity]
        return HadamardProduct(*matrices).doit()


class HadamardProduct(MatrixExpr):
    """
    Elementwise product of matrix expressions

    This is a symbolic object that simply stores its argument without
    evaluating it. To actually compute the product, use the function
    ``hadamard_product()``.

    >>> from sympy.matrices import hadamard_product, HadamardProduct, MatrixSymbol
    >>> A = MatrixSymbol('A', 5, 5)
    >>> B = MatrixSymbol('B', 5, 5)
    >>> isinstance(hadamard_product(A, B), HadamardProduct)
    True
    """
    is_HadamardProduct = True

    def __new__(cls, *args, **kwargs):
        args = list(map(sympify, args))
        check = kwargs.get('check', True)
        if check:
            validate(*args)
        return super(HadamardProduct, cls).__new__(cls, *args)

    @property
    def shape(self):
        return self.args[0].shape

    def _entry(self, i, j, **kwargs):
        return Mul(*[arg._entry(i, j, **kwargs) for arg in self.args])

    def _eval_transpose(self):
        from sympy.matrices.expressions.transpose import transpose
        return HadamardProduct(*list(map(transpose, self.args)))

    def doit(self, **ignored):
        return canonicalize(self)

    def _eval_derivative_matrix_lines(self, x):
        from sympy.core.expr import ExprBuilder
        from sympy.codegen.array_utils import CodegenArrayDiagonal, CodegenArrayTensorProduct
        from sympy.matrices.expressions.matexpr import _make_matrix

        with_x_ind = [i for i, arg in enumerate(self.args) if arg.has(x)]
        lines = []
        for ind in with_x_ind:
            left_args = self.args[:ind]
            right_args = self.args[ind+1:]

            d = self.args[ind]._eval_derivative_matrix_lines(x)
            hadam = hadamard_product(*(right_args + left_args))
            diagonal = [(0, 2), (3, 4)]
            diagonal = [e for i, e in enumerate(diagonal) if self.shape[i] != 1]
            for i in d:
                ptr1 = i.first_pointer
                ptr2 = i.second_pointer
                subexpr = ExprBuilder(
                    CodegenArrayDiagonal,
                    [
                        ExprBuilder(
                            CodegenArrayTensorProduct,
                            [
                                ExprBuilder(_make_matrix, [i._lines[0]]),
                                hadam,
                                ExprBuilder(_make_matrix, [i._lines[1]]),
                            ]
                        ),
                        *diagonal
                    ],

                )
                i._first_pointer_parent = subexpr.args[0].args[0].args
                i._first_pointer_index = 0
                i._second_pointer_parent = subexpr.args[0].args[2].args
                i._second_pointer_index = 0
                i._lines = [subexpr]
                lines.append(i)

        return lines


def validate(*args):
    if not all(arg.is_Matrix for arg in args):
        raise TypeError("Mix of Matrix and Scalar symbols")
    A = args[0]
    for B in args[1:]:
        if A.shape != B.shape:
            raise ShapeError("Matrices %s and %s are not aligned" % (A, B))

rules = (unpack,
         flatten)

canonicalize = exhaust(condition(lambda x: isinstance(x, HadamardProduct),
                                 do_one(*rules)))


def hadamard_power(base, exp):
    base = sympify(base)
    exp = sympify(exp)
    if exp == 1:
        return base
    if not base.is_Matrix:
        return base**exp
    if exp.is_Matrix:
        raise ValueError("cannot raise expression to a matrix")
    return HadamardPower(base, exp)


class HadamardPower(MatrixExpr):
    """
    Elementwise power of matrix expressions
    """

    def __new__(cls, base, exp):
        base = sympify(base)
        exp = sympify(exp)
        obj = super(HadamardPower, cls).__new__(cls, base, exp)
        return obj

    @property
    def base(self):
        return self._args[0]

    @property
    def exp(self):
        return self._args[1]

    @property
    def shape(self):
        return self.base.shape

    def _entry(self, i, j, **kwargs):
        return self.base._entry(i, j, **kwargs)**self.exp

    def _eval_transpose(self):
        from sympy.matrices.expressions.transpose import transpose
        return HadamardPower(transpose(self.base), self.exp)

    def _eval_derivative_matrix_lines(self, x):
        from sympy.codegen.array_utils import CodegenArrayTensorProduct
        from sympy.codegen.array_utils import CodegenArrayContraction, CodegenArrayDiagonal
        from sympy.core.expr import ExprBuilder
        from sympy.matrices.expressions.matexpr import _make_matrix

        lr = self.base._eval_derivative_matrix_lines(x)
        for i in lr:
            ptr1 = i.first_pointer
            ptr2 = i.second_pointer
            diagonal = [(1, 2), (3, 4)]
            diagonal = [e for i, e in enumerate(diagonal) if self.base.shape[i] != 1]
            subexpr = ExprBuilder(
                CodegenArrayDiagonal,
                [
                    ExprBuilder(
                        CodegenArrayTensorProduct,
                        [
                            ExprBuilder(_make_matrix, [ptr1]),
                            self.exp*hadamard_power(self.base, self.exp-1),
                            ExprBuilder(_make_matrix, [ptr2]),
                        ]
                    ),
                    *diagonal
                ],
                validator=CodegenArrayDiagonal._validate
            )
            i._first_pointer_parent = subexpr.args[0].args[0].args
            i._first_pointer_index = 0
            i._second_pointer_parent = subexpr.args[0].args[2].args
            i._second_pointer_index = 2
            i._lines = [subexpr]
        return lr
