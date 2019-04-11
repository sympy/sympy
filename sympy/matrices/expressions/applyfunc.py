from sympy.matrices.expressions import MatrixExpr
from sympy import MatrixBase, Dummy, Lambda, Function, FunctionClass
from sympy.matrices.expressions.diagonal import diagonalize_vector


class ElementwiseApplyFunction(MatrixExpr):
    r"""
    Apply function to a matrix elementwise without evaluating.

    Examples
    ========

    It can be created by calling ``.applyfunc(<function>)`` on a matrix
    expression:

    >>> from sympy.matrices.expressions import MatrixSymbol
    >>> from sympy.matrices.expressions.applyfunc import ElementwiseApplyFunction
    >>> from sympy import exp
    >>> X = MatrixSymbol("X", 3, 3)
    >>> X.applyfunc(exp)
    exp(X...)

    Otherwise using the class constructor:

    >>> from sympy import eye
    >>> expr = ElementwiseApplyFunction(exp, eye(3))
    >>> expr
    exp(Matrix([
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]])...)
    >>> expr.doit()
    Matrix([
    [E, 1, 1],
    [1, E, 1],
    [1, 1, E]])

    Notice the difference with the real mathematical functions:

    >>> exp(eye(3))
    Matrix([
    [E, 0, 0],
    [0, E, 0],
    [0, 0, E]])
    """

    def __new__(cls, function, expr):
        obj = MatrixExpr.__new__(cls, expr)
        if not isinstance(function, FunctionClass):
            d = Dummy("d")
            function = Lambda(d, function(d))
        obj._function = function
        obj._expr = expr
        return obj

    def _hashable_content(self):
        return (self.function, self.expr)

    @property
    def function(self):
        return self._function

    @property
    def expr(self):
        return self._expr

    @property
    def shape(self):
        return self.expr.shape

    @property
    def func(self):
        # This strange construction is required by the assumptions:
        # (.func needs to be a class)

        class ElementwiseApplyFunction2(ElementwiseApplyFunction):
            def __new__(obj, expr):
                return ElementwiseApplyFunction(self.function, expr)

        return ElementwiseApplyFunction2

    def doit(self, **kwargs):
        deep = kwargs.get("deep", True)
        expr = self.expr
        if deep:
            expr = expr.doit(**kwargs)
        if isinstance(expr, MatrixBase):
            return expr.applyfunc(self.function)
        else:
            return self

    def _entry(self, i, j, **kwargs):
        return self.function(self.expr._entry(i, j, **kwargs))

    def _eval_derivative_matrix_lines(self, x):
        from sympy import HadamardProduct, hadamard_product, Mul, MatMul, Identity, Transpose
        from sympy.matrices.expressions.diagonal import diagonalize_vector
        from sympy.matrices.expressions.matmul import validate as matmul_validate
        from sympy.codegen.array_utils import CodegenArrayContraction, CodegenArrayTensorProduct, CodegenArrayDiagonal
        from sympy.core.expr import ExprBuilder

        d = Dummy("d")
        function = self.function(d)
        fdiff = function.fdiff()
        if isinstance(fdiff, Function):
            fdiff = type(fdiff)
        else:
            fdiff = Lambda(d, fdiff)
        lr = self.expr._eval_derivative_matrix_lines(x)
        ewdiff = ElementwiseApplyFunction(fdiff, self.expr)
        if 1 in x.shape:
            # Vector:
            iscolumn = self.shape[1] == 1
            for i in lr:
                if iscolumn:
                    ptr1 = i.first_pointer
                    ptr2 = Identity(ewdiff.shape[1])
                else:
                    ptr1 = Identity(ewdiff.shape[0])
                    ptr2 = i.second_pointer

                subexpr = ExprBuilder(
                    CodegenArrayDiagonal,
                    [
                        ExprBuilder(
                            CodegenArrayTensorProduct,
                            [
                                ptr1,
                                ewdiff,
                                ptr2
                            ]
                        ),
                        (0, 2), (3, 5)
                    ],
                    validator=CodegenArrayDiagonal._validate
                )
                i._lines = [subexpr]
                i._first_pointer_parent = subexpr.args[0].args
                i._first_pointer_index = 0
                i._second_pointer_parent = subexpr.args[0].args
                i._second_pointer_index = 2
        else:
            # Matrix case:
            for i in lr:
                ptr1 = i.first_pointer
                ptr2 = i.second_pointer
                newptr1 = Identity(ptr1.shape[1])
                newptr2 = Identity(ptr2.shape[1])
                subexpr = ExprBuilder(
                    CodegenArrayContraction,
                    [
                        ExprBuilder(
                            CodegenArrayTensorProduct,
                            [ptr1, newptr1, ptr2, newptr2, ewdiff]
                        ),
                        (1, 2, 8),
                        (5, 6, 9),
                    ],
                    validator=CodegenArrayContraction._validate
                )
                i._first_pointer_parent = subexpr.args[0].args
                i._first_pointer_index = 1
                i._second_pointer_parent = subexpr.args[0].args
                i._second_pointer_index = 3
                i._lines = [subexpr]
        return lr
