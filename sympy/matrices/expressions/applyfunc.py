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
            ewdiff = diagonalize_vector(ewdiff)
            # TODO: check which axis is not 1
            for i in lr:
                if iscolumn:
                    ptr1 = i.first_pointer
                    ptr2 = Identity(ewdiff.shape[0])
                else:
                    ptr1 = Identity(ewdiff.shape[1])
                    ptr2 = i.second_pointer

                # TODO: check if pointers point to two different lines:

                def mul(*args):
                    return Mul.fromiter(args)

                def hadamard_or_mul(arg1, arg2):
                    if arg1.shape == arg2.shape:
                        return hadamard_product(arg1, arg2)
                    elif arg1.shape[1] == arg2.shape[0]:
                        return MatMul(arg1, arg2).doit()
                    elif arg1.shape[0] == arg2.shape[0]:
                        return MatMul(arg2.T, arg1).doit()
                    raise NotImplementedError

                i._lines = [[hadamard_or_mul, [[mul, [ewdiff, ptr1]], ptr2]]]
                i._first_pointer_parent = i._lines[0][1][0][1]
                i._first_pointer_index = 1
                i._second_pointer_parent = i._lines[0][1]
                i._second_pointer_index = 1

        else:
            # Matrix case:
            for i in lr:
                ptr1 = i.first_pointer
                ptr2 = i.second_pointer
                newptr1 = Identity(ptr1.shape[1])
                newptr2 = Identity(ptr2.shape[1])
                subexpr1 = ExprBuilder(
                    MatMul,
                    [ptr1, ExprBuilder(diagonalize_vector, [newptr1])],
                    validator=matmul_validate,
                )
                subexpr2 = ExprBuilder(
                    Transpose,
                    [ExprBuilder(
                        MatMul,
                        [
                            ptr2,
                            ExprBuilder(diagonalize_vector, [newptr2])
                            ,
                        ],
                    )],
                    validator=matmul_validate,
                )
                # TODO: replace the expressions above with the following:
                # subexpr = ExprBuilder(
                    # get_recognize([(1, 2, 8), (5, 6, 9)]),
                    # [ptr1, newptr1, ptr2, newptr2, ewdiff]
                # )
                # ... = recognize_matrix_expression(subexpr.build()).doit()
                i.first_pointer = subexpr1
                i.second_pointer = subexpr2
                i._first_pointer_parent = subexpr1.args[1].args
                i._first_pointer_index = 0
                i._second_pointer_parent = subexpr2.args[0].args[1].args
                i._second_pointer_index = 0
                # TODO: check if pointers point to two different lines:

                # Unify lines:
                l = i._lines
                # TODO: check nested fucntions, e.g. log(sin(...)), the second function should be a scalar one.
                i._lines = [ExprBuilder(MatMul, [l[0], ewdiff, l[1]], validator=matmul_validate)]
        return lr
