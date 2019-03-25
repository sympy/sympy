from sympy.matrices.expressions import MatrixExpr
from sympy import MatrixBase, Dummy, Lambda, Function, DiagonalMatrix, DiagonalOf
from sympy.matrices.expressions.diagonal import DiagonalizeVector, diagonalize_vector


class ElementwiseApplyFunction(MatrixExpr):
    r"""
    Apply function to a matrix elementwise without evaluating.

    Examples
    ========

    >>> from sympy.matrices.expressions import MatrixSymbol
    >>> from sympy.matrices.expressions.applyfunc import ElementwiseApplyFunction
    >>> from sympy import exp
    >>> X = MatrixSymbol("X", 3, 3)
    >>> X.applyfunc(exp)
    ElementwiseApplyFunction(exp, X)

    >>> from sympy import eye
    >>> expr = ElementwiseApplyFunction(exp, eye(3))
    >>> expr
    ElementwiseApplyFunction(exp, Matrix([
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]]))
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
        obj._function = function
        obj._expr = expr
        return obj

    @property
    def function(self):
        return self._function

    @property
    def expr(self):
        return self._expr

    @property
    def shape(self):
        return self.expr.shape

    def func(self, expr):
        return ElementwiseApplyFunction(self.function, expr)

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
        return self.function(self.expr[i, j])

    def _eval_derivative_matrix_lines(self, x):
        d = Dummy("d")
        function = self.function(d)
        fdiff = function.fdiff()
        if isinstance(fdiff, Function):
            fdiff = type(fdiff)
        else:
            fdiff = Lambda(d, function.diff(d))
        lr = self.expr._eval_derivative_matrix_lines(x)
        if 1 in self.shape:
            # Vector:
            ewdiff = ElementwiseApplyFunction(fdiff, self.expr)
            ewdiff = diagonalize_vector(ewdiff)
            # is it a vector or a matrix or a scalar?
            lr[0].first *= ewdiff
            return lr
        else:
            # Matrix case:
            raise NotImplementedError
