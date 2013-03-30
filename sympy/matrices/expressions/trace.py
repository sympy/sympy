from sympy import Basic, Expr
from matexpr import ShapeError


class Trace(Expr):
    """Matrix Trace

    Represents the trace of a matrix expression.

    >>> from sympy import MatrixSymbol, Trace, eye
    >>> A = MatrixSymbol('A', 3, 3)
    >>> Trace(A)
    Trace(A)

    >>> Trace(eye(3))
    3
    """
    is_Trace = True

    def __new__(cls, mat):
        if not mat.is_Matrix:
            raise TypeError("input to Trace, %s, is not a matrix" % str(mat))

        if not mat.is_square:
            raise ShapeError("Trace of a non-square matrix")

        try:
            return mat._eval_trace()
        except (AttributeError, NotImplementedError):
            return Basic.__new__(cls, mat)

    def _eval_transpose(self):
        return self

    @property
    def arg(self):
        return self.args[0]

    def doit(self):
        from sympy import Add
        return Add(*[self.arg[i, i] for i in range(self.arg.rows)])
