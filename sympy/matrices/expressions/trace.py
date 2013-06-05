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

    def doit(self, **kwargs):
        deep = kwargs.get('deep', False)
        from sympy import Sum, Dummy
        i = Dummy('i')
        rv = Sum(self.arg[i, i], (i, 0, self.arg.rows-1)).doit()
        if deep:
            return rv.doit(**kwargs)
        else:
            return rv

def trace(expr):
    """ Trace of a Matrix.  Sum of the diagonal elements

    >>> from sympy import trace, Symbol, MatrixSymbol, pprint, eye
    >>> n = Symbol('n')
    >>> X = MatrixSymbol('X', n, n)  # A square matrix
    >>> pprint(trace(X), use_unicode=False)
    n - 1
     __
     \ `
      )   X[i, i]
     /_,
    i = 0

    >>> trace(eye(3))
    3

    See Also:
        Trace
    """
    return Trace(expr).doit(deep=True)
