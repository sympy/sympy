from __future__ import print_function, division

from sympy import Basic, Q
from sympy.functions import adjoint, conjugate

from sympy.matrices.expressions.matexpr import MatrixExpr
from sympy.matrices import MatrixBase

class Transpose(MatrixExpr):
    """
    The transpose of a matrix expression.

    This is a symbolic object that simply stores its argument without
    evaluating it. To actually compute the transpose, use the ``transpose()``
    function, or the ``.T`` attribute of matrices.

    Examples
    ========

    >>> from sympy.matrices import MatrixSymbol, Transpose
    >>> from sympy.functions import transpose
    >>> A = MatrixSymbol('A', 3, 5)
    >>> B = MatrixSymbol('B', 5, 3)
    >>> Transpose(A)
    A'
    >>> A.T == transpose(A) == Transpose(A)
    True
    >>> Transpose(A*B)
    (A*B)'
    >>> transpose(A*B)
    B'*A'

    """
    is_Transpose = True

    def doit(self, **hints):
        arg = self.arg
        if hints.get('deep', True) and isinstance(arg, Basic):
            arg = arg.doit(**hints)
        try:
            result = arg._eval_transpose()
            return result if result is not None else Transpose(arg)
        except AttributeError:
            return Transpose(arg)

    @property
    def arg(self):
        return self.args[0]

    @property
    def shape(self):
        return self.arg.shape[::-1]

    def _entry(self, i, j):
        return self.arg._entry(j, i)

    def _eval_adjoint(self):
        return conjugate(self.arg)

    def _eval_conjugate(self):
        return adjoint(self.arg)

    def _eval_transpose(self):
        return self.arg

    def _eval_trace(self):
        from .trace import Trace
        return Trace(self.arg)  # Trace(X.T) => Trace(X)

    def _eval_determinant(self):
        from sympy.matrices.expressions.determinant import det
        return det(self.arg)

def transpose(expr):
    """ Matrix transpose """
    return Transpose(expr).doit()
