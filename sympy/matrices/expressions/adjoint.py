from sympy.matrices.expressions.matexpr import MatrixExpr
from sympy.core import Basic
from sympy.functions import conjugate, adjoint

class Adjoint(MatrixExpr):
    """
    The Hermitian adjoint of a matrix expression.

    Examples
    ========

    >>> from sympy.matrices import MatrixSymbol, Adjoint
    >>> from sympy.functions import adjoint
    >>> A = MatrixSymbol('A', 3, 5)
    >>> B = MatrixSymbol('B', 5, 3)
    >>> Adjoint(A*B)
    Adjoint(A*B)
    >>> adjoint(A*B)
    Adjoint(B)*Adjoint(A)
    >>> adjoint(A*B) == Adjoint(A*B)
    False
    >>> adjoint(A*B) == Adjoint(A*B).doit()
    True
    """
    is_Adjoint = True

    def doit(self, **hints):
        arg = self.arg
        if hints.get('deep', True) and isinstance(arg, Basic):
            return adjoint(arg.doit(**hints))
        else:
            return adjoint(self.arg)

    @property
    def arg(self):
        return self.args[0]

    @property
    def shape(self):
        return self.arg.shape[::-1]

    def _entry(self, i, j):
        return conjugate(self.arg._entry(j, i))

    def _eval_adjoint(self):
        return self.arg

    def _eval_trace(self):
        from sympy.matrices.expressions.trace import Trace
        return conjugate(Trace(self.arg))
