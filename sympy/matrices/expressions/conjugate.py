from __future__ import annotations

from sympy.core.sympify import _sympify
from sympy.matrices.expressions.matexpr import MatrixExpr


class Conjugate(MatrixExpr):
    """
    The elementwise complex conjugate of a matrix expression.

    This is a symbolic object that simply stores its argument without
    evaluating it. To actually compute the conjugate, use the ``conjugate()``
    function, or the ``.conjugate()`` method of matrices.

    Examples
    ========

    >>> from sympy import MatrixSymbol, Conjugate, conjugate
    >>> A = MatrixSymbol('A', 3, 5)
    >>> Conjugate(A)
    Conjugate(A)
    >>> Conjugate(A).shape
    (3, 5)
    >>> Conjugate(A).T
    Adjoint(A)
    >>> Conjugate(A).adjoint()
    A.T
    >>> Conjugate(Conjugate(A))
    A
    """
    is_Conjugate = True

    def __new__(cls, mat):
        mat = _sympify(mat)
        if not mat.is_Matrix:
            raise TypeError("argument to Conjugate must be a matrix")
        return super().__new__(cls, mat)

    def doit(self, **hints):
        from sympy.functions.elementary.complexes import conjugate
        arg = self.arg.doit(**hints) if hints.get('deep', True) else self.arg
        return conjugate(arg)

    @property
    def arg(self):
        return self.args[0]

    @property
    def shape(self):
        return self.arg.shape

    def _entry(self, i, j, **kwargs):
        from sympy.functions.elementary.complexes import conjugate
        return conjugate(self.arg._entry(i, j, **kwargs))

    def _eval_adjoint(self):
        return self.arg.transpose()

    def _eval_transpose(self):
        return self.arg.adjoint()

    def _eval_conjugate(self):
        return self.arg

    def _eval_trace(self):
        from sympy.matrices.expressions.trace import Trace
        from sympy.functions.elementary.complexes import conjugate
        return conjugate(Trace(self.arg))

    def _eval_determinant(self):
        from sympy.matrices.expressions.determinant import det
        from sympy.functions.elementary.complexes import conjugate
        return conjugate(det(self.arg))

    def _eval_derivative(self, x):
        return None


from sympy.assumptions.ask import ask, Q
from sympy.assumptions.refine import handlers_dict


def refine_Conjugate(expr, assumptions):
    """
    Refine a Conjugate matrix expression under assumptions.

    Examples
    ========

    >>> from sympy import MatrixSymbol, Q, assuming, refine
    >>> from sympy.matrices.expressions.conjugate import Conjugate
    >>> X = MatrixSymbol('X', 2, 2)
    >>> with assuming(Q.real_elements(X)):
    ...     print(refine(Conjugate(X)))
    X
    """
    if ask(Q.real_elements(expr.arg), assumptions):
        return expr.arg

    return expr


handlers_dict['Conjugate'] = refine_Conjugate

