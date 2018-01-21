from __future__ import print_function, division

from sympy import Basic, Expr, S, sympify
from .matexpr import ShapeError
from sympy.matrices.expressions.matexpr import MatrixExpr


class Determinant(Expr):
    """Matrix Determinant

    Represents the determinant of a matrix expression.

    >>> from sympy import MatrixSymbol, Determinant, eye
    >>> A = MatrixSymbol('A', 3, 3)
    >>> Determinant(A)
    Determinant(A)

    >>> Determinant(eye(3)).doit()
    1
    """

    def __new__(cls, mat):
        mat = sympify(mat)
        if not mat.is_Matrix:
            raise TypeError("Input to Determinant, %s, not a matrix" % str(mat))

        if not mat.is_square:
            raise ShapeError("Det of a non-square matrix")

        return Basic.__new__(cls, mat)

    @property
    def arg(self):
        return self.args[0]

    def doit(self, **kwargs):
        try:
            return self.arg._eval_determinant(**kwargs)
        except (AttributeError, NotImplementedError):
            return self

def det(matexpr, **kwargs):
    """ Matrix Determinant.

    Possible 'kwargs' are:
    method : string, optional (default='bareiss')
        specifies the algorithm to be used for computing the algorithm.
        It must be one of 'bareiss', 'lu' or 'berkowitz'.
        ``method`` is only applicable for a concrete matrix.

    >>> from sympy import MatrixSymbol, det, eye
    >>> A = MatrixSymbol('A', 3, 3)
    >>> det(A)
    Determinant(A)

    >>> det(eye(3))
    1
    """

    if isinstance(matexpr, MatrixExpr) and kwargs:
        raise ValueError(
            "Keyword arguments not supported for matrix expressions.")

    return Determinant(matexpr).doit(**kwargs)


from sympy.assumptions.ask import ask, Q
from sympy.assumptions.refine import handlers_dict


def refine_Determinant(expr, assumptions):
    """
    >>> from sympy import MatrixSymbol, Q, assuming, refine, det
    >>> X = MatrixSymbol('X', 2, 2)
    >>> det(X)
    Determinant(X)
    >>> with assuming(Q.orthogonal(X)):
    ...     print(refine(det(X)))
    1
    """
    if ask(Q.orthogonal(expr.arg), assumptions):
        return S.One
    elif ask(Q.singular(expr.arg), assumptions):
        return S.Zero
    elif ask(Q.unit_triangular(expr.arg), assumptions):
        return S.One

    return expr


handlers_dict['Determinant'] = refine_Determinant
