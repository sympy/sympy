from sympy import Basic, Expr, S, sympify
from sympy.matrices.common import NonSquareMatrixError


class Determinant(Expr):
    """Matrix Determinant

    Represents the determinant of a matrix expression.

    Examples
    ========

    >>> from sympy import MatrixSymbol, Determinant, eye
    >>> A = MatrixSymbol('A', 3, 3)
    >>> Determinant(A)
    Determinant(A)
    >>> Determinant(eye(3)).doit()
    1
    """
    is_commutative = True

    def __new__(cls, mat):
        mat = sympify(mat)
        if not mat.is_Matrix:
            raise TypeError("Input to Determinant, %s, not a matrix" % str(mat))

        if not mat.is_square:
            raise NonSquareMatrixError("Det of a non-square matrix")

        return Basic.__new__(cls, mat)

    @property
    def arg(self):
        return self.args[0]

    def doit(self, expand=False):
        try:
            return self.arg._eval_determinant()
        except (AttributeError, NotImplementedError):
            return self

def det(matexpr):
    """ Matrix Determinant

    Examples
    ========

    >>> from sympy import MatrixSymbol, det, eye
    >>> A = MatrixSymbol('A', 3, 3)
    >>> det(A)
    Determinant(A)
    >>> det(eye(3))
    1
    """

    return Determinant(matexpr).doit()

class Permanent(Expr):
    """Matrix Permanent

    Represents the permanent of a matrix expression.

    Examples
    ========

    >>> from sympy import MatrixSymbol, Permanent, ones
    >>> A = MatrixSymbol('A', 3, 3)
    >>> Permanent(A)
    Permanent(A)
    >>> Permanent(ones(3, 3)).doit()
    6
    """

    def __new__(cls, mat):
        mat = sympify(mat)
        if not mat.is_Matrix:
            raise TypeError("Input to Permanent, %s, not a matrix" % str(mat))

        return Basic.__new__(cls, mat)

    @property
    def arg(self):
        return self.args[0]

    def doit(self, expand=False):
        try:
            return self.arg.per()
        except (AttributeError, NotImplementedError):
            return self

def per(matexpr):
    """ Matrix Permanent

    Examples
    ========

    >>> from sympy import MatrixSymbol, Matrix, per, ones
    >>> A = MatrixSymbol('A', 3, 3)
    >>> per(A)
    Permanent(A)
    >>> per(ones(5, 5))
    120
    >>> M = Matrix([1, 2, 5])
    >>> per(M)
    8
    """

    return Permanent(matexpr).doit()

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
