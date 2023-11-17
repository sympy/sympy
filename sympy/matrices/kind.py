# sympy.matrices.kind

from sympy.core.kind import Kind, NumberKind

class MatrixKind(Kind):
    """
    Kind for all matrices in SymPy.

    Basic class for this kind is ``MatrixBase`` and ``MatrixExpr``,
    but any expression representing the matrix can have this.

    Parameters
    ==========

    element_kind : Kind
        Kind of the element. Default is
        :class:`sympy.core.kind.NumberKind`,
        which means that the matrix contains only numbers.

    Examples
    ========

    Any instance of matrix class has ``MatrixKind``:

    >>> from sympy import MatrixSymbol
    >>> A = MatrixSymbol('A', 2,2)
    >>> A.kind
    MatrixKind(NumberKind)

    Although expression representing a matrix may be not instance of
    matrix class, it will have ``MatrixKind`` as well:

    >>> from sympy import MatrixExpr, Integral
    >>> from sympy.abc import x
    >>> intM = Integral(A, x)
    >>> isinstance(intM, MatrixExpr)
    False
    >>> intM.kind
    MatrixKind(NumberKind)

    Use ``isinstance()`` to check for ``MatrixKind`` without specifying
    the element kind. Use ``is`` with specifying the element kind:

    >>> from sympy import Matrix
    >>> from sympy.core import NumberKind
    >>> from sympy.matrices import MatrixKind
    >>> M = Matrix([1, 2])
    >>> isinstance(M.kind, MatrixKind)
    True
    >>> M.kind is MatrixKind(NumberKind)
    True

    See Also
    ========

    sympy.core.kind.NumberKind
    sympy.core.kind.UndefinedKind
    sympy.core.containers.TupleKind
    sympy.sets.sets.SetKind

    """
    def __new__(cls, element_kind=NumberKind):
        obj = super().__new__(cls, element_kind)
        obj.element_kind = element_kind
        return obj

    def __repr__(self):
        return "MatrixKind(%s)" % self.element_kind
