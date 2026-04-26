"""
AST nodes for LAPACK function calls.

This module provides specialized
symbolic expression nodes that represents standard LAPACK operations
within matrix expressions.
"""

from __future__ import annotations

from .ast import Token
from sympy import Expr
from sympy.core.sympify import sympify

class Dgesv(Token, Expr):
    """
    Represents a symbolic node for the DGESV LAPACK operation.

    Parameters
    ==========

    matrix : MatrixSymbol

      Matrix representing the coefficients of variables in the linear
      equation. This matrix must be square and full-rank (i.e. all columns must
      be linearly independent) for the solving operation to be valid.

    vector : MatrixSymbol

      One-column matrix representing the solutions to the equations
      represented in ``matrix``.

    Examples
    ========

      >>> from sympy import MatrixSymbol
      >>> from sympy.codegen.lapack_nodes import Dgesv
      >>> A = MatrixSymbol('A', 3, 3)
      >>> b = MatrixSymbol('b', 3, 1)
      >>> Dgesv(A, b)
      Dgesv(A, vector=b)
      """

    __slots__ = _fields = ('matrix', 'vector')

    _construct_matrix = staticmethod(sympify)
    _construct_vector = staticmethod(sympify)

    @property
    def shape(self):
        return self.vector.shape
