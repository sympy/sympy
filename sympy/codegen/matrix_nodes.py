"""
AST nodes for operations on matrices
"""

from .ast import Token
from sympy.matrices import MatrixExpr
from sympy.core.sympify import sympify


class MatrixSolve(Token, MatrixExpr):
    __slots__ = ['matrix', 'vector']

    _construct_matrix = staticmethod(sympify)

    def __init__(self, *args, **kwargs):
        self.shape = self.vector.shape
