"""
AST nodes for LAPACK function calls.

This module provides specialized
symbolic expression nodes that represents standard LAPACK operations
within matrix expressions.
"""

from __future__ import annotations
from sympy.codegen.fnodes import dimension

from .ast import Token
from sympy.core.sympify import sympify
from sympy.codegen.ast import (Declaration, Variable, integer, Return,
                              real, FunctionPrototype,
                              FunctionDefinition, CodeBlock)
from sympy import symbols

class Dgesv(Token):
    """
    Represents a symbolic node for the DGESV LAPACK operation.

    Parameters
    ==========

    matrix : MatrixSymbol

      Matrix representing the coefficients of variables in the linear
      equation. This matrix must be square and full-rank (i.e. all
      columns must be linearly independent) for the solving operation
      to be valid.

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

def dgesv_function(matrix, vector, func_name="solve"):
      """
      Generates an AST for a function that calls LAPACK's dgesv.

      Parameters
      ==========

      matrix : MatrixSymbol
        Matrix representing the coefficients of variables in the linear
        equation. This matrix must be square and full-rank (i.e. all
        columns must be linearly independent) for the solving operation
        to be valid.

      vector : MatrixSymbol
        One-column matrix representing the solutions to the equations
        represented in ``matrix``.

      func_name : str
        Name of the generated function.

      Example
      =======

      >>> from sympy import MatrixSymbol, ccode
      >>> from sympy.codegen.lapack_nodes import dgesv_function
      >>> A = MatrixSymbol('A', 3, 3)
      >>> b = MatrixSymbol('b', 3, 1)
      >>> func = dgesv_function(A, b)
      >>> 'LAPACKE_dgesv' in ccode(func)
      True
      >>> from sympy.printing.c import C99CodePrinter
      >>> printer = C99CodePrinter()
      >>> _ = printer.doprint(func)
      >>> 'lapacke.h' in printer.headers
      True
      """
      n_sym = symbols('n')
      nrhs_sym = symbols('nrhs')
      ipiv_sym = symbols('ipiv')
      info_sym = symbols('info')
      n_dec = Declaration(Variable(n_sym, type=integer, value=matrix.shape[0]))
      nrhs_dec = Declaration(Variable(nrhs_sym, type=integer, value=vector.shape[1]))
      ipiv_dec = Declaration(Variable(ipiv_sym, type=integer, attrs=[dimension(matrix.shape[0])]))
      info_dec = Declaration(Variable(info_sym, type=integer))
      dgesv_call = Dgesv(matrix, vector)
      info_return = Return(info_sym)
      fp = FunctionPrototype(integer, func_name, [Variable(matrix.name, type=real), Variable(vector.name, type=real)])
      body = CodeBlock(n_dec, nrhs_dec, ipiv_dec, info_dec, dgesv_call, info_return)
      return FunctionDefinition.from_FunctionPrototype(fp, body)
