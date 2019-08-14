from __future__ import print_function, division

from .matexpr import MatrixExpr
from sympy.core.basic import Basic
from sympy.core.function import Lambda
from sympy.core.sympify import _sympify, sympify
from sympy.matrices import Matrix
from sympy.functions.elementary.complexes import re, im


class FunctionMatrix(MatrixExpr):
    """Represents a matrix using a function (``Lambda``) which gives
    outputs according to the coordinates of each matrix entries.

    Parameters
    ==========

    rows : nonnegative integer

    cols : nonnegative integer

    lamda : Lambda or str
        If it is a SymPy ``Lambda`` instance, it should have two
        arguments which represents the coordinate of the matrix.

        If it is a pure string containing python ``lambda`` semantics,
        it is interpreted by the SymPy parser and casted into a SymPy
        ``Lambda`` instance.

    Examples
    ========

    Creating a ``FunctionMatrix`` from ``Lambda``:

    >>> from sympy import FunctionMatrix, symbols, Lambda, MatPow, Matrix
    >>> i, j = symbols('i,j')
    >>> X = FunctionMatrix(3, 3, Lambda((i, j), i + j))

    Creating an explicit matrix from the ``FunctionMatrix``:

    >>> Matrix(X)
    Matrix([
    [0, 1, 2],
    [1, 2, 3],
    [2, 3, 4]])

    Creating a ``FunctionMatrix`` from python ``lambda``:

    >>> FunctionMatrix(3, 3, 'lambda i, j: i + j')
    FunctionMatrix(3, 3, Lambda((i, j), i + j))

    Example of lazy evaluation using the symbolic representation:

    >>> Y = FunctionMatrix(1000, 1000, Lambda((i, j), i + j))
    >>> isinstance(Y*Y, MatPow) # this is an expression object
    True
    >>> (Y**2)[10,10] # So this is evaluated lazily
    342923500

    Notes
    =====

    This class provides an alternative way to represent an extremely
    dense matrix with entries in some form of a sequence, in a most
    sparse way.
    """
    def __new__(cls, rows, cols, lamda):
        rows, cols = _sympify(rows), _sympify(cols)
        cls._check_dim(rows)
        cls._check_dim(cols)

        lamda = sympify(lamda)
        if not isinstance(lamda, Lambda):
            raise ValueError(
                "{} should be a SymPy Lambda instance.".format(lamda))
        if len(lamda.variables) != 2:
            raise ValueError(
                "{} should be a function of two variables.".format(lamda))

        return super(FunctionMatrix, cls).__new__(cls, rows, cols, lamda)

    @property
    def shape(self):
        return self.args[0:2]

    @property
    def lamda(self):
        return self.args[2]

    def _entry(self, i, j, **kwargs):
        return self.lamda(i, j)

    def _eval_trace(self):
        from sympy.matrices.expressions.trace import Trace
        from sympy import Sum
        return Trace(self).rewrite(Sum).doit()

    def as_real_imag(self):
        return (re(Matrix(self)), im(Matrix(self)))
