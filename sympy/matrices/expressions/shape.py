from sympy.core.relational import Eq
from sympy.core.expr import Expr
from sympy.logic.boolalg import Boolean, And
from sympy.matrices.expressions.matexpr import MatrixExpr
from typing import Union


def is_matadd_valid(*args: MatrixExpr) -> Boolean:
    rows, cols = zip(*(arg.shape for arg in args))
    return And(
        *(Eq(i, j) for i, j in zip(rows[:-1], rows[1:])),
        *(Eq(i, j) for i, j in zip(cols[:-1], cols[1:])),
    )


def is_matmul_valid(*args: Union[MatrixExpr, Expr]) -> Boolean:
    rows, cols = zip(*(arg.shape for arg in args if isinstance(arg, MatrixExpr)))
    return And(*(Eq(i, j) for i, j in zip(cols[:-1], rows[1:])))


def is_square(arg: MatrixExpr) -> Boolean:
    return Eq(arg.rows, arg.cols)
