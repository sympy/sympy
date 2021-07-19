from typing import Optional

from sympy import Derivative, Integer, Expr
from sympy.matrices.common import MatrixCommon
from .ndim_array import NDimArray
from .arrayop import derive_by_array
from sympy import MatrixExpr
from sympy import ZeroMatrix
from sympy.matrices.expressions.matexpr import _matrix_derivative


class ArrayDerivative(Derivative):

    is_scalar = False

    def __new__(cls, expr, *variables, **kwargs):
        obj = super().__new__(cls, expr, *variables, **kwargs)
        if isinstance(obj, ArrayDerivative):
            obj._shape = obj._get_shape()
        return obj

    def _get_shape(self):
        shape = ()
        for v, count in self.variable_count:
            if hasattr(v, "shape"):
                for i in range(count):
                    shape += v.shape
        if hasattr(self.expr, "shape"):
            shape += self.expr.shape
        return shape

    @property
    def shape(self):
        return self._shape

    @classmethod
    def _get_zero_with_shape_like(cls, expr):
        if isinstance(expr, (MatrixCommon, NDimArray)):
            return expr.zeros(*expr.shape)
        elif isinstance(expr, MatrixExpr):
            return ZeroMatrix(*expr.shape)
        else:
            raise RuntimeError("Unable to determine shape of array-derivative.")

    @staticmethod
    def _call_derive_scalar_by_matrix(expr, v):  # type: (Expr, MatrixCommon) -> Expr
        return v.applyfunc(lambda x: expr.diff(x))

    @staticmethod
    def _call_derive_scalar_by_matexpr(expr, v):  # type: (Expr, MatrixExpr) -> Expr
        if expr.has(v):
            return _matrix_derivative(expr, v)
        else:
            return ZeroMatrix(*v.shape)

    @staticmethod
    def _call_derive_scalar_by_array(expr, v):  # type: (Expr, NDimArray) -> Expr
        return v.applyfunc(lambda x: expr.diff(x))

    @staticmethod
    def _call_derive_matrix_by_scalar(expr, v):  # type: (MatrixCommon, Expr) -> Expr
        return _matrix_derivative(expr, v)

    @staticmethod
    def _call_derive_matexpr_by_scalar(expr, v):  # type: (MatrixExpr, Expr) -> Expr
        return expr._eval_derivative(v)

    @staticmethod
    def _call_derive_array_by_scalar(expr, v):  # type: (NDimArray, Expr) -> Expr
        return expr.applyfunc(lambda x: x.diff(v))

    @staticmethod
    def _call_derive_default(expr, v):  # type: (Expr, Expr) -> Optional[Expr]
        if expr.has(v):
            return _matrix_derivative(expr, v)
        else:
            return None

    @classmethod
    def _dispatch_eval_derivative_n_times(cls, expr, v, count):
        # Evaluate the derivative `n` times.  If
        # `_eval_derivative_n_times` is not overridden by the current
        # object, the default in `Basic` will call a loop over
        # `_eval_derivative`:

        if not isinstance(count, (int, Integer)) or ((count <= 0) == True):
            return None

        # TODO: this could be done with multiple-dispatching:
        if expr.is_scalar:
            if isinstance(v, MatrixCommon):
                result = cls._call_derive_scalar_by_matrix(expr, v)
            elif isinstance(v, MatrixExpr):
                result = cls._call_derive_scalar_by_matexpr(expr, v)
            elif isinstance(v, NDimArray):
                result = cls._call_derive_scalar_by_array(expr, v)
            elif v.is_scalar:
                # scalar by scalar has a special
                return super()._dispatch_eval_derivative_n_times(expr, v, count)
            else:
                return None
        elif v.is_scalar:
            if isinstance(expr, MatrixCommon):
                result = cls._call_derive_matrix_by_scalar(expr, v)
            elif isinstance(expr, MatrixExpr):
                result = cls._call_derive_matexpr_by_scalar(expr, v)
            elif isinstance(expr, NDimArray):
                result = cls._call_derive_array_by_scalar(expr, v)
            else:
                return None
        else:
            # Both `expr` and `v` are some array/matrix type:
            if isinstance(expr, MatrixCommon) or isinstance(expr, MatrixCommon):
                result = derive_by_array(expr, v)
            elif isinstance(expr, MatrixExpr) and isinstance(v, MatrixExpr):
                result = cls._call_derive_default(expr, v)
            elif isinstance(expr, MatrixExpr) or isinstance(v, MatrixExpr):
                # if one expression is a symbolic matrix expression while the other isn't, don't evaluate:
                return None
            else:
                result = derive_by_array(expr, v)
        if result is None:
            return None
        if count == 1:
            return result
        else:
            return cls._dispatch_eval_derivative_n_times(result, v, count - 1)
