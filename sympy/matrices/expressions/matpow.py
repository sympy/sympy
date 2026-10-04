from __future__ import annotations
from .matexpr import MatrixExpr
from .special import Identity
from sympy.core import S
from sympy.core.cache import cacheit
from sympy.core.sympify import _sympify
from sympy.matrices import MatrixBase
from sympy.matrices.exceptions import NonSquareMatrixError


class MatPow(MatrixExpr):
    def __new__(cls, base, exp, evaluate=False, **options) -> MatrixExpr: # type: ignore
        base = _sympify(base)
        if not base.is_Matrix:
            raise TypeError("MatPow base should be a matrix")

        if base.is_square is False:
            raise NonSquareMatrixError("Power of non-square matrix %s" % base)

        exp = _sympify(exp)
        obj = super().__new__(cls, base, exp)

        if evaluate:
            obj = obj.doit(deep=False)

        return obj

    @property
    def base(self):
        return self.args[0]

    @property
    def exp(self):
        return self.args[1]

    @property
    def shape(self):
        return self.base.shape

    @cacheit
    def _get_explicit_matrix(self):
        return self.base.as_explicit()**self.exp

    def _entry(self, i, j, **kwargs):
        from sympy.matrices.expressions import MatMul
        A = self.doit()
        if isinstance(A, MatPow):
            # We still have a MatPow, make an explicit MatMul out of it.
            if A.exp.is_Integer and A.exp.is_positive:
                A = MatMul(*[A.base for k in range(A.exp)])
            elif not self._is_shape_symbolic():
                return A._get_explicit_matrix()[i, j]
            else:
                # Leave the expression unevaluated:
                from sympy.matrices.expressions.matexpr import MatrixElement
                return MatrixElement(self, i, j)
        return A[i, j]

    def doit(self, **hints):
        if hints.get('deep', True):
            base, exp = (arg.doit(**hints) for arg in self.args)
        else:
            base, exp = self.args

        # combine all powers, e.g. (A ** 2) ** 3 -> A ** 6
        while isinstance(base, MatPow):
            exp *= base.args[1]
            base = base.args[0]

        if isinstance(base, MatrixBase):
            # Delegate
            return base ** exp

        # Handle simple cases so that _eval_power() in MatrixExpr sub-classes can ignore them
        if exp == S.One:
            return base
        if exp == S.Zero:
            return Identity(base.rows)
        if exp == S.NegativeOne:
            from sympy.matrices.expressions import Inverse
            return Inverse(base).doit(**hints)

        eval_power = getattr(base, '_eval_power', None)
        if eval_power is not None:
            return eval_power(exp)

        return MatPow(base, exp)

    def _eval_transpose(self):
        base, exp = self.args
        return MatPow(base.transpose(), exp)

    def _eval_adjoint(self):
        base, exp = self.args
        return MatPow(base.adjoint(), exp)

    def _eval_conjugate(self):
        base, exp = self.args
        return MatPow(base.conjugate(), exp)

    def _eval_derivative(self, x):
        assert not hasattr(x, "shape")
        # dbase = self.base.diff(x)
        dexp = self.exp.diff(x)
        if dexp != 0:
            return None
        if self.exp.is_Integer:
            if self.exp > 0:
                from sympy import MatMul
                return MatMul.fromiter(self.base for _ in range(self.exp))._eval_derivative(x)
            elif self.exp == 0:
                from .special import ZeroMatrix
                return ZeroMatrix(*self.base.shape)
        return None

    def _eval_inverse(self):
        return MatPow(self.base, -self.exp)
