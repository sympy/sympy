import itertools

from sympy.core.sympify import _sympify

from sympy import MatrixExpr, Expr, ShapeError, ZeroMatrix, Add, Mul, MatMul
from sympy.stats.rv import RandomSymbol
from sympy.stats.symbolic_probability import Variance, Covariance, Expectation, is_random


class ExpectationMatrix(Expectation, MatrixExpr):
    """
    Expectation of a matrix expression.
    """
    def __new__(cls, arg, condition=None):
        arg = _sympify(arg)

        if condition:
            obj = Expr.__new__(cls, arg, condition)
        else:
            obj = Expr.__new__(cls, arg)
        obj._shape = arg.shape
        obj._condition = condition
        return obj

    @property
    def shape(self):
        return self._shape

    def doit(self, **hints):
        expr = self.args[0]
        condition = self._condition

        if not is_random(expr):
            return expr

        if isinstance(expr, Add):
            return Add(*[Expectation(a, condition=condition).doit() for a in expr.args])
        elif isinstance(expr, (Mul, MatMul)):
            rv = []
            nonrv = []
            postnon = []

            for a in expr.args:
                if is_random(a):
                    if rv:
                        rv.extend(postnon)
                    else:
                        nonrv.extend(postnon)
                    postnon = []
                    rv.append(a)
                elif a.is_Matrix:
                    postnon.append(a)
                else:
                    nonrv.append(a)

            # In order to avoid infinite-looping (MatMul may call .doit() again),
            # do not rebuild
            if len(nonrv) == 0:
                return self
            return Mul.fromiter(nonrv)*Expectation(Mul.fromiter(rv), condition=condition)*Mul.fromiter(postnon)

        return self


class VarianceMatrix(Variance, MatrixExpr):
    """
    Variance of a matrix probability expression.
    """
    def __new__(cls, arg, condition=None):
        arg = _sympify(arg)

        if 1 not in arg.shape:
            raise ShapeError("Expression is not a vector")

        shape = (arg.shape[0], arg.shape[0]) if arg.shape[1] == 1 else (arg.shape[1], arg.shape[1])

        if condition:
            obj = Expr.__new__(cls, arg, condition)
        else:
            obj = Expr.__new__(cls, arg)

        obj._shape = shape
        obj._condition = condition
        return obj

    @property
    def shape(self):
        return self._shape

    def doit(self, **hints):
        arg = self.args[0]
        condition = self._condition

        if not is_random(arg):
            return ZeroMatrix(*self.shape)

        if isinstance(arg, RandomSymbol):
            return self
        elif isinstance(arg, Add):
            rv = []
            for a in arg.args:
                if a.has(RandomSymbol):
                    rv.append(a)
            variances = Add(*map(lambda xv: Variance(xv, condition).doit(), rv))
            map_to_covar = lambda x: 2*Covariance(*x, condition=condition).doit()
            covariances = Add(*map(map_to_covar, itertools.combinations(rv, 2)))
            return variances + covariances
        elif isinstance(arg, Mul):
            nonrv = []
            rv = []
            postnon = []
            for a in arg.args:
                if is_random(a):
                    rv.append(a)
                else:
                    nonrv.append(a)
                    postnon.append(a)
            if len(rv) == 0:
                return ZeroMatrix(*self.shape)
            # Avoid possible infinite loops with MatMul:
            if len(nonrv) == 0 and len(postnon) == 0:
                return self
            # Variance of many multiple matrix products is not implemented:
            if len(rv) > 1:
                return self
            return Mul.fromiter(nonrv)*Variance(Mul.fromiter(rv), condition)*(Mul.fromiter(postnon)).transpose()

        # this expression contains a RandomSymbol somehow:
        return self
