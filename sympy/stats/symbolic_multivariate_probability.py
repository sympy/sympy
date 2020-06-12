import itertools

from sympy import MatrixExpr, Expr, ShapeError, ZeroMatrix, Add, Mul, MatMul, S
from sympy.stats.rv import RandomSymbol
from sympy.core.sympify import _sympify
from sympy.stats.symbolic_probability import Variance, Covariance, Expectation, is_random

class ExpectationMatrix(Expectation, MatrixExpr):
    """
    Expectation of a random matrix expression.

    Examples
    ========

    >>> from sympy.stats import ExpectationMatrix
    >>> from sympy.stats.rv import RandomMatrixSymbol
    >>> from sympy import symbols, MatrixSymbol
    >>> k = symbols("k")
    >>> A, B = MatrixSymbol("A", k, k), MatrixSymbol("B", k, k)
    >>> X, Y = RandomMatrixSymbol("X", k, 1), RandomMatrixSymbol("Y", k, 1)
    >>> ExpectationMatrix(X)
    ExpectationMatrix(X)

    """
    def __new__(cls, expr, condition=None):
        expr = _sympify(expr)
        if condition is None:
            if not is_random(expr):
                return expr
            obj = Expr.__new__(cls, expr)
        else:
            condition = _sympify(condition)
            obj = Expr.__new__(cls, expr, condition)

        obj._shape = expr.shape
        obj._condition = condition
        return obj

    @property
    def shape(self):
        return self._shape

    def expand(self, **hints):
        expr = self.args[0]
        condition = self._condition
        if not is_random(expr):
            return expr

        if isinstance(expr, Add):
            return Add(*[Expectation(a, condition=condition).expand() for a in expr.args])
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
            return Mul.fromiter(nonrv)*Expectation(Mul.fromiter(rv),
                    condition=condition)*Mul.fromiter(postnon)

        return self

class VarianceMatrix(Variance, MatrixExpr):
    """
    Variance of a random matrix probability expression. Also known as
    Covariance matrix, auto-covariance matrix, dispersion matrix,
    or variance-covariance matrix

    Examples
    ========

    >>> from sympy.stats import VarianceMatrix
    >>> from sympy.stats.rv import RandomMatrixSymbol
    >>> from sympy import symbols, MatrixSymbol
    >>> k = symbols("k")
    >>> A, B = MatrixSymbol("A", k, k), MatrixSymbol("B", k, k)
    >>> X, Y = RandomMatrixSymbol("X", k, 1), RandomMatrixSymbol("Y", k, 1)
    >>> VarianceMatrix(X)
    VarianceMatrix(X)
    >>> VarianceMatrix(X).shape
    (k, k)

    To expand the variance in its expression, use ``expand()``:

    >>> VarianceMatrix(A*X).expand()
    A*VarianceMatrix(X)*A.T
    >>> VarianceMatrix(A*X + B*Y).expand()
    2*A*CrossCovarianceMatrix(X, Y)*B.T + A*VarianceMatrix(X)*A.T + B*VarianceMatrix(Y)*B.T
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

    def expand(self, **hints):
        arg = self.args[0]
        condition = self._condition

        if not is_random(arg):
            return ZeroMatrix(*self.shape)

        if isinstance(arg, RandomSymbol):
            return self
        elif isinstance(arg, Add):
            rv = []
            for a in arg.args:
                if is_random(a):
                    rv.append(a)
            variances = Add(*map(lambda xv: Variance(xv, condition).expand(), rv))
            map_to_covar = lambda x: 2*Covariance(*x, condition=condition).expand()
            covariances = Add(*map(map_to_covar, itertools.combinations(rv, 2)))
            return variances + covariances
        elif isinstance(arg, (Mul, MatMul)):
            nonrv = []
            rv = []
            for a in arg.args:
                if is_random(a):
                    rv.append(a)
                else:
                    nonrv.append(a)
            if len(rv) == 0:
                return ZeroMatrix(*self.shape)
            # Avoid possible infinite loops with MatMul:
            if len(nonrv) == 0:
                return self
            # Variance of many multiple matrix products is not implemented:
            if len(rv) > 1:
                return self
            return Mul.fromiter(nonrv)*Variance(Mul.fromiter(rv),
                            condition)*(Mul.fromiter(nonrv)).transpose()

        # this expression contains a RandomSymbol somehow:
        return self

class CrossCovarianceMatrix(Covariance, MatrixExpr):
    """
    Covariance of a random matrix probability expression.

    Examples
    ========

    >>> from sympy.stats import CrossCovarianceMatrix
    >>> from sympy.stats.rv import RandomMatrixSymbol
    >>> from sympy import symbols, MatrixSymbol
    >>> k = symbols("k")
    >>> A, B = MatrixSymbol("A", k, k), MatrixSymbol("B", k, k)
    >>> C, D = MatrixSymbol("C", k, k), MatrixSymbol("D", k, k)
    >>> X, Y = RandomMatrixSymbol("X", k, 1), RandomMatrixSymbol("Y", k, 1)
    >>> Z, W = RandomMatrixSymbol("Z", k, 1), RandomMatrixSymbol("W", k, 1)
    >>> CrossCovarianceMatrix(X, Y)
    CrossCovarianceMatrix(X, Y)
    >>> CrossCovarianceMatrix(X, Y).shape
    (k, k)

    To expand the variance in its expression, use ``expand()``:

    >>> CrossCovarianceMatrix(X + Y, Z).expand()
    CrossCovarianceMatrix(X, Z) + CrossCovarianceMatrix(Y, Z)
    >>> CrossCovarianceMatrix(A*X , Y).expand()
    A*CrossCovarianceMatrix(X, Y)
    >>> CrossCovarianceMatrix(A*X, B.T*Y).expand()
    A*CrossCovarianceMatrix(X, Y)*B
    >>> CrossCovarianceMatrix(A*X + B*Y, C.T*Z + D.T*W).expand()
    A*CrossCovarianceMatrix(X, W)*D + A*CrossCovarianceMatrix(X, Z)*C + B*CrossCovarianceMatrix(Y, W)*D + B*CrossCovarianceMatrix(Y, Z)*C

    """
    def __new__(cls, arg1, arg2, condition=None):
        arg1 = _sympify(arg1)
        arg2 = _sympify(arg2)

        if (1 not in arg1.shape) or (1 not in arg2.shape) or (arg1.shape[1] != arg2.shape[1]):
            raise ShapeError("Expression is not a vector")

        shape = (arg1.shape[0], arg2.shape[0]) if arg1.shape[1] == 1 and arg2.shape[1] == 1 \
                    else (1, 1)

        if condition:
            obj = Expr.__new__(cls, arg1, arg2, condition)
        else:
            obj = Expr.__new__(cls, arg1, arg2)

        obj._shape = shape
        obj._condition = condition
        return obj

    @property
    def shape(self):
        return self._shape

    def expand(self, **hints):
        arg1 = self.args[0]
        arg2 = self.args[1]
        condition = self._condition

        if arg1 == arg2:
            return VarianceMatrix(arg1, condition).expand()

        if not is_random(arg1) or not is_random(arg2):
            return ZeroMatrix(*self.shape)

        if isinstance(arg1, RandomSymbol) and isinstance(arg2, RandomSymbol):
            return CrossCovarianceMatrix(arg1, arg2, condition)

        coeff_rv_list1 = self._expand_single_argument(arg1.expand())
        coeff_rv_list2 = self._expand_single_argument(arg2.expand())

        addends = [a*CrossCovarianceMatrix(r1, r2, condition=condition)*b.transpose()
                   for (a, r1) in coeff_rv_list1 for (b, r2) in coeff_rv_list2]
        return Add(*addends)

    @classmethod
    def _expand_single_argument(cls, expr):
        # return (coefficient, random_symbol) pairs:
        if isinstance(expr, RandomSymbol):
            return [(S.One, expr)]
        elif isinstance(expr, Add):
            outval = []
            for a in expr.args:
                if isinstance(a, (Mul, MatMul)):
                    outval.append(cls._get_mul_nonrv_rv_tuple(a))
                elif is_random(a):
                    outval.append((S.One, a))

            return outval
        elif isinstance(expr, (Mul, MatMul)):
            return [cls._get_mul_nonrv_rv_tuple(expr)]
        elif is_random(expr):
            return [(S.One, expr)]

    @classmethod
    def _get_mul_nonrv_rv_tuple(cls, m):
        rv = []
        nonrv = []
        for a in m.args:
            if is_random(a):
                rv.append(a)
            else:
                nonrv.append(a)
        return (Mul(*nonrv), Mul(*rv))
