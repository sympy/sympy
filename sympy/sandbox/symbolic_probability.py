import itertools

from sympy.core.sympify import _sympify

from sympy.core.compatibility import default_sort_key

from sympy import Expr, Add, Mul, S
from sympy.core.evaluate import global_evaluate
from sympy.stats import variance, covariance
from sympy.stats.rv import RandomSymbol


class Probability(Expr):
    """
    Symbolic expression for the probability.
    """
    def __new__(cls, prob, condition=None, **kwargs):
        pass

class Variance(Expr):
    """
    Symbolic expression for the variance.
    """
    def __new__(cls, arg, condition=None, **kwargs):

        if not kwargs.pop('evaluate', global_evaluate[0]):
            if condition is None:
                obj = Expr.__new__(cls, arg, **kwargs)
            else:
                obj = Expr.__new__(cls, arg, _sympify(condition), **kwargs)
            obj._condition = condition
            return obj

        if not arg.has(RandomSymbol):
            return S.Zero

        if isinstance(arg, RandomSymbol):
            if condition is None:
                obj = Expr.__new__(cls, arg, **kwargs)
            else:
                obj = Expr.__new__(cls, arg, _sympify(condition), **kwargs)
            obj._condition = condition
            return obj
        elif isinstance(arg, Add):
            rv = []
            for a in arg.args:
                if a.has(RandomSymbol):
                    rv.append(a)
            variances = Add(*map(lambda xv: Variance(xv, condition), rv))
            map_to_covar = lambda x: 2*Covariance(*x, condition=condition)
            covariances = Add(*map(map_to_covar, itertools.combinations(rv, 2)))
            return variances + covariances
        elif isinstance(arg, Mul):
            nonrv = []
            rv = []
            for a in arg.args:
                if a.has(RandomSymbol):
                    rv.append(a)
                else:
                    nonrv.append(a**2)
            if len(rv) == 0:
                return S.Zero
            if condition is None:
                obj = Expr.__new__(cls, Mul(*rv), **kwargs)
            else:
                obj = Expr.__new__(cls, Mul(*rv), condition, **kwargs)
            obj._condition = condition
            return Mul(*nonrv)*obj

        raise NotImplementedError

    def doit(self, **kwargs):
        return variance(self.args[0], self._condition, **kwargs)


class Covariance(Expr):
    """
    Symbolic expression for the covariance.
    """
    def __new__(cls, arg1, arg2, condition=None, **kwargs):

        if not kwargs.pop('evaluate', global_evaluate[0]):
            if condition is None:
                obj = Expr.__new__(cls, arg1, arg2, **kwargs)
            else:
                obj = Expr.__new__(cls, arg1, arg2, _sympify(condition), **kwargs)
            obj._condition = condition
            return obj

        if arg1 == arg2:
            return Variance(arg1, condition)

        if not arg1.has(RandomSymbol):
            return S.Zero
        if not arg2.has(RandomSymbol):
            return S.Zero

        arg1, arg2 = sorted([arg1, arg2], key=default_sort_key)

        if isinstance(arg1, RandomSymbol) and isinstance(arg2, RandomSymbol):
            return Expr.__new__(cls, arg1, arg2, **kwargs)

        coeff_rv_list1 = cls._expand_single_argument(arg1.expand())
        coeff_rv_list2 = cls._expand_single_argument(arg2.expand())

        addends = [a*b*Covariance(r1, r2) for (a, r1) in coeff_rv_list1 for (b, r2) in coeff_rv_list2]
        return Add(*addends)

    @classmethod
    def _expand_single_argument(cls, expr):
        # return (coefficient, random_symbol) pairs:
        if isinstance(expr, RandomSymbol):
            return [(S.One, expr)]
        elif isinstance(expr, Add):
            outval = []
            for a in expr.args:
                rv = []
                nonrv = []
                if isinstance(a, Mul):
                    outval.append(cls._get_mul_nonrv_rv_tuple(a))
                elif isinstance(a, RandomSymbol):
                    outval.append((S.One, a))

                outval.append((Mul(*nonrv), Mul(*rv)))
            return outval
        elif isinstance(expr, Mul):
            return [cls._get_mul_nonrv_rv_tuple(expr)]

    @classmethod
    def _get_mul_nonrv_rv_tuple(cls, m):
        rv = []
        nonrv = []
        for a in m.args:
            if a.has(RandomSymbol):
                rv.append(a)
            else:
                nonrv.append(a)
        return (Mul(*nonrv), Mul(*rv))

    def doit(self, **kwargs):
        return covariance(self.args[0], self.args[1], self._condition, **kwargs)
