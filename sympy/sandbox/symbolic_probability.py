import itertools

from sympy.core.compatibility import default_sort_key

from sympy import Expr, Add, Mul, S
from sympy.core.evaluate import global_evaluate
from sympy.stats import variance, covariance
from sympy.stats.rv import RandomSymbol


class Variance(Expr):
    """
    Symbolic expression for the variance.
    """
    def __new__(cls, arg, **kwargs):

        if not kwargs.pop('evaluate', global_evaluate[0]):
            return Expr.__new__(cls, arg, **kwargs)

        if not arg.has(RandomSymbol):
            return S.Zero

        if isinstance(arg, RandomSymbol):
            return Expr.__new__(cls, arg, **kwargs)
        elif isinstance(arg, Add):
            rv = []
            for a in arg.args:
                if a.has(RandomSymbol):
                    rv.append(a)
            variances = Add(*map(Variance, rv))
            map_to_covar = lambda x: 2*Covariance(*x)
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
            obj = Expr.__new__(cls, Mul(*rv), **kwargs)
            return Mul(*nonrv)*obj

        raise NotImplementedError

    def doit(self, **kwargs):
        return variance(self.args[0], **kwargs)


class Covariance(Expr):
    """
    Symbolic expression for the covariance.
    """
    def __new__(cls, arg1, arg2, **kwargs):

        if not kwargs.pop('evaluate', global_evaluate[0]):
            return Expr.__new__(cls, arg1, arg2, **kwargs)

        if arg1 == arg2:
            return Variance(arg1)

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
        return covariance(self.args[0], self.args[1], **kwargs)
