import itertools

from sympy.core.sympify import _sympify

from sympy.core.compatibility import default_sort_key

from sympy import Expr, Add, Mul, S
from sympy.core.evaluate import global_evaluate
from sympy.stats import variance, covariance
from sympy.stats.rv import RandomSymbol, probability, expectation

__all__ = ['Probability', 'Expectation', 'Variance', 'Covariance']


class Probability(Expr):
    """
    Symbolic expression for the probability.

    Examples
    ========

    >>> from sympy.stats import Probability, Normal
    >>> X = Normal("X", 0, 1)
    >>> prob = Probability(X > 1)
    >>> prob
    Probability(X > 1)
    >>> prob.doit()
    sqrt(2)*(-sqrt(2)*sqrt(pi)*erf(sqrt(2)/2) + sqrt(2)*sqrt(pi))/(4*sqrt(pi))

    """
    def __new__(cls, prob, condition=None, **kwargs):
        prob = _sympify(prob)
        if condition is None:
            obj = Expr.__new__(cls, prob)
        else:
            condition = _sympify(condition)
            obj = Expr.__new__(cls, prob, condition)
        obj._condition = condition
        return obj

    def doit(self, **kwargs):
        return probability(self.args[0], given_condition=self._condition, **kwargs)


class Expectation(Expr):
    """
    Symbolic expression for the expectation.

    Examples
    ========

    >>> from sympy.stats import Expectation, Normal
    >>> from sympy import symbols
    >>> mu = symbols("mu")
    >>> sigma = symbols("sigma", positive=True)
    >>> X = Normal("X", mu, sigma)
    >>> Expectation(X)
    Expectation(X)
    >>> Expectation(X).doit().simplify()
    mu
    >>> Expectation(X).doit(evaluate=False)
    Integral(sqrt(2)*X*exp(-(X - mu)**2/(2*sigma**2))/(2*sqrt(pi)*sigma), (X, -oo, oo))

    This class is aware of some properties of the expectation:

    >>> from sympy.abc import a
    >>> Expectation(a*X)
    a*Expectation(X)
    >>> Y = Normal("Y", 0, 1)
    >>> Expectation(X + Y)
    Expectation(X) + Expectation(Y)

    To prevent this kind of automatic evaluation, use ``evaluate=False``:

    >>> Expectation(a*X + Y)
    a*Expectation(X) + Expectation(Y)
    >>> Expectation(a*X + Y, evaluate=False)
    Expectation(a*X + Y)
    """
    def __new__(cls, expr, condition=None, **kwargs):
        expr = _sympify(expr)
        if not kwargs.pop('evaluate', global_evaluate[0]):
            if condition is None:
                obj = Expr.__new__(cls, expr)
            else:
                condition = _sympify(condition)
                obj = Expr.__new__(cls, expr, condition)
            obj._condition = condition
            return obj

        if not expr.has(RandomSymbol):
            return expr

        if condition is not None:
            condition = _sympify(condition)

        if isinstance(expr, Add):
            return Add(*[Expectation(a, condition=condition) for a in expr.args])
        elif isinstance(expr, Mul):
            rv = []
            nonrv = []
            for a in expr.args:
                if isinstance(a, RandomSymbol) or a.has(RandomSymbol):
                    rv.append(a)
                else:
                    nonrv.append(a)
            return Mul(*nonrv)*Expectation(Mul(*rv), condition=condition, evaluate=False)
        else:
            if condition is None:
                obj = Expr.__new__(cls, expr)
            else:
                obj = Expr.__new__(cls, expr, condition)
            obj._condition = condition
            return obj

    def doit(self, **kwargs):
        return expectation(self.args[0], condition=self._condition, **kwargs)


class Variance(Expr):
    """
    Symbolic expression for the variance.

    Examples
    ========

    >>> from sympy import symbols
    >>> from sympy.stats import Normal, Expectation, Variance
    >>> mu = symbols("mu", positive=True)
    >>> sigma = symbols("sigma", positive=True)
    >>> X = Normal("X", mu, sigma)
    >>> Variance(X)
    Variance(X)
    >>> Variance(X).doit()
    sigma**2

    Rewrite the variance in terms of the expectation

    >>> Variance(X).rewrite(Expectation)
    -Expectation(X)**2 + Expectation(X**2)

    Some transformations based on the properties of the variance may happen:

    >>> from sympy.abc import a
    >>> Y = Normal("Y", 0, 1)
    >>> Variance(a*X)
    a**2*Variance(X)

    To prevent automatic evaluation, use ``evaluate=False``:

    >>> Variance(X + Y)
    2*Covariance(X, Y) + Variance(X) + Variance(Y)
    >>> Variance(X + Y, evaluate=False)
    Variance(X + Y)

    """
    def __new__(cls, arg, condition=None, **kwargs):
        arg = _sympify(arg)
        if not kwargs.pop('evaluate', global_evaluate[0]):
            if condition is None:
                obj = Expr.__new__(cls, arg)
            else:
                condition = _sympify(condition)
                obj = Expr.__new__(cls, arg, condition)
            obj._condition = condition
            return obj

        if not arg.has(RandomSymbol):
            return S.Zero

        if condition is not None:
            condition = _sympify(condition)

        if isinstance(arg, RandomSymbol):
            if condition is None:
                obj = Expr.__new__(cls, arg)
            else:
                obj = Expr.__new__(cls, arg, condition)
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
                obj = Expr.__new__(cls, Mul(*rv))
            else:
                obj = Expr.__new__(cls, Mul(*rv), condition)
            obj._condition = condition
            return Mul(*nonrv)*obj

        else:
            # this expression contains a RandomSymbol somehow:
            return Variance(arg, condition, evaluate=False)

    def _eval_rewrite_as_Expectation(self, arg, condition=None):
            e1 = Expectation(arg**2, condition)
            e2 = Expectation(arg, condition)**2
            return e1 - e2

    def doit(self, **kwargs):
        return variance(self.args[0], self._condition, **kwargs)


class Covariance(Expr):
    """
    Symbolic expression for the covariance.

    Examples
    ========

    >>> from sympy.stats import Covariance
    >>> from sympy.stats import Normal
    >>> X = Normal("X", 3, 2)
    >>> Y = Normal("Y", 0, 1)
    >>> Z = Normal("Z", 0, 1)
    >>> W = Normal("W", 0, 1)
    >>> cexpr = Covariance(X, Y)
    >>> cexpr
    Covariance(X, Y)

    Evaluate the covariance, `X` and `Y` are independent,
    therefore zero is the result:

    >>> cexpr.doit()
    0

    Rewrite the covariance expression in terms of expectations:

    >>> from sympy.stats import Expectation
    >>> cexpr.rewrite(Expectation)
    Expectation(X*Y) - Expectation(X)*Expectation(Y)

    This class is aware of some properties of the covariance:

    >>> Covariance(X, X)
    Variance(X)
    >>> from sympy.abc import a, b, c, d
    >>> Covariance(a*X, b*Y)
    a*b*Covariance(X, Y)

    In order to prevent such transformations, use ``evaluate=False``:

    >>> Covariance(a*X + b*Y, c*Z + d*W)
    a*c*Covariance(X, Z) + a*d*Covariance(W, X) + b*c*Covariance(Y, Z) + b*d*Covariance(W, Y)
    >>> Covariance(a*X + b*Y, c*Z + d*W, evaluate=False)
    Covariance(a*X + b*Y, c*Z + d*W)
    """
    def __new__(cls, arg1, arg2, condition=None, **kwargs):
        arg1 = _sympify(arg1)
        arg2 = _sympify(arg2)

        if condition is not None:
            condition = _sympify(condition)

        if not kwargs.pop('evaluate', global_evaluate[0]):
            return cls._call_super_constructor(arg1, arg2, condition)

        if arg1 == arg2:
            return Variance(arg1, condition)

        if not arg1.has(RandomSymbol):
            return S.Zero
        if not arg2.has(RandomSymbol):
            return S.Zero

        arg1, arg2 = sorted([arg1, arg2], key=default_sort_key)

        if isinstance(arg1, RandomSymbol) and isinstance(arg2, RandomSymbol):
            return cls._call_super_constructor(arg1, arg2, condition)

        coeff_rv_list1 = cls._expand_single_argument(arg1.expand())
        coeff_rv_list2 = cls._expand_single_argument(arg2.expand())

        addends = [a*b*Covariance(*sorted([r1, r2], key=default_sort_key), condition=condition, evaluate=False)
                   for (a, r1) in coeff_rv_list1 for (b, r2) in coeff_rv_list2]
        return Add(*addends)

    @classmethod
    def _call_super_constructor(cls, arg1, arg2, condition):
        if condition is not None:
            obj = Expr.__new__(cls, arg1, arg2, condition)
        else:
            obj = Expr.__new__(cls, arg1, arg2)
        obj._condition = condition
        return obj

    @classmethod
    def _expand_single_argument(cls, expr):
        # return (coefficient, random_symbol) pairs:
        if isinstance(expr, RandomSymbol):
            return [(S.One, expr)]
        elif isinstance(expr, Add):
            outval = []
            for a in expr.args:
                if isinstance(a, Mul):
                    outval.append(cls._get_mul_nonrv_rv_tuple(a))
                elif isinstance(a, RandomSymbol):
                    outval.append((S.One, a))

            return outval
        elif isinstance(expr, Mul):
            return [cls._get_mul_nonrv_rv_tuple(expr)]
        elif expr.has(RandomSymbol):
            return [(S.One, expr)]

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

    def _eval_rewrite_as_Expectation(self, arg1, arg2, condition=None):
        e1 = Expectation(arg1*arg2, condition)
        e2 = Expectation(arg1, condition)*Expectation(arg2, condition)
        return e1 - e2

    def doit(self, **kwargs):
        return covariance(self.args[0], self.args[1], self._condition, **kwargs)
