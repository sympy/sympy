from sympy import expand
from sympy import diff
from sympy import Sum

def finite_diff(expression, variable, increment=1):
    """
    Takes as input a polynomial expression and the variable used to construct it and returns the difference between function's value when the input is incremented to 1 and the original function value. The input expression can also be a polynomial summation. In that case it will return the difference between the summation with variable incremented by given increment value(1 if none is given) and the original function. More specifically, it takes as input a function or summation f(x) and x used in function expression as variable and returns f(x + increment) - f(x).If you want an increment other than 1, supply it as a third argument.

    Examples
    =========
    >>> from sympy.abc import x, y, z, k, n
    >>> from sympy.series.kauers import finite_diff
    >>> from  sympy import Sum
    >>> finite_diff(x**2, x)
    2*x + 1
    >>> finite_diff(y**3 + 2*y**2 + 3*y + 4, y)
    3*y**2 + 7*y + 6
    >>> finite_diff(x**2 + 3*x + 8, x, 2)
    4*x + 10
    >>> finite_diff(z**3 + 8*z, z, 3)
    9*z**2 + 27*z + 51
    >>> finite_diff(Sum(1/k, (k, 1, n)), k)
    1/(n + 1)
    >>> finite_diff(Sum(k, (k, 1, n)), k)
    n + 1

    """
    if isinstance(expression, Sum):
        function = expression.function
        limit = expression.limits
        limit_tuple = limit[0]
        var = limit_tuple[-1]
        return function.subs(variable, var + increment)
    else:
        expression = expression.expand()
        expression2 = expression.subs(variable, variable + increment)
        expression2 = expression2.expand()
        return expression2 - expression
