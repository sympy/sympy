from sympy import expand
from sympy import diff

def finite_diff(expression, variable, increment = 1):
    """
    Takes as input the expression and the variable used in constructing the expression
    and returns the diffrence between function's value when the input is incremented
    to 1 and the original function value. More specifically, it takes as input the
    function f(x) expression and x used in function expression as variable and returns
    f(x + 1)-f(x).If you want an increment other than 1, supply it as a
    third argument to the finite_diff

    Examples
    =========
    >>> from sympy.abc import x, y, z
    >>> from sympy.series.kauers import finite_diff
    >>> finite_diff(x**2, x)
    2*x + 1
    >>> finite_diff(y**3 + 2*y**2 + 3*y + 4, y)
    3*y**2 + 7*y + 6
    >>> finite_diff(x**2 + 3*x + 8, x, 2)
    4*x + 10
    >>> finite_diff(z**3 + 8*z, z, 3)
    9*z**2 + 27*z + 51
    """
    expression = expression.expand()
    expression2 = expression.subs(variable, variable + increment)
    expression2 = expression2.expand()
    return expression2 - expression

def dominant_term(expr, n):
     """
     Takes as input a polynomial expression along with it's variable and
     returns the dominant terms of the expression

     Examples
     =========
     >>> from sympy.abc import n, x, z
     >>> from sympy.series.kauers import dominant_term
     >>> dominant_term(3*n**2 + 4*n + 8, n)
     n**2
     >>> dominant_term(4*x**5 + 6*x**3 + 2*x**2 + 5, x)
     x**5
     >>> dominant_term(z**3 + 3*z**2 +4*z + 8, z)
     z**3
     """
     counter = 1
     while (expr != 0):
         expr = diff(expr, n)
         counter = counter + 1
     return n**(counter - 2)
