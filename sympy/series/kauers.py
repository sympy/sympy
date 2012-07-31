from sympy import expand
"""
Takes as input the expression and the variable used in constructing the expression and returns the diffrence between function's value when the input is incremented to 1 and the original function value

More specifically, it takes as input the function f(x) expression and x used in function expression and returns f(x + 1) - f(x) 

If you want an increment other than 1, supply it as a third argument to the Delta

Examples
=========

>>> Delta(x**2, x) 
2*x + 1
>>> Delta(y**3 + 2*y**2 + 3*y + 4, y) 
3*y**2 + 7*y + 6
>>> Delta(x**2 + 3*x + 8, x, 2) 
4*x + 10
>>>Delta(z**3 + 8*z, z, 3)
9*z**2 + 27*z + 51

"""

def Delta(expr, n, d = 1):
    expr =(expr).expand()
    expr2 = expr.subs(n, n + d)
    expr2 = (expr2).expand()
    expr = expr2 - expr
    expr = expr.expand()
    return expr




