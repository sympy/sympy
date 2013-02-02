from sympy.core import Mul, Add, Pow, sympify
from sympy.core.symbol import symbols, Symbol
from sympy import expand


def _give_variables(args, variable):
    if variable not in args:
        variable = args[0]
    args_list = list(args)
    args_list.remove(variable)
    return [variable, Add(*args_list)]


def binomial_expand(expr, variable = None):
    """
    Returns a valid binomial expansion of an expression, in terms of
    'variable', if there is any.
    Otherwise, returns the function itself.
    Output depends on the value of exponent.

    Examples
    ========

    >>> from sympy import binomial_expand, symbols
    >>> x,y,z = symbols('x,y,z')
    >>> binomial_expand((x + y + z) ** 2)
        x**2 + 2*x*y + 2*x*z + y**2 + 2*y*z + z**2
    >>> binomial_expand((1 + x) ** 0.3)
        1.0 + 0.3*x - 0.105*x**2 + 0.0595*x**3 - 0.0401625*x**4 +
        0.02972025*x**5 + O(x**6)

    """
    expr = sympify(expr).evalf()
    if expr.is_Pow:
        if expr.base.is_Add:
            if type(expr.exp) == int:
                return expand(expr)
            elif type(expr.exp) == float:
                if variable not in args:
                    variable = args[0]
                args_list = list(args)
                args_list.remove(variable)
                variable2 = Add(*args_list)
                xin, yin = symbols('xin,yin')
                return series((xin + yin) ** expr.exp, yin).subs(
                    xin, variable).subs(yin, variable2).expand().n(3).evalf()
            else:
                return expr
        else:
            return expr
    else:
        return expr
        
        
    
    
    
