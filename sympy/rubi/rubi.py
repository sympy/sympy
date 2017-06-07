from .patterns import rubi_object
from .operation import Int
import matchpy
from .sympy2matchpy import sympy2matchpy
from sympy.core.sympify import sympify

def rubi_integrate(expr, var):
    '''
    Main function for Rubi integeration.

    This function uses `eval`.
    '''

    if not isinstance(expr, matchpy.Expression):
        expr = Int(sympy2matchpy(expr), sympy2matchpy(var))

    rubi = rubi_object()
    result = rubi.replace(expr)

    result = sympify(str(result))

    return result
