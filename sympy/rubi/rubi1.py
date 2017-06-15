from .patterns1 import rubi_object
from .operation import Int
import matchpy
from .sympy2matchpy import sympy2matchpy
from sympy.core.sympify import sympify

rubi = rubi_object()

def rubi_integrate(expr, var):
    '''
    Main function for Rubi integeration.

    This function uses `eval`.
    '''

    if not isinstance(expr, matchpy.Expression):
        expr = Int(sympy2matchpy(expr), sympy2matchpy(var))

    result = rubi.replace(expr)

    if result == expr:
        print(('Unable to integrate: {}').format(expr))
        return None

    result = sympify(str(result))
    return result
