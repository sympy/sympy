from .patterns import rubi_object
from .operation import Int
import matchpy
from .sympy2matchpy import sympy2matchpy


def rubi_integrate(expr, var):
    '''
    Main function for Rubi integeration.

    This function uses `eval`
    '''

    if not isinstance(expr, matchpy.Expression):
        #expr = Int(sympy2matchpy(expr), var)
        expr = Int(sympy2matchpy(expr), sympy2matchpy(var))

    rubi = rubi_object()
    result = rubi.replace(expr)

    return result
