from patterns import rubi_object
import matchpy
from sympy2matchpy import sympy2matchpy

def rubi_integrate(expr):
    '''
    Main function for Rubi integeration.

    This function uses `eval`
    '''
    if not isinstance(expr, matchpy.Expression):
        expr = sympy2matchpy(expr)

    rubi = rubi_object()
    res = rubi.replace(expr)

    return res
