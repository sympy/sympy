from sympy.external import import_module
matchpy = import_module("matchpy")

if matchpy is None:
    raise ImportError('MatchPy could not be imported')

from sympy.integrals.rubi.patterns import rubi_object
from sympy.integrals.rubi.operation import Int
from sympy.integrals.rubi.sympy2matchpy import sympy2matchpy
from sympy.core.sympify import sympify
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core import S
from sympy.integrals.rubi.matchpy2sympy import matchpy2sympy

rubi = rubi_object()

def rubi_integrate(expr, var, showsteps=False):
    '''
    Main function for Rubi integeration.
    Raises NotImplementedError Error if expression is not integrated
    '''
    if expr == None:
        return None
    if isinstance(expr, int) or isinstance(expr, float):
        return expr*var
    elif not expr.has(var): # handle constant expressions
        return expr*var
    elif isinstance(expr, Add): # integrate each is_Add expression indivdually
        return Add(*[rubi_integrate(i, var) for i in expr.args])
    elif isinstance(expr, Mul): #seperate out constants in Mul expression
        e = 1
        c = 1
        for i in expr.args:
            if i.has(var):
                e = e*i
            else:
                c = c*i
        if c != 1:
            return c*rubi_integrate(e, var)

    if not isinstance(expr, matchpy.Expression):
        expr = Int(sympy2matchpy(expr), sympy2matchpy(var))

    result = rubi.replace(expr)

    if result == expr:
        raise NotImplementedError('Unable to intergate {}'.format(expr))

    return matchpy2sympy(result)
