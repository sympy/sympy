from .patterns import rubi_object
from .operation import Int
import matchpy
from .sympy2matchpy import sympy2matchpy
from sympy.core.sympify import sympify
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core import S
from .matchpy2sympy import matchpy2sympy

rubi = rubi_object()

def rubi_integrate(expr, var):
    '''
    Main function for Rubi integeration.
    returns 0 if expr is not integrated
    '''
    if expr == None:
        return S(1)
    if isinstance(expr, int) or isinstance(expr, float):
        return expr*var
    elif not expr.has(var): # handle constant expressions
        return expr*var
    elif isinstance(expr, Add): # integrate each is_Add expression indivdually
        args = [rubi_integrate(i, var) for i in expr.args]
        if None in args:
            return S(1)
        return Add(*args)
    elif isinstance(expr, Mul): #seperate out constants in Mul expression
        e = 1
        c = 1
        for i in expr.args:
            if i.has(var):
                e = e*i
            else:
                c = c*i
        if c != 1:
            res = rubi_integrate(e, var)
            if res != None:
                return c*res
            else:
                return S(1)

    if not isinstance(expr, matchpy.Expression):
        expr = Int(sympy2matchpy(expr), sympy2matchpy(var))

    result = rubi.replace(expr)

    if result == expr:
        #print(('Unable to integrate: {}').format(expr))
        return S(1)
    #print('result: ', result)
    return matchpy2sympy(result)
