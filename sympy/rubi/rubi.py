from .patterns import rubi_object
from .operation import Int
import matchpy
from .sympy2matchpy import sympy2matchpy
from sympy.core.sympify import sympify
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core import S
from sympy.rubi.matchpy2sympy import matchpy2sympy

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
        args = [rubi_integrate(i, var) for i in expr.args]
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
            return c*rubi_integrate(e, var)

    if not isinstance(expr, matchpy.Expression):
        expr = Int(sympy2matchpy(expr), sympy2matchpy(var))

    if showsteps:
        matches = rubi.matcher.match(expr)
        for _ in matches._match(rubi.matcher.root):
            for pattern_index in matches.patterns:
                renaming = rubi.matcher.pattern_vars[pattern_index]
                substitution = matches.substitution.rename({renamed: original for original, renamed in renaming.items()})
                pattern = rubi.matcher.patterns[pattern_index][0]
                print('{} matched with {}'.format(pattern , substitution))
                print()

    result = rubi.replace(expr)

    if result == expr:
        raise NotImplementedError('Unable to intergate {}'.format(expr))

    return matchpy2sympy(result)
