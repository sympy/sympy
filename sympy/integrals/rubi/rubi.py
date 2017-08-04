from sympy.external import import_module
matchpy = import_module("matchpy")
from sympy.utilities.decorator import doctest_depends_on

if matchpy:
    from sympy.core.add import Add
    from sympy.integrals import Integral
    from sympy.core import S

    from sympy.integrals.rubi.patterns import rubi_object
    rubi = rubi_object()

@doctest_depends_on(modules=('matchpy',))
def rubi_integrate(expr, var, showsteps=False):
    '''
    Main function for Rubi integeration.
    Returns Integral object if unable to integrate.
    '''
    if isinstance(expr, int):
        return S(expr)*var
    if expr.is_Add: # integrate each is_Add expression indivdually
        return Add(*[rubi_integrate(i, var) for i in expr.args])
    elif expr.is_Mul: # remove constants
        e, c = 1, 1
        for i in expr.args:
            if i.has(var):
                e = e*i
            else:
                c = c*i
        if c != 1:
            return c*rubi_integrate(e, var)

    result = rubi.replace(Integral(expr, var))

    return result
