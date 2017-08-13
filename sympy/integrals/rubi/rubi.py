from sympy.external import import_module
matchpy = import_module("matchpy")
from sympy.utilities.decorator import doctest_depends_on

if matchpy:
    from sympy.integrals import Integral
    from sympy.integrals.rubi.patterns import rubi_object
    from sympy import S
    rubi = rubi_object()

@doctest_depends_on(modules=('matchpy',))
def rubi_integrate(expr, var, showsteps=False):
    '''
    Function for Rubi integeration.
    Returns Integral object if unable to integrate.
    '''
    if isinstance(expr, int) or isinstance(expr, float):
        return S(expr)*var

    result = rubi.replace(Integral(expr, var))

    return result
