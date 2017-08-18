from sympy.external import import_module
matchpy = import_module("matchpy")
from sympy.utilities.decorator import doctest_depends_on

if matchpy:
    import inspect, re
    from sympy.integrals import Integral
    from sympy.integrals.rubi.patterns import rubi_object
    from sympy import S
    rubi = rubi_object()
    matcher = rubi.matcher

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

@doctest_depends_on(modules=('matchpy',))
def get_matching_rule_definition(expr, var):
    miter = matcher.match(Integral(expr, var))
    for fun, e in miter:
        print("Rule matching: ")
        print(inspect.getsourcefile(fun))
        code, lineno = inspect.getsourcelines(fun)
        print("On line: ", lineno)
        print("\n".join(code))
        print("Pattern matching: ")
        pattno = int(re.match(r"^\s*rule(\d+)", code[0]).group(1))
        print(matcher.patterns[pattno-1])
        print(e)
        print()
