from sympy.external import import_module
matchpy = import_module("matchpy")
from sympy.core import Integer, Float
import inspect, re
from sympy import powsimp
from sympy.utilities.decorator import doctest_depends_on

if matchpy:
    from sympy.integrals.rubi.matchpy_setup.operations import operation_register
    from sympy.integrals.rubi.matchpy_setup.rubi_object import get_rubi_object
    from sympy.integrals.rubi.utility_function import (rubi_unevaluated_expr, process_trig)
    from sympy import Basic, S, Function, E


    operation_register()
    _E = rubi_unevaluated_expr(E)
    Integrate = Function('Integrate')
    rubi, rules_applied, rules_applied = get_rubi_object()

def _has_cycle():
    if rules_applied.count(rules_applied[-1]) == 1:
        return False
    if rules_applied[-1] == rules_applied[-2] == rules_applied[-3] == rules_applied[-4] == rules_applied[-5]:
        return True

def process_final_integral(expr):
    '''
    When there is recursion for more than 10 rules or in total 20 rules have been applied
    rubi returns `Integrate` in order to stop any further matching. After complete integration,
    Integrate needs to be replaced back to Integral. Also rubi's `exp` need to be replaced back
    to sympy's general `exp`.

    Examples
    ========
    >>> from sympy import Function, E
    >>> from sympy.integrals.rubi.rubi import process_final_integral
    >>> from sympy.integrals.rubi.utility_function import rubi_unevaluated_expr
    >>> Integrate = Function("Integrate")
    >>> from sympy.abc import a, x
    >>> _E = rubi_unevaluated_expr(E)
    >>> process_final_integral(Integrate(a, x))
    Integral(a, x)
    >>> process_final_integral(_E**5)
    exp(5)

    '''
    if expr.has(Integrate):
        expr = expr.replace(Integrate, Integral)
    if expr.has(_E):
        expr = expr.replace(_E, E)
    return expr

def rubi_powsimp(expr):
    '''
    This function is needed to preprocess an expression as done in matchpy
    `x^a*x^b` in matchpy auotmatically transforms to `x^(a+b)`

    Examples
    ========

    >>> from sympy.integrals.rubi.rubi import rubi_powsimp
    >>> from sympy.abc import a, b, x
    >>> rubi_powsimp(x**a*x**b)
    x**(a+b)

    '''
    lst_pow =[]
    lst_non_pow = []
    if isinstance(expr, Mul):
        for i in expr.args:
            if isinstance(i, (Pow, exp, sym_exp)):
                lst_pow.append(i)
            else:
                lst_non_pow.append(i)
        return powsimp(Mul(*lst_pow))*Mul(*lst_non_pow)
    return expr

@doctest_depends_on(modules=('matchpy',))
def rubi_integrate(expr, var, showsteps=False):
    '''
    Rule based algorithm for integration. Integrates the expression by applying
    transformation rules to the expression.

    Returns `Integrate` if an expression cannot be integrated.

    Parameters
    ==========
    expr : integrand expression
    var : variable of integration

    Returns Integral object if unable to integrate.
    '''
    expr = expr.replace(sym_exp, exp)
    rules_applied[:] = []
    expr = process_trig(expr)
    expr = rubi_powsimp(expr)
    if isinstance(expr, (int, Integer)) or isinstance(expr, (float, Float)):
        return S(expr)*var
    if isinstance(expr, Add):
        results = 0
        for ex in expr.args:
            rules_applied[:] = []
            results += rubi.replace(Integral(ex, var))
            rules_applied[:] = []
        return process_final_integral(results)

    results = rubi.replace(Integral(expr, var), max_count = 10)
    return process_final_integral(results)

@doctest_depends_on(modules=('matchpy',))
def util_rubi_integrate(expr, var, showsteps=False):
    expr = process_trig(expr)
    expr = expr.replace(sym_exp, exp)
    if isinstance(expr, (int, Integer)) or isinstance(expr, (float, Float)):
        return S(expr)*var
    if isinstance(expr, Add):
        return rubi_integrate(expr, var)
    if len(rules_applied) > 10:
        if _has_cycle() or len(rules_applied) > 20:
            return Integrate(expr, var)
    results = rubi.replace(Integral(expr, var), max_count = 10)
    rules_applied[:] = []
    return results

@doctest_depends_on(modules=('matchpy',))
def get_matching_rule_definition(expr, var):
    '''
    Prints the list or rules which match to `expr`.

    Parameters
    ==========
    expr : integrand expression
    var : variable of integration
    '''
    matcher = rubi.matcher
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
