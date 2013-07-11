from sympy import Function, sympify, diff, Eq, S, Symbol
from sympy.utilities import numbered_symbols

def euler_equations(L, funcs, vars):
    """Find the Euler-Lagrange equations for a given Lagrangian.

    For the case of two functions (f(x,y), g(x,y)), two independent variables
    (x,y) the Lagrangian would have the form::

        L(f(x,y),g(x,y),df/dx,df/dy,dg/dx,dg/dy,x,y)

    Parameters
    ==========
    L : Expr
        The Lagrangian that should be a function of the functions listed
        in the second argument and their first derivatives.
    funcs : Function or list/tuple of Functions
        The functions that the Lagrangian depends on. The Euler equations
        are differential equations for each of these functions.
    vars : Symbol or list/tuple of Symbols
        The Symbols that are the independent variables of the functions.

    Returns
    =======
    eqns : set of Eq
        The set of differential equations, one for each function.

    Examples
    ========

        >>> from sympy import Symbol, Function
    """

    if not isinstance(funcs, (tuple, list)):
        funcs = (funcs,)
    if not isinstance(vars, (tuple, list)):
        vars = (vars,)

    vars = [sympify(var) for var in vars]
    constants = numbered_symbols(prefix='C', cls=Symbol, start=1)

    eqns = []
    for f in funcs:
        if not isinstance(f, Function):
            raise TypeError('Function expected, got: %r' % f)
        if not set(vars) == f.free_symbols:
            raise ValueError("Variables %r don't match function arguments: %r" % (vars, f))
        eq = diff(L, f)
        if eq == S.Zero and len(vars) == 1:
            eqns.append(Eq(diff(L, diff(f, var)), constants.next()))
        else:
            for var in vars:
                eq = eq - diff(L, diff(f, var), var)
            eqns.append(Eq(eq, 0))
    return set(eqns)
