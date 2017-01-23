"""
This module implements a method to find
Euler-Lagrange equations for given Lagrangian.
"""
from itertools import combinations_with_replacement
from sympy import Function, sympify, diff, Eq, S, Symbol, Derivative
from sympy.core.compatibility import iterable, range


def euler_equations(lagrangian, functions=(), variables=()):
    r"""
    Find the Euler-Lagrange equations [1]_ for a given Lagrangian.

    Parameters
    ==========

    lagrangian : Expr
        The Lagrangian that should be a function of the functions listed
        in the second argument and their derivatives.

        For example, in the case of two functions `f(x,y)`, `g(x,y)` and two
        independent variables `x`, `y` the Lagrangian L would have the form:

            .. math:: L\left(f(x,y),g(x,y),\frac{\partial f(x,y)}{\partial x},
                      \frac{\partial f(x,y)}{\partial y},
                      \frac{\partial g(x,y)}{\partial x},
                      \frac{\partial g(x,y)}{\partial y},x,y\right)

        In many cases it is not necessary to provide anything, except the
        Lagrangian, it will be auto-detected (and an error raised if this
        couldn't be done).

    functions : Function or an iterable of Functions
        The functions that the Lagrangian depends on. The Euler equations
        are differential equations for each of these functions.

    variables : Symbol or an iterable of Symbols
        The Symbols that are the independent variables of the functions.

    Returns
    =======

    equations : list of Eq
        The list of differential equations, one for each function.

    Examples
    ========

    >>> from sympy import Symbol, Function
    >>> from sympy.calculus.euler import euler_equations
    >>> x = Function('x')
    >>> t = Symbol('t')
    >>> L = (x(t).diff(t))**2/2 - x(t)**2/2
    >>> euler_equations(L, x(t), t)
    [Eq(-x(t) - Derivative(x(t), t, t), 0)]
    >>> u = Function('u')
    >>> x = Symbol('x')
    >>> L = (u(t, x).diff(t))**2/2 - (u(t, x).diff(x))**2/2
    >>> euler_equations(L, u(t, x), [t, x])
    [Eq(-Derivative(u(t, x), t, t) + Derivative(u(t, x), x, x), 0)]

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Euler%E2%80%93Lagrange_equation

    """

    functions = tuple(functions) if iterable(functions) else (functions,)

    if not functions:
        functions = tuple(lagrangian.atoms(Function))
    else:
        for function in functions:
            if not isinstance(function, Function):
                raise TypeError('Function expected, got: %s' % function)

    variables = tuple(variables) if iterable(variables) else (variables,)

    if not variables:
        variables = functions[0].args
    else:
        variables = tuple(sympify(var) for var in variables)

    if not all(isinstance(v, Symbol) for v in variables):
        raise TypeError('Variables are not symbols, got %s' % variables)

    for function in functions:
        if not variables == function.args:
            raise ValueError(
                "Variables %s don't match args: %s" % (variables, function)
            )

    order = max(len(d.variables) for d in lagrangian.atoms(Derivative)
                if d.expr in functions)

    equations = []
    for function in functions:
        equation = diff(lagrangian, function)
        for i in range(1, order + 1):
            for symbols in combinations_with_replacement(variables, i):
                equation += S.NegativeOne**i \
                    * diff(lagrangian, diff(function, *symbols), *symbols)
        equations.append(Eq(equation))

    return equations
