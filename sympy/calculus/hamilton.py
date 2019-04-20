"""
This module implements a method to find
Hamilton's Equations of motion for given Hamiltonian.
"""
from sympy import Function, sympify, diff, Eq, S, Symbol, Derivative
from sympy.core.compatibility import (iterable, range)


def hamilton_equations(H, coords, momenta, time):
    r"""
    Find the Hamilton's equations for a given Hamiltonian.

    Parameters
    ==========

    H : Expr
        The Hamiltonian that should be a function of the coordinates, momenta, and time.

    coords : Function or an iterable of Functions
        List of position coordinate that the Hamiltonian depends upon.

    momenta : Function or an iterable of Functions
        List of conjugate momentum coordinates that the Hamiltonian depends upon.

    time : Symbol
        The Symbol used for time, here the only free variable in the equation.

    Returns
    =======

    equations : list of Eq
        The list of differential equations, one for each function.

    Examples
    ========

    >>> from sympy import Symbol, Function
    >>> from sympy.calculus.hamilton import hamilton_equations
    >>> x = Function('x')
    >>> p = Function('p')
    >>> t = Symbol('t')
    >>> H = p(t)**2/2 + x(t)**2/2
    >>> hamilton_equations(H, x(t), p(t), t)
    [Eq(x(t) + Derivative(p(t), t), 0), Eq(p(t) - Derivative(x(t), t), 0)]

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Hamiltonian_mechanics#Deriving_Hamilton's_equations

    """

    coords = tuple(coords) if iterable(coords) else (coords,)
    momenta = tuple(momenta) if iterable(momenta) else (momenta,)

    for f in coords + momenta:
        if not isinstance(f, Function):
            raise TypeError('Function expected, got: %s' % f)

    if not isinstance(time, Symbol):
        raise TypeError('Symbol for Time is not a symbol, got %s' % time)

    for f in coords + momenta:
        if not tuple([time]) == f.args:
            raise ValueError("Arguments of Canonical Variable %s is not time" % f)

    equations = []
    funcs = list(zip(coords, momenta))
    for (coord, momentum) in funcs:
        eq = diff(H, coord) + diff(momentum, time)
        equations.append(Eq(eq, 0))
        eq = diff(H, momentum) - diff(coord, time)
        equations.append(Eq(eq, 0))

    return equations
