"""
This module implements a method to find
Euler-Lagrange Equations for given Lagrangian.
"""
from __future__ import annotations
from collections.abc import Iterable
from itertools import combinations_with_replacement
from typing import TYPE_CHECKING

from sympy.core.function import (Derivative, Function, diff)
from sympy.core.relational import Eq
from sympy.core.singleton import S
from sympy.core.symbol import Symbol
from sympy.core.sympify import sympify

if TYPE_CHECKING:
    from sympy.core.basic import Basic
    from sympy.core.expr import Expr


def euler_equations(
    L: Expr,
    funcs: Function | Iterable[Function] = (),
    vars: Symbol | Iterable[Symbol] = (),
) -> list[Eq]:
    r"""
    Find the Euler-Lagrange equations [1]_ for a given Lagrangian.

    Parameters
    ==========

    L : Expr
        The Lagrangian that should be a function of the functions listed
        in the second argument and their derivatives.

        For example, in the case of two functions $f(x,y)$, $g(x,y)$ and
        two independent variables $x$, $y$ the Lagrangian has the form:

            .. math:: L\left(f(x,y),g(x,y),\frac{\partial f(x,y)}{\partial x},
                      \frac{\partial f(x,y)}{\partial y},
                      \frac{\partial g(x,y)}{\partial x},
                      \frac{\partial g(x,y)}{\partial y},x,y\right)

        In many cases it is not necessary to provide anything, except the
        Lagrangian, it will be auto-detected (and an error raised if this
        cannot be done).

    funcs : Function or an iterable of Functions
        The functions that the Lagrangian depends on. The Euler equations
        are differential equations for each of these functions.

    vars : Symbol or an iterable of Symbols
        The Symbols that are the independent variables of the functions.

    Returns
    =======

    eqns : list of Eq
        The list of differential equations, one for each function.

    Examples
    ========

    >>> from sympy import euler_equations, Symbol, Function
    >>> x = Function('x')
    >>> t = Symbol('t')
    >>> L = (x(t).diff(t))**2/2 - x(t)**2/2
    >>> euler_equations(L, x(t), t)
    [Eq(-x(t) - Derivative(x(t), (t, 2)), 0)]
    >>> u = Function('u')
    >>> x = Symbol('x')
    >>> L = (u(t, x).diff(t))**2/2 - (u(t, x).diff(x))**2/2
    >>> euler_equations(L, u(t, x), [t, x])
    [Eq(-Derivative(u(t, x), (t, 2)) + Derivative(u(t, x), (x, 2)), 0)]

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Euler%E2%80%93Lagrange_equation

    """

    funcs_tuple: tuple[Function, ...]
    if isinstance(funcs, Iterable):
        funcs_tuple = tuple(funcs)
    else:
        funcs_tuple = (funcs,)

    if not funcs_tuple:
        funcs_tuple = tuple(L.atoms(Function))
    else:
        for f in funcs_tuple:
            if not isinstance(f, Function):
                raise TypeError('Function expected, got: %s' % f)

    vars_tuple: tuple[Basic, ...]
    if isinstance(vars, Iterable):
        vars_tuple = tuple(vars)
    else:
        vars_tuple = (vars,)

    if not vars_tuple:
        vars_tuple = funcs_tuple[0].args
    else:
        vars_tuple = tuple(sympify(var) for var in vars_tuple)

    if not all(isinstance(v, Symbol) for v in vars_tuple):
        raise TypeError('Variables are not symbols, got %s' % vars_tuple)

    for f in funcs_tuple:
        if not vars_tuple == f.args:
            raise ValueError("Variables %s do not match args: %s" % (vars_tuple, f))

    order = max([len(d.variables) for d in L.atoms(Derivative)
                        if d.expr in funcs_tuple] + [0])

    eqns: list[Eq] = []
    for f in funcs_tuple:
        eq = diff(L, f)
        for i in range(1, order + 1):
            for p in combinations_with_replacement(vars_tuple, i):
                eq = eq + S.NegativeOne**i*diff(L, diff(f, *p), *p)
        new_eq = Eq(eq, 0)
        if isinstance(new_eq, Eq):
            eqns.append(new_eq)

    return eqns
