"""
This module contains pdsolve() and different helper functions that it
uses. It is heavily inspired by the ode module and hence the basic
infrastructure remains the same.

**Functions in this module**

    These are the user functions in this module:

    - pdsolve()     - Solves PDE's
    - classify_pde() - Classifies PDEs into possible hints for dsolve().
    - pde_separate() - Separate variables in partial differential equation either by
                       additive or multiplicative separation approach.

    These are the helper functions in this module:

    - pde_separate_add() - Helper function for searching additive separable solutions.
    - pde_separate_mul() - Helper function for searching multiplicative
                           separable solutions.

**Currently implemented solver methods**

The following methods are implemented for solving partial differential
equations.  See the docstrings of the various pde_hint() functions for
more information on each (run help(pde)):

  - 1st order linear homogeneous partial differential equations
    with constant coefficients.
  - 1st order linear general partial differential equations
    with constant coefficients.
  - 1st order linear partial differential equations with
    variable coefficients.

"""
from __future__ import annotations
from functools import reduce
from abc import ABC, abstractmethod
from dataclasses import dataclass, field, replace

from itertools import combinations_with_replacement
from sympy.simplify import simplify  # type: ignore
from sympy.core import Add, S
from sympy.core.function import Function, expand, AppliedUndef, Subs, Derivative
from sympy.core.relational import Equality, Eq, Relational
from sympy.core.symbol import Symbol, Wild, Dummy, symbols
from sympy.core.numbers import pi, I
from sympy.core.power import Pow
from sympy.functions import exp
from sympy.integrals.integrals import Integral, integrate
from sympy.utilities.iterables import has_dups, is_sequence, sift
from sympy.utilities.misc import filldedent
from sympy.concrete.summations import Sum
from sympy.core.exprtools import factor_terms
from sympy.core.mul import Mul
from sympy.functions.elementary.complexes import conjugate, im, re
from sympy.functions.elementary.miscellaneous import Max, Min, sqrt
from sympy.functions.elementary.piecewise import Piecewise
from sympy.functions.elementary.trigonometric import cos, sin
from sympy.functions.elementary.hyperbolic import cosh, sinh
from sympy.sets.conditionset import ConditionSet
from sympy.sets.fancysets import ImageSet
from sympy.sets.sets import FiniteSet, Interval, Union, Intersection
from sympy.series.limits import Limit
from sympy.solvers.solveset import linear_eq_to_matrix, solveset
from sympy.solvers.diophantine.diophantine import diophantine
from sympy.logic.boolalg import And
from sympy.polys.polytools import Poly, cancel, factor
from sympy.polys.polyerrors import PolynomialError
from sympy.polys.domains import ZZ, QQ

from sympy.solvers.deutils import _preprocess, ode_order, _desolve
from sympy.solvers.solvers import solve
from sympy.simplify.radsimp import collect

import operator

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from sympy.core.expr import Expr

allhints = (
    "1st_linear_constant_coeff_homogeneous",
    "1st_linear_constant_coeff",
    "1st_linear_constant_coeff_Integral",
    "1st_linear_variable_coeff",
    "2nd_linear_constant_coeff_separable",
    )


def pdsolve(eq, func=None, hint='default', dict=False, solvefun=None, ics=None, **kwargs):
    """
    Solves any (supported) kind of partial differential equation.

    **Usage**

        pdsolve(eq, f(x,y), hint) -> Solve partial differential equation
        eq for function f(x,y), using method hint.

    **Details**

        ``eq`` can be any supported partial differential equation (see
            the pde docstring for supported methods).  This can either
            be an Equality, or an expression, which is assumed to be
            equal to 0.

        ``f(x,y)`` is a function of two variables whose derivatives in that
            variable make up the partial differential equation. In many
            cases it is not necessary to provide this; it will be autodetected
            (and an error raised if it could not be detected).

        ``hint`` is the solving method that you want pdsolve to use.  Use
            classify_pde(eq, f(x,y)) to get all of the possible hints for
            a PDE.  The default hint, 'default', will use whatever hint
            is returned first by classify_pde().  See Hints below for
            more options that you can use for hint.

        ``solvefun`` is the convention used for arbitrary functions returned
            by the PDE solver. If not set by the user, it is set by default
            to be F.

    **Hints**

        Aside from the various solving methods, there are also some
        meta-hints that you can pass to pdsolve():

        "default":
                This uses whatever hint is returned first by
                classify_pde(). This is the default argument to
                pdsolve().

        "all":
                To make pdsolve apply all relevant classification hints,
                use pdsolve(PDE, func, hint="all").  This will return a
                dictionary of hint:solution terms.  If a hint causes
                pdsolve to raise the NotImplementedError, value of that
                hint's key will be the exception object raised.  The
                dictionary will also include some special keys:

                - order: The order of the PDE.  See also ode_order() in
                  deutils.py
                - default: The solution that would be returned by
                  default.  This is the one produced by the hint that
                  appears first in the tuple returned by classify_pde().

        "all_Integral":
                This is the same as "all", except if a hint also has a
                corresponding "_Integral" hint, it only returns the
                "_Integral" hint.  This is useful if "all" causes
                pdsolve() to hang because of a difficult or impossible
                integral.  This meta-hint will also be much faster than
                "all", because integrate() is an expensive routine.

        See also the classify_pde() docstring for more info on hints,
        and the pde docstring for a list of all supported hints.

    **Tips**
        - You can declare the derivative of an unknown function this way:

            >>> from sympy import Function, Derivative
            >>> from sympy.abc import x, y # x and y are the independent variables
            >>> f = Function("f")(x, y) # f is a function of x and y
            >>> # fx will be the partial derivative of f with respect to x
            >>> fx = Derivative(f, x)
            >>> # fy will be the partial derivative of f with respect to y
            >>> fy = Derivative(f, y)

        - See test_pde.py for many tests, which serves also as a set of
          examples for how to use pdsolve().
        - pdsolve always returns an Equality class (except for the case
          when the hint is "all" or "all_Integral"). Note that it is not possible
          to get an explicit solution for f(x, y) as in the case of ODE's
        - Do help(pde.pde_hintname) to get help more information on a
          specific hint


    Examples
    ========

    >>> from sympy.solvers.pde import pdsolve
    >>> from sympy import Function, Eq
    >>> from sympy.abc import x, y
    >>> f = Function('f')
    >>> u = f(x, y)
    >>> ux = u.diff(x)
    >>> uy = u.diff(y)
    >>> eq = Eq(1 + (2*(ux/u)) + (3*(uy/u)), 0)
    >>> pdsolve(eq)
    Eq(f(x, y), F(3*x - 2*y)*exp(-2*x/13 - 3*y/13))

    """

    if not solvefun:
        solvefun = Function('F')

    # See the docstring of _desolve for more details.
    hints = _desolve(eq, func=func, hint=hint, ics=ics, simplify=True,
                     type='pde', solvefun=solvefun, **kwargs)
    eq = hints.pop('eq', False)
    all_ = hints.pop('all', False)

    if all_:
        # TODO : 'best' hint should be implemented when adequate
        # number of hints are added.
        pdedict = {}
        failed_hints = {}
        gethints = classify_pde(eq, dict=True, ics=ics)
        pdedict.update({'order': gethints['order'],
                        'default': gethints['default']})
        for hint in hints:
            try:
                rv = _helper_simplify(eq, hint, hints[hint]['func'],
                    hints[hint]['order'], hints[hint][hint], solvefun)
            except NotImplementedError as detail:
                failed_hints[hint] = detail
            else:
                pdedict[hint] = rv
        pdedict.update(failed_hints)
        return pdedict

    else:
        return _helper_simplify(eq, hints['hint'], hints['func'],
                                hints['order'], hints[hints['hint']], solvefun)


def _helper_simplify(eq, hint, func, order, match, solvefun):
    """Helper function of pdsolve that calls the respective
    pde functions to solve for the partial differential
    equations. This minimizes the computation in
    calling _desolve multiple times.
    """
    solvefunc = globals()["pde_" + hint.removesuffix("_Integral")]
    return _handle_Integral(solvefunc(eq, func, order,
        match, solvefun), func, order, hint)


def _handle_Integral(expr, func, order, hint):
    r"""
    Converts a solution with integrals in it into an actual solution.

    Simplifies the integral mainly using doit()
    """
    if hint.endswith("_Integral"):
        return expr

    elif hint == "1st_linear_constant_coeff":
        return simplify(expr.doit())

    else:
        return expr


def classify_pde(eq, func=None, dict=False, ics=None, *, prep=True, **kwargs):
    """
    Returns a tuple of possible pdsolve() classifications for a PDE.

    The tuple is ordered so that first item is the classification that
    pdsolve() uses to solve the PDE by default.  In general,
    classifications near the beginning of the list will produce
    better solutions faster than those near the end, though there are
    always exceptions.  To make pdsolve use a different classification,
    use pdsolve(PDE, func, hint=<classification>).  See also the pdsolve()
    docstring for different meta-hints you can use.

    If ``dict`` is true, classify_pde() will return a dictionary of
    hint:match expression terms. This is intended for internal use by
    pdsolve().  Note that because dictionaries are ordered arbitrarily,
    this will most likely not be in the same order as the tuple.

    You can get help on different hints by doing help(pde.pde_hintname),
    where hintname is the name of the hint without "_Integral".

    See sympy.pde.allhints or the sympy.pde docstring for a list of all
    supported hints that can be returned from classify_pde.


    Examples
    ========

    >>> from sympy.solvers.pde import classify_pde
    >>> from sympy import Function, Eq
    >>> from sympy.abc import x, y
    >>> f = Function('f')
    >>> u = f(x, y)
    >>> ux = u.diff(x)
    >>> uy = u.diff(y)
    >>> eq = Eq(1 + (2*(ux/u)) + (3*(uy/u)), 0)
    >>> classify_pde(eq)
    ('1st_linear_constant_coeff_homogeneous',)
    """

    if func and len(func.args) != 2:
        raise NotImplementedError("Right now only partial "
            "differential equations of two variables are supported")

    if prep or func is None:
        prep, func_ = _preprocess(eq, func)
        if func is None:
            func = func_

    if isinstance(eq, Equality):
        if eq.rhs != 0:
            return classify_pde(eq.lhs - eq.rhs, func, dict=dict, ics=ics, **kwargs)
        eq = eq.lhs

    f = func.func
    x = func.args[0]
    y = func.args[1]
    fx = f(x,y).diff(x)
    fy = f(x,y).diff(y)

    # TODO : For now pde.py uses support offered by the ode_order function
    # to find the order with respect to a multi-variable function. An
    # improvement could be to classify the order of the PDE on the basis of
    # individual variables.
    order = ode_order(eq, f(x,y))

    # hint:matchdict or hint:(tuple of matchdicts)
    # Also will contain "default":<default hint> and "order":order items.
    matching_hints = {'order': order}

    if not order:
        if dict:
            matching_hints["default"] = None
            return matching_hints
        return ()

    eq = expand(eq)

    a = Wild('a', exclude = [f(x,y)])
    b = Wild('b', exclude = [f(x,y), fx, fy, x, y])
    c = Wild('c', exclude = [f(x,y), fx, fy, x, y])
    d = Wild('d', exclude = [f(x,y), fx, fy, x, y])
    e = Wild('e', exclude = [f(x,y), fx, fy])
    n = Wild('n', exclude = [x, y])
    # Try removing the smallest power of f(x,y)
    # from the highest partial derivatives of f(x,y)
    reduced_eq = eq
    if eq.is_Add:
        power = None
        for i in set(combinations_with_replacement((x,y), order)):
            coeff = eq.coeff(f(x,y).diff(*i))
            if coeff == 1:
                continue
            match = coeff.match(a*f(x,y)**n)
            if match and match[a]:
                if power is None or match[n] < power:
                    power = match[n]
        if power:
            den = f(x,y)**power
            reduced_eq = Add(*[arg/den for arg in eq.args])

    if order == 1:
        reduced_eq = collect(reduced_eq, f(x, y))
        r = reduced_eq.match(b*fx + c*fy + d*f(x,y) + e)
        if r:
            if not r[e]:
                ## Linear first-order homogeneous partial-differential
                ## equation with constant coefficients
                r.update({'b': b, 'c': c, 'd': d})
                matching_hints["1st_linear_constant_coeff_homogeneous"] = r
            elif r[b]**2 + r[c]**2 != 0:
                ## Linear first-order general partial-differential
                ## equation with constant coefficients
                r.update({'b': b, 'c': c, 'd': d, 'e': e})
                matching_hints["1st_linear_constant_coeff"] = r
                matching_hints["1st_linear_constant_coeff_Integral"] = r

        else:
            b = Wild('b', exclude=[f(x, y), fx, fy])
            c = Wild('c', exclude=[f(x, y), fx, fy])
            d = Wild('d', exclude=[f(x, y), fx, fy])
            r = reduced_eq.match(b*fx + c*fy + d*f(x,y) + e)
            if r:
                r.update({'b': b, 'c': c, 'd': d, 'e': e})
                matching_hints["1st_linear_variable_coeff"] = r

    elif order == 2 and ics is not None:
        solvefun = kwargs.get("solvefun") or Function("F")
        solver = _SeparationOfVariables(_BVP.from_input(reduced_eq, func, ics, solvefun))
        if solver.matches():
            matching_hints["2nd_linear_constant_coeff_separable"] = solver

    # Order keys based on allhints.
    rettuple = tuple(i for i in allhints if i in matching_hints)

    if dict:
        # Dictionaries are ordered arbitrarily, so make note of which
        # hint would come first for pdsolve().  Use an ordered dict in Py 3.
        matching_hints["default"] = None
        matching_hints["ordered_hints"] = rettuple
        for i in allhints:
            if i in matching_hints:
                matching_hints["default"] = i
                break
        return matching_hints
    return rettuple


def checkpdesol(pde, sol, func=None, solve_for_func=True):
    """
    Checks if the given solution satisfies the partial differential
    equation.

    pde is the partial differential equation which can be given in the
    form of an equation or an expression. sol is the solution for which
    the pde is to be checked. This can also be given in an equation or
    an expression form. If the function is not provided, the helper
    function _preprocess from deutils is used to identify the function.

    If a sequence of solutions is passed, the same sort of container will be
    used to return the result for each solution.

    The following methods are currently being implemented to check if the
    solution satisfies the PDE:

        1. Directly substitute the solution in the PDE and check. If the
           solution has not been solved for f, then it will solve for f
           provided solve_for_func has not been set to False.

    If the solution satisfies the PDE, then a tuple (True, 0) is returned.
    Otherwise a tuple (False, expr) where expr is the value obtained
    after substituting the solution in the PDE. However if a known solution
    returns False, it may be due to the inability of doit() to simplify it to zero.

    Examples
    ========

    >>> from sympy import Function, symbols
    >>> from sympy.solvers.pde import checkpdesol, pdsolve
    >>> x, y = symbols('x y')
    >>> f = Function('f')
    >>> eq = 2*f(x,y) + 3*f(x,y).diff(x) + 4*f(x,y).diff(y)
    >>> sol = pdsolve(eq)
    >>> assert checkpdesol(eq, sol)[0]
    >>> eq = x*f(x,y) + f(x,y).diff(x)
    >>> checkpdesol(eq, sol)
    (False, (x*F(4*x - 3*y) - 6*F(4*x - 3*y)/25 + 4*Subs(Derivative(F(_xi_1), _xi_1), _xi_1, 4*x - 3*y))*exp(-6*x/25 - 8*y/25))
    """

    # Converting the pde into an equation
    if not isinstance(pde, Equality):
        pde = Eq(pde, 0)

    # If no function is given, try finding the function present.
    if func is None:
        try:
            _, func = _preprocess(pde.lhs)
        except ValueError:
            funcs = [s.atoms(AppliedUndef) for s in (
                sol if is_sequence(sol, set) else [sol])]
            funcs = set().union(funcs)
            if len(funcs) != 1:
                raise ValueError(
                    'must pass func arg to checkpdesol for this case.')
            func = funcs.pop()

    # If the given solution is in the form of a list or a set
    # then return a list or set of tuples.
    if is_sequence(sol, set):
        return type(sol)([checkpdesol(
            pde, i, func=func,
            solve_for_func=solve_for_func) for i in sol])

    # Convert solution into an equation
    if not isinstance(sol, Equality):
        sol = Eq(func, sol)
    elif sol.rhs == func:
        sol = sol.reversed

    # Try solving for the function
    solved = sol.lhs == func and not sol.rhs.has(func)
    if solve_for_func and not solved:
        solved = solve(sol, func)
        if solved:
            if len(solved) == 1:
                return checkpdesol(pde, Eq(func, solved[0]),
                    func=func, solve_for_func=False)
            else:
                return checkpdesol(pde, [Eq(func, t) for t in solved],
                    func=func, solve_for_func=False)

    # try direct substitution of the solution into the PDE and simplify
    if sol.lhs == func:
        pde = pde.lhs - pde.rhs
        s = simplify(pde.subs(func, sol.rhs).doit())
        return s is S.Zero, s

    raise NotImplementedError(filldedent('''
        Unable to test if %s is a solution to %s.''' % (sol, pde)))



def pde_1st_linear_constant_coeff_homogeneous(eq, func, order, match, solvefun):
    r"""
    Solves a first order linear homogeneous
    partial differential equation with constant coefficients.

    The general form of this partial differential equation is

    .. math:: a \frac{\partial f(x,y)}{\partial x}
              + b \frac{\partial f(x,y)}{\partial y} + c f(x,y) = 0

    where `a`, `b` and `c` are constants.

    The general solution is of the form:

    .. math::
        f(x, y) = F(- a y + b x ) e^{- \frac{c (a x + b y)}{a^2 + b^2}}

    and can be found in SymPy with ``pdsolve``::

        >>> from sympy.solvers import pdsolve
        >>> from sympy.abc import x, y, a, b, c
        >>> from sympy import Function, pprint
        >>> f = Function('f')
        >>> u = f(x,y)
        >>> ux = u.diff(x)
        >>> uy = u.diff(y)
        >>> genform = a*ux + b*uy + c*u
        >>> pprint(genform)
          d               d
        a*--(f(x, y)) + b*--(f(x, y)) + c*f(x, y)
          dx              dy

        >>> pprint(pdsolve(genform))
                                 -c*(a*x + b*y)
                                 ---------------
                                      2    2
                                     a  + b
        f(x, y) = F(-a*y + b*x)*e

    Examples
    ========

    >>> from sympy import pdsolve
    >>> from sympy import Function, pprint
    >>> from sympy.abc import x,y
    >>> f = Function('f')
    >>> pdsolve(f(x,y) + f(x,y).diff(x) + f(x,y).diff(y))
    Eq(f(x, y), F(x - y)*exp(-x/2 - y/2))
    >>> pprint(pdsolve(f(x,y) + f(x,y).diff(x) + f(x,y).diff(y)))
                          x   y
                        - - - -
                          2   2
    f(x, y) = F(x - y)*e

    References
    ==========

    - Viktor Grigoryan, "Partial Differential Equations"
      Math 124A - Fall 2010, pp.7

    """
    # TODO : For now homogeneous first order linear PDE's having
    # two variables are implemented. Once there is support for
    # solving systems of ODE's, this can be extended to n variables.

    f = func.func
    x = func.args[0]
    y = func.args[1]
    b = match[match['b']]
    c = match[match['c']]
    d = match[match['d']]
    return Eq(f(x,y), exp(-S(d)/(b**2 + c**2)*(b*x + c*y))*solvefun(c*x - b*y))


def pde_1st_linear_constant_coeff(eq, func, order, match, solvefun):
    r"""
    Solves a first order linear partial differential equation
    with constant coefficients.

    The general form of this partial differential equation is

    .. math:: a \frac{\partial f(x,y)}{\partial x}
              + b \frac{\partial f(x,y)}{\partial y}
              + c f(x,y) = G(x,y)

    where `a`, `b` and `c` are constants and `G(x, y)` can be an arbitrary
    function in `x` and `y`.

    The general solution of the PDE is:

    .. math::
        f(x, y) = \left. \left[F(\eta) + \frac{1}{a^2 + b^2}
        \int\limits^{a x + b y} G\left(\frac{a \xi + b \eta}{a^2 + b^2},
        \frac{- a \eta + b \xi}{a^2 + b^2} \right)
        e^{\frac{c \xi}{a^2 + b^2}}\, d\xi\right]
        e^{- \frac{c \xi}{a^2 + b^2}}
        \right|_{\substack{\eta=- a y + b x\\ \xi=a x + b y }}\, ,

    where `F(\eta)` is an arbitrary single-valued function. The solution
    can be found in SymPy with ``pdsolve``::

        >>> from sympy.solvers import pdsolve
        >>> from sympy.abc import x, y, a, b, c
        >>> from sympy import Function, pprint
        >>> f = Function('f')
        >>> G = Function('G')
        >>> u = f(x, y)
        >>> ux = u.diff(x)
        >>> uy = u.diff(y)
        >>> genform = a*ux + b*uy + c*u - G(x,y)
        >>> pprint(genform)
          d               d
        a*--(f(x, y)) + b*--(f(x, y)) + c*f(x, y) - G(x, y)
          dx              dy
        >>> pprint(pdsolve(genform, hint='1st_linear_constant_coeff_Integral'))
                  //          a*x + b*y                                             \         \|
                  ||              /                                                 |         ||
                  ||             |                                                  |         ||
                  ||             |                                      c*xi        |         ||
                  ||             |                                     -------      |         ||
                  ||             |                                      2    2      |         ||
                  ||             |      /a*xi + b*eta  -a*eta + b*xi\  a  + b       |         ||
                  ||             |     G|------------, -------------|*e        d(xi)|         ||
                  ||             |      |   2    2         2    2   |               |         ||
                  ||             |      \  a  + b         a  + b    /               |  -c*xi  ||
                  ||             |                                                  |  -------||
                  ||            /                                                   |   2    2||
                  ||                                                                |  a  + b ||
        f(x, y) = ||F(eta) + -------------------------------------------------------|*e       ||
                  ||                                  2    2                        |         ||
                  \\                                 a  + b                         /         /|eta=-a*y + b*x, xi=a*x + b*y

    Examples
    ========

    >>> from sympy.solvers.pde import pdsolve
    >>> from sympy import Function, pprint, exp
    >>> from sympy.abc import x,y
    >>> f = Function('f')
    >>> eq = -2*f(x,y).diff(x) + 4*f(x,y).diff(y) + 5*f(x,y) - exp(x + 3*y)
    >>> pdsolve(eq)
    Eq(f(x, y), (F(4*x + 2*y)*exp(x/2) + exp(x + 4*y)/15)*exp(-y))

    References
    ==========

    - Viktor Grigoryan, "Partial Differential Equations"
      Math 124A - Fall 2010, pp.7

    """

    # TODO : For now homogeneous first order linear PDE's having
    # two variables are implemented. Once there is support for
    # solving systems of ODE's, this can be extended to n variables.
    xi, eta = symbols("xi eta")
    f = func.func
    x = func.args[0]
    y = func.args[1]
    b = match[match['b']]
    c = match[match['c']]
    d = match[match['d']]
    e = -match[match['e']]
    expterm = exp(-S(d)/(b**2 + c**2)*xi)
    functerm = solvefun(eta)
    solvedict = solve((b*x + c*y - xi, c*x - b*y - eta), x, y)
    # Integral should remain as it is in terms of xi,
    # doit() should be done in _handle_Integral.
    genterm = (1/S(b**2 + c**2))*Integral(
        (1/expterm*e).subs(solvedict), (xi, b*x + c*y))
    return Eq(f(x,y), Subs(expterm*(functerm + genterm),
        (eta, xi), (c*x - b*y, b*x + c*y)))


def pde_1st_linear_variable_coeff(eq, func, order, match, solvefun):
    r"""
    Solves a first order linear partial differential equation
    with variable coefficients. The general form of this partial
    differential equation is

    .. math:: a(x, y) \frac{\partial f(x, y)}{\partial x}
                + b(x, y) \frac{\partial f(x, y)}{\partial y}
                + c(x, y) f(x, y) = G(x, y)

    where `a(x, y)`, `b(x, y)`, `c(x, y)` and `G(x, y)` are arbitrary
    functions in `x` and `y`. This PDE is converted into an ODE by
    making the following transformation:

    1. `\xi` as `x`

    2. `\eta` as the constant in the solution to the differential
       equation `\frac{dy}{dx} = -\frac{b}{a}`

    Making the previous substitutions reduces it to the linear ODE

    .. math:: a(\xi, \eta)\frac{du}{d\xi} + c(\xi, \eta)u - G(\xi, \eta) = 0

    which can be solved using ``dsolve``.

    >>> from sympy.abc import x, y
    >>> from sympy import Function, pprint
    >>> a, b, c, G, f= [Function(i) for i in ['a', 'b', 'c', 'G', 'f']]
    >>> u = f(x,y)
    >>> ux = u.diff(x)
    >>> uy = u.diff(y)
    >>> genform = a(x, y)*u + b(x, y)*ux + c(x, y)*uy - G(x,y)
    >>> pprint(genform)
                                         d                     d
    -G(x, y) + a(x, y)*f(x, y) + b(x, y)*--(f(x, y)) + c(x, y)*--(f(x, y))
                                         dx                    dy


    Examples
    ========

    >>> from sympy.solvers.pde import pdsolve
    >>> from sympy import Function, pprint
    >>> from sympy.abc import x,y
    >>> f = Function('f')
    >>> eq =  x*(u.diff(x)) - y*(u.diff(y)) + y**2*u - y**2
    >>> pdsolve(eq)
    Eq(f(x, y), F(x*y)*exp(y**2/2) + 1)

    References
    ==========

    - Viktor Grigoryan, "Partial Differential Equations"
      Math 124A - Fall 2010, pp.7

    """
    from sympy.solvers.ode import dsolve

    eta = symbols("eta")
    f = func.func
    x = func.args[0]
    y = func.args[1]
    b = match[match['b']]
    c = match[match['c']]
    d = match[match['d']]
    e = -match[match['e']]


    if not d:
         # To deal with cases like b*ux = e or c*uy = e
        if not (b and c):
            if c:
                try:
                    tsol = integrate(e/c, y)
                except NotImplementedError:
                    raise NotImplementedError("Unable to find a solution"
                        " due to inability of integrate")
                else:
                    return Eq(f(x,y), solvefun(x) + tsol)
            if b:
                try:
                    tsol = integrate(e/b, x)
                except NotImplementedError:
                    raise NotImplementedError("Unable to find a solution"
                        " due to inability of integrate")
                else:
                    return Eq(f(x,y), solvefun(y) + tsol)

    if not c:
        # To deal with cases when c is 0, a simpler method is used.
        # The PDE reduces to b*(u.diff(x)) + d*u = e, which is a linear ODE in x
        plode = f(x).diff(x)*b + d*f(x) - e
        sol = dsolve(plode, f(x))
        syms = sol.free_symbols - plode.free_symbols - {x, y}
        rhs = _simplify_variable_coeff(sol.rhs, syms, solvefun, y)
        return Eq(f(x, y), rhs)

    if not b:
        # To deal with cases when b is 0, a simpler method is used.
        # The PDE reduces to c*(u.diff(y)) + d*u = e, which is a linear ODE in y
        plode = f(y).diff(y)*c + d*f(y) - e
        sol = dsolve(plode, f(y))
        syms = sol.free_symbols - plode.free_symbols - {x, y}
        rhs = _simplify_variable_coeff(sol.rhs, syms, solvefun, x)
        return Eq(f(x, y), rhs)

    dummy = Function('d')
    h = (c/b).subs(y, dummy(x))
    sol = dsolve(dummy(x).diff(x) - h, dummy(x))
    if isinstance(sol, list):
        sol = sol[0]
    solsym = sol.free_symbols - h.free_symbols - {x, y}
    if len(solsym) == 1:
        solsym = solsym.pop()
        etat = (solve(sol, solsym)[0]).subs(dummy(x), y)
        ysub = solve(eta - etat, y)[0]
        deq = (b*(f(x).diff(x)) + d*f(x) - e).subs(y, ysub)
        final = (dsolve(deq, f(x), hint='1st_linear')).rhs
        if isinstance(final, list):
            final = final[0]
        finsyms = final.free_symbols - deq.free_symbols - {x, y}
        rhs = _simplify_variable_coeff(final, finsyms, solvefun, etat)
        return Eq(f(x, y), rhs)

    else:
        raise NotImplementedError("Cannot solve the partial differential equation due"
            " to inability of constantsimp")

def pde_2nd_linear_constant_coeff_separable(eq, func, order, match, solvefun):
    return Eq(func, match.solve())


def _simplify_variable_coeff(sol, syms, func, funcarg):
    r"""
    Helper function to replace constants by functions in 1st_linear_variable_coeff
    """
    eta = Symbol("eta")
    if len(syms) == 1:
        sym = syms.pop()
        final = sol.subs(sym, func(funcarg))

    else:
        for sym in syms:
            final = sol.subs(sym, func(funcarg))

    return simplify(final.subs(eta, funcarg))


def pde_separate(eq, fun, sep, strategy='mul'):
    """Separate variables in partial differential equation either by additive
    or multiplicative separation approach. It tries to rewrite an equation so
    that one of the specified variables occurs on a different side of the
    equation than the others.

    :param eq: Partial differential equation

    :param fun: Original function F(x, y, z)

    :param sep: List of separated functions [X(x), u(y, z)]

    :param strategy: Separation strategy. You can choose between additive
        separation ('add') and multiplicative separation ('mul') which is
        default.

    Examples
    ========

    >>> from sympy import E, Eq, Function, pde_separate, Derivative as D
    >>> from sympy.abc import x, t
    >>> u, X, T = map(Function, 'uXT')

    >>> eq = Eq(D(u(x, t), x), E**(u(x, t))*D(u(x, t), t))
    >>> pde_separate(eq, u(x, t), [X(x), T(t)], strategy='add')
    [exp(-X(x))*Derivative(X(x), x), exp(T(t))*Derivative(T(t), t)]

    >>> eq = Eq(D(u(x, t), x, 2), D(u(x, t), t, 2))
    >>> pde_separate(eq, u(x, t), [X(x), T(t)], strategy='mul')
    [Derivative(X(x), (x, 2))/X(x), Derivative(T(t), (t, 2))/T(t)]

    See Also
    ========
    pde_separate_add, pde_separate_mul
    """

    do_add = False
    if strategy == 'add':
        do_add = True
    elif strategy == 'mul':
        do_add = False
    else:
        raise ValueError('Unknown strategy: %s' % strategy)

    if isinstance(eq, Equality):
        if eq.rhs != 0:
            return pde_separate(Eq(eq.lhs - eq.rhs, 0), fun, sep, strategy)
    else:
        return pde_separate(Eq(eq, 0), fun, sep, strategy)

    if eq.rhs != 0:
        raise ValueError("Value should be 0")

    # Handle arguments
    orig_args = list(fun.args)
    subs_args = [arg for s in sep for arg in s.args]

    if do_add:
        functions = reduce(operator.add, sep)
    else:
        functions = reduce(operator.mul, sep)

    # Check whether variables match
    if len(subs_args) != len(orig_args):
        raise ValueError("Variable counts do not match")
    # Check for duplicate arguments like  [X(x), u(x, y)]
    if has_dups(subs_args):
        raise ValueError("Duplicate substitution arguments detected")
    # Check whether the variables match
    if set(orig_args) != set(subs_args):
        raise ValueError("Arguments do not match")

    # Substitute original function with separated...
    result = eq.lhs.subs(fun, functions).doit()

    # Divide by terms when doing multiplicative separation
    if not do_add:
        eq = 0
        for i in result.args:
            eq += i/functions
        result = eq

    svar = subs_args[0]
    dvar = subs_args[1:]
    return _separate(result, svar, dvar)


def pde_separate_add(eq, fun, sep):
    """
    Helper function for searching additive separable solutions.

    Consider an equation of two independent variables x, y and a dependent
    variable w, we look for the product of two functions depending on different
    arguments:

    `w(x, y, z) = X(x) + y(y, z)`

    Examples
    ========

    >>> from sympy import E, Eq, Function, pde_separate_add, Derivative as D
    >>> from sympy.abc import x, t
    >>> u, X, T = map(Function, 'uXT')

    >>> eq = Eq(D(u(x, t), x), E**(u(x, t))*D(u(x, t), t))
    >>> pde_separate_add(eq, u(x, t), [X(x), T(t)])
    [exp(-X(x))*Derivative(X(x), x), exp(T(t))*Derivative(T(t), t)]

    """
    return pde_separate(eq, fun, sep, strategy='add')


def pde_separate_mul(eq, fun, sep):
    """
    Helper function for searching multiplicative separable solutions.

    Consider an equation of two independent variables x, y and a dependent
    variable w, we look for the product of two functions depending on different
    arguments:

    `w(x, y, z) = X(x)*u(y, z)`

    Examples
    ========

    >>> from sympy import Function, Eq, pde_separate_mul, Derivative as D
    >>> from sympy.abc import x, y
    >>> u, X, Y = map(Function, 'uXY')

    >>> eq = Eq(D(u(x, y), x, 2), D(u(x, y), y, 2))
    >>> pde_separate_mul(eq, u(x, y), [X(x), Y(y)])
    [Derivative(X(x), (x, 2))/X(x), Derivative(Y(y), (y, 2))/Y(y)]

    """
    return pde_separate(eq, fun, sep, strategy='mul')


def _separate(eq, dep, others):
    """Separate expression into two parts based on dependencies of variables."""

    # FIRST PASS
    # Extract derivatives depending our separable variable...
    terms = set()
    for term in eq.args:
        if term.is_Mul:
            for i in term.args:
                if i.is_Derivative and not i.has(*others):
                    terms.add(term)
                    continue
        elif term.is_Derivative and not term.has(*others):
            terms.add(term)
    # Find the factor that we need to divide by
    div = set()
    for term in terms:
        ext, sep = term.expand().as_independent(dep)
        # Failed?
        if sep.has(*others):
            return None
        div.add(ext)
    # FIXME: Find lcm() of all the divisors and divide with it, instead of
    # current hack :(
    # https://github.com/sympy/sympy/issues/4597
    if len(div) > 0:
        # double sum required or some tests will fail
        eq = Add(*[simplify(Add(*[term/i for i in div])) for term in eq.args])
    # SECOND PASS - separate the derivatives
    div = set()
    lhs = rhs = 0
    for term in eq.args:
        # Check, whether we have already term with independent variable...
        if not term.has(*others):
            lhs += term
            continue
        # ...otherwise, try to separate
        temp, sep = term.expand().as_independent(dep)
        # Failed?
        if sep.has(*others):
            return None
        # Extract the divisors
        div.add(sep)
        rhs -= term.expand()
    # Do the division
    fulldiv = reduce(operator.add, div)
    lhs = simplify(lhs/fulldiv).expand()
    rhs = simplify(rhs/fulldiv).expand()
    # ...and check whether we were successful :)
    if lhs.has(*others) or rhs.has(dep):
        return None
    return [lhs, rhs]


@dataclass(frozen=True)
class _Cond:
    """TODO"""
    var: Symbol
    var_val: Expr
    unknown_val: Expr
    diff: int = 0

    @classmethod
    def parse(cls, cond_lhs, unknown_val, unknown):
        """TODO"""
        if unknown_val.has(unknown.func):
            raise NotImplementedError(f"{cond_lhs} = {unknown_val} relates the unknown to itself (eg. periodic), which is not yet supported")
        if isinstance(cond_lhs, Limit):
            raise NotImplementedError(f"limits are not currently supported; received {cond_lhs}")

        if isinstance(cond_lhs, Subs):
            if len(cond_lhs.variables) != 1:
                raise ValueError(f"{cond_lhs} must fix exactly one variable")
            (var,), (var_val,) = cond_lhs.variables, cond_lhs.point
            cond_derivative = cond_lhs.expr
            if not isinstance(cond_derivative, Derivative):
                raise ValueError(f"unrecognized condition: {cond_lhs}")
            base, diff = cond_derivative.expr, cond_derivative.derivative_count
            if set(cond_derivative.variables) != {var}:
                raise ValueError(f"{cond_lhs} differentiates in the wrong variable")
        else:
            diff_var = None
            if isinstance(cond_lhs, Derivative):
                base, diff = cond_lhs.expr, cond_lhs.derivative_count
                if len(set(cond_lhs.variables)) != 1:
                    raise ValueError(f"mixed derivatives not supported: {cond_lhs}")
                diff_var = cond_lhs.variables[0]
            else:
                base, diff = cond_lhs, 0
            fixed_vars = [(v, p) for v, p in zip(unknown.args, base.args) if v != p]
            if len(fixed_vars) != 1:
                raise ValueError(f"{cond_lhs} must fix exactly one variable")
            var, var_val = fixed_vars[0]
            if diff_var is not None and diff_var != var_val:
                raise ValueError(f"{cond_lhs} differentiates in the wrong variable")

        if var_val.is_infinite or unknown_val.is_infinite:
            raise NotImplementedError("conditions involving infinity are not currently supported")
        if base.func != unknown.func:
            raise ValueError(f"{cond_lhs} is not a condition on {unknown.func}")
        if unknown_val.has(var):
            raise ValueError(f"{cond_lhs} fixes {var}, so the resulting value {unknown_val} cannot depend on {var}")
        return cls(var, var_val, unknown_val, diff)

    def eval(self, expr):
        """TODO"""
        lhs = expr.diff(self.var, self.diff).subs(self.var, self.var_val)
        return Eq(lhs, self.unknown_val, evaluate=False)

@dataclass
class _BVP(ABC):
    """TODO"""
    unknown: Function
    pde: Expr
    coeffs: dict
    conds: dict
    domains: dict

    def validity_domain(self):
        """TODO"""
        return And(*[
            relation
            for var, interval in self.domains.items()
            for relation in And.make_args(interval.as_relational(var))
            if not relation.has(S.Infinity, S.NegativeInfinity)],
        )

    @classmethod
    def from_input(cls, pde, unknown, conditions=None, solvefun=None):
        """TODO"""
        cls.check_variables(unknown)
        coeffs = cls.find_coefficients(pde, unknown)
        orders = cls.check_derivative_orders(coeffs, unknown)
        pde_class = cls.classify(coeffs)

        solvefun = solvefun or Function("F")
        taken_function_names = {
            fun.func.__name__
            for expr in (pde, *map(S, (conditions or {}).values()))
            for fun in expr.atoms(AppliedUndef)
        }

        def new_func():
            """TODO"""
            index = 1
            while f"{solvefun.__name__}{index}" in taken_function_names:
                index += 1
            name = f"{solvefun.__name__}{index}"
            taken_function_names.add(name)
            return Function(name)

        conds, domains = pde_class.parse_conds(orders, coeffs, unknown, conditions, new_func)
        return pde_class(unknown, pde, coeffs, conds, domains)

    @staticmethod
    def check_variables(unknown):
        """TODO"""
        assert len(unknown.args) == 2, f"expected exactly 2 independent variables, received {unknown}"

        for var in unknown.args:
            if var.is_Symbol is False:
                raise ValueError(f"{unknown} must be applied to plain symbols, not {var}")
            if var.is_rational:
                raise ValueError(f"{var} must be a continuous variable")
            if var.is_real is False:
                raise NotImplementedError(f"{var} must be real, complex variables are not yet implemented")

    @staticmethod
    def find_coefficients(pde_expr, unknown):
        """TODO"""
        x, y = unknown.args
        terms = {
            "u_xx": unknown.diff(x, 2),
            "u_xy": unknown.diff(x, y),
            "u_yy": unknown.diff(y, 2),
            "u_x": unknown.diff(x),
            "u_y": unknown.diff(y),
            "u": unknown,
        }

        coeffs = {name: expand(pde_expr).coeff(term) for name, term in terms.items()}
        coeffs["c"] = expand(
            pde_expr - sum(coeffs[n] * term for n, term in terms.items())
        )

        if any(coeff.has(unknown.func) for coeff in coeffs.values()):
            raise NotImplementedError(f"received {pde_expr} = 0, which is not linear in {unknown}")
        if any(coeff.has(x, y) for name, coeff in coeffs.items() if name != "c"):
            raise NotImplementedError(f"received {pde_expr} = 0, which does not have constant coefficients")

        leading_terms = ["u_xx", "u_yy", "u_xy"]
        complex_leading_terms = [coeffs[term] for term in leading_terms if coeffs[term].is_real is False]
        if complex_leading_terms:
            coeffs = {name: simplify(coeff * conjugate(complex_leading_terms[0])) for name, coeff in coeffs.items()}
            if any(coeffs[term].is_real is False for term in leading_terms):
                raise NotImplementedError("second-order coefficients must be real, up to a common factor")
        return coeffs

    @staticmethod
    def check_derivative_orders(coeffs, unknown):
        """TODO"""
        x, y = unknown.args
        x_order = (
            2 if coeffs["u_xx"] != 0
            else 1 if (coeffs["u_x"] != 0 or coeffs["u_xy"] != 0)
            else 0
        )
        y_order = (
            2 if coeffs["u_yy"] != 0
            else 1 if (coeffs["u_y"] != 0 or coeffs["u_xy"] != 0)
            else 0
        )

        if x_order == 0 or y_order == 0:
            raise ValueError("expected a PDE, received an ODE")
        if x_order < 2 and y_order < 2:
            raise NotImplementedError("expected a second derivative in at least one variable")
        return {x: x_order, y: y_order}

    @staticmethod
    def classify(coeffs):
        """TODO"""
        discriminant = simplify(coeffs["u_xy"] ** 2 - 4 * coeffs["u_xx"] * coeffs["u_yy"])
        if discriminant.is_zero:
            return _ParabolicBVP
        if discriminant.is_nonnegative:
            return _HyperbolicBVP
        if discriminant.is_nonpositive:
            return _EllipticBVP
        raise ValueError(f"failed to classify the PDE due to ambiguous sign of the discriminant {discriminant}; declare the coefficients' signs")

    @staticmethod
    @abstractmethod
    def check_conditions(conds, orders, coeffs, unknown, new_func):
        """TODO"""

    @classmethod
    def parse_conds(cls, orders, coeffs, unknown, conditions, new_func):
        """TODO"""
        parsed = [_Cond.parse(lhs, S(val), unknown) for lhs, val in (conditions or {}).items()]
        conds = {var: tuple(c for c in parsed if c.var == var) for var in unknown.args}

        for var in unknown.args:
            if len(conds[var]) > orders[var]:
                raise ValueError(f"too many conditions given for {var}; expected at most {orders[var]} conditions, received {conds[var]}")
            if any(cond.diff >= orders[var] for cond in conds[var]):
                raise ValueError(f"conditions in {var} may only fix derivatives of order below {orders[var]}, as higher ones are already determined by the PDE")

        return cls.check_conditions(conds, orders, coeffs, unknown, new_func)

    @staticmethod
    def fill_bcs(var_conds, var, other_var, new_func):
        """TODO"""
        if len(var_conds) == 2:
            return var_conds

        L = Dummy("L", positive=True)  # arbitrary length

        if len(var_conds) == 1:
            x0 = var_conds[0].var_val
            return var_conds + (_Cond(var, x0 + L, new_func()(other_var)),)

        x0 = S.Zero
        return (
            _Cond(var, x0, new_func()(other_var)),
            _Cond(var, x0 + L, new_func()(other_var)),
        )

    @staticmethod
    def fill_ics(var_conds, var_order, var, other_var, new_func):
        """TODO"""
        diffs = [cond.diff for cond in var_conds]
        t0 = var_conds[0].var_val if var_conds else S.Zero

        return var_conds + tuple(
            _Cond(var, t0, new_func()(other_var), k)
            for k in range(var_order)
            if k not in diffs
        )

    @staticmethod
    def can_be_spatial(var_conds):
        """TODO"""
        return len(var_conds) == len({cond.var_val for cond in var_conds})

    @staticmethod
    def parse_interval(var_conds, direction=S.Zero):
        """TODO"""
        fixed_points = {cond.var_val for cond in var_conds}
        if len(fixed_points) == 2:
            start, end = Min(*fixed_points), Max(*fixed_points)
            if isinstance(start, Min):
                raise ValueError(f"cannot order the points {fixed_points}")
            return Interval(start, end)

        (start,) = fixed_points
        if direction.is_zero:
            return Interval(-S.Infinity, S.Infinity)
        if direction.is_nonnegative:
            return Interval(start, S.Infinity)
        if direction.is_nonpositive:
            return Interval(-S.Infinity, start)
        raise ValueError(f"cannot tell the sign of {direction}, which decides whether the solution holds above or below {start}; declare the coefficients' signs")

class _EllipticBVP(_BVP):
    """TODO"""
    @staticmethod
    def check_conditions(conds, orders, coeffs, unknown, new_func):
        """TODO"""
        x, y = unknown.args
        domains = {}
        for var, other_var in ((x, y), (y, x)):
            if not _BVP.can_be_spatial(conds[var]):
                raise ValueError(f"elliptic PDEs cannot accept multiple conditions fixing {var} at a single point; received {conds[var]}")

            conds[var] = _BVP.fill_bcs(conds[var], var, other_var, new_func)
            domains[var] = _BVP.parse_interval(conds[var])

        if coeffs["u"].is_zero and all(cond.diff > 0 for var_conds in conds.values() for cond in var_conds):
            raise ValueError("every condition fixes a derivative and the PDE has no u term, so any constant can be added to a solution, and one exists only if the data's net flux is zero")
        return conds, domains


class _ParabolicBVP(_BVP):
    """TODO"""
    @staticmethod
    def check_conditions(conds, orders, coeffs, unknown, new_func):
        """TODO"""
        if coeffs["u_xy"] != 0:
            raise NotImplementedError("parabolic PDEs with a mixed derivative are not yet supported")

        (space_var,) = [var for var in conds if orders[var] == 2]
        (time_var,) = [var for var in conds if orders[var] == 1]

        if not _BVP.can_be_spatial(conds[space_var]):
            raise ValueError(f"parabolic PDEs cannot accept multiple conditions fixing the spatial variable {space_var} at a single point")

        conds[space_var] = _BVP.fill_bcs(conds[space_var], space_var, time_var, new_func)
        conds[time_var] = _BVP.fill_ics(conds[time_var], orders[time_var], time_var, space_var, new_func,)

        x, y = unknown.args
        time_coeff, space_coeff = (
            (coeffs["u_y"], coeffs["u_xx"]) if space_var == x
            else (coeffs["u_x"], coeffs["u_yy"])
        )
        direction = re(simplify(-space_coeff / time_coeff))
        return conds, {space_var: _BVP.parse_interval(conds[space_var]), time_var: _BVP.parse_interval(conds[time_var], direction)}

class _HyperbolicBVP(_BVP):
    """TODO"""
    @staticmethod
    def check_conditions(conds, orders, coeffs, unknown, new_func):
        """TODO"""
        x, y = unknown.args
        space_vars = [var for var in conds if len({c.var_val for c in conds[var]}) == 2]

        # both wave speeds share a sign, so x behaves like a second time variable
        if (coeffs["u_xx"] * coeffs["u_yy"]).is_nonnegative:
            if space_vars:
                raise ValueError(f"conditions at two points in {space_vars[0]} overlap: when both wave speeds have the same sign, each variable must behave like a time, taking both its conditions at a single point")

            for var, other_var in ((x, y), (y, x)):
                conds[var] = _BVP.fill_ics(conds[var], orders[var], var, other_var, new_func)
            return conds, {
                x: _BVP.parse_interval(conds[x], coeffs["u_xy"] * (coeffs["u_xx"] + coeffs["u_yy"])),
                y: _BVP.parse_interval(conds[y], S.One),  # convention: y >= y0; the mirrored quadrant is equally valid
            }

        if len(space_vars) == 2:
            raise ValueError(f"conditions at two points in both {x} and {y} can resonate, leaving no solutions or infinitely many")
        if not space_vars:
            space_vars = [var for var, var_conds in conds.items() if _BVP.can_be_spatial(var_conds)]
            if not space_vars:
                raise ValueError(f"conditions at a single point in both {x} and {y} overdetermine the problem")
            if len(space_vars) == 2:
                raise ValueError(f"cannot tell which variable is spatial; give two conditions in either {x} or {y}")

        (space_var,) = space_vars
        (time_var,) = (var for var in conds if var != space_var)

        conds[space_var] = _BVP.fill_bcs(conds[space_var], space_var, time_var, new_func)
        conds[time_var] = _BVP.fill_ics(conds[time_var], orders[time_var], time_var, space_var, new_func)

        return conds, {space_var: _BVP.parse_interval(conds[space_var]), time_var: _BVP.parse_interval(conds[time_var])}

@dataclass
class _BVPSolver(ABC):
    """TODO"""
    bvp: _BVP
    _matched: bool | None = field(default=None, init=False)

    def matches(self):
        """TODO"""
        if self._matched is None:
            self._matched = self._matches()
        return self._matched

    def solve(self):
        """TODO"""
        if not self.matches():
            raise ValueError(f"{type(self).__name__} does not apply to this BVP")
        return self._solve()

    @abstractmethod
    def _matches(self):
        """Without solving, checks if the current BVP meets this solution method's requirements."""

    @abstractmethod
    def _solve(self):
        """Solves the current BVP, to be called only if _matched is true."""


@dataclass
class _SeparationOfVariables(_BVPSolver):
    """TODO"""
    def _matches(self):
        """TODO"""
        if not self.bvp.coeffs["u_xy"].is_zero:
            return False
        if not self.bvp.coeffs["c"].is_zero:
            return False

        homogeneous_spatial = [
            var
            for var, var_conds in self.bvp.conds.items()
            if len(var_conds) == 2
            and _BVP.can_be_spatial(var_conds)
            and all(cond.unknown_val.is_zero for cond in var_conds)
        ]
        if not homogeneous_spatial:
            return False

        self.space_var = homogeneous_spatial[0]
        (self.other_var,) = [var for var in self.bvp.conds if var != self.space_var]
        assert (len(self.bvp.unknown.args) == 2), "This solution method expects an unknown with exactly two independent variables."
        if self.space_var == self.bvp.unknown.args[0]:
            self.coeffs = self.bvp.coeffs
        else:
            swap = str.maketrans("xy", "yx")
            self.coeffs = {name.translate(swap): coeff for name, coeff in self.bvp.coeffs.items()}
        return True

    def _solve(self):
        """TODO"""
        reverse_subs = self.normalize() # 0 <= x <= L and (0 <= y <= H or (t >= 0 or t <= 0 or -oo < t < oo))
        self.n = Dummy("n", integer=True, positive=True)  # summation index

        terms = []
        for eigenvalue, eigenfunction in self.find_eigenpairs():
            if len({cond.var_val for cond in self.conds[self.y]}) == 2:
                self.check_resonance(eigenvalue)
            terms.append(self.sol_eigenpair(eigenvalue, eigenfunction))

        return Piecewise((Add(*terms).subs(reverse_subs), self.bvp.validity_domain()))

    def normalize(self):
        """TODO"""
        self.x = Dummy(str(self.space_var), real=True)
        self.y = Dummy(str(self.other_var), real=True)

        shifts = {var: Min(*{cond.var_val for cond in var_conds}) for var, var_conds in self.bvp.conds.items()}
        dummies = {self.space_var: self.x, self.other_var: self.y}
        shifted = {var: dummies[var] + shifts[var] for var in dummies}

        self.conds = {
            dummies[var]: tuple(
                _Cond(dummies[var], cond.var_val - shifts[var], cond.unknown_val.subs(shifted), cond.diff)
                for cond in var_conds
            )
            for var, var_conds in self.bvp.conds.items()
        }

        self.domains = {
            dummies[var]: Interval(*(
                    p if p.is_infinite
                    else p - shifts[var]
                    for p in (interval.start, interval.end)
            ))
            for var, interval in self.bvp.domains.items()
        }

        reverse_subs = {dummies[var]: var - shifts[var] for var in dummies}
        return reverse_subs

    def find_eigenpairs(self):
        """TODO"""
        eigenpairs = []
        a, b = self.coeffs["u_xx"], self.coeffs["u_x"]
        x = self.x

        C1, C2 = symbols("C1 C2", cls=Dummy)
        mu = Dummy("mu")
        lamda = self.coeffs["u"] - b**2 / (4 * a) - a * mu**2

        # Special case: mu = 0
        X_sol = exp(-b * x / (2 * a)) * (C1 + C2 * x)
        BC_matrix, rhs = linear_eq_to_matrix([cond.eval(X_sol) for cond in self.conds[x]], (C1, C2))
        assert (rhs.is_zero_matrix), f"Expected homogeneous boundary conditions, received{rhs}."
        det = expand(BC_matrix.det())
        if det.is_zero is None:
            det = simplify(det)
        if det.is_zero:
            eigenvalue = lamda.subs(mu, S.Zero)
            eigenfunction = self.find_eigenfunction(BC_matrix, X_sol, (C1, C2))
            eigenpairs.append((eigenvalue, eigenfunction))

        # Regular case: mu != 0
        X_sol = exp(-b * x / (2 * a)) * (C1 * cos(mu * x) + C2 * sin(mu * x))
        BC_matrix, _ = linear_eq_to_matrix([cond.eval(X_sol) for cond in self.conds[x]], (C1, C2))
        roots = solveset(Eq(self.strip_nonzero_factors(BC_matrix.det().rewrite(exp)), S.Zero), mu, S.Complexes)

        implicit_roots, explicit_roots = sift(sorted(Union.make_args(roots), key=str), lambda s: isinstance(s, ConditionSet), binary=True)
        if implicit_roots:
            raise NotImplementedError(f"implicit eigenvalue roots are not yet implemented; found {implicit_roots}")

        for mu_n in self.find_mu_values(Union(*explicit_roots)):
            eigenvalue = lamda.subs(mu, mu_n)
            eigenfunction = self.find_eigenfunction(expand(BC_matrix.subs(mu, mu_n)), X_sol.subs(mu, mu_n), (C1, C2))
            eigenpairs.append((eigenvalue, eigenfunction))
        return eigenpairs

    @staticmethod
    def find_eigenfunction(BC_matrix_n, X_n, X_consts):
        """TODO"""
        p, q = BC_matrix_n.row(0)
        return X_n.subs(dict(zip(X_consts, (q, -p))))

    @staticmethod
    def strip_nonzero_factors(expr):
        """TODO"""
        return Mul(*[
                factor
                for factor in Mul.make_args(factor_terms(expand(expr)))
                if not (factor.func == exp or (factor.is_number and factor.is_nonzero))
        ])

    def find_mu_values(self, roots):
        """TODO"""
        mu_values = []

        for s in sorted(Union.make_args(roots), key=str):
            if isinstance(s, FiniteSet):
                for mu_n in s:
                    # Deduplicate roots, since mu**2 = (-mu)**2
                    if (not re(mu_n).is_nonpositive or (re(mu_n).is_zero and not im(mu_n).is_nonpositive)) and -mu_n not in mu_values:
                        mu_values.append(mu_n)

            elif s == S.Integers or s == S.Naturals:
                mu_n = self.n
                mu_values.append(mu_n)

            elif isinstance(s, ImageSet):
                indexes = s.lamda.variables
                assert (len(indexes) == 1), f"Expected exactly one index for a family of eigenvalue roots, received {s}"

                (index,) = indexes
                mu_n = s.lamda.expr

                positive_re = solveset(re(mu_n) > 0, index, s.base_set)
                zero_re = solveset(Eq(re(mu_n), 0), index, s.base_set)
                positive_im = solveset(im(mu_n) > 0, index, s.base_set)
                index_range = Union(positive_re, Intersection(zero_re, positive_im))
                if index_range.is_empty:
                    continue

                if index_range.inf.is_finite:
                    mu_n = mu_n.subs(index, index_range.inf + self.n - 1)
                else:
                    mu_n = mu_n.subs(index, index_range.sup - self.n + 1)
                mu_values.append(mu_n)

            else:
                raise ValueError(f"received unexpected object {s} for an eigenvalue root; expected FiniteSet, ImageSet, or S.Integers")
        return mu_values

    def sol_eigenpair(self, lamda_n, X_n):
        """TODO"""
        generic_term = self.sol_Y(lamda_n, X_n)
        if not lamda_n.has(self.n): # isolated eigenvalue
            return generic_term

        a, b = self.coeffs["u_yy"], self.coeffs["u_y"]
        characteristic_discriminant = b**2 - 4*a*lamda_n
        special_indices = solveset(characteristic_discriminant, self.n, S.Naturals)
        if not isinstance(special_indices, FiniteSet):
            special_indices = S.EmptySet # can't list them (symbolic coefficients), so treat every n as generic

        relation_roots = [solveset(relation.lhs - relation.rhs, self.n, S.Naturals) for relation in generic_term.atoms(Relational)]
        exceptional_indices = Union(special_indices, *[roots for roots in relation_roots if isinstance(roots, FiniteSet)])
        exceptional_terms = [
            self.sol_Y(lamda_n.subs(self.n, index), X_n.subs(self.n, index)) if index in special_indices
            else generic_term.subs(self.n, index)
            for index in exceptional_indices
        ]

        range_bounds = sorted(exceptional_indices + {0, S.Infinity})
        range_sums = [self.range_sum(generic_term, range_bounds[i] + 1, range_bounds[i + 1] - 1) for i in range(len(range_bounds) - 1)]
        return Add(*exceptional_terms, *range_sums)

    def range_sum(self, term, start, end):
        """TODO"""
        j = Dummy("j", integer=True, nonnegative=True)
        shifts = [(start + j, self.n - start)]
        if end.is_finite:
            shifts.append((end - j, end - self.n))

        for shift, unshift in shifts:
            shifted_term = term.subs(self.n, shift)
            imaginary_sqrts = {root: (I*sqrt(-root.base))**(2*root.exp) for root in shifted_term.atoms(Pow) if root.exp.is_Rational and root.exp.q == 2 and root.base.is_negative}
            term = shifted_term.xreplace(imaginary_sqrts).subs(j, unshift)

        if term == 0:
            return S.Zero

        max_explicit_terms = 2
        if end - start + 1 <= max_explicit_terms:
            return Add(*[term.subs(self.n, index) for index in range(start, end + 1)])

        return Sum(term, (self.n, start, end))

    def sol_Y(self, lamda_n, X_n):
        """TODO"""
        Y_n, Y_consts = self.Y_basis(lamda_n)
        norm = self.weighted_inner_product(X_n, X_n)
        projected_conds = [replace(cond, unknown_val = factor(cancel(self.weighted_inner_product(cond.unknown_val, X_n) / norm))) for cond in self.conds[self.y]]
        M, data = linear_eq_to_matrix([cond.eval(Y_n) for cond in projected_conds], Y_consts)
        if M.det().is_zero:
            raise ValueError("the conditions are resonant, so there is no solution or infinitely many")
        Y_n = Y_n.subs(dict(zip(Y_consts, M.LUsolve(data))))
        return factor_terms(X_n * Y_n)

    def Y_basis(self, lam):
        """TODO"""
        a, b = self.coeffs["u_yy"], self.coeffs["u_y"]
        y = self.y
        C3, C4 = symbols("C3 C4", cls=Dummy)

        if a == 0:
            return C3 * exp(-lam * y/b), (C3,)
        assert len(self.conds[y]) == 2, f"Expected 2 conditions in Y, received {self.conds[y]}"

        E = exp(-b*y / (2*a))
        d = b**2 - 4*a*lam
        k = sqrt(d) / (2*a)
        Y1 = E*y if d.is_zero else E*sinh(k*y) / k # Y1(0) = 0, Y1'(0) = 1
        Y2 = E*cosh(k*y) + (b/(2*a)) * Y1 # Y2(0) = 1, Y2'(0) = 0

        basis = [] # B(0) = B'(0) = 0
        for cond in self.conds[y]:
            s = y - cond.var_val
            Y1_s, Y2_s = Y1.subs(y, s), Y2.subs(y, s)
            p, q = cond.eval(Y1_s).lhs, cond.eval(Y2_s).lhs
            B = q * Y1_s - p * Y2_s # cond.eval(B) = 0 at s = 0
            basis.append(B)

        return C3*basis[0] + C4*basis[1], (C3, C4)

    def weighted_inner_product(self, g, h):
        """TODO"""
        x = self.x
        limits = (x, self.domains[x].start, self.domains[x].end)
        weight = exp(self.coeffs["u_x"] * x / self.coeffs["u_xx"])
        return integrate(g * h * weight, limits)

    def check_resonance(self, lamda_n):
        """TODO"""
        if not lamda_n.has(self.n):
            return

        a, b = self.coeffs["u_yy"], self.coeffs["u_y"]
        H = self.domains[self.y].end
        m = Dummy("m", integer=True, positive=True)
        cond_diffs = {cond.diff for cond in self.conds[self.y]}

        resonance_conditions = []
        if len(cond_diffs) == 1:
            resonance_conditions.append(lamda_n - (b**2 / (4*a) + a * (m*pi/H)**2))
        elif b.is_zero:
            resonance_conditions.append(lamda_n - a * ((m - S.Half)*pi/H)**2)
        else:
            return

        if cond_diffs == {1}: # all Neumann
            zero_root_indices = solveset(lamda_n, self.n, S.Naturals)
            if zero_root_indices != S.EmptySet and isinstance(zero_root_indices, FiniteSet):
                raise ValueError("the conditions are resonant, so there is no solution or infinitely many")

        if any(self.has_positive_integer_solutions(condition, self.n, m) for condition in resonance_conditions):
            raise ValueError("the conditions are resonant, so there is no solution or infinitely many")

    @staticmethod
    def has_positive_integer_solutions(expr, n, m):
        """TODO"""
        try:
            poly = Poly(expand(expr), n, m)
        except PolynomialError:
            return None

        poly = Poly(expand(expr / poly.LC()), n, m)
        if poly.domain not in (ZZ, QQ):
            return None  # some coefficient is irrational or symbolic beyond a shared factor

        _, integer_poly = poly.clear_denoms()
        return bool(diophantine(integer_poly.as_expr(), syms=(n, m)))
