"""
This module contains pdesolve() and different helper functions that it
uses.It is heavily inspired by the ode module and hence the basic
infrastructure remains the same.

**Functions in this module**

    These are the user functions in this module:

    - pdesolve() - Solves PDEs.
    - classify_pde() - Classifies PDEs into possible hints for dsolve().
    - pde_separate() - Separate variables in partial differential equation either by
                       additive or multiplicative separation approach.

    These are the user functions in the ode module that are used here.
    - ode_order() - Returns the order (degree) of an ODE. However this can
      be used to find the order of a PDE as well.

    These are the helper functions in this module
    - pde_separate_add() - Helper function for searching additive separable solutions.
    - pde_separate_mul() - Helper function for searching multiplicative
                         separable solutions.

    These are the helper functions in the ode module that are used here.
    - preprocess - prepare the equation and detect function to solve for.

**Currently implemented solver methods**

The following methods are implemented for solving partial differential
equations.  See the docstrings of the various pde_hint() functions for
more information on each (run help(pde)):

  - 1st order linear homogeneous partial differential equations
    with constant coefficients.

"""
from copy import deepcopy
from itertools import combinations_with_replacement

from sympy import Eq, Equality
from sympy.simplify import simplify
from sympy.core import Add, C, S, Mul, Pow, oo
from sympy.core.compatibility import reduce
from sympy.core.function import Function, Derivative, expand, diff
from sympy.core.numbers import Rational
from sympy.core.symbol import Symbol, Wild, Dummy, symbols
from sympy.functions import exp
from sympy.utilities.iterables import has_dups

from sympy.solvers.ode import preprocess, ode_order
import operator

allhints = (
    "1st_linear_constant_coeff_homo",
    )

def pdesolve(eq, func=None, hint='default', dict=False, **kwargs):
    """
    Solves any (supported) kind of partial differential equation.

    **Usage**

        pdesolve(eq, f(x,y), hint) -> Solve partial differential equation
        eq for function f(x,y), using method hint.

    **Details**

        ``eq`` can be any supported partial differential equation (see
            the pde docstring for supported methods).  This can either
            be an Equality, or an expression, which is assumed to be
            equal to 0.

        ``f(x,y)`` is a function of two variables whose derivatives in that
            variable make up the partial differential equation. In many
            cases it is not necessary to provide this; it will be autodetected
            (and an error raised if it couldn't be detected).

        ``hint`` is the solving method that you want pdesolve to use.  Use
            classify_pde(eq, f(x,y)) to get all of the possible hints for
            a PDE.  The default hint, 'default', will use whatever hint
            is returned first by classify_pde().  See Hints below for
            more options that you can use for hint.

    **Hints**

        Aside from the various solving methods, there are also some
        meta-hints that you can pass to pdesolve():

        "default":
                This uses whatever hint is returned first by
                classify_pde(). This is the default argument to
                pdesolve().

        See also the classify_pde() docstring for more info on hints,
        and the ode docstring for a list of all supported hints.

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
          examples for how to use pdesolve().
        - pdesolve always returns an Equality class (except for the case
          when the hint is "all" or "all_Integral"). Note that it is not possible
          to get an explicit solution for f(x, y) as in the case of ODE's
        - Do help(ode.ode_hintname) to get help more information on a
          specific hint


    Examples
    ========

    >>> from sympy.solvers.pde import pdesolve
    >>> from sympy import Function, diff, Eq
    >>> from sympy.abc import x, y
    >>> f = Function('f')
    >>> u = f(x, y)
    >>> ux = u.diff(x)
    >>> uy = u.diff(y)
    >>> eq = Eq(1 + (2*(ux/u)) + (3*(uy/u)))
    >>> pdesolve(eq)
    f(x, y) == g(3*x - 2*y)*exp(-2*x/13 - 3*y/13)
    """
    prep = kwargs.get('prep', True)

    if isinstance(eq, Equality):
        eq = eq.lhs - eq.rhs

    # preprocess the equation and find func if not given
    if prep or func is None:
        eq, func = preprocess(eq, func)
        prep = False

    # Magic that should only be used internally.  Prevents classify_ode from
    # being called more than it needs to be by passing its results through
    # recursive calls.
    if kwargs.get('classify', True):
        hints = classify_pde(eq, func = func, dict = True, prep = prep)
    else:
        # Here is what all this means:
        #
        # hint:    The hint method given to pdesolve() by the user.
        # hints:   The dictionary of hints that match the PDE, along with other
        #          information (including the internal pass-through magic).
        # default: The default hint to return, the first hint from allhints
        #          that matches the hint; obtained from classify_pde().
        # match:   Dictionary containing the match dictionary for each hint
        #          (the parts of the PDE for solving).  When going through the
        #          hints in "all", this holds the match string for the current
        #          hint.
        # order:   The order of the PDE, as determined by ode_order().
        hints = kwargs.get('hint',
                           {'default': hint,
                            hint: kwargs['match'],
                            'order': kwargs['order']})
    if hints['order'] == 0:
        raise ValueError(
            str(eq) + " is not a differential equation in " + str(func))

    if not hints['default']:
        # classify_pde will set hints['default'] to None if no hints match
        if hint not in allhints and hint != 'default':
            raise ValueError("Hint not recognized: " + hint)
        elif hint not in hints['ordered_hints'] and hint != 'default':
            raise ValueError("PDE " + str(eq) + " does not match hint " + hint)
        else:
            raise NotImplementedError("pdesolve: Cannot solve " + str(eq))

    if hint == 'default':
        return pdesolve(eq, func, hint=hints['default'], simplify=simplify,
                      prep=prep, classify=False, order=hints['order'],
                      match=hints[hints['default']])
    elif hint not in allhints:  # and hint not in ('default', 'ordered_hints'):
        raise ValueError("Hint not recognized: " + hint)
    elif hint not in hints:
        raise ValueError("PDE " + str(eq) + " does not match hint " + hint)

    else:
        # convert the string into a function
        solvefunc = globals()['pde_' + hint]
        return solvefunc(eq, func, order=hints['order'], match=hints[hint])

def classify_pde(eq, func=None, dict=False, **kwargs):
    """
    Returns a tuple of possible pdesolve() classifications for a PDE.

    The tuple is ordered so that first item is the classification that
    pdesolve() uses to solve the PDE by default.  In general,
    classifications at the near the beginning of the list will produce
    better solutions faster than those near the end, thought there are
    always exceptions.  To make dsolve use a different classification,
    use pdesolve(PDE, func, hint=<classification>).  See also the pdesolve()
    docstring for different meta-hints you can use.

    If ``dict`` is true, classify_pde() will return a dictionary of
    hint:match expression terms. This is intended for internal use by
    pdesolve().  Note that because dictionaries are ordered arbitrarily,
    this will most likely not be in the same order as the tuple.

    You can get help on different hints by doing help(pde.pde_hintname),
    where hintname is the name of the hint without "_Integral".

    See sympy.pde.allhints or the sympy.pde docstring for a list of all
    supported hints that can be returned from classify_pde.


    Examples
    ========
    >>> from sympy.solvers.pde import classify_pde
    >>> from sympy import Function, diff, Eq
    >>> from sympy.abc import x, y
    >>> f = Function('f')
    >>> u = f(x, y)
    >>> ux = u.diff(x)
    >>> uy = u.diff(y)
    >>> eq = Eq(1 + (2*(ux/u)) + (3*(uy/u)))
    >>> classify_pde(eq)
    ('1st_linear_constant_coeff_homo',)
    """

    prep = kwargs.pop('prep', True)

    if func and len(func.args) != 2:
        raise NotImplementedError("Right now only partial"
        "differential equations of two variables are supported")

    if prep or func is None:
        prep, func_ = preprocess(eq, func)
        if func is None:
            func = func_

    if isinstance(eq, Equality):
        if eq.rhs != 0:
            return classify_pde(eq.lhs - eq.rhs, func)
        eq = eq.lhs

    f = func.func
    x = func.args[0]
    y = func.args[1]
    fx = f(x,y).diff(x)
    fy = f(x,y).diff(y)

    order = ode_order(eq, f(x,y))

    # hint:matchdict or hint:(tuple of matchdicts)
    # Also will contain "default":<default hint> and "order":order items.
    matching_hints = {'order': order}

    if not order:
        if dict:
            matching_hints["default"] = None
            return matching_hints
        else:
            return ()

    eq = expand(eq)

    a = Wild('a', exclude = [f(x,y)])
    b = Wild('b', exclude = [f(x,y), fx, fy])
    c = Wild('c', exclude = [f(x,y), fx, fy])
    d = Wild('d', exclude = [f(x,y), fx, fy])
    n = Wild('n', exclude = [x, y])
    # Try removing the smallest power of f(x,y)
    # from the highest partial derivatives of f(x,y)
    reduced_eq = None
    if eq.is_Add:
        var = set(combinations_with_replacement((x,y), order))
        dummyvar = deepcopy(var)
        power = None
        for i in var:
            coeff = eq.coeff(f(x,y).diff(*i))
            if coeff != 1:
                match = coeff.match(a*f(x,y)**n)
                if match and match[a]:
                    power = match[n]
                    dummyvar.remove(i)
                    break
            dummyvar.remove(i)
        for i in dummyvar:
            coeff = eq.coeff(f(x,y).diff(*i))
            if coeff != 1:
                match = coeff.match(a*f(x,y)**n)
                if match and match[a] and match[n] < power:
                    power = match[n]
        if power:
            den = f(x,y)**power
            reduced_eq = Add(*[arg/den for arg in eq.args])
        if not reduced_eq:
            reduced_eq = eq

    if order == 1:
        ## Linear first-order homogeneous partial-differential
        ## equation with constant coefficients
        r = reduced_eq.match(b*fx + c*fy + d*f(x,y))
        if r and all(isinstance(r[var], Rational) for var in r):
            r.update({'b': b, 'c': c, 'd': d})
            matching_hints["1st_linear_constant_coeff_homo"] = r

    # Order keys based on allhints.
    retlist = []
    for i in allhints:
        if i in matching_hints:
            retlist.append(i)

    if dict:
        # Dictionaries are ordered arbitrarily, so make note of which
        # hint would come first for pdesolve().  Use an ordered dict in Py 3.
        matching_hints["default"] = None
        matching_hints["ordered_hints"] = tuple(retlist)
        for i in allhints:
            if i in matching_hints:
                matching_hints["default"] = i
                break
        return matching_hints
    else:
        return tuple(retlist)

def pde_1st_linear_constant_coeff_homo(eq, func, order, match):
    r"""
    Solves a first order linear homogeneous
    partial differential equation with constant coefficients.

    The general form of this partial differential equation is of
    the form a*f(x,y).diff(x) + b*f(x,y).diff(y) + f(x,y) = 0
    where a, b and c are Rationals.

    The general solution of the differential equation, can be found
    put by the method of characteristics. It is given by
    f(x,y) = G(b*x - a*y)*exp(-c/(a**2 + b**2)*(a*x + b*y))

    Examples
    ========

    >>> from sympy.solvers.pde import (pde_1st_linear_constant_coeff_homo,
    ... pdesolve)
    >>> from sympy import Function, diff, pprint
    >>> from sympy.abc import x,y
    >>> f = Function('f')
    >>> pdesolve(f(x,y) + f(x,y).diff(x) + f(x,y).diff(y))
    f(x, y) == g(x - y)*exp(-x/2 - y/2)
    >>> pprint(pdesolve(f(x,y) + f(x,y).diff(x) + f(x,y).diff(y)))
                          x   y
                        - - - -
                          2   2
    f(x, y) = g(x - y)*e

    References
    ==========

    - Viktor Grigoryan, "Partial Differential Equations"
      Math 124A - Fall 2010, pp.7

    """
    f = func.func
    x = func.args[0]
    y = func.args[1]
    g = Function('g')
    b = match[match['b']]
    c = match[match['c']]
    d = match[match['d']]
    return Eq(f(x,y), exp(-S(d)/(b**2 + c**2)*(b*x + c*y))*g(c*x - b*y))

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
    [Derivative(X(x), x, x)/X(x), Derivative(T(t), t, t)/T(t)]

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
        assert ValueError('Unknown strategy: %s' % strategy)

    if isinstance(eq, Equality):
        if eq.rhs != 0:
            return pde_separate(Eq(eq.lhs - eq.rhs), fun, sep, strategy)
    assert eq.rhs == 0

    # Handle arguments
    orig_args = list(fun.args)
    subs_args = []
    for s in sep:
        for j in range(0, len(s.args)):
            subs_args.append(s.args[j])

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
    [Derivative(X(x), x, x)/X(x), Derivative(Y(y), y, y)/Y(y)]

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
    # http://code.google.com/p/sympy/issues/detail?id=1498
    if len(div) > 0:
        final = 0
        for term in eq.args:
            eqn = 0
            for i in div:
                eqn += term / i
            final += simplify(eqn)
        eq = final

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
