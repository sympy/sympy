"""This module contains various helpers for working with Partial Differential
Equations (PDEs)

"""

from sympy import Derivative, diff, Eq, Equality, Mul
from sympy.simplify import simplify

def pde_separate_add(eqn, func, sep_list):
    """
    Helper function for searching additive separable solutions.

    Consider an equation of two independent variables x, y and a dependent
    variable w, we look for the product of two functions depending on different
    arguments:

    w(x,t) = f(x) + g(t)

    Examples

    >>> from sympy import *
    >>> x, t = symbols('xt')
    >>> u, X, T = map(Function, 'uXT')

    >>> eq = Eq(Derivative(u(x, t), x), E**(u(x, t))*Derivative(u(x, t), t))
    >>> pde_separate_add(eq, u(x, t), [X(x), T(t)])
    [D(X(x), x)*exp(-X(x)), D(T(t), t)*exp(T(t))]

    """
    if isinstance(eqn, Equality):
        if eqn.rhs != 0:
            return pde_separate_add(Eq(eqn.lhs - eqn.rhs), func, sep_list)
    assert eqn.rhs == 0

    functions = 0
    subs_args = []
    orig_args = list(func.args)

    # FIXME: Check whether we are dealing with functions in func and sep_list
    for sep in sep_list:
        functions += sep
        for j in range(0, len(sep.args)):
            subs_args.append(sep.args[j])
    # FIXME!
    # Check whether first depends on only one variable

    # Check whether variables match
    if len(subs_args) != len(orig_args):
        raise ValueError("Variable counts do not match")
    # Check for duplicate arguments like  [X(x), u(x, y)]
    if len(subs_args) != len(set(subs_args)):
        raise ValueError("Duplicate substitution arguments detected")
    # Check whether the variables match
    if set(orig_args) != set(subs_args):
        raise ValueError("Arguments do not match")

    # Substitute...
    result = eqn.lhs.subs(func, functions)

    # Variable for which we separate

    svar = subs_args[0]
    dvar = subs_args[1:]
    return _separate(result, svar, dvar)

def pde_separate_mul(eq, fun, sep):
    """
    Helper function for searching multiplicative separable solutions.

    Consider an equation of two independent variables x, y and a dependent
    variable w, we look for the product of two functions depending on different
    arguments:

    w(x,t) = f(x)*g(t)

    Examples

    >>> from sympy import *
    >>> x, y = symbols('xy')
    >>> u, X, Y = map(Function, 'uXY')

    >>> eq = Eq(Derivative(u(x, y), x, 2), Derivative(u(x, y), y, 2))
    >>> pde_separate_mul(eq, u(x, y), [X(x), Y(y)])
    [D(X(x), x, x)/X(x), D(Y(y), y, y)/Y(y)]

    """

    if isinstance(eq, Equality):
        if eq.rhs != 0:
            return pde_separate_mul(Eq(eq.lhs - eq.rhs), fun, sep)
    assert eq.rhs == 0

    # Deal with function arguments
    subs_args = []
    orig_args = list(fun.args)

    # FIXME: Check whether there are functions in fun and sep
    functions = 1
    for s in sep:
        functions = functions * s
        for j in range(0, len(s.args)):
            subs_args.append(s.args[j])

    # Check whether variables match
    if len(subs_args) != len(orig_args):
        raise ValueError("Variable counts do not match")
    # Check for duplicate arguments like  [X(x), u(x, y)]
    if len(subs_args) != len(set(subs_args)):
        raise ValueError("Duplicate substitution arguments detected")
    # Check whether the variables match
    if set(orig_args) != set(subs_args):
        raise ValueError("Arguments do not match")

    # Substitute original function with separated...
    eq = eq.lhs.subs(fun, functions)

    result = 0
    for i in eq.args:
        result += i/functions

    svar = subs_args[0]
    dvar = subs_args[1:]
    return _separate(result, svar, dvar)

def _separate(eq, dep, others):
    """Separate expression into two parts based on dependencies of variables."""

    # FIRST PASS
    # Extract derivatives depending our separable variable...
    terms = set()
    for term in eq.args:
        if term.is_Mul:
            for i in term.args:
                if i.is_Derivative and not i.has_any_symbols(*others):
                    terms.add(term)
                    continue
        elif term.is_Derivative and not term.has_any_symbols(*others):
            terms.add(term)
    # Find the factor that we need to divide by
    div = set()
    for term in terms:
        ext, sep = term.expand().as_independent(dep)
        # Failed?
        if sep.has_any_symbols(*others):
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
        if not term.has_any_symbols(*others):
            lhs += term
            continue
        # ...otherwise, try to separate
        temp, sep = term.expand().as_independent(dep)
        # Failed?
        if sep.has_any_symbols(*others):
            return None
        # Extract the divisors
        div.add(sep)
        rhs -= term.expand()
    # Do the division
    fulldiv = 0
    for d in div:
        fulldiv += d
    lhs = simplify(lhs/fulldiv).expand()
    rhs = simplify(rhs/fulldiv).expand()
    # ...and check whether we were successful :)
    if lhs.has_any_symbols(*others) or rhs.has_any_symbols(dep):
        return None
    return [lhs, rhs]
