"""
Analytical methods for solving Partial Differential Equations
Currently implemented methods:
    - separation of variables - pde_separate

"""

from sympy import Derivative, diff, Eq, Equality, Mul
from sympy.simplify import simplify
import operator

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
    if len(subs_args) != len(set(subs_args)):
        raise ValueError("Duplicate substitution arguments detected")
    # Check whether the variables match
    if set(orig_args) != set(subs_args):
        raise ValueError("Arguments do not match")

    # Substitute original function with separated...
    result = eq.lhs.subs(fun, functions)

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

    Examples:

    >>> from sympy import *
    >>> from sympy.abc import x, t
    >>> u, X, T = map(Function, 'uXT')

    >>> eq = Eq(Derivative(u(x, t), x), E**(u(x, t))*Derivative(u(x, t), t))
    >>> pde_separate_add(eq, u(x, t), [X(x), T(t)])
    [D(X(x), x)*exp(-X(x)), D(T(t), t)*exp(T(t))]

    """
    return pde_separate(eq, fun, sep, strategy='add')

def pde_separate_mul(eq, fun, sep):
    """
    Helper function for searching multiplicative separable solutions.

    Consider an equation of two independent variables x, y and a dependent
    variable w, we look for the product of two functions depending on different
    arguments:

    `w(x, y, z) = X(x)*u(y, z)`

    Examples:

    >>> from sympy import *
    >>> from sympy.abc import x, y
    >>> u, X, Y = map(Function, 'uXY')

    >>> eq = Eq(Derivative(u(x, y), x, 2), Derivative(u(x, y), y, 2))
    >>> pde_separate_mul(eq, u(x, y), [X(x), Y(y)])
    [D(X(x), x, x)/X(x), D(Y(y), y, y)/Y(y)]

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
    fulldiv = reduce(operator.add, div)
    lhs = simplify(lhs/fulldiv).expand()
    rhs = simplify(rhs/fulldiv).expand()
    # ...and check whether we were successful :)
    if lhs.has_any_symbols(*others) or rhs.has_any_symbols(dep):
        return None
    return [lhs, rhs]
