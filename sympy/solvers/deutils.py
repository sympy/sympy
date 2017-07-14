"""Utility functions for classifying and solving
ordinary and partial differential equations.

Contains
========
_preprocess
ode_order
_desolve

"""
from __future__ import print_function, division

from sympy.core.function import Derivative, AppliedUndef
from sympy.core.relational import Equality
from sympy.core.symbol import Wild

def _find_func(expr):
    """
    Autodetect the function to be solved for in expr.
    """
    derivs = expr.atoms(Derivative)
    funcs = set().union(*[d.atoms(AppliedUndef) for d in derivs])
    if len(funcs) != 1:
        raise ValueError('The function cannot be '
            'automatically detected for %s.' % expr)
    func, = funcs
    return func


def _preprocess(expr, func):
    """Prepare expr for solving by making sure that differentiation
    is done so that only func remains in unevaluated derivatives.

    >>> from sympy.solvers.deutils import _preprocess
    >>> from sympy import Derivative, Function, Integral, sin
    >>> from sympy.abc import x, y, z
    >>> f, g = map(Function, 'fg')

    Apply doit to derivatives that contain more than the function
    of interest:

    >>> _preprocess(Derivative(f(x) + x, x), f(x))
    Derivative(f(x), x) + 1

    Do others if the differentiation variable(s) intersect with those
    of the function of interest or contain the function of interest:

    >>> _preprocess(Derivative(g(x), y, z), f(y))
    0
    >>> _preprocess(Derivative(f(y), z), f(y))
    0
    """
    derivs = expr.atoms(Derivative)
    fvars = set(func.args)
    reps = [(d, d.doit()) for d in derivs if d.has(func) or
            set(d.variables) & fvars]
    return expr.subs(reps)


def ode_order(expr, func):
    """
    Returns the order of a given differential
    equation with respect to func.

    This function is implemented recursively.

    Examples
    ========

    >>> from sympy import Function
    >>> from sympy.solvers.deutils import ode_order
    >>> from sympy.abc import x
    >>> f, g = map(Function, ['f', 'g'])
    >>> ode_order(f(x).diff(x, 2) + f(x).diff(x)**2 +
    ... f(x).diff(x), f(x))
    2
    >>> ode_order(f(x).diff(x, 2) + g(x).diff(x, 3), f(x))
    2
    >>> ode_order(f(x).diff(x, 2) + g(x).diff(x, 3), g(x))
    3

    """
    a = Wild('a', exclude=[func])
    if expr.match(a):
        return 0

    if isinstance(expr, Derivative):
        if expr.args[0] == func:
            return len(expr.variables)
        else:
            order = 0
            for arg in expr.args[0].args:
                order = max(order, ode_order(arg, func) + len(expr.variables))
            return order
    else:
        order = 0
        for arg in expr.args:
            order = max(order, ode_order(arg, func))
        return order

def _desolve(eq, func, hint="default", ics=None, simplify=True,
    classifier=None, **kwargs):
    """This is a helper function to dsolve and pdsolve in the ode
    and pde modules.

    If the hint provided to the function is "default", then a dict with
    the following keys are returned

    'default' - The default key as returned by classifier functions in ode
                and pde.py

    'hint'    - The hint given by the user for which the differential equation
                is to be solved. If the hint given by the user is 'default',
                then the value of 'hint' and 'default' is the same.

    'order'   - The order of the function as returned by ode_order

    'match'   - It returns the match as given by the classifier functions, for
                the default hint.

    If the hint provided to the function is not "default" and is not in
    ('all', 'all_Integral', 'best'), then a dict with the above mentioned keys
    is returned along with the keys which are returned when dict in
    classify_ode or classify_pde is set True

    If the hint given is in ('all', 'all_Integral', 'best'), then this function
    returns a nested dict, with the keys, being the set of classified hints
    returned by classifier functions, and the values being the dict of form
    as mentioned above.

    See Also
    ========
    classify_ode(ode.py)
    classify_pde(pde.py)
    """
    xi = kwargs.get('xi')
    eta = kwargs.get('eta')
    x0 = kwargs.get('x0', 0)
    terms = kwargs.get('n')

    hints = classifier(eq, func, dict=True, ics=ics, xi=xi, eta=eta,
        n=terms, x0=x0, prep=False)

    if hints['order'] == 0:
        raise ValueError(
            str(eq) + " is not a differential equation in " + str(func))

    if not hints['default']:
        # classify_ode will set hints['default'] to None if no hints match.
        if hint not in classifier.allhints and hint != 'default':
            raise ValueError("Hint not recognized: " + hint)
        elif hint not in hints['ordered_hints'] and hint != 'default':
            raise ValueError(
                "%s %s does not match hint %s" % (classifier.kind, eq, hint))
        else:
            raise NotImplementedError(
                "%s(): Cannot solve %s" % (classifier.solve_func, eq))

    if hint == 'default':
        hint = hints['default']
        hints = _filter_hints(hints, hint)
    elif hint in ('all', 'all_Integral', 'best'):
        retdict = {}
        gethints = set(hints) - set(['order', 'default', 'ordered_hints'])
        if hint == 'all_Integral':
            for i in hints:
                if i.endswith('_Integral'):
                    gethints.remove(i[:-len('_Integral')])
            # special cases
            for k in ["1st_homogeneous_coeff_best", "1st_power_series",
                    "lie_group", "2nd_power_series_ordinary",
                    "2nd_power_series_regular"]:
                if k in gethints:
                    gethints.remove(k)
        for i in gethints:
            newhints = _filter_hints(hints, i)
            newhints['hint'] = i
            retdict[i] = newhints
        retdict['all'] = True
        return retdict
    elif hint not in classifier.allhints:
        raise ValueError("Hint not recognized: " + hint)
    elif hint not in hints:
        raise ValueError(
            "%s %s does not match hint %s" % (classifier.kind, eq, hint))
    # Key added to identify the hint needed to solve the equation
    hints['hint'] = hint
    return hints

def _filter_hints(hints, hint):
    return {
        'default': hint, hint: hints[hint], 'order': hints['order']}
