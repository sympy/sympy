"""Utility functions for classifying and solving
ordinary and partial differential equations.

Contains
========
_preprocess
de_order
_desolve

"""
from sympy.core.compatibility import set_union
from sympy.core.function import Function, Derivative, AppliedUndef
from sympy.core.relational import Equality, Eq
from sympy.core.symbol import Wild
#from sympy.solvers.ode import allhints as allhints_ode
#from sympy.solvers.pde import allhints as allhints_pde

def _preprocess(expr, func=None, hint='_Integral'):
    """Prepare expr for solving by making sure that differentiation
    is done so that only func remains in unevaluated derivatives and
    (if hint doesn't end with _Integral) that doit is applied to all
    other derivatives. If hint is None, don't do any differentiation.
    (Currently this may cause some simple differential equations to
    fail.)

    In case func is None, an attempt will be made to autodetect the
    function to be solved for.

    >>> from sympy.solvers.util import _preprocess
    >>> from sympy import Derivative, Function, Integral, sin
    >>> from sympy.abc import x, y, z
    >>> f, g = map(Function, 'fg')

    Apply doit to derivatives that contain more than the function
    of interest:

    >>> _preprocess(Derivative(f(x) + x, x))
    (Derivative(f(x), x) + 1, f(x))

    Do others if the differentiation variable(s) intersect with those
    of the function of interest or contain the function of interest:

    >>> _preprocess(Derivative(g(x), y, z), f(y))
    (0, f(y))
    >>> _preprocess(Derivative(f(y), z), f(y))
    (0, f(y))

    Do others if the hint doesn't end in '_Integral' (the default
    assumes that it does):

    >>> _preprocess(Derivative(g(x), y), f(x))
    (Derivative(g(x), y), f(x))
    >>> _preprocess(Derivative(f(x), y), f(x), hint='')
    (0, f(x))

    Don't do any derivatives if hint is None:

    >>> eq = Derivative(f(x) + 1, x) + Derivative(f(x), y)
    >>> _preprocess(eq, f(x), hint=None)
    (Derivative(f(x) + 1, x) + Derivative(f(x), y), f(x))

    If it's not clear what the function of interest is, it must be given:

    >>> eq = Derivative(f(x) + g(x), x)
    >>> _preprocess(eq, g(x))
    (Derivative(f(x), x) + Derivative(g(x), x), g(x))
    >>> try: _preprocess(eq)
    ... except ValueError: print "A ValueError was raised."
    A ValueError was raised.

    """

    derivs = expr.atoms(Derivative)
    if not func:
        funcs = set_union(*[d.atoms(AppliedUndef) for d in derivs])
        if len(funcs) != 1:
            raise ValueError('The function cannot be '
                'automatically detected for %s.' % expr)
        func = funcs.pop()
    fvars = set(func.args)
    if hint is None:
        return expr, func
    reps = [(d, d.doit()) for d in derivs if not hint.endswith('_Integral') or
            d.has(func) or set(d.variables) & fvars]
    eq = expr.subs(reps)
    return eq, func

def de_order(expr, func):
    """
    Returns the order of a given differential
    equation with respect to func.

    This function is implemented recursively.

    Examples
    ========

    >>> from sympy import Function
    >>> from sympy.solvers.util import de_order
    >>> from sympy.abc import x
    >>> f, g = map(Function, ['f', 'g'])
    >>> de_order(f(x).diff(x, 2) + f(x).diff(x)**2 +
    ... f(x).diff(x), f(x))
    2
    >>> de_order(f(x).diff(x, 2) + g(x).diff(x, 3), f(x))
    2
    >>> de_order(f(x).diff(x, 2) + g(x).diff(x, 3), g(x))
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
                order = max(order, de_order(arg, func) + len(expr.variables))
            return order
    else:
        order = 0
        for arg in expr.args:
            order = max(order, de_order(arg, func))
        return order

def _desolve(eq, func=None, hint="default", simplify=True, **kwargs):
    """This is a helper function to dsolve and pdsolve in the ode
    and pde modules.

    If the hint is given in ('all', 'all_Integral', 'best'), then this function
    returns a nested dict, with the keys, being the
    set of identified hints in classify_ode(). The value of each key is a dict
    consisting of a 'default' key, a 'hint' key which helps in
    identifying the hint needed to be passed to the solver functions
    in ode and pde.py, an 'order' key that identifies the order of the
    differential equation given by 'de_order', a 'func' key that
    tells what function, the differential equation is to be solved for,
    an 'all' key that tells if the hint given is in all or not,
    and the match for the corresponding hint as returned by classify_ode

    If the hint given is not present in ('all', 'all_Integral', 'best')
    then a single dict corresponding to the key value mentioned
    above that has a 'default', 'order', 'func' and 'match' key is
    returned.

    See Also
    ========
    classify_ode(ode.py)
    classify_pde(pde.py)
    """
    prep = kwargs.pop('prep', True)
    if isinstance(eq, Equality):
        eq = eq.lhs - eq.rhs

    # preprocess the equation and find func if not given
    if prep or func is None:
        eq, func = _preprocess(eq, func)
        prep = False

    # type is an argument passed by the solve functions in ode and pde.py
    # that identifies whether the function caller is an ordinary
    # or partial differential equation. Accordingly corresponding
    # changes are made in the function.
    type = kwargs.get('type', None)
    if type == 'ode':
        from sympy.solvers.ode import classify_ode, allhints
        classifier = classify_ode
        string = 'ODE '
        dummy = ''

    elif type == 'pde':
        from sympy.solvers.pde import classify_pde, allhints
        classifier = classify_pde
        allhints = allhints_pde
        string = 'PDE '
        dummy = 'p'

    # Magic that should only be used internally.  Prevents classify_ode from
    # being called more than it needs to be by passing its results through
    # recursive calls.
    if kwargs.get('classify', True):
        hints = classifier(eq, func, dict=True, prep=prep)

    else:
        # Here is what all this means:
        #
        # hint:    The hint method given to _desolve() by the user.
        # hints:   The dictionary of hints that match the DE, along with other
        #          information (including the internal pass-through magic).
        # default: The default hint to return, the first hint from allhints
        #          that matches the hint; obtained from classify_ode().
        # match:   Dictionary containing the match dictionary for each hint
        #          (the parts of the DE for solving).  When going through the
        #          hints in "all", this holds the match string for the current
        #          hint.
        # order:   The order of the DE, as determined by de_order().
        hints = kwargs.get('hint',
                           {'default': hint,
                            hint: kwargs['match'],
                            'order': kwargs['order']})
    if hints['order'] == 0:
        raise ValueError(
            str(eq) + " is not a differential equation in " + str(func))

    if not hints['default']:
        # classify_ode will set hints['default'] to None if no hints match.
        if hint not in allhints and hint != 'default':
            raise ValueError("Hint not recognized: " + hint)
        elif hint not in hints['ordered_hints'] and hint != 'default':
            raise ValueError(string + str(eq) + " does not match hint " + hint)
        else:
            raise NotImplementedError(dummy + "solve" + ": Cannot solve " + str(eq))
    if hint == 'default':
        return _desolve(eq, func, hint=hints['default'], simplify=simplify,
                      prep=prep, classify=False, order=hints['order'],
                      match=hints[hints['default']], type=type)
    elif hint in ('all', 'all_Integral', 'best'):
        retdict = {}
        failedhints = {}
        gethints = set(hints) - set(['order', 'default', 'ordered_hints'])
        if hint == 'all_Integral':
            for i in hints:
                if i.endswith('_Integral'):
                    gethints.remove(i[:-len('_Integral')])
            # special case
            if "1st_homogeneous_coeff_best" in gethints:
                gethints.remove("1st_homogeneous_coeff_best")
        for i in gethints:
            sol = _desolve(eq, func, hint=i, simplify=simplify, prep=prep,
                classify=False, order=hints['order'], match=hints[i], type=type)
            retdict[i] = sol
        retdict['all'] = True
        return retdict
    elif hint not in allhints:  # and hint not in ('default', 'ordered_hints'):
        raise ValueError("Hint not recognized: " + hint)
    elif hint not in hints:
        raise ValueError(string + str(eq) + " does not match hint " + hint)
    else:
        # Key added to identify the hint needed to solve the equation
        hints['hint'] = hint
    hints['func'] = func
    return hints
