"""Integration method that uses pattern matching.

The ``pattern_integrate`` function tries to compute an integral by matching a
set of user-supplied patterns to it. The functions ``int_patterns``,
``add_int_pattern`` and ``clear_int_patterns`` can be used to manage the set of
patterns.

"""

from sympy import ask, Symbol
from sympy.core.relational import Relational
from sympy.integrals.integrals import Integral

_patterns = list()

def int_patterns():
    """Return the list of currently defined integration patterns.
    """

    global _patterns
    return _patterns

def add_int_pattern(vars, pat, assumption, res):
    """Add a new pattern to the list of integration patterns.

    The integration variables ``vars`` should be a list of tuples (V, L, U),
    where V is a dummy variable, L is a set of allowed values for the lower
    bound of the domain of integration and U is a set of the allowed values for
    the corresponding upper bound. The pattern ``pat`` should be an expression
    which may contain the dummy variables from ``vars`` plus additional wildcard
    symbols. Make sure to ``exclude`` the integration variables from the
    wildcards to avoid unintentional matching.

    It is possible to further restrict the pattern by giving an assumption which
    involves the wildcard symbols. This assumption can optionally be checked
    when it is decided whether or not the pattern applies.

    The result ``res`` of the integral should be an expression which may involve
    any of the wildcard symbols or integration variables. Wildcards are
    substituted using the mapping that is defined by matching the pattern to the
    actual integrand. If integration variables appear, their respective upper
    and lower bounds are substituted and the difference between the two
    expressions is returned as the integral. Otherwise, it is assumed that the
    given expression already accounts for this difference and it is taken as is.
    """

    global _patterns
    _patterns = _patterns + [(vars, pat, assumption, res)]

def clear_int_patterns():
    """Clears all previously defined integration patterns.
    """

    global _patterns
    _patterns = list()

def pattern_integrate(integral, strict=True):
    r"""Try to solve an integral using pattern matching.

    This integration method is intended for cases where SymPy can do most
    but not all integrals in a large computation and the user is able to solve
    the few remaining ones by hand or can look them up. Then, instead of
    extracting the problematic integrals from a large expression and handling
    them separately, their solutions can be supplied as patterns and the
    integrals will be solved automatically.

    If ``strict`` is ``True``, only apply a pattern if the assumption needed for
    its validity holds. If it is ``None``, apply a pattern also if it is unknown
    or undecided whether the assumption holds. If it is ``False``, ignore
    assumptions for pattern validity and always apply matching patterns.
    """

    global _patterns

    def bounds_allowed(dummies, actuals):
        for dum, act in zip(dummies, actuals):
            for i in 1, 2:
                condition = dum[i].contains(act[i])
                if isinstance(condition, Relational):
                    # condition undecided, assume it is not met
                    return False
                elif not condition:
                    return False
        return True

    for pattern in _patterns:
        if len(integral.limits) < len(pattern[0]):
            # The integral has too few integration variables, so the pattern
            # cannot match.
            continue

        inner_vars = integral.limits[-len(pattern[0]):] # pattern will be matched to these variables
        outer_vars = integral.limits[:-len(pattern[0])]

        if not bounds_allowed(pattern[0], inner_vars):
            continue

        # Generate substitution rules for switching back and forth between
        # dummy variables from patterns and actual integration variables.
        varsubs = [(dum[0], act[0]) for dum, act in zip(pattern[0], inner_vars)]
        inv_varsubs = [(v2, v1) for v1, v2 in varsubs]

        mapping = integral.function.subs(inv_varsubs).match(pattern[1], old=True)
        # Use dummy variables for matching so that exclusions for wildcard
        # symbols are honored.

        if not mapping:
            continue

        # Check assumptions.
        if strict != False and pattern[2] != None:
            assumption_res = ask(pattern[2].subs(mapping))
            if not assumption_res:
                if assumption_res == False or strict == True:
                    continue

        # Pattern matches, use it.
        result = pattern[3].subs(mapping).subs(varsubs)
        for var in inner_vars:
            if var[0] in result.atoms(Symbol):
                result = (result.subs(var[0], var[2]) - result.subs(var[0], var[1]))

        if len(outer_vars):
            return pattern_integrate(Integral(result, *outer_vars), strict=strict)
        else:
            return result

    return integral
