from sympy import Function, Symbol, diff, symbols
from sympy.core.compatibility import reduce


def dynamicsymbols(names, level=0):
    """Uses symbols and Function for functions of time.

    Creates a SymPy UndefinedFunction, which is then initialized as a function
    of a variable, the default being Symbol('t').

    Parameters
    ==========

    names : str
        Names of the dynamic symbols you want to create; works the same way as
        inputs to symbols
    level : int
        Level of differentiation of the returned function; d/dt once of t,
        twice of t, etc.

    Examples
    ========

    >>> from sympy.physics.vector import dynamicsymbols
    >>> from sympy import diff, Symbol
    >>> q1 = dynamicsymbols('q1')
    >>> q1
    q1(t)
    >>> diff(q1, Symbol('t'))
    Derivative(q1(t), t)

    """

    esses = symbols(names, cls=Function)
    t = dynamicsymbols._t
    if hasattr(esses, '__iter__'):
        esses = [reduce(diff, [t]*level, e(t)) for e in esses]
        return esses
    else:
        return reduce(diff, [t]*level, esses(t))


dynamicsymbols._t = Symbol('t')
dynamicsymbols._str = '\''
