from sympy.simplify import simplify

def is_even(f, x):
    """
    Determines if a function is even or not.

    **Usage**

        is_even(function, variable)

    **Details**

        ``function`` can be any function defined in sympy or by the user
        ``variable`` is the variable with respect to which evenness has to be determined

    Examples
    ========

    >>> from sympy.functions import *
    >>> from sympy import is_even
    >>> from sympy.core.symbols import symbols
    >>> x = symbols('x')
    >>> is_even(sin(x), x)
    False
    >>> is_even(cos(x), x)
    True
    >>> is_even(x**2, x)
    True

    References
    ==========

    http://en.wikipedia.org/wiki/Even_and_odd_functions
    """
    
    return simplify(f - f.subs(x, -x)) == 0

def is_odd(f, x):
    """
    Determines if a function is odd or not.
    A function f is odd iff f(-x) = -f(x).

    **Usage**

        is_odd(function, variable)

    **Details**

        ``function`` can be any function defined in sympy or by the user
        ``variable`` is the variable with respect to which oddness has to be determined

    Examples
    ========

    >>> from sympy.functions import *
    >>> from sympy import is_odd
    >>> from sympy.core.symbols import symbols
    >>> x = symbols('x')
    >>> is_even(sin(x), x)
    True
    >>> is_even(cos(x), x)
    False
    >>> is_even(x**3, x)
    True

    References
    ==========

    http://en.wikipedia.org/wiki/Even_and_odd_functions
    """
    return simplify(f.subs(x, -x) + f) == 0
