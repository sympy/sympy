"""Integration functions that integrates a sympy expression.

    Examples
    --------
    >>> from sympy import integrate, sin
    >>> from sympy.abc import x
    >>> integrate(1/x,x)
    log(x)
    >>> integrate(sin(x),x)
    -cos(x)
"""
from integrals import integrate, Integral, line_integrate
