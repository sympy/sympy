"""Integration functions that integrates a sympy expression.

>>> from sympy import integrate, log, cos
>>> from sympy.abc import x
>>> integrate(1/x,x)
log(x)
>>> integrate(sin(x),x)
-cos(x)
"""
from integrals import integrate, Integral, line_integrate
