Basic Usage
-----------

The nsolve function is used to numerically solve equations or systems of equations::

    >>> from sympy import nsolve, Symbol, sin, cos
    >>> x = Symbol('x')
    >>> nsolve(sin(x), x, 2)
    3.14159265358979
    >>> nsolve(sin(x), x, 0)
    0.0

Solving Single Equations
------------------------

For single equations, nsolve can find roots numerically::

    >>> from sympy import exp, log
    >>> nsolve(exp(x) - 1, x, 0)
    0.0
    >>> nsolve(log(x) + x, x, 2)
    0.567143290409784

Systems of Equations
-------------------

For systems of equations, provide a list of equations and a list of symbols::

    >>> y = Symbol('y')
    >>> eq1 = x**2 + y**2 - 1
    >>> eq2 = x - y
    >>> nsolve((eq1, eq2), (x, y), (1, 1))
    Matrix([
    [0.707106781186548],
    [0.707106781186547]])

Solving Transcendental Equations
--------------------------------

nsolve is particularly useful for transcendental equations::

    >>> nsolve(sin(x) - x/2, x, 2)
    1.89549426703398 