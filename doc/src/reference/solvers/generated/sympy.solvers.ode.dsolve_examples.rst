Basic Usage
-----------

The dsolve function is used to solve ordinary differential equations (ODEs)::

    >>> from sympy import Function, dsolve, Eq, Derivative, sin, cos
    >>> from sympy.abc import x
    >>> f = Function('f')
    >>> dsolve(Derivative(f(x), x, x) + 9*f(x), f(x))
    f(x) = C1*sin(3*x) + C2*cos(3*x)

First Order ODEs
----------------

First-order ODEs can be solved using dsolve::

    >>> dsolve(Derivative(f(x), x) - f(x), f(x))
    f(x) = C1*exp(x)
    
    >>> dsolve(Derivative(f(x), x) + f(x), f(x))
    f(x) = C1*exp(-x)

Initial Value Problems
---------------------

You can provide initial conditions to get a particular solution::

    >>> from sympy import symbols
    >>> C1 = symbols('C1')
    >>> eq = Eq(f(x), C1*exp(x))
    >>> subs = {C1: 5}
    >>> eq.subs(subs)
    Eq(f(x), 5*exp(x))

Solving with a Different Method
-------------------------------

You can specify the method to use when solving the ODE::

    >>> dsolve(Derivative(f(x), x, x) + 9*f(x), f(x), hint='nth_linear_constant_coeff_homogeneous')
    f(x) = C1*sin(3*x) + C2*cos(3*x)

Higher Order Equations
----------------------

Higher-order ODEs can also be solved::

    >>> dsolve(Derivative(f(x), x, x, x) + Derivative(f(x), x), f(x))
    f(x) = C1 + C2*sin(x) + C3*cos(x)
