Basic Usage
-----------

The pdsolve function is used to solve partial differential equations (PDEs)::

    >>> from sympy import Function, pdsolve, Eq, Derivative
    >>> from sympy.abc import x, y, t
    >>> f = Function('f')
    >>> u = f(x, y)
    >>> ux = u.diff(x)
    >>> uy = u.diff(y)
    >>> eq = Eq(ux + uy, u)
    >>> pdsolve(eq)
    f(x, y) = F(x - y)*exp(x)

Heat Equation
-------------

Solving the 1D heat equation::

    >>> u = Function('u')
    >>> eq = Eq(Derivative(u(x, t), t), Derivative(u(x, t), x, x))
    >>> pdsolve(eq)
    u(x, t) = F(x - 2*sqrt(t)*y)*exp(-y**2) + G(x + 2*sqrt(t)*y)*exp(-y**2)

Wave Equation
-------------

Solving the 1D wave equation::

    >>> u = Function('u')
    >>> eq = Eq(Derivative(u(x, t), t, t), Derivative(u(x, t), x, x))
    >>> pdsolve(eq)
    u(x, t) = F(x - t) + G(x + t)

Laplace Equation
----------------

Solving the 2D Laplace equation::

    >>> u = Function('u')
    >>> eq = Eq(Derivative(u(x, y), x, x) + Derivative(u(x, y), y, y), 0)
    >>> pdsolve(eq)
    u(x, y) = F(x + I*y) + G(x - I*y)

