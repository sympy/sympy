Represnting a Holonomic Function in SymPy
=========================================

Let's understand this with an example:

Take :math:`sin(x)` for instance, the differential equation satisfied by it
is :math:`y''(x) + y(x) = 0`. So by definition we conclude it is a Holonomic
Function. The general solution of this ODE is
:math:`C_{1}sin(x) + C_{2}cos(x)` but to get :math:`sin(x)` we need to
provide initial conditions i.e. :math:`y(0) = 0, y'(0) = 1`.

To represent the same in this module one needs to provide the differential
equation in the form of annihilator. Basically a differential operator is an
operator on functions that differentiates them. So :math:`D^{n}.y(x) = y^{n}(x)`
where :math:`y^{n}(x)` denotes ``n`` times differentiation of :math:`y(x)` with
respect to ``x``.

So the differential equation can also be written as
:math:`D^{2}y(x) + y(x) = 0` or :math:`(D^{2} + 1)y(x) = 0`.
The part left of :math:`y(x)` is called the annihilator i.e. :math:`D^{2}+1`.

So this is how one will represent ``sin(x)`` as a Holonomic Function:

    >>> from sympy.holonomic import *
    >>> from sympy.abc import x
    >>> from sympy import ZZ
    >>> R, D = DifferentialOperators(ZZ.old_poly_ring(x), 'D')
    >>> HolonomicFunction(D**2 + 1, x, 0, [0, 1])
    HolonomicFunction((1) + (1)*D**2, x, 0, [0, 1])

.. module:: sympy.holonomic.holonomic

.. autofunction:: DifferentialOperators

