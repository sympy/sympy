Holonomic Functions
*******************

This text aims to explain what a Holonomic Function is and our
implemetations on it. We will also mention its uses in our scope.

Overview
========
The holonomic module is intended to deal with Holonomic Functions along
with various operations on them like addition, multiplication, composition,
integration and differentiation. The module also implements various kinds of
conversions like converting a Holonomic Function to an other form and the
other way around.

Definition of Holonomic Function
---------------------------------

Holonomic Function is a very general type of function and includes a lot of
known functions in it. In fact more known Hypergeometric Functon and
Meijer G-function are also a special case of it.

The mathematical definition of this function is pretty simple:

A Holonomic Function is a solution to an ordinary differential equation having
polynomial coefficients only. In this module the differential equation will be
represented by a Differential Operator annihilating the function.

Since solution of a Differential equation often consists of a family of
functions rather than a unique function it is generally advisable to give
initial conditions also when defining a Holonomic Function.

Take :math:`sin(x)` for instance, differential equation satisfied by the
function is :math:`y''(x) + y(x) = 0`. The general solution of this ODE
is :math:`C_{1}sin(x) + C_{2}cos(x)` but to get :math:`sin(x)` we need to
provide initial conditions on the function i.e. :math:`y(0) = 0, y'(0) = 1`.

To represent the same in this module you need to give the annihilator.
Basically a differential operator is an operator on functions that
differentiates them. So :math:`D^{n}.y(x) = y^{n}(x)` where :math:`y^{n}(x)`
denotes ``n`` times differentiation of :math:`y(x)` with respect to ``x``.

So the differential equation can also be written as
:math:`D^{2}y(x) + y(x) = 0` or :math:`(D^{2} + 1)y(x) = 0`.
The part left of :math:`y(x)` is called the annihilator i.e. :math:`D^{2}+1`.

So this is how one will represent ``sin(x)`` as a Holonomic Function:

    >>> from sympy.holonomic import *
    >>> from sympy.abc import x
    >>> from sympy import ZZ
    >>> R, D = DifferentialOperators(ZZ.old_poly_ring(x), 'D')
    >>> HolonomicFunction(D**2 + 1, x, 0, [0, 1])

