About Holonomic Functions
=========================

This text aims to explain holonomic functions. We assume you
have a basic idea of Differential equations and Abstract algebra.

Definition
----------

Holonomic function is a very general type of special function that includes
lots of simple known functions as its special cases. In fact the more known
hypergeometric function and Meijer G-function are also a special case of it.

A function is called holonomic if it's a solution to an ordinary differential
equation having polynomial coefficients only.
Since the general solution of a differential equation consists of a family of
functions rather than a single function, holonomic functions are usually defined
by a set of initial conditions along with the differential equation.

Let :math:`K` be a field of characteristic ``0``. For example, :math:`K` can be
``QQ`` or ``RR``.
A function :math:`f(x)` will be holonomic if there exists polynomials
:math:`p_0, p_1, p_2, ... p_r \in K[x]` such that

.. math::

    p_0 \cdot f(x) + p_1 \cdot f^{(1)}(x) + p_2 \cdot f^{(2)}(x) + ... + p_r \cdot f^{(r)}(x) = 0

This differential equation can also be written as :math:`L \cdot f(x) = 0` where

.. math::

    L = p_0 + p_1 \cdot D + p_2 \cdot D^2 + ... p_r \cdot D^r

Here `D` is the Differential Operator and `L` is called the annihilator
of the function.

A unique holonomic function can be defined from the annihilator and a set of
initial conditions.
For instance:

.. math::

    f(x) = \exp(x): L = D - 1,\: f(0) = 1

    f(x) = \sin(x): L = D^2 + 1,\: f(0) = 0, f'(0) = 1

Other fundamental functions such as `\cos(x)`, `\log(x)`, bessel functions etc. are also holonomic.

The family of holonomic functions is closed under addition, multiplication, integration,
composition. This means if two functions are given are holonomic, then the
function resulting on applying these operation on them will also be holonomic.

References
----------
https://en.wikipedia.org/wiki/Holonomic_function
