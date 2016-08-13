About Holonomic Function
========================

This text aims to explain Holonomic Functions. We assume you
have a basic idea of Differential equations and Abstract algebra.

Definition
----------

Holonomic Function is a very general type of special function that includes
lots of simple known functions as its special cases. In fact the more known
Hypergeometric Function and Meijer G-function are also a special case of it.

A function is called Holonomic if it's a solution to an ordinary differential
equation having polynomial coefficients only.
Since the general solution of a Differential equation consists of a family of
functions rather than a single function, Holonomic Functions are usually defined
by a set of initial conditions along with the differential equation.

Let :math:`K` be a field of characteristic ``0``. For example, :math:`K` can be
``QQ``, ``ZZ`` or ``RR``.
A function :math:`f(x)` will be holonomic if there exists polynomials
:math:`p_0, p_1, p_2, ... p_r \in K[x]` such that

.. math::

    p_0 \cdot f(x) + p_1 \cdot f^{(1)}(x) + p_2 \cdot f^{(2)}(x) + ... + p_r \cdot f^{(r)}(x) = 0

This differential equation can also be written as :math:`L \cdot f(x) = 0` where

.. math::

    L = p_0 + p_1 \cdot D + p_2 \cdot D^2 + ... p_r \cdot D^r

Here ``D`` is the Differential Operator and :math:`L` is called the annihilator of the function.

A unique holonomic function can be defined from the annihilator and a set of initial conditons.
For instance:

.. math::

    f(x) = exp(x): L = D - 1,\: f(0) = 1

Other fundamental functions such as :math:`sin(x), cos(x), log(x)` etc. are also holonomic.
