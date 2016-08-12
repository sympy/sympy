Holonomic Functions
*******************

.. module:: sympy.holonomic

This text aims to explain what a Holonomic Function is, its significance and SymPy's implemetations
on it.

Overview
========

A Holonomic Function is a solution to an ordinary differential equation with polynomial coefficients.
In this module the differential equation will be represented by a Differential Operator, also called
annihilator of the function.
Since solution of a differential equation often consists of a family of functions rather than a
unique function it is generally advisable to give initial conditions when defining a Holonomic Function.

.. autoclass:: sympy.holonomic.holonomic.HolonomicFunction
    :members:
