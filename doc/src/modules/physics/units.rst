.. _physics-units:

=====
Units
=====

Introduction
============

This module provides around 200 predefined units that are commonly used in the
sciences.  Additionally, it provides the :class:`Unit` class which allows you
to define your own units.

Examples
========

All examples in this tutorial are computable, so one can just copy and
paste them into a Python shell and do something useful with them. All
computations were done using the following setup::

    >>> from sympy.physics.units import *

Dimensionless quantities
------------------------

Provides variables for commonly written dimensionless (unit-less) quantities.

    >>> 2*ten
    20
    >>> 20*percent
    1/5
    >>> 300*kilo*20*percent
    60000
    >>> nano*deg
    pi/180000000000

Base units
----------

The SI base units are defined variable name that are commonly used in written
and verbal communication.  The singular abbreviated versions are what is used
for display purposes, but the plural non-abbreviated versions are defined in
case it helps readability.

    >>> 5*meters
    5*m
    >>> milli*kilogram
    kg/1000
    >>> gram
    kg/1000

Note that British Imperial and U.S. customary units are not included.
We strongly urge the use of SI units; only Myanmar (Burma), Liberia, and the
United States have not officially accepted the SI system.


Derived units
-------------

Common SI derived units.

    >>> joule
    kg*m**2/s**2

Docstring
=========

.. automodule:: sympy.physics.units
   :members:
