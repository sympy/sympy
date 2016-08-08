===========================================================
Solving Beam Bending Problems using Discontinuity Functions
===========================================================

To make this document easier to read, we are going to enable pretty printing.

    >>> from sympy import *
    >>> x, y, z = symbols('x y z')
    >>> init_printing(use_unicode=True)

Beam
====

A Beam is a structural element that is capable of withstanding load
primarily by resisting against bending. Beams are characterized by
their second moment of area, their length and their elastic modulus.

In Sympy, we can constuct a beam object with the following properties :

- Length
- Elastic Modulus
- Second Moment of Area
- A symbol that can be used as a variable along the length. By default,
  this is set to ``Symbol(x)``

For example :

    >>> from sympy.physics.continuum_mechanics.beam import Beam
    >>> from sympy import Symbol
    >>> E = Symbol('E')
    >>> I = Symbol('I')
    >>> b = Beam(4, E, I)
