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
- Variable : A symbol that can be used as a variable along the length. By default,
  this is set to ``Symbol(x)``.
- Boundary Conditions
    - bc_moment : Boundary conditions for moment.
    - bc_slope : Boundary conditions for slope.
    - bc_deflection : Boundary conditions for deflection.
- Load Distribution

We have following methods under the beam class:

- apply_load
- shear_force
- bending_moment
- slope
- deflection


Examples
========

Let us solve some beam bending problems using this module :

**Example 1**

**Example 2**

**Example 3**
