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

In Sympy, we can constuct a 2D beam objects with the following properties :

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
- evaluate_reaction_forces
- shear_force
- bending_moment
- slope
- deflection


Examples
========

Let us solve some beam bending problems using this module :

**Example 1**

A beam of length 9 meters is having a fixed support at the start.
A distributed constant load of 8 kN/m is applied downward from the starting
point till 5 meters away from the start. A clockwise moment of 50 kN.m is
applied at 5 meters away from the start of the beam. A downward point load
of 12 kN is applied at the end.

.. note::

    Since a user is free to choose its own sign convention we are considering
    the upward forces and clockwise bending moment being positive.


>>> from sympy.physics.continuum_mechanics.beam import Beam
>>> from sympy import symbols
>>> E, I = symbols('E, I')
>>> R1, M1 = symbols('R1, M1')
>>> b = Beam(9, E, I)
>>> b.apply_load(R1, 0, -1)
>>> b.apply_load(M1, 0, -2)
>>> b.apply_load(-8, 0, 0, end=5)
>>> b.apply_load(50, 5, -2)
>>> b.apply_load(-12, 9, -1)
>>> b.bc_slope = [(0, 0)]
>>> b.bc_deflection = [(0, 0)]
>>> b.evaluate_reaction_forces(R1, M1)
>>> b.reaction_forces
{M₁: -258, R₁: 52}
>>> b.load
         -2         -1        0             -2            0             -1
- 258⋅<x>   + 52⋅<x>   - 8⋅<x>  + 50⋅<x - 5>   + 8⋅<x - 5>  - 12⋅<x - 9>  
>>> b.shear_force()
         -1         0        1             -1            1             0
- 258⋅<x>   + 52⋅<x>  - 8⋅<x>  + 50⋅<x - 5>   + 8⋅<x - 5>  - 12⋅<x - 9> 
>>> b.bending_moment()
         0         1        2             0            2             1
- 258⋅<x>  + 52⋅<x>  - 4⋅<x>  + 50⋅<x - 5>  + 4⋅<x - 5>  - 12⋅<x - 9> 
>>> b.slope()
                            3                          3             
         1         2   4⋅<x>              1   4⋅<x - 5>             2
- 258⋅<x>  + 26⋅<x>  - ────── + 50⋅<x - 5>  + ────────── - 6⋅<x - 9> 
                         3                        3                  
─────────────────────────────────────────────────────────────────────
                                 E⋅I                                 
>>> b.deflection()
                   3      4                        4             
         2   26⋅<x>    <x>              2   <x - 5>             3
- 129⋅<x>  + ─────── - ──── + 25⋅<x - 5>  + ──────── - 2⋅<x - 9> 
                3       3                      3                 
─────────────────────────────────────────────────────────────────
                               E⋅I                               

**Example 2**

**Example 3**
