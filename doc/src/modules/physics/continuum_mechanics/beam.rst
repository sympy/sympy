=================
Beam
====

The `Beam` class represents a prismatic beam that can have loads,
supports, and boundary conditions applied to it. It is used to perform
structural analysis such as shear force and bending moment diagrams.

Example
-------

.. code-block:: python

    from sympy.physics.continuum_mechanics import Beam
    from sympy import symbols

    E, I = symbols('E I')
    b = Beam(3, E, I)
    b.apply_load(10, 1, -1)  # Point load of 10N at x=1m
    b.solve_for_sole()
=================

.. automodule:: sympy.physics.continuum_mechanics.beam
   :members:
