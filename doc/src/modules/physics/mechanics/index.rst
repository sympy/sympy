===================
Classical Mechanics
===================

.. topic:: Abstract

   In this documentation many components of the physics/mechanics module will
   be discussed. :mod:`mechanics` has been written to allow for creation of
   symbolic equations of motion for complicated multibody systems.

Vector
======

This module derives the vector-related abilities and related functionalities
from :mod:`physics.vector`. Please have a look at the documentation of
:mod:`physics.vector` and its necessary API to understand the vector capabilities
of :mod:`mechanics`.

Mechanics
=========

In physics, mechanics describes conditions of rest (statics) or motion
(dynamics). There are a few common steps to all mechanics problems. First, an
idealized representation of a system is described. Next, we use physical laws
to generate equations that define the system's behavior. Then, we solve these
equations, sometimes analytically but usually numerically. Finally, we extract
information from these equations and solutions. The current scope of the module
is multi-body dynamics: the motion of systems of multiple particles and/or
rigid bodies. For example, this module could be used to understand the motion
of a double pendulum, planets, robotic manipulators, bicycles, and any
other system of rigid bodies that may fascinate us.

Often, the objective in multi-body dynamics is to obtain the trajectory of a
system of rigid bodies through time. The challenge for this task is to first
formulate the equations of motion of the system. Once they are formulated, they
must be solved, that is, integrated forward in time. When digital computers
came around, solving became the easy part of the problem. Now, we can
tackle more complicated problems, which leaves the challenge of formulating the
equations.

The term "equations of motion" is used to describe the application of Newton's
second law to multi-body systems. The form of the equations of motion depends
on the method used to generate them. This package implements two of these
methods: Kane's method and Lagrange's method. This module facilitates the
formulation of equations of motion, which can then be solved (integrated) using
generic ordinary differential equation (ODE) solvers.

The approach to a particular class of dynamics problems, that of forward
dynamics, has the following steps:

1.  describing the system's geometry and configuration,
2.  specifying the way the system can move, including constraints on its motion
3.  describing the external forces and moments on the system,
4.  combining the above information according to Newton's second law
    (:math:`\mathbf{F}=m\mathbf{a}`), and
5.  organizing the resulting equations so that they can be integrated to obtain
    the system's trajectory through time.

Together with the rest of :mod:`SymPy`, this module performs steps 4 and 5,
provided that the user can perform 1 through 3 for the module. That is to say,
the user must provide a complete representation of the free
body diagrams that themselves represent the system, with which this code can
provide equations of motion in a form amenable to numerical integration. Step
5 above amounts to arduous algebra for even fairly simple multi-body systems.
Thus, it is desirable to use a symbolic math package, such as Sympy, to
perform this step. It is for this reason that this module is a part of Sympy.
Step 4 amounts to this specific module, sympy.physics.mechanics.


Guide to Mechanics
==================

.. toctree::
    :maxdepth: 2

    masses.rst
    kane.rst
    lagrange.rst
    linearize.rst
    examples.rst
    advanced.rst
    reference.rst

Mechanics API
=============

.. toctree::
    :maxdepth: 2

    api/part_bod.rst
    api/kane_lagrange.rst
    api/linearize.rst
    api/expr_manip.rst
    api/printing.rst
