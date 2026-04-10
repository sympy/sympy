======================================
Lagrange's Method in Physics/Mechanics
======================================

:mod:`sympy.physics.mechanics` provides functionality for deriving equations of motion
using `Lagrange's method <https://en.wikipedia.org/wiki/Lagrangian_mechanics>`_.
This document will describe Lagrange's method as used in this module, but not
how the equations are actually derived.

Structure of Equations
======================

In :mod:`sympy.physics.mechanics` we are assuming there are 3 basic sets of equations needed
to describe a system; the constraint equations, the time differentiated
constraint equations and the dynamic equations.

.. math::
  \mathbf{m_{c}}(q, t) \dot{q} + \mathbf{f_{c}}(q, t) &= 0\\
  \mathbf{m_{dc}}(\dot{q}, q, t) \ddot{q} + \mathbf{f_{dc}}(\dot{q}, q, t) &= 0\\
  \mathbf{m_d}(\dot{q}, q, t) \ddot{q} + \mathbf{\Lambda_c}(q, t)
  \lambda + \mathbf{f_d}(\dot{q}, q, t) &= 0\\

In this module, the expressions formed by using Lagrange's equations of the
second kind are rearranged into the following form:

 :math:`\mathbf{M}(q, t) x = \mathbf{f}(q, \dot{q}, t)`

where in the case of a system without constraints:

 :math:`x = \ddot{q}`

For a constrained system with `n` generalized speeds and `m` constraints, we
will get n - m equations. The mass-matrix/forcing equations are then augmented
in the following fashion:

.. math::
  x = \begin{bmatrix} \ddot{q} \\ \lambda \end{bmatrix} \\
  \mathbf{M}(q, t) &= \begin{bmatrix} \mathbf{m_d}(q, t) &
  \mathbf{\Lambda_c}(q, t) \end{bmatrix}\\
  \mathbf{F}(\dot{q}, q, t) &= \begin{bmatrix} \mathbf{f_d}(q, \dot{q}, t)
  \end{bmatrix}\\


Lagrange's Method in Physics/Mechanics
======================================

The formulation of the equations of motion in :mod:`sympy.physics.mechanics` using
Lagrange's Method starts with the creation of generalized coordinates and a
Lagrangian. The Lagrangian can either be created with the ``Lagrangian``
function or can be a user supplied function. In this case we will supply the
Lagrangian. ::

  >>> from sympy.physics.mechanics import *
  >>> q1, q2 = dynamicsymbols('q1 q2')
  >>> q1d, q2d = dynamicsymbols('q1 q2', 1)
  >>> L = q1d**2 + q2d**2

To formulate the equations of motion we create a ``LagrangesMethod``
object. The Lagrangian and generalized coordinates need to be supplied upon
initialization. ::

  >>> LM = LagrangesMethod(L, [q1, q2])

With that the equations of motion can be formed. ::

  >>> mechanics_printing(pretty_print=False)
  >>> LM.form_lagranges_equations()
  Matrix([
  [2*q1''],
  [2*q2'']])

It is possible to obtain the mass matrix and the forcing vector. ::

  >>> LM.mass_matrix
  Matrix([
  [2, 0],
  [0, 2]])

  >>> LM.forcing
  Matrix([
  [0],
  [0]])

If there are any holonomic or non-holonomic constraints, they must be supplied
as keyword arguments (``hol_coneqs`` and ``nonhol_coneqs`` respectively) in a
list of expressions which are equal to zero. Modifying the example above, the
equations of motion can then be generated: ::

  >>> LM = LagrangesMethod(L, [q1, q2], hol_coneqs=[q1 - q2])

When the equations of motion are generated in this case, the Lagrange
multipliers are introduced; they are represented by ``lam1`` in this case. In
general, there will be as many multipliers as there are constraint equations. ::

  >>> LM.form_lagranges_equations()
  Matrix([
  [ lam1 + 2*q1''],
  [-lam1 + 2*q2'']])

Also in the case of systems with constraints, the 'full' mass matrix is
augmented by the :math:`k_{dc}(q, t)` matrix, and the forcing vector by the
:math:`f_{dc}(q, \dot{q}, t)` vector. The 'full' mass matrix is of size
(2n + o) x (2n + o), i.e. it's a square matrix. ::

  >>> LM.mass_matrix_full
  Matrix([
  [1, 0, 0,  0,  0],
  [0, 1, 0,  0,  0],
  [0, 0, 2,  0, -1],
  [0, 0, 0,  2,  1],
  [0, 0, 1, -1,  0]])
  >>> LM.forcing_full
  Matrix([
  [q1'],
  [q2'],
  [  0],
  [  0],
  [  0]])

If there are any non-conservative forces or moments acting on the system,
they must also be supplied as keyword arguments in a list of 2-tuples of the
form ``(Point, Vector)`` or ``(ReferenceFrame, Vector)`` where the ``Vector``
represents the non-conservative forces and torques. Along with this 2-tuple,
the inertial frame must also be specified as a keyword argument. This is shown
below by modifying the example above: ::

  >>> N = ReferenceFrame('N')
  >>> P = Point('P')
  >>> P.set_vel(N, q1d * N.x)
  >>> FL = [(P, 7 * N.x)]
  >>> LM = LagrangesMethod(L, [q1, q2], forcelist=FL, frame=N)
  >>> LM.form_lagranges_equations()
  Matrix([
  [2*q1'' - 7],
  [    2*q2'']])

Exploration of the provided examples is encouraged in order to gain more
understanding of the ``LagrangesMethod`` object.
