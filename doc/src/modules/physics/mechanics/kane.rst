==================================
Kane's Method in Physics/Mechanics
==================================

:mod:`mechanics` provides functionality for deriving equations of motion
using Kane's method [Kane1985]_. This document will describe Kane's method
as used in this module, but not how the equations are actually derived.

Structure of Equations
======================

In :mod:`mechanics` we are assuming there are 5 basic sets of equations needed
to describe a system. They are: holonomic constraints, non-holonomic
constraints, kinematic differential equations, dynamic equations, and
differentiated non-holonomic equations.

.. math::
  \mathbf{f_h}(q, t) &= 0\\
  \mathbf{k_{nh}}(q, t) u + \mathbf{f_{nh}}(q, t) &= 0\\
  \mathbf{k_{k\dot{q}}}(q, t) \dot{q} + \mathbf{k_{ku}}(q, t) u +
  \mathbf{f_k}(q, t) &= 0\\
  \mathbf{k_d}(q, t) \dot{u} + \mathbf{f_d}(q, \dot{q}, u, t) &= 0\\
  \mathbf{k_{dnh}}(q, t) \dot{u} + \mathbf{f_{dnh}}(q, \dot{q}, u, t) &= 0\\

In :mod:`mechanics` holonomic constraints are only used for the linearization
process; it is assumed that they will be too complicated to solve for the
dependent coordinate(s).  If you are able to easily solve a holonomic
constraint, you should consider redefining your problem in terms of a smaller
set of coordinates. Alternatively, the time-differentiated holonomic
constraints can be supplied.

Kane's method forms two expressions, :math:`F_r` and :math:`F_r^*`, whose sum
is zero. In this module, these expressions are rearranged into the following
form:

 :math:`\mathbf{M}(q, t) \dot{u} = \mathbf{f}(q, \dot{q}, u, t)`

For a non-holonomic system with `o` total speeds and `m` motion constraints, we
will get o - m equations. The mass-matrix/forcing equations are then augmented
in the following fashion:

.. math::
  \mathbf{M}(q, t) &= \begin{bmatrix} \mathbf{k_d}(q, t) \\
  \mathbf{k_{dnh}}(q, t) \end{bmatrix}\\
  \mathbf{_{(forcing)}}(q, \dot{q}, u, t) &= \begin{bmatrix}
  - \mathbf{f_d}(q, \dot{q}, u, t) \\ - \mathbf{f_{dnh}}(q, \dot{q}, u, t)
  \end{bmatrix}\\


Kane's Method in Physics/Mechanics
==================================

The formulation of the equations of motion in :mod:`mechanics` starts with
creation of a ``KanesMethod`` object. Upon initialization of the
``KanesMethod`` object, an inertial reference frame needs to be supplied. along
with some basic system information, suchs as coordinates and speeds ::

  >>> from sympy.physics.mechanics import *
  >>> N = ReferenceFrame('N')
  >>> q1, q2, u1, u2 = dynamicsymbols('q1 q2 u1 u2')
  >>> q1d, q2d, u1d, u2d = dynamicsymbols('q1 q2 u1 u2', 1)
  >>> KM = KanesMethod(N, [q1, q2], [u1, u2])

It is also important to supply the order of coordinates and speeds properly if
there are dependent coordinates and speeds. They must be supplied after
independent coordinates and speeds or as a keyword argument; this is shown
later. ::

  >>> q1, q2, q3, q4 = dynamicsymbols('q1 q2 q3 q4')
  >>> u1, u2, u3, u4 = dynamicsymbols('u1 u2 u3 u4')
  >>> # Here we will assume q2 is dependent, and u2 and u3 are dependent
  >>> # We need the constraint equations to enter them though
  >>> KM = KanesMethod(N, [q1, q3, q4], [u1, u4])

Additionally, if there are auxiliary speeds, they need to be identified here.
See the examples for more information on this. In this example u4 is the
auxiliary speed. ::

  >>> KM = KanesMethod(N, [q1, q3, q4], [u1, u2, u3], u_auxiliary=[u4])

Kinematic differential equations must also be supplied; there are to be
provided as a list of expressions which are each equal to zero. A trivial
example follows: ::

  >>> kd = [q1d - u1, q2d - u2]

Turning on ``mechanics_printing()`` makes the expressions significantly
shorter and is recommended. Alternatively, the ``mprint`` and ``mpprint``
commands can be used.

If there are non-holonomic constraints, dependent speeds need to be specified
(and so do dependent coordinates, but they only come into play when linearizing
the system). The constraints need to be supplied in a list of expressions which
are equal to zero, trivial motion and configuration constraints are shown
below: ::

  >>> N = ReferenceFrame('N')
  >>> q1, q2, q3, q4 = dynamicsymbols('q1 q2 q3 q4')
  >>> q1d, q2d, q3d, q4d = dynamicsymbols('q1 q2 q3 q4', 1)
  >>> u1, u2, u3, u4 = dynamicsymbols('u1 u2 u3 u4')
  >>> #Here we will assume q2 is dependent, and u2 and u3 are dependent
  >>> speed_cons = [u2 - u1, u3 - u1 - u4]
  >>> coord_cons = [q2 - q1]
  >>> q_ind = [q1, q3, q4]
  >>> q_dep = [q2]
  >>> u_ind = [u1, u4]
  >>> u_dep = [u2, u3]
  >>> kd = [q1d - u1, q2d - u2, q3d - u3, q4d - u4]
  >>> KM = KanesMethod(N, q_ind, u_ind, kd,
  ...           q_dependent=q_dep,
  ...           configuration_constraints=coord_cons,
  ...           u_dependent=u_dep,
  ...           velocity_constraints=speed_cons)

A dictionary returning the solved :math:`\dot{q}`'s can also be solved for: ::

  >>> mechanics_printing(pretty_print=False)
  >>> KM.kindiffdict()
  {q1': u1, q2': u2, q3': u3, q4': u4}

The final step in forming the equations of motion is supplying a list of
bodies and particles, and a list of 2-tuples of the form ``(Point, Vector)``
or ``(ReferenceFrame, Vector)`` to represent applied forces and torques. ::

  >>> N = ReferenceFrame('N')
  >>> q, u = dynamicsymbols('q u')
  >>> qd, ud = dynamicsymbols('q u', 1)
  >>> P = Point('P')
  >>> P.set_vel(N, u * N.x)
  >>> Pa = Particle('Pa', P, 5)
  >>> BL = [Pa]
  >>> FL = [(P, 7 * N.x)]
  >>> KM = KanesMethod(N, [q], [u], [qd - u])
  >>> (fr, frstar) = KM.kanes_equations(BL, FL)
  >>> KM.mass_matrix
  Matrix([[5]])
  >>> KM.forcing
  Matrix([[7]])

When there are motion constraints, the mass matrix is augmented by the
:math:`k_{dnh}(q, t)` matrix, and the forcing vector by the :math:`f_{dnh}(q,
\dot{q}, u, t)` vector.

There are also the "full" mass matrix and "full" forcing vector terms, these
include the kinematic differential equations; the mass matrix is of size (n +
o) x (n + o), or square and the size of all coordinates and speeds. ::

  >>> KM.mass_matrix_full
  Matrix([
  [1, 0],
  [0, 5]])
  >>> KM.forcing_full
  Matrix([
  [u],
  [7]])

Exploration of the provided examples is encouraged in order to gain more
understanding of the ``KanesMethod`` object.
