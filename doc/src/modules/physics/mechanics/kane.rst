==================================
Kane's Method in Physics/Mechanics
==================================

:mod:`mechanics` has been written for use with Kane's method of forming
equations of motion [Kane1985]_. This document will describe Kane's Method
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
  >>> (fr, frstar) = KM.kanes_equations(FL, BL)
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

The forcing vector can be linearized as well; its Jacobian is taken only with
respect to the independent coordinates and speeds. The linearized forcing
vector is of size (n + o) x (n - l + o - m), where l is the number of
configuration constraints and m is the number of motion constraints. Two
matrices are returned; the first is an "A" matrix, or the Jacobian with respect
to the independent states, the second is a "B" matrix, or the Jacobian with
respect to 'forces'; this can be an empty matrix if there are no 'forces'.
Forces here are undefined functions of time (dynamic symbols); they are only
allowed to be in the forcing vector and their derivatives are not allowed to be
present. If dynamic symbols appear in the mass matrix or kinematic differential
equations, an error with be raised. ::

  >>> KM.linearize()[0]
  Matrix([
  [0, 1],
  [0, 0]])

Exploration of the provided examples is encouraged in order to gain more
understanding of the ``KanesMethod`` object.

======================================
Lagrange's Method in Physics/Mechanics
======================================

Structure of Equations
======================

In :mod:`mechanics` we are assuming there are 3 basic sets of equations needed
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

The formulation of the equations of motion in :mod:`mechanics` using
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
as keyword arguments in a list of expressions which are equal to zero. It
should be noted that :mod:`mechanics` requires that the holonomic constraint
equations must be supplied as velocity level constraint equations i.e. the
holonomic constraint equations must be supplied after they have been
differentiated with respect to time. Modifying the example above, the equations
of motion can then be generated: ::

  >>> LM = LagrangesMethod(L, [q1, q2], coneqs = [q1d - q2d])

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
  >>> LM = LagrangesMethod(L, [q1, q2], forcelist = FL, frame = N)
  >>> LM.form_lagranges_equations()
  Matrix([
  [2*q1'' - 7],
  [    2*q2'']])

Exploration of the provided examples is encouraged in order to gain more
understanding of the ``LagrangesMethod`` object.
