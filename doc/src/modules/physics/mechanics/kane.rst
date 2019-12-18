.. _kane_method:

=============
Kane's Method
=============

The :mod:`~sympy.physics.mechanics.kane` module provides functionality for
deriving equations of motion using Kane's method [Kane1985]_. This document
describes Kane's method as used in this module, please refer to [Kane1985]_ and
the source code for the mathematical details of the underlying algorithms.

Structure of Equations
======================

In :mod:`~sympy.physics.mechanics.kane` we describe a multi-body system with 5
general sets of equations given the:

- :math:`M` : number of holonomic constraint equations
- :math:`n` : number of generalized coordinates and number of generalized speeds, i.e. :math:`n=N-M`
- :math:`m` : number of nonholonomic constraint equations
- :math:`p` : number of independent generalized speeds, i.e. :math:`p=n-m`
- :math:`\mathbf{q}` : vector of independent generalized coordinates where :math:`\mathbf{q} \in \mathbb{R}^n`
- :math:`\mathbf{u}` : vector of generalized speeds where :math:`\mathbf{u} \in \mathbb{R}^n`
- :math:`\mathbf{u}_r` : dependent generalized speeds where :math:`\mathbf{u}_r \in \mathbb{R}^m`
- :math:`\mathbf{u}_s` : independent generalized speeds where :math:`\mathbf{u}_s \in \mathbb{R}^p`

The equations that result from applying Kane's Method are then as follows:

1. Holonomic (configuration) constraints

   .. math::
      \mathbf{f}_h(\mathbf{q}, t) = 0 \quad
      \mathrm{where} \quad
      \mathbf{f}_h \in \mathbb{R}^M

2. Kinematic differential equations

   .. math::
      \mathbf{M}_{k\dot{q}}(\mathbf{q}, t) \dot{\mathbf{q}} + \mathbf{M}_{ku}(\mathbf{q}, t) \mathbf{u} + \mathbf{f}_{k}(\mathbf{q}, t) = 0 \quad
      \mathrm{where} \quad
      \mathbf{M}_{k\dot{q}} \in \mathbb{R}^{n \times n}
      \mathrm{,\ }
      \mathbf{M}_{ku} \in \mathbb{R}^{n \times n}
      \mathrm{,\ }
      \mathbf{f}_{k} \in \mathbb{R}^n

3. Nonholonomic (velocity) constraints

   .. math::
      \mathbf{M}_{n}(\mathbf{q}, t) \mathbf{u} + \mathbf{f}_{n}(\mathbf{q}, t) = 0 \quad
      \mathrm{where} \quad
      \mathbf{M}_{n} \in \mathbb{R}^{m \times n}
      \mathrm{,\ }
      \mathbf{f}_{n} \in \mathbb{R}^m

   or with the dependent and independent speeds explicit:

   .. math::
      \mathbf{M}_{nr}(\mathbf{q}, t) \mathbf{u}_r + \mathbf{M}_{ns}(\mathbf{q}, t) \mathbf{u}_s + \mathbf{f}_{n}(\mathbf{q}, t) = 0 \quad
      \mathrm{where} \quad
      \mathbf{M}_{nr} \in \mathbb{R}^{m \times m}
      \mathrm{,\ }
      \mathbf{M}_{ns} \in \mathbb{R}^{m \times p}
      \mathrm{,\ }
      \mathbf{f}_{n} \in \mathbb{R}^m

4. Dynamic differential equations

   .. math::
      \mathbf{M}_{d}(\mathbf{q}, t) \dot{\mathbf{u}}_s + \mathbf{f}_{d}(\mathbf{q}, \dot{\mathbf{q}}, \mathbf{u}, t) = 0 \quad
      \mathrm{where} \quad
      \mathbf{M}_{d} \in \mathbb{R}^{p \times p}
      \mathrm{,\ }
      \mathbf{f}_{d} \in \mathbb{R}^p

5. Differentiated non-holonomic equations

   .. math::
      \mathbf{M}_{\dot{n}}(\mathbf{q}, \dot{\mathbf{q}}, t) \dot{\mathbf{u}} + \mathbf{f}_{\dot{n}}(\mathbf{q}, \dot{\mathbf{q}}, \mathbf{u}, t) = 0 \quad
      \mathrm{where} \quad
      \mathbf{M}_{\dot{n}} \in \mathbb{R}^{m \times n}
      \mathrm{,\ }
      \mathbf{f}_{\dot{n}} \in \mathbb{R}^m

   or with the dependent and independent speeds explicit:

   .. math::
      \mathbf{M}_{\dot{n}r}(\mathbf{q}, \dot{\mathbf{q}}, t) \dot{\mathbf{u}}_r + \mathbf{M}_{\dot{n}s}(\mathbf{q}, \dot{\mathbf{q}}, t) \dot{\mathbf{u}}_s + \mathbf{f}_{\dot{n}}(\mathbf{q}, \dot{\mathbf{q}}, \mathbf{u}, t) = 0 \quad
      \mathrm{where} \quad
      \mathbf{M}_{\dot{n}r} \in \mathbb{R}^{m \times m}
      \mathrm{,\ }
      \mathbf{M}_{\dot{n}s} \in \mathbb{R}^{m \times p}
      \mathrm{,\ }
      \mathbf{f}_{\dot{n}} \in \mathbb{R}^m

Equation sets 1 through 3 are provided by the user, whereas the sets 4 and
5 are computed by the :class:`~sympy.physics.mechanics.kane.KanesMethod` class.

Holonomic Constraint Equations
------------------------------

The first set of equations describes the configuration constraints which are,
in general, non-linear in the generalized coordinates. For the purposes of
forming the non-linear equations of motion, these are additional algebraic equations that
must be satisfied along with the other equations. Note that
:class:`~sympy.physics.mechanics.kane.KanesMethod` only makes use of the holonomic
constraints during linearization of the equations. If provided, it is assumed that the
equations are not easily solvable for the dependent coordinate(s) [1]_.

The expressions that make up the non-zero side of the equations should be
passed to :class:`~sympy.physics.mechanics.kane.KanesMethod` on initialization,
e.g.::

   >>> KanesMethod(..., configuration_constraints=(expr_0, expr_1, ...), ...)

where the constraint expressions are equivalent to zero.

.. [1] The equations can be either linear or non-linear in the coordinates. It
   is recommended that the user analytically solve for the independent
   generalized coordinates if possible before passing them to
   :py:class:`KanesMethod`. If you are able to easily solve a holonomic
   constraint, you should consider redefining your problem in terms of a
   smaller set of coordinates. Alternatively, the time-differentiated holonomic
   constraints can be supplied.

Kinematic Differential Equations
--------------------------------

The second set of equations are the kinematic differential equations and they
describe the relationship between the generalized speeds and the derivatives of
the generalized coordinates. These are defined by the analyst and can reduce
the length of the dynamic equations of motion if chosen carefully
[Mitiguy1996]_. The simplest and always valid choice is :math:`\mathbf{u} =
\dot{\mathbf{q}}`. These equations are needed to transform the second order
equations of motion into first order form.

These are passed into :class:`~sympy.physics.mechanics.kane.KanesMethod` class
as such::

   >>> KanesMethod(..., kd_eqs=(expr_0, expr_1), ...)

where each expression the left hand side of the above equations.

The :meth:`~sympy.physics.mechanics.kane.KanesMethod.kindiff`` method returns a
dictionary with expressions for derivatives of the generalized coordinates,
i.e.:

.. math::

   \dot{\mathbf{q}} = -\mathbf{M}_{k\dot{q}}(\mathbf{q}, t)^{-1}\left[\mathbf{M}_{ku}(\mathbf{q}, t) \mathbf{u} + \mathbf{f}_{k}(\mathbf{q}, t)\right]

Nonholonomic Constraint Equations
---------------------------------

The third set of equations describe the nonholonomic constraints, otherwise
known as velocity constraints. These must be linear in the generalized speeds.
There are fewer equations than generalized speeds, and thus describe the
relationship between the dependent and independent generalized speeds. There
are :math:`m` dependent speeds and :math:`p=n-m` independent speeds.

To solve these for the dependent speeds :math:`\mathbf{u}` can be broken into
the dependent and independent speeds:

.. math::
   \mathbf{u} = [\mathbf{u}_r, \mathbf{u}_s]^T \\
   \mathbf{u}_r = -\mathbf{M}_{nr}(\mathbf{q}, t)^{-1} \mathbf{f}_{ns}(\mathbf{q}, \mathbf{u}_s, t)

These nonholonomic constraint expressions should be passed directly to the
:class:`~sympy.physics.mechanics.kane.KanesMethod` class as such::

   >>> KanesMethod(..., velocity_constraints=(expr_0, expr_1), ...)

where each expression is one entry of the left hand side of the second set of
equations above.

Dynamic Differential Equations
------------------------------

The fourth set of equations are the dynamical differential equations. These equations are
linear in the derivatives of the generalized speeds and are initially provided as Kane's
:math:`\mathbf{F}_r + \mathbf{F}_r^* = 0`. These equations are the primary
result from executing the
:meth:`~sympy.physics.mechanics.kane.KanesMethod.kanes_equation` method::

   >>> kane = KanesMethod(...)
   >>> fr, frstar = kane.kanes_equations(bodies, loads)

Kane's equations are linear in :math:`\dot{\mathbf{u}}_s`. The coefficient mass matrix, :math:`\mathbf{M}_{d}`, can be accessed with::

   >>> kane.mass_matrix

and the forcing function, :math:`-\mathbf{f}_{d}`, can be accessed with::

   >>> kane.forcing

Note the negative sign.

Derivative of the Nonholonomic Constraint Equations
---------------------------------------------------

The fifth set of equations are the derivatives of the nonholonomic constraints. These can
be used to augment the independent dynamic equations if it is desired to solve for the dependent generalized speeds.

These can be optionally passed into
:class:`~sympy.physics.mechanics.kane.KanesMethod` as::

   >>> KanesMethod(..., acceleration_constraints=(expr_0, expr_1), ...)

where each expression is equal to zero, but otherwise they are automatically
computed from the provided velocity constraints.

Accessing the Variables and the Equations
-----------------------------------------

For a nonholonomic system with :math:`n` speeds and :math:`m` motion
constraints, we will get :math:`p = n - m` dynamic equations. The
:class:`~sympy.physics.mechanics.KanesMethod` class organizes the equations in
the following fashion:

.. math::
   \mathbf{M}(\mathbf{q}, t) &=
   \begin{bmatrix}
     \mathbf{M}_{d}(\mathbf{q}, t) & \mathbf{0}_{m \times p} \\
     \mathbf{0}_{p \times m} & \mathbf{M}_{\dot{n}r}(\mathbf{q}, t) \end{bmatrix}\\
=======
In :mod:`sympy.physics.mechanics` we are assuming there are 5 basic sets of equations needed
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

In :mod:`sympy.physics.mechanics` holonomic constraints are only used for the linearization
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
>>>>>>> master

.. math::

  \mathbf{f}(\mathbf{q}, \dot{\mathbf{q}}, \mathbf{u}, t) &=
   \begin{bmatrix}
  - \mathbf{f}_{d}(\mathbf{q}, \dot{\mathbf{q}}, \mathbf{u}, t) \\
  - \mathbf{f}_{\dot{n}s}(\mathbf{q}, \dot{\mathbf{q}}, \mathbf{u}, t)
  \end{bmatrix}\\

such that

.. math::

   \mathbf{M}(\mathbf{q}, t) \dot{\mathbf{u}} = \mathbf{f}(\mathbf{q}, \dot{\mathbf{q}}, \mathbf{u}, t)

Each component is accessed as such::

   >>> kane = KanesMethod(...)
   >>> kane.kanes_equations(bodies, loads)
   >>> kane.mass_matrix
   >>> kane.u.diff()
   >>> kane.forcing

where the total equation is::

   >>> Equality(kane.mass_matrix * kane.u.diff(), kane.forcing)

Additionally, :class:`~sympy.physics.mechanics.KanesMethod` provides the
combined dynamic and kinematic equations:

.. math::
  \tilde{\mathbf{M}}(\mathbf{q}, t) &=
   \begin{bmatrix}
     \mathbf{M}(\mathbf{q}, t) & \mathbf{0}_{n \times n} \\
     \mathbf{0}_{n \times n} & \mathbf{M}_{k\dot{q}}(\mathbf{q}, t) \end{bmatrix}\\

.. math::

  \tilde{\mathbf{f}}(\mathbf{q}, \dot{\mathbf{q}}, \mathbf{u}, t) &=
   \begin{bmatrix}
     \mathbf{f}(\mathbf{q}, \dot{\mathbf{q}}, \mathbf{u}, t) \\
     - \mathbf{M}_{ku}(\mathbf{q}, t) \mathbf{u} - \mathbf{f}_{k}(\mathbf{q}, t)
  \end{bmatrix}\\

Each component is accessed as such::

   >>> kane = KanesMethod(...)
   >>> kane.kanes_equations(bodies, loads)
   >>> kane.mass_matrix_full
   >>> kane.u.diff()
   >>> kane.q.diff()
   >>> kane.forcing_full

where the total equation is::

   >>> Equality(kane.mass_matrix_full * kane.u.col_join(kane.q).diff(), kane.forcing_full)

Simple Example
==============

The formulation of the equations of motion in
:mod:`~sympy.physics.mechanics` starts with creation of a ``KanesMethod``
object. Upon initialization of the ``KanesMethod`` object, an inertial
reference frame needs to be supplied. along with some basic system information,
such as coordinates and speeds::

The formulation of the equations of motion in :mod:`sympy.physics.mechanics` starts with
creation of a ``KanesMethod`` object. Upon initialization of the
``KanesMethod`` object, an inertial reference frame needs to be supplied. along
with some basic system information, such as coordinates and speeds ::

  >>> from sympy.physics.mechanics import *
  >>> N = ReferenceFrame('N')
  >>> q1, q2, u1, u2 = dynamicsymbols('q1 q2 u1 u2')
  >>> q1d, q2d, u1d, u2d = dynamicsymbols('q1 q2 u1 u2', 1)
  >>> KM = KanesMethod(N, [q1, q2], [u1, u2])

It is also important to supply the order of coordinates and speeds properly if
there are dependent coordinates and speeds. They must be supplied after
independent coordinates and speeds or as a keyword argument; this is shown
later.::

  >>> q1, q2, q3, q4 = dynamicsymbols('q1 q2 q3 q4')
  >>> u1, u2, u3, u4 = dynamicsymbols('u1 u2 u3 u4')
  >>> # Here we will assume q2 is dependent, and u2 and u3 are dependent
  >>> # We need the constraint equations to enter them though
  >>> KM = KanesMethod(N, [q1, q3, q4], [u1, u4])

Additionally, if there are auxiliary speeds, they need to be identified here.
See the examples for more information on this. In this example ``u4`` is the
auxiliary speed.::

  >>> KM = KanesMethod(N, [q1, q3, q4], [u1, u2, u3], u_auxiliary=[u4])

Kinematic differential equations must also be supplied; there are to be
provided as a list of expressions which are each equal to zero. A trivial
example follows::

  >>> kd = [q1d - u1, q2d - u2]

Turning on ``mechanics_printing()`` makes the expressions significantly shorter
and is recommended. Alternatively, the ``mprint`` and ``mpprint`` commands can
be used.

If there are non-holonomic constraints, dependent speeds need to be specified
(and so do dependent coordinates, but they only come into play when linearizing
the system). The constraints need to be supplied in a list of expressions which
are equal to zero, trivial motion and configuration constraints are shown
below::

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

A dictionary returning the solved :math:`\dot{q}`'s can also be solved for::

  >>> mechanics_printing(pretty_print=False)
  >>> KM.kindiffdict()
  {q1': u1, q2': u2, q3': u3, q4': u4}

The final step in forming the equations of motion is supplying a list of bodies
and particles, and a list of 2-tuples of the form ``(Point, Vector)`` or
``(ReferenceFrame, Vector)`` to represent applied forces and torques.::

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
:math:`\mathbf{M}_{\dot{n}}(\mathbf{q}, t)` matrix, and the forcing vector by
the :math:`\mathbf{f}_{\dot{n}}(\mathbf{q}, \dot{\mathbf{q}}, \mathbf{u}, t)`
vector.

There are also the "full" mass matrix and "full" forcing vector terms, these
include the kinematic differential equations; the mass matrix is of size (m +
p) x (m + p), or square and the size of all coordinates and speeds.::

  >>> KM.mass_matrix_full
  Matrix([
  [1, 0],
  [0, 5]])
  >>> KM.forcing_full
  Matrix([
  [u],
  [7]])

Exploration of the provided examples is encouraged in order to gain more
understanding of the :py:class:`KanesMethod` object.
