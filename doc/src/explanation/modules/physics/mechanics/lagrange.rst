======================================
Lagrange's Method in Physics/Mechanics
======================================

:mod:`sympy.physics.mechanics` provides functionality for deriving equations of
motion using `Lagrange's method
<https://en.wikipedia.org/wiki/Lagrangian_mechanics>`_. This page describes
Lagrange's method as used in this module, but not how Lagrange's equations are
derived.

Structure of Equations
======================

In :mod:`sympy.physics.mechanics` there are three sets of equations needed to
describe a system: the holonomic constraint equations, the nonholonomic
constraint equations and the dynamical differential equations.

The :math:`M` holonomic constraints are:

.. math::

   \mathbf{f}_h(\mathbf{q}, t) = \mathbf{0}

The :math:`m` nonholonomic constraints are:

.. math::

   \mathbf{m}_n(\mathbf{q}, t) \dot{\mathbf{q}} + \mathbf{f}_n(\mathbf{q}, t) = \mathbf{0}

These can be combined into :math:`M + m` velocity constraints (time
differentiated holonomic constraints with the nonholonomic constraints):

.. math::

   \dot{\mathbf{f}}_h(\mathbf{q}, t) =
   \mathbf{m}_{hd} \dot{\mathbf{q}} +
   \mathbf{f}_{hd} = \mathbf{0}

.. math::

   \begin{bmatrix}
     \mathbf{m}_{hd} \\
     \mathbf{m}_n
   \end{bmatrix}
   \dot{\mathbf{q}} +
   \begin{bmatrix}
     \mathbf{f}_{hd} \\
     \mathbf{f}_n
   \end{bmatrix} =
   \mathbf{m}_v(\mathbf{q}, t) \dot{\mathbf{q}} + \mathbf{f}_v(\mathbf{q}, t) = 0

Time differentiating the velocity constraints then gives the acceleration level
constraints:

.. math::

   \mathbf{m}_v(\mathbf{q}, t)\ddot{\mathbf{q}} + \mathbf{f}_a(\dot{\mathbf{q}}, \mathbf{q}, t) = 0

:math:`\mathbf{m}_v` is called "the Jacobian of the constraints" and can be
used to augment Lagrange's :math:`n` dynamical differential equations with
constraint forces:

.. math::

  \mathbf{m}_d(\mathbf{q}, t) \ddot{\mathbf{q}} + \mathbf{\Lambda}(\mathbf{q}, t) \mathbf{\lambda} + \mathbf{f}_d(\dot{\mathbf{q}}, \mathbf{q}, t) = 0

where :math:`\mathbf{\Lambda}(\mathbf{q}, t) = \mathbf{m}_v^T` and :math:`\mathbf{\lambda}` are
the :math:`M + m` unknown Lagrange multipliers.

The complete equations of motion take the final form:

.. math::

   \mathbf{M}(\mathbf{q}, t) x = \mathbf{F}(\dot{\mathbf{q}}, \mathbf{q}, t)

where:

.. math::

   \mathbf{x} =
   \begin{bmatrix}
     \ddot{\mathbf{q}} \\
     \mathbf{\lambda}
   \end{bmatrix}

For a constrained system with `n` generalized coordinates and `m + M`
constraints, there are :math:`n` equations of the form:

.. math::

   \mathbf{M}(\mathbf{q}, t)\mathbf{x} =
   \begin{bmatrix}
     \mathbf{m}_d(\mathbf{q}, t) & \mathbf{\Lambda}(\mathbf{q}, t)
   \end{bmatrix}\mathbf{x} =
   -\mathbf{f}_d(\dot{\mathbf{q}}, \mathbf{q}, t) =
   \mathbf{F}(\dot{\mathbf{q}}, \mathbf{q}, t)

Lagrange's Method in Physics/Mechanics
======================================

The formulation of the equations of motion in :mod:`sympy.physics.mechanics`
using Lagrange's Method starts with the creation of generalized coordinates and
a Lagrangian. The Lagrangian can either be created with the
:py:func:`~.Lagrangian` function or can be a user supplied expression. In this
case we will supply the Lagrangian.::

  >>> from sympy.physics.mechanics import (dynamicsymbols, LagrangesMethod,
  ...     mechanics_printing)
  >>> q1, q2 = dynamicsymbols('q1 q2')
  >>> q1d, q2d = dynamicsymbols('q1 q2', 1)
  >>> lagrangian = q1d**2 + q2d**2 - 3*q1 - 4*q2

To formulate the equations of motion we create a ``LagrangesMethod``
object. The Lagrangian and generalized coordinates need to be supplied upon
initialization. ::

  >>> lm = LagrangesMethod(lagrangian, (q1, q2))

With that, the equations of motion can be formed::

  >>> mechanics_printing(pretty_print=False)
  >>> lm.form_lagranges_equations()
  Matrix([
  [2*q1'' + 3],
  [2*q2'' + 4]])

The mass matrix :math:`\mathbf{M}` and the forcing vector :math:`\mathbf{F}`
are now accessible::

  >>> lm.mass_matrix
  Matrix([
  [2, 0],
  [0, 2]])
  >>> lm.forcing
  Matrix([
  [-3],
  [-4]])

If there are any holonomic or non-holonomic constraints, they must be supplied
as keyword arguments (``hol_coneqs`` and ``nonhol_coneqs`` respectively) in a
list of expressions which are equal to zero. Modifying the example above, the
equations of motion can then be generated::

  >>> lm = LagrangesMethod(lagrangian, (q1, q2), hol_coneqs=[q1 - q2])

When the equations of motion are generated in this case, the Lagrange
multipliers are introduced; they are represented by ``lam1`` here. In general,
there will be as many multipliers as there are constraint equations.

::

  >>> lm.form_lagranges_equations()
  Matrix([
  [ lam1 + 2*q1'' + 3],
  [-lam1 + 2*q2'' + 4]])

The Lagrange multipliers :math:`\mathbf{\lambda}` are accessed with::

  >>> lm.lam_vec
  Matrix([[lam1]])

The :math:`-\mathbf{\Lambda}` matrix is accessed with::

  >>> lm.lam_coeffs
  Matrix([[-1, 1]])

The augmented mass matrix :math:`\left[ \mathbf{m}_d \quad
-\mathbf{\Lambda}\right]` is then::

  >>> lm.mass_matrix
  Matrix([
  [2, 0, -1],
  [0, 2,  1]])
  >>> lm.forcing
  Matrix([
  [-3],
  [-4]])

In the case of systems with constraints, the 'full' mass matrix includes the
kinematical differential equations. The 'full' mass matrix is of size (2n + M +
m) x (2n + M + m), i.e. it's a square matrix. ::

  >>> lm.mass_matrix_full
  Matrix([
  [1, 0, 0,  0,  0],
  [0, 1, 0,  0,  0],
  [0, 0, 2,  0, -1],
  [0, 0, 0,  2,  1],
  [0, 0, 1, -1,  0]])
  >>> lm.forcing_full
  Matrix([
  [q1'],
  [q2'],
  [ -3],
  [ -4],
  [  0]])

If there are any non-conservative forces or moments acting on the system, they
must also be supplied as keyword arguments in a list of 2-tuples of the form
(:py:obj:`~.vector.point.Point`, :py:obj:`~.physics.vector.vector.Vector`) or
(:py:obj:`~.ReferenceFrame`, :py:obj:`~.physics.vector.vector.Vector`) where
the ``Vector`` represents the non-conservative forces and torques. Along with
this 2-tuple, the inertial frame must also be specified as a keyword argument.
This is shown below by modifying the example above::

  >>> from sympy.physics.mechanics import ReferenceFrame, Point
  >>> N = ReferenceFrame('N')
  >>> P = Point('P')
  >>> P.set_vel(N, q1d*N.x)
  >>> loads = [(P, 7*N.x)]
  >>> lm = LagrangesMethod(lagrangian, (q1, q2), forcelist=loads, frame=N)
  >>> lm.form_lagranges_equations()
  Matrix([
  [2*q1'' - 4],
  [2*q2'' + 4]])

Exploration of the provided :ref:`examples <mechanics_tutorial>` is encouraged
in order to gain more understanding of the :py:class:`~.LagrangesMethod` class.
