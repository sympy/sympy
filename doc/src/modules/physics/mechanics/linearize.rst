==================================
Linearization in Physics/Mechanics
==================================

:mod:`mechanics` includes methods for linearizing the generated equations of
motion (EOM) about an operating point (also known as the trim condition).
Note that this operating point doesn't have to be an equilibrium position, it
just needs to satisfy the equations of motion.

Linearization is accomplished by taking the first order Taylor expansion of
the EOM about the operating point. When there are no dependent coordinates
or speeds this is simply the jacobian of the right hand side about `q` and `u`.
However, in the presence of constraints more care needs to be taken. The
linearization methods provided here handle these constraints correctly.

Background
==========

In :mod:`mechanics` we assume all systems can be represented in the following
general form:

.. math::
  f_{c}(q, t) &= 0_{l \times 1}\\
  f_{v}(q, u, t) &= 0_{m \times 1}\\
  f_{a}(q, \dot{q}, u, \dot{u}, t) &= 0_{m \times 1}\\
  f_{0}(q, \dot{q}, t) + f_{1}(q, u, t) &= 0_{n \times 1}\\
  f_{2}(q, u, \dot{u}, t) + f_{3}(q, \dot{q}, u, r, t) +
  f_{4}(q, \lambda, t) &= 0_{(o-m+k) \times 1}

where

.. math::
  q, \dot{q} & \in \mathbb{R}^n\\
  u, \dot{u} & \in \mathbb{R}^o\\
  r & \in \mathbb{R}^s\\
  \lambda & \in \mathbb{R}^k

In this form,

- :math:`f_{c}` represents the configuration constraint equations
- :math:`f_{v}` represents the velocity constraint equations
- :math:`f_{a}` represents the acceleration constraint equations
- :math:`f_{0}` and :math:`f_{1}` form the kinematic differential equations
- :math:`f_{2}`, :math:`f_{3}`, and :math:`f_{4}` form the dynamic differential equations
- :math:`q` and :math:`\dot{q}` are the generalized coordinates and their derivatives
- :math:`u` and :math:`\dot{u}` are the generalized speeds and their derivatives
- :math:`r` is the system inputs
- :math:`\lambda` is the Lagrange multipliers

This generalized form is held inside the ``Linearizer`` class, which
performs the actual linearization. Both ``KanesMethod`` and
``LagrangesMethod`` objects have methods for forming the linearizer using
the ``to_linearizer`` class method.

.. topic:: A Note on Dependent Coordinates and Speeds

  If the system being linearized contains constraint equations, this results in
  not all generalized coordinates being independent (i.e. `q_1` may depend on
  `q_2`). With `l` configuration constraints, and `m` velocity constraints,
  there are `l` dependent coordinates and `m` dependent speeds.

  In general, you may pick any of the coordinates and speeds to be dependent,
  but in practice some choices may result in undesirable singularites. Methods
  for deciding which coordinates/speeds to make dependent is behind the scope of
  this guide. For more information, please see [Blajer1994]_.

Once the system is coerced into the generalized form, the linearized EOM can be
solved for. The methods provided in :mod:`mechanics` allow for two different
forms of the linearized EOM:

`M`, `A`, and `B`
  In this form, the forcing matrix is linearized into two separate matrices `A`
  and `B`. This is the default form of the linearized EOM. The resulting
  equations are:

  .. math::
    M \begin{bmatrix} \delta \dot{q} \\ \delta \dot{u} \\ \delta \lambda \end{bmatrix} =
    A \begin{bmatrix} \delta q_i \\ \delta u_i \end{bmatrix} + B \begin{bmatrix} \delta r \end{bmatrix}

  where

  .. math::
    M &\in \mathbb{R}^{(n+o+k) \times (n+o+k)}\\
    A &\in \mathbb{R}^{(n+o+k) \times (n-l+o-m)}\\
    B &\in \mathbb{R}^{(n+o+k) \times s}

  Note that `q_i` and `u_i` are just the independent coordinates and speeds,
  while `q` and `u` contains both the independent and dependent coordinates
  and speeds.

`A` and `B`
  In this form, the linearized EOM are brought into explicit first order form,
  in terms of just the independent coordinates and speeds. This form is often
  used in stability analysis or control theory. The resulting equations
  are:

  .. math::
    \begin{bmatrix} \delta \dot{q_i} \\ \delta \dot{u_i} \end{bmatrix} =
    A \begin{bmatrix} \delta q_i \\ \delta u_i \end{bmatrix} + B \begin{bmatrix} \delta r \end{bmatrix}

  where

  .. math::
    A &\in \mathbb{R}^{(n-l+o-m) \times (n-l+o-m)}\\
    B &\in \mathbb{R}^{(n-l+o-m) \times s}

  To use this form set ``A_and_B=True`` in the ``linearize`` class method.

Linearizing Kane's Equations
============================

After initializing the ``KanesMethod`` object and forming `F_r` and `F_r^*`
using the ``kanes_equations`` class method, linearization can be accomplished
in a couple ways. The different methods will be demonstrated with a simple
pendulum system: ::

  >>> from sympy import symbols, Matrix
  >>> from sympy.physics.mechanics import *
  >>> q1 = dynamicsymbols('q1')                     # Angle of pendulum
  >>> u1 = dynamicsymbols('u1')                     # Angular velocity
  >>> q1d = dynamicsymbols('q1', 1)
  >>> L, m, t, g = symbols('L, m, t, g')

  >>> # Compose world frame
  >>> N = ReferenceFrame('N')
  >>> pN = Point('N*')
  >>> pN.set_vel(N, 0)

  >>> # A.x is along the pendulum
  >>> A = N.orientnew('A', 'axis', [q1, N.z])
  >>> A.set_ang_vel(N, u1*N.z)

  >>> # Locate point P relative to the origin N*
  >>> P = pN.locatenew('P', L*A.x)
  >>> vel_P = P.v2pt_theory(pN, N, A)
  >>> pP = Particle('pP', P, m)

  >>> # Create Kinematic Differential Equations
  >>> kde = Matrix([q1d - u1])

  >>> # Input the force resultant at P
  >>> R = m*g*N.x

  >>> # Solve for eom with kanes method
  >>> KM = KanesMethod(N, q_ind=[q1], u_ind=[u1], kd_eqs=kde)
  >>> fr, frstar = KM.kanes_equations([pP], [(P, R)])

1. Using the ``Linearizer`` class directly:
-------------------------------------------

A linearizer object can be created using the ``to_linearizer`` class method.
This coerces the representation found in the ``KanesMethod`` object into the
generalized form described above. As the independent and dependent
coordinates and speeds are specified upon creation of the KanesMethod object,
there is no need to specifiy them here. ::

  >>> linearizer = KM.to_linearizer()

The linearized EOM can then be formed with the ``linearize`` method of the
``Linearizer`` object: ::

  >>> M, A, B = linearizer.linearize()
  >>> M
  Matrix([
  [1,       0],
  [0, -L**2*m]])
  >>> A
  Matrix([
  [                 0, 1],
  [L*g*m*cos(q1(t)), 0]])
  >>> B
  Matrix(0, 0, [])

Alternatively, the `A` and `B` form can be generated instead by specifying
``A_and_B=True``: ::

  >>> A, B = linearizer.linearize(A_and_B=True)
  >>> A
  Matrix([
  [                0, 1],
  [-g*cos(q1(t))/L, 0]])
  >>> B
  Matrix(0, 0, [])

An operating point can also be specified as a dictionary or an iterable of
dictionaries. This will evaluate the linearized form at the specified
point before returning the matrices: ::

  >>> op_point = {q1: 0, u1: 0}
  >>> A_op, B_op = linearizer.linearize(A_and_B=True, op_point=op_point)
  >>> A_op
  Matrix([
  [     0, 1],
  [-g/L, 0]])

Note that the same effect can be had by applying ``msubs`` to the matrices
generated without the ``op_point`` kwarg: ::

  >>> assert msubs(A, op_point) == A_op

Sometimes the returned matrices may not be in the most simplified form.
Simplification can be performed after the fact, or the ``Linearizer`` object
can be made to perform simplification internally by setting the ``simplify``
kwarg to ``True``.

2. Using the ``linearize`` class method:
----------------------------------------

The ``linearize`` method of the ``KanesMethod`` class is provided as a nice
wrapper that calls ``to_linearizer`` internally, performs the linearization,
and returns the result. Note that all the kwargs available in the
``linearize`` method described above are also available here: ::

  >>> A, B, inp_vec = KM.linearize(A_and_B=True, op_point=op_point, new_method=True)
  >>> A
  Matrix([
  [     0, 1],
  [-g/L, 0]])

The additional output ``inp_vec`` is a vector containing all found
``dynamicsymbols`` not included in the generalized coordinate or speed
vectors. These are assumed to be inputs to the system, forming the `r` vector
described in the background above. In this example there are no inputs, so
the vector is empty: ::

  >>> inp_vec
  Matrix(0, 0, [])

.. topic:: What's with the ``new_method`` kwarg?

  Previous releases of :mod:`SymPy` contained a linearization method for
  `KanesMethod`` objects. This method is deprecated, and will be removed
  from future releases. Until then, you must set ``new_method=True`` in all
  calls to ``KanesMethod.linearize``. After the old method is removed, this
  kwarg will no longer be needed.

Linearizing Lagrange's Equations
================================

Linearization of Lagrange's equations proceeds much the same as that of
Kane's equations. As before, the process will be demonstrated with a simple
pendulum system: ::

  >>> # Redefine A and P in terms of q1d, not u1
  >>> A = N.orientnew('A', 'axis', [q1, N.z])
  >>> A.set_ang_vel(N, q1d*N.z)
  >>> P = pN.locatenew('P', L*A.x)
  >>> vel_P = P.v2pt_theory(pN, N, A)
  >>> pP = Particle('pP', P, m)

  >>> # Solve for eom with Lagranges method
  >>> Lag = Lagrangian(N, pP)
  >>> LM = LagrangesMethod(Lag, [q1], forcelist=[(P, R)], frame=N)
  >>> lag_eqs = LM.form_lagranges_equations()

1. Using the ``Linearizer`` class directly:
-------------------------------------------

A ``Linearizer`` object can be formed from a ``LagrangesMethod`` object using
the ``to_linearizer`` class method. The only difference between this process
and that of the ``KanesMethod`` class is that the ``LagrangesMethod`` object
doesn't already have its independent and dependent coordinates and speeds
specified internally. These must be specified in the call to
``to_linearizer``. In this example there are no dependent coordinates and
speeds, but if there were they would be included in the ``q_dep`` and
``qd_dep`` kwargs: ::

  >>> linearizer = LM.to_linearizer(q_ind=[q1], qd_ind=[q1d])

Once in this form, everything is the same as it was before with the
``KanesMethod`` example: ::

  >>> A, B = linearizer.linearize(A_and_B=True, op_point=op_point)
  >>> A
  Matrix([
  [     0, 1],
  [-g/L, 0]])

2. Using the ``linearize`` class method:
----------------------------------------

Similar to ``KanesMethod``, the ``LagrangesMethod`` class also provides a
``linearize`` method as a nice wrapper that calls ``to_linearizer``
internally, performs the linearization, and returns the result. As before, the
only difference is that the independent and dependent coordinates and speeds
must be specified in the call as well: ::

  >>> A, B, inp_vec = LM.linearize(q_ind=[q1], qd_ind=[q1d], A_and_B=True, op_point=op_point)
  >>> A
  Matrix([
  [     0, 1],
  [-g/L, 0]])

Potential Issues
================

While the ``Linearizer`` class *should* be able to linearize all systems,
there are some potential issues that could occur. These are discussed below,
along with some troubleshooting tips for solving them.

1. Symbolic linearization with ``A_and_B=True`` is slow
-------------------------------------------------------
This could be due to a number of things, but the most likely one is that
solving a large linear system symbolically is an expensive operation.
Specifying an operating point will reduce the expression size and speed
this up. If a purely symbolic solution is desired though (for application
of many operating points at a later period, for example) a way to get
around this is to evaluate with ``A_and_B=False``, and then solve
manually after applying the operating point: ::

  >>> M, A, B = linearizer.linearize()
  >>> M_op = msubs(M, op_point)
  >>> A_op = msubs(A, op_point)
  >>> perm_mat = linearizer.perm_mat
  >>> A_lin = perm_mat.T * M_op.LUsolve(A_op)
  >>> A_lin
  Matrix([
  [     0, 1],
  [-g/L, 0]])

The fewer symbols in ``A`` and ``M`` before solving, the faster this
solution will be. Thus, for large expressions, it may be to your benefit
to delay conversion to the `A` and `B` form until most symbols are subbed
in for their numeric values.

2. The linearized form has ``nan``, ``zoo``, or ``oo`` as matrix elements
-------------------------------------------------------------------------
There are two potential causes for this. The first (and the one you
should check first) is that some choices of dependent coordinates
will result in singularities at certain operating points. Coordinate
partitioning in a systemic manner to avoid this is beyond the scope
of this guide; see [Blajer1994]_ for more information.

The other potential cause for this is that the matrices may not have
been in the most reduced form before the operating point was substituted
in. A simple example of this behavior is: ::

  >>> from sympy import sin, tan
  >>> expr = sin(q1)/tan(q1)
  >>> op_point = {q1: 0}
  >>> expr.subs(op_point)
  nan

Note that if this expression was simplified before substitution, the
correct value results: ::

  >>> expr.simplify().subs(op_point)
  1

A good way of avoiding this hasn't been found yet. For expressions of
reasonable size, using ``msubs`` with ``smart=True`` will apply an
algorithm that tries to avoid these conditions. For large expressions
though this is extremely time consuming. ::

  >>> msubs(expr, op_point, smart=True)
  1

Further Examples
================

The pendulum example used above was simple, but didn't include any dependent
coordinates or speeds. For a more thorough example, the same pendulum
was linearized with dependent coordinates using both Kane's and Lagrange's
methods:

.. toctree::

    examples/lin_pend_nonmin_example.rst
