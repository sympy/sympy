===============================
Nonminimal Coordinates Pendulum
===============================

In this example we demonstrate the use of the functionality provided in
:mod:`mechanics` for deriving the equations of motion (EOM) for a pendulum
with a nonminimal set of coordinates. As the pendulum is a one degree of
freedom system, it can be described using one coordinate and one speed (the
pendulum angle, and the angular velocity respectively). Choosing instead to
describe the system using the `x` and `y` coordinates of the mass results in
a need for constraints. The system is shown below:

.. image:: pendulum_nonmin.*
   :align: center

The system will be modeled using both Kane's and Lagrange's methods, and the
resulting EOM linearized. While this is a simple problem, it should illustrate
the use of the linearization methods in the presence of constraints.

Kane's Method
=============

First we need to create the ``dynamicsymbols`` needed to describe the system as
shown in the above diagram. In this case, the generalized coordinates `q_1` and
`q_2` represent the mass `x` and `y` coordinates in the inertial `N` frame.
Likewise, the generalized speeds `u_1` and `u_2` represent the velocities in
these directions. We also create some ``symbols`` to represent the length and
mass of the pendulum, as well as gravity and time. ::

  >>> from sympy.physics.mechanics import *
  >>> from sympy import symbols, atan, Matrix, solve
  >>> # Create generalized coordinates and speeds for this non-minimal realization
  >>> # q1, q2 = N.x and N.y coordinates of pendulum
  >>> # u1, u2 = N.x and N.y velocities of pendulum
  >>> q1, q2 = dynamicsymbols('q1:3')
  >>> q1d, q2d = dynamicsymbols('q1:3', level=1)
  >>> u1, u2 = dynamicsymbols('u1:3')
  >>> u1d, u2d = dynamicsymbols('u1:3', level=1)
  >>> L, m, g, t = symbols('L, m, g, t')

Next, we create a world coordinate frame `N`, and its origin point `N^*`. The
velocity of the origin is set to 0. A second coordinate frame `A` is oriented
such that its x-axis is along the pendulum (as shown in the diagram above). ::

  >>> # Compose world frame
  >>> N = ReferenceFrame('N')
  >>> pN = Point('N*')
  >>> pN.set_vel(N, 0)

  >>> # A.x is along the pendulum
  >>> theta1 = atan(q2/q1)
  >>> A = N.orientnew('A', 'axis', [theta1, N.z])

Locating the pendulum mass is then as easy as specifying its location with in
terms of its x and y coordinates in the world frame. A ``Particle`` object is
then created to represent the mass at this location. ::

  >>> # Locate the pendulum mass
  >>> P = pN.locatenew('P1', q1*N.x + q2*N.y)
  >>> pP = Particle('pP', P, m)

The kinematic differential equations (KDEs) relate the derivatives of the
generalized coordinates to the generalized speeds. In this case the speeds are
the derivatives, so these are simple. A dictionary is also created to map
`\dot{q}` to `u`: ::

  >>> # Calculate the kinematic differential equations
  >>> kde = Matrix([q1d - u1,
  ...               q2d - u2])
  >>> dq_dict = solve(kde, [q1d, q2d])

The velocity of the mass is then the time derivative of the position from the
origin `N^*`: ::

  >>> # Set velocity of point P
  >>> P.set_vel(N, P.pos_from(pN).dt(N).subs(dq_dict))

As this system has more coordinates than degrees of freedom, constraints are
needed. The configuration constraints relate the coordinates to each other. In
this case the constraint is that the distance from the origin to the mass is
always the length `L` (the pendulum doesn't get longer). Likewise, the velocity
constraint is that the mass velocity in the ``A.x`` direction is always 0 (no
radial velocity). ::

  >>> f_c = Matrix([P.pos_from(pN).magnitude() - L])
  >>> f_v = Matrix([P.vel(N).express(A).dot(A.x)])
  >>> f_v.simplify()

The force on the system is just gravity, at point ``P``. ::

  >>> # Input the force resultant at P
  >>> R = m*g*N.x

With the problem setup, the equations of motion can be generated using the
``KanesMethod`` class. As there are constraints, dependent and independent
coordinates need to be provided to the class. In this case we'll use `q_2` and
`u_2` as the independent coordinates and speeds: ::

  >>> # Derive the equations of motion using the KanesMethod class.
  >>> KM = KanesMethod(N, q_ind=[q2], u_ind=[u2], q_dependent=[q1],
  ...                  u_dependent=[u1], configuration_constraints=f_c,
  ...                  velocity_constraints=f_v, kd_eqs=kde)
  >>> (fr, frstar) = KM.kanes_equations([pP],[(P, R)])

For linearization, operating points can be specified on the call, or be
substituted in afterwards. In this case we'll provide them in the call,
supplied in a list.  The ``A_and_B=True`` kwarg indicates to solve invert the
`M` matrix and solve for just the explicit linearized `A` and `B` matrices. The
``simplify=True`` kwarg indicates to simplify inside the linearize call, and
return the presimplified matrices. The cost of doing this is small for simple
systems, but for larger systems this can be a costly operation, and should be
avoided. ::

  >>> # Set the operating point to be straight down, and non-moving
  >>> q_op = {q1: L, q2: 0}
  >>> u_op = {u1: 0, u2: 0}
  >>> ud_op = {u1d: 0, u2d: 0}
  >>> # Perform the linearization
  >>> A, B, inp_vec = KM.linearize(op_point=[q_op, u_op, ud_op], A_and_B=True,
  ...                              new_method=True, simplify=True)
  >>> A
  Matrix([
  [   0, 1],
  [-g/L, 0]])
  >>> B
  Matrix(0, 0, [])

The resulting `A` matrix has dimensions 2 x 2, while the number of total states
is ``len(q) + len(u) = 2 + 2 = 4``. This is because for constrained systems the
resulting ``A_and_B`` form has a partitioned state vector only containing
the independent coordinates and speeds. Written out mathematically, the system
linearized about this point would be written as:

.. math::
  \begin{bmatrix} \dot{q_2} \\ \dot{u_2} \end{bmatrix} =
  \begin{bmatrix} 0 & 1 \\ \frac{-g}{L} & 0 \end{bmatrix}
  \begin{bmatrix} q_2 \\ u_2 \end{bmatrix}


Lagrange's Method
=================

The derivation using Lagrange's method is very similar to the approach using
Kane's method described above. As before, we first create the
``dynamicsymbols`` needed to describe the system. In this case, the generalized
coordinates `q_1` and `q_2` represent the mass `x` and `y` coordinates in the
inertial `N` frame.  This results in the time derivatives `\dot{q_1}` and
`\dot{q_2}` representing the velocities in these directions. We also create some
``symbols`` to represent the length and mass of the pendulum, as well as
gravity and time. ::

  >>> from sympy.physics.mechanics import *
  >>> from sympy import symbols, atan, Matrix
  >>> q1, q2 = dynamicsymbols('q1:3')
  >>> q1d, q2d = dynamicsymbols('q1:3', level=1)
  >>> L, m, g, t = symbols('L, m, g, t')

Next, we create a world coordinate frame `N`, and its origin point `N^*`. The
velocity of the origin is set to 0. A second coordinate frame `A` is oriented
such that its x-axis is along the pendulum (as shown in the diagram above). ::

  >>> # Compose World Frame
  >>> N = ReferenceFrame('N')
  >>> pN = Point('N*')
  >>> pN.set_vel(N, 0)
  >>> # A.x is along the pendulum
  >>> theta1 = atan(q2/q1)
  >>> A = N.orientnew('A', 'axis', [theta1, N.z])

Locating the pendulum mass is then as easy as specifying its location with in
terms of its x and y coordinates in the world frame. A ``Particle`` object is
then created to represent the mass at this location. ::

  >>> # Create point P, the pendulum mass
  >>> P = pN.locatenew('P1', q1*N.x + q2*N.y)
  >>> P.set_vel(N, P.pos_from(pN).dt(N))
  >>> pP = Particle('pP', P, m)

As this system has more coordinates than degrees of freedom, constraints are
needed. In this case only a single holonomic constraints is needed: the
distance from the origin to the mass is always the length `L` (the pendulum
doesn't get longer). ::

  >>> # Holonomic Constraint Equations
  >>> f_c = Matrix([q1**2 + q2**2 - L**2])

The force on the system is just gravity, at point ``P``. ::

  >>> # Input the force resultant at P
  >>> R = m*g*N.x

With the problem setup, the Lagrangian can be calculated, and the equations of
motion formed. Note that the call to ``LagrangesMethod`` includes the
Lagrangian, the generalized coordinates, the constraints (specified by
``hol_coneqs`` or ``nonhol_coneqs``), the list of (body, force) pairs, and the
inertial frame. In contrast to the ``KanesMethod`` initializer, independent and
dependent coordinates are not partitioned inside the ``LagrangesMethod``
object. Such a partition is supplied later. ::

  >>> # Calculate the lagrangian, and form the equations of motion
  >>> Lag = Lagrangian(N, pP)
  >>> LM = LagrangesMethod(Lag, [q1, q2], hol_coneqs=f_c, forcelist=[(P, R)], frame=N)
  >>> lag_eqs = LM.form_lagranges_equations()

Next, we compose the operating point dictionary, set in the hanging at rest
position: ::

  >>> # Compose operating point
  >>> op_point = {q1: L, q2: 0, q1d: 0, q2d: 0, q1d.diff(t): 0, q2d.diff(t): 0}

As there are constraints in the formulation, there will be corresponding
Lagrange Multipliers. These may appear inside the linearized form as well, and
thus should also be included inside the operating point dictionary.
Fortunately, the ``LagrangesMethod`` class provides an easy way of solving
for the multipliers at a given operating point using the ``solve_multipliers``
method. ::

  >>> # Solve for multiplier operating point
  >>> lam_op = LM.solve_multipliers(op_point=op_point)

With this solution, linearization can be completed. Note that in contrast to
the ``KanesMethod`` approach, the ``LagrangesMethod.linearize`` method also
requires the partitioning of the generalized coordinates and their time
derivatives into independent and dependent vectors.  This is the same as what
was passed into the ``KanesMethod`` constructor above:

  >>> op_point.update(lam_op)
  >>> # Perform the Linearization
  >>> A, B, inp_vec = LM.linearize([q2], [q2d], [q1], [q1d],
  ...                             op_point=op_point, A_and_B=True)
  >>> A
  Matrix([
  [     0, 1],
  [-g/L, 0]])
  >>> B
  Matrix(0, 0, [])

The resulting `A` matrix has dimensions 2 x 2, while the number of total states
is ``2*len(q) = 4``. This is because for constrained systems the resulting
``A_and_B`` form has a partitioned state vector only containing the independent
coordinates and their derivatives. Written out mathematically, the system
linearized about this point would be written as:

.. math::
  \begin{bmatrix} \dot{q_2} \\ \ddot{q_2} \end{bmatrix} =
  \begin{bmatrix} 0 & 1 \\ \frac{-g}{L} & 0 \end{bmatrix}
  \begin{bmatrix} q_2 \\ \dot{q_2} \end{bmatrix}
