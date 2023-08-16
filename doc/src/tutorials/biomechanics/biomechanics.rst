======================================
Introduction to Biomechanical Modeling
======================================

:obj:`~sympy.physics._biomechanics` provides features to enhance models created
with :obj:`~sympy.physics.mechanics` with force producing elements that model
muscles and other biomechanical elements. In this tutorial, we will introduce
the features of this package by adding muscles to a simple model of a human arm
that moves a lever.

Model Description
=================

.. _fig-biomechanics-steerer:
.. figure:: biomechanics-steerer.svg

   TODO: add caption

The lever :math:`A` can rotate about :math:`\hat{n}_z` through angle
:math:`q_1`. The shoulder is located at :math:`P_2` and the upper arm :math:`C`
can extend about :math:`\hat{n}_y` through angle :math:`q_2` and rotate about
:math:`\hat{b}_z` through angle :math:`q_3`. The elbow is located at point
:math:`P_3`.  The lower arm can extend about :math:`\hat{c}_y` through angle
:math:`q_4`. The hand is located at point :math:`P_4`. The hand will be
constrained to the lever by :math:`\mathbf{r}^{P_4/O} = \mathbf{r}^{P_1/O}`.
The lever, upper arm, and lower arm will be modeled as thin cylinders for
inertial simplicity.

To begin, we will introduce two muscles that represent the bicep and the
tricep. A circular arc of radius :math:`r` is defined with its center at
:math:`P_3` and normal to :math:`\hat{c}_y`. Two muscle attachment points
:math:`C_m` and :math:`D_m` are fixed on the upper arm and lower arm,
respectively. The bicep muscle will act along a linear path from :math:`C_m` to
:math:`D_m`. The tricep will wrap around the circular arc and also attach at
the same points as the bicep.

::

   import sympy as sm
   import sympy.physics.mechanics as me
   import sympy.physics.biomechanics as bm

Define variables
----------------

Introduce the four coordinates :math:`\mathbf{q} = [q_1, q_2, q_3, q_4]^T` for
the lever angle, shoulder extension, shoulder rotation, and elbow extension. We
will also need generalized speeds :math:`\mathbf{u} = [u_1,u_2,u_3,u_4]^T`
which we define as :math:`\mathbf{u} = \dot{\mathbf{q}}`.

::

   q1, q2, q3, q4 = me.dynamicsymbols('q1, q2, q3, q4')
   u1, u2, u3, u4 = me.dynamicsymbols('u1, u2, u3, u4')
   u1d, u2d, u3d, u4d = me.dynamicsymbols('u1d, u2d, u3d, u4d')

The necessary constant parameters are:

- :math:`d_x, l_A`: locates :math:`P_1` from :math:`O` along the
  :math:`\hat{n}_x` and :math:`\hat{a}_y`, respectively
- :math:`d_y, d_z`: locates :math:`P_2` from :math:`O` along the :math:`N` unit
  vector directions
- :math:`l_C,l_D` : length of upper and lower arm
- :math:`m_A,m_C,m_D` : mass of lever, upper arm, and lower arm
- :math:`g` : acceleration due to gravity
- :math:`k` : lever linear rotational spring coefficient
- :math:`c` : leverl linear rotational damper coefficient

::

   dx, dy, dz, lA, lC, lD = sm.symbols('dx, dy, dz, lA, lC, lD', real=True,
                                       nonnegative=True)
   mA, mC, mD = sm.symbols('mA, mC, mD', real=True, positive=True)
   g, k, c, r = sm.symbols('g, k, c, r', real=True, positive=True)

::

   q = sm.Matrix([q1, q2, q3, q4])
   u = sm.Matrix([u1, u2, u3, u4])
   ud = u.diff(me.dynamicsymbols._t)
   ud_zerod = {udi: 0 for udi in ud}
   p = sm.Matrix([
       dx,
       dy,
       dz,
       lA,
       lC,
       lD,
       mA,
       mC,
       mD,
       g,
       kA,
       cA,
       r,
   ])

Define kinematics
-----------------

Define all the reference frames and points show in
:numref:`fig-biomechanics-steerer`. :math:`C_o` and :math:`D_o` are the mass
centers of the upper and lower arm, respectively.

::

   N, A, B, C, D = sm.symbols('N, A, B, C, D', cls=me.ReferenceFrame)
   O, P1, P2, P3, P4 = sm.symbols('O, P1, P2, P3, P4 ', cls=me.Point)
   Ao, Co, Cm, Dm, Do = sm.symbols('Ao, Co, Cm, Dm, Do', cls=me.Point)

The orientations and angular velocities of the reference frames are::

   A.orient_axis(N, q1, N.z)
   B.orient_axis(N, q2, N.y)
   C.orient_axis(B, q3, B.z)
   D.orient_axis(C, q4, C.y)
   A.set_ang_vel(N, u1*N.z)
   B.set_ang_vel(N, u2*N.y)
   C.set_ang_vel(B, u3*B.z)
   D.set_ang_vel(C, u4*C.y)

All of the points locations and velocities are::

   P1.set_pos(O, dx*N.x + lA*A.y)
   P2.set_pos(O, dy*N.y + dz*N.z)
   Co.set_pos(P2, lC/2*C.z)
   Cm.set_pos(P2, 2*lC/3*C.z)
   P3.set_pos(P2, lC*C.z)
   Dm.set_pos(P3, 1*lD/3*D.z)
   Do.set_pos(P3, lD/2*D.z)
   P4.set_pos(P3, lD*D.z)

   O.set_vel(N, 0)
   Ao.set_vel(N, 0)
   P1.v2pt_theory(O, N, A)
   P2.set_vel(N, 0)
   Co.v2pt_theory(P2, N, C)
   Cm.v2pt_theory(P2, N, C)
   P3.v2pt_theory(P2, N, C)
   Dm.v2pt_theory(P3, N, D)
   Do.v2pt_theory(P3, N, D)
   P4.v2pt_theory(P3, N, D)

There are three holonomic constrain equations needed to keep the hand
:math:`P_4` on the lever :math:`P_1`::

   holonomic = (P4.pos_from(O) - P1.pos_from(O)).to_matrix(N)

Define inertia
--------------

The inertia dyadics can be formed assuming the lever, upper arm, and lower arm
are thin cylinders::

   IA = me.Inertia(me.inertia(A, mA/12*lA**2, mA/2*lA**2, mA/12*lA**2), Ao)
   IC = me.Inertia(me.inertia(C, mC/12*lC**2, mC/12*lC**2, mC/2*lC**2), Co)
   ID = me.Inertia(me.inertia(D, mD/12*lD**2, mD/12*lD**2, mD/2*lD**2), Do)

   lever = me.RigidBody('lever', masscenter=Ao, frame=A, mass=mA, inertia=IA)
   u_arm = me.RigidBody('upper arm', masscenter=Co, frame=C, mass=mC, inertia=IC)
   l_arm = me.RigidBody('lower arm', masscenter=Do, frame=D, mass=mD, inertia=ID)

Define forces
-------------

::

   lever_resistance = me.Torque(A, (-kA*q1 - cA*u2)*N.z)
   lever_resistance = me.Torque(A, 0*N.z)

   gravC = me.Force(u_arm, mC*g*N.z)
   gravD = me.Force(l_arm, mD*g*N.z)

Bicep
~~~~~

::

   bicep_pathway = LinearPathway(Cm, Dm)
   bicep_activation = FirstOrderActivationDeGroote2016.with_default_constants('bicep')
   bicep = MusculotendonDeGroote2016('bicep', bicep_pathway, activation_dynamics=bicep_activation)
   bicep_constants = {
       bicep._F_M_max: 500.0,
       bicep._l_M_opt: 0.6 * 0.3,
       bicep._l_T_slack: 0.55 * 0.3,
       bicep._v_M_max: 10.0,
       bicep._alpha_opt: 0.0,
       bicep._beta: 0.1,
   }

Tricep
~~~~~~

::

   tricep_pathway = ExtensorPathway(C.y, P3, -C.z, D.z, Cm, Dm, r, q4)
   tricep_activation = FirstOrderActivationDeGroote2016.with_default_constants('tricep')
   tricep = MusculotendonDeGroote2016('tricep', tricep_pathway, activation_dynamics=tricep_activation)
   tricep_constants = {
       tricep._F_M_max: 500.0,
       tricep._l_M_opt: 0.6 * 0.3,
       tricep._l_T_slack: 0.65 * 0.3,
       tricep._v_M_max: 10.0,
       tricep._alpha_opt: 0.0,
       tricep._beta: 0.1,
   }

::

   loads = (
       bicep.to_loads() +
       tricep.to_loads() +
       [steer_resistance, gravA, gravB]
   )

Muscle Differential Equations
=============================

::

   musculotendon_constants = {**bicep_constants, **tricep_constants}
   mt = sm.Matrix(list(musculotendon_constants.keys()))

   a = list(bicep.activation_dynamics.state_variables) + list(tricep.activation_dynamics.state_variables)
   e = list(bicep.activation_dynamics.control_variables) + list(tricep.activation_dynamics.control_variables)
   da = list(bicep.activation_dynamics.state_equations.values()) + list(tricep.activation_dynamics.state_equations.values())
   eval_da = sm.lambdify((e, a), da, cse=True)


Equations of Motion
===================

::

   kane = me.KanesMethod(
       N,
       (q1,),
       (u1,),
       kd_eqs=(
           u1 - q1.diff(),
           u2 - q2.diff(),
           u3 - q3.diff(),
           u4 - q4.diff(),
       ),
       q_dependent=(q2, q3, q4),
       configuration_constraints=holonomic,
       velocity_constraints=holonomic.diff(t),
       u_dependent=(u2, u3, u4),
       bodies=(lever, u_arm, l_arm),
       forcelist=loads,
   )
