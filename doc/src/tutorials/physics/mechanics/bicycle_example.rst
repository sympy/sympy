=========================================
Linearized Carvallo-Whipple Bicycle Model
=========================================

The bicycle is an interesting system in that it can be modeled with multiple
rigid bodies, non-holonomic constraints, and a holonomic constraint. The
linearized equations of motion of the Carvallo-Whipple bicycle model are
presented and benchmarked in [Meijaard2007]_. This example will construct the
same linear equations of motion using :mod:`sympy.physics.mechanics`. ::

  >>> import sympy as sm
  >>> import sympy.physics.mechanics as me
  >>> me.mechanics_printing(pretty_print=False)

Declaration of Coordinates & Speeds
===================================

The simple definition of :math:`\mathbf{u} = \dot{\mathbf{q}}` is used in this model. The generalized speeds are:

- yaw frame angular rate :math:`u_1`,
- roll frame angular rate :math:`u_2`,
- rear wheel frame angular rate (spinning motion) :math:`u_3`,
- frame angular rate (pitching motion) :math:`u_4`,
- steering frame angular rate :math:`u_5`, and
- front wheel angular rate (spinning motion) :math:`u_6`.

Wheel positions are ignorable coordinates, so they are not introduced. ::

  >>> q1, q2, q3, q4, q5 = me.dynamicsymbols('q1 q2 q3 q4 q5')
  >>> q1d, q2d, q4d, q5d = me.dynamicsymbols('q1 q2 q4 q5', 1)
  >>> u1, u2, u3, u4, u5, u6 = me.dynamicsymbols('u1 u2 u3 u4 u5 u6')
  >>> u1d, u2d, u3d, u4d, u5d, u6d = me.dynamicsymbols('u1 u2 u3 u4 u5 u6', 1)

Declaration of System's Parameters
==================================

The constant parameters of the model are::

  >>> WFrad, WRrad, htangle, forkoffset = sm.symbols('WFrad WRrad htangle forkoffset')
  >>> forklength, framelength, forkcg1 = sm.symbols('forklength framelength forkcg1')
  >>> forkcg3, framecg1, framecg3, Iwr11 = sm.symbols('forkcg3 framecg1 framecg3 Iwr11')
  >>> Iwr22, Iwf11, Iwf22, Iframe11 = sm.symbols('Iwr22 Iwf11 Iwf22 Iframe11')
  >>> Iframe22, Iframe33, Iframe31, Ifork11 = sm.symbols('Iframe22 Iframe33 Iframe31 Ifork11')
  >>> Ifork22, Ifork33, Ifork31, g = sm.symbols('Ifork22 Ifork33 Ifork31 g')
  >>> mframe, mfork, mwf, mwr = sm.symbols('mframe mfork mwf mwr')

Kinematics of the Bicycle
=========================

Set up reference frames for the system
--------------------------------------

- ``N`` - inertial
- ``Y`` - yaw
- ``R`` - roll
- ``WR`` - rear wheel, rotation angle is ignorable coordinate so not oriented
- ``Frame`` - bicycle frame
- ``TempFrame`` - statically rotated frame for easier reference inertia definition
- ``Fork`` - bicycle fork
- ``TempFork`` - statically rotated frame for easier reference inertia definition
- ``WF`` - front wheel, again posses an ignorable coordinate

::

  >>> N = me.ReferenceFrame('N')
  >>> Y = N.orientnew('Y', 'Axis', [q1, N.z])
  >>> R = Y.orientnew('R', 'Axis', [q2, Y.x])
  >>> Frame = R.orientnew('Frame', 'Axis', [q4 + htangle, R.y])
  >>> WR = me.ReferenceFrame('WR')
  >>> TempFrame = Frame.orientnew('TempFrame', 'Axis', [-htangle, Frame.y])
  >>> Fork = Frame.orientnew('Fork', 'Axis', [q5, Frame.x])
  >>> TempFork = Fork.orientnew('TempFork', 'Axis', [-htangle, Fork.y])
  >>> WF = me.ReferenceFrame('WF')

Define relevant points for the system
-------------------------------------

``WR_cont`` - rear wheel contact
``WR_mc``- rear wheel's center of mass
``Steer`` - frame/fork connection
``Frame_mc`` - frame's center of mass
``Fork_mc`` - fork's center of mass
``WF_mc`` - front wheel's center of mass
``WF_cont`` - front wheel contact point

  >>> WR_cont = me.Point('WR_cont')
  >>> WR_mc = WR_cont.locatenew('WR_mc', WRrad*R.z)
  >>> Steer = WR_mc.locatenew('Steer', framelength*Frame.z)
  >>> Frame_mc = WR_mc.locatenew('Frame_mc', -framecg1*Frame.x + framecg3*Frame.z)
  >>> Fork_mc = Steer.locatenew('Fork_mc', -forkcg1*Fork.x + forkcg3*Fork.z)
  >>> WF_mc = Steer.locatenew('WF_mc', forklength*Fork.x + forkoffset*Fork.z)
  >>> WF_cont = WF_mc.locatenew('WF_cont', WFrad*(me.dot(Fork.y, Y.z)*Fork.y - Y.z).normalize())

Set the angular velocity of each frame
--------------------------------------

Angular accelerations end up being calculated automatically by differentiating
the angular velocities when first needed.

- ``u1`` is yaw rate
- ``u2`` is roll rate
- ``u3`` is rear wheel rate
- ``u4`` is frame pitch rate
- ``u5`` is fork steer rate
- ``u6`` is front wheel rate

::

  >>> Y.set_ang_vel(N, u1 * Y.z)
  >>> R.set_ang_vel(Y, u2 * R.x)
  >>> WR.set_ang_vel(Frame, u3 * Frame.y)
  >>> Frame.set_ang_vel(R, u4 * Frame.y)
  >>> Fork.set_ang_vel(Frame, u5 * Fork.x)
  >>> WF.set_ang_vel(Fork, u6 * Fork.y)

Form the velocities of the points, using the 2-point theorem. Accelerations
again are calculated automatically when first needed. ::

  >>> WR_cont.set_vel(N, 0)
  >>> WR_mc.v2pt_theory(WR_cont, N, WR)
  WRrad*(u1*sin(q2) + u3 + u4)*R.x - WRrad*u2*R.y
  >>> Steer.v2pt_theory(WR_mc, N, Frame)
  WRrad*(u1*sin(q2) + u3 + u4)*R.x - WRrad*u2*R.y + framelength*(u1*sin(q2) + u4)*Frame.x - framelength*(-u1*sin(htangle + q4)*cos(q2) + u2*cos(htangle + q4))*Frame.y
  >>> Frame_mc.v2pt_theory(WR_mc, N, Frame)
  WRrad*(u1*sin(q2) + u3 + u4)*R.x - WRrad*u2*R.y + framecg3*(u1*sin(q2) + u4)*Frame.x + (-framecg1*(u1*cos(htangle + q4)*cos(q2) + u2*sin(htangle + q4)) - framecg3*(-u1*sin(htangle + q4)*cos(q2) + u2*cos(htangle + q4)))*Frame.y + framecg1*(u1*sin(q2) + u4)*Frame.z
  >>> Fork_mc.v2pt_theory(Steer, N, Fork)
  WRrad*(u1*sin(q2) + u3 + u4)*R.x - WRrad*u2*R.y + framelength*(u1*sin(q2) + u4)*Frame.x - framelength*(-u1*sin(htangle + q4)*cos(q2) + u2*cos(htangle + q4))*Frame.y + forkcg3*((sin(q2)*cos(q5) + sin(q5)*cos(htangle + q4)*cos(q2))*u1 + u2*sin(htangle + q4)*sin(q5) + u4*cos(q5))*Fork.x + (-forkcg1*((-sin(q2)*sin(q5) + cos(htangle + q4)*cos(q2)*cos(q5))*u1 + u2*sin(htangle + q4)*cos(q5) - u4*sin(q5)) - forkcg3*(-u1*sin(htangle + q4)*cos(q2) + u2*cos(htangle + q4) + u5))*Fork.y + forkcg1*((sin(q2)*cos(q5) + sin(q5)*cos(htangle + q4)*cos(q2))*u1 + u2*sin(htangle + q4)*sin(q5) + u4*cos(q5))*Fork.z
  >>> WF_mc.v2pt_theory(Steer, N, Fork)
  WRrad*(u1*sin(q2) + u3 + u4)*R.x - WRrad*u2*R.y + framelength*(u1*sin(q2) + u4)*Frame.x - framelength*(-u1*sin(htangle + q4)*cos(q2) + u2*cos(htangle + q4))*Frame.y + forkoffset*((sin(q2)*cos(q5) + sin(q5)*cos(htangle + q4)*cos(q2))*u1 + u2*sin(htangle + q4)*sin(q5) + u4*cos(q5))*Fork.x + (forklength*((-sin(q2)*sin(q5) + cos(htangle + q4)*cos(q2)*cos(q5))*u1 + u2*sin(htangle + q4)*cos(q5) - u4*sin(q5)) - forkoffset*(-u1*sin(htangle + q4)*cos(q2) + u2*cos(htangle + q4) + u5))*Fork.y - forklength*((sin(q2)*cos(q5) + sin(q5)*cos(htangle + q4)*cos(q2))*u1 + u2*sin(htangle + q4)*sin(q5) + u4*cos(q5))*Fork.z
  >>> WF_cont.v2pt_theory(WF_mc, N, WF)
  - WFrad*((-sin(q2)*sin(q5)*cos(htangle + q4) + cos(q2)*cos(q5))*u6 + u4*cos(q2) + u5*sin(htangle + q4)*sin(q2))/sqrt((-sin(q2)*cos(q5) - sin(q5)*cos(htangle + q4)*cos(q2))*(sin(q2)*cos(q5) + sin(q5)*cos(htangle + q4)*cos(q2)) + 1)*Y.x + WFrad*(u2 + u5*cos(htangle + q4) + u6*sin(htangle + q4)*sin(q5))/sqrt((-sin(q2)*cos(q5) - sin(q5)*cos(htangle + q4)*cos(q2))*(sin(q2)*cos(q5) + sin(q5)*cos(htangle + q4)*cos(q2)) + 1)*Y.y + WRrad*(u1*sin(q2) + u3 + u4)*R.x - WRrad*u2*R.y + framelength*(u1*sin(q2) + u4)*Frame.x - framelength*(-u1*sin(htangle + q4)*cos(q2) + u2*cos(htangle + q4))*Frame.y + (-WFrad*(sin(q2)*cos(q5) + sin(q5)*cos(htangle + q4)*cos(q2))*((-sin(q2)*sin(q5) + cos(htangle + q4)*cos(q2)*cos(q5))*u1 + u2*sin(htangle + q4)*cos(q5) - u4*sin(q5))/sqrt((-sin(q2)*cos(q5) - sin(q5)*cos(htangle + q4)*cos(q2))*(sin(q2)*cos(q5) + sin(q5)*cos(htangle + q4)*cos(q2)) + 1) + forkoffset*((sin(q2)*cos(q5) + sin(q5)*cos(htangle + q4)*cos(q2))*u1 + u2*sin(htangle + q4)*sin(q5) + u4*cos(q5)))*Fork.x + (forklength*((-sin(q2)*sin(q5) + cos(htangle + q4)*cos(q2)*cos(q5))*u1 + u2*sin(htangle + q4)*cos(q5) - u4*sin(q5)) - forkoffset*(-u1*sin(htangle + q4)*cos(q2) + u2*cos(htangle + q4) + u5))*Fork.y + (WFrad*(sin(q2)*cos(q5) + sin(q5)*cos(htangle + q4)*cos(q2))*(-u1*sin(htangle + q4)*cos(q2) + u2*cos(htangle + q4) + u5)/sqrt((-sin(q2)*cos(q5) - sin(q5)*cos(htangle + q4)*cos(q2))*(sin(q2)*cos(q5) + sin(q5)*cos(htangle + q4)*cos(q2)) + 1) - forklength*((sin(q2)*cos(q5) + sin(q5)*cos(htangle + q4)*cos(q2))*u1 + u2*sin(htangle + q4)*sin(q5) + u4*cos(q5)))*Fork.z

The kinematic differential equations are as follows. Each entry in this list is
equal to zero. ::

  >>> kd = [q1d - u1, q2d - u2, q4d - u4, q5d - u5]

Setup the constraints
---------------------

The nonholonomic constraints are the velocity of the front wheel contact point
dotted into the X, Y, and Z directions; the yaw frame is used as it is "closer"
to the front wheel (one fewer direction cosine matrix connecting them). These
constraints force the velocity of the front wheel contact point to be zero in
the inertial frame; the X and Y direction constraints enforce a "no-slip"
condition, and the Z direction constraint forces the front wheel contact point
to not move away from the ground frame, essentially replicating the holonomic
constraint which does not allow the frame pitch to change in an invalid
fashion. ::

  >>> conlist_speed = [me.dot(WF_cont.vel(N), Y.x),
  ...                  me.dot(WF_cont.vel(N), Y.y),
  ...                  me.dot(WF_cont.vel(N), Y.z)]

The holonomic constraint is that the position from the rear wheel contact point
to the front wheel contact point when dotted into the normal-to-ground plane
direction must be zero; effectively that the front and rear wheel contact
points are always touching the ground plane. This is actually not part of the
dynamical differential equations, but is necessary for the linearization
process. ::

  >>> conlist_coord = [me.dot(WF_cont.pos_from(WR_cont), Y.z)]

Inertia and Rigid Bodies
========================

Sets the inertias of each body. Uses the inertia frame to construct the inertia
dyadics. Wheel inertias are only defined by principal moments of inertia, and
are in fact constant in the frame and fork reference frames; it is for this
reason that the orientations of the wheels does not need to be defined. The
frame and fork inertias are defined in the 'Temp' frames which are fixed to the
appropriate body frames; this is to allow easier input of the reference values
of the benchmark paper. Note that due to slightly different orientations, the
products of inertia need to have their signs flipped; this is done later when
entering the numerical value. ::

  >>> Frame_I = (me.inertia(TempFrame, Iframe11, Iframe22, Iframe33, 0, 0,
  ...                       Iframe31), Frame_mc)
  >>> Fork_I = (me.inertia(TempFork, Ifork11, Ifork22, Ifork33, 0, 0, Ifork31), Fork_mc)
  >>> WR_I = (me.inertia(Frame, Iwr11, Iwr22, Iwr11), WR_mc)
  >>> WF_I = (me.inertia(Fork, Iwf11, Iwf22, Iwf11), WF_mc)

Declaration of the ``RigidBody`` containers. ::

  >>> BodyFrame = me.RigidBody('BodyFrame', Frame_mc, Frame, mframe, Frame_I)
  >>> BodyFork = me.RigidBody('BodyFork', Fork_mc, Fork, mfork, Fork_I)
  >>> BodyWR = me.RigidBody('BodyWR', WR_mc, WR, mwr, WR_I)
  >>> BodyWF = me.RigidBody('BodyWF', WF_mc, WF, mwf, WF_I)
  >>> bodies = [BodyFrame, BodyFork, BodyWR, BodyWF]

Gravitational Loads
===================

The force list; each body has the appropriate gravitational force applied at
its center of mass. ::

  >>> forces = [(Frame_mc, -mframe * g * Y.z),
  ...           (Fork_mc, -mfork * g * Y.z),
  ...           (WF_mc, -mwf * g * Y.z),
  ...           (WR_mc, -mwr * g * Y.z)]
  ...

Nonlinear Equations of Motion
=============================

The ``N`` frame is the inertial frame, coordinates are supplied in the order of
independent, dependent coordinates. The kinematic differential equations are
also entered here. Here the independent speeds are specified, followed by the
dependent speeds, along with the non-holonomic constraints. The dependent
coordinate is also provided, with the holonomic constraint. Again, this is only
comes into play in the linearization process, but is necessary for the
linearization to correctly work. ::

  >>> kane = me.KanesMethod(
  ...     N,
  ...     q_ind=[q1, q2, q5],
  ...     q_dependent=[q4],
  ...     configuration_constraints=conlist_coord,
  ...     u_ind=[u2, u3, u5],
  ...     u_dependent=[u1, u4, u6],
  ...     velocity_constraints=conlist_speed,
  ...     kd_eqs=kd,
  ...     constraint_solver='CRAMER')
  >>> fr, frstar = kane.kanes_equations(bodies, loads=forces)

Linearized Equations of Motion
==============================

This is the start of entering in the numerical values from the benchmark paper
to validate the eigenvalues of the linearized equations from this model to the
reference eigenvalues. Look at the aforementioned paper for more information.
Some of these are intermediate values, used to transform values from the paper
into the coordinate systems used in this model. ::

  >>> PaperRadRear  =  0.3
  >>> PaperRadFront =  0.35
  >>> HTA           =  sm.evalf.N(sm.pi/2 - sm.pi/10)
  >>> TrailPaper    =  0.08
  >>> rake          =  sm.evalf.N(-(TrailPaper*sm.sin(HTA) - (PaperRadFront*sm.cos(HTA))))
  >>> PaperWb       =  1.02
  >>> PaperFrameCgX =  0.3
  >>> PaperFrameCgZ =  0.9
  >>> PaperForkCgX  =  0.9
  >>> PaperForkCgZ  =  0.7
  >>> FrameLength   =  sm.evalf.N(PaperWb*sm.sin(HTA) - (rake -
  ...                             (PaperRadFront - PaperRadRear)*sm.cos(HTA)))
  >>> FrameCGNorm   =  sm.evalf.N((PaperFrameCgZ - PaperRadRear -
  ...                             (PaperFrameCgX/sm.sin(HTA))*sm.cos(HTA))*sm.sin(HTA))
  >>> FrameCGPar    =  sm.evalf.N((PaperFrameCgX / sm.sin(HTA) +
  ...                             (PaperFrameCgZ - PaperRadRear -
  ...                              PaperFrameCgX / sm.sin(HTA)*sm.cos(HTA))*sm.cos(HTA)))
  >>> tempa         =  sm.evalf.N((PaperForkCgZ - PaperRadFront))
  >>> tempb         =  sm.evalf.N((PaperWb-PaperForkCgX))
  >>> tempc         =  sm.evalf.N(sm.sqrt(tempa**2 + tempb**2))
  >>> PaperForkL    =  sm.evalf.N((PaperWb*sm.cos(HTA) -
  ...                             (PaperRadFront - PaperRadRear)*sm.sin(HTA)))
  >>> ForkCGNorm    =  sm.evalf.N(rake + (tempc*sm.sin(sm.pi/2 -
  ...                             HTA - sm.acos(tempa/tempc))))
  >>> ForkCGPar     =  sm.evalf.N(tempc*sm.cos((sm.pi/2 - HTA) -
  ...                             sm.acos(tempa/tempc)) - PaperForkL)

Here is the final assembly of the numerical values. The symbol 'v' is the
forward speed of the bicycle (a concept which only makes sense in the upright,
static equilibrium case?). These are in a dictionary which will later be
substituted in. Again the sign on the *product* of inertia values is flipped
here, due to different orientations of coordinate systems. ::

  >>> v = sm.Symbol('v')
  >>> val_dict = {
  ...     WFrad: PaperRadFront,
  ...     WRrad: PaperRadRear,
  ...     htangle: HTA,
  ...     forkoffset: rake,
  ...     forklength: PaperForkL,
  ...     framelength: FrameLength,
  ...     forkcg1: ForkCGPar,
  ...     forkcg3: ForkCGNorm,
  ...     framecg1: FrameCGNorm,
  ...     framecg3: FrameCGPar,
  ...     Iwr11: 0.0603,
  ...     Iwr22: 0.12,
  ...     Iwf11: 0.1405,
  ...     Iwf22: 0.28,
  ...     Ifork11: 0.05892,
  ...     Ifork22: 0.06,
  ...     Ifork33: 0.00708,
  ...     Ifork31: 0.00756,
  ...     Iframe11: 9.2,
  ...     Iframe22: 11,
  ...     Iframe33: 2.8,
  ...     Iframe31: -2.4,
  ...     mfork: 4,
  ...     mframe: 85,
  ...     mwf: 3,
  ...     mwr: 2,
  ...     g: 9.81,
  ... }
  ...

Linearize the equations of motion about the equilibrium point::

  >>> eq_point = {
  ...     u1d: 0,
  ...     u2d: 0,
  ...     u3d: 0,
  ...     u4d: 0,
  ...     u5d: 0,
  ...     u6d: 0,
  ...     q1: 0,
  ...     q2: 0,
  ...     q4: 0,
  ...     q5: 0,
  ...     u1: 0,
  ...     u2: 0,
  ...     u3: v/PaperRadRear,
  ...     u4: 0,
  ...     u5: 0,
  ...     u6: v/PaperRadFront,
  ... }
  ...
  >>> Amat, _, _ = kane.linearize(A_and_B=True, op_point=eq_point, linear_solver='CRAMER')
  >>> Amat = me.msubs(Amat, val_dict)

Calculate the Eigenvalues
-------------------------

Finally, we construct an "A" matrix for the form :math:`\dot{\mathbf{x}} =
\mathbf{A} \mathbf{x}` (:math:`\mathbf{x}` being the state vector, although in
this case, the sizes are a little off). The following line extracts only the
minimum entries required for eigenvalue analysis, which correspond to rows and
columns for lean, steer, lean rate, and steer rate.

::

  >>> A = Amat.extract([1, 2, 3, 5], [1, 2, 3, 5])
  >>> A
  Matrix([
  [               0,                                           0,                    1,                    0],
  [               0,                                           0,                    0,                    1],
  [9.48977444677355, -0.891197738059089*v**2 - 0.571523173729245, -0.105522449805691*v, -0.330515398992311*v],
  [11.7194768719633,    30.9087533932407 - 1.97171508499972*v**2,   3.67680523332152*v,  -3.08486552743311*v]])
  >>> print('v = 1')
  v = 1
  >>> print(A.subs(v, 1).eigenvals())
  {-3.13423125066578 - 1.05503732448615e-65*I: 1, 3.52696170990069 - 0.807740275199311*I: 1, 3.52696170990069 + 0.807740275199311*I: 1, -7.11008014637441: 1}
  >>> print('v = 2')
  v = 2
  >>> print(A.subs(v, 2).eigenvals())
  {2.68234517512745 - 1.68066296590676*I: 1, 2.68234517512745 + 1.68066296590676*I: 1, -3.07158645641514: 1, -8.67387984831737: 1}
  >>> print('v = 3')
  v = 3
  >>> print(A.subs(v, 3).eigenvals())
  {1.70675605663973 - 2.31582447384324*I: 1, 1.70675605663973 + 2.31582447384324*I: 1, -2.63366137253665: 1, -10.3510146724592: 1}
  >>> print('v = 4')
  v = 4
  >>> print(A.subs(v, 4).eigenvals())
  {0.413253315211239 - 3.07910818603205*I: 1, 0.413253315211239 + 3.07910818603205*I: 1, -1.42944427361326 + 1.65070329233125e-64*I: 1, -12.1586142657644: 1}
  >>> print('v = 5')
  v = 5
  >>> print(A.subs(v, 5).eigenvals())
  {-0.775341882195845 - 4.46486771378823*I: 1, -0.322866429004087 + 3.32140410564766e-64*I: 1, -0.775341882195845 + 4.46486771378823*I: 1, -14.0783896927982: 1}

The eigenvalues shown above match those in Table 2 on pg. 1971 of
[Meijaard2007]_. This concludes the bicycle example.

References
==========

.. [Meijaard2007] Meijaard, J. P., Papadopoulos, J. M., Ruina, A., & Schwab, A.
   L. (2007). Linearized dynamics equations for the balance and steer of a
   bicycle: A benchmark and review. Proceedings of the Royal Society A:
   Mathematical, Physical and Engineering Sciences, 463(2084), 1955–1982.
   https://doi.org/10.1098/rspa.2007.1857
