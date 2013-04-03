==============================
Examples for Physics.Mechanics
==============================

This module has been developed to aid in formulation of equations of motion;
here, examples follow to help display how to use this module.

The Rolling Disc
================

The first example is that of a rolling disc. The disc is assumed to be
infinitely thin, in contact with the ground at only 1 point, and it is rolling
without slip on the ground. See the image below.

.. image:: ex_rd.*
   :height: 550
   :width: 550
   :align: center

Here the definition of the rolling disc's kinematics is formed from the contact
point up, removing the need to introduce generalized speeds. Only 3
configuration and three speed variables are need to describe this system, along
with the disc's mass and radius, and the local gravity (note that mass will
drop out). ::

  >>> from sympy import symbols, sin, cos, tan
  >>> from sympy.physics.mechanics import *
  >>> q1, q2, q3, u1, u2, u3  = dynamicsymbols('q1 q2 q3 u1 u2 u3')
  >>> q1d, q2d, q3d, u1d, u2d, u3d = dynamicsymbols('q1 q2 q3 u1 u2 u3', 1)
  >>> r, m, g = symbols('r m g')
  >>> mechanics_printing()

The kinematics are formed by a series of simple rotations. Each simple rotation
creates a new frame, and the next rotation is defined by the new frame's basis
vectors. This example uses a 3-1-2 series of rotations, or Z, X, Y series of
rotations. Angular velocity for this is defined using the second frame's basis
(the lean frame); it is for this reason that we defined intermediate frames,
rather than using a body-three orientation. ::

  >>> N = ReferenceFrame('N')
  >>> Y = N.orientnew('Y', 'Axis', [q1, N.z])
  >>> L = Y.orientnew('L', 'Axis', [q2, Y.x])
  >>> R = L.orientnew('R', 'Axis', [q3, L.y])
  >>> R.set_ang_vel(N, u1 * L.x + u2 * L.y + u3 * L.z)
  >>> R.set_ang_acc(N, R.ang_vel_in(N).dt(R) + (R.ang_vel_in(N) ^ R.ang_vel_in(N)))

This is the translational kinematics. We create a point with no velocity
in N; this is the contact point between the disc and ground. Next we form
the position vector from the contact point to the disc's center of mass.
Finally we form the velocity and acceleration of the disc. ::

  >>> C = Point('C')
  >>> C.set_vel(N, 0)
  >>> Dmc = C.locatenew('Dmc', r * L.z)
  >>> Dmc.v2pt_theory(C, N, R)
  r*u2*L.x - r*u1*L.y
  >>> Dmc.a2pt_theory(C, N, R)
  (r*u1*u3 + r*u2')*L.x + (-r*(-u3*q3' + u1') + r*u2*u3)*L.y + (-r*u1**2 - r*u2**2)*L.z

This is a simple way to form the inertia dyadic. The inertia of the disc does
not change within the lean frame as the disc rolls; this will make for simpler
equations in the end. ::

  >>> I = inertia(L, m / 4 * r**2, m / 2 * r**2, m / 4 * r**2)
  >>> mprint(I)
  m*r**2/4*(L.x|L.x) + m*r**2/2*(L.y|L.y) + m*r**2/4*(L.z|L.z)

Kinematic differential equations; how the generalized coordinate time
derivatives relate to generalized speeds. Here these were computed by hand. ::

  >>> kd = [q1d - u3/cos(q2), q2d - u1, q3d - u2 + u3 * tan(q2)]

Creation of the force list; it is the gravitational force at the center of mass of
the disc. Then we create the disc by assigning a Point to the center of mass
attribute, a ReferenceFrame to the frame attribute, and mass and inertia. Then
we form the body list. ::

  >>> ForceList = [(Dmc, - m * g * Y.z)]
  >>> BodyD = RigidBody('BodyD', Dmc, R, m, (I, Dmc))
  >>> BodyList = [BodyD]

Finally we form the equations of motion, using the same steps we did before.
Specify inertial frame, supply generalized coordinates and speeds, supply
kinematic differential equation dictionary, compute Fr from the force list and
Fr* from the body list, compute the mass matrix and forcing terms, then solve
for the u dots (time derivatives of the generalized speeds). ::

  >>> KM = KanesMethod(N, q_ind=[q1, q2, q3], u_ind=[u1, u2, u3], kd_eqs=kd)
  >>> (fr, frstar) = KM.kanes_equations(ForceList, BodyList)
  >>> MM = KM.mass_matrix
  >>> forcing = KM.forcing
  >>> rhs = MM.inv() * forcing
  >>> kdd = KM.kindiffdict()
  >>> rhs = rhs.subs(kdd)
  >>> rhs.simplify()
  >>> mprint(rhs)
  [(4*g*sin(q2)/5 + 2*r*u2*u3 - r*u3**2*tan(q2))/r]
  [                                     -2*u1*u3/3]
  [                        (-2*u2 + u3*tan(q2))*u1]

This concludes the rolling disc example.

The Rolling Disc, again
=======================

We will now revisit the rolling disc example, except this time we are bringing
the non-contributing (constraint) forces into evidence. See [Kane1985]_ for a
more thorough explanation of this. Here, we will turn on the automatic
simplifcation done when doing vector operations. It makes the outputs nicer for
small problems, but can cause larger vector operations to hang. ::

  >>> from sympy.physics.mechanics import *
  >>> mechanics_printing()
  >>> Vector.simp = True
  >>> q1, q2, q3, u1, u2, u3  = dynamicsymbols('q1 q2 q3 u1 u2 u3')
  >>> q1d, q2d, q3d, u1d, u2d, u3d = dynamicsymbols('q1 q2 q3 u1 u2 u3', 1)
  >>> r, m, g = symbols('r m g')

These two lines introduce the extra quantities needed to find the constraint
forces. ::

  >>> u4, u5, u6, f1, f2, f3 = dynamicsymbols('u4 u5 u6 f1 f2 f3')

Most of the main code is the same as before. ::

  >>> N = ReferenceFrame('N')
  >>> Y = N.orientnew('Y', 'Axis', [q1, N.z])
  >>> L = Y.orientnew('L', 'Axis', [q2, Y.x])
  >>> R = L.orientnew('R', 'Axis', [q3, L.y])
  >>> R.set_ang_vel(N, u1 * L.x + u2 * L.y + u3 * L.z)
  >>> R.set_ang_acc(N, R.ang_vel_in(N).dt(R) + (R.ang_vel_in(N) ^ R.ang_vel_in(N)))

The definition of rolling without slip necessitates that the velocity of the
contact point is zero; as part of bringing the constraint forces into evidence,
we have to introduce speeds at this point, which will by definition always be
zero. They are normal to the ground, along the path which the disc is rolling,
and along the ground in an perpendicular direction. ::

  >>> C = Point('C')
  >>> C.set_vel(N, u4 * L.x + u5 * (Y.z ^ L.x) + u6 * Y.z)
  >>> Dmc = C.locatenew('Dmc', r * L.z)
  >>> vel = Dmc.v2pt_theory(C, N, R)
  >>> acc = Dmc.a2pt_theory(C, N, R)
  >>> I = inertia(L, m / 4 * r**2, m / 2 * r**2, m / 4 * r**2)
  >>> kd = [q1d - u3/cos(q2), q2d - u1, q3d - u2 + u3 * tan(q2)]

Just as we previously introduced three speeds as part of this process, we also
introduce three forces; they are in the same direction as the speeds, and
represent the constraint forces in those directions. ::

  >>> ForceList = [(Dmc, - m * g * Y.z), (C, f1 * L.x + f2 * (Y.z ^ L.x) + f3 * Y.z)]
  >>> BodyD = RigidBody('BodyD', Dmc, R, m, (I, Dmc))
  >>> BodyList = [BodyD]

  >>> KM = KanesMethod(N, q_ind=[q1, q2, q3], u_ind=[u1, u2, u3], kd_eqs=kd,
  ...           u_auxiliary=[u4, u5, u6])
  >>> (fr, frstar) = KM.kanes_equations(ForceList, BodyList)
  >>> MM = KM.mass_matrix
  >>> forcing = KM.forcing
  >>> rhs = MM.inv() * forcing
  >>> kdd = KM.kindiffdict()
  >>> rhs = rhs.subs(kdd)
  >>> rhs.simplify()
  >>> mprint(rhs)
  [(4*g*sin(q2)/5 + 2*r*u2*u3 - r*u3**2*tan(q2))/r]
  [                                     -2*u1*u3/3]
  [                        (-2*u2 + u3*tan(q2))*u1]
  >>> from sympy import signsimp, collect, factor_terms
  >>> mprint(KM.auxiliary_eqs.applyfunc(lambda w: signsimp(collect(factor_terms(w), m*r))))
  [                                                   m*r*(u1*u3 + u2') - f1]
  [      m*r*((u1**2 + u2**2)*sin(q2) + (u2*u3 + u3*q3' - u1')*cos(q2)) - f2]
  [g*m - m*r*((u1**2 + u2**2)*cos(q2) - (u2*u3 + u3*q3' - u1')*sin(q2)) - f3]

The Bicycle
===========

The bicycle is an interesting system in that it has multiple rigid bodies,
non-holonomic constraints, and a holonomic constraint. The linearized equations
of motion are presented in [Meijaard2007]_. This example will go through
construction of the equations of motion in :mod:`mechanics`. ::

  >>> from sympy import *
  >>> from sympy.physics.mechanics import *
  >>> print('Calculation of Linearized Bicycle \"A\" Matrix, '
  ...       'with States: Roll, Steer, Roll Rate, Steer Rate')
  Calculation of Linearized Bicycle "A" Matrix, with States: Roll, Steer, Roll Rate, Steer Rate


Note that this code has been crudely ported from Autolev, which is the reason
for some of the unusual naming conventions. It was purposefully as similar as
possible in order to aid initial porting & debugging. We also turn off
Vector.simp (turned on in the last example) to avoid hangups when doing
computations in this problem. ::

  >>> Vector.simp = False
  >>> mechanics_printing()

Declaration of Coordinates & Speeds:
A Simple definitions for qdots: (qd = u) is used in this code.  Speeds are: yaw
frame ang. rate, roll frame ang. rate, rear wheel frame ang.  rate (spinning
motion), frame ang. rate (pitching motion), steering frame ang. rate, and front
wheel ang. rate (spinning motion).  Wheel positions are ignorable coordinates,
so they are not introduced. ::

  >>> q1, q2, q4, q5 = dynamicsymbols('q1 q2 q4 q5')
  >>> q1d, q2d, q4d, q5d = dynamicsymbols('q1 q2 q4 q5', 1)
  >>> u1, u2, u3, u4, u5, u6 = dynamicsymbols('u1 u2 u3 u4 u5 u6')
  >>> u1d, u2d, u3d, u4d, u5d, u6d = dynamicsymbols('u1 u2 u3 u4 u5 u6', 1)

Declaration of System's Parameters:
The below symbols should be fairly self-explanatory. ::

  >>> WFrad, WRrad, htangle, forkoffset = symbols('WFrad WRrad htangle forkoffset')
  >>> forklength, framelength, forkcg1 = symbols('forklength framelength forkcg1')
  >>> forkcg3, framecg1, framecg3, Iwr11 = symbols('forkcg3 framecg1 framecg3 Iwr11')
  >>> Iwr22, Iwf11, Iwf22, Iframe11 = symbols('Iwr22 Iwf11 Iwf22 Iframe11')
  >>> Iframe22, Iframe33, Iframe31, Ifork11 = \
  ...     symbols('Iframe22 Iframe33 Iframe31 Ifork11')
  >>> Ifork22, Ifork33, Ifork31, g = symbols('Ifork22 Ifork33 Ifork31 g')
  >>> mframe, mfork, mwf, mwr = symbols('mframe mfork mwf mwr')

Set up reference frames for the system:
N - inertial
Y - yaw
R - roll
WR - rear wheel, rotation angle is ignorable coordinate so not oriented
Frame - bicycle frame
TempFrame - statically rotated frame for easier reference inertia definition
Fork - bicycle fork
TempFork - statically rotated frame for easier reference inertia definition
WF - front wheel, again posses a ignorable coordinate ::

  >>> N = ReferenceFrame('N')
  >>> Y = N.orientnew('Y', 'Axis', [q1, N.z])
  >>> R = Y.orientnew('R', 'Axis', [q2, Y.x])
  >>> Frame = R.orientnew('Frame', 'Axis', [q4 + htangle, R.y])
  >>> WR = ReferenceFrame('WR')
  >>> TempFrame = Frame.orientnew('TempFrame', 'Axis', [-htangle, Frame.y])
  >>> Fork = Frame.orientnew('Fork', 'Axis', [q5, Frame.x])
  >>> TempFork = Fork.orientnew('TempFork', 'Axis', [-htangle, Fork.y])
  >>> WF = ReferenceFrame('WF')


Kinematics of the Bicycle:
First block of code is forming the positions of the relevant points rear wheel
contact -> rear wheel's center of mass -> frame's center of mass + frame/fork connection
-> fork's center of mass + front wheel's center of mass -> front wheel contact point. ::

  >>> WR_cont = Point('WR_cont')
  >>> WR_mc = WR_cont.locatenew('WR_mc', WRrad * R.z)
  >>> Steer = WR_mc.locatenew('Steer', framelength * Frame.z)
  >>> Frame_mc = WR_mc.locatenew('Frame_mc', -framecg1 * Frame.x + framecg3 * Frame.z)
  >>> Fork_mc = Steer.locatenew('Fork_mc', -forkcg1 * Fork.x + forkcg3 * Fork.z)
  >>> WF_mc = Steer.locatenew('WF_mc', forklength * Fork.x + forkoffset * Fork.z)
  >>> WF_cont = WF_mc.locatenew('WF_cont', WFrad*(dot(Fork.y, Y.z)*Fork.y - \
  ...                                             Y.z).normalize())

Set the angular velocity of each frame:
Angular accelerations end up being calculated automatically by differentiating
the angular velocities when first needed. ::
u1 is yaw rate
u2 is roll rate
u3 is rear wheel rate
u4 is frame pitch rate
u5 is fork steer rate
u6 is front wheel rate ::

  >>> Y.set_ang_vel(N, u1 * Y.z)
  >>> R.set_ang_vel(Y, u2 * R.x)
  >>> WR.set_ang_vel(Frame, u3 * Frame.y)
  >>> Frame.set_ang_vel(R, u4 * Frame.y)
  >>> Fork.set_ang_vel(Frame, u5 * Fork.x)
  >>> WF.set_ang_vel(Fork, u6 * Fork.y)

Form the velocities of the points, using the 2-point theorem.  Accelerations
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
  WRrad*(u1*sin(q2) + u3 + u4)*R.x - WRrad*u2*R.y + framelength*(u1*sin(q2) + u4)*Frame.x - framelength*(-u1*sin(htangle + q4)*cos(q2) + u2*cos(htangle + q4))*Frame.y + (-WFrad*(sin(q2)*cos(q5) + sin(q5)*cos(htangle + q4)*cos(q2))*((-sin(q2)*sin(q5) + cos(htangle + q4)*cos(q2)*cos(q5))*u1 + u2*sin(htangle + q4)*cos(q5) - u4*sin(q5))/sqrt((-sin(q2)*cos(q5) - sin(q5)*cos(htangle + q4)*cos(q2))*(sin(q2)*cos(q5) + sin(q5)*cos(htangle + q4)*cos(q2)) + 1) + forkoffset*((sin(q2)*cos(q5) + sin(q5)*cos(htangle + q4)*cos(q2))*u1 + u2*sin(htangle + q4)*sin(q5) + u4*cos(q5)))*Fork.x + (forklength*((-sin(q2)*sin(q5) + cos(htangle + q4)*cos(q2)*cos(q5))*u1 + u2*sin(htangle + q4)*cos(q5) - u4*sin(q5)) - forkoffset*(-u1*sin(htangle + q4)*cos(q2) + u2*cos(htangle + q4) + u5))*Fork.y + (WFrad*(sin(q2)*cos(q5) + sin(q5)*cos(htangle + q4)*cos(q2))*(-u1*sin(htangle + q4)*cos(q2) + u2*cos(htangle + q4) + u5)/sqrt((-sin(q2)*cos(q5) - sin(q5)*cos(htangle + q4)*cos(q2))*(sin(q2)*cos(q5) + sin(q5)*cos(htangle + q4)*cos(q2)) + 1) - forklength*((sin(q2)*cos(q5) + sin(q5)*cos(htangle + q4)*cos(q2))*u1 + u2*sin(htangle + q4)*sin(q5) + u4*cos(q5)))*Fork.z - WFrad*((-sin(q2)*sin(q5)*cos(htangle + q4) + cos(q2)*cos(q5))*u6 + u4*cos(q2) + u5*sin(htangle + q4)*sin(q2))/sqrt((-sin(q2)*cos(q5) - sin(q5)*cos(htangle + q4)*cos(q2))*(sin(q2)*cos(q5) + sin(q5)*cos(htangle + q4)*cos(q2)) + 1)*Y.x + WFrad*(u2 + u5*cos(htangle + q4) + u6*sin(htangle + q4)*sin(q5))/sqrt((-sin(q2)*cos(q5) - sin(q5)*cos(htangle + q4)*cos(q2))*(sin(q2)*cos(q5) + sin(q5)*cos(htangle + q4)*cos(q2)) + 1)*Y.y


Sets the inertias of each body. Uses the inertia frame to construct the inertia
dyadics. Wheel inertias are only defined by principal moments of inertia, and
are in fact constant in the frame and fork reference frames; it is for this
reason that the orientations of the wheels does not need to be defined. The
frame and fork inertias are defined in the 'Temp' frames which are fixed to the
appropriate body frames; this is to allow easier input of the reference values
of the benchmark paper. Note that due to slightly different orientations, the
products of inertia need to have their signs flipped; this is done later when
entering the numerical value. ::

  >>> Frame_I = (inertia(TempFrame, Iframe11, Iframe22, Iframe33, 0, 0,
  ...                                                   Iframe31), Frame_mc)
  >>> Fork_I = (inertia(TempFork, Ifork11, Ifork22, Ifork33, 0, 0, Ifork31), Fork_mc)
  >>> WR_I = (inertia(Frame, Iwr11, Iwr22, Iwr11), WR_mc)
  >>> WF_I = (inertia(Fork, Iwf11, Iwf22, Iwf11), WF_mc)

Declaration of the RigidBody containers. ::

  >>> BodyFrame = RigidBody('BodyFrame', Frame_mc, Frame, mframe, Frame_I)
  >>> BodyFork = RigidBody('BodyFork', Fork_mc, Fork, mfork, Fork_I)
  >>> BodyWR = RigidBody('BodyWR', WR_mc, WR, mwr, WR_I)
  >>> BodyWF = RigidBody('BodyWF', WF_mc, WF, mwf, WF_I)

  >>> print('Before Forming the List of Nonholonomic Constraints.')
  Before Forming the List of Nonholonomic Constraints.

The kinematic differential equations; they are defined quite simply. Each entry
in this list is equal to zero. ::

  >>> kd = [q1d - u1, q2d - u2, q4d - u4, q5d - u5]

The nonholonomic constraints are the velocity of the front wheel contact point
dotted into the X, Y, and Z directions; the yaw frame is used as it is "closer"
to the front wheel (1 less DCM connecting them). These constraints force the
velocity of the front wheel contact point to be 0 in the inertial frame; the X
and Y direction constraints enforce a "no-slip" condition, and the Z direction
constraint forces the front wheel contact point to not move away from the
ground frame, essentially replicating the holonomic constraint which does not
allow the frame pitch to change in an invalid fashion. ::

  >>> conlist_speed = [WF_cont.vel(N) & Y.x,
  ...                  WF_cont.vel(N) & Y.y,
  ...                  WF_cont.vel(N) & Y.z]

The holonomic constraint is that the position from the rear wheel contact point
to the front wheel contact point when dotted into the normal-to-ground plane
direction must be zero; effectively that the front and rear wheel contact
points are always touching the ground plane. This is actually not part of the
dynamic equations, but instead is necessary for the linearization process. ::

  >>> conlist_coord = [WF_cont.pos_from(WR_cont) & Y.z]

The force list; each body has the appropriate gravitational force applied at
its center of mass. ::

  >>> FL = [(Frame_mc, -mframe * g * Y.z), (Fork_mc, -mfork * g * Y.z),
  ...       (WF_mc, -mwf * g * Y.z), (WR_mc, -mwr * g * Y.z)]
  >>> BL = [BodyFrame, BodyFork, BodyWR, BodyWF]

The N frame is the inertial frame, coordinates are supplied in the order of
independent, dependent coordinates. The kinematic differential equations are
also entered here. Here the independent speeds are specified, followed by the
dependent speeds, along with the non-holonomic constraints. The dependent
coordinate is also provided, with the holonomic constraint. Again, this is only
comes into play in the linearization process, but is necessary for the
linearization to correctly work. ::

  >>> KM = KanesMethod(N, q_ind=[q1, q2, q3],
  ...           q_dependent=[q4], configuration_constraints=conlist_coord,
  ...           u_ind=[u2, u3, u5],
  ...           u_dependent=[u1, u4, u6], velocity_constraints=conlist_speed,
  ...           kd_eqs=kd)
  >>> print('Before Forming Generalized Active and Inertia Forces, Fr and Fr*')
  Before Forming Generalized Active and Inertia Forces, Fr and Fr*
  >>> (fr, frstar) = KM.kanes_equations(FL, BL)
  >>> print('Base Equations of Motion Computed')
  Base Equations of Motion Computed

This is the start of entering in the numerical values from the benchmark paper
to validate the eigenvalues of the linearized equations from this model to the
reference eigenvalues. Look at the aforementioned paper for more information.
Some of these are intermediate values, used to transform values from the paper
into the coordinate systems used in this model. ::

  >>> PaperRadRear  =  0.3
  >>> PaperRadFront =  0.35
  >>> HTA           =  evalf.N(pi/2-pi/10)
  >>> TrailPaper    =  0.08
  >>> rake          =  evalf.N(-(TrailPaper*sin(HTA)-(PaperRadFront*cos(HTA))))
  >>> PaperWb       =  1.02
  >>> PaperFrameCgX =  0.3
  >>> PaperFrameCgZ =  0.9
  >>> PaperForkCgX  =  0.9
  >>> PaperForkCgZ  =  0.7
  >>> FrameLength   =  evalf.N(PaperWb*sin(HTA) - (rake - \
  ...                         (PaperRadFront - PaperRadRear)*cos(HTA)))
  >>> FrameCGNorm   =  evalf.N((PaperFrameCgZ - PaperRadRear - \
  ...                          (PaperFrameCgX/sin(HTA))*cos(HTA))*sin(HTA))
  >>> FrameCGPar    =  evalf.N((PaperFrameCgX / sin(HTA) + \
  ...                          (PaperFrameCgZ - PaperRadRear - \
  ...                           PaperFrameCgX / sin(HTA) * cos(HTA)) * cos(HTA)))
  >>> tempa         =  evalf.N((PaperForkCgZ - PaperRadFront))
  >>> tempb         =  evalf.N((PaperWb-PaperForkCgX))
  >>> tempc         =  evalf.N(sqrt(tempa**2 + tempb**2))
  >>> PaperForkL    =  evalf.N((PaperWb*cos(HTA) - \
  ...                          (PaperRadFront - PaperRadRear)*sin(HTA)))
  >>> ForkCGNorm    =  evalf.N(rake + (tempc * sin(pi/2 - \
  ...                          HTA - acos(tempa/tempc))))
  >>> ForkCGPar     =  evalf.N(tempc * cos((pi/2 - HTA) - \
  ...                          acos(tempa/tempc)) - PaperForkL)

Here is the final assembly of the numerical values. The symbol 'v' is the
forward speed of the bicycle (a concept which only makes sense in the upright,
static equilibrium case?). These are in a dictionary which will later be
substituted in. Again the sign on the *product* of inertia values is flipped
here, due to different orientations of coordinate systems. ::

  >>> v = Symbol('v')
  >>> val_dict = {
  ...       WFrad: PaperRadFront,
  ...       WRrad: PaperRadRear,
  ...       htangle: HTA,
  ...       forkoffset: rake,
  ...       forklength: PaperForkL,
  ...       framelength: FrameLength,
  ...       forkcg1: ForkCGPar,
  ...       forkcg3: ForkCGNorm,
  ...       framecg1: FrameCGNorm,
  ...       framecg3: FrameCGPar,
  ...       Iwr11: 0.0603,
  ...       Iwr22: 0.12,
  ...       Iwf11: 0.1405,
  ...       Iwf22: 0.28,
  ...       Ifork11: 0.05892,
  ...       Ifork22: 0.06,
  ...       Ifork33: 0.00708,
  ...       Ifork31: 0.00756,
  ...       Iframe11: 9.2,
  ...       Iframe22: 11,
  ...       Iframe33: 2.8,
  ...       Iframe31: -2.4,
  ...       mfork: 4,
  ...       mframe: 85,
  ...       mwf: 3,
  ...       mwr: 2,
  ...       g: 9.81,
  ...       q1: 0,
  ...       q2: 0,
  ...       q4: 0,
  ...       q5: 0,
  ...       u1: 0,
  ...       u2: 0,
  ...       u3: v/PaperRadRear,
  ...       u4: 0,
  ...       u5: 0,
  ...       u6: v/PaperRadFront}
  >>> kdd = KM.kindiffdict()
  >>> print('Before Linearization of the \"Forcing\" Term')
  Before Linearization of the "Forcing" Term

Linearizes the forcing vector; the equations are set up as MM udot = forcing,
where MM is the mass matrix, udot is the vector representing the time
derivatives of the generalized speeds, and forcing is a vector which contains
both external forcing terms and internal forcing terms, such as centripetal or
Coriolis forces.  This actually returns a matrix with as many rows as *total*
coordinates and speeds, but only as many columns as independent coordinates and
speeds. (Note that below this is commented out, as it takes a few minutes to
run, which is not good when performing the doctests) ::

  >>> # forcing_lin = KM.linearize()[0].subs(sub_dict)

As mentioned above, the size of the linearized forcing terms is expanded to
include both q's and u's, so the mass matrix must have this done as well.  This
will likely be changed to be part of the linearized process, for future
reference. ::

  >>> MM_full = (KM._k_kqdot).row_join(zeros(4, 6)).col_join(
  ...           (zeros(6, 4)).row_join(KM.mass_matrix))
  >>> print('Before Substitution of Numerical Values')
  Before Substitution of Numerical Values

I think this is pretty self explanatory. It takes a really long time though.
I've experimented with using evalf with substitution, this failed due to
maximum recursion depth being exceeded; I also tried lambdifying this, and it
is also not successful. (again commented out due to speed) ::

  >>> # MM_full = MM_full.subs(val_dict)
  >>> # forcing_lin = forcing_lin.subs(val_dict)
  >>> # print('Before .evalf() call')

  >>> # MM_full = MM_full.evalf()
  >>> # forcing_lin = forcing_lin.evalf()

Finally, we construct an "A" matrix for the form xdot = A x (x being the state
vector, although in this case, the sizes are a little off). The following line
extracts only the minimum entries required for eigenvalue analysis, which
correspond to rows and columns for lean, steer, lean rate, and steer rate.
(this is all commented out due to being dependent on the above code, which is
also commented out)::

  >>> # Amat = MM_full.inv() * forcing_lin
  >>> # A = Amat.extract([1,2,4,6],[1,2,3,5])
  >>> # print(A)
  >>> # print('v = 1')
  >>> # print(A.subs(v, 1).eigenvals())
  >>> # print('v = 2')
  >>> # print(A.subs(v, 2).eigenvals())
  >>> # print('v = 3')
  >>> # print(A.subs(v, 3).eigenvals())
  >>> # print('v = 4')
  >>> # print(A.subs(v, 4).eigenvals())
  >>> # print('v = 5')
  >>> # print(A.subs(v, 5).eigenvals())

Upon running the above code yourself, enabling the commented out lines, compare
the computed eigenvalues to those is the referenced paper. This concludes the
bicycle example.

The Rolling Disc Example using Lagranges Method
===============================================

Here the rolling disc is formed from the contact point up, removing the
need to introduce generalized speeds. Only 3 configuration and 3
speed variables are needed to describe this system, along with the
disc's mass and radius, and the local gravity. ::

  >>> from sympy import symbols, cos, sin
  >>> from sympy.physics.mechanics import *
  >>> mechanics_printing()
  >>> q1, q2, q3 = dynamicsymbols('q1 q2 q3')
  >>> q1d, q2d, q3d = dynamicsymbols('q1 q2 q3', 1)
  >>> r, m, g = symbols('r m g')

The kinematics are formed by a series of simple rotations. Each simple
rotation creates a new frame, and the next rotation is defined by the new
frame's basis vectors. This example uses a 3-1-2 series of rotations, or
Z, X, Y series of rotations. Angular velocity for this is defined using
the second frame's basis (the lean frame). ::

  >>> N = ReferenceFrame('N')
  >>> Y = N.orientnew('Y', 'Axis', [q1, N.z])
  >>> L = Y.orientnew('L', 'Axis', [q2, Y.x])
  >>> R = L.orientnew('R', 'Axis', [q3, L.y])

This is the translational kinematics. We create a point with no velocity
in N; this is the contact point between the disc and ground. Next we form
the position vector from the contact point to the disc's center of mass.
Finally we form the velocity and acceleration of the disc. ::

  >>> C = Point('C')
  >>> C.set_vel(N, 0)
  >>> Dmc = C.locatenew('Dmc', r * L.z)
  >>> Dmc.v2pt_theory(C, N, R)
  r*(sin(q2)*q1' + q3')*L.x - r*q2'*L.y

Forming the inertia dyadic. ::

  >>> I = inertia(L, m / 4 * r**2, m / 2 * r**2, m / 4 * r**2)
  >>> mprint(I)
  m*r**2/4*(L.x|L.x) + m*r**2/2*(L.y|L.y) + m*r**2/4*(L.z|L.z)
  >>> BodyD = RigidBody('BodyD', Dmc, R, m, (I, Dmc))

We then set the potential energy and determine the Lagrangian of the rolling
disc. ::

  >>> BodyD.set_potential_energy(- m * g * r * cos(q2))
  >>> Lag = Lagrangian(N, BodyD)

Then the equations of motion are generated by initializing the
``LagrangesMethod`` object. Finally we solve for the generalized
accelerations(q double dots) with the ``rhs`` method. ::

  >>> q = [q1, q2, q3]
  >>> l = LagrangesMethod(Lag, q)
  >>> le = l.form_lagranges_equations()
  >>> le.simplify(); le
  [m*r**2*(12*sin(q2)*q3'' + 10*sin(2*q2)*q1'*q2' + 12*cos(q2)*q2'*q3' - 5*cos(2*q2)*q1'' + 7*q1'')/8]
  [                     m*r*(8*g*sin(q2) - 5*r*sin(2*q2)*q1'**2 - 12*r*cos(q2)*q1'*q3' + 10*r*q2'')/8]
  [                                                3*m*r**2*(sin(q2)*q1'' + cos(q2)*q1'*q2' + q3'')/2]
  >>> lrhs = l.rhs(); lrhs.simplify(); lrhs
  [                                                                q1']
  [                                                                q2']
  [                                                                q3']
  [                     -2*(2*cos(q2)*tan(q2)*q1' + 3*q3')*q2'/cos(q2)]
  [(-8*g*sin(q2) + 5*r*sin(2*q2)*q1'**2 + 12*r*cos(q2)*q1'*q3')/(10*r)]
  [               (-5*cos(q2)*q1' + 6*tan(q2)*q3' + 4*q1'/cos(q2))*q2']