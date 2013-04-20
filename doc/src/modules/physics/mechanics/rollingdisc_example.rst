==============
A rolling disc
==============

The disc is assumed to be infinitely thin, in contact with the ground at only 1
point, and it is rolling without slip on the ground. See the image below.

.. image:: rollingdisc.svg
   :height: 350
   :width: 350
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



A rolling disc, with constraint forces
======================================

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

The rolling disc using Lagrange's Method
========================================

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
