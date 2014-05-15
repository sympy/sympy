==================================
A rolling disc, with Kane's method
==================================

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
  >>> mechanics_printing(pretty_print=False)

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
  >>> w_R_N_qd = R.ang_vel_in(N)
  >>> R.set_ang_vel(N, u1 * L.x + u2 * L.y + u3 * L.z)

This is the translational kinematics. We create a point with no velocity
in N; this is the contact point between the disc and ground. Next we form
the position vector from the contact point to the disc's center of mass.
Finally we form the velocity and acceleration of the disc. ::

  >>> C = Point('C')
  >>> C.set_vel(N, 0)
  >>> Dmc = C.locatenew('Dmc', r * L.z)
  >>> Dmc.v2pt_theory(C, N, R)
  r*u2*L.x - r*u1*L.y

This is a simple way to form the inertia dyadic. The inertia of the disc does
not change within the lean frame as the disc rolls; this will make for simpler
equations in the end. ::

  >>> I = inertia(L, m / 4 * r**2, m / 2 * r**2, m / 4 * r**2)
  >>> mprint(I)
  m*r**2/4*(L.x|L.x) + m*r**2/2*(L.y|L.y) + m*r**2/4*(L.z|L.z)

Kinematic differential equations; how the generalized coordinate time
derivatives relate to generalized speeds. ::

  >>> kd = [dot(R.ang_vel_in(N) - w_R_N_qd, uv) for uv in L]

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
  Matrix([
  [(4*g*sin(q2) + 6*r*u2*u3 - r*u3**2*tan(q2))/(5*r)],
  [                                       -2*u1*u3/3],
  [                          (-2*u2 + u3*tan(q2))*u1]])
