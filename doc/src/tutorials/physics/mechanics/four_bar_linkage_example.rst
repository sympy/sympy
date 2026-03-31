==================
A four bar linkage
==================

The four bar linkage is a common example used in mechanics, which can be
formulated with only two holonomic constraints. This example will make use of
joints functionality provided in :mod:`sympy.physics.mechanics`. In summary we
will use bodies and joints to define the open loop system. Next, we define the
configuration constraints to close the loop. :class:`~.System` will be used to
do the "book-keeping" of the entire system with :class:`~.KanesMethod` as the
backend.

.. raw:: html
   :file: four_bar_linkage.svg

First we need to create the :func:`~.dynamicsymbols` needed to describe the
system as shown in the above diagram. In this case, the generalized coordinates
`q_1`, `q_2` and `q_3` represent the angles between the links. Likewise, the
generalized speeds `u_1`, `u_2` and `u_3` represent the angular velocities
between the links. We also create some :func:`~.symbols` to represent the
lengths and density of the links. ::

   >>> from sympy import Matrix, linear_eq_to_matrix, pi, simplify, symbols
   >>> from sympy.physics.mechanics import *
   >>> mechanics_printing(pretty_print=False)
   >>> q1, q2, q3, u1, u2, u3 = dynamicsymbols('q1:4, u1:4')
   >>> l1, l2, l3, l4, rho = symbols('l1:5, rho')

With all symbols defined, we can now define the bodies and initialize our
instance of :class:`~.System`. ::

   >>> N = ReferenceFrame('N')
   >>> mass_centers = [Point(f'mc{i}') for i in range(1, 5)]
   >>> inertias = [Inertia.from_inertia_scalars(P, N, 0, 0, rho*l**3/12)
   ...             for P, l in zip(mass_centers, (l1, l2, l3, l4))]
   >>> link1 = RigidBody('Link1', frame=N, mass=rho*l1,
   ...                   masscenter=mass_centers[0], inertia=inertias[0])
   >>> link2 = RigidBody('Link2', mass=rho*l2, masscenter=mass_centers[1],
   ...                   inertia=inertias[1])
   >>> link3 = RigidBody('Link3', mass=rho*l3, masscenter=mass_centers[2],
   ...                   inertia=inertias[2])
   >>> link4 = RigidBody('Link4', mass=rho*l4, masscenter=mass_centers[3],
   ...                   inertia=inertias[3])
   >>> system = System.from_newtonian(link1)

Next, we also define the first three joints, which create the open loop pendulum, and
add them to the system. ::

   >>> joint1 = PinJoint('J1', link1, link2, coordinates=q1, speeds=u1,
   ...                   parent_point=l1/2*link1.x,
   ...                   child_point=-l2/2*link2.x, joint_axis=link1.z)
   >>> joint2 = PinJoint('J2', link2, link3, coordinates=q2, speeds=u2,
   ...                   parent_point=l2/2*link2.x,
   ...                   child_point=-l3/2*link3.x, joint_axis=link2.z)
   >>> joint3 = PinJoint('J3', link3, link4, coordinates=q3, speeds=u3,
   ...                   parent_point=l3/2*link3.x,
   ...                   child_point=-l4/2*link4.x, joint_axis=link3.z)
   >>> system.add_joints(joint1, joint2, joint3)

Now we can formulate the holonomic constraint that will close the kinematic
loop. ::

   >>> start_point = link1.masscenter.locatenew('start_point', -l1/2*link1.x)
   >>> end_point = link4.masscenter.locatenew('end_point', l4/2*link4.x)
   >>> loop = end_point.pos_from(start_point)
   >>> system.add_holonomic_constraints(loop.dot(link1.x), loop.dot(link1.y))

Before generating the equations of motion we need to specify which generalized
coordinates and speeds are independent and which are dependent. After which we
can run :meth:`~.System.validate_system` to do some basic consistency checks. ::

   >>> system.q_ind = [q1]
   >>> system.u_ind = [u1]
   >>> system.q_dep = [q2, q3]
   >>> system.u_dep = [u2, u3]
   >>> system.validate_system()

As we have the entire system ready, we can now form the equations of motion
using :class:`~.KanesMethod` as the backend. ::

   >>> simplify(system.form_eoms())
    Matrix([[l2*rho*(-2*l2**2*sin(q3)*u1' + 3*l2*l3*u1**2*sin(q2 + q3)*sin(q2) + 3*l2*l3*sin(q2)*cos(q2 + q3)*u1' - 3*l2*l3*sin(q3)*u1' + 3*l2*l4*u1**2*sin(q2 + q3)*sin(q2) + 3*l2*l4*sin(q2)*cos(q2 + q3)*u1' + 3*l3**2*u1**2*sin(q2)*sin(q3) + 6*l3**2*u1*u2*sin(q2)*sin(q3) + 3*l3**2*u2**2*sin(q2)*sin(q3) + 2*l3**2*sin(q2)*cos(q3)*u1' + 2*l3**2*sin(q2)*cos(q3)*u2' - l3**2*sin(q3)*cos(q2)*u1' - l3**2*sin(q3)*cos(q2)*u2' + 3*l3*l4*u1**2*sin(q2)*sin(q3) + 6*l3*l4*u1*u2*sin(q2)*sin(q3) + 3*l3*l4*u2**2*sin(q2)*sin(q3) + 3*l3*l4*sin(q2)*cos(q3)*u1' + 3*l3*l4*sin(q2)*cos(q3)*u2' + l4**2*sin(q2)*u1' + l4**2*sin(q2)*u2' + l4**2*sin(q2)*u3')/(6*sin(q3))]])

Revealing noncontributing forces
--------------------------------

To reveal the noncontributing forces at the closing joint, we must introduce
auxiliary speeds in the x and y-direction at the endpoint.

   >>> uaux1, uaux2 = dynamicsymbols('uaux1:3')
   >>> end_point_aux = end_point.locatenew('end_point_aux', 0)
   >>> end_point_aux.set_vel(N, end_point.vel(N) + uaux1*N.x + uaux2*N.y)

To ensure that speeds are included in the velocity constraints, we must manually
overwrite the velocity constraints because those are by default specified as the
time derivatives of the holonomic constraints.

   >>> system.velocity_constraints = [
   ...    end_point_aux.vel(N).dot(N.x), end_point_aux.vel(N).dot(N.y)]

When adding the noncontributing forces we need them to depend only on the
auxiliary velocity and not the velocity that is eliminated by the constraints.
This can be achieved by applying an equal and opposite force to the
non-auxiliary endpoint.

   >>> faux1, faux2 = dynamicsymbols('faux1:3')
   >>> noncontributing_forces = [
   ...   Force(end_point_aux, faux1*N.x + faux2*N.y),
   ...   Force(end_point, -(faux1*N.x + faux2*N.y)),
   ... ]

Alternatively, we can specify a new point that already subtracts the velocity
eliminated by the constraints.

   >>> end_point_forces = end_point.locatenew('end_point_forces', 0)
   >>> end_point_forces.set_vel(N, uaux1*N.x + uaux2*N.y)
   >>> noncontributing_forces = [Force(end_point_forces, faux1*N.x + faux2*N.y)]

Next, we can add the auxiliary speeds and noncontributing forces to the system.

   >>> system.add_loads(*noncontributing_forces)
   >>> system.u_aux = [uaux1, uaux2]

To include gravity we can use :meth:`~.System.apply_uniform_gravity` before
validating the system and forming the equations of motion.

   >>> g = symbols('g')
   >>> system.apply_uniform_gravity(-g*N.y)
   >>> system.validate_system()
   >>> eoms = system.form_eoms()

With the equations of motion formed we can solve the auxiliary equations for the
noncontributing forces and compute their values for a simple configuration.

   >>> auxiliary_eqs = system.eom_method.auxiliary_eqs
   >>> forces_eqs = Matrix.LUsolve(
   ...   *linear_eq_to_matrix(auxiliary_eqs, [faux1, faux2]))
   >>> subs = {
   ...   l1: 2, l2: 1, l3: 2, l4: 1,
   ...   rho: 5, g: 9.81,
   ...   q1: pi/2, q2: pi/2, q3: pi/2,
   ...   u1: 0, u2: 0, u3: 0, u1.diff(): 0, u2.diff(): 0, u3.diff(): 0,
   ... }
   >>> forces_eqs.xreplace(subs)
   Matrix([
   [    0],
   [-98.1]])
