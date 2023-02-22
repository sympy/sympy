=========================================
Multi Degree of Freedom Holonomic System
=========================================

In this example we demonstrate the use of the functionality provided in
:mod:`sympy.physics.mechanics` for deriving the equations of motion (EOMs) of a
holonomic system that includes both particles and rigid bodies with contributing
forces and torques, some of which are specified forces and torques. The system
is shown below:

.. image:: multidof-holonomic.*
   :align: center

The system will be modeled using ``System``. First we need to create the
``dynamicsymbols`` needed to describe the system as shown in the above diagram.
In this case, the generalized coordinate :math:`q_1` represents the lateral
distance of the block from the wall, :math:`q_2` represents the angle of the
compound pendulum from the vertical, and :math:`q_3` represents the angle of
the simple pendulum relative to the compound pendulum. The generalized speed
:math:`u_1` represents the lateral speed of the block, :math:`u_2` represents
the angular velocity of the compound pendulum, and :math:`u_3` represents the
angular velocity of the simple pendulum relative to the compound pendulum.

We also create some ``symbols`` to represent the length and mass of the
pendulum, as well as gravity and others. ::

   >>> from sympy import zeros, symbols
   >>> from sympy.physics.mechanics import *
   >>> q1, q2, q3, u1, u2, u3 = dynamicsymbols('q1, q2, q3, u1, u2, u3')
   >>> l, k, c, g, kT = symbols('l, k, c, g, kT')
   >>> ma, mb, mc, IBzz = symbols('ma, mb, mc, IBzz')

Next, we create the bodies and a system to do the book-keeping.

   >>> wall = RigidBody('N')
   >>> system = System.from_newtonian(wall)
   >>> block = Particle('A', mass=ma)
   >>> IB = Inertia.from_inertia_scalars(Point('B_masscenter'),
   ...                                   wall.frame, 0, 0, IBzz)
   >>> compound_pend = RigidBody('B', IB.point, mass=mb, inertia=IB)
   >>> simple_pend = Particle('C', mass=mc)

Now we can connect the bodies using joints to establish the kinematics.
Particles do not have frames themselves. Therefore we can either create a frame
beforehand or reuse a frame from the previous joint. After creating the joints
we can add them to the system. ::

   >>> slider = PrismaticJoint('J1', wall, block, coordinates=q1, speeds=u1)
   >>> rev1 = PinJoint('J2', block, compound_pend, coordinates=q2, speeds=u2,
   ...                 parent_interframe=slider.child_interframe,
   ...                 joint_axis=wall.z, child_point=l*2/3*compound_pend.y)
   >>> C_frame = ReferenceFrame('C_frame')
   >>> rev2 = PinJoint('J3', compound_pend, simple_pend, coordinates=q3, speeds=u3,
   ...                 joint_axis=compound_pend.z, parent_point=-l/3*compound_pend.y,
   ...                 child_interframe=C_frame, child_point=l*C_frame.y)
   >>> system.add_joints(slider, rev1, rev2)

Now we can apply loads (forces and torques) to the bodies: gravity acts on all
bodies, a linear spring and damper act on the block and the wall, a rotational
linear spring acts between the compound_pend and the simple_pend, a specified
torque T acts on the compound_pend and the block, and a specified force F acts
on the block. ::

    >>> F, T = dynamicsymbols('F, T')
    >>> system.apply_gravity(-g * wall.y)
    >>> system.apply_force(block, F * wall.x)
    >>> system.apply_force(block, -k*q1*wall.x, wall)
    >>> system.apply_force(block, -c*u1*wall.x, wall)
    >>> system.apply_torque(compound_pend, T * compound_pend.z, slider.child_interframe)
    >>> system.apply_torque(C_frame, -kT*q3*compound_pend.z, compound_pend)

With the problem set up, the equations of motion can be generated using the
``System`` class with :class:`~.KanesMethod` in the backend. ::

    >>> system.form_eoms()
    Matrix([
    [                                -c*u1(t) - k*q1(t) + 2*l*mb*u2(t)**2*sin(q2(t))/3 - l*mc*(-sin(q2(t))*sin(q3(t)) + cos(q2(t))*cos(q3(t)))*Derivative(u3(t), t) - l*mc*(-sin(q2(t))*cos(q3(t)) - sin(q3(t))*cos(q2(t)))*(u2(t) + u3(t))**2 + l*mc*u2(t)**2*sin(q2(t)) - (2*l*mb*cos(q2(t))/3 + mc*(l*(-sin(q2(t))*sin(q3(t)) + cos(q2(t))*cos(q3(t))) + l*cos(q2(t))))*Derivative(u2(t), t) - (ma + mb + mc)*Derivative(u1(t), t) + F(t)],
    [-2*g*l*mb*sin(q2(t))/3 - g*l*mc*(sin(q2(t))*cos(q3(t)) + sin(q3(t))*cos(q2(t))) - g*l*mc*sin(q2(t)) + l**2*mc*(u2(t) + u3(t))**2*sin(q3(t)) - l**2*mc*u2(t)**2*sin(q3(t)) - mc*(l**2*cos(q3(t)) + l**2)*Derivative(u3(t), t) - (2*l*mb*cos(q2(t))/3 + mc*(l*(-sin(q2(t))*sin(q3(t)) + cos(q2(t))*cos(q3(t))) + l*cos(q2(t))))*Derivative(u1(t), t) - (IBzz + 4*l**2*mb/9 + mc*(2*l**2*cos(q3(t)) + 2*l**2))*Derivative(u2(t), t) + T(t)],
    [                                                                                                                                                                        -g*l*mc*(sin(q2(t))*cos(q3(t)) + sin(q3(t))*cos(q2(t))) - kT*q3(t) - l**2*mc*u2(t)**2*sin(q3(t)) - l**2*mc*Derivative(u3(t), t) - l*mc*(-sin(q2(t))*sin(q3(t)) + cos(q2(t))*cos(q3(t)))*Derivative(u1(t), t) - mc*(l**2*cos(q3(t)) + l**2)*Derivative(u2(t), t)]])

    >>> system.mass_matrix_full
    Matrix([
    [-1,  0,  0,                                                                                            0,                                                                                            0,                         0],
    [ 0, -1,  0,                                                                                            0,                                                                                            0,                         0],
    [ 0,  0, -1,                                                                                            0,                                                                                            0,                         0],
    [ 0,  0,  0,                                                                                 ma + mb + mc, 2*l*mb*cos(q2(t))/3 + mc*(l*(-sin(q2(t))*sin(q3(t)) + cos(q2(t))*cos(q3(t))) + l*cos(q2(t))), l*mc*(-sin(q2(t))*sin(q3(t)) + cos(q2(t))*cos(q3(t)))],
    [ 0,  0,  0, 2*l*mb*cos(q2(t))/3 + mc*(l*(-sin(q2(t))*sin(q3(t)) + cos(q2(t))*cos(q3(t))) + l*cos(q2(t))),                                         IBzz + 4*l**2*mb/9 + mc*(2*l**2*cos(q3(t)) + 2*l**2),                           mc*(l**2*cos(q3(t)) + l**2)],
    [ 0,  0,  0,                                        l*mc*(-sin(q2(t))*sin(q3(t)) + cos(q2(t))*cos(q3(t))),                                                                  mc*(l**2*cos(q3(t)) + l**2),                   l**2*mc]])


    >>> system.forcing_full
    Matrix([
    [                                                                                                                                                                          -u1(t)],
    [                                                                                                                                                                          -u2(t)],
    [                                                                                                                                                                          -u3(t)],
    [                  -c*u1(t) - k*q1(t) + 2*l*mb*u2(t)**2*sin(q2(t))/3 - l*mc*(-sin(q2(t))*cos(q3(t)) - sin(q3(t))*cos(q2(t)))*(u2(t) + u3(t))**2 + l*mc*u2(t)**2*sin(q2(t)) + F(t)],
    [-2*g*l*mb*sin(q2(t))/3 - g*l*mc*(sin(q2(t))*cos(q3(t)) + sin(q3(t))*cos(q2(t))) - g*l*mc*sin(q2(t)) + l**2*mc*(u2(t) + u3(t))**2*sin(q3(t)) - l**2*mc*u2(t)**2*sin(q3(t)) + T(t)],
    [                                                                                -g*l*mc*(sin(q2(t))*cos(q3(t)) + sin(q3(t))*cos(q2(t))) - kT*q3(t) - l**2*mc*u2(t)**2*sin(q3(t))]])

