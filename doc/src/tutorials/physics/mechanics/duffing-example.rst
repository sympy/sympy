.. _duffing-example:

==================================
Duffing Oscillator with a Pendulum
==================================

In this example we demonstrate the use of functionality provided in
:obj:`sympy.physics.mechanics` for deriving the equations of motion for a system
consisting of a Duffing oscillator with a pendulum. This example is inspired by the
paper [P.Brzeskia2012]_ section 2.

.. raw:: html
   :file: duffing.svg

The system will be modeled using Lagrange equations. `M` is mass of the Duffing oscillator,
`m` is mass of the pendulum, `l` is length of the pendulum. `k_1` and `k_2` are linear and
non-linear parts of spring stiffness, and `c_1` is a viscous damping coefficient of the Duffing oscillator.

   >>> import sympy as sm
   >>> import sympy.physics.mechanics as me
   >>> me.init_vprinting()

Define Variables
================

   >>> M, m, l, k1, k2, c1, g, h, w, d, r = sm.symbols('M, m, l, k1, k2, c1, g, h, w, d, r')
   >>> q1, q2 = me.dynamicsymbols('q1 q2')
   >>> q1d = me.dynamicsymbols('q1', 1)

- :math:`h`: Height of the Duffing oscillator
- :math:`w`: Width of the Duffing oscillator
- :math:`d`: Depth of the Duffing oscillator
- :math:`r`: Radius of the massive bob of the pendulum
- :math:`q_1`: Generalized coordinate representing the position of the Duffing oscillator
- :math:`q_2`: Generalized coordinate representing the angle of the pendulum

Define Kinematics
=================

Define all the reference frames and points.

   >>> N = me.ReferenceFrame('N')
   >>> B = N.orientnew('B', 'axis', (q2, N.z))

The angular velocity of the pendulum in the reference frame is:

   >>> B.ang_vel_in(N)
   q2'(t) n_z

Locations and velocities of the Duffing Oscillator block and the pendulum are:

   >>> O = me.Point('O')
   >>> block_point = O.locatenew('block', q1 * N.y)
   >>> pendulum_point = block_point.locatenew('pendulum', l * B.y)

O is a fixed point in the inertial reference frame.

   >>> O.set_vel(N, 0)
   >>> block_point.set_vel(N, q1d * N.y)
   >>> pendulum_point.v2pt_theory(block_point, N, B)
   q1'(t) n_y + -l*q2'(t) b_x

Define inertia and rigid bodies.
Here, we assume a simple pendulum which consists of a bob of mass m hanging from a massless string of length l
and fixed at a pivot point (Duffing Oscillator Block).

   >>> I_block = M / 12 * me.inertia(N, h**2 + d**2, w**2 + d**2, w**2 + h**2)
   >>> I_pendulum = 2*m*r**2/5*me.inertia(B, 1, 0, 1)

   >>> block_body = me.RigidBody('block', block_point, N, M, (I_block, block_point))
   >>> pendulum_body = me.RigidBody('pendulum', pendulum_point, B, m, (I_pendulum, pendulum_point))

Define Forces
=============

We calculate the forces acting on the system.
In this example, we set the potential energy to zero in the Lagrangian, and include the
conservative forces (gravity and the Duffing spring) in the loads.

   >>> path = me.LinearPathway(O, block_point)
   >>> spring = me.DuffingSpring(k1, k2, path, 0)
   >>> damper = me.LinearDamper(c1, path)

   >>> loads = spring.to_loads() + damper.to_loads()

   >>> bodies = [block_body, pendulum_body]

   >>> for body in bodies:
   ...     loads.append(me.Force(body, body.mass * g * N.y))

   >>> loads
               /      _____           3/2\                  /        _____           3/2\
         |     /   2       /  2\   |                  |       /   2       /  2\   |
         \k1*\/  q1   + k2*\q1 /   /*q1               \- k1*\/  q1   - k2*\q1 /   /*q1
    [(O, ------------------------------ n_y), (block, -------------------------------- n_y), (O, c1*q1'(t) n_y), (block, -c1*q1'(t) n_y), (block, M*g n_y), (pendulum, g*m n_y)]
                       _____                                         _____
                      /   2                                         /   2
                    \/  q1                                        \/  q1

Lagrange's Method
=================

With the problem setup, the Lagrangian can be calculated, and the equations of motion formed.

   >>> L = me.Lagrangian(N, block_body, pendulum_body)
   >>> L
            2      2       2     / 2       2                                     2\
    M*q1'(t)    m*r *q2'(t)    m*\l *q2'(t)  - 2*l*sin(q2)*q1'(t)*q2'(t) + q1'(t) /
    --------- + ------------ + ----------------------------------------------------
        2            5                                  2

   >>> LM = me.LagrangesMethod(L, [q1, q2], bodies=bodies, forcelist=loads, frame=N)
   >>> sm.simplify(LM.form_lagranges_equations())
    [                                       /                                    2          \   /          2\   ]
    [-M*g + M*q1''(t) + c1*q1'(t) - g*m - m*\l*sin(q2)*q2''(t) + l*cos(q2)*q2'(t)  - q1''(t)/ + \k1 + k2*q1 /*q1]
    [                                                                                                           ]
    [                     /                   2                                    2        \                   ]
    [                   m*\5*g*l*sin(q2) + 5*l *q2''(t) - 5*l*sin(q2)*q1''(t) + 2*r *q2''(t)/                   ]
    [                   ---------------------------------------------------------------------                   ]
    [                                                     5                                                     ]

Equations of motion in [P.Brzeskia2012]_:

.. math::

   (M + m)y'' - ml \phi'' \sin(\phi) - ml(\phi')^2 \cos(\phi) + k_1 y + k_2 y^3 + c_1 y' = F_0 \cos(\nu t)

   ml^2 \phi'' - mly'' \sin(\phi) + mlg \sin(\phi) + c_2 \phi' = 0

Equations of motion in this example:

.. math::

   (M + m)q_1'' - mlq_2'' \sin(q_2) - ml(q_2')^2 \cos(q_2) + k_1 q_1 + k_2 q_1^3 + c_1 q_1' - (M + m)g = 0

   ml^2 q_2'' - mlq_1'' \sin(q_2) + mlg \sin(q_2) + \frac{2r^2q_2''}{5} = 0

The differences in the equations of motion are attributed to several factors: the gravitational force, a damping torque characterized by the damping coefficient `c_2`,
and a periodically varying excitation :math:`F_0 \cos(\nu t)`.

References
==========

.. [P.Brzeskia2012] P. Brzeskia, P. Perlikowskia, S. Yanchukb, T. Kapitaniaka,
   The dynamics of the pendulum suspended on the forced Duffing oscillator,
   Journal of Sound and Vibration, 2012, https://doi.org/10.48550/arXiv.1202.5937
