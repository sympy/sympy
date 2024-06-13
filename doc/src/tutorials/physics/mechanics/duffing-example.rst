.. _duffing-example:

==================================
Duffing Oscillator with a Pendulum
==================================

In this example we demonstrate the use of functionality provided in
:obj:`sympy.physics.mechanics` for deriving the quations of motion (EOM) for a system
consisting of a Duffing oscillator with a pendulum.

.. _fig-duffing-oscillator-pendulum:
.. figure:: duffing-oscillator-pendulum.svg

The system will be modeled using Lagrange equations. `M` is mass of the Duffing oscillator,
`m` is mass of the pendulum, `l` is length of the pendulum. `k_1` and `k_2` are linear and
non-linear parts of spring stiffness, and `c_1` is a viscous damping coefficient of the Duffing oscillator.

   >>> import sympy as sm
   >>> import sympy.physics.mechanics as me

Define Variables
================

   >>> M, m, l, k1, k2, c1, g, h, w, d, r = sm.symbols('M, m, l, k1, k2, c1, g, h, w, d, r')
   >>> q1, q2 = me.dynamicsymbols('q1 q2')
   >>> q1d, q2d = me.dynamicsymbols('q1 q2', 1)
   >>> u1, u2 = me.dynamicsymbols('u1 u2')
   >>> u1d, u2d = me.dynamicsymbols('u1 u2', 1)

- :math:`h`: Height of the Duffing oscillator
- :math:`w`: Width of the Duffing oscillator
- :math:`d`: Depth of the Duffing oscillator
- :math:`r`: Radius of the massive bob of the pendulum
- :math:`q_1`: Generalized coordinate representing the position of the Duffing oscillator
- :math:`q_2`: Generalized coordinate representing the angle of the pendulum
- :math:`\dot{q}_1`: First time derivative of `q1`, representing the velocity of the Duffing oscillator
- :math:`\dot{q}_2`: First time derivative of `q2`, representing the angular velocity of the pendulum
- :math:`u_1`: Generalized speed associated with the Duffing oscillator
- :math:`u_2`: Generalized speed associated with the pendulum
- :math:`\dot{u}_1`: First time derivative of `u1`, representing the acceleration of the Duffing oscillator
- :math:`\dot{u}_2`: First time derivative of `u2`, representing the angular acceleration of the pendulum

Define Kinematics
=================

Define all the reference frames and points.

   >>> N = me.ReferenceFrame('N')
   >>> B = N.orientnew('B', 'axis', (q2, N.z))

The angular velocity of the pendulum in the reference frame is:

   >>> B.set_ang_vel(N, q2.diff() * N.z)
   >>> B.ang_vel_in(N)
   Derivative(q2(t), t)*N.z

Locations and velocities of the Duffing Oscillator block and the pendulum are:

   >>> O = me.Point('O')
   >>> Block = O.locatenew('Block', q1 * N.y)
   >>> Pendulum = Block.locatenew('Pendulum', l * B.y)

O is a fixed point in the inertial reference frame.

   >>> O.set_vel(N, 0)
   >>> Block.set_vel(N, q1d * N.y)
   >>> Pendulum.v2pt_theory(Block, N, B)
   Derivative(q1(t), t)*N.y - l*Derivative(q2(t), t)*B.x

Define inertia and rigid bodies.
Here, we assume a simple pendulum which consists of a bob of mass m hanging from a massless string of length l
and fixed at a pivot point (Duffing Oscillator Block).

   >>> I_block = me.inertia(N, M*(h**2 + d**2)/12, M*(w**2 + h**2)/12, M*(w**2 + d**2)/12)
   >>> I_pendulum = me.inertia(B, 0, 0, 2*m*r**2/5 + m*l**2/3)

   >>> par_block = me.RigidBody('block', Block, N, M, (I_block, Block))
   >>> par_pendulum = me.RigidBody('pendulum', Pendulum, B, m, (I_pendulum, Pendulum))

Define Force, Energy
====================

We obtain the Duffing force using the `DuffingSpring` actuator from ``sympy.physics.mechanics.actuator``.
This force will be used to calculate the potential energy of the Duffing Oscillator block.

Define the Duffing spring force.

   >>> pathway = me.LinearPathway(O, Block)
   >>> duffing_spring = me.DuffingSpring(k1, k2, pathway, 0)
   >>> duffing_force = -duffing_spring.force
   >>> duffing_force
   k1*sqrt(q1(t)**2) + k2*(q1(t)**2)**(3/2)

Define Rayleigh dissipation.

   >>> D = (1/2) * c1 * q1d**2

In relation to Lagrange's Method, we derive both the kinetic and potential energies of the system.

Kinetic Energy

   >>> Kinetic = par_block.kinetic_energy(N) + par_pendulum.kinetic_energy(N)
   >>> Kinetic
   M*Derivative(q1(t), t)**2/2 + m*r**2*Derivative(q2(t), t)**2/5 + m*(l**2*Derivative(q2(t), t)**2 - 2*l*sin(q2(t))*Derivative(q1(t), t)*Derivative(q2(t), t) + Derivative(q1(t), t)**2)/2

Potential Energy

   >>> par_block.potential_energy = sm.integrate(duffing_force, q1)
   >>> par_pendulum.potential_energy = m * g * (l - l/(sm.sqrt(1+q2**2/q1**2)))
   >>> me.potential_energy(par_block, par_pendulum)
   g*m*(l - l/sqrt(1 + q2(t)**2/q1(t)**2)) + k1*sqrt(q1(t)**2)*q1(t)/2 + k2*(q1(t)**2)**(3/2)*q1(t)/4

Lagrange's Method
=================

With the problem setup, the Lagrangian can be calculated, and the equations of motion formed.
In the force list FL, we specify forces in the format (point, the force acting on the particle).

   >>> L = me.Lagrangian(N, par_block, par_pendulum)
   >>> me.Lagrangian(N, par_block, par_pendulum)
   M*Derivative(q1(t), t)**2/2 - g*m*(l - l/sqrt(1 + q2(t)**2/q1(t)**2)) - k1*sqrt(q1(t)**2)*q1(t)/2 - k2*(q1(t)**2)**(3/2)*q1(t)/4 + m*r**2*Derivative(q2(t), t)**2/5 + m*(l**2*Derivative(q2(t), t)**2 - 2*l*sin(q2(t))*Derivative(q1(t), t)*Derivative(q2(t), t) + Derivative(q1(t), t)**2)/2

   >>> FL = [(Block, duffing_force * N.y + D * N.y), (Pendulum, - m * g * N.y)]
   >>> LM = me.LagrangesMethod(L, [q1, q2], forcelist = FL, frame = N)
   >>> LM.form_lagranges_equations()
   Matrix([
    [                                                                    M*Derivative(q1(t), (t, 2)) - 0.5*c1*Derivative(q1(t), t)**2 - g*l*m*q2(t)**2/((1 + q2(t)**2/q1(t)**2)**(3/2)*q1(t)**3) + g*m + m*(-2*l*sin(q2(t))*Derivative(q2(t), (t, 2)) - 2*l*cos(q2(t))*Derivative(q2(t), t)**2 + 2*Derivative(q1(t), (t, 2)))/2],
    [-g*l*m*sin(q2(t)) + g*l*m*q2(t)/((1 + q2(t)**2/q1(t)**2)**(3/2)*q1(t)**2) + l*m*cos(q2(t))*Derivative(q1(t), t)*Derivative(q2(t), t) + 2*m*r**2*Derivative(q2(t), (t, 2))/5 + m*(2*l**2*Derivative(q2(t), (t, 2)) - 2*l*sin(q2(t))*Derivative(q1(t), (t, 2)) - 2*l*cos(q2(t))*Derivative(q1(t), t)*Derivative(q2(t), t))/2]])

References
==========

.. [P.Brzeskia2012] P. Brzeskia, P. Perlikowskia, S. Yanchukb, T. Kapitaniaka,
   The dynamics of the pendulum suspended on the forced Duffing oscillator,
   Journal of Sound and Vibration, 2012, https://doi.org/10.48550/arXiv.1202.5937
