.. _duffing-example:

==================================
Duffing Oscillator with a Pendulum
==================================

In this example we demonstrate the use of functionality provided in
:obj:`sympy.physics.mechanics` for deriving the quations of motion (EOM) for a system
consisting of a Duffing oscillator with a pendulum.

(schematic diagram)

The system will be modeled using Lagrange equations. `M` is mass of the Duffing oscillator,
`m` is mass of the pendulum, `l` is length of the pendulum. `k_1` and `k_2` are linear and
non-linear parts of spring stiffness, and `c_1` is a viscous damping coefficient of the Duffing oscillator.

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs
   :nofigs:

   >>> import sympy as sm
   >>> import sympy.physics.mechanics as me
   >>> from sympy.physics.mechanics.actuator import DuffingSpring
   >>> from sympy.physics.mechanics.lagrange import LagrangesMethod
   >>> from sympy.physics.mechanics.pathway import LinearPathway

Define Variables
================

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs
   :nofigs:

   >>> M, m, l, k1, k2, c1, g = sm.symbols('M, m, l, k1, k2, c1, g')
   >>> q1, q2 = me.dynamicsymbols('q1 q2')
   >>> q1d, q2d = me.dynamicsymbols('q1 q2', 1)
   >>> u1, u2 = me.dynamicsymbols('u1 u2')
   >>> u1d, u2d = me.dynamicsymbols('u1 u2', 1)

- :math:`q1`: Generalized coordinate representing the position of the Duffing oscillator
- :math:`q2`: Generalized coordinate representing the angle of the pendulum
- :math:`q1d`: First time derivative of `q1`, representing the velocity of the Duffing oscillator
- :math:`q2d`: First time derivative of `q2`, representing the angular velocity of the pendulum
- :math:`u1`: Generalized speed associated with the Duffing oscillator
- :math:`u2`: Generalized speed associated with the pendulum
- :math:`u1d`: First time derivative of `u1`, representing the acceleration of the Duffing oscillator
- :math:`u2d`: First time derivative of `u2`, representing the angular acceleration of the pendulum

Define Kinematics
=================

Define all the reference frames and points.

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs
   :nofigs:

   >>> # Define reference frames
   >>> N = me.ReferenceFrame('N')
   >>> A = N.orientnew('A', 'axis', (q1, N.y)) # Duffing Oscillator Block
   >>> B = N.orientnew('B', 'axis', (q2, N.z)) # Pendulum

The angular velocity of the pendulum in the reference frame is:

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs
   :nofigs:

   >>> B.set_ang_vel(N, u2 * N.z)

Locations and velocities of the Duffing Oscillator block and the pendulum are:

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs
   :nofigs:

   >>> O = me.Point('O') # fixed point in the inertial reference frame
   >>> Block = O.locatenew('Block', q1 * N.y)
   >>> Pendulum = Block.locatenew('Pendulum', l * B.y)

   >>> O.set_vel(N, 0)
   >>> Block.set_vel(N, u1 * N.y)
   >>> Pendulum.v2pt_theory(Block, N, B)
   u1(t)*N.y - l*u2(t)*B.x

   >>> ParBlock = me.Particle('ParBlock', Block, M)
   >>> ParPendulum = me.Particle('ParPendulum', Pendulum, m)

Define Force, Energy
====================

We obtain the Duffing force using the `DuffingSpring` actuator from ``sympy.physics.mechanics.actuator``.
This force will be used to calculate the potential energy of the Duffing Oscillator block.

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs
   :nofigs:

   >>> # Define the Duffing spring force
   >>> pathway = LinearPathway(O, Block)
   >>> duffing_spring = DuffingSpring(k1, k2, pathway, 0)
   >>> duffing_force = -duffing_spring.force
   >>> duffing_force
   k1*sqrt(q1(t)**2) + k2*(q1(t)**2)**(3/2)

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs
   :nofigs:

   >>> # Define Rayleigh dissipation
   >>> D = (1/2) * c1 * q1d**2

In relation to Lagrange's Method, we derive both the kinetic and potential energies of the system.

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs
   :nofigs:

   >>> # Kinetic Energy
   >>> KBlock = (1/2) * M * Block.vel(N).dot(Block.vel(N))
   >>> KPendulum = (1/2) * m * Pendulum.vel(N).dot(Pendulum.vel(N))
   >>> Kinetic = KBlock + KPendulum
   >>> Kinetic
   0.5*M*u1(t)**2 + 0.5*m*(l**2*u2(t)**2 - 2*l*u1(t)*u2(t)*sin(q2(t)) + u1(t)**2)

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs
   :nofigs:

   >>> # Potential Energy
   >>> PBlock = sm.integrate(duffing_force, q1)
   >>> PPendulum = m * g * (l - l/(sm.sqrt(1+q2**2/q1**2)))
   >>> Potential = PBlock + PPendulum
   >>> Potential
   g*m*(l - l/sqrt(1 + q2(t)**2/q1(t)**2)) + k1*sqrt(q1(t)**2)*q1(t)/2 + k2*(q1(t)**2)**(3/2)*q1(t)/4

Lagrange's Method
=================

With the problem setup, the Lagrangian can be calculated, and the equations of motion formed.

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs
   :nofigs:

   >>> L = Kinetic - Potential
   >>> L
   0.5*M*u1(t)**2 - g*m*(l - l/sqrt(1 + q2(t)**2/q1(t)**2)) - k1*sqrt(q1(t)**2)*q1(t)/2 - k2*(q1(t)**2)**(3/2)*q1(t)/4 + 0.5*m*(l**2*u2(t)**2 - 2*l*u1(t)*u2(t)*sin(q2(t)) + u1(t)**2)

   >>> FL = [(Block, duffing_force * N.y + D * N.y), (Pendulum, - m * g * N.y)] # [(point, the force acting on the particle)]
   >>> LM = LagrangesMethod(L, [q1, q2], forcelist = FL, frame = N)
   >>> LM.form_lagranges_equations()
   Matrix([
    [-g*l*m*q2(t)**2/((1 + q2(t)**2/q1(t)**2)**(3/2)*q1(t)**3) + k1*sqrt(q1(t)**2) + k2*(q1(t)**2)**(3/2)],
    [              g*l*m*q2(t)/((1 + q2(t)**2/q1(t)**2)**(3/2)*q1(t)**2) + 1.0*l*m*u1(t)*u2(t)*cos(q2(t))]])

References
==========

P. Brzeskia, P. Perlikowskia, S. Yanchukb, T. Kapitaniaka,
The dynamics of the pendulum suspended on the forced Duffing oscillator,
Journal of Sound and Vibration, 2012, https://doi.org/10.48550/arXiv.1202.5937.
