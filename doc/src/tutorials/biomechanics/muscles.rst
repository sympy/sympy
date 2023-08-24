.. _biomechanics-tutorial:

=====================================
Introduction to Biomechanics Modeling
=====================================

:obj:`~sympy.physics._biomechanics` provides features to enhance models created
with :obj:`~sympy.physics.mechanics` with force producing elements that model
muscles and other biomechanical components. In this tutorial, we will introduce
the features of this package.

A Simple Musculotendon Model
============================

To demonstrate a muscle's effect on a simple system, we can model a particle of
mass :math:`m` under the influence of gravity with a muscle pulling the mass
against gravity. The mass :math:`m` has a single generalized coordinate
:math:`q` and generalized speed :math:`u` to describe its position and motion.
The following code establishes the kinematics and gravitational force and an
associated particle::

   >>> import sympy as sm
   >>> import sympy.physics.mechanics as me

   >>> q, u = me.dynamicsymbols('q, u')
   >>> m, g = sm.symbols('m, g')

   >>> N = me.ReferenceFrame('N')
   >>> O, P = sm.symbols('O, P', cls=me.Point)

   >>> P.set_pos(O, q*N.x)
   >>> O.set_vel(N, 0)
   >>> P.set_vel(N, u*N.x)

   >>> gravity = me.Force(P, m*g*N.x)

   >>> block = me.Particle('block', P, m)

SymPy Biomechanics includes musculotendon actuator models. Here we will use a
specific musculotendon model implementation. A musculotendon actuator is
instantiated with two input components, the pathway and the activation dynamics
model. The actuator must act along a pathway that connects the origin and
insertion points of the muscle. Our origin will attach to the fixed point
:math:`O` and insert on the moving particle :math:`P`.

::

   >>> from sympy.physics.mechanics.pathway import LinearPathway

   >>> muscle_pathway = LinearPathway(O, P)

A pathway has attachment points::

   >>> muscle_pathway.attachments
   (O, P)

TODO : note the sign conventions

and knows the length between the end attachment points as well as the relative
speed between the two attachment points::

   >>> muscle_pathway.length
   sqrt(q(t)**2)
   >>> muscle_pathway.extension_velocity
   sqrt(q(t)**2)*Derivative(q(t), t)/q(t)

Finally, the pathway can determine the forces acting on the two attachment
points give a force magnitude::

   >>> muscle_pathway.to_loads(m*g)
   [(O, - g*m*q(t)/sqrt(q(t)**2)*N.x), (P, g*m*q(t)/sqrt(q(t)**2)*N.x)]

The activation dynamics model represents a set of algebraic or ordinary
differential equations that relate the muscle excitation to the muscle
activation. In our case, we will use a first order ordinary differential
equation that gives a smooth, but delayed activation :math:`a(t)` from the
excitation :math:`e(t)`.

TODO : We could plot dadt as a function of a for different e from 0 to 1.

::

   >>> from sympy.physics._biomechanics import FirstOrderActivationDeGroote2016
   >>> muscle_activation = FirstOrderActivationDeGroote2016.with_defaults('muscle')

The activation model has a state variable, input variable, and some constant
parameters::

   >>> muscle_activation.x
   Matrix([[a_muscle(t)]])
   >>> muscle_activation.r
   Matrix([[e_muscle(t)]])
   >>> muscle_activation.p
   Matrix([
   [0.015],
   [ 0.06],
   [   10]])

These are associated with its first order differential equation::

   >>> muscle_activation.rhs()
   Matrix([[((1/2 - tanh(10*a_muscle(t) - 10*e_muscle(t))/2)/(0.0225*a_muscle(t) + 0.0075) + 16.6666666666667*(3*a_muscle(t)/2 + 1/2)*(tanh(10*a_muscle(t) - 10*e_muscle(t))/2 + 1/2))*(-a_muscle(t) + e_muscle(t))]])

With the pathway and activation dynamics, the musculotendon model created using
them both and needs some parameters to define the muscle and tendon specific
properties. You need to specify the tendon slack length, peak isometric force,
optimal fiber length, maximal fiber velocity, optimal pennation angle, and
fiber damping coefficients.

TODO : How do we know this is a rigid tendon model?

::

   >>> from sympy.physics._biomechanics import MusculotendonDeGroote2016

   >>> F_M_max, l_M_opt, l_T_slack = sm.symbols('F_M_max, l_M_opt, l_T_slack')
   >>> v_M_max, alpha_opt, beta = sm.symbols('v_M_max, alpha_opt, beta')

   >>> muscle = MusculotendonDeGroote2016(
   ...     'muscle',
   ...     muscle_pathway,
   ...     muscle_activation,
   ...     tendon_slack_length=l_T_slack,
   ...     peak_isometric_force=F_M_max,
   ...     optimal_fiber_length=l_M_opt,
   ...     maximal_fiber_velocity=v_M_max,
   ...     optimal_pennation_angle=alpha_opt,
   ...     fiber_damping_coefficient=beta,
   ... )
   ...

TODO : Explain why the rhs() is different for the muscle than the activation.
TODO : Needs explanation about rigid tendon

Because this musculotendon actuator has a rigid tendon model, it has the same
state and ordinary differential equation as the activation model::

   >>> muscle.musculotendon_dynamics
   MusculotendonFormulation.RIGID_TENDON
   >>> muscle.x
   Matrix([[a_muscle(t)]])
   >>> muscle.r
   Matrix([[e_muscle(t)]])
   >>> muscle.p
   Matrix([
   [l_T_slack],
   [  F_M_max],
   [  l_M_opt],
   [  v_M_max],
   [alpha_opt],
   [     beta],
   [    0.015],
   [     0.06],
   [       10]])
   >>> muscle.rhs()
   Matrix([[(-0.5625*a_muscle(t)**3*tanh(10*a_muscle(t) - 10*e_muscle(t)) - 0.5625*a_muscle(t)**3 + 0.5625*a_muscle(t)**2*e_muscle(t)*tanh(10*a_muscle(t) - 10*e_muscle(t)) + 0.5625*a_muscle(t)**2*e_muscle(t) - 0.375*a_muscle(t)**2*tanh(10*a_muscle(t) - 10*e_muscle(t)) - 0.375*a_muscle(t)**2 + 0.375*a_muscle(t)*e_muscle(t)*tanh(10*a_muscle(t) - 10*e_muscle(t)) + 0.375*a_muscle(t)*e_muscle(t) + 0.9375*a_muscle(t)*tanh(10*a_muscle(t) - 10*e_muscle(t)) - 1.0625*a_muscle(t) - 0.9375*e_muscle(t)*tanh(10*a_muscle(t) - 10*e_muscle(t)) + 1.0625*e_muscle(t))/(0.045*a_muscle(t) + 0.015)]])

The musculotendon provides the extra ordinary differential equations as well as
the muscle specific forces applied to the pathway::

   >>> muscle_loads = muscle.to_loads()
   >>> muscle_loads[0]
   (O, F_M_max*(beta*(-l_T_slack + sqrt(q(t)**2))*sqrt(q(t)**2)*Derivative(q(t), t)/(v_M_max*sqrt(l_M_opt**2*sin(alpha_opt)**2 + (-l_T_slack + sqrt(q(t)**2))**2)*q(t)) + a_muscle(t)*FiberForceLengthActiveDeGroote2016(sqrt(l_M_opt**2*sin(alpha_opt)**2 + (-l_T_slack + sqrt(q(t)**2))**2)/l_M_opt, 0.814, 1.06, 0.162, 0.0633, 0.433, 0.717, -0.0299, 1/5, 1/10, 1, 0.354, 0)*FiberForceVelocityDeGroote2016((-l_T_slack + sqrt(q(t)**2))*sqrt(q(t)**2)*Derivative(q(t), t)/(v_M_max*sqrt(l_M_opt**2*sin(alpha_opt)**2 + (-l_T_slack + sqrt(q(t)**2))**2)*q(t)), -0.318, -8.149, -0.374, 0.886) + FiberForceLengthPassiveDeGroote2016(sqrt(l_M_opt**2*sin(alpha_opt)**2 + (-l_T_slack + sqrt(q(t)**2))**2)/l_M_opt, 3/5, 4))*q(t)/sqrt(q(t)**2)*N.x)
   >>> muscle_loads[1]
   (P, - F_M_max*(beta*(-l_T_slack + sqrt(q(t)**2))*sqrt(q(t)**2)*Derivative(q(t), t)/(v_M_max*sqrt(l_M_opt**2*sin(alpha_opt)**2 + (-l_T_slack + sqrt(q(t)**2))**2)*q(t)) + a_muscle(t)*FiberForceLengthActiveDeGroote2016(sqrt(l_M_opt**2*sin(alpha_opt)**2 + (-l_T_slack + sqrt(q(t)**2))**2)/l_M_opt, 0.814, 1.06, 0.162, 0.0633, 0.433, 0.717, -0.0299, 1/5, 1/10, 1, 0.354, 0)*FiberForceVelocityDeGroote2016((-l_T_slack + sqrt(q(t)**2))*sqrt(q(t)**2)*Derivative(q(t), t)/(v_M_max*sqrt(l_M_opt**2*sin(alpha_opt)**2 + (-l_T_slack + sqrt(q(t)**2))**2)*q(t)), -0.318, -8.149, -0.374, 0.886) + FiberForceLengthPassiveDeGroote2016(sqrt(l_M_opt**2*sin(alpha_opt)**2 + (-l_T_slack + sqrt(q(t)**2))**2)/l_M_opt, 3/5, 4))*q(t)/sqrt(q(t)**2)*N.x)

These loads are made up of various functions that describe the length and
velocity relationships to the fiber force.

Now that we have the forces that the muscles and tendons produce the equations
of motion of the system can be formed with, for example, Kanes Method::

   >>> kane = me.KanesMethod(N, (q,), (u,), kd_eqs=(u - q.diff(),))
   >>> Fr, Frs = kane.kanes_equations((block,), (muscle_loads + [gravity]))

The equations of motion are made up of the kinematical differential equation,
the dynamical differential equation (Newton's Second Law), and the muscle
activation differential equation. The explicit form of each can be formed like
so::

   >>> dqdt = u
   >>> dudt = kane.forcing[0]/m
   >>> dadt = muscle.rhs()[0]

We can now create a numerical function that evaluates the equations of motion
given the state, inputs, and constant parameters. Start by listing each
symbolically::

   >>> a = muscle.a
   >>> e = muscle.e
   >>> state = [q, u, a]
   >>> inputs = [e]
   >>> constants = [m, g, F_M_max, l_M_opt, l_T_slack, v_M_max, alpha_opt, beta]

Then the numerical function is::

   >>> eval_eom = sm.lambdify((state, inputs, constants), (dqdt, dudt, dadt))

It will additionally be interesting to numerically evaluate the muscle force,
so create a function for it too::

   >>> force = muscle.force.xreplace({q.diff(): u})
   >>> eval_force = sm.lambdify((state, constants), force)

To test these functions we need some suitable numerical values. This muscle
will be able to produce a maximum force of 10 N to lift a mass of 0.5 kg::

   >>> import numpy as np
   >>> p_vals = np.array([
   ...     0.5,  # m [kg]
   ...     9.81,  # g [m/s/s]
   ...     10.0,  # F_M_max
   ...     0.18,  # l_M_opt, length of muscle at which max force is produced
   ...     0.17,  # l_T_slack, always fixed (rigid tendon)
   ...     10.0,  # v_M_max
   ...     0.0,  # alpha_opt
   ...     0.1,  # beta
   ... ])
   ...

Our tendon is rigid, so the length of the muscle will be :math:`q-l_T_slack`
and we want to give an initial muscle length near its force producing peak, so
we choose :math:`q_0=l_M_opt + l_T_slack`::

   >>> x_vals = np.array([
   ...     p_vals[3] + p_vals[4],  # q [m]
   ...     0.0,  # u [m/s]
   ...     0.0,  # a [?]
   ... ])
   ...

We can set the excitation to zero to test the numerical functions::

   >>> r_vals = np.array([
   ...     0.0,  # e
   ... ])
   ...
   >>> eval_eom(x_vals, r_vals, p_vals)
   (0.0, 9.81, 0.0)
   >>> eval_force(x_vals, p_vals)
   1.4499681738213515e-16

The two functions work so we can now simulate this system to see if and how the
muscle lifts the mass::

   >>> def eval_rhs(t, x):
   ...
   ...     r = np.array([1.0])
   ...
   ...     return eval_eom(x, r, p_vals)
   ...

   >>> from scipy.integrate import solve_ivp
   >>> t0, tf = 0.0, 10.0
   >>> times = np.linspace(t0, tf, num=1001)
   >>> sol = solve_ivp(eval_rhs,
   ...                 (t0, tf),
   ...                 x_vals, t_eval=times)
   ...
   >>> import matplotlib.pyplot as plt
   >>> fig, axes = plt.subplots(4, 1, sharex=True)
   >>> axes[0].plot(sol.t, sol.y[0] - p_vals[4], label='length of muscle')
   >>> axes[1].plot(sol.t, sol.y[1], label=state[1])
   >>> axes[2].plot(sol.t, sol.y[2], label=state[2])
   >>> axes[3].plot(sol.t, eval_force(sol.y, p_vals).T, label='force')
   >>> axes[0].legend(), axes[1].legend(), axes[2].legend(), axes[3].legend()
