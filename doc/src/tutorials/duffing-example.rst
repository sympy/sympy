.. _duffing-example:

===========================
Duffing Spring Example
===========================

:obj:`sympy.physics.mechanics` provides features to model
mechanical systems and simulate their dynamics.
In this tutorial, we will demonstrate the use of this package
by modeling and simulating the behavior of a Duffing Spring system,
a type of nonlinear oscillator

Model Description
=================

The Duffing Spring system is a type of nonlinear oscillator
with a restoring force that includes both linear and cubic terms.
Assuming mass m = 1 for simplicity, the equation of motion is given by:

.. math::
   \ddot{x} + \beta x + \alpha x^3 = 0

Here, :math:`x` is the displacement from the equilibrium position,
:math:`\beta` is the linear stiffness coefficient,
and :math:`\alpha` is the coefficient of the nonlinear cubic term.

.. plot::
   :format: doctest
   :include-source: True
   :context: reset
   :nofigs:

   >>> import numpy as np
   >>> from scipy.integrate import solve_ivp
   >>> import sympy as sm
   >>> import sympy.physics.mechanics as me
   >>> from sympy.physics.mechanics.actuator import DuffingSpring
   >>> from sympy.physics.mechanics.pathway import LinearPathway

Define variables and parameters for the Duffing Spring
======================================================

We define the necessary variables and parameters for
the Duffing Spring problem using SymPy.

.. plot::
   :format: doctest
   :include-source: True
   :context:
   :nofigs:

   >>> t = sm.symbols('t')
   >>> q = me.dynamicsymbols('q')
   >>> l = sm.symbols('l')
   >>> alpha, beta = sm.symbols('alpha beta')
   >>> N = me.ReferenceFrame('N')
   >>> O = me.Point('O')
   >>> P = me.Point('P')
   >>> P.set_pos(O, q * N.x)
   >>> P.set_vel(N, q.diff(t) * N.x)
   >>> pathway = LinearPathway(O, P)

Now, let's calculate the force and equation of motion using DuffingSpring.

.. plot::
   :format: doctest
   :include-source: True
   :context:
   :nofigs:

   >>> duffing_spring = DuffingSpring(linear_stiffness=beta, nonlinear_stiffness=alpha, pathway=pathway, equilibrium_length=l)
   >>> force_expr = duffing_spring.force
   >>> force_func = sm.lambdify((q, alpha, beta, l), force_expr)

We also numerically solve the Duffing Spring equation using solve_ivp from
scipy.integrate, which provides up with the system's dynamics
over a specific time span.

.. plot::
   :format: doctest
   :include-source: True
   :context:
   :nofigs:

   >>> def duffing_oscillator(t, y, alpha_val, beta_val, l_val):
   ...     q_val, qdot_val = y
   ...     force = force_func(q_val, alpha_val, beta_val, l_val)
   ...     mass = 1
   ...     qddot_val = force / mass
   ...     return [qdot_val, qddot_val]

   >>> # Parameters for the simulation
   >>> alpha_val = -0.1
   >>> beta_val = 1.0
   >>> l_val = 10
   >>> initial_conditions = [l_val, 1]  # [initial displacement, initial velocity]
   >>> t_span = (0, 100)
   >>> t_eval = np.linspace(t_span[0], t_span[1], 1000)

   >>> # Solve the differential equation
   >>> solution = solve_ivp(duffing_oscillator, t_span, initial_conditions, args=(alpha_val, beta_val, l_val), t_eval=t_eval, method='RK45')

   >>> # Extract the time and displacement from the solution
   >>> time = solution.t
   >>> displacements = solution.y[0]
   >>> velocities = solution.y[1]

Visualize the System
====================

We can plot the displacement and velocity over time.

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs

   >>> import matplotlib.pyplot as plt

   >>> fig, ax = plt.subplots()
   >>> _ = ax.plot(time, displacements, label='Time vs Displacement')
   >>> _ = ax.set_xlabel('Time (s)')
   >>> _ = ax.set_ylabel('Displacement (m)')
   >>> plt.show()

   >>> fig, ax = plt.subplots()
   >>> _ = ax.plot(time, velocities, label='Time vs Velocity')
   >>> _ = ax.set_xlabel('Time (s)')
   >>> _ = ax.set_ylabel('Velocity (m/s)')
   >>> plt.show()

Phase Space Plot
================

We can also create a phase space plot which is a plot of velocity vs displacement.
This phase space plot graphs velocity against displacement, visually representing
the system's state over time.

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs

   >>> import matplotlib.pyplot as plt

   >>> fig, ax = plt.subplots()
   >>> _ = ax.plot(displacements, velocities, label='Phase Space')
   >>> _ = ax.set_xlabel('Displacement (m)')
   >>> _ = ax.set_ylabel('Velocity (m/s)')

Parameter Exploration
=====================

To further understand the dynamics, let's vary parameters like alpha and beta
and observe how the system's behavior changes.
Multiple subplots explore variations in displacement over time for different values of alpha (:math:`\alpha`) and beta (:math:`\beta`),
demonstrating how the system's response varies with changes in stiffness and nonlinearity parameters.

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs

   >>> import matplotlib.pyplot as plt

   >>> alpha_values = [-1, 0, 1]
   >>> beta_values = [0.5, 1, 1.5]
   >>> fig, axs = plt.subplots(len(alpha_values), len(beta_values), figsize=(15, 10))
   >>> for i, alpha_val in enumerate(alpha_values):
   ...     for j, beta_val in enumerate(beta_values):
   ...         solution = solve_ivp(duffing_oscillator, t_span, initial_conditions, args=(alpha_val, beta_val, l_val), t_eval=t_eval)
   ...         axs[i, j].plot(solution.t, solution.y[0])
   ...         axs[i, j].set_title(f'alpha = {alpha_val}, beta = {beta_val}')
   ...         axs[i, j].set_xlabel('Time (s)')
   ...         axs[i, j].set_ylabel('Displacement (m)')
   >>> plt.tight_layout()

Energy Plot
===========

Let's add a plot for the total mechanical energy (kinetic + potential) over time
to check energy conservation.

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs

   >>> import matplotlib.pyplot as plt

   >>> def total_energy(q, qdot, alpha, beta):
   ...     kinetic = 0.5 * qdot ** 2
   ...     potential = 0.5 * beta * q**2 + 0.25 * alpha * q**4
   ...     return kinetic + potential
   >>> energy = [total_energy(q, qdot, alpha_chaos, beta_chaos) for q, qdot in zip(long_solution.y[0], long_solution.y[1])]

   >>> fig, ax = plt.subplots()
   >>> _ = ax.plot(long_solution.t, energy, label='Total Energy')
   >>> _ = ax.set_xlabel('Time (s)')
   >>> _ = ax.set_ylabel('Energy')

This plot tracks the total mechanical energy of the system over time.

Analytical Solutions
====================

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs

   >>> import matplotlib.pyplot as plt

   >>> beta_val = 1
   >>> l_val = 10
   >>> alpha_values = [-1, 0, 1]
   >>> x_vals = np.linspace(0, 20, 400)

   >>> fig, ax = plt.subplots(figsize=(8, 5))
   >>> for alpha_val in alpha_values:
   ...     force_vals = [force_func(x, alpha_val, beta_val, l_val) for x in x_vals]
   ...     ax.plot(x_vals, force_vals, label=f'α = {alpha_val}', linewidth=2)
