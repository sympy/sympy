.. _duffing-spring-example:

===========================
Duffing Spring Example
===========================

:obj:`sympy.physics.mechanics` provides features to model 
mechanical systems and simulate their dynamics. 
In this tutorial, we will demonstrate the use of this package 
by modeling and simulating the behavior of a Duffing Spring system, 
a type of nonlinear oscillator.

The Duffing Spring system is characterized by a restoring force 
that includes both linear and cubic terms. 
The equation of motion for the Duffing Spring.

Model Description
=================

The Duffing Spring system is a type of nonlinear oscillator 
with a restoring force that includes both linear and cubic terms. 
Assuming mass m = 1 for simplicity, the equation of motion is given by:

::math::`\ddot{x} + \beta x + \alpha x^3 = 0`

Here, :math:`x` is the displacement from the equilibrium position, 
:math:`\beta` is the linear stiffness coefficient, 
and :math:`\alpha` is the coefficient of the nonlinear cubic term.

.. plot::
   :format: doctest
   :include-source: True
   :context: reset
   :nofigs:

   >>> import sympy as sm
   >>> from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point

Define Variables
================

We define the necessary variables and parameters for 
the Duffing Spring problem using SymPy:

.. plot::
   :format: doctest
   :include-source: True
   :context: reset
   :nofigs:

   >>> import sympy as sm
   >>> from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
   >>> t = sm.symbols('t')
   >>> x = dynamicsymbols('x')
   >>> xdot = x.diff(t)
   >>> xddot = xdot.diff(t)
   >>> alpha, beta = sm.symbols('alpha beta')
   >>> N = ReferenceFrame('N')
   >>> O = Point('O')
   >>> P = Point('P')
   >>> P.set_vel(N, xdot * N.x)
   >>> force_expression = -beta * x - alpha * x**3
   >>> equation_of_motion = sm.Eq(xddot, force_expression)

Simulate the System
===================

Let's numerically solve the Duffing Spring equation using 
`solve_ivp` from `scipy.integrate`, which provides us with 
the system's dynamics over a specific time span.

.. plot::
   :include-source: True
   :context: reset
   :nofigs:

   >>> import sympy as sm
   >>> from sympy.physics.mechanics import dynamicsymbols
   >>> from scipy.integrate import solve_ivp
   >>> import numpy as np
   >>> def duffing_oscillator(t, y, alpha_val, beta_val):
   ...     x, xdot = y
   ...     xddot = -beta_val * x - alpha_val * x**3
   ...     return [xdot, xddot]
   >>> alpha_val = -0.5
   >>> beta_val = 1.0
   >>> initial_conditions = [0, 1]
   >>> t_span = (0, 20)
   >>> t_eval = np.linspace(t_span[0], t_span[1], 400)
   >>> solution = solve_ivp(duffing_oscillator, t_span, initial_conditions, args=(alpha_val, beta_val), t_eval=t_eval, method='RK45')
   >>> time = solution.t
   >>> displacements = solution.y[0]
   >>> velocities = solution.y[1]

Visualize the System
====================

We can plot the displacement and velocity over time and also 
create a phase space plot which is a plot of velocity vs displacement.

.. plot::
   :include-source: True
   :context: reset

   >>> import matplotlib.pyplot as plt
   >>> import numpy as np
   >>> from scipy.integrate import solve_ivp
   >>> def duffing_oscillator(t, y, alpha_val, beta_val):
   ...     x, xdot = y
   ...     xddot = -beta_val * x - alpha_val * x**3
   ...     return [xdot, xddot]
   >>> alpha_val = -0.5
   >>> beta_val = 1.0
   >>> initial_conditions = [0, 1]
   >>> t_span = (0, 20)
   >>> t_eval = np.linspace(t_span[0], t_span[1], 400)
   >>> solution = solve_ivp(duffing_oscillator, t_span, initial_conditions, args=(alpha_val, beta_val), t_eval=t_eval, method='RK45')
   >>> time = solution.t
   >>> displacements = solution.y[0]
   >>> velocities = solution.y[1]
   >>> plt.figure(figsize=(12, 6))
   >>> plt.subplot(2, 1, 1)
   >>> plt.plot(time, displacements, label='Displacement')
   >>> plt.title('Displacement Over Time')
   >>> plt.xlabel('Time (seconds)')
   >>> plt.ylabel('Displacement (meters)')
   >>> plt.grid(True)
   >>> plt.legend()
   >>> plt.subplot(2, 1, 2)
   >>> plt.plot(time, velocities, color='r', label='Velocity')
   >>> plt.title('Velocity Over Time')
   >>> plt.xlabel('Time (seconds)')
   >>> plt.ylabel('Velocity (meters/second)')
   >>> plt.grid(True)
   >>> plt.legend()
   >>> plt.tight_layout()
   >>> plt.show()

Phase Space Plot
================
.. plot::
   :include-source: True
   :context: reset

   >>> import matplotlib.pyplot as plt
   >>> import numpy as np
   >>> from scipy.integrate import solve_ivp
   >>> def duffing_oscillator(t, y, alpha_val, beta_val):
   ...     x, xdot = y
   ...     xddot = -beta_val * x - alpha_val * x**3
   ...     return [xdot, xddot]
   >>> alpha_val = -0.5
   >>> beta_val = 1.0
   >>> initial_conditions = [0, 1]
   >>> t_span = (0, 20)
   >>> t_eval = np.linspace(t_span[0], t_span[1], 400)
   >>> solution = solve_ivp(duffing_oscillator, t_span, initial_conditions, args=(alpha_val, beta_val), t_eval=t_eval, method='RK45')
   >>> time = solution.t
   >>> displacements = solution.y[0]
   >>> velocities = solution.y[1]
   >>> plt.figure(figsize=(6, 6))
   >>> plt.plot(displacements, velocities, label='Phase Space')
   >>> plt.title('Phase Space Plot')
   >>> plt.xlabel('Displacement (meters)')
   >>> plt.ylabel('Velocity (meters/second)')
   >>> plt.grid(True)
   >>> plt.legend()
   >>> plt.show()

Parameter Exploration
=====================

To further understand the dynamics, let's vary parameters like alpha and beta 
and observe how the system's behavior changes.

.. plot::
   :include-source: True
   :context: reset

   >>> import matplotlib.pyplot as plt
   >>> import numpy as np
   >>> from scipy.integrate import solve_ivp
   >>> def duffing_oscillator(t, y, alpha_val, beta_val):
   ...     x, xdot = y
   ...     xddot = -beta_val * x - alpha_val * x**3
   ...     return [xdot, xddot]
   >>> alpha_values = [-1, 0, 1]
   >>> beta_values = [0.5, 1, 1.5]
   >>> initial_conditions = [0, 1]
   >>> t_span = np.linspace(0, 20, 400)
   >>> fig, axs = plt.subplots(len(alpha_values), len(beta_values), figsize=(15, 10))
   >>> for i, alpha_val in enumerate(alpha_values):
   ...     for j, beta_val in enumerate(beta_values):
   ...         solution = solve_ivp(duffing_oscillator, [t_span[0], t_span[-1]], initial_conditions, args=(alpha_val, beta_val), t_eval=t_span)
   ...         axs[i, j].plot(solution.t, solution.y[0])
   ...         axs[i, j].set_title(f'alpha = {alpha_val}, beta = {beta_val}')
   ...         axs[i, j].set_xlabel('Time (s)')
   ...         axs[i, j].set_ylabel('Displacement (m)')
   >>> plt.tight_layout()
   >>> plt.show()

Longer Simulation
=================

Let's extend the time span to observe long-term behavior.

.. plot::
   :include-source: True
   :context: reset

   >>> import matplotlib.pyplot as plt
   >>> import numpy as np
   >>> from scipy.integrate import solve_ivp
   >>> def duffing_oscillator(t, y, alpha_val, beta_val):
   ...     x, xdot = y
   ...     xddot = -beta_val * x - alpha_val * x**3
   ...     return [xdot, xddot]
   >>> alpha_chaos = -1
   >>> beta_chaos = 1.5
   >>> initial_conditions = [0, 1]
   >>> t_span_long = np.linspace(0, 100, 1000)
   >>> long_solution = solve_ivp(duffing_oscillator, [t_span_long[0], t_span_long[-1]], initial_conditions, args=(alpha_chaos, beta_chaos), t_eval=t_span_long)
   >>> plt.figure(figsize=(12, 6))
   >>> plt.plot(long_solution.t, long_solution.y[0], label='Displacement')
   >>> plt.plot(long_solution.t, long_solution.y[1], label='Velocity')
   >>> plt.title('Long-term Behavior')
   >>> plt.xlabel('Time (s)')
   >>> plt.ylabel('Displacement and Velocity')
   >>> plt.legend()
   >>> plt.grid(True)
   >>> plt.show()

Energy Plot
===========

Let's add a plot for the total mechanical energy (kinetic + potential) over time 
to check energy conservation.

.. plot::
   :include-source: True
   :context: reset

   >>> import matplotlib.pyplot as plt
   >>> def total_energy(t, y, alpha, beta):
   ...     kinetic = 0.5 * y[1]**2
   ...     potential = 0.5 * beta * y[0]**2 + (alpha / 4) * y[0]**4
   ...     return kinetic + potential
   >>> alpha_chaos = -1
   >>> beta_chaos = 1.5
   >>> initial_conditions = [0, 1]
   >>> t_span_long = np.linspace(0, 100, 1000)
   >>> long_solution = solve_ivp(duffing_oscillator, [t_span_long[0], t_span_long[-1]], initial_conditions, args=(alpha_chaos, beta_chaos), t_eval=t_span_long)
   >>> energy = [total_energy(t, [y0, y1], alpha_chaos, beta_chaos) for t, y0, y1 in zip(long_solution.t, long_solution.y[0], long_solution.y[1])]
   >>> plt.figure(figsize=(12, 6))
   >>> plt.plot(long_solution.t, energy, label='Total Energy')
   >>> plt.title('Total Mechanical Energy Over Time')
   >>> plt.xlabel('Time (s)')
   >>> plt.ylabel('Energy')
   >>> plt.legend()
   >>> plt.grid(True)
   >>> plt.show()

Analytical Solutions
====================

Let's compare with analytical solutions/results from the literature for validation.

.. plot::
   :include-source: True
   :context: reset

   >>> import matplotlib.pyplot as plt
   >>> import numpy as np
   >>> beta = 1
   >>> alpha_values = [-1, 0, 1]
   >>> x = np.linspace(-2, 2, 400)
   >>> plt.figure(figsize=(6, 6))
   >>> for alpha in alpha_values:
   ...     F = -beta * x - alpha * x**3
   ...     plt.plot(x, F, label=f'α = {alpha}', linewidth=2)
   >>> plt.title('Duffing Oscillator Restoring Force')
   >>> plt.xlabel('Displacement (x)')
   >>> plt.ylabel('Force (F)')
   >>> plt.axhline(0, color='black', linewidth=0.5)
   >>> plt.axvline(0, color='black', linewidth=0.5)
   >>> plt.grid(True)
   >>> plt.legend(title='Parameter α')
   >>> plt.show()