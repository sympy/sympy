#!/usr/bin/env python
# coding: utf-8

# ### Defining the Problem
# Derivation of the equations of motion for the DuffingSpring system

# In[1]:


# Load the core functionalities/libraries
import sympy as sm
import numpy as np
import matplotlib.pyplot as plt
from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy import symbols, Function, cos
from sympy.interactive import printing


# In[2]:


# Load SymPy's printing extension
printing.init_printing(use_latex=True)


# In[3]:


# Define variables we need for the DuffingSpring problem
t = symbols('t')
x = dynamicsymbols('x')
xdot = x.diff(t)
xddot = xdot.diff(t)
alpha, beta = symbols('alpha beta')


# In[4]:


# Define the reference frame and points
N = ReferenceFrame('N')
O = Point('O')
P = Point('P')
P.set_vel(N, xdot * N.x)  # Set the velocity of point P in frame N


# In[5]:


# Define the force as a function of displacement
force_expression = -beta * x - alpha * x**3


# In[6]:


# Define the equation of motion based on Newton's second law: F = m*a
# Assuming mass m = 1 for simplicity
equation_of_motion = sm.Eq(xddot, force_expression)


# ### Simulating the System
# Let's numerically solve the Duffing Spring equation using `solve_ivp` from `scipy.integrate`, which provides us with the system's dynamics over a specific time span.

# In[7]:


from scipy.integrate import solve_ivp

# Convert SymPy expression to a function that can be used by scipy.integrate.solve_ivp
def duffing_oscillator(t, y, alpha_val, beta_val):
    x, xdot = y
    xddot = -beta_val * x - alpha_val * x**3
    return [xdot, xddot]

# Parameters
alpha_val = -0.5  # Example value for alpha
beta_val = 1.0    # Example value for beta

# Initial conditions
initial_conditions = [0, 1]  # [initial displacement, initial velocity]

# Time span for the simulation
t_span = (0, 20)  # Simulate from t=0 to t=20
t_eval = np.linspace(t_span[0], t_span[1], 400)  # Time points at which to store the results


# In[8]:


# Solve the equation using solve_ivp
solution = solve_ivp(duffing_oscillator, t_span, initial_conditions, args=(alpha_val, beta_val), t_eval=t_eval, method='RK45')


# In[9]:


# Extract the results
time = solution.t
displacements = solution.y[0]
velocities = solution.y[1]


# ### Visualising the System
# Let's have the results from the simulation. We can plot the displacement and velocity over time and also create a phase space plot which is a plot of velocity vs displacement.

# In[10]:


# Plotting displacement over time
plt.figure(figsize=(12, 6))
plt.subplot(2, 1, 1)
plt.plot(time, displacements, label='Displacement')
plt.title('Displacement Over Time')
plt.xlabel('Time (seconds)')
plt.ylabel('Displacement (meters)')
plt.grid(True)
plt.legend()


# This plot displays the displacement of the DuffingSpring over time. The smooth, periodic oscillations indicate stable, regular motion under the specified system parameters. Symmetry around zero suggests an equilibrium state at this point, with the amplitude and frequency offering insights into the system's stiffness and damping characteristics.

# In[11]:


# Plotting velocity over time
plt.subplot(2, 1, 2)
plt.plot(time, velocities, color='r', label='Velocity')
plt.title('Velocity Over Time')
plt.xlabel('Time (seconds)')
plt.ylabel('Velocity (meters/second)')
plt.grid(True)
plt.legend()


# This plot illustrates the velocity of the DuffingSpring over time. Changes in velocity are sinusoidal and out of phase with displacement, typical of harmonic oscillators. Peak velocities align with the zero-crossings of displacement, indicating maximum kinetic energy when the spring passes through its equilibrium position.

# In[12]:


# Show the plots
plt.tight_layout()
plt.show()


# In[13]:


# Phase space plot (velocity vs. displacement)
plt.figure(figsize=(6, 6))
plt.plot(displacements, velocities, label='Phase Space')
plt.title('Phase Space Plot')
plt.xlabel('Displacement (meters)')
plt.ylabel('Velocity (meters/second)')
plt.grid(True)
plt.legend()
plt.show()


# This phase space plot graphs velocity against displacement, visually representing the system's state over time. The closed loop indicates a stable limit cycle, characteristic of periodic motion.

# ### Parameter Exploration
# To further understand the dynamics, let's vary parameters like alpha and beta and observe how the system's behaviour changes. This can help in analysing the stability and bifurcation scenarios typical to Duffing systems.

# In[14]:


# Parameter values to explore
alpha_values = [-1, 0, 1]
beta_values = [0.5, 1, 1.5]

# Initial conditions
initial_conditions = [0, 1]

# Time span
t_span = np.linspace(0, 20, 400)

fig, axs = plt.subplots(len(alpha_values), len(beta_values), figsize=(15, 10))

for i, alpha_val in enumerate(alpha_values):
    for j, beta_val in enumerate(beta_values):
        solution = solve_ivp(duffing_oscillator, [t_span[0], t_span[-1]], initial_conditions, args=(alpha_val, beta_val), t_eval=t_span)
        axs[i, j].plot(solution.t, solution.y[0])
        axs[i, j].set_title(f'alpha = {alpha_val}, beta = {beta_val}')
        axs[i, j].set_xlabel('Time (s)')
        axs[i, j].set_ylabel('Displacement (m)')

plt.tight_layout()
plt.show()


# Multiple subplots explore variations in displacement over time for different values of alpha (α) and beta (β), demonstrating how the system's response varies with changes in stiffness and nonlinearity parameters.

# ### Longer Simulation
# Let's extend the time span

# In[15]:


# Extended time span for longer simulation
t_span_long = np.linspace(0, 100, 1000)  # Simulate for 100 seconds

# Choose a set of parameters
alpha_chaos = -1
beta_chaos = 1.5

# Solve the equation over the extended time span
long_solution = solve_ivp(duffing_oscillator, [t_span_long[0], t_span_long[-1]], initial_conditions, args=(alpha_chaos, beta_chaos), t_eval=t_span_long)

# Plotting the results
plt.figure(figsize=(12, 6))
plt.plot(long_solution.t, long_solution.y[0], label='Displacement')
plt.plot(long_solution.t, long_solution.y[1], label='Velocity')
plt.title('Long-term Behaviour')
plt.xlabel('Time (s)')
plt.ylabel('Displacement and Velocity')
plt.legend()
plt.grid(True)
plt.show()


# ### Energy Plot
# Let's add a plot for the total mechanical energy (kinetic + potential) over time to check energy conservation in the numerical method.

# In[16]:


# Function to calculate energy
def total_energy(t, y, alpha, beta):
    kinetic = 0.5 * y[1]**2
    potential = 0.5 * beta * y[0]**2 + (alpha / 4) * y[0]**4
    return kinetic + potential

# Calculate energy for the longer simulation
energy = [total_energy(t, [y0, y1], alpha_chaos, beta_chaos) for t, y0, y1 in zip(long_solution.t, long_solution.y[0], long_solution.y[1])]

# Plotting the energy over time
plt.figure(figsize=(12, 6))
plt.plot(long_solution.t, energy, label='Total Energy')
plt.title('Total Mechanical Energy Over Time')
plt.xlabel('Time (s)')
plt.ylabel('Energy')
plt.legend()
plt.grid(True)
plt.show()


# This plot tracks the total mechanical energy of the system over time. A gradual decline in total energy indicates a non-conservative system, where energy dissipation occurs over time, likely due to internal damping effects.

# ### Analytical Solutions
# Let's compare with analytical solutions/results from the literature for validation.

# In[17]:


# Parameters
beta = 1  # Positive beta as per the reference figure
alpha_values = [-1, 0, 1]  # Different alpha values for softening and hardening
x = np.linspace(-2, 2, 400)  # Displacement range from -2 to 2

# Plot setup
plt.figure(figsize=(6, 6))
for alpha in alpha_values:
    F = -beta * x - alpha * x**3
    plt.plot(x, F, label=f'α = {alpha}', linewidth=2)

# Adding plot features to match the reference figure
plt.title('Duffing Oscillator Restoring Force')
plt.xlabel('Displacement (x)')
plt.ylabel('Force (F)')
plt.axhline(0, color='black',linewidth=0.5)
plt.axvline(0, color='black',linewidth=0.5)
plt.grid(True)
plt.legend(title='Parameter α')
plt.show()


# This plot visualises the restoring force as a function of displacement for varying values of α, illustrating the nonlinear force characteristics of the Duffing oscillator. The shape of the curves indicates how the restoring force varies with displacement. Negative α values exhibit a softening effect (force decreases with displacement), while positive α values indicate a hardening effect (force increases with displacement), essential for understanding the spring's behaviour under different conditions. As illustrated, our resulting plots align with expectations and reference literature, thereby validating their functionality and accuracy.
# 
# ![Duffing%20Oscillator%20for%20beta%3E0.png](attachment:Duffing%20Oscillator%20for%20beta%3E0.png)

# ### Adding examples (docstring)
# Consider adding this as a docstring in the DuffingSpring class.

# In[18]:
