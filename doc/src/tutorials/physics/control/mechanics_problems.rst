.. _mechanics_problems-physics:

=============================================
Mechanics Problems using StateSpace
=============================================

Below are some Mechanics problems that can be solved using
StateSpace.

Example 1
---------

        .. image:: Mechanics_Problems_Q1.svg
           :align: center

A spring-mass-damping system can be modeled using a mass (m), a spring with a constant (k), and a damper with a damping coefficient (b). The spring force is proportional to the displacement of the mass, and the damping force is proportional to the velocity of the mass.
Find the frequency response of the system.
The free-body diagram for this system is shown below:

        .. image:: Mechanics_Problems_Q1_FBD.svg
           :align: center

The equation of motion for the mass-spring-damper system is given by:

.. math::

   m\ddot{x} + b\dot{x} + kx = F(t)

where:

* :math:`x` is the displacement of the mass,
* :math:`\dot{x}` is the velocity of the mass,
* :math:`\ddot{x}` is the acceleration of the mass,
* :math:`F(t)` is the external force applied to the system.

To determine the state-space representation of the mass-spring-damper system, we reduce the second-order differential equation to a set of two first-order differential equations. We choose the position and velocity as our state variables:

.. math::

   x_1 = x \quad \text{and} \quad x_2 = \dot{x}

The state equations become:

.. math::

   \dot{x}_1 = x_2 \\

   \dot{x}_2 = -\frac{k}{m}x_1 - \frac{b}{m}x_2 + \frac{1}{m}F(t)\\

The state-space can be represented by:

.. math::

   \mathbf{A} = \begin{bmatrix} 0 & 1 \\ -\frac{k}{m} & -\frac{b}{m} \end{bmatrix}, \quad
   \mathbf{B} = \begin{bmatrix} 0 \\ \frac{1}{m} \end{bmatrix}, \quad
   \mathbf{C} = \begin{bmatrix} 1 & 0 \end{bmatrix}

The state equation can be written as

.. math::

   \dot{x} =
   \begin{bmatrix}
   \dot{x} \\
   \ddot{x}
   \end{bmatrix}
   =
   \begin{bmatrix}
   0 & 1 \\
   -\frac{k}{m} & -\frac{b}{m}
   \end{bmatrix}
   \begin{bmatrix}
   x \\
   \dot{x}
   \end{bmatrix}
   +
   \begin{bmatrix}
   0 \\
   \frac{1}{m}
   \end{bmatrix}
   F(t)

Using SymPy's Control Systems Toolbox (CST), we can define the state-space representation and convert it to the transfer function.

Solution
^^^^^^^^

The following code demonstrates how to define the state-space representation of the spring-mass-damper system and convert it to a transfer function using SymPy:

    >>> from sympy import symbols, Matrix
    >>> from sympy.physics.control import *

    Define the variables

    >>> m, k, b = symbols('m k b')

    Define the state-space matrices

    >>> A = Matrix([[0, 1], [-k/m, -b/m]])
    >>> B = Matrix([[0], [1/m]])
    >>> C = Matrix([[1, 0]])
    >>> D = Matrix([[0]])

    Create the StateSpace model

    >>> ss = StateSpace(A, B, C, D)
    >>> ss
    StateSpace(Matrix([
    [   0,    1],
    [-k/m, -b/m]]), Matrix([
    [  0],
    [1/m]]), Matrix([[1, 0]]), Matrix([[0]]))

    Converting StateSpace to TransferFunction by rewrite method.

    >>> tf = ss.rewrite(TransferFunction)[0][0]
    >>> tf
    TransferFunction(1, b*s + k + m*s**2, s)

References
^^^^^^^^^^
1. `ctms.engin.umich.edu <https://ctms.engin.umich.edu/CTMS/index.php?example=Introduction&section=SystemModeling>`_

Example 2
---------

        .. image:: Mechanics_Problems_Q2.svg
           :align: center

This problem explains how to model a rotaional system to state-space model. The system has input torque `τ_a` and damping effects `B_{r1}` and `B_{r2}`. The system consists of two flywheels connected by a spring, with the angular positions denoted by `θ_1` and `θ_2`.

The energy variables for the rotating system are potential energy stored in springs `1/2 K_r \theta^ 2` and kinetic energy stored in inertial elments `1/2 J \omega ^ 2`.

The **State Variables:** can be written as:

.. math::

   x_1 = \theta_1 \quad \text{(angular position of the first flywheel)}

   x_2 = \dot{\theta}_1 \quad \text{(angular velocity of the first flywheel)}

   x_3 = \dot{\theta_2} \quad \text{(angular velocity of the second flywheel)}

The goal is to find a set of first-order differential equations that describe the system in terms of these state variables.

First, we write the equations of motion for the two flywheels, including the effects of damping.

   1. For the first flywheel (`J_1`):

        .. image:: Mechanics_Problems_Q2_FBD1.svg
           :align: center

      .. math::

         J_1 \ddot{\theta}_1 + B_{r1} \dot{\theta}_1 + K_r\theta_1 - B_{r1} \dot{\theta}_2 = - \tau_a

   2. For the second flywheel (`J_2`):

        .. image:: Mechanics_Problems_Q2_FBD2.svg
           :align: center

      .. math::

         J_2 \ddot{\theta}_2 + (B_{r2} + B_{r1}) \dot{\theta}_2 - B_{r1} \dot{\theta}_1 = 0

   Now we want the equations for the derivates of state variables.

      .. math::

         \dot{x}_1 = \dot{\theta}_1 = x_2

      .. math::

         \dot{x}_2 = \ddot{\theta_1} = \frac{1}{J_1} \left(-\tau_a - B_{r1} \dot{\theta_1} - K_r \theta_1 + B_{r1}\dot{\theta_2} \right)

         \dot{x}_2 = \frac{1}{J_1} \left(-\tau_a - B_{r1} x_2 - K_r x_1 + B_{r1} x_3 \right)

      .. math::

         \dot{x}_3 = \ddot{\theta}_2 = \frac{1}{J_2} \left(- (B_{r2} + B_{r1}) \dot{\theta_2} + B_{r1} \dot{\theta_1} \right)

         \dot{x}_3 = \frac{1}{J_2} \left( - (B_{r2} + B_{r1}) x_3 + B_{r1} x_2 \right)

The state-space model of the system can be expressed in the standard form:

.. math::

   \dot{x} = A x + B u

   y = C x + D u

Where:

- **x** is the state vector:

  .. math::

     x = \begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} = \begin{bmatrix} \theta_1 \\ \dot{\theta}_1 \\ \theta_2 \end{bmatrix}

- **u** is the input torque (`τ_a`).
- **y** is the output angular position (`θ_1`).

The matrices **A**, **B**, **C**, and **D** are defined as follows:

   The **A matrix** represents the relationship between the state variables. It is defined as:

   .. math::

      A = \begin{bmatrix}
      0 & 1 & 0 \\
      -\frac{K_r}{J_1} & -\frac{B_{r1}}{J_1} & \frac{B_{r1}}{J_1} \\
      0 & \frac{B_{r1}}{J_2} x_2 & -\frac{B_{r2} + B_{r1}}{J_2}
      \end{bmatrix}

   The **B matrix** represents the influence of the input torque on the system. It is defined as:

   .. math::

      B = \begin{bmatrix}
      0 \\
      \frac{-1}{J_1} \\
      0
      \end{bmatrix}

   The **C matrix** defines the relationship between the output (y) and the state variables (x). Since we are only interested in the angular position θ₁, the **C matrix** is:

   .. math::

      C = \begin{bmatrix} 1 & 0 & 0 \end{bmatrix}

   The **D matrix** is the direct transmission matrix. Since there is no direct transmission frothe input to the output, **D** is zero:

   .. math::

      D = 0

Solution
^^^^^^^^

   >>> from sympy import symbols, Matrix
   >>> from sympy.physics.control import StateSpace
   >>> K_r, J1, J2, B_r1, B_r2, x2 = symbols('K_r J1 J2 B_r1 B_r2 x2')
   >>> A = Matrix([[0, 1, 0], [-K_r/J1, -B_r1/J1, B_r1/J1], [0, B_r1/J2 * x2, - (B_r2 + B_r1)/J2]])
   >>> B = Matrix([[0], [-1/J1], [0]])
   >>> C = Matrix([[1, 0, 0]])
   >>> ss = StateSpace(A, B, C)
   >>> ss
   StateSpace(Matrix([
   [      0,          1,                 0],
   [-K_r/J1,   -B_r1/J1,           B_r1/J1],
   [      0, B_r1*x2/J2, (-B_r1 - B_r2)/J2]]), Matrix([
   [    0],
   [-1/J1],
   [    0]]), Matrix([[1, 0, 0]]), Matrix([[0]]))

References
^^^^^^^^^^
1. `https://lpsa.swarthmore.edu/ <https://lpsa.swarthmore.edu/Representations/SysRepSS.html#:~:text=Example%3A%20Direct%20Derivation%20of%20State%20Space%20Model%20(Electrical),-Derive%20a%20state&text=The%20input%20is%20ia%20and%20the%20output%20is%20e2.&text=space%20representation%20becomes-,This%20technique%20does%20not%20always%20easily%20yield%20a%20set%20of,Transfer%20functions%20are%20discussed%20elsewhere.>`_
