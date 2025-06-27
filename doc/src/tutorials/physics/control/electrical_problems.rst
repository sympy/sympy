.. _electrical_problems-physics:

=============================================
Electrical Problems using StateSpace
=============================================

The state-space approach is a powerful method used to model and analyze systems in control
theory. Instead of focusing solely on the input-output relationships like the transfer function
approach, the state-space approach represents systems as a set of first-order differential
equations.

The state-space representation of a system can be written as:

        .. math::

            \dot{x}(t) = A x(t) + B u(t) \\
            y(t) = C x(t) + D u(t)


Where :math:`x(t)` is the state vector, :math:`u(t)` is the input vector, :math:`y(t)` is the output vector,
:math:`A`, :math:`B`, :math:`C`, and :math:`D` are matrices that define the system dynamics.

Below are some examples to demonstrate the use of StateSpace to solve Electrical problems.

Example 1
---------

        .. image:: Electrical_Problems_Q1.svg
           :align: center

In a series RLC circuit, we have a resistor :math:`R`, an inductor :math:`L`, and a capacitor :math:`C`
connected in series with an input voltage :math:`v_{in}(t)`. The state variables are the current through
the inductor :math:`i(t)` and the voltage across the capacitor :math:`v_C(t)`.

Applying **Kirchhoff's Voltage Law** (KVL) around the loop in the above diagram gives:

        .. math::

            v_{in}(t) = R \cdot i(t) + L \frac{di(t)}{dt} +  V_C(t)

Where: :math:`v_{in}(t)` is the input voltage, :math:`i(t)` is the current through the inductor and
:math:`v_C(t)` is the voltage across the capacitor.

This equation relates the input voltage to the elements of the RLC circuit.

**Capacitor Voltage Equation**

The voltage across the capacitor can be related to the current by:

        .. math::

            V_C(t) = \frac{1}{C} \int i(t) \, dt

Taking the time derivative of both sides, we obtain the rate of change of the capacitor voltage:

        .. math::

            \dot{v}_C(t) = \frac{d v_C(t)}{dt} = \frac{i(t)}{C}

This equation shows that the rate of change of the capacitor voltage is proportional to the current through the circuit.

From the KVL equation, solving for the derivative of the current gives:

        .. math::

            \frac{di(t)}{dt} = -\frac{R}{L} i(t) - \frac{1}{L} v_C(t) + \frac{1}{L} v_{in}(t)

This is the first-order differential equation that describes the rate of change of the current in terms of the circuit's components and input voltage.

The state-space representation expresses the system in terms of state variables, which are typically the variables that describe the energy stored in the circuit elements (such as current and voltage).

We define the state vector `X(t)` as:

        .. math::

            X(t) = \begin{bmatrix} x_1(t) \\ x_2(t) \end{bmatrix} = \begin{bmatrix} i(t) \\ v_C(t) \end{bmatrix}

Here `x_1(t) = i(t)` is the current through the inductor and `x_2(t) = v_C(t)` is the voltage across the capacitor.

The input vector `U(t)` is the input voltage:

        .. math::

            U(t) = v_{in}(t)

The system of differential equations in terms of the state variables becomes:

1. The derivative of the current:

        .. math::

            \dot{x}_1(t) = -\frac{R}{L} x_1(t) - \frac{1}{L} x_2(t) + \frac{1}{L} v_{in}(t)

2. The derivative of the capacitor voltage:

        .. math::

            \dot{x}_2(t) = \frac{x_1(t)}{C}


The matrices for the series RLC circuit are:

        .. math::

            A = \begin{bmatrix}
            -\frac{R}{L} & -\frac{1}{L} \\
            \frac{1}{C} & 0
            \end{bmatrix},
            B = \begin{bmatrix}
            \frac{1}{L} \\
            0
            \end{bmatrix},
            C = \begin{bmatrix} 0 & 1 \end{bmatrix},
            D = \begin{bmatrix} 0 \end{bmatrix}


Thus, the state-space representation of the series RLC circuit is:

        .. math::

            \dot{X}(t) = \begin{bmatrix}
            -\frac{R}{L} & -\frac{1}{L} \\
            \frac{1}{C} & 0
            \end{bmatrix}
            \begin{bmatrix} x_1(t) \\ x_2(t) \end{bmatrix}
            + \begin{bmatrix}
            \frac{1}{L} \\
            0
            \end{bmatrix} V_{in}(t)

            Y(t) = \begin{bmatrix} 0 & 1 \end{bmatrix}
            \begin{bmatrix} x_1(t) \\ x_2(t) \end{bmatrix}
            + \begin{bmatrix} 0 \end{bmatrix} V_{in}(t)


The state-space representation provides a compact way of modeling
the series RLC circuit by using matrices to describe the system's
dynamics. The matrices :math:`A`, :math:`B`, :math:`C`, and :math:`D`
capture the relationships between the circuit's state variables,
input, and output. This representation is particularly useful for
analyzing the system's behavior in the time domain and for designing
control systems.

Solution

    >>> from sympy import Matrix, symbols, pprint
    >>> from sympy.physics.control import *
    >>> R, L, C = symbols('R L C')
    >>> A = Matrix([[-R/L, -1/L], [1/C, 0]])
    >>> B = Matrix([[1/L], [0]])
    >>> C = Matrix([[0, 1]])
    >>> D = Matrix([[0]])
    >>> ss = StateSpace(A, B, C, D)
    >>> ss
    StateSpace(Matrix([
    [-R/L, -1/L],
    [ 1/C,    0]]), Matrix([
    [1/L],
    [  0]]), Matrix([[0, 1]]), Matrix([[0]]))

    We can convert the StateSpace to TransferFunction by rewrite method.

    >>> tf = ss.rewrite(TransferFunction)[0][0]
    >>> tf
    TransferFunction(1, C*L*s**2 + C*R*s + 1, s)


Example 2
---------

        .. image:: Electrical_Problems_Q2.svg
           :align: center

Obtain the state model for a system represented by an electrical
system as shown in figure

The system is modeled with two state variables,
`x_1(t)` and `x_2(t)`, which are related to the physical voltages at the nodes
`v_1(t)` and `v_2(t)` respectively.

Let the two state variables be defined as:

        .. math::

           v_1(t) = x_1(t)

           v_2(t) = x_2(t)

The governing equations are derived by applying Kirchhoff's Current Law (KCL) at the nodes `v_1(t)` and `v_2(t)`.

Applying KCL at node `v_1(t)`:

        .. math::

           \frac{v_1(t) - u(t)}{R} + C \frac{d v_1(t)}{dt} + \frac{v_1(t) - v_2(t)}{R} = 0

Substituting the state variables:

        .. math::

           \frac{x_1(t) - u(t)}{R} + C \frac{dx_1(t)}{dt} + \frac{x_1(t) - x_2(t)}{R} = 0

Simplifying:

        .. math::

           C \dot{x_1}(t) = -\frac{2x_1(t)}{R} + \frac{x_2(t)}{R} + \frac{u(t)}{R}

Thus, the state equation for `x_1(t)` becomes:

        .. math::

           \dot{x_1}(t) = -\frac{2x_1(t)}{RC} + \frac{x_2(t)}{RC} + \frac{u(t)}{RC}


Applying KCL at node `v_2(t)`:

        .. math::

            C \frac{d v_2(t)}{dt} + \frac{v_2(t) - v_1(t)}{R} = 0

Substituting the state variables:

        .. math::

           C \frac{d x_2(t)}{dt} + \frac{x_2(t) - x_1(t)}{R} = 0

Simplifying:

        .. math::

           C \dot{x_2}(t) = \frac{x_1(t)}{R} - \frac{x_2(t)}{R}

Thus, the state equation for `x_2(t)` becomes:

        .. math::

           \dot{x_2}(t) = \frac{x_1(t)}{RC} - \frac{x_2(t)}{RC}

The state-space representation is given by the following matrix equation:

        .. math::

           \begin{bmatrix}
           \dot{x_1}(t) \\
           \dot{x_2}(t)
           \end{bmatrix}
           =
           \begin{bmatrix}
           -\frac{2}{RC} & \frac{1}{RC} \\
           \frac{1}{RC} & -\frac{1}{RC}
           \end{bmatrix}
           \begin{bmatrix}
           x_1(t) \\
           x_2(t)
           \end{bmatrix}
           +
           \begin{bmatrix}
           \frac{1}{RC} \\
           0
           \end{bmatrix}
           u(t)

The output of the circuit is defined as:

        .. math::

           y(t) = v_2(t) = x_2(t)

Thus, the output equation can be written as:

        .. math::

           y(t) = \begin{bmatrix} 0 & 1 \end{bmatrix}
           \begin{bmatrix}
           x_1(t) \\
           x_2(t)
           \end{bmatrix}


Solution

    >>> from sympy import symbols, Matrix
    >>> from sympy.physics.control import *
    >>> R, C = symbols('R C')
    >>> A = Matrix([[-2/(R*C), 1/(R*C)], [1/(R*C), -1/(R*C)]])
    >>> B = Matrix([[1/(R*C)], [0]])
    >>> C = Matrix([[0, 1]])
    >>> ss = StateSpace(A, B, C)
    >>> ss
    StateSpace(Matrix([
    [-2/(C*R),  1/(C*R)],
    [ 1/(C*R), -1/(C*R)]]), Matrix([
    [1/(C*R)],
    [      0]]), Matrix([[0, 1]]), Matrix([[0]]))


References
----------
1. `bmsce.ac.in <https://bmsce.ac.in/Content/TE/STATE_SPACE_ANALYSIS.pdf>`_
