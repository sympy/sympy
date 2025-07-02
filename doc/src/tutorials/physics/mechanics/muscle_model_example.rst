.. _muscle_model_example:

=====================================
Muscle Wrapping Over a Spherical Bone
=====================================

.. _fig-muscle-model:
.. figure:: muscle_model.svg

The wrapping sphere muscle model consists of a muscle that runs from an origin point $A$ to an insertion point $B$, wrapping smoothly over a fixed spherical bone of radius $R$ centered at the origin. The muscle contacts the sphere at a moving point $P(t)$ that varies with time based on the muscle's activation and the mechanical constraints of the system.

The complete muscle path is composed of three segments:

1. **Straight line segment** from origin $A$ to tangent point $T_A$
2. **Geodesic segment** on the sphere surface from $T_A$ to $T_B$  
3. **Straight line segment** from tangent point $T_B$ to insertion $B$

Here, $T_A$ and $T_B$ are the tangent points on the sphere where the straight lines from $A$ and $B$ touch the spherical surface. This model is particularly useful in biomechanics for analyzing how muscles wrap around bones and joints during movement.

:obj:`sympy.physics.mechanics` provides wrapping geometry classes including ``WrappingSphere`` that are specifically designed to handle the complex geodesic calculations required for this type of analysis. This example demonstrates how to set up the complete kinematic and dynamic model, compute the total muscle path length, and derive the equations of motion using Lagrangian mechanics.

Define Symbols and Import Modules
=================================

First, we import the necessary symbols, frames, points, and classes from SymPy and its mechanics module. We'll use spherical coordinates $\theta(t)$ and $\phi(t)$ to describe the contact point on the sphere surface.

    >>> from sympy import (
    ...      symbols, sin, cos, sqrt,
    ...      simplify, init_printing,
    ...      pprint, Rational
    ... )
    >>> from sympy.physics.mechanics import (
    ...      dynamicsymbols, ReferenceFrame,
    ...      Point, dot, LagrangesMethod,
    ...      Lagrangian, WrappingSphere,
    ...      Particle
    ... )

    >>> init_printing(use_unicode=True)

The ``init_printing`` function enables pretty mathematical output formatting, making the equations more readable. The ``WrappingSphere`` class will handle the complex geodesic length calculations on the spherical surface.

Generalized Coordinates and Physical Parameters
===============================================

We define the time variable and establish our generalized coordinates. The spherical coordinates $\theta(t)$ and $\phi(t)$ completely describe the position of the contact point $P$ on the sphere surface, where $\theta$ is the polar angle (measured from the positive z-axis) and $\phi$ is the azimuthal angle (measured from the positive x-axis in the xy-plane).

    >>> # Time variable
    >>> t = symbols('t')

    >>> # Generalized coordinates and their time derivatives
    >>> theta, phi = dynamicsymbols('theta phi')
    >>> dtheta, dphi = dynamicsymbols('theta phi', 1)

    >>> # Physical parameters
    >>> R, m, k, L0 = symbols('R m k L0', positive=True)
    >>> Ax, Ay, Az, Bx, By, Bz = symbols('Ax Ay Az Bx By Bz')

The physical parameters represent:
- $R$: radius of the spherical bone
- $m$: effective mass associated with the contact point (representing muscle inertia)
- $k$: elastic stiffness constant of the muscle
- $L_0$: natural (rest) length of the muscle
- $(A_x, A_y, A_z)$ and $(B_x, B_y, B_z)$: fixed spatial coordinates of the muscle's origin and insertion points

Inertial Frame and WrappingSphere Setup
=======================================

We establish an inertial reference frame $N$ and place the sphere's center $O$ at the origin with zero velocity. The ``WrappingSphere`` object will automatically handle the geodesic computations required for calculating the shortest path along the spherical surface.

    >>> # Inertial frame and fixed origin
    >>> N = ReferenceFrame('N')
    >>> O = Point('O')
    >>> O.set_vel(N, 0)

    >>> # Define wrapping sphere centered at origin
    >>> sphere = WrappingSphere(R, O)

The sphere is assumed to be fixed in space (representing a bone), so its center has zero velocity in the inertial frame.

Position and Velocity of Contact Point
======================================

The contact point $P$ on the sphere surface is parameterized using spherical coordinates. In the standard spherical coordinate system, the position vector is:

$$\mathbf{r}_P = R\sin\theta\cos\phi\,\hat{\mathbf{N}}_x + R\sin\theta\sin\phi\,\hat{\mathbf{N}}_y + R\cos\theta\,\hat{\mathbf{N}}_z$$

The velocity of point $P$ is obtained by differentiating its position vector with respect to time in the inertial frame $N$.

    >>> # Contact point P in spherical coordinates
    >>> P = Point('P')
    >>> x = R*sin(theta)*cos(phi)
    >>> y = R*sin(theta)*sin(phi)
    >>> z = R*cos(theta)

    >>> P.set_pos(O, x*N.x + y*N.y + z*N.z)
    >>> P.set_vel(N, P.pos_from(O).diff(t, N))

The resulting velocity will be a function of $\dot{\theta}$ and $\dot{\phi}$, representing how the contact point moves along the sphere surface.

Kinetic Energy Calculation
==========================

We model the contact point $P$ as a particle with mass $m$ to account for the inertial effects of the muscle tissue. The kinetic energy of this particle is computed using the standard formula $T = \frac{1}{2}mv^2$.

    >>> # Model P as a particle with mass m
    >>> P_part = Particle('P_part', P, m)

    >>> # Kinetic energy T = (1/2)mv²
    >>> T = simplify(P_part.kinetic_energy(N))
    >>> pprint(T)

The kinetic energy will be expressed in terms of $\dot{\theta}$ and $\dot{\phi}$, showing how the energy depends on the rates of change of the spherical coordinates.

Fixed Points and Tangent Point Determination
============================================

The muscle's origin point $A$ and insertion point $B$ are fixed in space. To find where the muscle would naturally contact the sphere (the tangent points), we project rays from the sphere center $O$ toward points $A$ and $B$ onto the sphere surface.

    >>> # Fixed muscle origin and insertion points
    >>> A = Point('A'); B = Point('B')
    >>> A.set_pos(O, Ax*N.x + Ay*N.y + Az*N.z)
    >>> B.set_pos(O, Bx*N.x + By*N.y + Bz*N.z)

    >>> # Unit direction vectors from sphere center
    >>> uA = A.pos_from(O).normalize()
    >>> uB = B.pos_from(O).normalize()

    >>> # Tangent points on sphere surface
    >>> TA = Point('TA'); TB = Point('TB')
    >>> TA.set_pos(O, R*uA)
    >>> TB.set_pos(O, R*uB)

The tangent points $T_A$ and $T_B$ represent the "natural" contact points where the muscle would touch the sphere if it followed the shortest possible path from $A$ to $B$ while wrapping around the spherical obstacle.

Total Muscle Path Length Computation
====================================

The total muscle path consists of three segments whose lengths must be computed and summed:

1. **Straight segment** from $A$ to the current contact point $P$
2. **Geodesic segment** on the sphere from $T_A$ to $T_B$ (shortest path on sphere surface)
3. **Straight segment** from $P$ to insertion point $B$

The ``WrappingSphere`` class provides the ``geodesic_length`` method to compute the arc length of the shortest path between two points on the sphere surface.

    >>> # Geodesic segment length on sphere surface
    >>> L_wrap = simplify(sphere.geodesic_length(TA, TB))

    >>> # Straight-line segment vectors and lengths
    >>> TA_vec = P.pos_from(A)
    >>> TB_vec = P.pos_from(B)

    >>> # Total muscle path length
    >>> L_tot = simplify(
    ...     sqrt(dot(TA_vec, TA_vec)) +
    ...     L_wrap +
    ...     sqrt(dot(TB_vec, TB_vec))
    ... )
    >>> pprint(L_tot)

The total length $L_{total}$ will be a function of the generalized coordinates $\theta$ and $\phi$, representing how the muscle length changes as the contact point moves around the sphere.

Potential Energy and Lagrangian Formulation
===========================================

The muscle is modeled as an elastic element with spring constant $k$ and natural length $L_0$. The elastic potential energy is: $V = \frac{1}{2}k(L_{total} - L_0)^2$, where the energy increases quadratically with the deviation from the natural length.

    >>> # Elastic potential energy of the muscle
    >>> P_part.potential_energy = (
    ...     Rational(1, 2)*k*(L_tot - L0)**2
    ... )

    >>> # Lagrangian L = T - V (kinetic minus potential energy)
    >>> Lag = Lagrangian(N, P_part)
    >>> pprint(Lag)

The ``Lagrangian`` class automatically computes $L = T - V$ using the kinetic energy of the particle and the potential energy we've assigned to it. This Lagrangian captures the complete dynamics of the system.

Equations of Motion Derivation
==============================

Using Lagrange's method, we derive the equations of motion for the system. The ``LagrangesMethod`` class automatically applies the Euler-Lagrange equations:

$$\frac{d}{dt}\left(\frac{\partial L}{\partial \dot{q}_i}\right) - \frac{\partial L}{\partial q_i} = 0$$

for each generalized coordinate $q_i$ (in our case, $\theta$ and $\phi$).

    >>> # Form Lagrange's equations for θ and φ
    >>> LM = LagrangesMethod(Lag, [theta, phi])
    >>> eqns = LM.form_lagranges_equations()

    >>> # Display the resulting ordinary differential equations
    >>> for i, eq in enumerate(eqns, 1):
    ...     print(f"Equation {i}:")
    ...     pprint(simplify(eq))

The resulting equations will be coupled, nonlinear second-order ordinary differential equations in $\theta(t)$ and $\phi(t)$. These equations describe how the contact point moves on the sphere surface under the influence of the muscle's elastic forces and the geometric constraints imposed by the spherical wrapping.

Conclusion
==========

This tutorial has demonstrated how to construct a comprehensive dynamic model of muscle wrapping around a spherical bone using SymPy's mechanics framework. The ``WrappingSphere`` class automatically handles the complex geodesic calculations, while the Lagrangian formulation provides a systematic approach to deriving the equations of motion.

This framework can be extended to more complex scenarios involving multiple muscles, non-spherical bones, and time-varying activation patterns.
