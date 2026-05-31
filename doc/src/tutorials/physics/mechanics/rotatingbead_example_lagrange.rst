============================================
Bead Motion on a Rotating Parabolic Path
============================================

This example analyzes the dynamics of a particle constrained to move along a rotating
parabolic path using Lagrangian mechanics.

System Configuration
--------------------
A particle of mass *M* slides frictionlessly along a parabolic guide defined by
:math:`z = \kappa \rho^2`, where :math:`\rho` is the horizontal displacement from
the rotation axis and :math:`\kappa` is the curvature parameter. The guide rotates
with constant angular rate :math:`\varpi` about the vertical axis. We determine the
required curvature :math:`\kappa` for sustained circular motion at fixed radius
:math:`\lambda`.

Kinematic Relationships
-----------------------
The position vector is expressed in cylindrical coordinates with angular orientation
:math:`\phi = \varpi t` due to the guide's rotation::

  >>> from sympy import symbols, cos, sin
  >>> from sympy.physics.mechanics import *
  >>> mechanics_printing(pretty_print=False)

  # System parameters and coordinates
  >>> M, grav, varpi, kappa, lam = symbols('M grav varpi kappa λ')
  >>> time = dynamicsymbols._t
  >>> rho = dynamicsymbols('ρ')  # Radial displacement
  
  # Define inertial reference frame
  >>> Frame = ReferenceFrame('Frame')
  
  # Position vector components
  >>> phi = varpi * time
  >>> pos_vec = rho * (cos(phi)*Frame.x + sin(phi)*Frame.y) + kappa*rho**2 * Frame.z

Velocity computation in the inertial frame::

  >>> vel_vec = pos_vec.dt(Frame)
  >>> mprint(vel_vec.simplify())
  (-varpi*ρ*sin(varpi*t) + cos(varpi*t)*ρ')*Frame.x + 
  (varpi*ρ*cos(varpi*t) + sin(varpi*t)*ρ')*Frame.y + 
  2*kappa*ρ*ρ'*Frame.z

Energy Formulation
------------------
The kinetic energy (T) and gravitational potential energy (V)::

  >>> T_kin = M/2 * vel_vec.dot(vel_vec)
  >>> V_pot = M * grav * kappa * rho**2
  >>> Lagrangian = T_kin - V_pot

Equations of Motion
-------------------
Applying Lagrange's method with generalized coordinate :math:`\rho`::

  >>> LM = LagrangesMethod(Lagrangian, [rho])
  >>> LM.form_lagranges_equations().simplify()
  
  # Display simplified equations
  >>> LM.eom
  Matrix([[M*(2*grav*kappa*ρ + 4*kappa**2*ρ**2*ρ'' + 4*kappa**2*ρ*ρ'**2 - varpi**2*ρ + ρ'')]])

Steady Circular Motion Condition
--------------------------------
For sustained circular motion at radius :math:`\lambda`, substitute:
:math:`\rho = \lambda`, :math:`\dot{\rho} = 0`, :math:`\ddot{\rho} = 0`::

  >>> steady_cond = {rho: lam, rho.diff(time): 0, rho.diff(time, 2): 0}
  >>> balance_eq = LM.eom.subs(steady_cond)
  
  # Solve for curvature parameter
  >>> from sympy import solve
  >>> kappa_sol = solve(balance_eq, kappa)
  >>> mprint(kappa_sol[kappa])
  varpi**2/(2*grav)

Fundamental Relationship
------------------------
The curvature parameter for stable circular motion is:

.. math::
  \kappa = \frac{\varpi^2}{2g}

This result is confirmed by substituting back into the equations of motion::

  >>> balance_eq.subs(kappa, kappa_sol[kappa])
  Matrix([[0]])

The analysis demonstrates how the required path curvature scales quadratically with
the rotation rate to maintain equilibrium motion at fixed radius.

For more details on the problem, check page 245 on [Thornton2003]_.

References
==========

.. [Thornton2003] Thornton, S., & Marion, J. (2003). Classical Dynamics of Particles and Systems (5th edition). Cengage Learning.
