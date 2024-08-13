.. _rotating_disc_block_friction:

==========================
A Block on a Rotating Disc
==========================

This page demonstrates how to use the functionalities in :obj:`sympy.physics.mechanics`
and the **CoulombKineticFriction** actuator to model a sliding block on a rotating disc.
The block is located at a point on the disc surface and is subjected to forces, including
friction force with additional Stribeck and viscous effects.

.. raw:: html
   :file: rotating_disc_block.svg

This system will be modeled using Kane's method. Here, :math:`F_f` is the friction force.
We assume that the block is initially sliding (relative motion is non-zero), moving
in the direction opposite to the disc's rotation due to the relative motion, and the
kinetic friction acts in the direction opposite to this relative motion.
The disc has an infinite radius, preventing the block from flying off the disc at any
point in this model. The normal force is a positive scalar, which allows the valid usage
of the **CoulombKineticFriction** actuator.

   >>> import sympy as sm
   >>> import sympy.physics.mechanics as me
   >>> me.init_vprinting()

Let's define the necessary variables and coordinates.

   >>> m, g, mu_k, mu_s, v_s, sigma = sm.symbols('m, g, mu_k, mu_s, v_s, sigma')
   >>> q1, q2, q3 = me.dynamicsymbols('q1:4')
   >>> q1d, q2d, q3d = me.dynamicsymbols('q1:4', level=1)
   >>> u1, u2, u3 = me.dynamicsymbols('u1:4')
   >>> u1d, u2d, u3d = me.dynamicsymbols('u1:4', level=1)

- :math:`q1`: Generalized coordinate representing the angular velocity of the disc w.r.t. the inertial frame
- :math:`q2, q3`: Generalized coordinates of the block w.r.t the disc
- :math:`mu_k`: coefficient of kinetic friction
- :math:`mu_s`: coefficient of static friction
- :math:`v_s`: Stribeck friction coefficient
- :math:`sigma`: viscous friction coefficient
- :math:`m`: mass of the block
- :math:`g`: gravitational constant

We will also define the necessary reference frames, points, and velocities.
:math:`N` is the inertial reference frame and :math:`A` is the rotating reference frame.
:math:`O` is the center of the disc, :math:`P` is the actual contact point between
the block and the disc.

   >>> N = me.ReferenceFrame('N')
   >>> A = me.ReferenceFrame('A')
   >>> A.set_ang_vel(N, u1 * N.x)
   >>> A.ang_vel_in(N)
   u1 n_x

   >>> O = me.Point('O')
   >>> P = O.locatenew('P', q2 * N.x + q3 * N.y)
   >>> O.set_vel(N, 0)
   >>> P.v2pt_theory(P, N, A)
   q2'(t) n_x + q3'(t) n_y

We define the particle, pathway, and forces using **CoulombKineticFriction**.
A unique custom pathway, SlidingPathway, represent the motion of the block along a plane.

   >>> class SlidingPathway(me.PathwayBase):
   ...
   ...     def __init__(self, fixed_point, contact_point, frame):
   ...         """A custom pathway that moves along a rotating disc.
   ...
   ...         Two points, ``fixed_point`` and ``contact_point``, start from the same
   ...         location on the disc. ``contact_point``, the actual sliding point of the block,
   ...         moves along the ``SlidingPathway`` and is subjected to the friction force.
   ...         ``SlidingPathway`` is valid only as long as the two points lie within the
   ...         same plane, and it ensures the direction of the force is accurate.
   ...
   ...         Parameters
   ...         ==========
   ...         fixed_point : Point
   ...             Point fixed w.r.t the provided reference frame; the point on the disc.
   ...         contact_point : Point
   ...             Point with a time-dependent position w.r.t the fixed point, where the
   ...             actuator forces will be applied. This point is assumed to be fixed to
   ...             another body; the point on the block.
   ...         frame : ReferenceFrame
   ...             Reference frame in which the ``fixed_point`` has a zero velocity.
   ...
   ...         """
   ...
   ...         self.fixed_point = fixed_point
   ...         self.contact_point = contact_point
   ...         self.frame = frame
   ...
   ...     @property
   ...     def length(self):
   ...         return self.contact_point.pos_from(self.fixed_point).magnitude()
   ...
   ...     @property
   ...     def extension_velocity(self):
   ...         """Extension velocity of the pathway.
   ...
   ...         The extension velocity of the pathway is the magnitude of the velocity of
   ...         the ``contact_point`` relative to the frame in which the ``fixed_point``
   ...         is stationary.
   ...
   ...         """
   ...
   ...         return self.contact_point.vel(self.frame).magnitude()
   ...
   ...     def to_loads(self, force):
   ...         """Loads in the correct format to be supplied to `KanesMethod`.
   ...
   ...         Forces applied to the ``contact_point`` and ``fixed_point``
   ...         based on the friction force.
   ...
   ...         """
   ...
   ...         direction = -self.contact_point.vel(self.frame).normalize()
   ...         force = mu_k * m * g
   ...
   ...         return [
   ...             me.Force(self.fixed_point, -force * direction),
   ...             me.Force(self.contact_point, force * direction),
   ...             ]

   >>> block = me.Particle('block', P, m)
   >>> normal_force = m * g
   >>> pathway = SlidingPathway(O, P, N)
   >>> friction = me.CoulombKineticFriction(mu_k, normal_force, pathway, v_s=v_s, sigma=sigma, mu_s=mu_k)
   >>> loads = friction.to_loads()
   >>> loads
           g*m*mu_k*q2'(t)              g*m*mu_k*q3'(t)                 -g*m*mu_k*q2'(t)             -g*m*mu_k*q3'(t)
    [(O, ---------------------- n_x + ---------------------- n_y), (P, ---------------------- n_x + ---------------------- n_y)]
            ___________________          ___________________              ___________________          ___________________
           /       2         2          /       2         2              /       2         2          /       2         2
         \/  q2'(t)  + q3'(t)         \/  q2'(t)  + q3'(t)             \/  q2'(t)  + q3'(t)         \/  q2'(t)  + q3'(t)

Now, we're ready to use Kane's method to obtain the equations of motion.

   >>> BL = [block]
   >>> kane = me.KanesMethod(
   ...     N,
   ...     q_ind=[q2, q3],
   ...     u_ind=[u2, u3],
   ...     kd_eqs=[q2d - u2, q3d - u3],
   ...     bodies=BL
   ...     )

   >>> fr, frstar = kane.kanes_equations(BL, loads)
   >>> eom = fr + frstar
   >>> eom
   [   g*m*mu_k*u2             ]
   [- -------------- - m*u2'(t)]
   [     ___________           ]
   [    /   2     2            ]
   [  \/  u2  + u3             ]
   [                           ]
   [   g*m*mu_k*u3             ]
   [- -------------- - m*u3'(t)]
   [     ___________           ]
   [    /   2     2            ]
   [  \/  u2  + u3             ]
