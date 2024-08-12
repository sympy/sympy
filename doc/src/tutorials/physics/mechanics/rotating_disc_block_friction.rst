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

This system will be modeled using Kane's method. Here, :math:`m` is the mass of the block and the disc,
:math:`r` represents the location of the block from the center of the disc, :math:`g` is the
gravitational constant, :math:`F_f` is the friction force, :math:`F_{cp}` is the centripetal force,
and :math:`F_{cf}` is the centrifugal force.

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

   >>> m, g, r, mu_k, mu_s, v_s, sigma = sm.symbols('m, g, r, mu_k, mu_s, v_s, sigma')
   >>> q1, q2, q3 = me.dynamicsymbols('q1 q2 q3')
   >>> u1, u2, u3 = me.dynamicsymbols('q1 q2 u3', 1)
   >>> t = me.dynamicsymbols._t

- :math:`q1`: Generalized coordinate representing the angular velocity of the disc w.r.t. the inertial frame
- :math:`q2, q3`: Generalized coordinates of the block w.r.t the disc
- :math:`mu_k`: coefficient of kinetic friction
- :math:`mu_s`: coefficient of static friction
- :math:`v_s`: Stribeck friction coefficient
- :math:`sigma`: viscous friction coefficient

We will also define the necessary reference frames, points, and velocities.
:math:`N` is the inertial reference frame, :math:`A1` is the rotating reference frame,
and :math:`A2` is the frame fixed to the block.

:math:`O` is the center of the disc, :math:`P` is a fixed point on the disc at a
distance :math:`r` from the center and in contact with the block, and :math:`Q` is
the actual contact point between the block and the disc.

   >>> N, A1, A2 = sm.symbols('N, A1, A2', cls=me.ReferenceFrame)
   >>> O, P, Q = sm.symbols('O, P, Q', cls=me.Point)
   >>> O.set_vel(N, 0)

   >>> A1.orient_axis(N, q1, N.z)
   >>> A1.set_ang_vel(N, u1 * N.z)

   >>> A2.set_ang_vel(A1, u2 * A1.x + u3 * A2.y)

   >>> P.set_pos(O, q1 * N.z)
   >>> P.set_vel(A2, u1 * N.z)

   >>> Q.set_pos(O, q2 * A1.x + q3 * A1.y)
   >>> Q.set_vel(A1, u2 * A1.x + u3 * A1.y)

We define the particle, pathway, and forces using **CoulombKineticFriction**.
The pathway should be able to represent the motion of the block along a plane, so we use
a unique custom pathway, SlidingPathway.

   >>> class SlidingPathway(me.PathwayBase):
   ...
   ...     def __init__(self, fixed_point, contact_point, frame, m):
   ...         """A custom pathway that moves along a rotating disc.
   ...
   ...         We assume that two points, which represent a block and a disc respectively,
   ...         lie in the same plane, and ``SlidingPathway`` ensures that the forces and velocities
   ...         are correctly related in this plane.
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
   ...             Reference frame in which the ``fixed_point`` has a zero velocity. The pathway
   ...             is only valid as long as two points lie within the same plane.
   ...         m : sympifiable
   ...             Mass of the block.
   ...
   ...         """
   ...
   ...         self.fixed_point = fixed_point
   ...         self.contact_point = contact_point
   ...         self.frame = frame
   ...         self.m = m
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
   ...         Forces and torques applied to the ``contact_point`` and ``fixed_point``
   ...         based on the friction force.
   ...
   ...         """
   ...
   ...         direction = -self.contact_point.vel(self.frame).normalize()
   ...         acc = self.contact_point.vel(A1).diff(t, A1).magnitude()
   ...         force = acc * self.m
   ...
   ...         return [
   ...             me.Force(self.contact_point, force * direction),
   ...             me.Force(self.fixed_point, -force * direction),
   ...             me.Torque(self.frame, me.cross(self.contact_point.pos_from(self.fixed_point), -force * direction)),
   ...         ]

   >>> block = me.Particle('block', Q, m)
   >>> disc = me.Particle('disc', P, m)
   >>> normal_force = m * g
   >>> pathway = SlidingPathway(P, Q, N, m)
   >>> friction = me.CoulombKineticFriction(mu_k, normal_force, pathway, v_s=v_s, sigma=sigma, mu_s=mu_k)
   >>> loads = friction.to_loads()

Now, we're ready to use Kane's method to obtain the equations of motion.

   >>> kane = me.KanesMethod(
   ...     N,
   ...     q_ind=[q2, q3],
   ...     u_ind=[u2, u3],
   ...     kd_eqs=[q2.diff() - u2, q3.diff() - u3],
   ...     bodies=[block]
   ...     )

   >>> fr, frstar = kane.kanes_equations([block], loads)
   >>> eom = fr + frstar
