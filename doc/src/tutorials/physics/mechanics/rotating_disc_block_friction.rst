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

   >>> m, g, r, mu_k, mu_s, v_s, sigma, v0 = sm.symbols('m, g, r, mu_k, mu_s, v_s, sigma, v0')
   >>> q1, q2, q3 = me.dynamicsymbols('q1 q2 q3')
   >>> u1, u2, u3 = me.dynamicsymbols('q1 q2 u3', 1)

- :math:`q1`: Generalized coordinate representing the angular velocity of the disc w.r.t. the inertial frame
- :math:`q2, q3`: Generalized coordinates of the block w.r.t the disc
- :math:`v0`: initial relative velocity of the block
- :math:`mu_k`: coefficient of kinetic friction
- :math:`mu_s`: coefficient of static friction
- :math:`v_s`: Stribeck friction coefficient
- :math:`sigma`: viscous friction coefficient

We will also define the inertial reference frame :math:`N` and the rotating reference frame :math:`A`.

   >>> N, A = sm.symbols('N, A', cls=me.ReferenceFrame)
   >>> A.set_ang_vel(N, u1 * A.z)

Next, we define the essential points and velocities.
:math:`O` is the center of the disc, :math:`P` is a fixed point on the disc at a
distance :math:`r` from the center and in contact with the block, and :math:`Q` is 
the actual contact point between the block and the disc.

   >>> O = me.Point('O')
   >>> P = O.locatenew('P', r * q1 * A.x)
   >>> Q = P.locatenew('Q', r * q2 * A.x + r * q3 * A.y)

The point :math:`P` moves with the disk, having a tangential velocity due to the
disc's rotation, and the point :math:`Q`has a velocity that includes the tangential
velocity of :math:`P` (from the disc's rotation) and the block's sliding velocity
relative to the disc.

   >>> O.set_vel(N, 0)
   >>> P.set_vel(N, r * u1 * A.y)
   >>> Q.set_vel(N, P.vel(N) + v0 * A.x + v0 * A.y)

We define the particle, pathway, and forces using **CoulombKineticFriction**.
The pathway should be able to represent the motion of the block along a plane, so we use
a unique custom pathway, SlidingPathway.

   >>> class SlidingPathway(me.PathwayBase):
   ...
   ...     def __init__(self, fixed_point, contact_point, frame):
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
   ...         Forces and torques applied to the ``contact_point`` and ``fixed_point``
   ...         based on the friction force.
   ...
   ...         """
   ...
   ...         F = me.dynamicsymbols('F')
   ...         direction = -self.contact_point.vel(self.frame).normalize()
   ...         force = F
   ...         return [
   ...             me.Force(self.contact_point, force * direction),
   ...             me.Force(self.fixed_point, -force * direction),
   ...             me.Torque(self.frame, me.cross(self.contact_point.pos_from(self.fixed_point), -force*direction)),
   ...         ]

   >>> block = me.Particle('block', Q, m)
   >>> disc = me.Particle('disc', P, m)
   >>> normal_force = m * g
   >>> pathway = SlidingPathway(P, Q, N)
   >>> friction = me.CoulombKineticFriction(mu_k, normal_force, pathway, v_s=v_s, sigma=sigma, mu_s=mu_k)
   >>> loads = friction.to_loads() + pathway.to_loads(force=None)

Now, we're ready to use Kane's method to obtain the equations of motion.

   >>> kane = me.KanesMethod(
   ...     N,
   ...     q_ind=[q1, q2, q3],
   ...     u_ind=[u1, u2, u3],
   ...     kd_eqs=[q1.diff() - u1, q2.diff() - u2, q3.diff() - u3],
   ...     bodies=[block, disc]
   ...     )

   >>> fr, frstar = kane.kanes_equations([block, disc], loads)
   >>> eom = fr + frstar
