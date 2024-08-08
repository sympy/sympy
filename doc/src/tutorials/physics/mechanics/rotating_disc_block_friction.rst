.. _rotating_disc_block_friction:

==========================
A Block on a Rotating Disc
==========================

This page demonstrates how to use the functionalities in :obj:`sympy.physics.mechanics`
and the **CoulombKineticFriction** actuator to model a sliding block on a rotating disc.
The particle is located at a point on the disc surface and is subjected to forces, including
friction force with additional Stribeck and viscous effects.

.. raw:: html
   :file: rotating_disc_block.svg

This system will be modeled using Kane's method. Here, :math:`m` is the mass of the block,
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
   >>> q1, q2 = me.dynamicsymbols('q1 q2')
   >>> u1, u2 = me.dynamicsymbols('q1 q2', 1)

- :math:`q1, q2`: generalized coordinates of the block w.r.t the disc [rad]
- :math:`u1, u2`: generalized angual velocities of the block [rad/sec]
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
   >>> P = O.locatenew('P', r * sm.cos(q1) * A.x + r * sm.sin(q1) * A.y)
   >>> Q = P.locatenew('Q', 0)

The point :math:`P` moves with the disk, having a tangential velocity due to the
disc's rotation, and the initial velocity of the block at point :math:`Q` includes
the disc's tangential velocity and the sliding velocity :math:`v0` relative to the disc.

   >>> O.set_vel(N, 0)
   >>> P.set_vel(N, r * u1 * A.y)
   >>> Q.set_vel(N, P.vel(N) + v0 * A.y)

The position, velocity, and acceleration of point :math:`P` will give the inertia force,
which includes the Coriolis force, centrifugal force, and centripetal force.

   >>> Q_pos = Q.pos_from(P)
   >>> Q_vel = Q.vel(N)
   >>> Q_acc = Q.acc(N)
   >>> Q_acc
   -(r*q1'(t) + v0)*q1'(t) a_x + r*q1''(t) a_y

We define the particle, pathway, and forces using **CoulombKineticFriction**.
The pathway should be able to represent the motion of the block along a plane, so we use
a unique custom pathway.

   >>> class CustomPathway(me.PathwayBase):
   ...     
   ...     def __init__(self,)
   ...         self.
   ...     
   ...     @property
   ...     def extension_velocity(self):
   ...         return

   >>> block = me.Particle('block', Q, m)
   >>> normal_force = m * g
   >>> pathway = CustomPathway()
   >>> friction = me.CoulombKineticFriction(mu_k, normal_force, pathway, v_s=v_s, sigma=sigma, mu_s=mu_k)

Now, we're ready to use Kane's method to obtain the equations of motion.

   >>> kanes_method = me.KanesMethod(
   ...     N,
   ...     q_ind=[q1, q2],
   ...     u_ind=[u1, u2],
   ...     kd_eqs=[q1.diff() - u1, q2.diff() - u2],
   ...     bodies=[block]
   ...     )
   >>> bodies = [block]

   >>> loads = friction.to_loads()
   >>> fr, frstar = kane.kanes_equations([block], loads)
   >>> eom = fr + frstar
