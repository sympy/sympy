.. _rotating_disc_block_friction:

==========================
A Block on a Rotating Disc
==========================

This page demonstrates how to use the functionalities in :obj:`sympy.physics.mechanics`
and the :class:`~.CoulombKineticFriction` actuator to model a sliding block on a rotating disc.
This is a nonholonomic system affected by friction, including Coulomb friction with the
Stribeck effect, viscous friction, and other constraints.

.. raw:: html
   :file: block_on_rotating_disc.svg

This system will be modeled using Lagrange's method. We treat the block as a point mass that
slides under friction on top of a rotating disc. The disc rotates with constant angular speed, and
the block is allowed to move radially and tangentially relative to the disc.

Friction forces are modeled using the :class:`~.CoulombKineticFriction` actuator, which supports regularized
static and kinetic Coulomb friction along with viscous damping.

   >>> import sympy as sm
   >>> import sympy.physics.mechanics as me
   >>> me.init_vprinting()

Let's define the necessary variables and coordinates.

   >>> t = sm.Symbol('t')
   >>> r, theta = me.dynamicsymbols('r, theta') # generalized coordinates
   >>> r_dot, theta_dot = me.dynamicsymbols('r, theta', 1)

   >>> m, omega = sm.symbols('m omega')
   >>> mu_kr, mu_sr, vr_s, br = sm.symbols('mu_kr mu_sr v_sr b_r')
   >>> mu_kt, mu_st, vt_s, bt = sm.symbols('mu_kt mu_st v_st b_t')
   >>> F_nr, F_nt = sm.symbols('F_nr F_nt', positive=True)

- :math:`r(t)`: Generalized coordinate representing the radial position of the block from the disc center
- :math:`\\theta(t)`: Generalized coordinate representing the angular position of the block in the disc’s rotating frame
- :math:`m`: Mass of the block
- :math:`\\omega`: Angular velocity of the rotating disc (assumed constant)
- :math:`\\mu_{kr}, \\mu_{kt}`: Coefficients of kinetic friction (radial / tangential)
- :math:`\\mu_{sr}, \\mu_{st}`: Coefficients of static friction (radial / tangential)
- :math:`v_{sr}, v_{st}`: Velocity thresholds for Stribeck effect smoothing (radial / tangential)
- :math:`b_r, b_t`: Smoothing steepness parameters in exponential regularization (radial / tangential)
- :math:`F_{nr}, F_{nt}`: Positive normal forces associated with radial / tangential friction

Reference frames, points, and velocities are defined next.
:math:`N` is the inertial reference frame and :math:`A` is the rotating reference frame.
:math:`O` is the center of the disc, :math:`P` is the actual contact point between
the block and the disc, and :math:`Q` is a nearby point on the disc, offset by a small angle
:math:`dtheta` from :math:`P`, at the same radius :math:`r`.

   >>> N = me.ReferenceFrame('N')
   >>> O = me.Point('O')
   >>> O.set_vel(N, 0)
   
   >>> A = N.orientnew('A', 'Axis', [omega * t, N.z])
   >>> A.set_ang_vel(N, omega * N.z)

   >>> P = O.locatenew('P', r * sm.cos(theta) * A.x + r * sm.sin(theta) * A.y)
   >>> P.set_vel(N, P.pos_from(O).dt(N)) # velocity in inertial frame
   >>> P.set_vel(A, P.vel(N).express(A))

We define the block as a particle with mass :math:`m` located at :math:`P`.

   >>> block = me.Particle('Block', P, m)

We now define a custom ``SlidingPathway`` that describes how the block slides on the disc surface.

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
   ...         direction = self.contact_point.vel(self.frame).normalize()
   ...         force = self.extension_velocity * force
   ...
   ...         return [
   ...             me.Force(self.fixed_point, -force * direction),
   ...             me.Force(self.contact_point, force * direction),
   ...             ]

The ``SlidingPathway`` allows friction to act along the relative motion direction between
two points in a given reference frame. In our model, we define two such directions:

- **Radial**: along the line from center :math:`O` to block :math:`P`
- **Tangential**: perpendicular to the radial, around the disc's edge

We now define both friction forces.

**Radial Friction**:

   >>> radial_pathway = SlidingPathway(O, P, A)
   >>> radial_friction = me.CoulombKineticFriction(mu_kr, F_nr, radial_pathway,
   ...                                             mu_s=mu_sr, v_s=vr_s, sigma=br)

**Tangential Friction**: To define a direction of motion in the tangential direction,
we construct a nearby point :math:`Q` slightly ahead of :math:`P` in the angular direction.

   >>> dtheta = sm.Symbol('dtheta') # small angle increment
   >>> Q = O.locatenew('Q', r * sm.cos(theta + dtheta) * A.x + r * sm.sin(theta + dtheta) * A.y)
   >>> Q.set_vel(N, Q.pos_from(O).dt(N))
   >>> Q.set_vel(A, Q.vel(N).express(A))

   >>> tangent_pathway = SlidingPathway(Q, P, A)
   >>> tangent_friction = me.CoulombKineticFriction(mu_kt, F_nt, tangent_pathway,
   ...                                              mu_s=mu_st, v_s=vt_s, sigma=bt)

We now collect the forces applied by both friction actuators.

   >>> loads = radial_friction.to_loads() + tangent_friction.to_loads()

Now, we use Lagrange's method to obtain the equations of motion.

   >>> coordinates = [r, theta]
   >>> speeds = [r.diff(t), theta.diff(t)]

   >>> L = me.Lagrangian(N, block)
   >>> L
   ⎛                                                                  2                                                                    2⎞
   ⎜⎛                                   d                    d       ⎞    ⎛                                  d                    d       ⎞ ⎟
   m⋅⎜⎜-ω⋅r(t)⋅sin(θ(t)) - r(t)⋅sin(θ(t))⋅──(θ(t)) + cos(θ(t))⋅──(r(t))⎟  + ⎜ω⋅r(t)⋅cos(θ(t)) + r(t)⋅cos(θ(t))⋅──(θ(t)) + sin(θ(t))⋅──(r(t))⎟ ⎟
   ⎝⎝                                   dt                   dt      ⎠    ⎝                                  dt                   dt      ⎠ ⎠
   ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
                                                                        2

   >>> lagranges_method = me.LagrangesMethod(L, coordinates, forcelist=loads, frame=N)
   >>> lagranges_method.form_lagranges_equations()

The obtained Lagrangian matches our theoretical Lagrangian:

.. math::

   L = \frac{1}{2} m \left( \dot{r}^2 + r^2 (\dot{\theta} + \omega)^2 \right)

The resulting equations of motion are nonlinear and coupled, due to friction and the interaction
between radial and tangential components.

**Equation 1**: Radial

.. math::

   m \ddot{r} = m r (\omega + \dot{\theta})^2 + F_r

**Equation 2**: Angular

.. math::

   m \left( 2 r \dot{r} (\omega + \dot{\theta}) + r^2 \ddot{\theta} \right) = r F_\theta
