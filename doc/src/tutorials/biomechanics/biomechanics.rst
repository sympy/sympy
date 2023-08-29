.. _biomechanics-tutorial:

======================================
Introduction to Biomechanical Modeling
======================================

:obj:`~sympy.physics._biomechanics` provides features to enhance models created
with :obj:`~sympy.physics.mechanics` with force producing elements that model
muscles and tendons. In this tutorial, we will introduce the features of this
package.

The initial primary purpose of the biomechanics package is to introduce tools
for modeling the forces produced by `Hill-type muscle models`_. These models
generate forces applied to the skeletal structure of an organism based on the
contraction state of the muscle coupled with the passive stretch of tendons. In
this tutorial, we introduce the elements that make up a musculotendon model and
then demonstrate it in operation with a specific implementation, the
:obj:`~sympy.physics._biomechanics.MusculotendonDeGroote2016` model.

.. _Hill-type muscle models: https://en.wikipedia.org/wiki/Hill%27s_muscle_model

Loads
=====

:obj:`~sympy.physics.mechanics.loads` includes two types of loads:
:obj:`~sympy.physics.mechanics.loads.Force` and
:obj:`~sympy.physics.mechanics.loads.Torque`. Forces represent bound vector
quantities that act directed along a line of action and torques are unbound
vectors which represent the resulting torque of a couple from a set of forces.

An example of very common force model is a linear spring and linear damper in
parallel. The force acting on a particle of mass :math:`m` with 1D motion
described by generalized coordinate :math:`x(t)`  with linear spring and damper
coefficients :math:`k` and :math:`c` has the familiar equation of motion:

.. math::

   m \ddot{x} = \sum F = -kx - c\dot{x}

In SymPy, we can formulate the force acting on the particle :math:`P` that has
motion in reference frame :math:`N` and position relative to point :math:`O`
fixed in :math:`N` like so:

.. plot::
   :format: doctest
   :include-source: True
   :context: reset
   :nofigs:

   >>> import sympy as sm
   >>> import sympy.physics.mechanics as me

   >>> k, c = sm.symbols('k, c', real=True, nonnegative=True)
   >>> x = me.dynamicsymbols('x', real=True)

   >>> N = me.ReferenceFrame('N')
   >>> O, P = me.Point('O'), me.Point('P')

   >>> P.set_pos(O, x*N.x)
   >>> P.set_vel(N, x.diff()*N.x)

   >>> force_on_P = me.Force(P, -k*P.pos_from(O) - c*P.vel(N))
   >>> force_on_P
   (P, (-c*Derivative(x(t), t) - k*x(t))*N.x)

and there would be an equal and opposite force acting on :math:`O`:

.. plot::
   :format: doctest
   :include-source: True
   :context:
   :nofigs:

   >>> force_on_O = me.Force(O, k*P.pos_from(O) + c*P.vel(N))
   >>> force_on_O
   (O, (c*Derivative(x(t), t) + k*x(t))*N.x)

Forces that a single muscle and tendon applies to a set of rigid bodies will
be the primary output of the musculotendon models developed further in this
tutorial.

Pathways
========

Muscles and their associated tendons wrap around the moving skeletal system and
other muscles and organs. This imposes the challenge of determining the lines
of action of the forces that the muscle and tendon produce on the skeleton and
organs it touches. We have introduced the
:obj:`~sympy.physics.mechanics.pathway` module to help manage the specification
of the geometric relationships to the forces' lines of action.

The spring-damper example above has the simplest line of action definition so
we can use a :obj:`~sympy.physics.mechanics.pathway.LinearPathway` to establish
that line of action. First provide the two endpoints where the force will have
equal and opposite application to and the distance between the points and the
relative speed between the two points are calculated by the pathway with
:obj:`~sympy.physics.mechanics.pathway.LinearPathway.length` and
:obj:`~sympy.physics.mechanics.pathway.LinearPathway.extension_velocity`. Note
that a positive speed implies the points are moving away from each other. Also
note that the formulation handles the case where :math:`x` is positive or
negative.

.. plot::
   :format: doctest
   :include-source: True
   :context:
   :nofigs:

   >>> lpathway = me.LinearPathway(O, P)
   >>> lpathway
   LinearPathway(O, P)
   >>> lpathway.length
   Abs(x(t))
   >>> lpathway.extension_velocity
   sign(x(t))*Derivative(x(t), t)

The :obj:`~sympy.physics.mechanics.pathway.LinearPathway.to_loads` method then
takes the magnitude of a force with a sign convention that positive magnitudes
push the two points away from each other and returns a list of all forces
acting on the two points.

.. plot::
   :format: doctest
   :include-source: True
   :context:
   :nofigs:

   >>> import pprint
   >>> pprint.pprint(lpathway.to_loads(-k*x - k*x.diff()))
   [Force(point=O, force=(k*x(t) + k*Derivative(x(t), t))*x(t)/Abs(x(t))*N.x),
    Force(point=P, force=(-k*x(t) - k*Derivative(x(t), t))*x(t)/Abs(x(t))*N.x)]

Pathways can be constructed with any arbitrary geometry and any number of
interconnected particles and rigid bodies. An example, a more complicated
pathway is an :obj:`~sympy.physics.mechanics.pathway.ObstacleSetPathway`. You
can specify any number of intermediate points between the two pathway endpoints
which the actuation path of the forces will follow along. For example, if we
introduce two points fixed in :math:`N` then the force will act along a set of
linear segments connecting :math:`O` to :math:`Q` to :math:`R`: then to
:math:`P`. Each of the four points will experience resultant forces. For
simplicity we show the effect of only the spring force.

.. plot::
   :format: doctest
   :include-source: True
   :context:
   :nofigs:

   >>> Q, R = me.Point('Q'), me.Point('R')
   >>> Q.set_pos(O, 1*N.y)
   >>> R.set_pos(O, 1*N.x + 1*N.y)
   >>> opathway = me.ObstacleSetPathway(O, Q, R, P)
   >>> opathway.length
   sqrt((x(t) - 1)**2 + 1) + 2
   >>> opathway.extension_velocity
   (x(t) - 1)*Derivative(x(t), t)/sqrt((x(t) - 1)**2 + 1)
   >>> pprint.pprint(opathway.to_loads(-k*opathway.length))
   [Force(point=O, force=k*(sqrt((x(t) - 1)**2 + 1) + 2)*N.y),
    Force(point=Q, force=- k*(sqrt((x(t) - 1)**2 + 1) + 2)*N.y),
    Force(point=Q, force=k*(sqrt((x(t) - 1)**2 + 1) + 2)*N.x),
    Force(point=R, force=- k*(sqrt((x(t) - 1)**2 + 1) + 2)*N.x),
    Force(point=R, force=k*(sqrt((x(t) - 1)**2 + 1) + 2)*(x(t) - 1)/sqrt((x(t) - 1)**2 + 1)*N.x - k*(sqrt((x(t) - 1)**2 + 1) + 2)/sqrt((x(t) - 1)**2 + 1)*N.y),
    Force(point=P, force=- k*(sqrt((x(t) - 1)**2 + 1) + 2)*(x(t) - 1)/sqrt((x(t) - 1)**2 + 1)*N.x + k*(sqrt((x(t) - 1)**2 + 1) + 2)/sqrt((x(t) - 1)**2 + 1)*N.y)]

If you set :math:`x=1`, it is a bit easier to see that the collection of forces
are correct:

.. plot::
   :format: doctest
   :include-source: True
   :context:
   :nofigs:

   >>> for load in opathway.to_loads(-k*opathway.length):
   ...     pprint.pprint(me.Force(load[0], load[1].subs({x: 1})))
   Force(point=O, force=3*k*N.y)
   Force(point=Q, force=- 3*k*N.y)
   Force(point=Q, force=3*k*N.x)
   Force(point=R, force=- 3*k*N.x)
   Force(point=R, force=- 3*k*N.y)
   Force(point=P, force=3*k*N.y)

You can create your own pathways by sub-classing
:obj:`~sympy.physics.mechanics.pathway.PathwayBase`.

Wrapping Geometries
===================

It is common for muscles to wrap over bones, tissue, or organs. We have
introduced wrapping geometries and associated wrapping pathways to help manage
their complexities. For example, if two pathway endpoints lie on the surface of
a cylinder the forces act along lines that are tangent to the geodesic
connecting the two points at the endpoints. The
:obj:`~sympy.physics.mechanics.wrapping_geometry.WrappingCylinder` object
calculates the complex geometry for the pathway. A
:obj:`~sympy.physics.mechanics.pathway.WrappingPathway` then uses the geometry
to construct the forces. A spring force along this pathway is constructed
below:

.. plot::
   :format: doctest
   :include-source: True
   :context:
   :nofigs:

   >>> r = sm.symbols('r', real=True, nonegative=True)
   >>> theta = me.dynamicsymbols('theta', real=True)
   >>> O, P, Q = sm.symbols('O, P, Q', cls=me.Point)
   >>> A = me.ReferenceFrame('A')

   >>> A.orient_axis(N, theta, N.z)

   >>> P.set_pos(O, r*N.x)
   >>> Q.set_pos(O, N.z + r*A.x)

   >>> cyl = me.WrappingCylinder(r, O, N.z)
   >>> wpathway = me.WrappingPathway(P, Q, cyl)
   >>> pprint.pprint(wpathway.to_loads(-k*wpathway.length))
   [Force(point=P, force=- k*r*Abs(theta(t))*N.y - k*N.z),
    Force(point=Q, force=k*N.z + k*r*Abs(theta(t))*A.y),
    Force(point=O, force=k*r*Abs(theta(t))*N.y - k*r*Abs(theta(t))*A.y)]

Actuators
=========

Models of multibody systems commonly have time varying inputs in the form of
the magnitudes of forces or torques. In many cases, these specified inputs may
be derived from the state of the system or even from the output of another
dynamic system. The :obj:`sympy.physics.mechanics.actuator` module includes
classes to help manage the creation of such models of force and torque inputs.
An actuator is intended to represent a real physical component. For example,
the spring-damper force from above can be created by sub-classing
:obj:`sympy.physics.mechanics.ActuatorBase` and developing a method that
generates the loads associated with that spring-damper actuator.

.. plot::
   :format: doctest
   :include-source: True
   :context:
   :nofigs:

   >>> N = me.ReferenceFrame('N')
   >>> O, P = me.Point('O'), me.Point('P')
   >>> P.set_pos(O, x*N.x)

   >>> class SpringDamper(me.ActuatorBase):
   ...     # positive x spring is in tension
   ...     # negative x spring is in compression
   ...     def __init__(self, P1, P2, spring_constant, damper_constant):
   ...         self.P1 = P1
   ...         self.P2 = P2
   ...         self.k = spring_constant
   ...         self.c = damper_constant
   ...     def to_loads(self):
   ...         x = self.P2.pos_from(self.P1).magnitude()
   ...         v = x.diff(me.dynamicsymbols._t)
   ...         dir_vec = self.P2.pos_from(self.P1).normalize()
   ...         force_P1 = me.Force(self.P1,
   ...                             self.k*x*dir_vec + self.c*v*dir_vec)
   ...         force_P2 = me.Force(self.P2,
   ...                             -self.k*x*dir_vec - self.c*v*dir_vec)
   ...         return [force_P1, force_P2]
   ...

   >>> spring_damper = SpringDamper(O, P, k, c)
   >>> pprint.pprint(spring_damper.to_loads())
   [Force(point=O, force=(c*x(t)*sign(x(t))*Derivative(x(t), t)/Abs(x(t)) + k*x(t))*N.x),
    Force(point=P, force=(-c*x(t)*sign(x(t))*Derivative(x(t), t)/Abs(x(t)) - k*x(t))*N.x)]

There is also a :obj:`sympy.physics.mechanics.actuator.ForceActuator` that
allows seamless integration with pathway objects. You only need to set the
``.force`` attribute in initialization in the sub-class.

.. plot::
   :format: doctest
   :include-source: True
   :context:
   :nofigs:

   >>> class SpringDamper(me.ForceActuator):
   ...     # positive x spring is in tension
   ...     # negative x spring is in compression
   ...     def __init__(self, pathway, spring_constant, damper_constant):
   ...         self.pathway = pathway
   ...         self.force = (-spring_constant*pathway.length -
   ...                       damper_constant*pathway.extension_velocity)
   ...
   >>> spring_damper2 = SpringDamper(lpathway, k, c)
   >>> pprint.pprint(spring_damper2.to_loads())
   [Force(point=O, force=(c*sign(x(t))*Derivative(x(t), t) + k*Abs(x(t)))*x(t)/Abs(x(t))*N.x),
    Force(point=P, force=(-c*sign(x(t))*Derivative(x(t), t) - k*Abs(x(t)))*x(t)/Abs(x(t))*N.x)]

This then makes it easy to apply a spring-damper force to other pathways, e.g.:

.. plot::
   :format: doctest
   :include-source: True
   :context:
   :nofigs:

   >>> spring_damper3 = SpringDamper(wpathway, k, c)
   >>> pprint.pprint(spring_damper3.to_loads())
   [Force(point=P, force=r*(-c*r**2*theta(t)*Derivative(theta(t), t)/sqrt(r**2*theta(t)**2 + 1) - k*sqrt(r**2*theta(t)**2 + 1))*Abs(theta(t))/sqrt(r**2*theta(t)**2 + 1)*N.y + (-c*r**2*theta(t)*Derivative(theta(t), t)/sqrt(r**2*theta(t)**2 + 1) - k*sqrt(r**2*theta(t)**2 + 1))/sqrt(r**2*theta(t)**2 + 1)*N.z),
    Force(point=Q, force=- (-c*r**2*theta(t)*Derivative(theta(t), t)/sqrt(r**2*theta(t)**2 + 1) - k*sqrt(r**2*theta(t)**2 + 1))/sqrt(r**2*theta(t)**2 + 1)*N.z - r*(-c*r**2*theta(t)*Derivative(theta(t), t)/sqrt(r**2*theta(t)**2 + 1) - k*sqrt(r**2*theta(t)**2 + 1))*Abs(theta(t))/sqrt(r**2*theta(t)**2 + 1)*A.y),
    Force(point=O, force=- r*(-c*r**2*theta(t)*Derivative(theta(t), t)/sqrt(r**2*theta(t)**2 + 1) - k*sqrt(r**2*theta(t)**2 + 1))*Abs(theta(t))/sqrt(r**2*theta(t)**2 + 1)*N.y + r*(-c*r**2*theta(t)*Derivative(theta(t), t)/sqrt(r**2*theta(t)**2 + 1) - k*sqrt(r**2*theta(t)**2 + 1))*Abs(theta(t))/sqrt(r**2*theta(t)**2 + 1)*A.y)]

Activation Dynamics
===================

Musculotendon Curves
====================

Musculotendon Actuators
=======================

A Simple Musculotendon Model
============================

To demonstrate a muscle's effect on a simple system, we can model a particle of
mass :math:`m` under the influence of gravity with a muscle pulling the mass
against gravity. The mass :math:`m` has a single generalized coordinate
:math:`q` and generalized speed :math:`u` to describe its position and motion.
The following code establishes the kinematics and gravitational force and an
associated particle:

.. plot::
   :format: doctest
   :include-source: True
   :context: reset
   :nofigs:

   >>> import sympy as sm
   >>> import sympy.physics.mechanics as me

   >>> q, u = me.dynamicsymbols('q, u', real=True)
   >>> m, g = sm.symbols('m, g', real=True, positive=True)

   >>> N = me.ReferenceFrame('N')
   >>> O, P = sm.symbols('O, P', cls=me.Point)

   >>> P.set_pos(O, q*N.x)
   >>> O.set_vel(N, 0)
   >>> P.set_vel(N, u*N.x)

   >>> gravity = me.Force(P, m*g*N.x)

   >>> block = me.Particle('block', P, m)

SymPy Biomechanics includes musculotendon actuator models. Here we will use a
specific musculotendon model implementation. A musculotendon actuator is
instantiated with two input components, the pathway and the activation dynamics
model. The actuator must act along a pathway that connects the origin and
insertion points of the muscle. Our origin will attach to the fixed point
:math:`O` and insert on the moving particle :math:`P`.

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs
   :nofigs:

   >>> from sympy.physics.mechanics.pathway import LinearPathway

   >>> muscle_pathway = LinearPathway(O, P)

A pathway has attachment points:

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs
   :nofigs:

   >>> muscle_pathway.attachments
   (O, P)

TODO : note the sign conventions

and knows the length between the end attachment points as well as the relative
speed between the two attachment points:

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs
   :nofigs:

   >>> muscle_pathway.length
   Abs(q(t))
   >>> muscle_pathway.extension_velocity
   sign(q(t))*Derivative(q(t), t)

The sign convention is that a positive speed

Finally, the pathway can determine the forces acting on the two attachment
points give a force magnitude:

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs
   :nofigs:

   >>> muscle_pathway.to_loads(m*g)
   [(O, - g*m*q(t)/Abs(q(t))*N.x), (P, g*m*q(t)/Abs(q(t))*N.x)]

The activation dynamics model represents a set of algebraic or ordinary
differential equations that relate the muscle excitation to the muscle
activation. In our case, we will use a first order ordinary differential
equation that gives a smooth, but delayed activation :math:`a(t)` from the
excitation :math:`e(t)`.

TODO : We could plot dadt as a function of a for different e from 0 to 1.

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs
   :nofigs:

   >>> from sympy.physics._biomechanics import FirstOrderActivationDeGroote2016
   >>> muscle_activation = FirstOrderActivationDeGroote2016.with_defaults('muscle')

The activation model has a state variable, input variable, and some constant
parameters:

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs
   :nofigs:

   >>> muscle_activation.x
   Matrix([[a_muscle(t)]])
   >>> muscle_activation.r
   Matrix([[e_muscle(t)]])
   >>> muscle_activation.p
   Matrix([
   [0.015],
   [ 0.06],
   [   10]])

These are associated with its first order differential equation:

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs
   :nofigs:

   >>> muscle_activation.rhs()
   Matrix([[((1/2 - tanh(10*a_muscle(t) - 10*e_muscle(t))/2)/(0.0225*a_muscle(t) + 0.0075) + 16.6666666666667*(3*a_muscle(t)/2 + 1/2)*(tanh(10*a_muscle(t) - 10*e_muscle(t))/2 + 1/2))*(-a_muscle(t) + e_muscle(t))]])

With the pathway and activation dynamics, the musculotendon model created using
them both and needs some parameters to define the muscle and tendon specific
properties. You need to specify the tendon slack length, peak isometric force,
optimal fiber length, maximal fiber velocity, optimal pennation angle, and
fiber damping coefficients.

TODO : How do we know this is a rigid tendon model?

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs
   :nofigs:

   >>> from sympy.physics._biomechanics import MusculotendonDeGroote2016

   >>> F_M_max, l_M_opt, l_T_slack = sm.symbols('F_M_max, l_M_opt, l_T_slack', real=True)
   >>> v_M_max, alpha_opt, beta = sm.symbols('v_M_max, alpha_opt, beta', real=True)

   >>> muscle = MusculotendonDeGroote2016(
   ...     'muscle',
   ...     muscle_pathway,
   ...     muscle_activation,
   ...     tendon_slack_length=l_T_slack,
   ...     peak_isometric_force=F_M_max,
   ...     optimal_fiber_length=l_M_opt,
   ...     maximal_fiber_velocity=v_M_max,
   ...     optimal_pennation_angle=alpha_opt,
   ...     fiber_damping_coefficient=beta,
   ... )
   ...

TODO : Explain why the rhs() is different for the muscle than the activation.
TODO : Needs explanation about rigid tendon

Because this musculotendon actuator has a rigid tendon model, it has the same
state and ordinary differential equation as the activation model:

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs
   :nofigs:

   >>> muscle.musculotendon_dynamics
   MusculotendonFormulation.RIGID_TENDON
   >>> muscle.x
   Matrix([[a_muscle(t)]])
   >>> muscle.r
   Matrix([[e_muscle(t)]])
   >>> muscle.p
   Matrix([
   [l_T_slack],
   [  F_M_max],
   [  l_M_opt],
   [  v_M_max],
   [alpha_opt],
   [     beta],
   [    0.015],
   [     0.06],
   [       10]])
   >>> muscle.rhs()
   Matrix([[(-0.5625*a_muscle(t)**3*tanh(10*a_muscle(t) - 10*e_muscle(t)) - 0.5625*a_muscle(t)**3 + 0.5625*a_muscle(t)**2*e_muscle(t)*tanh(10*a_muscle(t) - 10*e_muscle(t)) + 0.5625*a_muscle(t)**2*e_muscle(t) - 0.375*a_muscle(t)**2*tanh(10*a_muscle(t) - 10*e_muscle(t)) - 0.375*a_muscle(t)**2 + 0.375*a_muscle(t)*e_muscle(t)*tanh(10*a_muscle(t) - 10*e_muscle(t)) + 0.375*a_muscle(t)*e_muscle(t) + 0.9375*a_muscle(t)*tanh(10*a_muscle(t) - 10*e_muscle(t)) - 1.0625*a_muscle(t) - 0.9375*e_muscle(t)*tanh(10*a_muscle(t) - 10*e_muscle(t)) + 1.0625*e_muscle(t))/(0.045*a_muscle(t) + 0.015)]])

The musculotendon provides the extra ordinary differential equations as well as
the muscle specific forces applied to the pathway:

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs
   :nofigs:

   >>> muscle_loads = muscle.to_loads()
   >>> pprint.pprint(muscle_loads)
   [Force(point=O, force=F_M_max*(beta*(-l_T_slack + Abs(q(t)))*sign(q(t))*Derivative(q(t), t)/(v_M_max*sqrt(l_M_opt**2*sin(alpha_opt)**2 + (-l_T_slack + Abs(q(t)))**2)) + a_muscle(t)*FiberForceLengthActiveDeGroote2016(sqrt(l_M_opt**2*sin(alpha_opt)**2 + (-l_T_slack + Abs(q(t)))**2)/l_M_opt, 0.814, 1.06, 0.162, 0.0633, 0.433, 0.717, -0.0299, 1/5, 1/10, 1, 0.354, 0)*FiberForceVelocityDeGroote2016((-l_T_slack + Abs(q(t)))*sign(q(t))*Derivative(q(t), t)/(v_M_max*sqrt(l_M_opt**2*sin(alpha_opt)**2 + (-l_T_slack + Abs(q(t)))**2)), -0.318, -8.149, -0.374, 0.886) + FiberForceLengthPassiveDeGroote2016(sqrt(l_M_opt**2*sin(alpha_opt)**2 + (-l_T_slack + Abs(q(t)))**2)/l_M_opt, 3/5, 4))*q(t)/Abs(q(t))*N.x),
    Force(point=P, force=- F_M_max*(beta*(-l_T_slack + Abs(q(t)))*sign(q(t))*Derivative(q(t), t)/(v_M_max*sqrt(l_M_opt**2*sin(alpha_opt)**2 + (-l_T_slack + Abs(q(t)))**2)) + a_muscle(t)*FiberForceLengthActiveDeGroote2016(sqrt(l_M_opt**2*sin(alpha_opt)**2 + (-l_T_slack + Abs(q(t)))**2)/l_M_opt, 0.814, 1.06, 0.162, 0.0633, 0.433, 0.717, -0.0299, 1/5, 1/10, 1, 0.354, 0)*FiberForceVelocityDeGroote2016((-l_T_slack + Abs(q(t)))*sign(q(t))*Derivative(q(t), t)/(v_M_max*sqrt(l_M_opt**2*sin(alpha_opt)**2 + (-l_T_slack + Abs(q(t)))**2)), -0.318, -8.149, -0.374, 0.886) + FiberForceLengthPassiveDeGroote2016(sqrt(l_M_opt**2*sin(alpha_opt)**2 + (-l_T_slack + Abs(q(t)))**2)/l_M_opt, 3/5, 4))*q(t)/Abs(q(t))*N.x)]

These loads are made up of various functions that describe the length and
velocity relationships to the fiber force.

Now that we have the forces that the muscles and tendons produce the equations
of motion of the system can be formed with, for example, Kanes Method:

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs
   :nofigs:

   >>> kane = me.KanesMethod(N, (q,), (u,), kd_eqs=(u - q.diff(),))
   >>> Fr, Frs = kane.kanes_equations((block,), (muscle_loads + [gravity]))

The equations of motion are made up of the kinematical differential equation,
the dynamical differential equation (Newton's Second Law), and the muscle
activation differential equation. The explicit form of each can be formed like
so:

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs
   :nofigs:

   >>> dqdt = u
   >>> dudt = kane.forcing[0]/m
   >>> dadt = muscle.rhs()[0]

We can now create a numerical function that evaluates the equations of motion
given the state, inputs, and constant parameters. Start by listing each
symbolically:

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs
   :nofigs:

   >>> a = muscle.a
   >>> e = muscle.e
   >>> state = [q, u, a]
   >>> inputs = [e]
   >>> constants = [m, g, F_M_max, l_M_opt, l_T_slack, v_M_max, alpha_opt, beta]

Then the numerical function is:

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs
   :nofigs:

   >>> eval_eom = sm.lambdify((state, inputs, constants), (dqdt, dudt, dadt))

It will additionally be interesting to numerically evaluate the muscle force,
so create a function for it too:

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs
   :nofigs:

   >>> force = muscle.force.xreplace({q.diff(): u})
   >>> eval_force = sm.lambdify((state, constants), force)

To test these functions we need some suitable numerical values. This muscle
will be able to produce a maximum force of 10 N to lift a mass of 0.5 kg:

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs
   :nofigs:

   >>> import numpy as np
   >>> p_vals = np.array([
   ...     0.5,  # m [kg]
   ...     9.81,  # g [m/s/s]
   ...     10.0,  # F_M_max
   ...     0.18,  # l_M_opt, length of muscle at which max force is produced
   ...     0.17,  # l_T_slack, always fixed (rigid tendon)
   ...     10.0,  # v_M_max
   ...     0.0,  # alpha_opt
   ...     0.1,  # beta
   ... ])
   ...

Our tendon is rigid, so the length of the muscle will be :math:`q-l_T_slack`
and we want to give an initial muscle length near its force producing peak, so
we choose :math:`q_0=l_M_opt + l_T_slack`:

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs
   :nofigs:

   >>> x_vals = np.array([
   ...     p_vals[3] + p_vals[4],  # q [m]
   ...     0.0,  # u [m/s]
   ...     0.0,  # a [?]
   ... ])
   ...

We can set the excitation to zero to test the numerical functions:

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs
   :nofigs:

   >>> r_vals = np.array([
   ...     1.0,  # e
   ... ])
   ...
   >>> eval_eom(x_vals, r_vals, p_vals)
   (0.0, 9.81, 133.33333307568913)
   >>> eval_force(x_vals, p_vals)
   1.4499681738213515e-16

The two functions work so we can now simulate this system to see if and how the
muscle lifts the mass:

.. plot::
   :format: doctest
   :include-source: True
   :context: close-figs

   >>> def eval_rhs(t, x):
   ...
   ...     r = np.array([1.0])
   ...
   ...     return eval_eom(x, r, p_vals)
   ...

   >>> from scipy.integrate import solve_ivp
   >>> t0, tf = 0.0, 10.0
   >>> times = np.linspace(t0, tf, num=1001)
   >>> sol = solve_ivp(eval_rhs,
   ...                 (t0, tf),
   ...                 x_vals, t_eval=times)
   ...
   >>> import matplotlib.pyplot as plt
   >>> fig, axes = plt.subplots(4, 1, sharex=True)
   >>> axes[0].plot(sol.t, sol.y[0] - p_vals[4], label='length of muscle')
   >>> axes[0].set_ylabel('Distance [m]')
   >>> axes[1].plot(sol.t, sol.y[1], label=state[1])
   >>> axes[1].set_ylabel('Speed [m/s]')
   >>> axes[2].plot(sol.t, sol.y[2], label=state[2])
   >>> axes[2].set_ylabel('Activation')
   >>> axes[3].plot(sol.t, eval_force(sol.y, p_vals).T, label='force')
   >>> axes[3].set_ylabel('Force [N]')
   >>> axes[3].set_xlabel('Time [s]')
   >>> axes[0].legend(), axes[1].legend(), axes[2].legend(), axes[3].legend()
