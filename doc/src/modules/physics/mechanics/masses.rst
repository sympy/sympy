.. _masses:

=================================================================
Masses, Inertias, Particles and Rigid Bodies in Physics/Mechanics
=================================================================

This document will describe how to represent masses and inertias in
:mod:`sympy.physics.mechanics` and use of the :class:`~.RigidBody` and
:class:`~.Particle` classes.

It is assumed that the reader is familiar with the basics of these topics, such
as finding the center of mass for a system of particles, how to manipulate an
inertia tensor, and the definition of a particle and rigid body. Any advanced
dynamics text can provide a reference for these details.

Mass
====

The only requirement for a mass is that it needs to be a ``sympify``-able
expression. Keep in mind that masses can be time varying.

Particle
========

Particles are created with the class :class:`~.Particle` in
:mod:`sympy.physics.mechanics`. A :class:`~.Particle` object has an associated
point and an associated mass which are the only two attributes of the object.::

  >>> from sympy.physics.mechanics import Particle, Point
  >>> from sympy import Symbol
  >>> m = Symbol('m')
  >>> po = Point('po')
  >>> # create a particle container
  >>> pa = Particle('pa', po, m)

The associated point contains the position, velocity and acceleration of the
particle. :mod:`sympy.physics.mechanics` allows one to perform kinematic
analysis of points separate from their association with masses.

Inertia
=======

Inertia consists out of two parts: a quantity and a reference. The quantity is
expressed as a :class:`Dyadic<sympy.physics.vector.dyadic.Dyadic>` and the
reference is a :class:`Point<sympy.physics.vector.point.Point>`. The
:class:`Dyadic<sympy.physics.vector.dyadic.Dyadic>` can be defined as the outer
product between two vectors, which returns the juxtaposition of these vectors.
For further information, please refer to the :ref:`Dyadic` section in the
advanced documentation of the :mod:`sympy.physics.vector` module. Another more
intuitive method to define the
:class:`Dyadic<sympy.physics.vector.dyadic.Dyadic>` is to use the
:func:`~.inertia` function as described below in the section
'Inertia (Dyadics)'. The :class:`Point<sympy.physics.vector.point.Point>` about
which the :class:`Dyadic<sympy.physics.vector.dyadic.Dyadic>` is specified can
be any point, as long as it is defined with respect to the center of mass. The
most common reference point is of course the center of mass itself.

The inertia of a body can be specified using either an :class:`~.Inertia` object
or a ``tuple``. If a ``tuple`` is used, then it should have a length of two,
with the first entry being a :class:`Dyadic<sympy.physics.vector.dyadic.Dyadic>`
and the second entry being a :class:`Point<sympy.physics.vector.point.Point>`
about which the inertia dyadic is defined. Internally this ``tuple`` gets
converted to an :class:`~.Inertia` object. An example of using a ``tuple`` about
the center of mass is given below in the 'Rigid Body' section. The
:class:`~.Inertia` object can be created as follows.::

   >>> from sympy.physics.mechanics import ReferenceFrame, Point, outer, Inertia
   >>> A = ReferenceFrame('A')
   >>> P = Point('P')
   >>> Inertia(P, outer(A.x, A.x))
   ((A.x|A.x), P)


Inertia (Dyadics)
=================

A dyadic tensor is a second order tensor formed by the juxtaposition of a pair
of vectors. There are various operations defined with respect to dyadics,
which have been implemented in :obj:`~.sympy.physics.vector` in the form of
class :class:`Dyadic<sympy.physics.vector.dyadic.Dyadic>`. To know more, refer
to the :obj:`sympy.physics.vector.dyadic.Dyadic` and
:obj:`sympy.physics.vector.vector.Vector` class APIs. Dyadics are used to
define the inertia of bodies within :mod:`sympy.physics.mechanics`. Inertia
dyadics can be defined explicitly using the outer product, but the
:func:`~.inertia` function is typically much more convenient for the user.::

  >>> from sympy.physics.mechanics import ReferenceFrame, inertia
  >>> N = ReferenceFrame('N')

  Supply a reference frame and the moments of inertia if the object
  is symmetrical:

  >>> inertia(N, 1, 2, 3)
  (N.x|N.x) + 2*(N.y|N.y) + 3*(N.z|N.z)

  Supply a reference frame along with the products and moments of inertia
  for a general object:

  >>> inertia(N, 1, 2, 3, 4, 5, 6)
  (N.x|N.x) + 4*(N.x|N.y) + 6*(N.x|N.z) + 4*(N.y|N.x) + 2*(N.y|N.y) + 5*(N.y|N.z) + 6*(N.z|N.x) + 5*(N.z|N.y) + 3*(N.z|N.z)

Notice that the :func:`~.inertia` function returns a dyadic with each component
represented as two unit vectors separated by a ``|`` (outer product). Refer to
the :obj:`sympy.physics.vector.dyadic.Dyadic` section for more information about
dyadics.

Inertia is often expressed in a matrix, or tensor, form, especially for
numerical purposes. Since the matrix form does not contain any information
about the reference frame(s) the inertia dyadic is defined in, you must provide
one or two reference frames to extract the measure numbers from the dyadic.
There is a convenience function to do this::

  >>> inertia(N, 1, 2, 3, 4, 5, 6).to_matrix(N)
  Matrix([
  [1, 4, 6],
  [4, 2, 5],
  [6, 5, 3]])

Rigid Body
==========

Rigid bodies are created in a similar fashion as particles. The
:class:`~.RigidBody` class generates objects with four attributes: mass, center
of mass, a reference frame, and an :class:`~.Inertia` (a ``tuple`` can be passed
as well).::

  >>> from sympy import Symbol
  >>> from sympy.physics.mechanics import ReferenceFrame, Point, RigidBody
  >>> from sympy.physics.mechanics import outer
  >>> m = Symbol('m')
  >>> A = ReferenceFrame('A')
  >>> P = Point('P')
  >>> I = outer(A.x, A.x)
  >>> # create a rigid body
  >>> B = RigidBody('B', P, A, m, (I, P))

The mass is specified exactly as is in a particle. Similar to the
:class:`~.Particle`'s ``.point``, the :class:`~.RigidBody`'s center of mass,
``.masscenter`` must be specified. The reference frame is stored in an analogous
fashion and holds information about the body's orientation and angular velocity.

Loads
=====

In :mod:`sympy.physics.mechanics` loads can either be represented with tuples or
with the dedicated classes :class:`~.Force` and :class:`~.Torque`. Generally the
first argument (or item in the case of a tuple) is the location of the load. The
second argument is the vector. In the case of a force the first argument is a
point and the second a vector.

   >>> from sympy.physics.mechanics import Point, ReferenceFrame, Force
   >>> N = ReferenceFrame('N')
   >>> Po = Point('Po')
   >>> Force(Po, N.x)
   (Po, N.x)

The location of a torque, on the other hand, is a frame.

   >>> from sympy.physics.mechanics import Torque
   >>> Torque(N, 2 * N.x)
   (N, 2*N.x)

Optionally, one can also pass the body when using dedicated classes. If so,
the force will use the center of mass and the torque will use the associated
frame.

   >>> from sympy.physics.mechanics import RigidBody
   >>> rb = RigidBody('rb')
   >>> Force(rb, 3 * N.x)
   (rb_masscenter, 3*N.x)
   >>> Torque(rb, 4 * N.x)
   (rb_frame, 4*N.x)

Linear Momentum
===============

The linear momentum of a particle P is defined as:

.. math::
  L_P = m\mathbf{v}

where :math:`m` is the mass of the particle P and :math:`\mathbf{v}` is the
velocity of the particle in the inertial frame.[Likins1973]_.

Similarly the linear momentum of a rigid body is defined as:

.. math::
  L_B = m\mathbf{v^*}

where :math:`m` is the mass of the rigid body, B, and :math:`\mathbf{v^*}` is
the velocity of the mass center of B in the inertial frame.

Angular Momentum
================

The angular momentum of a particle P about an arbitrary point O in an inertial
frame N is defined as:

.. math::
  ^N \mathbf{H} ^ {P/O} = \mathbf{r} \times m\mathbf{v}

where :math:`\mathbf{r}` is a position vector from point O to the particle of
mass :math:`m` and :math:`\mathbf{v}` is the velocity of the particle in the
inertial frame.

Similarly the angular momentum of a rigid body B about a point O in an inertial
frame N is defined as:

.. math::
  ^N \mathbf{H} ^ {B/O} = ^N \mathbf{H} ^ {B/B^*} + ^N \mathbf{H} ^ {B^*/O}

where the angular momentum of the body about it's mass center is:

.. math::
  ^N \mathbf{H} ^ {B/B^*} = \mathbf{I^*} \cdot \omega

and the angular momentum of the mass center about O is:

.. math::
  ^N \mathbf{H} ^ {B^*/O} = \mathbf{r^*} \times m \mathbf{v^*}

where :math:`\mathbf{I^*}` is the central inertia dyadic of rigid body B,
:math:`\omega` is the inertial angular velocity of B, :math:`\mathbf{r^*}` is a
position vector from point O to the mass center of B, :math:`m` is the mass of
B and :math:`\mathbf{v^*}` is the velocity of the mass center in the inertial
frame.

Using momenta functions in Mechanics
====================================

The following example shows how to use the momenta functions in
:mod:`sympy.physics.mechanics`.

One begins by creating the requisite symbols to describe the system. Then
the reference frame is created and the kinematics are done. ::

  >>> from sympy import symbols
  >>> from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame
  >>> from sympy.physics.mechanics import RigidBody, Particle, Point, outer
  >>> from sympy.physics.mechanics import linear_momentum, angular_momentum
  >>> from sympy.physics.vector import init_vprinting
  >>> init_vprinting(pretty_print=False)
  >>> m, M, l1 = symbols('m M l1')
  >>> q1d = dynamicsymbols('q1d')
  >>> N = ReferenceFrame('N')
  >>> O = Point('O')
  >>> O.set_vel(N, 0 * N.x)
  >>> Ac = O.locatenew('Ac', l1 * N.x)
  >>> P = Ac.locatenew('P', l1 * N.x)
  >>> a = ReferenceFrame('a')
  >>> a.set_ang_vel(N, q1d * N.z)
  >>> Ac.v2pt_theory(O, N, a)
  l1*q1d*N.y
  >>> P.v2pt_theory(O, N, a)
  2*l1*q1d*N.y

Finally, the bodies that make up the system are created. In this case the
system consists of a particle Pa and a RigidBody A. ::

  >>> Pa = Particle('Pa', P, m)
  >>> I = outer(N.z, N.z)
  >>> A = RigidBody('A', Ac, a, M, (I, Ac))

Then one can either choose to evaluate the momenta of individual components
of the system or of the entire system itself. ::

  >>> linear_momentum(N,A)
  M*l1*q1d*N.y
  >>> angular_momentum(O, N, Pa)
  4*l1**2*m*q1d*N.z
  >>> linear_momentum(N, A, Pa)
  (M*l1*q1d + 2*l1*m*q1d)*N.y
  >>> angular_momentum(O, N, A, Pa)
  (M*l1**2*q1d + 4*l1**2*m*q1d + q1d)*N.z

It should be noted that the user can determine either momenta in any frame
in :mod:`sympy.physics.mechanics` as the user is allowed to specify the reference frame when
calling the function. In other words the user is not limited to determining
just inertial linear and angular momenta. Please refer to the docstrings on
each function to learn more about how each function works precisely.

Kinetic Energy
==============

The kinetic energy of a particle P is defined as

.. math::
  T_P = \frac{1}{2} m \mathbf{v^2}

where :math:`m` is the mass of the particle P and :math:`\mathbf{v}`
is the velocity of the particle in the inertial frame.

Similarly the kinetic energy of a rigid body B is defined as

.. math::
  T_B = T_t + T_r

where the translational kinetic energy is given by:

.. math::
  T_t = \frac{1}{2} m \mathbf{v^*} \cdot \mathbf{v^*}

and the rotational kinetic energy is given by:

.. math::
  T_r = \frac{1}{2} \omega \cdot \mathbf{I^*} \cdot \omega

where :math:`m` is the mass of the rigid body, :math:`\mathbf{v^*}` is the
velocity of the mass center in the inertial frame, :math:`\omega` is the
inertial angular velocity of the body and :math:`\mathbf{I^*}` is the central
inertia dyadic.

Potential Energy
================

Potential energy is defined as the energy possessed by a body or system by
virtue of its position or arrangement.

Since there are a variety of definitions for potential energy, this is not
discussed further here. One can learn more about this in any elementary text
book on dynamics.

Lagrangian
==========

The Lagrangian of a body or a system of bodies is defined as:

.. math::
   \mathcal{L} = T - V

where :math:`T` and :math:`V` are the kinetic and potential energies
respectively.

Using energy functions in Mechanics
===================================

The following example shows how to use the energy functions in
:mod:`sympy.physics.mechanics`.

As was discussed above in the momenta functions, one first creates the system
by going through an identical procedure. ::

  >>> from sympy import symbols
  >>> from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, outer
  >>> from sympy.physics.mechanics import RigidBody, Particle
  >>> from sympy.physics.mechanics import kinetic_energy, potential_energy, Point
  >>> from sympy.physics.vector import init_vprinting
  >>> init_vprinting(pretty_print=False)
  >>> m, M, l1, g, h, H = symbols('m M l1 g h H')
  >>> omega = dynamicsymbols('omega')
  >>> N = ReferenceFrame('N')
  >>> O = Point('O')
  >>> O.set_vel(N, 0 * N.x)
  >>> Ac = O.locatenew('Ac', l1 * N.x)
  >>> P = Ac.locatenew('P', l1 * N.x)
  >>> a = ReferenceFrame('a')
  >>> a.set_ang_vel(N, omega * N.z)
  >>> Ac.v2pt_theory(O, N, a)
  l1*omega*N.y
  >>> P.v2pt_theory(O, N, a)
  2*l1*omega*N.y
  >>> Pa = Particle('Pa', P, m)
  >>> I = outer(N.z, N.z)
  >>> A = RigidBody('A', Ac, a, M, (I, Ac))

The user can then determine the kinetic energy of any number of entities of the
system: ::

  >>> kinetic_energy(N, Pa)
  2*l1**2*m*omega**2
  >>> kinetic_energy(N, Pa, A)
  M*l1**2*omega**2/2 + 2*l1**2*m*omega**2 + omega**2/2

It should be noted that the user can determine either kinetic energy relative
to any frame in :mod:`sympy.physics.mechanics` as the user is allowed to specify the
reference frame when calling the function. In other words the user is not
limited to determining just inertial kinetic energy.

For potential energies, the user must first specify the potential energy of
every entity of the system using the
:obj:`sympy.physics.mechanics.rigidbody.RigidBody.potential_energy` property.
The potential energy of any number of entities comprising the system can then
be determined: ::

  >>> Pa.potential_energy = m * g * h
  >>> A.potential_energy = M * g * H
  >>> potential_energy(A, Pa)
  H*M*g + g*h*m

One can also determine the Lagrangian for this system: ::

  >>> from sympy.physics.mechanics import Lagrangian
  >>> from sympy.physics.vector import init_vprinting
  >>> init_vprinting(pretty_print=False)
  >>> Lagrangian(N, Pa, A)
  -H*M*g + M*l1**2*omega**2/2 - g*h*m + 2*l1**2*m*omega**2 + omega**2/2

Please refer to the docstrings to learn more about each function.
