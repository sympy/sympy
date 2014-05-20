=======================================================
Mechanics: Masses, Inertias, Particles and Rigid Bodies
=======================================================

This document will describe how to represent masses and inertias in
:mod:`mechanics` and use of the ``RigidBody`` and ``Particle`` classes.

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

Particles are created with the class ``Particle`` in :mod:`mechanics`.
A ``Particle`` object has an associated point and an associated mass which are
the only two attributes of the object.::

  >>> from sympy.physics.mechanics import Particle, Point
  >>> from sympy import Symbol
  >>> m = Symbol('m')
  >>> po = Point('po')
  >>> # create a particle container
  >>> pa = Particle('pa', po, m)

The associated point contains the position, velocity and acceleration of the
particle. :mod:`mechanics` allows one to perform kinematic analysis of points
separate from their association with masses.

Inertia
=======

See the Inertia (Dyadics) section in 'Advanced Topics' part of
:mod:`physics/vector` docs.

Rigid Body
==========

Rigid bodies are created in a similar fashion as particles. The ``RigidBody``
class generates objects with four attributes: mass, center of mass, a reference
frame, and an inertia tuple::

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
``Particle``'s ``.point``, the ``RigidBody``'s center of mass, ``.masscenter``
must be specified. The reference frame is stored in an analogous fashion and
holds information about the body's orientation and angular velocity. Finally,
the inertia for a rigid body needs to be specified about a point. In
:mod:`mechanics`, you are allowed to specify any point for this. The most
common is the center of mass, as shown in the above code. If a point is selected
which is not the center of mass, ensure that the position between the point and
the center of mass has been defined. The inertia is specified as a tuple of length
two with the first entry being a ``Dyadic`` and the second entry being a
``Point`` of which the inertia dyadic is defined about.

.. _Dyadic:

Dyadic
======

In :mod:`mechanics`, dyadics are used to represent inertia ([Kane1985]_,
[WikiDyadics]_, [WikiDyadicProducts]_). A dyadic is a linear polynomial of
component unit dyadics, similar to a vector being a linear polynomial of
component unit vectors. A dyadic is the outer product between two vectors which
returns a new quantity representing the juxtaposition of these two vectors. For
example:

.. math::
  \mathbf{\hat{a}_x} \otimes \mathbf{\hat{a}_x} &= \mathbf{\hat{a}_x}
  \mathbf{\hat{a}_x}\\
  \mathbf{\hat{a}_x} \otimes \mathbf{\hat{a}_y} &= \mathbf{\hat{a}_x}
  \mathbf{\hat{a}_y}\\

Where :math:`\mathbf{\hat{a}_x}\mathbf{\hat{a}_x}` and
`\mathbf{\hat{a}_x}\mathbf{\hat{a}_y}` are the outer products obtained by
multiplying the left side as a column vector by the right side as a row vector.
Note that the order is significant.

Some additional properties of a dyadic are:

.. math::
  (x \mathbf{v}) \otimes \mathbf{w} &= \mathbf{v} \otimes (x \mathbf{w}) = x
  (\mathbf{v} \otimes \mathbf{w})\\
  \mathbf{v} \otimes (\mathbf{w} + \mathbf{u}) &= \mathbf{v} \otimes \mathbf{w}
  + \mathbf{v} \otimes \mathbf{u}\\
  (\mathbf{v} + \mathbf{w}) \otimes \mathbf{u} &= \mathbf{v} \otimes \mathbf{u}
  + \mathbf{w} \otimes \mathbf{u}\\

A vector in a reference frame can be represented as
:math:`\begin{bmatrix}a\\b\\c\end{bmatrix}` or :math:`a \mathbf{\hat{i}} + b
\mathbf{\hat{j}} + c \mathbf{\hat{k}}`. Similarly, a dyadic can be represented
in tensor form:

.. math::
  \begin{bmatrix}
  a_{11} & a_{12} & a_{13} \\
  a_{21} & a_{22} & a_{23} \\
  a_{31} & a_{32} & a_{33}
  \end{bmatrix}\\

or in dyadic form:

.. math::
  a_{11} \mathbf{\hat{a}_x}\mathbf{\hat{a}_x} +
  a_{12} \mathbf{\hat{a}_x}\mathbf{\hat{a}_y} +
  a_{13} \mathbf{\hat{a}_x}\mathbf{\hat{a}_z} +
  a_{21} \mathbf{\hat{a}_y}\mathbf{\hat{a}_x} +
  a_{22} \mathbf{\hat{a}_y}\mathbf{\hat{a}_y} +
  a_{23} \mathbf{\hat{a}_y}\mathbf{\hat{a}_z} +
  a_{31} \mathbf{\hat{a}_z}\mathbf{\hat{a}_x} +
  a_{32} \mathbf{\hat{a}_z}\mathbf{\hat{a}_y} +
  a_{33} \mathbf{\hat{a}_z}\mathbf{\hat{a}_z}\\

Just as with vectors, the later representation makes it possible to keep track
of which frames the dyadic is defined with respect to. Also, the two
components of each term in the dyadic need not be in the same frame. The
following is valid:

.. math::
  \mathbf{\hat{a}_x} \otimes \mathbf{\hat{b}_y} = \mathbf{\hat{a}_x}
  \mathbf{\hat{b}_y}

Dyadics can also be crossed and dotted with vectors; again, order matters:

.. math::
  \mathbf{\hat{a}_x}\mathbf{\hat{a}_x} \cdot \mathbf{\hat{a}_x} &=
  \mathbf{\hat{a}_x}\\
  \mathbf{\hat{a}_y}\mathbf{\hat{a}_x} \cdot \mathbf{\hat{a}_x} &=
  \mathbf{\hat{a}_y}\\
  \mathbf{\hat{a}_x}\mathbf{\hat{a}_y} \cdot \mathbf{\hat{a}_x} &= 0\\
  \mathbf{\hat{a}_x} \cdot \mathbf{\hat{a}_x}\mathbf{\hat{a}_x} &=
  \mathbf{\hat{a}_x}\\
  \mathbf{\hat{a}_x} \cdot \mathbf{\hat{a}_x}\mathbf{\hat{a}_y} &=
  \mathbf{\hat{a}_y}\\
  \mathbf{\hat{a}_x} \cdot \mathbf{\hat{a}_y}\mathbf{\hat{a}_x} &= 0\\
  \mathbf{\hat{a}_x} \times \mathbf{\hat{a}_y}\mathbf{\hat{a}_x} &=
  \mathbf{\hat{a}_z}\mathbf{\hat{a}_x}\\
  \mathbf{\hat{a}_x} \times \mathbf{\hat{a}_x}\mathbf{\hat{a}_x} &= 0\\
  \mathbf{\hat{a}_y}\mathbf{\hat{a}_x} \times \mathbf{\hat{a}_z} &=
  - \mathbf{\hat{a}_y}\mathbf{\hat{a}_y}\\

One can also take the time derivative of dyadics or express them in different
frames, just like with vectors.

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
:mod:`mechanics`.

One begins by creating the requisite symbols to describe the system. Then
the reference frame is created and the kinematics are done. ::

  >> from sympy import symbols
  >> from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame
  >> from sympy.physics.mechanics import RigidBody, Particle, Point, outer
  >> from symp.physics.mechanics import linear_momentum, angular_momentum
  >> m, M, l1 = symbols('m M l1')
  >> q1d = dynamicsymbols('q1d')
  >> N = ReferenceFrame('N')
  >> O = Point('O')
  >> O.set_vel(N, 0 * N.x)
  >> Ac = O.locatenew('Ac', l1 * N.x)
  >> P = Ac.locatenew('P', l1 * N.x)
  >> a = ReferenceFrame('a')
  >> a.set_ang_vel(N, q1d * N.z)
  >> Ac.v2pt_theory(O, N, a)
  >> P.v2pt_theory(O, N, a)

Finally, the bodies that make up the system are created. In this case the
system consists of a particle Pa and a RigidBody A. ::

  >> Pa = Particle('Pa', P, m)
  >> I = outer(N.z, N.z)
  >> A = RigidBody('A', Ac, a, M, (I, Ac))

Then one can either choose to evaluate the the momenta of individual components
of the system or of the entire system itself. ::

  >> linear_momentum(N,A)
  M*l1*q1d*N.y
  >> angular_momentum(O, N, Pa)
  4*l1**2*m*q1d*N.z
  >> linear_momentum(N, A, Pa)
  (M*l1*q1d + 2*l1*m*q1d)*N.y
  >> angular_momentum(O, N, A, Pa)
  (4*l1**2*m*q1d + q1d)*N.z

It should be noted that the user can determine either momenta in any frame
in :mod:`mechanics` as the user is allowed to specify the reference frame when
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
:mod:`mechanics`.

As was discussed above in the momenta functions, one first creates the system
by going through an identical procedure. ::

  >> from sympy import symbols
  >> from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, outer
  >> from sympy.physics.mechanics import RigidBody, Particle, mechanics_printing
  >> from symp.physics.mechanics import kinetic_energy, potential_energy, Point
  >> mechanics_printing()
  >> m, M, l1, g, h, H = symbols('m M l1 g h H')
  >> omega = dynamicsymbols('omega')
  >> N = ReferenceFrame('N')
  >> O = Point('O')
  >> O.set_vel(N, 0 * N.x)
  >> Ac = O.locatenew('Ac', l1 * N.x)
  >> P = Ac.locatenew('P', l1 * N.x)
  >> a = ReferenceFrame('a')
  >> a.set_ang_vel(N, omega * N.z)
  >> Ac.v2pt_theory(O, N, a)
  >> P.v2pt_theory(O, N, a)
  >> Pa = Particle('Pa', P, m)
  >> I = outer(N.z, N.z)
  >> A = RigidBody('A', Ac, a, M, (I, Ac))

The user can then determine the kinetic energy of any number of entities of the
system: ::

  >> kinetic_energy(N, Pa)
  2*l1**2*m*q1d**2
  >> kinetic_energy(N, Pa, A)
  M*l1**2*q1d**2/2 + 2*l1**2*m*q1d**2 + q1d**2/2

It should be noted that the user can determine either kinetic energy relative
to any frame in :mod:`mechanics` as the user is allowed to specify the
reference frame when calling the function. In other words the user is not
limited to determining just inertial kinetic energy.

For potential energies, the user must first specify the potential energy of
every entity of the system using the :mod:`set_potential_energy` method. The
potential energy of any number of entities comprising the system can then be
determined: ::

  >> Pa.set_potential_energy(m * g * h)
  >> A.set_potential_energy(M * g * H)
  >> potential_energy(A, Pa)
  H*M*g + g*h*m

One can also determine the Lagrangian for this system: ::

  >> Lagrangian(Pa, A)
  -H*M*g + M*l1**2*q1d**2/2 - g*h*m + 2*l1**2*m*q1d**2 + q1d**2/2

Please refer to the docstrings to learn more about each function.
