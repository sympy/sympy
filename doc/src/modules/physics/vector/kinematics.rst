=====================
Vector: Kinematics
=====================

This document will give some mathematical background
to describing a system's kinematics as well as how to represent the kinematics
in :mod:`physics.vector`.

Introduction to Kinematics
==========================

The first topic is rigid motion kinematics. A rigid body is an idealized
representation of a physical object which has mass and rotational inertia.
Rigid bodies are obviously not flexible. We can break down rigid body motion
into translational motion, and rotational motion (when dealing with particles, we
only have translational motion). Rotational motion can further be broken down
into simple rotations and general rotations.

Translation of a rigid body is defined as a motion where the orientation of the
body does not change during the motion; or during the motion any line segment
would be parallel to itself at the start of the motion.

Simple rotations are rotations in which the orientation of the body may change,
but there is always one line which remains parallel to itself at the start of
the motion.

General rotations are rotations which there is not always one line parallel to
itself at the start of the motion.

Angular Velocity
----------------

The angular velocity of a rigid body refers to the rate of change of its
orientation. The angular velocity of a body is written down as:
:math:`^{\mathbf{N}}\mathbf{\omega}^{\mathbf{B}}`, or the angular velocity of
:math:`\mathbf{B}` in :math:`\mathbf{N}`, which is a vector. Note that here,
the term rigid body was used, but reference frames can also have angular
velocities. Further discussion of the distinction between a rigid body and a
reference frame will occur later when describing the code representation.

Angular velocity is defined as being positive in the direction which causes the
orientation angles to increase (for simple rotations, or series of simple
rotations).

.. image:: kin_angvel1.*
   :height: 350
   :width: 250
   :align: center

The angular velocity vector represents the time derivative of the orientation.
As a time derivative vector quantity, like those covered in the Vector &
ReferenceFrame documentation, this quantity (angular velocity) needs to be
defined in a reference frame. That is what the :math:`\mathbf{N}` is in the
above definition of angular velocity; the frame in which the angular velocity
is defined in.

The angular velocity of :math:`\mathbf{B}` in :math:`\mathbf{N}` can also be
defined by:

.. math::
  ^{\mathbf{N}}\mathbf{\omega}^{\mathbf{B}} =
  (\frac{^{\mathbf{N}}d \mathbf{\hat{b}_y}}{dt}\cdot\mathbf{\hat{b}_z}
  )\mathbf{\hat{b}_x} + (\frac{^{\mathbf{N}}d \mathbf{\hat{b}_z}}{dt}\cdot
  \mathbf{\hat{b}_x})\mathbf{\hat{b}_y} + (\frac{^{\mathbf{N}}d
  \mathbf{\hat{b}_x}}{dt}\cdot\mathbf{\hat{b}_y})\mathbf{\hat{b}_z}

It is also common for a body's angular velocity to be written as:

.. math::
  ^{\mathbf{N}}\mathbf{\omega}^{\mathbf{B}} = w_x \mathbf{\hat{b}_x} +
  w_y \mathbf{\hat{b}_y} + w_z \mathbf{\hat{b}_z}

There are a few additional important points relating to angular velocity. The
first is the addition theorem for angular velocities, a way of relating the
angular velocities of multiple bodies and frames. The theorem follows:

.. math::
  ^{\mathbf{N}}\mathbf{\omega}^{\mathbf{D}} =
  ^{\mathbf{N}}\mathbf{\omega}^{\mathbf{A}} +
  ^{\mathbf{A}}\mathbf{\omega}^{\mathbf{B}} +
  ^{\mathbf{B}}\mathbf{\omega}^{\mathbf{C}} +
  ^{\mathbf{C}}\mathbf{\omega}^{\mathbf{D}}

This is also shown in the following example:

.. image:: kin_angvel2.*
   :height: 300
   :width: 450
   :align: center

.. math::
  ^{\mathbf{N}}\mathbf{\omega}^{\mathbf{A}} &= 0\\
  ^{\mathbf{A}}\mathbf{\omega}^{\mathbf{B}} &= \dot{q_1} \mathbf{\hat{a}_x}\\
  ^{\mathbf{B}}\mathbf{\omega}^{\mathbf{C}} &= - \dot{q_2} \mathbf{\hat{b}_z}\\
  ^{\mathbf{C}}\mathbf{\omega}^{\mathbf{D}} &= \dot{q_3} \mathbf{\hat{c}_y}\\
  ^{\mathbf{N}}\mathbf{\omega}^{\mathbf{D}} &= \dot{q_1} \mathbf{\hat{a}_x}
  - \dot{q_2} \mathbf{\hat{b}_z} + \dot{q_3} \mathbf{\hat{c}_y}\\

Note the signs used in the angular velocity definitions, which are related to
how the displacement angle is defined in this case.


This theorem makes defining angular velocities of multibody systems much
easier, as the angular velocity of a body in a chain needs to only be defined
to the previous body in order to be fully defined (and the first body needs
to be defined in the desired reference frame). The following figure shows an
example of when using this theorem can make things easier.

.. image:: kin_angvel3.*
   :height: 250
   :width: 400
   :align: center

Here we can easily write the angular velocity of the body
:math:`\mathbf{D}` in the reference frame of the first body :math:`\mathbf{A}`:

.. math::
  ^\mathbf{A}\mathbf{\omega}^\mathbf{D} = w_1 \mathbf{\hat{p_1}} +
  w_2 \mathbf{\hat{p_2}} + w_3 \mathbf{\hat{p_3}}\\

It is very important to remember to only use this with angular velocities; you
cannot use this theorem with the velocities of points.

There is another theorem commonly used: the derivative theorem. It provides an
alternative method (which can be easier) to calculate the time derivative of a
vector in a reference frame:

.. math::
  \frac{^{\mathbf{N}} d \mathbf{v}}{dt} = \frac{^{\mathbf{B}} d \mathbf{v}}{dt}
  + ^{\mathbf{N}}\mathbf{\omega}^{\mathbf{B}} \times \mathbf{v}

The vector :math:`\mathbf{v}` can be any vector quantity: a position vector,
a velocity vector, angular velocity vector, etc. Instead of taking the time
derivative of the vector in :math:`\mathbf{N}`, we take it in
:math:`\mathbf{B}`, where :math:`\mathbf{B}` can be any reference frame or
body, usually one in which it is easy to take the derivative on
:math:`\mathbf{v}` in (:math:`\mathbf{v}` is usually composed only of the basis
vector set belonging to :math:`\mathbf{B}`). Then we add the cross product of
the angular velocity of our newer frame,
:math:`^{\mathbf{N}}\mathbf{\omega}^{\mathbf{B}}` and our vector quantity
:math:`\mathbf{v}`. Again, you can choose any alternative frame for this.
Examples follow:

.. % need multiple examples here showing the derivative theorem


Angular Acceleration
--------------------
Angular acceleration refers to the time rate of change of the angular velocity
vector. Just as the angular velocity vector is for a body and is specified in a
frame, the angular acceleration vector is for a body and is specified in a
frame: :math:`^{\mathbf{N}}\mathbf{\alpha}^{\mathbf{B}}`, or the angular
acceleration of :math:`\mathbf{B}` in :math:`\mathbf{N}`, which is a vector.

Calculating the angular acceleration is relatively straight forward:

.. math::
  ^{\mathbf{N}}\mathbf{\alpha}^{\mathbf{B}} =
  \frac{^{\mathbf{N}} d ^{\mathbf{N}}\mathbf{\omega}^{\mathbf{B}}}{dt}

Note that this can be calculated with the derivative theorem, and when the
angular velocity is defined in a body fixed frame, becomes quite simple:

.. math::

  ^{\mathbf{N}}\mathbf{\alpha}^{\mathbf{B}} &=
  \frac{^{\mathbf{N}} d ^{\mathbf{N}}\mathbf{\omega}^{\mathbf{B}}}{dt}\\

  ^{\mathbf{N}}\mathbf{\alpha}^{\mathbf{B}} &=
  \frac{^{\mathbf{B}} d ^{\mathbf{N}}\mathbf{\omega}^{\mathbf{B}}}{dt}
  + ^{\mathbf{N}}\mathbf{\omega}^{\mathbf{B}} \times
  ^{\mathbf{N}}\mathbf{\omega}^{\mathbf{B}}\\

  \textrm{if } ^{\mathbf{N}}\mathbf{\omega}^{\mathbf{B}} &=
  w_x \mathbf{\hat{b}_x} + w_y \mathbf{\hat{b}_y} + w_z \mathbf{\hat{b}_z}\\

  \textrm{then } ^{\mathbf{N}}\mathbf{\alpha}^{\mathbf{B}} &=
  \frac{^{\mathbf{B}} d ^{\mathbf{N}}\mathbf{\omega}^{\mathbf{B}}}{dt}
  + \underbrace{^{\mathbf{N}}\mathbf{\omega}^{\mathbf{B}} \times
  ^{\mathbf{N}}\mathbf{\omega}^{\mathbf{B}}}_{
  \textrm{this is 0 by definition}}\\

  ^{\mathbf{N}}\mathbf{\alpha}^{\mathbf{B}}&=\frac{d w_x}{dt}\mathbf{\hat{b}_x}
  + \frac{d w_y}{dt}\mathbf{\hat{b}_y} + \frac{d w_z}{dt}\mathbf{\hat{b}_z}\\

  ^{\mathbf{N}}\mathbf{\alpha}^{\mathbf{B}}&= \dot{w_x}\mathbf{\hat{b}_x} +
  \dot{w_y}\mathbf{\hat{b}_y} + \dot{w_z}\mathbf{\hat{b}_z}\\

Again, this is only for the case in which the angular velocity of the body is
defined in body fixed components.



Point Velocity & Acceleration
-----------------------------

Consider a point, :math:`P`: we can define some characteristics of the point.
First, we can define a position vector from some other point to :math:`P`.
Second, we can define the velocity vector of :math:`P` in a reference frame of
our choice. Third, we can define the acceleration vector of :math:`P` in a
reference frame of our choice.

These three quantities are read as:

.. math::
  \mathbf{r}^{OP} \textrm{, the position vector from } O
  \textrm{ to }P\\
  ^{\mathbf{N}}\mathbf{v}^P \textrm{, the velocity of } P
  \textrm{ in the reference frame } \mathbf{N}\\
  ^{\mathbf{N}}\mathbf{a}^P \textrm{, the acceleration of } P
  \textrm{ in the reference frame } \mathbf{N}\\

Note that the position vector does not have a frame associated with it; this is
because there is no time derivative involved, unlike the velocity and
acceleration vectors.

We can find these quantities for a simple example easily:

.. image:: kin_1.*
   :height: 300
   :width: 300
   :align: center

.. math::
  \textrm{Let's define: }
  \mathbf{r}^{OP} &= q_x \mathbf{\hat{n}_x} + q_y \mathbf{\hat{n}_y}\\
  ^{\mathbf{N}}\mathbf{v}^P &= \frac{^{\mathbf{N}} d \mathbf{r}^{OP}}{dt}\\
  \textrm{then we can calculate: }
  ^{\mathbf{N}}\mathbf{v}^P &= \dot{q}_x\mathbf{\hat{n}_x} +
  \dot{q}_y\mathbf{\hat{n}_y}\\
  \textrm{and :}
  ^{\mathbf{N}}\mathbf{a}^P &= \frac{^{\mathbf{N}} d
  ^{\mathbf{N}}\mathbf{v}^P}{dt}\\
  ^{\mathbf{N}}\mathbf{a}^P &= \ddot{q}_x\mathbf{\hat{n}_x} +
  \ddot{q}_y\mathbf{\hat{n}_y}\\

It is critical to understand in the above example that the point :math:`O` is
fixed in the reference frame :math:`\mathbf{N}`. There is no addition theorem
for translational velocities; alternatives will be discussed later though.
Also note that the position of every point might not
always need to be defined to form the dynamic equations of motion.
When you don't want to define the position vector of a point, you can start by
just defining the velocity vector. For the above example:

.. math::
  \textrm{Let us instead define the velocity vector as: }
  ^{\mathbf{N}}\mathbf{v}^P &= u_x \mathbf{\hat{n}_x} +
  u_y \mathbf{\hat{n}_y}\\
  \textrm{then acceleration can be written as: }
  ^{\mathbf{N}}\mathbf{a}^P &= \dot{u}_x \mathbf{\hat{n}_x} +
  \dot{u}_y \mathbf{\hat{n}_y}\\


There will often be cases when the velocity of a point is desired and a related
point's velocity is known. For the cases in which we have two points fixed on a
rigid body, we use the 2-Point Theorem:

.. image:: kin_2pt.*
   :height: 300
   :width: 300
   :align: center

Let's say we know the velocity of the point :math:`S` and the angular
velocity of the body :math:`\mathbf{B}`, both defined in the reference frame
:math:`\mathbf{N}`. We can calculate the velocity and acceleration
of the point :math:`P` in :math:`\mathbf{N}` as follows:

.. math::
  ^{\mathbf{N}}\mathbf{v}^P &= ^\mathbf{N}\mathbf{v}^S +
  ^\mathbf{N}\mathbf{\omega}^\mathbf{B} \times \mathbf{r}^{SP}\\
  ^{\mathbf{N}}\mathbf{a}^P &= ^\mathbf{N}\mathbf{a}^S +
  ^\mathbf{N}\mathbf{\alpha}^\mathbf{B} \times \mathbf{r}^{SP} +
  ^\mathbf{N}\mathbf{\omega}^\mathbf{B} \times
  (^\mathbf{N}\mathbf{\omega}^\mathbf{B} \times \mathbf{r}^{SP})\\

When only one of the two points is fixed on a body, the 1 point theorem is used
instead.

.. image:: kin_1pt.*
   :height: 400
   :width: 400
   :align: center

Here, the velocity of point :math:`S` is known in the frame :math:`\mathbf{N}`,
the angular velocity of :math:`\mathbf{B}` is known in :math:`\mathbf{N}`, and
the velocity of the point :math:`P` is known in the frame associated with body
:math:`\mathbf{B}`. We can then write the velocity and acceleration of
:math:`P` in :math:`\mathbf{N}` as:

.. math::
  ^{\mathbf{N}}\mathbf{v}^P &= ^\mathbf{B}\mathbf{v}^P +
  ^\mathbf{N}\mathbf{v}^S + ^\mathbf{N}\mathbf{\omega}^\mathbf{B} \times
  \mathbf{r}^{SP}\\

  ^{\mathbf{N}}\mathbf{a}^P &= ^\mathbf{B}\mathbf{a}^S +
  ^\mathbf{N}\mathbf{a}^O + ^\mathbf{N}\mathbf{\alpha}^\mathbf{B}
  \times \mathbf{r}^{SP} + ^\mathbf{N}\mathbf{\omega}^\mathbf{B} \times
  (^\mathbf{N}\mathbf{\omega}^\mathbf{B} \times \mathbf{r}^{SP}) +
  2 ^\mathbf{N}\mathbf{\omega}^\mathbf{B} \times ^\mathbf{B} \mathbf{v}^P \\


Examples of applications of the 1 point and 2 point theorem follow.

.. image:: kin_2.*
   :height: 300
   :width: 400
   :align: center

This example has a disc translating and rotating in a plane. We can easily
define the angular velocity of the body :math:`\mathbf{B}` and velocity of the
point :math:`O`:

.. math::
  ^\mathbf{N}\mathbf{\omega}^\mathbf{B} &= u_3 \mathbf{\hat{n}_z} = u_3
  \mathbf{\hat{b}_z}\\
  ^\mathbf{N}\mathbf{v}^O &= u_1 \mathbf{\hat{n}_x} + u_2 \mathbf{\hat{n}_y}\\

and accelerations can be written as:

.. math::
  ^\mathbf{N}\mathbf{\alpha}^\mathbf{B} &= \dot{u_3} \mathbf{\hat{n}_z} =
  \dot{u_3} \mathbf{\hat{b}_z}\\
  ^\mathbf{N}\mathbf{a}^O &= \dot{u_1} \mathbf{\hat{n}_x} + \dot{u_2}
  \mathbf{\hat{n}_y}\\

We can use the 2 point theorem to calculate the velocity and acceleration of
point :math:`P` now.

.. math::
  \mathbf{r}^{OP} &= R \mathbf{\hat{b}_x}\\
  ^\mathbf{N}\mathbf{v}^P &= ^\mathbf{N}\mathbf{v}^O +
  ^\mathbf{N}\mathbf{\omega}^\mathbf{B} \times \mathbf{r}^{OP}\\
  ^\mathbf{N}\mathbf{v}^P &= u_1 \mathbf{\hat{n}_x} + u_2 \mathbf{\hat{n}_y}
  + u_3 \mathbf{\hat{b}_z} \times R \mathbf{\hat{b}_x} = u_1
  \mathbf{\hat{n}_x} + u_2 \mathbf{\hat{n}_y} + u_3 R \mathbf{\hat{b}_y}\\
  ^{\mathbf{N}}\mathbf{a}^P &= ^\mathbf{N}\mathbf{a}^O +
  ^\mathbf{N}\mathbf{\alpha}^\mathbf{B} \times \mathbf{r}^{OP} +
  ^\mathbf{N}\mathbf{\omega}^\mathbf{B} \times
  (^\mathbf{N}\mathbf{\omega}^\mathbf{B} \times \mathbf{r}^{OP})\\
  ^{\mathbf{N}}\mathbf{a}^P &= \dot{u_1} \mathbf{\hat{n}_x} + \dot{u_2}
  \mathbf{\hat{n}_y} + \dot{u_3}\mathbf{\hat{b}_z}\times R \mathbf{\hat{b}_x}
  +u_3\mathbf{\hat{b}_z}\times(u_3\mathbf{\hat{b}_z}\times
  R\mathbf{\hat{b}_x})\\
  ^{\mathbf{N}}\mathbf{a}^P &= \dot{u_1} \mathbf{\hat{n}_x} + \dot{u_2}
  \mathbf{\hat{n}_y} + R\dot{u_3}\mathbf{\hat{b}_y} - R u_3^2
  \mathbf{\hat{b}_x}\\

.. image:: kin_3.*
   :height: 200
   :width: 200
   :align: center


In this example we have a double pendulum. We can use the two point theorem
twice here in order to find the velocity of points :math:`Q` and :math:`P`;
point :math:`O`'s velocity is zero in :math:`\mathbf{N}`.

.. math::
  \mathbf{r}^{OQ} &= l \mathbf{\hat{b}_x}\\
  \mathbf{r}^{QP} &= l \mathbf{\hat{c}_x}\\
  ^\mathbf{N}\mathbf{\omega}^\mathbf{B} &= u_1 \mathbf{\hat{b}_z}\\
  ^\mathbf{N}\mathbf{\omega}^\mathbf{C} &= u_2 \mathbf{\hat{c}_z}\\
  ^\mathbf{N}\mathbf{v}^Q &= ^\mathbf{N}\mathbf{v}^O +
  ^\mathbf{N}\mathbf{\omega}^\mathbf{B} \times \mathbf{r}^{OQ}\\
  ^\mathbf{N}\mathbf{v}^Q &= u_1 l \mathbf{\hat{b}_y}\\
  ^\mathbf{N}\mathbf{v}^P &= ^\mathbf{N}\mathbf{v}^Q +
  ^\mathbf{N}\mathbf{\omega}^\mathbf{C} \times \mathbf{r}^{QP}\\
  ^\mathbf{N}\mathbf{v}^Q &= u_1 l \mathbf{\hat{b}_y} +u_2 \mathbf{\hat{c}_z}
  \times l \mathbf{\hat{c}_x}\\
  ^\mathbf{N}\mathbf{v}^Q &= u_1 l\mathbf{\hat{b}_y}+u_2 l\mathbf{\hat{c}_y}\\

.. image:: kin_4.*
   :height: 400
   :width: 300
   :align: center

In this example we have a particle moving on a ring; the ring is supported by a
rod which can rotate about the :math:`\mathbf{\hat{n}_x}` axis. First we use
the two point theorem to find the velocity of the center point of the ring,
:math:`Q`, then use the 1 point theorem to find the velocity of the particle on
the ring.

.. math::
  ^\mathbf{N}\mathbf{\omega}^\mathbf{C} &= u_1 \mathbf{\hat{n}_x}\\
  \mathbf{r}^{OQ} &= -l \mathbf{\hat{c}_z}\\
  ^\mathbf{N}\mathbf{v}^Q &= u_1 l \mathbf{\hat{c}_y}\\
  \mathbf{r}^{QP} &= R(cos(q_2) \mathbf{\hat{c}_x}
  + sin(q_2) \mathbf{\hat{c}_y} )\\
  ^\mathbf{C}\mathbf{v}^P &= R u_2 (-sin(q_2) \mathbf{\hat{c}_x}
  + cos(q_2) \mathbf{\hat{c}_y} )\\
  ^\mathbf{N}\mathbf{v}^P &= ^\mathbf{C}\mathbf{v}^P +^\mathbf{N}\mathbf{v}^Q
  + ^\mathbf{N}\mathbf{\omega}^\mathbf{C} \times \mathbf{r}^{QP}\\
  ^\mathbf{N}\mathbf{v}^P &= R u_2 (-sin(q_2) \mathbf{\hat{c}_x}
  + cos(q_2) \mathbf{\hat{c}_y} ) + u_1 l \mathbf{\hat{c}_y} +
  u_1 \mathbf{\hat{c}_x} \times R(cos(q_2) \mathbf{\hat{c}_x}
  + sin(q_2) \mathbf{\hat{c}_y}\\
  ^\mathbf{N}\mathbf{v}^P &= - R u_2 sin(q_2) \mathbf{\hat{c}_x}
  + (R u_2 cos(q_2)+u_1 l)\mathbf{\hat{c}_y} + R u_1 sin(q_2)
  \mathbf{\hat{c}_z}\\

A final topic in the description of velocities of points is that of rolling, or
rather, rolling without slip. Two bodies are said to be rolling without slip if
and only if the point of contact on each body has the same velocity in another
frame. See the following figure:

.. image:: kin_rolling.*
   :height: 250
   :width: 450
   :align: center

This is commonly used to form the velocity of a point on one object rolling on
another fixed object, such as in the following example:

.. % rolling disc kinematics here


Kinematics in physics.vector
============================

It should be clear by now that the topic of kinematics here has been mostly
describing the correct way to manipulate vectors into representing the
velocities of points. Within :mod:`vector` there are convenient methods for
storing these velocities associated with frames and points. We'll now revisit
the above examples and show how to represent them in :mod:`sympy`.

The topic of reference frame creation has already been covered. When a
``ReferenceFrame`` is created though, it automatically calculates the angular
velocity of the frame using the time derivative of the DCM and the angular
velocity definition. ::

  >>> from sympy import Symbol, sin, cos
  >>> from sympy.physics.vector import *
  >>> N = ReferenceFrame('N')
  >>> q1 = dynamicsymbols('q1')
  >>> A = N.orientnew('A', 'Axis', [q1, N.x])
  >>> A.ang_vel_in(N)
  q1'*N.x

Note that the angular velocity can be defined in an alternate way: ::

  >>> B = ReferenceFrame('B')
  >>> u1 = dynamicsymbols('u1')
  >>> B.set_ang_vel(N, u1 * B.y)
  >>> B.ang_vel_in(N)
  u1*B.y
  >>> N.ang_vel_in(B)
  - u1*B.y

Both upon frame creation during ``orientnew`` and when calling ``set_ang_vel``,
the angular velocity is set in both frames involved, as seen above.

.. image:: kin_angvel2.*
   :height: 300
   :width: 450
   :align: center

Here we have multiple bodies with angular velocities defined relative to each
other. This is coded as: ::

  >>> N = ReferenceFrame('N')
  >>> A = ReferenceFrame('A')
  >>> B = ReferenceFrame('B')
  >>> C = ReferenceFrame('C')
  >>> D = ReferenceFrame('D')
  >>> u1, u2, u3 = dynamicsymbols('u1 u2 u3')
  >>> A.set_ang_vel(N, 0)
  >>> B.set_ang_vel(A, u1 * A.x)
  >>> C.set_ang_vel(B, -u2 * B.z)
  >>> D.set_ang_vel(C, u3 * C.y)
  >>> D.ang_vel_in(N)
  u1*A.x - u2*B.z + u3*C.y

In :mod:`vector` the shortest path between two frames is used when finding
the angular velocity. That would mean if we went back and set: ::

  >>> D.set_ang_vel(N, 0)
  >>> D.ang_vel_in(N)
  0

The path that was just defined is what is used.
This can cause problems though, as now the angular
velocity definitions are inconsistent. It is recommended that you avoid
doing this.

.. % put some stuff to go with derivative theorem here

Points are a translational analog to the rotational ``ReferenceFrame``.
Creating a ``Point`` can be done in two ways, like ``ReferenceFrame``: ::

  >>> O = Point('O')
  >>> P = O.locatenew('P', 3 * N.x + N.y)
  >>> P.pos_from(O)
  3*N.x + N.y
  >>> Q = Point('Q')
  >>> Q.set_pos(P, N.z)
  >>> Q.pos_from(P)
  N.z
  >>> Q.pos_from(O)
  3*N.x + N.y + N.z

Similar to ``ReferenceFrame``, the position vector between two points is found
by the shortest path (number of intermediate points) between them. Unlike
rotational motion, there is no addition theorem for the velocity of points. In
order to have the velocity of a ``Point`` in a ``ReferenceFrame``, you have to
set the value. ::

  >>> O = Point('O')
  >>> O.set_vel(N, u1*N.x)
  >>> O.vel(N)
  u1*N.x

For both translational and rotational accelerations, the value is computed by
taking the time derivative of the appropriate velocity, unless the user sets it
otherwise.

  >>> O.acc(N)
  u1'*N.x
  >>> O.set_acc(N, u2*u1*N.y)
  >>> O.acc(N)
  u1*u2*N.y


Next is a description of the 2 point and 1 point theorems, as used in
``sympy``.

.. image:: kin_2.*
   :height: 300
   :width: 400
   :align: center

First is the translating, rotating disc. ::

  >>> N = ReferenceFrame('N')
  >>> u1, u2, u3 = dynamicsymbols('u1 u2 u3')
  >>> R = Symbol('R')
  >>> B = ReferenceFrame('B')
  >>> O = Point('O')
  >>> O.set_vel(N, u1 * N.x + u2 * N.y)
  >>> P = O.locatenew('P', R * B.x)
  >>> B.set_ang_vel(N, u3 * B.z)
  >>> P.v2pt_theory(O, N, B)
  u1*N.x + u2*N.y + R*u3*B.y
  >>> P.a2pt_theory(O, N, B)
  u1'*N.x + u2'*N.y - R*u3**2*B.x + R*u3'*B.y

We will also cover implementation of the 1 point theorem.

.. image:: kin_4.*
   :height: 400
   :width: 300
   :align: center

This is the particle moving on a ring, again. ::

  >>> N = ReferenceFrame('N')
  >>> u1, u2 = dynamicsymbols('u1 u2')
  >>> q1, q2 = dynamicsymbols('q1 q2')
  >>> l = Symbol('l')
  >>> R = Symbol('R')
  >>> C = N.orientnew('C', 'Axis', [q1, N.x])
  >>> C.set_ang_vel(N, u1 * N.x)
  >>> O = Point('O')
  >>> O.set_vel(N, 0)
  >>> Q = O.locatenew('Q', -l * C.z)
  >>> P = Q.locatenew('P', R * (cos(q2) * C.x + sin(q2) * C.y))
  >>> P.set_vel(C, R * u2 * (-sin(q2) * C.x + cos(q2) * C.y))
  >>> Q.v2pt_theory(O, N, C)
  l*u1*C.y
  >>> P.v1pt_theory(Q, N, C)
  - R*u2*sin(q2)*C.x + (R*u2*cos(q2) + l*u1)*C.y + R*u1*sin(q2)*C.z
