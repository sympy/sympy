.. _joints_framework:

=====================================
Joints Framework in Physics/Mechanics
=====================================

:mod:`sympy.physics.mechanics` provides a joints framework. This system consists
of two parts. The first are the :obj:`joints<sympy.physics.mechanics.joint>`
themselves, which are used to create connections between
:class:`bodies<sympy.physics.mechanics.body.Body>`. The second part is the
:class:`~.JointsMethod`, which is used to form the equations of motion. Both of
these parts are doing what we can call "book-keeping": keeping track of the
relationships between :class:`bodies<sympy.physics.mechanics.body.Body>`.

Joints in Physics/Mechanics
===========================

The general task of the :mod:`joints<sympy.physics.mechanics.joint>` is creating
kinematic relationships between
:class:`bodies<sympy.physics.mechanics.body.Body>`. A joint is generally
described as shown in the image below.

.. image:: api/joint_explanation.svg
   :align: center
   :width: 600

As can be seen in this image, each joint needs several objects in order to
define the relationships. First off it needs two bodies: the parent body (shown
in green) and the child body (shown in blue). The attachment of the joint to the
body consists of a point and a frame. In the parent body the point is called
``parent_point`` and the frame ``parent_interframe``. For the child body these
are called ``child_point`` and ``child_interframe``. The arguments for of the
attachment points are ``parent_joint_pos`` and ``child_joint_pos``.

For describing the joint transformation the joint generally needs
:func:`~.dynamicsymbols` for the generalized coordinates and speeds. Some joints
like the :class:`~.PinJoint`, :class:`~.PrismaticJoint` also require a
``joint_axis``, which is an axis that does not change between the
``parent_interframe`` and ``child_interframe``. So in case of the
:class:`~.PinJoint`, also shown below, this means that the ``joint_axis`` is the
axis of rotation. With the generalized coordinate :math:`\theta` as the angle of
rotation and the generalized speed :math:`\omega` as the angular velocity.

.. image:: api/PinJoint.svg
   :align: center
   :width: 600

With the information listed above, the joint defines the following
relationships. It first defines the kinematic differential equations, which
relate the generalized coordinates to the generalized speeds. Next, it orients
the parent and child body with respect to each other. After which it also
defines their velocity relationships.

The code below shows the creation of a :class:`~.PinJoint` as shown above
with arbitrary linked position vectors. In this code the attachment points are
set using vectors, which define the attachment point with respect to the body's
mass center. The intermediate frames are not set, so those are the same as the
body's frame. ::

   >>> from sympy.physics.mechanics import *
   >>> theta, omega = dynamicsymbols('theta, omega')
   >>> parent = Body('parent')
   >>> child = Body('child')
   >>> joint = PinJoint(
   ...     'hinge', parent, child, theta, omega,
   ...     parent_joint_pos=3 * parent.frame.x,
   ...     child_joint_pos=-3 * child.frame.x,
   ...     joint_axis=parent.frame.z)
   >>> joint.kdes
   [omega - theta']
   >>> joint.parent_point.pos_from(parent.masscenter)
   3*parent_frame.x
   >>> joint.parent_interframe
   parent_frame
   >>> joint.joint_axis.express(child.frame)
   child_frame.z
   >>> child.masscenter.pos_from(parent.masscenter)
   3*parent_frame.x + 3*child_frame.x
   >>> child.masscenter.vel(parent.frame)
   3*omega*child_frame.y

JointsMethod in Physics/Mechanics
=================================
After defining the entire system you can use the :class:`~.JointsMethod` to
parse the system and form the equations of motion. In this process the
:class:`~.JointsMethod` only does the "book-keeping" of the joints. It uses
another method, like the :class:`~.KanesMethod`, as its backend for forming the
equations of motion.

In the code below we form the equations of motion of the single
:class:`~.PinJoint` shown previously. ::

   >>> method = JointsMethod(parent, joint)
   >>> method.form_eoms()
   Matrix([[-(child_izz + 9*child_mass)*omega']])
   >>> type(method.method)  # The method working in the backend
   <class 'sympy.physics.mechanics.kane.KanesMethod'>
