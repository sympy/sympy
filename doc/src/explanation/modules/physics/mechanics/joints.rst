.. _joints_framework:

=====================================
Joints Framework in Physics/Mechanics
=====================================

:mod:`sympy.physics.mechanics` provides a joints framework. This system consists
of two parts. The first are the :obj:`joints<sympy.physics.mechanics.joint>`
themselves, which are used to create connections between
:class:`bodies<sympy.physics.mechanics.rigidbody.RigidBody>`. The second part is the
:class:`~.System`, which is used to form the equations of motion. Both of these
parts are doing what we can call "book-keeping": keeping track of the
relationships between
:class:`bodies<sympy.physics.mechanics.rigidbody.RigidBody>`.

Joints in Physics/Mechanics
===========================

The general task of the :mod:`joints<sympy.physics.mechanics.joint>` is creating
kinematic relationships between
:class:`bodies<sympy.physics.mechanics.rigidbody.RigidBody>`. A joint is
generally described as shown in the image below.

.. raw:: html
   :file: joint_explanation.svg

As can be seen in this image, each joint needs several objects in order to
define the relationships. First off it needs two bodies: the parent body (shown
in green) and the child body (shown in blue). The transformation made by the
joint is defined between the joint attachments of both bodies. A joint
attachment of a body consists of a point and a body-fixed frame. In the parent
body the point is called ``parent_point`` and the frame ``parent_interframe``.
For the child body these are called ``child_point`` and ``child_interframe``.
For most joints it is the case that when the generalized coordinates are zero,
that there is no rotation or translation between the parent and child joint
attachments. So the ``child_point`` is at the same location as the
``parent_point`` and the ``child_interframe`` is in the same orientation as the
``parent_interframe``.

For describing the joint transformation the joint generally needs
:func:`~.dynamicsymbols` for the generalized coordinates and speeds. Some joints
like the :class:`~.PinJoint`, :class:`~.PrismaticJoint` also require a
``joint_axis``, which consists of the same components in the
``parent_interframe`` and ``child_interframe``. This means that if for example
the joint axis is defined in the ``parent_interframe`` as $2\hat{p}_x +
4\hat{p}_y + 3\hat{p}_z$, then this will also be $2\hat{c}_x + 4\hat{c}_y +
3\hat{c}_z$ in the ``child_interframe``. Practically this means that in the case
of the :class:`~.PinJoint`, also shown below, the ``joint_axis`` is the axis of
rotation, with the generalized coordinate :math:`q` as the angle of
rotation and the generalized speed :math:`u` as the angular velocity.

.. raw:: html
   :file: PinJoint.svg

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
   >>> mechanics_printing(pretty_print=False)
   >>> q, u = dynamicsymbols('q, u')
   >>> parent = RigidBody('parent')
   >>> child = RigidBody('child')
   >>> joint = PinJoint(
   ...     'hinge', parent, child, coordinates=q, speeds=u,
   ...     parent_point=3 * parent.frame.x,
   ...     child_point=-3 * child.frame.x,
   ...     joint_axis=parent.frame.z)
   >>> joint.kdes
   Matrix([[u - q']])
   >>> joint.parent_point.pos_from(parent.masscenter)
   3*parent_frame.x
   >>> joint.parent_interframe
   parent_frame
   >>> joint.joint_axis.express(child.frame)
   child_frame.z
   >>> child.masscenter.pos_from(parent.masscenter)
   3*parent_frame.x + 3*child_frame.x
   >>> child.masscenter.vel(parent.frame)
   3*u*child_frame.y

System in Physics/Mechanics
===========================
After defining the entire system you can use the :class:`~.System` to parse the
system and form the equations of motion. In this process the :class:`~.System`
only does the "book-keeping" of the joints. It uses another method, like the
:class:`~.KanesMethod`, as its backend for forming the equations of motion.

In the code below we form the equations of motion of the single
:class:`~.PinJoint` shown previously. ::

   >>> system = System.from_newtonian(parent)
   >>> system.add_joints(joint)
   >>> system.form_eoms()
   Matrix([[-(child_izz + 9*child_mass)*u']])
   >>> type(system.eom_method)  # The method working in the backend
   <class 'sympy.physics.mechanics.kane.KanesMethod'>
