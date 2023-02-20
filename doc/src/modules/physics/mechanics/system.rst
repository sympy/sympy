.. _system:

===========================
System in Physics/Mechanics
===========================

:mod:`sympy.physics.mechanics` provides a :class:`~.System` class. In a
:class:`~.System` instance you store information about the model, e.g. bodies,
joints, constraints, etc. With all relations of the :class:`~.System` being
defined it can form the equations of motion using a backend of your choosing.
In general :class:`~.System` has been designed to be compatible with third party
libraries, like `PyDy`.

Usage example
=============

:class:`~.System` requires at least an inertial frame and global origin on
initialization. If those are not provided, then they will be created for you. An
instance of :class:`~.System` can be created as follows:

   >>> from sympy.physics.mechanics import *
   >>> mechanics_printing(pretty_print=False)
   >>> N = ReferenceFrame('N')
   >>> O = Point('O')
   >>> system = System(O, N)

Another option is to formulate the :class:`~.System` with respect to a Newtonian
body. This can be done as follows:

   >>> ceiling = RigidBody('ceiling')
   >>> system = System.from_newtonian(ceiling)
   >>> system.origin, system.frame
   (ceiling_masscenter, ceiling_frame)

Next you can start defining the system by for example adding two
:class:`PinJoints<sympy.physics.mechanics.joint.PinJoint>`, making it a double
pendulum:

   >>> from sympy import symbols
   >>> l1, l2 = symbols('l1:3')
   >>> link1 = RigidBody('link1')
   >>> link2 = RigidBody('link2')
   >>> system.add_joints(PinJoint('hinge1', ceiling, link1,
   ...                            child_point=l1 / 2 * link1.y,
   ...                            joint_axis=ceiling.z))
   >>> system.add_joints(PinJoint('hinge2', link1, link2,
   ...                            parent_point=-l1 / 2 * link1.y,
   ...                            child_point=l2 / 2 * link2.y,
   ...                            joint_axis=link1.z))
   >>> system.joints
   (PinJoint: hinge1  parent: ceiling  child: link1, PinJoint: hinge2  parent: link1  child: link2)
   >>> system.q
   Matrix([
   [q_hinge1],
   [q_hinge2]])
   >>> [body.name for body in system.bodies]
   ['ceiling', 'link1', 'link2']

For more information of the joints check out the joints documentation.
Constraints can be added in the following way. The constraint below restricts
the tip of the ``link2`` to only move along the ``y-axis`` of the ``ceiling``
frame.

   >>> r_tip = link2.masscenter.pos_from(system.origin) - l2 / 2 * link2.y
   >>> system.add_holonomic_constraints(r_tip.dot(ceiling.y))
   >>> system.holonomic_constraints
   Matrix([[-l1*cos(q_hinge1) - l2*(-sin(q_hinge1)*sin(q_hinge2) + cos(q_hinge1)*cos(q_hinge2))]])

Note that the corresponding velocity constraint is added by default as well.
However you will still have to choose which generalized coordinates and speeds
become dependent. This can also be seen when using the
:obj:`~.System.validate_system` method:

   >>> try:
   ...     system.validate_system()
   ... except ValueError as e:
   ...     print(e)
   The number of dependent generalized coordinates 0 should be equal to the number of holonomic constraints 1.
   The number of dependent generalized speeds 0 should be equal to the number of velocity constraints 1.

This can be done by manually specifying which should become dependent.

   >>> q_hinge1, q_hinge2 = system.q
   >>> system.q_ind = q_hinge1
   >>> system.q_dep = q_hinge2
   >>> u_hinge1, u_hinge2 = system.u
   >>> system.u_ind = u_hinge1
   >>> system.u_dep = u_hinge2

Now the :obj:`~.System.validate_system` method does not give any errors, we can
form the equations of motion.

   >>> system.validate_system()
   >>> eoms = system.form_eoms()

