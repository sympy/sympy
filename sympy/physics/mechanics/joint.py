from sympy import acos
from sympy.physics.vector import cross, Vector, dot, dynamicsymbols
from sympy.physics.mechanics.functions import convert_tuple_to_vector

__all__ = ['Joint', 'PinJoint', 'SlidingJoint']


class Joint(object):
    """Abstract Base class for all specific joints.

    A Joint connects two bodies (a parent and child) by adding different degrees
    of freedom to the child object with respect to the parent object.

    This is the base class for all specific joints and holds all common methods
    acting as an interface for all joints. Custom joints can be created by
    creating a subclass of this class and overriding 'apply_joint()' method to
    add the dynamics of the custom joint.

    Parameters
    ==========

    name : String
        Name of the joint which makes it unique. Should be different from other
        joints.
    parent : Body
        Instance of Body class serving as the parent in the joint.
    child : Body
        Instance of Body class serving as the child in the joint.
    parent_point_pos : 3 Tuple, optional
        Defines the joint's point where the parent will be connected to child.
        3 Tuple defines the values of x, y and z directions w.r.t parent's
        frame. If it is not supplied, center of mass is added as default.
    child_point_pos : 3 Tuple, optional
        Defines the joint's point where the child will be connected to parent.
        3 Tuple defines the values of x, y and z directions w.r.t child's frame.
        If it is not supplied, center of mass is added as default.
    parent_axis : Vector, optional
        Defines the orientation as a vector which must be aligned with child's
        axis before adding joint. If it is not passed, default is x axis in
        parent's frame.
    child_axis : Vector, optional
        Defines the orientation as a vector which must be aligned with parent's
        axis before adding joint. If it is not passed, default is x axis in
        child's frame.
    """

    def __init__(self, name, parent, child, parent_point_pos=None,
                 child_point_pos=None):
        self.name = name
        self.parent = parent
        self.child = child
        self.coordinates = []
        self.speeds = []
        self.kds = []

        if parent_point_pos is None:
            parent_point_pos = (0, 0, 0)  # Parent's Center of mass

        if child_point_pos is None:
            child_point_pos = (0, 0, 0)  # Child's Center of mass

        self._parent_joint_location = convert_tuple_to_vector(
            parent.frame, parent_point_pos)
        self._child_joint_location = convert_tuple_to_vector(
            child.frame, child_point_pos)

        self._locate_joint_point()
        self.apply_joint()

    def _locate_joint_point(self):
        self.parent_joint_point = self.parent.masscenter.locatenew(
            self.name + '_parent_joint',
            self._parent_joint_location)
        self.child_joint_point = self.child.masscenter.locatenew(
            self.name + '_child_joint',
            self._child_joint_location)

    def _align_axes(self, parent_axis, child_axis):
        """Rotates child_frame so that child_axis is aligned to parent_axis."""
        mag1 = parent_axis.magnitude()
        mag2 = child_axis.magnitude()
        angle = acos(dot(parent_axis, child_axis)/(mag1 * mag2))
        axis = cross(child_axis, parent_axis)
        if axis != Vector(0):
            self.child.frame.orient(
                self.parent.frame, 'Axis', [angle, axis])

    def apply_joint(self):
        """To create a custom joint, this method should add degrees of freedom
        to the bodies in subclass of Joint"""
        raise NotImplementedError("To define a custom pydy.Joint, you need to" +
                                  " override apply_joint method in Joint's" +
                                  " subclass.")


class PinJoint(Joint):
    """
    Pin (Revolute) Joint.

    It is defined such that the child body rotates with respect to the parent
    body about the body fixed parent axis through the angle theta. The point of
    joint (pin joint) is defined by two points in each body which coincides
    together.  They are supplied as parent_point_pos and child_point_pos.  In
    addition, the child's reference frame can be arbitrarily rotated a constant
    amount with respect to the parent axis by passing an axis in child body
    which must align with parent axis after the rotation.

    Attributes
    ==========

    name : string
        This is the name of the joint
    parent : Body
        This returns the parent body object of the joint
    child : Body
        This returns the child body object of the joint
    coordinates : list
        This returns the list of coordinates that define the joint. For a
        PinJoint this is the angle between the two body link frames.
    speeds : list
        This returns the list of generalized speeds for the joint. For a pin
        joint this is equivalent to the derivative of the generalized
        coordinate.
    kds : list
        This returns a list of the kinematic differential equations.
    parent_axis : Vector
        This is the axis in the parent body at the joint location that aligns
        with the given axis of the frame in the child body located at the joint.
    child_axis : Vector
        This is the axis in the child body at the joint location that aligns
        with the given axis of the frame in the parent body located at the joint.

    Parameters
    ==========

    name : String
        Name of the joint which makes it unique. Should be different from other
        joints.
    parent : Body
        Instance of Body class serving as the parent in the joint.
    child : Body
        Instance of Body class serving as the child in the joint.
    parent_point_pos : 3 Tuple, optional
        Defines the pin joint's point where the parent will be connected to
        child.  3 Tuple defines the values of x, y and z directions w.r.t
        parent's frame. If it is not supplied, center of mass is added as
        default.
    child_point_pos : 3 Tuple, optional
        Defines the joint's point where the child will be connected to parent.
        3 Tuple defines the values of x, y and z directions w.r.t child's frame.
        If it is not supplied, center of mass is added as default.
    parent_axis : 3 Tuple or string, optional
        Defines the orientation as a vector which must be aligned with child's
        axis before adding joint. If it is not passed, default is x axis in
        parent's frame. Can accept 'X', 'Y', 'Z' or a 3 tuple describing the
        vector.
    child_axis : 3 Tuple or string, optional
        Defines the orientation as a vector which must be aligned with parent's
        axis before adding joint. If it is not passed, default is x axis in
        child's frame. Can accept 'X', 'Y', 'Z' or a 3 tuple describing the
        vector.

    Examples
    ========

    Adding a Pin Joint which connects center of mass of parent to a point
    pointed by child.frame.x + child.frame.y w.r.t. child's center of mass.
    Gravity is along y axis of parent. ::

        >>> from sympy import Symbol
        >>> from sympy.physics.mechanics import Body, PinJoint
        >>> parent = Body('parent')
        >>> child = Body('child')
        >>> gravity = Symbol('gravity')
        >>> l = Symbol('l')
        >>> child.apply_force(child.mass * gravity * parent.frame.y,
        ...                   child.masscenter)
        >>> pin_joint = PinJoint('pinjoint', parent, \
        ...                      child, child_point_pos=(l, l, 0))
        >>> pin_joint.coordinates
        [pinjoint_theta(t)]
        >>> pin_joint.speeds
        [pinjoint_omega(t)]

    """
    def __init__(self, name, parent, child, parent_point_pos=None,
                 child_point_pos=None, parent_axis=None, child_axis=None):

        parent_axes_str = {'X': parent.frame.x, 'Y': parent.frame.y,
                           'Z': parent.frame.z}
        child_axes_str = {'X': child.frame.x, 'Y': child.frame.y,
                          'Z': child.frame.z}

        if parent_axis is None:
            self.parent_axis = parent.frame.x
        elif isinstance(parent_axis, tuple):
            self.parent_axis = convert_tuple_to_vector(self.parent.frame,
                                                       parent_axis)
        elif isinstance(parent_axis, str) and parent_axis.upper() \
                in parent_axes_str.keys():
            self.parent_axis = parent_axes_str[parent_axis.upper()]
        else:
            raise TypeError("Parent Axis must either be one of 'X', 'Y', 'Z'" +
                            " or a 3-Tuple.")

        if child_axis is None:
            self.child_axis = child.frame.x
        elif isinstance(child_axis, tuple):
            self.child_axis = convert_tuple_to_vector(self.child.frame,
                                                      child_axis)
        elif isinstance(child_axis, str) and child_axis.upper() \
                in child_axes_str.keys():
            self.child_axis = child_axes_str[child_axis.upper()]
        else:
            raise TypeError("Child Axis must either be one of 'X', 'Y', 'Z'" +
                            " or a 3 Tuple")

        super(PinJoint, self).__init__(name, parent, child, parent_point_pos,
                                       child_point_pos)

    def apply_joint(self):
        theta = dynamicsymbols(self.name + '_theta')
        thetad = dynamicsymbols(self.name + '_theta', 1)
        omega = dynamicsymbols(self.name + '_omega')

        self.coordinates.append(theta)
        self.speeds.append(omega)
        self.kds.append(thetad - omega)

        self.child.frame.orient(self.parent.frame, 'Axis',
                                [theta, self.parent_axis])
        self.child.frame.set_ang_vel(self.parent.frame,
                                     omega * self.parent_axis)
        self._align_axes(self.parent_axis, self.child_axis)
        self.parent_joint_point.set_vel(self.parent.frame, 0)
        self.child_joint_point.set_vel(self.parent.frame, 0)
        self.child_joint_point.set_pos(self.parent_joint_point, 0)
        self.child.masscenter.v2pt_theory(self.parent.masscenter,
                                          self.parent.frame, self.child.frame)


class SlidingJoint(Joint):
    """
    Sliding (Prismatic) Joint.

    It is defined such that the child body translates with respect to the parent
    body along the body fixed parent axis by the displacement x. The point of
    joint (sliding joint) is defined by two points in each body which coincides
    initially. They are supplied as parent_point_pos and child_point_pos.
    In addition, the child's reference frame can be arbitrarily rotated a
    constant amount with respect to the parent axis by passing an axis in child
    body which must align with parent axis after the rotation.

    Attributes
    ==========

    name : string
        This is the name of the joint
    parent : Body
        This returns the parent body object of the joint
    child : Body
        This returns the child body object of the joint
    coordinates : list
        This returns the list of coordinates that define the joint. For a
        SlidingJoint this is the distance between the two body link points.
    speeds : list
        This returns the list of generalized speeds for the joint. For a sliding
        joint this is equivalent to the derivative of the generalized
        coordinate.
    kds : list
        This returns a list of the kinematic differential equations.
    parent_axis : Vector
        This is the axis in the parent body at the joint location that aligns
        with the given axis of the frame in the child body located at the joint.
    child_axis : Vector
        This is the axis in the child body at the joint location that aligns
        with the given axis of the frame in the parent body located at the joint.

    Parameters
    ==========

    name : String
        Name of the joint which makes it unique. Should be different from other
        joints.
    parent : Body
        Instance of Body class serving as the parent in the joint.
    child : Body
        Instance of Body class serving as the child in the joint.
    parent_point_pos : 3 Tuple, optional
        Defines the joint's point where the parent will be connected to child.
        3 Tuple defines the values of x, y and z directions w.r.t parent's
        frame. If it is not supplied, center of mass is added as default.
    child_point_pos : 3 Tuple, optional
        Defines the joint's point where the child will be connected to parent.
        3 Tuple defines the values of x, y and z directions w.r.t child's frame.
        If it is not supplied, center of mass is added as default.
    parent_axis : 3 Tuple or string, optional
        Defines the orientation as a vector which must be aligned with child's
        axis before adding joint. If it is not passed, default is x axis in
        parent's frame. Can pass in 'X', 'Y', 'Z' or a 3 tuple describing the
        vector.
    child_axis : 3 Tuple or string, optional
        Defines the orientation as a vector which must be aligned with parent's
        axis before adding joint. If it is not passed, default is x axis in
        child's frame. Can pass in 'X', 'Y', 'Z' or a 3 tuple describing the
        vector.

    Examples
    ========

    Adds sliding Joint between parent's masscenter and a point located at unit
    distance in x axis of child. ::

        >>> from sympy import Symbol
        >>> from sympy.physics.mechanics import Body, SlidingJoint
        >>> parent = Body('parent')
        >>> child = Body('child')
        >>> l = Symbol('l')
        >>> sliding_joint = SlidingJoint('slidingjoint', parent, child, \
        ...                             child_point_pos=(l, 0, 0))
        >>> sliding_joint.coordinates
        [slidingjoint_x(t)]
        >>> sliding_joint.speeds
        [slidingjoint_v(t)]

    """
    def __init__(self, name, parent, child, parent_point_pos=None,
                 child_point_pos=None, parent_axis=None, child_axis=None):

        parent_axes_str = {'X': parent.frame.x, 'Y': parent.frame.y,
                           'Z': parent.frame.z}
        child_axes_str = {'X': child.frame.x, 'Y': child.frame.y,
                          'Z': child.frame.z}

        if parent_axis is None:
            self.parent_axis = parent.frame.x
        elif isinstance(parent_axis, tuple):
            self.parent_axis = convert_tuple_to_vector(self.parent.frame,
                                                       parent_axis)
        elif isinstance(parent_axis, str) and parent_axis.upper() \
                in parent_axes_str.keys():
            self.parent_axis = parent_axes_str[parent_axis.upper()]
        else:
            raise TypeError("Parent Axis must either be one of 'X', 'Y', 'Z'" +
                            " or a 3-Tuple.")

        if child_axis is None:
            self.child_axis = child.frame.x
        elif isinstance(child_axis, tuple):
            self.child_axis = convert_tuple_to_vector(self.child.frame,
                                                      child_axis)
        elif isinstance(child_axis, str) and child_axis.upper() \
                in child_axes_str.keys():
            self.child_axis = child_axes_str[child_axis.upper()]
        else:
            raise TypeError("Child Axis must either be one of 'X', 'Y', 'Z'" +
                            " or a 3 Tuple")

        super(SlidingJoint, self).__init__(name, parent, child,
                                           parent_point_pos, child_point_pos)

    def apply_joint(self):
        x = dynamicsymbols(self.name + '_x')
        xd = dynamicsymbols(self.name + '_x', 1)
        v = dynamicsymbols(self.name + '_v')

        self.coordinates.append(x)
        self.speeds.append(v)
        self.kds.append(xd - v)

        self.child.frame.orient(self.parent.frame, 'Axis',
                                [0, self.parent.frame.x])
        self._align_axes(self.parent_axis, self.child_axis)
        self._locate_joint_point()

        self.parent_joint_point.set_vel(self.parent.frame, 0)
        self.child_joint_point.set_vel(self.child.frame, 0)

        self.child_joint_point.set_pos(self.parent_joint_point,
                                       x * self.parent_axis)
        self.child_joint_point.set_vel(self.parent.frame,
                                       v * self.parent_axis)
        self.child.masscenter.set_vel(self.parent.frame,
                                      v * self.parent_axis)
