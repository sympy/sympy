from sympy import acos, cos, Matrix, sin
from sympy.physics.vector import cross, Vector, dot, dynamicsymbols
from sympy.physics.vector.spatial import rot, xlt
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

    def __init__(self, name, parent, child, coordinates, speeds,
                 parent_point_pos=(0, 0, 0), child_point_pos=(0, 0, 0),
                 parent_axis=None, child_axis=None):
        self.name = name
        self.parent = parent
        self.parent.child_joints.append(self)
        self.child = child
        self.child.parent_joint = self
        self.coordinates = Matrix([i for i in coordinates])
        self.speeds = Matrix([i for i in speeds])
        self.kin_diff = []

        # Process the parent and child joint access input possibilities
        parent_axes_str = {'X': parent.frame.x, 'Y': parent.frame.y,
                           'Z': parent.frame.z}
        child_axes_str = {'X': child.frame.x, 'Y': child.frame.y,
                          'Z': child.frame.z}

        if parent_axis is None:
            self.parent_axis = parent.frame.z
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
            self.child_axis = child.frame.z
        elif isinstance(child_axis, tuple):
            self.child_axis = convert_tuple_to_vector(self.child.frame,
                                                      child_axis)
        elif isinstance(child_axis, str) and child_axis.upper() \
                in child_axes_str.keys():
            self.child_axis = child_axes_str[child_axis.upper()]
        else:
            raise TypeError("Child Axis must either be one of 'X', 'Y', 'Z'" +
                            " or a 3 Tuple.")

        # Orient the frame at the joint in the parent body such that the joint
        # axis is the z axis of the new frame
        mag = self.parent_axis.magnitude()
        angle = acos(dot(self.parent.frame.z, self.parent_axis)/(mag))
        axis = cross(self.parent.frame.z, self.parent_axis)
        if axis == 0:
            axis = self.parent.frame.z
        temp = self.parent.frame.orientnew(name + "_parent_joint_frame", 'Axis',
                                           [angle, axis])
        self.parent_joint_frame = temp
        self.parent_joint_frame.set_ang_vel(self.parent.frame, 0)

        # Create the Point objects for the parent and child bodies at the joint
        # location
        parent_joint_location = convert_tuple_to_vector(parent.frame,
                                                        parent_point_pos)
        self.parent_joint_point = self.parent.masscenter.locatenew(
            self.name + '_parent_joint_point', parent_joint_location)
        self.parent_joint_point.set_vel(self.parent.frame, 0)
        self.parent_joint_point.set_vel(self.parent_joint_frame, 0)

        child_joint_location = convert_tuple_to_vector(child.frame,
                                                       child_point_pos)
        self.child_joint_point = self.child.masscenter.locatenew(
            self.name + '_child_joint_point', child_joint_location)
        self.child_joint_point.set_vel(child.frame, 0)

        # Run the joint specific code
        self.apply_joint()

        # Orient the frame at the joint in the child body such that the joint
        # axis is the z axis of the new frame. Note this has to be done after
        # apply_joint has been run due to each reference frame only being able
        # to be oriented once.
        mag = self.child_axis.magnitude()
        angle = -1 * acos(dot(self.child.frame.z, self.child_axis)/(mag))
        axis = cross(self.child.frame.z, self.child_axis)
        if axis == 0:
            axis = self.child.frame.z

        # Get child_axis defined in the child_joint_frame
        temp_frame = self.child.frame.orientnew("temp_frame", "Axis",
                                                [angle, axis])
        temp_axis = temp_frame.dcm(self.child.frame) * \
                    axis.to_matrix(self.child.frame)
        axis = temp_axis[0]*self.child_joint_frame.x + \
               temp_axis[1]*self.child_joint_frame.y + \
               temp_axis[2]*self.child_joint_frame.z


        self.child.frame.orient(self.child_joint_frame, 'Axis', [angle, axis])

        self.child_joint_frame.set_ang_vel(self.child.frame, 0)
        self.child_joint_point.set_vel(self.child_joint_frame, 0)
        self.child.masscenter.v2pt_theory(self.child_joint_point,
                                          self.parent_joint_frame,
                                          self.child_joint_frame)

    def apply_joint(self):
        """To create a custom joint, this method should add degrees of freedom
        to the bodies in subclass of Joint"""
        raise NotImplementedError("To define a custom pydy.Joint, you need to" +
                                  " override apply_joint method in Joint's" +
                                  " subclass.")

    def spatial_info(self):
        try:
            spat_velocity = self.motion_subspace * self.coordinates.transpose()
            return self.joint_transform, self.motion_subspace, spat_velocity
        except AttributeError:
            raise NotImplementedError("To use this method you need to define" +
                                      " motion subspace and joint transform" +
                                      " for your joint.")

    def XT_child(self):
        if self.parent.parent_joint is None:
            rotation = rot(self.parent_joint_frame.dcm(self.parent.frame))
            vector = self.parent_joint_point.pos_from(self.parent.masscenter)
            translation = xlt(vector.to_matrix(self.parent.frame))
        else:
            temp_frame = self.parent.parent_joint.child_joint_frame
            dcm1 = self.parent.frame.dcm(temp_frame)
            rotation1 = rot(dcm1)
            rotation2 = rot(self.parent_joint_frame.dcm(self.parent.frame))

            temp_point = self.parent.parent_joint.child_joint_point
            vector1 = self.parent.masscenter.pos_from(temp_point)
            vector2 = self.parent_joint_point.pos_from(self.parent.masscenter)
            translation1 = xlt(vector1.to_matrix(temp_frame))
            translation2 = xlt(vector2.to_matrix(self.parent.frame))

            rotation = rotation1*translation1
            translation = rotation2 * translation2

        return rotation*translation


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
    kin_diff : list
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
    coord : function of time
        This is the generalized coordinate to associate with the PinJoint.
    speed : function of time
        This is the generalized speed to associate with the PinJoint.
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
        axis before adding joint. If it is not passed, default is z axis in
        parent's frame. Can accept 'X', 'Y', 'Z' or a 3 tuple describing the
        vector.
    child_axis : 3 Tuple or string, optional
        Defines the orientation as a vector which must be aligned with parent's
        axis before adding joint. If it is not passed, default is z axis in
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
        >>> pin_joint = PinJoint('pinjoint', parent, child,
        ...                      child_point_pos=(l, l, 0))
        >>> pin_joint.coordinates
        [pinjoint_theta(t)]
        >>> pin_joint.speeds
        [pinjoint_omega(t)]

    """
    def __init__(self, name, parent, child, coord, speed,
                 parent_point_pos=(0, 0, 0), child_point_pos=(0, 0, 0),
                 parent_axis=None, child_axis=None):


        super(PinJoint, self).__init__(name, parent, child, [coord], [speed],
                                       parent_point_pos, child_point_pos,
                                       parent_axis, child_axis)

        self.kin_diff.append(coord - speed)

        # Define spatial information
        self.joint_transform = Matrix([[cos(coord), sin(coord),  0, 0, 0, 0],
                                       [-sin(coord), cos(coord), 0, 0, 0, 0],
                                       [0,           0,          1, 0, 0, 0],
                                       [0, 0, 0,  cos(coord), sin(coord), 0],
                                       [0, 0, 0, -sin(coord), cos(coord), 0],
                                       [0, 0, 0,  0,          0,          1]])

        self.motion_subspace = Matrix([0, 0, 1, 0, 0, 0])

    def apply_joint(self):
        # Orient the two joint frames in the two bodies
        temp = self.parent_joint_frame.orientnew(self.name +
                                                 "_child_joint_frame", 'Axis',
                                                 [self.coordinates[0],
                                                  self.parent_joint_frame.z])
        self.child_joint_frame = temp
        self.child_joint_frame.set_ang_vel(self.parent_joint_frame,
                                           self.speeds[0] *
                                           self.parent_joint_frame.z)

        # Set the velocities of the joint's points for the two bodies
        self.child_joint_point.set_pos(self.parent_joint_point, 0)
        self.child_joint_point.set_vel(self.parent_joint_frame, 0)


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
    kin_diff : list
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
    coord : function of time
        This is the generalized coordinate to associate with the SlidingJoint.
    speed : function of time
        This is the generalized speed to associate with the SlidingJoint.
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
        axis before adding joint. If it is not passed, default is z axis in
        parent's frame. Can pass in 'X', 'Y', 'Z' or a 3 tuple describing the
        vector.
    child_axis : 3 Tuple or string, optional
        Defines the orientation as a vector which must be aligned with parent's
        axis before adding joint. If it is not passed, default is z axis in
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
        ...                              child_point_pos=(l, 0, 0))
        >>> sliding_joint.coordinates
        [slidingjoint_x(t)]
        >>> sliding_joint.speeds
        [slidingjoint_v(t)]

    """
    def __init__(self, name, parent, child, coord, speed, parent_point_pos=None,
                 child_point_pos=None, parent_axis=None, child_axis=None):


        super(SlidingJoint, self).__init__(name, parent, child, [coord],
                                           [speed], parent_point_pos,
                                           child_point_pos, parent_axis,
                                           child_axis)

        self.kin_diff.append(coord - speed)

        # Define spatial information
        self.joint_transform = Matrix([[1,     0,      0, 0, 0, 0],
                                       [0,     1,      0, 0, 0, 0],
                                       [0,     0,      1, 0, 0, 0],
                                       [0,     -coord, 0, 1, 0, 0],
                                       [coord, 0,      0, 0, 1, 0],
                                       [0,     0,      0, 0, 0, 1]])

        self.motion_subspace = Matrix([0, 0, 0, 0, 0, 1])

    def apply_joint(self):
        # Orient the two joint frames in the two bodies
        temp = self.parent_joint_frame.orientnew(self.name +
                                                 "_child_joint_frame", 'Axis',
                                                 [0, self.parent_joint_frame.z])
        self.child_joint_frame = temp
        self.child_joint_frame.set_ang_vel(self.parent_joint_frame,
                                           0*self.parent_joint_frame.z)

        # Set he velocities of the joint's points for the two bodies
        self.child_joint_point.set_pos(self.parent_joint_point,
                                       self.coordinates[0] *
                                       self.child_joint_frame.z)
        self.child_joint_point.set_vel(self.parent_joint_frame,
                                       self.speeds[0]*self.parent_joint_frame.z)
        self.child.masscenter.v2pt_theory(self.child_joint_point,
                                          self.parent_joint_frame,
                                          self.child_joint_frame)
