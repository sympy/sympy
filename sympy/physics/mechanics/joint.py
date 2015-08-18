from sympy import acos, Symbol
from sympy.physics.vector import cross, Vector, dot, dynamicsymbols
from sympy.physics.mechanics.functions import convert_tuple_to_vector

__all__ = ['Joint', 'PinJoint', 'SlidingJoint', 'CylindricalJoint',
           'SphericalJoint', 'PlanarJoint']


class Joint(object):
    """Abstract Base class for all specific joints.

    A Joint connects two bodies (a parent and child) by adding different degrees
    of freedom to child with respect to parent.
    This is the base class for all specific joints and holds all common methods
    acting as an interface for all joints. Custom joint can be created by
    creating a subclass of this class and overriding 'apply_joint()' method to
    add the dynamics of the joint.

    Parameters
    ----------
    name: String
        Name of the joint which makes it unique. Should be different from other
        joints.
    parent: Body
        Instance of Body class serving as the parent in the joint.
    child: Body
        Instance of Body class serving as the child in the joint.
    parent_point_pos: 3 Tuple, Optional
        Defines the joint's point where the parent will be connected to child.
        3 Tuple defines the values of x, y and z directions w.r.t parent's
        frame. If it is not supplied, center of mass is added as default.
    child_point_pos: 3 Tuple, Optional
        Defines the joint's point where the child will be connected to parent.
        3 Tuple defines the values of x, y and z directions w.r.t child's frame.
        If it is not supplied, center of mass is added as default.
    parent_axis: Vector, Optional
        Defines the orientation as a vector which must be aligned with child's
        axis before adding joint. If it is not passed, default is x axis in
        parent's frame.
    child_axis: Vector, Optional
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

    It is defined such that the child body rotates with respect to the parent body
    about the body fixed parent axis through the angle theta. The point of joint
    (pin joint) is defined by two points in each body which coincides together.
    They are supplied as parent_point_pos and child_point_pos.
    In addition, the child's reference frame can be arbitrarily rotated a constant
    amount with respect to the parent axis by passing an axis in child body which
    must align with parent axis after the rotation.

    Parameters
    ----------
    name: String
        Name of the joint which makes it unique. Should be different from other
        joints.
    parent: Body
        Instance of Body class serving as the parent in the joint.
    child: Body
        Instance of Body class serving as the child in the joint.
    parent_point_pos: 3 Tuple, Optional
        Defines the pin joint's point where the parent will be connected to child.
        3 Tuple defines the values of x, y and z directions w.r.t parent's
        frame. If it is not supplied, center of mass is added as default.
    child_point_pos: 3 Tuple, Optional
        Defines the joint's point where the child will be connected to parent.
        3 Tuple defines the values of x, y and z directions w.r.t child's frame.
        If it is not supplied, center of mass is added as default.
    parent_axis: Vector, Optional
        Defines the orientation as a vector which must be aligned with child's
        axis before adding joint. If it is not passed, default is x axis in
        parent's frame.
    child_axis: Vector, Optional
        Defines the orientation as a vector which must be aligned with parent's
        axis before adding joint. If it is not passed, default is x axis in
        child's frame.

    Examples
    --------
    Adding a Pin Joint which connects center of mass of parent to a point
    pointed by child.frame.x + child.frame.y w.r.t. child's center of mass.
    Gravity is along y axis of parent.

    >>> from sympy import Symbol
    >>> from sympy.physics.mechanics import Body, PinJoint
    >>> parent = Body('parent')
    >>> child = Body('child')
    >>> gravity = Symbol('gravity')
    >>> l = Symbol('l')
    >>> child.apply_force(child.mass * gravity * parent.frame.y, child.masscenter)
    >>> pin_joint = PinJoint('pinjoint', parent, \
                             child, child_point_pos=(l, l, 0))
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
            self.parent_axis = convert_tuple_to_vector(self.parent.frame, parent_axis)
        elif isinstance(parent_axis, str) and parent_axis.upper() \
                in parent_axes_str.keys():
            self.parent_axis = parent_axes_str[parent_axis.upper()]
        else:
            raise TypeError("Parent Axis must either be one of 'X', 'Y', 'Z'" +
                            " or a 3-Tuple.")

        if child_axis is None:
            self.child_axis = child.frame.x
        elif isinstance(child_axis, tuple):
            self.child_axis = convert_tuple_to_vector(self.child.frame, child_axis)
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

    Parameters
    ----------
    name: String
        Name of the joint which makes it unique. Should be different from other
        joints.
    parent: Body
        Instance of Body class serving as the parent in the joint.
    child: Body
        Instance of Body class serving as the child in the joint.
    parent_point_pos: 3 Tuple, Optional
        Defines the joint's point where the parent will be connected to child.
        3 Tuple defines the values of x, y and z directions w.r.t parent's
        frame. If it is not supplied, center of mass is added as default.
    child_point_pos: 3 Tuple, Optional
        Defines the joint's point where the child will be connected to parent.
        3 Tuple defines the values of x, y and z directions w.r.t child's frame.
        If it is not supplied, center of mass is added as default.
    parent_axis: Vector, Optional
        Defines the orientation as a vector which must be aligned with child's
        axis before adding joint. If it is not passed, default is x axis in
        parent's frame.
    child_axis: Vector, Optional
        Defines the orientation as a vector which must be aligned with parent's
        axis before adding joint. If it is not passed, default is x axis in
        child's frame.

    Examples
    --------
    Adds sliding Joint between parent's masscenter and a point located at unit
    distance in x axis of child.

    >>> from sympy import Symbol
    >>> from sympy.physics.mechanics import Body, SlidingJoint
    >>> parent = Body('parent')
    >>> child = Body('child')
    >>> l = Symbol('l')
    >>> sliding_joint = SlidingJoint('slidingjoint', parent, child, \
                                     child_point_pos=(l, 0, 0))
    >>> sliding_joint.coordinates
    [slidingjoint_dis(t)]
    >>> sliding_joint.speeds
    [slidingjoint_vel(t)]

    """
    def __init__(self, name, parent, child, parent_point_pos=None, child_point_pos=None,
                 parent_axis=None, child_axis=None):

        parent_axes_str = {'X': parent.frame.x, 'Y': parent.frame.y,
                           'Z': parent.frame.z}
        child_axes_str = {'X': child.frame.x, 'Y': child.frame.y,
                          'Z': child.frame.z}

        if parent_axis is None:
            self.parent_axis = parent.frame.x
        elif isinstance(parent_axis, tuple):
            self.parent_axis = convert_tuple_to_vector(self.parent.frame, parent_axis)
        elif isinstance(parent_axis, str) and parent_axis.upper() \
                in parent_axes_str.keys():
            self.parent_axis = parent_axes_str[parent_axis.upper()]
        else:
            raise TypeError("Parent Axis must either be one of 'X', 'Y', 'Z'" +
                            " or a 3-Tuple.")

        if child_axis is None:
            self.child_axis = child.frame.x
        elif isinstance(child_axis, tuple):
            self.child_axis = convert_tuple_to_vector(self.child.frame, child_axis)
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


class CylindricalJoint(Joint):
    """
    Cylindrical Joint.

    It is defined such that the child body rotates about as well as translates
    along the fixed parent axis with respect to parent body through the angle
    theta and displacement x. The point of joint (cylindrical joint) is defined
    by the two points in each body which coincides initially. They are supplied
    as parent_point_pos and child_point_pos.
    In addition, the child's reference frame can be arbitrarily rotated a
    constant amount with respect to the parent axis by passing an axis in child
    body which must align with parent axis after rotation.

    Parameters
    ----------
    name: String
        Name of the joint which makes it unique. Should be different from other
        joints.
    parent: Body
        Instance of Body class serving as the parent in the joint.
    child: Body
        Instance of Body class serving as the child in the joint.
    parent_point_pos: 3 Tuple, Optional
        Defines the joint's point where the parent will be connected to child.
        3 Tuple defines the values of x, y and z directions w.r.t parent's
        frame. If it is not supplied, center of mass is added as default.
    child_point_pos: 3 Tuple, Optional
        Defines the joint's point where the child will be connected to parent.
        3 Tuple defines the values of x, y and z directions w.r.t child's frame.
        If it is not supplied, center of mass is added as default.
    parent_axis: Vector, Optional
        Defines the orientation as a vector which must be aligned with child's
        axis before adding joint. If it is not passed, default is x axis in
        parent's frame.
    child_axis: Vector, Optional
        Defines the orientation as a vector which must be aligned with parent's
        axis before adding joint. If it is not passed, default is x axis in
        child's frame.

    Examples
    --------
    Adds cylindrical Joint between parent's masscenter and a point located by
    child.frame.x + child.frame.y in child.

    >>> from sympy.physics.mechanics import Body, CylindricalJoint
    >>> parent = Body('parent')
    >>> child = Body('child')
    >>> cylindrical_joint = CylindricalJoint('cylindricaljoint', parent, child, \
                                             child_point_pos=(1, 1, 0))
    >>> cylindrical_joint.coordinates
    [cylindricaljoint_dis(t), cylindricaljoint_theta(t)]
    >>> cylindrical_joint.speeds
    [cylindricaljoint_vel(t), cylindricaljoint_omega(t)]

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
            self.parent_axis = convert_tuple_to_vector(self.parent.frame, parent_axis)
        elif isinstance(parent_axis, str) and parent_axis.upper() \
                in parent_axes_str.keys():
            self.parent_axis = parent_axes_str[parent_axis.upper()]
        else:
            raise TypeError("Parent Axis must either be one of 'X', 'Y', 'Z'" +
                            " or a 3-Tuple.")

        if child_axis is None:
            self.child_axis = child.frame.x
        elif isinstance(child_axis, tuple):
            self.child_axis = convert_tuple_to_vector(self.child.frame, child_axis)
        elif isinstance(child_axis, str) and child_axis.upper() \
                in child_axes_str.keys():
            self.child_axis = child_axes_str[child_axis.upper()]
        else:
            raise TypeError("Child Axis must either be one of 'X', 'Y', 'Z'" +
                            " or a 3 Tuple")

        super(CylindricalJoint, self).__init__(name, parent, child,
                                               parent_point_pos, child_point_pos)

    def apply_joint(self):
        x = dynamicsymbols(self.name + '_x')
        xd = dynamicsymbols(self.name + '_x', 1)
        v = dynamicsymbols(self.name + '_v')

        theta = dynamicsymbols(self.name + '_theta')
        thetad = dynamicsymbols(self.name + '_theta', 1)
        omega = dynamicsymbols(self.name + '_omega')

        self.coordinates.append(x)
        self.speeds.append(v)
        self.kds.append(xd - v)

        self.coordinates.append(theta)
        self.speeds.append(omega)
        self.kds.append(thetad - omega)

        self.child.frame.orient(self.parent.frame, 'Axis',
                                [0, self.parent.frame.x])
        self._align_axes(self.parent_axis, self.child_axis)
        self._locate_joint_point()

        self.child.frame.orient(self.parent.frame, 'Axis',
                                [theta, self.parent_axis])
        self.child.frame.set_ang_vel(self.parent.frame,
                                     omega * self.parent_axis)

        self.parent_joint_point.set_vel(self.parent.frame, 0)
        self.child_joint_point.set_vel(self.child.frame, 0)

        self.child_joint_point.set_pos(self.parent_joint_point,
                                       x * self.parent_axis)
        self.child_joint_point.set_vel(self.parent.frame,
                                       v * self.parent_axis)
        self.child.masscenter.set_vel(self.parent.frame,
                                      v * self.parent_axis)
        self.child.masscenter.v2pt_theory(self.parent.masscenter,
                                          self.parent.frame,
                                          self.child.frame)


class PlanarJoint(Joint):
    """
    Planar Joint.

    It is defined such that child body can rotates about the body fixed parent
    axis through the angle theta and translates along the plane perpendicular
    to parent axis (parent axis is the normal to the plane) through displacement
    as x_x and x_y. along two perpendicular vectors in the plane.The plane
    passes though a point which is at a distance from the parent's center of
    mass and parent axis passes though it.
    To define translation along the plane, two perpendicular vector in the plane
    are calculated using parent axis and distance of plane from parent's center
    of mass, arbitrarily named vecx and vecy.
    In addition, the child's reference frame can be arbitrarily rotated a constant
    amount with respect to the parent axis by passing an axis in child body which
    must align with parent axis after the rotation.

    Parameters
    ----------
    name: String
        Name of the joint which makes it unique. Should be different from other
        joints.
    parent: Body
        Instance of Body class serving as the parent in the joint.
    child: Body
        Instance of Body class serving as the child in the joint.
    parent_point_pos: 3 Tuple, Optional
        Defines the joint's point where the parent will be connected to child.
        3 Tuple defines the values of x, y and z directions w.r.t parent's
        frame. If it is not supplied, center of mass is added as default.
    child_point_pos: 3 Tuple, Optional
        Defines the joint's point where the child will be connected to parent.
        3 Tuple defines the values of x, y and z directions w.r.t child's frame.
        If it is not supplied, center of mass is added as default.
    parent_axis: Vector, Optional
        Defines the orientation as a vector which must be aligned with child's
        axis before adding joint. If it is not passed, default is x axis in
        parent's frame.
    child_axis: Vector, Optional
        Defines the orientation as a vector which must be aligned with parent's
        axis before adding joint. If it is not passed, default is x axis in
        child's frame.

    Note: Both child_point_pos and parent_point_pos cannot be None, atleast
    one of them must be supplied.

    Examples
    --------
    Adds planar Joint between parent's masscenter and a point at unit distance
    along child's y axis.

    >>> from sympy import Symbol
    >>> from sympy.physics.mechanics import Body, PlanarJoint
    >>> parent = Body('parent')
    >>> child = Body('child')
    >>> l = Symbol('l')
    >>> planar_joint = PlanarJoint('planarjoint',parent,child,parent_distance=l)
    >>> planar_joint.coordinates
    [planarjoint_theta(t), planarjoint_disx(t), planarjoint_disy(t)]
    >>> planar_joint.speeds
    [planarjoint_omega(t), planarjoint_velx(t), planarjoint_vely(t)]

    """
    def __init__(self, name, parent, child, parent_axis=None, child_axis=None, parent_distance=None,
                 child_distance=None):

        parent_axes_str = {'X': parent.frame.x, 'Y': parent.frame.y,
                           'Z': parent.frame.z}
        child_axes_str = {'X': child.frame.x, 'Y': child.frame.y,
                          'Z': child.frame.z}

        if parent_axis is None:
            self.parent_axis = parent.frame.x
        elif isinstance(parent_axis, tuple):
            self.parent_axis = convert_tuple_to_vector(self.parent.frame, parent_axis)
        elif isinstance(parent_axis, str) and parent_axis.upper() \
                in parent_axes_str.keys():
            self.parent_axis = parent_axes_str[parent_axis.upper()]
        else:
            raise TypeError("Parent Axis must either be one of 'X', 'Y', 'Z'" +
                            " or a 3-Tuple.")

        if child_axis is None:
            self.child_axis = child.frame.x
        elif isinstance(child_axis, tuple):
            self.child_axis = convert_tuple_to_vector(self.child.frame, child_axis)
        elif isinstance(child_axis, str) and child_axis.upper() \
                in child_axes_str.keys():
            self.child_axis = child_axes_str[child_axis.upper()]
        else:
            raise TypeError("Child Axis must either be one of 'X', 'Y', 'Z'" +
                            " or a 3 Tuple")

        if parent_distance is None and child_distance is None:
            raise TypeError("Atleast one of parent_distance or chid_distance " +
                            "must be passed.")
        if parent_distance is not None:
            self.parent_distance = parent_distance * self.parent_axis
            _parent_distance_mat = self.parent_distance.to_matrix(parent.frame)
            parent_point_pos = (_parent_distance_mat[0], _parent_distance_mat[1],
                            _parent_distance_mat[2])
        else:
            parent_point_pos = None

        if child_distance is not None:
            self.child_distance = child_distance * self.child_axis
            _child_distance_mat = self.child_distance.to_matrix(child.frame)
            child_point_pos = (_child_distance_mat[0], _child_distance_mat[1],
                               _child_distance_mat[2])
        else:
            child_point_pos = None

        super(PlanarJoint, self).__init__(name, parent, child,
                                          parent_point_pos=parent_point_pos,
                                          child_point_pos=child_point_pos)

    def apply_joint(self):
        # generalized coordinates in specific order.
        theta = dynamicsymbols(self.name + '_theta')
        thetad = dynamicsymbols(self.name + '_theta', 1)
        omega = dynamicsymbols(self.name + '_omega')
        x_x = dynamicsymbols(self.name + '_x_x')
        x_xd = dynamicsymbols(self.name + '_x_x', 1)
        v_x = dynamicsymbols(self.name + '_v_x')
        x_y = dynamicsymbols(self.name + '_x_y')
        x_yd = dynamicsymbols(self.name + '_x_y', 1)
        v_y = dynamicsymbols(self.name + '_v_y')

        self.coordinates.append(theta)
        self.speeds.append(omega)
        self.kds.append(thetad - omega)

        self.coordinates.append(x_x)
        self.speeds.append(v_x)
        self.kds.append(x_xd - v_x)

        self.coordinates.append(x_y)
        self.speeds.append(v_y)
        self.kds.append(x_yd - v_y)

        # getting vectors in the joint plane, atleast one of them is not 0 vec.
        if self._parent_joint_location is not Vector(0):
            self.vecx = cross(self.parent_axis, self._parent_joint_location)
        else:
            self.vecx = cross(self.parent_axis, self._child_joint_location)

        self.vecy = cross(self.parent_axis, self.vecx)
        self.vecz = self.parent_axis  # To maintain symmetry, not used anywhere.

        # Adding rotation
        self.child.frame.orient(self.parent.frame, 'Axis',
                                [theta, self.parent_axis])
        self.child.frame.set_ang_vel(self.parent.frame,
                                     omega * self.parent_axis)
        self.child_joint_point.set_pos(self.parent_joint_point, 0)
        self._align_axes(self.parent_axis, self.child_axis)
        self._locate_joint_point()

        # Adding translation along self.vecx.
        self.child_joint_point.set_pos(self.parent_joint_point, x_x * self.vecx)
        self.child_joint_point.set_vel(self.parent.frame, v_x * self.vecx)

        # Adding translation along y axis
        self.child_joint_point.set_pos(self.parent_joint_point, x_y * self.vecy)
        self.child_joint_point.set_vel(self.parent.frame, v_y * self.vecy)

        self.child.masscenter.v2pt_theory(self.parent.masscenter,
                                          self.parent.frame, self.child.frame)


class SphericalJoint(Joint):
    """
    Spherical (Ball) Joint.

    Provides three rotational degrees of freedom to the child w.r.t parent.
    The three generalized coordinates (thetax, thetay and thetax) are always
    the three angles along three perpendicular vectors and the three
    generalized speeds (omegax, omegay and omegaz) are the angular velocities
    along those vectors respectively.

    Parameters
    ----------
    name: String
        Name of the joint which makes it unique. Should be different from other
        joints.
    parent: Body
        Instance of Body class serving as the parent in the joint.
    child: Body
        Instance of Body class serving as the child in the joint.
    parent_point_pos: 3 Tuple, Optional
        Defines the joint's point where the parent will be connected to child.
        3 Tuple defines the values of x, y and z directions w.r.t parent's
        frame. If it is not supplied, center of mass is added as default.
    child_point_pos: 3 Tuple, Optional
        Defines the joint's point where the child will be connected to parent.
        3 Tuple defines the values of x, y and z directions w.r.t child's frame.
        If it is not supplied, center of mass is added as default.
    parent_axis: Vector, Optional
        Defines the orientation as a vector which must be aligned with child's
        axis before adding joint. If it is not passed, default is x axis in
        parent's frame.
    child_axis: Vector, Optional
        Defines the orientation as a vector which must be aligned with parent's
        axis before adding joint. If it is not passed, default is x axis in
        child's frame.

    Examples
    --------
    Adds spherical Joint between parent's masscenter and a point at unit distance
    along child's y axis.

    >>> from sympy import Symbol
    >>> from sympy.physics.mechanics import Body, SphericalJoint
    >>> parent = Body('parent')
    >>> child = Body('child')
    >>> l = Symbol('l')
    >>> spherical_joint = SphericalJoint('sphericaljoint', parent, child, \
                                         child_point_pos=(0, l, 0))
    >>> spherical_joint.coordinates
    [sphericaljoint_thetax(t), sphericaljoint_thetay(t), sphericaljoint_thetaz(t)]
    >>> spherical_joint.speeds
    [sphericaljoint_omegax(t), sphericaljoint_omegay(t), sphericaljoint_omegaz(t)]

    """
    def __init__(self, name, parent, child, parent_point_pos=None,
                 child_point_pos=None):

        super(SphericalJoint, self).__init__(name, parent, child,
                                             parent_point_pos, child_point_pos)

    def apply_joint(self):
        theta_x = dynamicsymbols(self.name + '_theta_x')
        theta_y = dynamicsymbols(self.name + '_theta_y')
        theta_z = dynamicsymbols(self.name + '_theta_z')
        theta_xd = dynamicsymbols(self.name + '_theta_x', 1)
        theta_yd = dynamicsymbols(self.name + '_theta_y', 1)
        theta_zd = dynamicsymbols(self.name + '_theta_z', 1)
        omega_x = dynamicsymbols(self.name + '_omega_x')
        omega_y = dynamicsymbols(self.name + '_omega_y')
        omega_z = dynamicsymbols(self.name + '_omega_z')

        self.coordinates.append(theta_x)
        self.speeds.append(omega_x)
        self.kds.append(theta_xd - omega_x)

        self.coordinates.append(theta_y)
        self.speeds.append(omega_y)
        self.kds.append(theta_yd - omega_y)

        self.coordinates.append(theta_z)
        self.speeds.append(omega_z)
        self.kds.append(theta_zd - omega_z)

        self._locate_joint_point()
        self.child_joint_point.set_pos(self.parent_joint_point, 0)

        self.child.frame.orient(self.parent.frame, 'Body',
                                [theta_x, theta_y, theta_z], 'XYZ')
        self.child.frame.set_ang_vel(self.parent.frame,
                                     omega_x * self.parent.frame.x +
                                     omega_y * self.parent.frame.y +
                                     omega_z * self.parent.frame.z)
        self.child.masscenter.v2pt_theory(self.parent.masscenter,
                                          self.parent.frame,
                                          self.child.frame)
