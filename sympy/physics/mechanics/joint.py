from sympy.physics.mechanics.body import Body
from sympy import acos, Symbol
from sympy.physics.vector import cross, Vector, dot, dynamicsymbols
from sympy.physics.mechanics.functions import convert_tuple_to_vector
from abc import ABC, abstractmethod

__all__ = ['Joint', 'PinJoint', 'SlidingJoint', 'CylindricalJoint',
           'SphericalJoint', 'PlanarJoint']


class Joint(ABC):
    """Abstract Base class for all specific joints.

    A Joint connects two bodies (a parent and child) by adding different degrees
    of freedom to child with respect to parent.
    This is the base class for all specific joints and holds all common methods
    acting as an interface for all joints. Custom joint can be created by
    inheriting Joint class and defining all abstract functions.

    Parameters
    ==========

    name : string
        The Joint's name which makes it unique.
    parent : Body
        The parent body of joint.
    child : Body
        The child body of joint.
    coordinates: List, optional
        Coordinates of joint.
    speeds : List, optional
        Speeds of joint.
    parent_joint_pos : Vector, optional
        Defines the joint's point where the parent will be connected to child.
        Default value is masscenter of Parent Body.
    child_joint_pos : Vector, optional
        Defines the joint's point where the child will be connected to parent.
        Default value is masscenter of Child Body. 
    parent_axis : Vector, optional
        Axis of parent frame which would be be aligned with child's
        axis. Default is x axis in parent's frame.
    child_axis : Vector, optional
        Axis of Child frame which would be be aligned with parent's
        axis. Default is x axis in child's frame.

    """

    def __init__(self, name, parent, child, coordinates=None, speeds=None, parent_joint_pos=None, 
        child_joint_pos=None, parent_axis = None, child_axis=None):

        if not isinstance(name, str):
            raise TypeError('Supply a valid name.')
        self._name = name

        if not isinstance(parent, Body):
            raise TypeError('Parent must be an instance of Body.')
        self._parent = parent

        if not isinstance(child, Body):
            raise TypeError('Parent must be an instance of Body.')
        self._child = child

        self._coordinates = self._generate_coordinates(coordinates)
        self._speeds = self._generate_speeds(speeds)
        self._kdes = self._generate_kdes()

        self.child_axis = self._axis(child, child_axis)
        self.parent_axis = self._axis(parent, parent_axis)

        self.parent_joint = self._locate_joint_pos(parent,parent_joint_pos)
        self.child_joint = self._locate_joint_pos(child, child_joint_pos)

        self._orient_frames()
        self._align_axes(self.parent_axis, self.child_axis)
        self._set_angular_velocity()
        self._set_linear_velocity()

    def __str__(self):
        return self._name

    def __repr__(self):
        return self.__str__()

    def parent(self):
        """ Parent body of Joint."""
        return self._parent

    def child(self):
        """ Child body of Joint."""
        return self._child

    def coordinates(self):
        """ List of Coordinates of Joint."""
        return self._coordinates

    def speeds(self):
        """ List of speeds of Joint."""
        return self._speeds

    def kdes(self):
        """ KDE of Joint."""
        return self._kdes

    @abstractmethod
    def _generate_coordinates(self, coordinates):
        """ Generate list of coordinates."""
        pass

    @abstractmethod
    def _generate_speeds(self, speeds):
        """ Generate list of speeds. """
        pass

    @abstractmethod
    def _orient_frames(self):
        """Orient frames as per the joint"""
        pass

    @abstractmethod
    def _set_angular_velocity(self):
        pass

    @abstractmethod
    def _set_linear_velocity(self):
        pass

    def _generate_kdes(self):
        kdes = []
        t = dynamicsymbols._t
        for i in range(len(self._coordinates)):
            kdes.append(-self._coordinates[0].diff(t) + self._speeds[0])

    def _axis(self, body, ax):
        if ax is None:
            ax = body.frame.x
            return ax
        if not isinstance(ax, Vector):
            raise TypeError("Axis must be of type Vector. Example-> A.x wehere 'A' is ReferenceFrame.")
        return ax

    def _locate_joint_pos(self, body, joint_pos):
        if joint_pos is None:
            joint_pos = Vector(0)
        if not isinstance(joint_pos, Vector):
            raise ValueError('Joint Position must be supplied as Vector.')
        
        return body.masscenter.locatenew(self.name + '_' + body.name + '_joint', joint_pos)

    def _align_axes(self, parent_axis, child_axis):
        """Rotates child_frame so that child_axis is aligned to parent_axis."""
        mag1 = parent_axis.magnitude()
        mag2 = child_axis.magnitude()
        angle = acos(dot(parent_axis, child_axis)/(mag1 * mag2))
        axis = cross(child_axis, parent_axis)
        if axis != Vector(0):
            self.child.frame.orient(
                self.parent.frame, 'Axis', [angle, axis])



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


class SphericalJoint(Joint):
    """
    Spherical (Ball) Joint.

    It is defined such that the child body rotates in any direction with respect
    to the parent body about the point of spherical joint through the angle
    theta_x, theta_y, theta_z along x, y and z directions of the parent's frame.
    The point of joint (pin joint) is defined by two points in each body which
    coincides together. They are supplied as parent_point_pos and child_point_pos.
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
    [sphericaljoint_alpha(t), sphericaljoint_beta(t), sphericaljoint_gamma(t)]
    >>> spherical_joint.speeds
    [sphericaljoint_omega_alpha(t), sphericaljoint_omega_beta(t), sphericaljoint_omega_gamma(t)]

    """
    def __init__(self, name, parent, child, parent_point_pos=None,
                 child_point_pos=None, rot_order=None):

        if rot_order is None:
            self.rot_order = 'XYZ'
        else:
            self.rot_order = rot_order

        super(SphericalJoint, self).__init__(name, parent, child,
                                             parent_point_pos, child_point_pos)

    def apply_joint(self):
        alpha = dynamicsymbols(self.name + '_alpha')
        beta = dynamicsymbols(self.name + '_beta')
        gamma = dynamicsymbols(self.name + '_gamma')
        alphad = dynamicsymbols(self.name + '_alpha', 1)
        betad = dynamicsymbols(self.name + '_beta', 1)
        gammad = dynamicsymbols(self.name + '_gamma', 1)
        omega_alpha = dynamicsymbols(self.name + '_omega_alpha')
        omega_beta = dynamicsymbols(self.name + '_omega_beta')
        omega_gamma = dynamicsymbols(self.name + '_omega_gamma')

        self.coordinates.append(alpha)
        self.speeds.append(omega_alpha)
        self.kds.append(alphad - omega_alpha)

        self.coordinates.append(beta)
        self.speeds.append(omega_beta)
        self.kds.append(betad - omega_beta)

        self.coordinates.append(gamma)
        self.speeds.append(omega_gamma)
        self.kds.append(gammad - omega_gamma)

        self._locate_joint_point()
        self.child_joint_point.set_pos(self.parent_joint_point, 0)
        self.child_joint_point.set_vel(self.child.frame, 0)
        self.child_joint_point.set_vel(self.parent.frame, 0)
        self.parent_joint_point.set_vel(self.parent.frame, 0)

        self.child.frame.orient(self.parent.frame, 'Body',
                                [alpha, beta, gamma], self.rot_order)
        self.child.frame.set_ang_vel(self.parent.frame,
                                     omega_alpha * self.parent.frame.x +
                                     omega_beta * self.parent.frame.y +
                                     omega_gamma * self.parent.frame.z)
        self.child.masscenter.v2pt_theory(self.parent.masscenter,
                                          self.parent.frame,
                                          self.child.frame)
