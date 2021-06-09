from sympy.physics.vector.frame import ReferenceFrame
from sympy import acos
from sympy.physics.mechanics.body import Body
from sympy.physics.vector import Vector, dynamicsymbols, dot
from abc import ABC, abstractmethod

__all__ = ['Joint', 'PinJoint']


class Joint(ABC):
    """Abstract Base class for all specific joints.

    Explanation
    ===========

    A joint subtracts degrees of freedom from a body.
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
    coordinates: List of dynamicsymbols/dynamicsymbol, optional
        Coordinates of joint.
    speeds : List of dynamicsymbols/dynamicsymbol, optional
        Speeds of joint.
    parent_joint_pos : Vector, optional
        Defines the joint's point where the parent will be connected to child.
        Default value is masscenter of Parent Body.
    child_joint_pos : Vector, optional
        Defines the joint's point where the child will be connected to parent.
        Default value is masscenter of Child Body.
    parent_axis : Vector, optional
        Axis of parent frame which would be aligned with child's
        axis. Default is x axis in parent's frame.
    child_axis : Vector, optional
        Axis of child frame which would be aligned with parent's
        axis. Default is x axis in child's frame.

    """

    def __init__(self, name, parent, child, coordinates=None, speeds=None, parent_joint_pos=None,
        child_joint_pos=None, parent_axis=None, child_axis=None):

        if not isinstance(name, str):
            raise TypeError('Supply a valid name.')
        self._name = name

        if not isinstance(parent, Body):
            raise TypeError('Parent must be an instance of Body.')
        self._parent = parent

        if not isinstance(child, Body):
            raise TypeError('Child must be an instance of Body.')
        self._child = child

        self._coordinates = self._generate_coordinates(coordinates)
        self._speeds = self._generate_speeds(speeds)
        self._kdes = self._generate_kdes()

        self._parent_axis = self._axis(parent, parent_axis)
        self._child_axis = self._axis(child, child_axis)

        self._parent_joint = self._locate_joint_pos(parent,parent_joint_pos)
        self._child_joint = self._locate_joint_pos(child, child_joint_pos)

        self._orient_frames()
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
        """ List of coordinates of Joint."""
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
            kdes.append(-self._coordinates[i].diff(t) + self._speeds[i])
        return kdes

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
        return body.masscenter.locatenew(self._name + '_' + body.name + '_joint', joint_pos)

    def _align_axis(self):
        axis1 = self._parent_axis
        axis2 = self._child_axis
        angle = acos(dot(axis1, axis2)/(axis1.magnitude() * axis2.magnitude()))
        if angle!=0:
            I = ReferenceFrame('I')
            I.orient_axis(self._parent.frame, axis1, angle)
            self._child.frame.orient_axis(I, axis1, self._coordinates[0])


class PinJoint(Joint):
    """Pin (Revolute) Joint.

    Explanation
    ===========

    It is defined such that the child body rotates with respect to the parent body
    about the body fixed parent axis through the angle theta. The point of joint
    (pin joint) is defined by two points in each body which coincides together.
    In addition, the child's reference frame can be arbitrarily rotated a constant
    amount with respect to the parent axis by passing an axis in child body which
    must align with parent axis after the rotation.

    Parameters
    ==========

    name : string
        The Joint's name which makes it unique.
    parent : Body
        The parent body of joint.
    child : Body
        The child body of joint.
    coordinates: dynamicsymbol, optional
        Coordinates of joint.
    speeds : dynamicsymbol, optional
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
        Axis of child frame which would be aligned with parent's
        axis. Default is x axis in child's frame.

    Examples
    =========

    >>> from sympy.physics.mechanics import Body, PinJoint
    >>> parent = Body('parent')
    >>> child = Body('child')
    >>> P = PinJoint('P', parent, child)
    >>> P.coordinates()
    [P_theta(t)]
    >>> P.speeds()
    [P_omega(t)]

    """

    def __init__(self, name, parent, child, coordinates=None, speeds=None, parent_joint_pos=None,
        child_joint_pos=None, parent_axis=None, child_axis=None):

        super().__init__(name, parent, child, coordinates, speeds, parent_joint_pos, child_joint_pos,
            parent_axis, child_axis)

    def _generate_coordinates(self, coordinate):
        coordinates = []
        if coordinate is None:
            theta = dynamicsymbols(self._name + '_theta')
            coordinate = theta
        coordinates.append(coordinate)
        return coordinates

    def _generate_speeds(self, speed):
        speeds = []
        if speed is None:
            omega = dynamicsymbols(self._name + '_omega')
            speed = omega
        speeds.append(speed)
        return speeds

    def _orient_frames(self):
        self._child.frame.orient_axis(self._parent.frame, self._parent_axis, self._coordinates[0])
        self._align_axis()

    def _set_angular_velocity(self):
        self._child.frame.set_ang_vel(self._parent.frame, self._speeds[0] * self._parent_axis)

    def _set_linear_velocity(self):
        self._parent_joint.set_vel(self._parent.frame, 0)
        self._child_joint.set_vel(self._parent.frame, 0)
        self._child_joint.set_pos(self._parent_joint, 0)
        self._child.masscenter.v2pt_theory(self._parent.masscenter, self._parent.frame, self._child.frame)
