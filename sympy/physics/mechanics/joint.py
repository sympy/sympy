from sympy.physics.vector.frame import ReferenceFrame
from sympy.physics.mechanics.body import Body
from sympy.physics.vector import Vector, dynamicsymbols, cross
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
    Abstract methods are:

    - ``'_generate_coordinates'``
    - ``'_generate_speeds'``
    - ``'_orient_frames'``
    - ``'_set_angular_velocity'``
    - ``'_set_linar_velocity'``

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
        """Parent body of Joint."""
        return self._parent

    def child(self):
        """Child body of Joint."""
        return self._child

    def coordinates(self):
        """ List of coordinates of Joint."""
        return self._coordinates

    def speeds(self):
        """ List of speeds of Joint."""
        return self._speeds

    def kdes(self):
        """Kinematical differential equations of the joint."""
        return self._kdes

    @abstractmethod
    def _generate_coordinates(self, coordinates):
        """List generalized coordinates of the joint."""
        pass

    @abstractmethod
    def _generate_speeds(self, speeds):
        """List generalized speeds of the joint."""
        pass

    @abstractmethod
    def _orient_frames(self):
        """Orient frames as per the joint."""
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
            raise TypeError("Axis must be of type Vector.")
        if not ax.dt(body.frame) == 0:
            raise ValueError('Axis cannot be time-varying when viewed from the associated body.')
        return ax

    def _locate_joint_pos(self, body, joint_pos):
        if joint_pos is None:
            joint_pos = Vector(0)
        if not isinstance(joint_pos, Vector):
            raise ValueError('Joint Position must be supplied as Vector.')
        if not joint_pos.dt(body.frame) == 0:
            raise ValueError('Position Vector cannot be time-varying when viewed from the associated body.')
        return body.masscenter.locatenew(self._name + '_' + body.name + '_joint', joint_pos)

    def _align_axis(self, parent, child):
        # Returns the axis and angle between two axis(vectors).
        angle = parent.angle_between(child)
        axis = cross(child, parent).normalize()
        return angle, axis


class PinJoint(Joint):
    """Pin (Revolute) Joint.

    Explanation
    ===========

    It is defined such that the pin joint axis is ``body fixed`` in both child and parent.
    The location of the joint is defined relative to the mass center of each body which
    coincides together.

    Parameters
    ==========

    name : string
        The joint's name which makes it unique.
    parent : Body
        The parent body of joint.
    child : Body
        The child body of joint.
    coordinate: dynamicsymbol, optional
        Coordinates of joint.
    speed : dynamicsymbol, optional
        Speeds of joint.
    parent_joint_pos : Vector, optional
        Locates the joint relative to the parent's mass center in the parent body.
        Default value is masscenter of Parent Body.
    child_joint_pos : Vector, optional
        locates the joint relative to the child's mass center in the child body.
        Default value is masscenter of Child Body.
    parent_axis : Vector, optional
        Axis of parent frame which would be be aligned with child's
        axis. Default is x axis in parent's frame.
    child_axis : Vector, optional
        Axis of child frame which would be aligned with parent's
        axis. Default is x axis in child's frame.

    Examples
    =========

    This is an example of simple double pendulum where we will do all the kinematics
    using PinJoints.

    >>> from sympy.physics.mechanics import PinJoint, Body
    >>> from sympy import symbols
    >>> l1, l2 = symbols('l1 l2')
    >>> ceiling = Body('C')
    >>> upper_bob = Body('U')
    >>> lower_bob = Body('L')
    >>> ceiling_joint = PinJoint('P1', ceiling, upper_bob, child_joint_pos=l1*upper_bob.frame.x,parent_axis=ceiling.frame.z, child_axis=upper_bob.frame.z)
    >>> pendulum_joint = PinJoint('P2', upper_bob, lower_bob, child_joint_pos=l2*lower_bob.frame.x, parent_axis=upper_bob.frame.z, child_axis=lower_bob.frame.z)
    >>> upper_bob.masscenter.vel(ceiling.frame)
    - l1*omega_P1(t)*U_frame.y
    >>> lower_bob.masscenter.vel(ceiling.frame)
    - l1*omega_P1(t)*U_frame.y - l2*(omega_P1(t) + omega_P2(t))*L_frame.y
    >>> upper_bob.frame.ang_vel_in(ceiling.frame)
    omega_P1(t)*C_frame.z
    >>> lower_bob.frame.ang_vel_in(ceiling.frame)
    omega_P1(t)*C_frame.z + omega_P2(t)*U_frame.z
    >>> upper_bob.frame.dcm(ceiling.frame)
    Matrix([
    [ cos(theta_P1(t)), sin(theta_P1(t)), 0],
    [-sin(theta_P1(t)), cos(theta_P1(t)), 0],
    [                0,                0, 1]])
    >>> lower_bob.frame.dcm(upper_bob.frame)
    Matrix([
    [ cos(theta_P2(t)), sin(theta_P2(t)), 0],
    [-sin(theta_P2(t)), cos(theta_P2(t)), 0],
    [                0,                0, 1]])

    Notes
    ======

    All dynamicsymbols used in `Joints` are `real ` by default. It is adviced to use
    `real=True` assumption when passing dynamicsymbols to `Joints` else it may give
    incorrect result in some cases.

    """

    def __init__(self, name, parent, child, coordinates=None, speeds=None, parent_joint_pos=None,
        child_joint_pos=None, parent_axis=None, child_axis=None):

        super().__init__(name, parent, child, coordinates, speeds, parent_joint_pos, child_joint_pos,
            parent_axis, child_axis)

    def _generate_coordinates(self, coordinate):
        coordinates = []
        if coordinate is None:
            theta = dynamicsymbols('theta' + '_' + self._name, real=True)
            coordinate = theta
        coordinates.append(coordinate)
        return coordinates

    def _generate_speeds(self, speed):
        speeds = []
        if speed is None:
            omega = dynamicsymbols('omega' + '_' + self._name, real=True)
            speed = omega
        speeds.append(speed)
        return speeds

    def _orient_frames(self):
        self._child.frame.orient_axis(self._parent.frame, self._parent_axis, 0)
        angle, axis = self._align_axis(self._parent_axis, self._child_axis)
        if axis != Vector(0):
            I = ReferenceFrame('I')
            I.orient_axis(self._child.frame, self._child_axis, 0)
            I.orient_axis(self._parent.frame, axis, angle)
            self._child.frame.orient_axis(I, self._parent_axis, self._coordinates[0])
        else:
            self._child.frame.orient_axis(self._parent.frame, self._parent_axis, self._coordinates[0])

    def _set_angular_velocity(self):
        self._child.frame.set_ang_vel(self._parent.frame, self._speeds[0] * self._parent_axis.normalize())

    def _set_linear_velocity(self):
        self._parent_joint.set_vel(self._parent.frame, 0)
        self._child_joint.set_vel(self._parent.frame, 0)
        self._child_joint.set_pos(self._parent_joint, 0)
        self._child.masscenter.v2pt_theory(self._parent.masscenter, self._parent.frame, self._child.frame)
