from sympy.physics.vector.frame import ReferenceFrame
from sympy import acos
from sympy.physics.mechanics.body import Body
from sympy.physics.vector import Vector, dynamicsymbols, dot, cross
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

    def _align_axis(self, parent, child):
        # Returns the axis and angle between two axis(vectors).
        axis1 = parent.normalize()
        axis2 = child.normalize()
        angle = acos(dot(axis2, axis1))
        axis = cross(child, parent)
        return angle, axis


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

    This is minimal working example where parent body and child body is
    connected via PinJoint through their masscenters.

        >>> from sympy.physics.mechanics import Body, PinJoint
        >>> parent = Body('parent')
        >>> child = Body('child')
        >>> P = PinJoint('P', parent, child)
        >>> P.coordinates()
        [P_theta(t)]
        >>> P.speeds()
        [P_omega(t)]
        >>> P.kdes()
        [P_omega(t) - Derivative(P_theta(t), t)]

    This is an example of chaos pendulum where we will do all kinematics
    using PinJoints.

        >>> from sympy import symbols
        >>> from sympy.physics.mechanics import PinJoint, Body
        >>> from sympy.physics.vector import dynamicsymbols, ReferenceFrame

    Declaring the mass, frame, speeds, coordinates etc isn't necessary as `Body` and
    `PinJoint` can declare it themselves. But we would be declaring them for better understanding.

        >>> mA, mB, lA, lB, h = symbols('mA, mB, lA, lB, h')
        >>> theta, phi, omega, alpha = dynamicsymbols('theta phi omega alpha')
        >>> N = ReferenceFrame('N')
        >>> A = ReferenceFrame('A')
        >>> B = ReferenceFrame('B')

    Declaring the bodies.

        >>> rod = Body('rod', frame=A, mass=mA)
        >>> plate = Body('plate', mass=mB, frame=B)
        >>> C = Body('C', frame=N) #Ceiling

    Declaring the joint position wrt rod's masscenter, other bodies are connected through their
    masscenters.

        >>> lA = (lB - h / 2) / 2 * A.z
        >>> lC = (lB/2 + h/4) * A.z

    Rod's y axis is aligned to ceiling's(C) y axis.
    Rod's z axis is aligned to plate's z axis.

        >>> J1 = PinJoint('J1', C, rod, coordinates=theta, speeds=omega, child_joint_pos=lA, parent_axis=N.y, child_axis=A.y)
        >>> J2 = PinJoint('J2', rod, plate, coordinates=phi, speeds=alpha, parent_joint_pos=lC, parent_axis=A.z, child_axis=B.z)

    We can check the kinematics now.

        >>> J1.kdes()
        [omega(t) - Derivative(theta(t), t)]
        >>> J2.kdes()
        [alpha(t) - Derivative(phi(t), t)]

    This example can also be done in another way, i.e, without explicitly defining constants,
    dyamicsymbols, frames etc.

        >>> rod = Body('rod')
        >>> plate = Body('plate')
        >>> C = Body('C') #Ceiling
        >>> lA = (lB - h / 2) / 2 * rod.frame.z
        >>> lC = (lB/2 + h/4) * rod.frame.z
        >>> J1 = PinJoint('J1', C, rod, child_joint_pos=lA, parent_axis=C.frame.y, child_axis=rod.frame.y)
        >>> J2 = PinJoint('J2', rod, plate, parent_joint_pos=lC, parent_axis=rod.frame.z, child_axis=plate.frame.z)
        >>> J1.kdes()
        [J1_omega(t) - Derivative(J1_theta(t), t)]
        >>> J2.kdes()
        [J2_omega(t) - Derivative(J2_theta(t), t)]

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
        self._child.frame.set_ang_vel(self._parent.frame, self._speeds[0] * self._parent_axis)

    def _set_linear_velocity(self):
        self._parent_joint.set_vel(self._parent.frame, 0)
        self._child_joint.set_vel(self._parent.frame, 0)
        self._child_joint.set_pos(self._parent_joint, 0)
        self._child.masscenter.v2pt_theory(self._parent.masscenter, self._parent.frame, self._child.frame)
