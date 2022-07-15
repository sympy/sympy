# coding=utf-8

from abc import ABC, abstractmethod

from sympy.core.numbers import pi
from sympy.physics.mechanics.body import Body
from sympy.physics.vector import Vector, dynamicsymbols, cross
from sympy.physics.vector.frame import ReferenceFrame

import warnings

__all__ = ['Joint', 'PinJoint', 'PrismaticJoint']


class Joint(ABC):
    """Abstract base class for all specific joints.

    Explanation
    ===========

    A joint subtracts degrees of freedom from a body. This is the base class
    for all specific joints and holds all common methods acting as an interface
    for all joints. Custom joint can be created by inheriting Joint class and
    defining all abstract functions.

    The abstract methods are:

    - ``_generate_coordinates``
    - ``_generate_speeds``
    - ``_orient_frames``
    - ``_set_angular_velocity``
    - ``_set_linar_velocity``

    Parameters
    ==========

    name : string
        A unique name for the joint.
    parent : Body
        The parent body of joint.
    child : Body
        The child body of joint.
    coordinates: List of dynamicsymbols, optional
        Generalized coordinates of the joint.
    speeds : List of dynamicsymbols, optional
        Generalized speeds of joint.
    parent_joint_pos : Vector, optional
        Vector from the parent body's mass center to the point where the parent
        and child are connected. The default value is the zero vector.
    child_joint_pos : Vector, optional
        Vector from the child body's mass center to the point where the parent
        and child are connected. The default value is the zero vector.
    parent_axis : Vector, optional
        Axis fixed in the parent body which aligns with an axis fixed in the
        child body. The default is x axis in parent's reference frame.
    child_axis : Vector, optional
        Axis fixed in the child body which aligns with an axis fixed in the
        parent body. The default is x axis in child's reference frame.

    Attributes
    ==========

    name : string
        The joint's name.
    parent : Body
        The joint's parent body.
    child : Body
        The joint's child body.
    coordinates : list
        List of the joint's generalized coordinates.
    speeds : list
        List of the joint's generalized speeds.
    parent_point : Point
        The point fixed in the parent body that represents the joint.
    child_point : Point
        The point fixed in the child body that represents the joint.
    parent_axis : Vector
        The axis fixed in the parent frame that represents the joint.
    child_axis : Vector
        The axis fixed in the child frame that represents the joint.
    kdes : list
        Kinematical differential equations of the joint.

    Notes
    =====

    The direction cosine matrix between the child and parent is formed using a
    simple rotation about an axis that is normal to both ``child_axis`` and
    ``parent_axis``. In general, the normal axis is formed by crossing the
    ``child_axis`` into the ``parent_axis`` except if the child and parent axes
    are in exactly opposite directions. In that case the rotation vector is chosen
    using the rules in the following table where ``C`` is the child reference
    frame and ``P`` is the parent reference frame:

    .. list-table::
       :header-rows: 1

       * - ``child_axis``
         - ``parent_axis``
         - ``rotation_axis``
       * - ``-C.x``
         - ``P.x``
         - ``P.z``
       * - ``-C.y``
         - ``P.y``
         - ``P.x``
       * - ``-C.z``
         - ``P.z``
         - ``P.y``
       * - ``-C.x-C.y``
         - ``P.x+P.y``
         - ``P.z``
       * - ``-C.y-C.z``
         - ``P.y+P.z``
         - ``P.x``
       * - ``-C.x-C.z``
         - ``P.x+P.z``
         - ``P.y``
       * - ``-C.x-C.y-C.z``
         - ``P.x+P.y+P.z``
         - ``(P.x+P.y+P.z) × P.x``

    """

    def __init__(self, name, parent, child, coordinates=None, speeds=None,
                 parent_joint_pos=None, child_joint_pos=None, parent_axis=None,
                 child_axis=None):

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

        self._parent_axis = self._axis(parent, parent_axis)
        self._child_axis = self._axis(child, child_axis)

        self._parent_point = self._locate_joint_pos(parent, parent_joint_pos)
        self._child_point = self._locate_joint_pos(child, child_joint_pos)

        self._orient_frames()
        self._set_angular_velocity()
        self._set_linear_velocity()

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.__str__()

    @property
    def name(self):
        return self._name

    @property
    def parent(self):
        """Parent body of Joint."""
        return self._parent

    @property
    def child(self):
        """Child body of Joint."""
        return self._child

    @property
    def coordinates(self):
        """List generalized coordinates of the joint."""
        return self._coordinates

    @property
    def speeds(self):
        """List generalized coordinates of the joint.."""
        return self._speeds

    @property
    def kdes(self):
        """Kinematical differential equations of the joint."""
        return self._kdes

    @property
    def parent_axis(self):
        """The axis of parent frame."""
        return self._parent_axis

    @property
    def child_axis(self):
        """The axis of child frame."""
        return self._child_axis

    @property
    def parent_point(self):
        """The joint's point where parent body is connected to the joint."""
        return self._parent_point

    @property
    def child_point(self):
        """The joint's point where child body is connected to the joint."""
        return self._child_point

    @abstractmethod
    def _generate_coordinates(self, coordinates):
        """Generate list generalized coordinates of the joint."""
        pass

    @abstractmethod
    def _generate_speeds(self, speeds):
        """Generate list generalized speeds of the joint."""
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
        for i in range(len(self.coordinates)):
            kdes.append(-self.coordinates[i].diff(t) + self.speeds[i])
        return kdes

    def _axis(self, body, ax):
        if ax is None:
            ax = body.frame.x
            return ax
        if not isinstance(ax, Vector):
            raise TypeError("Axis must be of type Vector.")
        if not ax.dt(body.frame) == 0:
            msg = ('Axis cannot be time-varying when viewed from the '
                   'associated body.')
            raise ValueError(msg)
        return ax

    def _locate_joint_pos(self, body, joint_pos):
        if joint_pos is None:
            joint_pos = Vector(0)
        if not isinstance(joint_pos, Vector):
            raise ValueError('Joint Position must be supplied as Vector.')
        if not joint_pos.dt(body.frame) == 0:
            msg = ('Position Vector cannot be time-varying when viewed from '
                   'the associated body.')
            raise ValueError(msg)
        point_name = self._name + '_' + body.name + '_joint'
        return body.masscenter.locatenew(point_name, joint_pos)

    def _alignment_rotation(self, parent, child):
        # Returns the axis and angle between two axis(vectors).
        angle = parent.angle_between(child)
        axis = cross(child, parent).normalize()
        return angle, axis

    def _generate_vector(self):
        parent_frame = self.parent.frame
        components = self.parent_axis.to_matrix(parent_frame)
        x, y, z = components[0], components[1], components[2]

        if x != 0:
            if y!=0:
                if z!=0:
                    return cross(self.parent_axis,
                                parent_frame.x)
            if z!=0:
                return parent_frame.y
            return parent_frame.z

        if x == 0:
            if y!=0:
                if z!=0:
                    return parent_frame.x
                return parent_frame.x
            return parent_frame.y

    def _set_orientation(self):
        #Helper method for `orient_axis()`
        self.child.frame.orient_axis(self.parent.frame, self.parent_axis, 0)
        angle, axis = self._alignment_rotation(self.parent_axis,
                                                self.child_axis)

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)

            if axis != Vector(0) or angle == pi:

                if angle == pi:
                    axis = self._generate_vector()

                int_frame = ReferenceFrame('int_frame')
                int_frame.orient_axis(self.child.frame, self.child_axis, 0)
                int_frame.orient_axis(self.parent.frame, axis, angle)
                return int_frame
        return self.parent.frame


class PinJoint(Joint):
    """Pin (Revolute) Joint.

    .. image:: PinJoint.png

    Explanation
    ===========

    A pin joint is defined such that the joint rotation axis is fixed in both
    the child and parent and the location of the joint is relative to the mass
    center of each body. The child rotates an angle, θ, from the parent about
    the rotation axis and has a simple angular speed, ω, relative to the
    parent. The direction cosine matrix between the child and parent is formed
    using a simple rotation about an axis that is normal to both ``child_axis``
    and ``parent_axis``, see the Notes section for a detailed explanation of
    this.

    Parameters
    ==========

    name : string
        A unique name for the joint.
    parent : Body
        The parent body of joint.
    child : Body
        The child body of joint.
    coordinates: dynamicsymbol, optional
        Generalized coordinates of the joint.
    speeds : dynamicsymbol, optional
        Generalized speeds of joint.
    parent_joint_pos : Vector, optional
        Vector from the parent body's mass center to the point where the parent
        and child are connected. The default value is the zero vector.
    child_joint_pos : Vector, optional
        Vector from the child body's mass center to the point where the parent
        and child are connected. The default value is the zero vector.
    parent_axis : Vector, optional
        Axis fixed in the parent body which aligns with an axis fixed in the
        child body. The default is x axis in parent's reference frame.
    child_axis : Vector, optional
        Axis fixed in the child body which aligns with an axis fixed in the
        parent body. The default is x axis in child's reference frame.

    Attributes
    ==========

    name : string
        The joint's name.
    parent : Body
        The joint's parent body.
    child : Body
        The joint's child body.
    coordinates : list
        List of the joint's generalized coordinates.
    speeds : list
        List of the joint's generalized speeds.
    parent_point : Point
        The point fixed in the parent body that represents the joint.
    child_point : Point
        The point fixed in the child body that represents the joint.
    parent_axis : Vector
        The axis fixed in the parent frame that represents the joint.
    child_axis : Vector
        The axis fixed in the child frame that represents the joint.
    kdes : list
        Kinematical differential equations of the joint.

    Examples
    =========

    A single pin joint is created from two bodies and has the following basic
    attributes:

    >>> from sympy.physics.mechanics import Body, PinJoint
    >>> parent = Body('P')
    >>> parent
    P
    >>> child = Body('C')
    >>> child
    C
    >>> joint = PinJoint('PC', parent, child)
    >>> joint
    PinJoint: PC  parent: P  child: C
    >>> joint.name
    'PC'
    >>> joint.parent
    P
    >>> joint.child
    C
    >>> joint.parent_point
    PC_P_joint
    >>> joint.child_point
    PC_C_joint
    >>> joint.parent_axis
    P_frame.x
    >>> joint.child_axis
    C_frame.x
    >>> joint.coordinates
    [theta_PC(t)]
    >>> joint.speeds
    [omega_PC(t)]
    >>> joint.child.frame.ang_vel_in(joint.parent.frame)
    omega_PC(t)*P_frame.x
    >>> joint.child.frame.dcm(joint.parent.frame)
    Matrix([
    [1,                 0,                0],
    [0,  cos(theta_PC(t)), sin(theta_PC(t))],
    [0, -sin(theta_PC(t)), cos(theta_PC(t))]])
    >>> joint.child_point.pos_from(joint.parent_point)
    0

    To further demonstrate the use of the pin joint, the kinematics of simple
    double pendulum that rotates about the Z axis of each connected body can be
    created as follows.

    >>> from sympy import symbols, trigsimp
    >>> from sympy.physics.mechanics import Body, PinJoint
    >>> l1, l2 = symbols('l1 l2')

    First create bodies to represent the fixed ceiling and one to represent
    each pendulum bob.

    >>> ceiling = Body('C')
    >>> upper_bob = Body('U')
    >>> lower_bob = Body('L')

    The first joint will connect the upper bob to the ceiling by a distance of
    ``l1`` and the joint axis will be about the Z axis for each body.

    >>> ceiling_joint = PinJoint('P1', ceiling, upper_bob,
    ... child_joint_pos=-l1*upper_bob.frame.x,
    ... parent_axis=ceiling.frame.z,
    ... child_axis=upper_bob.frame.z)

    The second joint will connect the lower bob to the upper bob by a distance
    of ``l2`` and the joint axis will also be about the Z axis for each body.

    >>> pendulum_joint = PinJoint('P2', upper_bob, lower_bob,
    ... child_joint_pos=-l2*lower_bob.frame.x,
    ... parent_axis=upper_bob.frame.z,
    ... child_axis=lower_bob.frame.z)

    Once the joints are established the kinematics of the connected bodies can
    be accessed. First the direction cosine matrices of pendulum link relative
    to the ceiling are found:

    >>> upper_bob.frame.dcm(ceiling.frame)
    Matrix([
    [ cos(theta_P1(t)), sin(theta_P1(t)), 0],
    [-sin(theta_P1(t)), cos(theta_P1(t)), 0],
    [                0,                0, 1]])
    >>> trigsimp(lower_bob.frame.dcm(ceiling.frame))
    Matrix([
    [ cos(theta_P1(t) + theta_P2(t)), sin(theta_P1(t) + theta_P2(t)), 0],
    [-sin(theta_P1(t) + theta_P2(t)), cos(theta_P1(t) + theta_P2(t)), 0],
    [                              0,                              0, 1]])

    The position of the lower bob's masscenter is found with:

    >>> lower_bob.masscenter.pos_from(ceiling.masscenter)
    l1*U_frame.x + l2*L_frame.x

    The angular velocities of the two pendulum links can be computed with
    respect to the ceiling.

    >>> upper_bob.frame.ang_vel_in(ceiling.frame)
    omega_P1(t)*C_frame.z
    >>> lower_bob.frame.ang_vel_in(ceiling.frame)
    omega_P1(t)*C_frame.z + omega_P2(t)*U_frame.z

    And finally, the linear velocities of the two pendulum bobs can be computed
    with respect to the ceiling.

    >>> upper_bob.masscenter.vel(ceiling.frame)
    l1*omega_P1(t)*U_frame.y
    >>> lower_bob.masscenter.vel(ceiling.frame)
    l1*omega_P1(t)*U_frame.y + l2*(omega_P1(t) + omega_P2(t))*L_frame.y

    """

    def __init__(self, name, parent, child, coordinates=None, speeds=None,
                 parent_joint_pos=None, child_joint_pos=None, parent_axis=None,
                 child_axis=None):

        super().__init__(name, parent, child, coordinates, speeds,
                         parent_joint_pos, child_joint_pos, parent_axis,
                         child_axis)

    def __str__(self):
        return (f'PinJoint: {self.name}  parent: {self.parent}  '
                f'child: {self.child}')

    def _generate_coordinates(self, coordinate):
        coordinates = []
        if coordinate is None:
            theta = dynamicsymbols('theta' + '_' + self._name)
            coordinate = theta
        coordinates.append(coordinate)
        return coordinates

    def _generate_speeds(self, speed):
        speeds = []
        if speed is None:
            omega = dynamicsymbols('omega' + '_' + self._name)
            speed = omega
        speeds.append(speed)
        return speeds

    def _orient_frames(self):
        frame = self._set_orientation()
        self.child.frame.orient_axis(frame, self.parent_axis,
                                    self.coordinates[0])

    def _set_angular_velocity(self):
        self.child.frame.set_ang_vel(self.parent.frame, self.speeds[0] *
                                     self.parent_axis.normalize())

    def _set_linear_velocity(self):
        self.parent_point.set_vel(self.parent.frame, 0)
        self.child_point.set_vel(self.child.frame, 0)
        self.child_point.set_pos(self.parent_point, 0)
        self.child.masscenter.v2pt_theory(self.parent_point,
                                          self.parent.frame, self.child.frame)


class PrismaticJoint(Joint):
    """Prismatic (Sliding) Joint.

    .. image:: PrismaticJoint.png

    Explanation
    ===========

    It is defined such that the child body translates with respect to the parent
    body along the body fixed parent axis. The location of the joint is defined
    by two points in each body which coincides when the generalized coordinate is zero. The direction cosine matrix between
    the child and parent is formed using a simple rotation about an axis that is normal to
    both ``child_axis`` and ``parent_axis``, see the Notes section for a detailed explanation of
    this.

    Parameters
    ==========

    name : string
        A unique name for the joint.
    parent : Body
        The parent body of joint.
    child : Body
        The child body of joint.
    coordinates: dynamicsymbol, optional
        Generalized coordinates of the joint.
    speeds : dynamicsymbol, optional
        Generalized speeds of joint.
    parent_joint_pos : Vector, optional
        Vector from the parent body's mass center to the point where the parent
        and child are connected. The default value is the zero vector.
    child_joint_pos : Vector, optional
        Vector from the child body's mass center to the point where the parent
        and child are connected. The default value is the zero vector.
    parent_axis : Vector, optional
        Axis fixed in the parent body which aligns with an axis fixed in the
        child body. The default is x axis in parent's reference frame.
    child_axis : Vector, optional
        Axis fixed in the child body which aligns with an axis fixed in the
        parent body. The default is x axis in child's reference frame.

    Attributes
    ==========

    name : string
        The joint's name.
    parent : Body
        The joint's parent body.
    child : Body
        The joint's child body.
    coordinates : list
        List of the joint's generalized coordinates.
    speeds : list
        List of the joint's generalized speeds.
    parent_point : Point
        The point fixed in the parent body that represents the joint.
    child_point : Point
        The point fixed in the child body that represents the joint.
    parent_axis : Vector
        The axis fixed in the parent frame that represents the joint.
    child_axis : Vector
        The axis fixed in the child frame that represents the joint.
    kdes : list
        Kinematical differential equations of the joint.

    Examples
    =========

    A single prismatic joint is created from two bodies and has the following basic
    attributes:

    >>> from sympy.physics.mechanics import Body, PrismaticJoint
    >>> parent = Body('P')
    >>> parent
    P
    >>> child = Body('C')
    >>> child
    C
    >>> joint = PrismaticJoint('PC', parent, child)
    >>> joint
    PrismaticJoint: PC  parent: P  child: C
    >>> joint.name
    'PC'
    >>> joint.parent
    P
    >>> joint.child
    C
    >>> joint.parent_point
    PC_P_joint
    >>> joint.child_point
    PC_C_joint
    >>> joint.parent_axis
    P_frame.x
    >>> joint.child_axis
    C_frame.x
    >>> joint.coordinates
    [x_PC(t)]
    >>> joint.speeds
    [v_PC(t)]
    >>> joint.child.frame.ang_vel_in(joint.parent.frame)
    0
    >>> joint.child.frame.dcm(joint.parent.frame)
    Matrix([
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]])
    >>> joint.child_point.pos_from(joint.parent_point)
    x_PC(t)*P_frame.x

    To further demonstrate the use of the prismatic joint, the kinematics of
    two masses sliding, one moving relative to a fixed body and the other relative to the
    moving body. about the X axis of each connected body can be created as follows.

    >>> from sympy.physics.mechanics import PrismaticJoint, Body

    First create bodies to represent the fixed ceiling and one to represent
    a particle.

    >>> wall = Body('W')
    >>> Part1 = Body('P1')
    >>> Part2 = Body('P2')

    The first joint will connect the particle to the ceiling and the
    joint axis will be about the X axis for each body.

    >>> J1 = PrismaticJoint('J1', wall, Part1)

    The second joint will connect the second particle to the first particle
    and the joint axis will also be about the X axis for each body.

    >>> J2 = PrismaticJoint('J2', Part1, Part2)

    Once the joint is established the kinematics of the connected bodies can
    be accessed. First the direction cosine matrices of Part relative
    to the ceiling are found:

    >>> Part1.dcm(wall)
    Matrix([
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]])

    >>> Part2.dcm(wall)
    Matrix([
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]])

    The position of the particles' masscenter is found with:

    >>> Part1.masscenter.pos_from(wall.masscenter)
    x_J1(t)*W_frame.x

    >>> Part2.masscenter.pos_from(wall.masscenter)
    x_J1(t)*W_frame.x + x_J2(t)*P1_frame.x

    The angular velocities of the two particle links can be computed with
    respect to the ceiling.

    >>> Part1.ang_vel_in(wall)
    0

    >>> Part2.ang_vel_in(wall)
    0

    And finally, the linear velocities of the two particles can be computed
    with respect to the ceiling.

    >>> Part1.masscenter_vel(wall)
    v_J1(t)*W_frame.x

    >>> Part2.masscenter.vel(wall.frame)
    v_J1(t)*W_frame.x + Derivative(x_J2(t), t)*P1_frame.x

    """

    def __init__(self, name, parent, child, coordinates=None, speeds=None, parent_joint_pos=None,
        child_joint_pos=None, parent_axis=None, child_axis=None):

        super().__init__(name, parent, child, coordinates, speeds, parent_joint_pos,
                        child_joint_pos, parent_axis, child_axis)

    def __str__(self):
        return (f'PrismaticJoint: {self.name}  parent: {self.parent}  '
                f'child: {self.child}')

    def _generate_coordinates(self, coordinate):
        coordinates = []
        if coordinate is None:
            x = dynamicsymbols('x' + '_' + self._name)
            coordinate = x
        coordinates.append(coordinate)
        return coordinates

    def _generate_speeds(self, speed):
        speeds = []
        if speed is None:
            y = dynamicsymbols('v' + '_' + self._name)
            speed = y
        speeds.append(speed)
        return speeds

    def _orient_frames(self):
        frame = self._set_orientation()
        self.child.frame.orient_axis(frame, self.parent_axis, 0)

    def _set_angular_velocity(self):
        self.child.frame.set_ang_vel(self.parent.frame, 0)

    def _set_linear_velocity(self):
        self.parent_point.set_vel(self.parent.frame, 0)
        self.child_point.set_vel(self.child.frame, 0)
        self.child_point.set_pos(self.parent_point, self.coordinates[0] * self.parent_axis.normalize())
        self.child_point.set_vel(self.parent.frame, self.speeds[0] * self.parent_axis.normalize())
        self.child.masscenter.set_vel(self.parent.frame, self.speeds[0] * self.parent_axis.normalize())
