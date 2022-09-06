# coding=utf-8

from abc import ABC, abstractmethod

from sympy.core.numbers import pi
from sympy.physics.mechanics.body import Body
from sympy.physics.vector import (Vector, dynamicsymbols, cross, Point,
                                  ReferenceFrame)

from sympy.utilities.exceptions import sympy_deprecation_warning
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
    - ``_set_linear_velocity``

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
    parent_point : Point or Vector, optional
        Attachment point where the joint is fixed to the parent body. If a
        vector is provided, then the attachment point is computed by adding the
        vector to the body's mass center. The default value is the parent's mass
        center.
    child_point : Point or Vector, optional
        Attachment point where the joint is fixed to the child body. If a
        vector is provided, then the attachment point is computed by adding the
        vector to the body's mass center. The default value is the child's mass
        center.
    parent_axis : Vector, optional
        .. deprecated:: 1.12
            Axis fixed in the parent body which aligns with an axis fixed in the
            child body. The default is the x axis of parent's reference frame.
            For more information on this deprecation, see
            :ref:`deprecated-mechanics-joint-axis`.
    child_axis : Vector, optional
        .. deprecated:: 1.12
            Axis fixed in the child body which aligns with an axis fixed in the
            parent body. The default is the x axis of child's reference frame.
            For more information on this deprecation, see
            :ref:`deprecated-mechanics-joint-axis`.
    parent_interframe : ReferenceFrame, optional
        Intermediate frame of the parent body with respect to which the joint
        transformation is formulated. If a Vector is provided then an interframe
        is created which aligns its X axis with the given vector. The default
        value is the parent's own frame.
    child_interframe : ReferenceFrame, optional
        Intermediate frame of the child body with respect to which the joint
        transformation is formulated. If a Vector is provided then an interframe
        is created which aligns its X axis with the given vector. The default
        value is the child's own frame.
    parent_joint_pos : Point or Vector, optional
        .. deprecated:: 1.12
            This argument is replaced by parent_point and will be removed in a
            future version.
            See :ref:`deprecated-mechanics-joint-pos` for more information.
    child_joint_pos : Point or Vector, optional
        .. deprecated:: 1.12
            This argument is replaced by child_point and will be removed in a
            future version.
            See :ref:`deprecated-mechanics-joint-pos` for more information.

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
        Attachment point where the joint is fixed to the parent body.
    child_point : Point
        Attachment point where the joint is fixed to the child body.
    parent_axis : Vector
        The axis fixed in the parent frame that represents the joint.
    child_axis : Vector
        The axis fixed in the child frame that represents the joint.
    parent_interframe : ReferenceFrame
        Intermediate frame of the parent body with respect to which the joint
        transformation is formulated.
    child_interframe : ReferenceFrame
        Intermediate frame of the child body with respect to which the joint
        transformation is formulated.
    kdes : list
        Kinematical differential equations of the joint.

    Notes
    =====

    When providing a vector as the intermediate frame, a new intermediate frame
    is created which aligns its X axis with the provided vector. This is done
    with a single fixed rotation around a rotation axis. This rotation axis is
    determined by taking the cross product of the ``body.x`` axis with the
    provided vector. In the case where the provided vector is in the ``-body.x``
    direction, the rotation is done around the ``body.y`` axis.

    """

    def __init__(self, name, parent, child, coordinates=None, speeds=None,
                 parent_point=None, child_point=None, parent_axis=None,
                 child_axis=None, parent_interframe=None, child_interframe=None,
                 parent_joint_pos=None, child_joint_pos=None):

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

        self._parent_axis = self._axis(parent_axis, parent.frame)
        self._child_axis = self._axis(child_axis, child.frame)

        if parent_joint_pos is not None or child_joint_pos is not None:
            sympy_deprecation_warning(
                """
                The parent_joint_pos and child_joint_pos arguments for the Joint
                classes are deprecated. Instead use parent_point and child_point.
                """,
                deprecated_since_version="1.12",
                active_deprecations_target="deprecated-mechanics-joint-pos",
                stacklevel=4
            )
            if parent_point is None:
                parent_point = parent_joint_pos
            if child_point is None:
                child_point = child_joint_pos
        self._parent_point = self._locate_joint_pos(parent, parent_point)
        self._child_point = self._locate_joint_pos(child, child_point)
        if parent_axis is not None or child_axis is not None:
            sympy_deprecation_warning(
                """
                The parent_axis and child_axis arguments for the Joint classes
                are deprecated. Instead use parent_interframe, child_interframe.
                """,
                deprecated_since_version="1.12",
                active_deprecations_target="deprecated-mechanics-joint-axis",
                stacklevel=4
            )
            if parent_interframe is None:
                parent_interframe = parent_axis
            if child_interframe is None:
                child_interframe = child_axis
        self._parent_interframe = self._locate_joint_frame(parent,
                                                           parent_interframe)
        self._child_interframe = self._locate_joint_frame(child,
                                                          child_interframe)

        self._orient_frames()
        self._set_angular_velocity()
        self._set_linear_velocity()

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.__str__()

    @property
    def name(self):
        """Name of the joint."""
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
        # Will be removed with `deprecated-mechanics-joint-axis`
        return self._parent_axis

    @property
    def child_axis(self):
        """The axis of child frame."""
        # Will be removed with `deprecated-mechanics-joint-axis`
        return self._child_axis

    @property
    def parent_point(self):
        """Attachment point where the joint is fixed to the parent body."""
        return self._parent_point

    @property
    def child_point(self):
        """Attachment point where the joint is fixed to the child body."""
        return self._child_point

    @property
    def parent_interframe(self):
        return self._parent_interframe

    @property
    def child_interframe(self):
        return self._child_interframe

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
        """Set angular velocity of the joint related frames."""
        pass

    @abstractmethod
    def _set_linear_velocity(self):
        """Set velocity of related points to the joint."""
        pass

    @staticmethod
    def _to_vector(matrix, frame):
        """Converts a matrix to a vector in the given frame."""
        return Vector([(matrix, frame)])

    @staticmethod
    def _axis(ax, *frames):
        """Check whether an axis is fixed in one of the frames."""
        if ax is None:
            ax = frames[0].x
            return ax
        if not isinstance(ax, Vector):
            raise TypeError("Axis must be a Vector.")
        ref_frame = None  # Find a body in which the axis can be expressed
        for frame in frames:
            try:
                ax.to_matrix(frame)
                ref_frame = frame
                break
            except ValueError:
                pass
        if ref_frame is None:
            raise ValueError("Axis cannot be expressed in one of the body's "
                             "frames.")
        if not ax.dt(ref_frame) == 0:
            raise ValueError('Axis cannot be time-varying when viewed from the '
                             'associated body.')
        return ax

    @staticmethod
    def _choose_rotation_axis(frame, axis):
        components = axis.to_matrix(frame)
        x, y, z = components[0], components[1], components[2]

        if x != 0:
            if y != 0:
                if z != 0:
                    return cross(axis, frame.x)
            if z != 0:
                return frame.y
            return frame.z

        if x == 0:
            if y != 0:
                if z != 0:
                    return frame.x
                return frame.x
            return frame.y

    @staticmethod
    def _create_aligned_interframe(frame, align_axis, frame_axis=None,
                                   frame_name=None):
        """
        Returns an intermediate frame, where the ``frame_axis`` defined in
        ``frame`` is aligned with ``axis``. By default this means that the X
        axis will be aligned with ``axis``.

        Parameters
        ==========

        frame: Body or ReferenceFrame
            The body or reference frame with respect to which the intermediate
            frame is oriented
        align_axis: Vector
            The vector with respect to which the intermediate frame will be
            aligned.
        frame_axis: Vector
            The vector of the frame which should get aligned with ``axis``. The
            default is the X axis of the frame.
        frame_name: string
            Name of the to be created intermediate frame. The default adds
            "_int_frame" to the name of ``frame``.

        Example
        =======

        An intermediate frame, where the X axis of the parent becomes aligned
        with ``parent.y + parent.z`` can be created as follows:

        >>> from sympy.physics.mechanics.joint import Joint
        >>> from sympy.physics.mechanics import Body
        >>> parent = Body('parent')
        >>> parent_interframe = Joint._create_aligned_interframe(
        ...     parent, parent.y + parent.z)
        >>> parent_interframe
        parent_int_frame
        >>> parent.dcm(parent_interframe)
        Matrix([
        [        0, -sqrt(2)/2, -sqrt(2)/2],
        [sqrt(2)/2,        1/2,       -1/2],
        [sqrt(2)/2,       -1/2,        1/2]])
        >>> (parent.y + parent.z).express(parent_interframe)
        sqrt(2)*parent_int_frame.x

        Notes
        =====

        The direction cosine matrix between the given frame and intermediate
        frame is formed using a simple rotation about an axis that is normal to
        both ``align_axis`` and ``frame_axis``. In general, the normal axis is
        formed by crossing the ``frame_axis`` with the ``align_axis`` except if
        the axes are in exactly the opposite direction. In that case the
        rotation vector is chosen using the rules in the following table with
        the vectors expressed in the given frame:

        .. list-table::
           :header-rows: 1

           * - ``align_axis``
             - ``frame_axis``
             - ``rotation_axis``
           * - ``-x``
             - ``x``
             - ``z``
           * - ``-y``
             - ``y``
             - ``x``
           * - ``-z``
             - ``z``
             - ``y``
           * - ``-x-y``
             - ``x+y``
             - ``z``
           * - ``-y-z``
             - ``y+z``
             - ``x``
           * - ``-x-z``
             - ``x+z``
             - ``y``
           * - ``-x-y-z``
             - ``x+y+z``
             - ``(x+y+z) × x``

        """
        if isinstance(frame, Body):
            frame = frame.frame
        if frame_axis is None:
            frame_axis = frame.x
        if frame_name is None:
            if frame.name[-6:] == '_frame':
                frame_name = f'{frame.name[:-6]}_int_frame'
            else:
                frame_name = f'{frame.name}_int_frame'
        angle = frame_axis.angle_between(align_axis)
        rotation_axis = cross(frame_axis, align_axis)
        if rotation_axis != Vector(0) or angle == pi:
            if angle == pi:
                rotation_axis = Joint._choose_rotation_axis(frame, align_axis)

            int_frame = ReferenceFrame(frame_name)
            int_frame.orient_axis(frame, rotation_axis, angle)
            int_frame.set_ang_vel(frame, 0 * rotation_axis)
            return int_frame
        return frame

    def _generate_kdes(self):
        """Generate kinematical differential equations."""
        kdes = []
        t = dynamicsymbols._t
        for i in range(len(self.coordinates)):
            kdes.append(-self.coordinates[i].diff(t) + self.speeds[i])
        return kdes

    def _locate_joint_pos(self, body, joint_pos):
        """Returns the attachment point of a body."""
        if joint_pos is None:
            return body.masscenter
        if not isinstance(joint_pos, (Point, Vector)):
            raise TypeError('Attachment point must be a Point or Vector.')
        if isinstance(joint_pos, Vector):
            point_name = self._name + '_' + body.name + '_joint'
            joint_pos = body.masscenter.locatenew(point_name, joint_pos)
        if not joint_pos.pos_from(body.masscenter).dt(body.frame) == 0:
            raise ValueError('Attachment point must be fixed to the associated '
                             'body.')
        return joint_pos

    def _locate_joint_frame(self, body, interframe):
        """Returns the attachment frame of a body."""
        if interframe is None:
            return body.frame
        if isinstance(interframe, Vector):
            interframe = Joint._create_aligned_interframe(
                body, interframe,
                frame_name=f'{self.name}_{body.name}_int_frame')
        elif not isinstance(interframe, ReferenceFrame):
            raise TypeError('Interframe must be a ReferenceFrame.')
        if not interframe.ang_vel_in(body.frame) == 0:
            raise ValueError(f'Interframe {interframe} is not fixed to body '
                             f'{body}.')
        body.masscenter.set_vel(interframe, 0)  # Fixate interframe to body
        return interframe

    def _alignment_rotation(self, parent, child):
        # Returns the axis and angle between two axis(vectors), as if the parent
        # and child body are not rotated with respect to each other.
        # Will be removed with `deprecated-mechanics-joint-axis`
        child = self._to_vector(child.to_matrix(self.child.frame),
                                self.parent.frame)
        angle = parent.angle_between(child)
        axis = cross(child, parent).normalize()
        return angle, axis

    def _generate_vector(self):
        # Will be removed with `deprecated-mechanics-joint-axis`
        return self._choose_rotation_axis(self.parent.frame, self.parent_axis)

    def _set_orientation(self):
        # Helper method to determine parent_interframe based on parent_axis and
        # child_axis.
        # Will be removed with `deprecated-mechanics-joint-axis`
        angle, axis = self._alignment_rotation(self.parent_axis,
                                               self.child_axis)

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)

            if axis != Vector(0) or angle == pi:

                if angle == pi:
                    axis = self._generate_vector()

                int_frame = ReferenceFrame('int_frame')
                int_frame.orient_axis(self.parent.frame, axis, angle)
                return int_frame
        return self.parent.frame


class _JointAxisMixin:
    """Mixin class, which provides joint axis support methods.

    This class is strictly only to be inherited by joint types that have a joint
    axis. This is because the methods of this class assume that the object has
    a joint_axis attribute.
    """

    __slots__ = ()

    def _express_joint_axis(self, frame):
        """Helper method to express the joint axis in a specified frame."""
        try:
            ax_mat = self.joint_axis.to_matrix(self.parent_interframe)
        except ValueError:
            ax_mat = self.joint_axis.to_matrix(self.child_interframe)
        try:
            self.parent_interframe.dcm(frame)  # Check if connected
        except ValueError:
            self.child_interframe.dcm(frame)  # Should be connected
            int_frame = self.child_interframe
        else:
            int_frame = self.parent_interframe
        return self._to_vector(ax_mat, int_frame).express(frame)


class PinJoint(Joint, _JointAxisMixin):
    """Pin (Revolute) Joint.

    .. image:: PinJoint.svg

    Explanation
    ===========

    A pin joint is defined such that the joint rotation axis is fixed in both
    the child and parent and the location of the joint is relative to the mass
    center of each body. The child rotates an angle, θ, from the parent about
    the rotation axis and has a simple angular speed, ω, relative to the
    parent. The direction cosine matrix between the childinterframe and
    parentinterframe is formed using a simple rotation about the joint axis.
    The page on the joints framework gives a more detailed explanation on the
    intermediate frames.

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
    parent_point : Point or Vector, optional
        Attachment point where the joint is fixed to the parent body. If a
        vector is provided, then the attachment point is computed by adding the
        vector to the body's mass center. The default value is the parent's mass
        center.
    child_point : Point or Vector, optional
        Attachment point where the joint is fixed to the child body. If a
        vector is provided, then the attachment point is computed by adding the
        vector to the body's mass center. The default value is the child's mass
        center.
    parent_axis : Vector, optional
        .. deprecated:: 1.12
            Axis fixed in the parent body which aligns with an axis fixed in the
            child body. The default is the x axis of parent's reference frame.
            For more information on this deprecation, see
            :ref:`deprecated-mechanics-joint-axis`.
    child_axis : Vector, optional
        .. deprecated:: 1.12
            Axis fixed in the child body which aligns with an axis fixed in the
            parent body. The default is the x axis of child's reference frame.
            For more information on this deprecation, see
            :ref:`deprecated-mechanics-joint-axis`.
    parent_interframe : ReferenceFrame, optional
        Intermediate frame of the parent body with respect to which the joint
        transformation is formulated. If a Vector is provided then an interframe
        is created which aligns its X axis with the given vector. The default
        value is the parent's own frame.
    child_interframe : ReferenceFrame, optional
        Intermediate frame of the child body with respect to which the joint
        transformation is formulated. If a Vector is provided then an interframe
        is created which aligns its X axis with the given vector. The default
        value is the child's own frame.
    joint_axis : Vector
        The axis around which the rotation occurs. Note that the components
        of this axis are the same in the parent_interframe and child_interframe.
    parent_joint_pos : Point or Vector, optional
        .. deprecated:: 1.12
            This argument is replaced by parent_point and will be removed in a
            future version.
            See :ref:`deprecated-mechanics-joint-pos` for more information.
    child_joint_pos : Point or Vector, optional
        .. deprecated:: 1.12
            This argument is replaced by child_point and will be removed in a
            future version.
            See :ref:`deprecated-mechanics-joint-pos` for more information.

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
        Attachment point where the joint is fixed to the parent body.
    child_point : Point
        Attachment point where the joint is fixed to the child body.
    parent_axis : Vector
        The axis fixed in the parent frame that represents the joint.
    child_axis : Vector
        The axis fixed in the child frame that represents the joint.
    parent_interframe : ReferenceFrame
        Intermediate frame of the parent body with respect to which the joint
        transformation is formulated.
    child_interframe : ReferenceFrame
        Intermediate frame of the child body with respect to which the joint
        transformation is formulated.
    joint_axis : Vector
        The axis around which the rotation occurs. Note that the components of
        this axis are the same in the parent_interframe and child_interframe.
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
    P_masscenter
    >>> joint.child_point
    C_masscenter
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
    ... child_point=-l1*upper_bob.frame.x,
    ... joint_axis=ceiling.frame.z)

    The second joint will connect the lower bob to the upper bob by a distance
    of ``l2`` and the joint axis will also be about the Z axis for each body.

    >>> pendulum_joint = PinJoint('P2', upper_bob, lower_bob,
    ... child_point=-l2*lower_bob.frame.x,
    ... joint_axis=upper_bob.frame.z)

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
                 parent_point=None, child_point=None, parent_axis=None,
                 child_axis=None, parent_interframe=None, child_interframe=None,
                 joint_axis=None, parent_joint_pos=None, child_joint_pos=None):

        self._joint_axis = joint_axis
        super().__init__(name, parent, child, coordinates, speeds, parent_point,
                         child_point, parent_axis, child_axis,
                         parent_interframe, child_interframe, parent_joint_pos,
                         child_joint_pos)

    def __str__(self):
        return (f'PinJoint: {self.name}  parent: {self.parent}  '
                f'child: {self.child}')

    @property
    def joint_axis(self):
        """Axis along which the child rotates with respect to the parent."""
        return self._joint_axis

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
        self._joint_axis = self._axis(
            self.joint_axis, self.parent_interframe, self.child_interframe)
        axis = self._express_joint_axis(self.parent_interframe)
        self.child_interframe.orient_axis(
            self.parent_interframe, axis, self.coordinates[0])

    def _set_angular_velocity(self):
        self.child_interframe.set_ang_vel(self.parent_interframe, self.speeds[
            0] * self.joint_axis.normalize())

    def _set_linear_velocity(self):
        self.child_point.set_pos(self.parent_point, 0)
        self.parent_point.set_vel(self.parent.frame, 0)
        self.child_point.set_vel(self.child.frame, 0)
        self.child.masscenter.v2pt_theory(self.parent_point,
                                          self.parent.frame, self.child.frame)


class PrismaticJoint(Joint, _JointAxisMixin):
    """Prismatic (Sliding) Joint.

    .. image:: PrismaticJoint.svg

    Explanation
    ===========

    It is defined such that the child body translates with respect to the parent
    body along the body fixed joint axis. The location of the joint is defined
    by two points, one in each body, which coincides when the generalized
    coordinate is zero. The direction cosine matrix between the
    parent_interframe and child_interframe is the identity matrix. So the
    direction cosine matrix between the parent and child frames, is fully
    defined by the definition of the intermediate frames. For a detailed
    explanation on the intermediate frames, see the page on the joints
    framework.

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
    parent_point : Point or Vector, optional
        Attachment point where the joint is fixed to the parent body. If a
        vector is provided, then the attachment point is computed by adding the
        vector to the body's mass center. The default value is the parent's mass
        center.
    child_point : Point or Vector, optional
        Attachment point where the joint is fixed to the child body. If a
        vector is provided, then the attachment point is computed by adding the
        vector to the body's mass center. The default value is the child's mass
        center.
    parent_axis : Vector, optional
        .. deprecated:: 1.12
            Axis fixed in the parent body which aligns with an axis fixed in the
            child body. The default is the x axis of parent's reference frame.
            For more information on this deprecation, see
            :ref:`deprecated-mechanics-joint-axis`.
    child_axis : Vector, optional
        .. deprecated:: 1.12
            Axis fixed in the child body which aligns with an axis fixed in the
            parent body. The default is the x axis of child's reference frame.
            For more information on this deprecation, see
            :ref:`deprecated-mechanics-joint-axis`.
    parent_interframe : ReferenceFrame, optional
        Intermediate frame of the parent body with respect to which the joint
        transformation is formulated. If a Vector is provided then an interframe
        is created which aligns its X axis with the given vector. The default
        value is the parent's own frame.
    child_interframe : ReferenceFrame, optional
        Intermediate frame of the child body with respect to which the joint
        transformation is formulated. If a Vector is provided then an interframe
        is created which aligns its X axis with the given vector. The default
        value is the child's own frame.
    joint_axis : Vector
        The axis across which the translation occurs. Note that the components
        of this axis are the same in the parent_interframe and child_interframe.
    parent_joint_pos : Point or Vector, optional
        .. deprecated:: 1.12
            This argument is replaced by parent_point and will be removed in a
            future version.
            See :ref:`deprecated-mechanics-joint-pos` for more information.
    child_joint_pos : Point or Vector, optional
        .. deprecated:: 1.12
            This argument is replaced by child_point and will be removed in a
            future version.
            See :ref:`deprecated-mechanics-joint-pos` for more information.

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
        Attachment point where the joint is fixed to the parent body.
    child_point : Point
        Attachment point where the joint is fixed to the child body.
    parent_axis : Vector
        The axis fixed in the parent frame that represents the joint.
    child_axis : Vector
        The axis fixed in the child frame that represents the joint.
    parent_interframe : ReferenceFrame
        Intermediate frame of the parent body with respect to which the joint
        transformation is formulated.
    child_interframe : ReferenceFrame
        Intermediate frame of the child body with respect to which the joint
        transformation is formulated.
    kdes : list
        Kinematical differential equations of the joint.

    Examples
    =========

    A single prismatic joint is created from two bodies and has the following
    basic attributes:

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
    P_masscenter
    >>> joint.child_point
    C_masscenter
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

    To further demonstrate the use of the prismatic joint, the kinematics of two
    masses sliding, one moving relative to a fixed body and the other relative
    to the moving body. about the X axis of each connected body can be created
    as follows.

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

    def __init__(self, name, parent, child, coordinates=None, speeds=None,
                 parent_point=None, child_point=None, parent_axis=None,
                 child_axis=None, parent_interframe=None, child_interframe=None,
                 joint_axis=None, parent_joint_pos=None, child_joint_pos=None):

        self._joint_axis = joint_axis
        super().__init__(name, parent, child, coordinates, speeds, parent_point,
                         child_point, parent_axis, child_axis,
                         parent_interframe, child_interframe, parent_joint_pos,
                         child_joint_pos)

    def __str__(self):
        return (f'PrismaticJoint: {self.name}  parent: {self.parent}  '
                f'child: {self.child}')

    @property
    def joint_axis(self):
        """Axis along which the child translates with respect to the parent."""
        return self._joint_axis

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
        self._joint_axis = self._axis(
            self.joint_axis, self.parent_interframe, self.child_interframe)
        axis = self._express_joint_axis(self.parent_interframe)
        self.child_interframe.orient_axis(self.parent_interframe, axis, 0)

    def _set_angular_velocity(self):
        self.child_interframe.set_ang_vel(self.parent_interframe, 0)

    def _set_linear_velocity(self):
        self.child_point.set_pos(self.parent_point, self.coordinates[
            0] * self.joint_axis.normalize())
        self.parent_point.set_vel(self.parent.frame, 0)
        self.child_point.set_vel(self.child.frame, 0)
        self.child_point.set_vel(self.parent.frame, self.speeds[
            0] * self.joint_axis.normalize())
        self.child.masscenter.set_vel(self.parent.frame, self.speeds[
            0] * self.joint_axis.normalize())
